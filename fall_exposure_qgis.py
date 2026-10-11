# -*- coding: utf-8 -*-
"""
Absturzrisiko Wanderwege - QGIS Processing-Skript (v3: Kacheln + Threads)

Installation:
  Processing-Werkzeugkiste -> Python-Symbol -> "Skript zur Werkzeugkiste hinzufuegen..."
  diese Datei waehlen. Danach unter "Skripte > Wanderwege > Absturzrisiko Wanderwege (v3)".

Benoetigt nur numpy und GDAL (beides in QGIS enthalten).

Modell (angelehnt an Skitourenguru, Wiki a0077, eigene Umsetzung):
    risk = max(slope, k1 * fold, p_min) * v_max * forest_mod
  slope       Hangneigung am Weg [Grad]
  fold        Faltigkeit (Vector Ruggedness Measure, 0..1), mit k1 auf Grad skaliert
  v_max       max. Geschwindigkeit [m/s] eines Gleitblocks, der vom Wegrand der
              Falllinie folgt: h_kin += dz - mu*dl - h_kin*dl/drag_len
  forest_mod  fd < 60 % -> 1, sonst 1 - (fd-60)*0.5/40 (Walddichte entlang Sturzbahn)
  Zusatzregel: Hangneigung > 35 Grad -> mindestens Klasse 2.
  Klassen: 1 gering < T1 <= 2 maessig < T2 <= 3 hoch < T3 <= 4 sehr hoch
"""

import csv
import datetime
import math
import os
import platform
import threading
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from osgeo import gdal, osr

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsFeature,
    QgsFeatureSink,
    QgsField,
    QgsFieldConstraints,
    QgsFields,
    QgsGeometry,
    QgsPointXY,
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterString,
    QgsProcessingOutputNumber,
    QgsWkbTypes,
)

G = 9.81
CLASS_LABEL = {1: "gering", 2: "maessig", 3: "hoch", 4: "sehr hoch"}

# Name, Beschreibung, Standard, Typ (0=int, 1=double), Minimum
NUM_PARAMS = [
    ("RES", "Arbeitsaufloesung DEM [m]", 2.0, 1, 0.25),
    ("SEG_LEN", "Segmentlaenge Ausgabe [m]", 25.0, 1, 1.0),
    ("STEP", "Abstand Startpunkte entlang Weg [m]", 5.0, 1, 0.5),
    ("HALFWIDTH", "Wegtrasse beidseitig [m]", 1.5, 1, 0.0),
    ("MU", "Gleitreibung mu", 0.6, 1, 0.0),
    ("DRAG_LEN", "Bremslaenge v^2-Term [m] (0 = aus)", 100.0, 1, 0.0),
    ("MAX_PATH", "max. Laenge Sturzbahn [m]", 400.0, 1, 10.0),
    ("K1", "Skalierung Faltigkeit k1", 300.0, 1, 0.0),
    ("FOLD_WIN", "Fenster Faltigkeit [Zellen]", 5, 0, 3),
    ("P_MIN", "Untergrenze Wahrscheinlichkeit [Grad]", 10.0, 1, 0.0),
    ("FOREST_LEN", "Walddichte gemittelt ueber erste [m] der Sturzbahn", 50.0, 1, 0.0),
    ("STEEP_DEG", "Zusatzregel steiler Weg [Grad]", 35.0, 1, 0.0),
    ("CLIFF_DEG", "Absturzkante ab [Grad]", 60.0, 1, 0.0),
    ("NEAR_DIST", "'direkt neben dem Weg' [m]", 5.0, 1, 0.5),
    ("CHUNK_LEN", "lange Wege teilen alle [m]", 500.0, 1, 50.0),
    ("TILE_SIZE", "Kachelgroesse fuer gemeinsames DEM-Lesen [m]", 1000.0, 1, 100.0),
]


# --------------------------------------------------------------------------
# numpy-Helfer (Ersatz fuer scipy/shapely)
# --------------------------------------------------------------------------
def bilinear(arr, rows, cols):
    """Bilineare Interpolation; rows/cols in Pixelkoordinaten (Zellzentrum = i)."""
    h, w = arr.shape
    r = np.clip(rows, 0, h - 1)
    c = np.clip(cols, 0, w - 1)
    r0 = np.floor(r).astype(np.int64)
    c0 = np.floor(c).astype(np.int64)
    r1 = np.minimum(r0 + 1, h - 1)
    c1 = np.minimum(c0 + 1, w - 1)
    fr = r - r0
    fc = c - c0
    return ((arr[r0, c0] * (1 - fc) + arr[r0, c1] * fc) * (1 - fr)
            + (arr[r1, c0] * (1 - fc) + arr[r1, c1] * fc) * fr)


def box_filter(a, k):
    """Gleitender Mittelwert k x k (Rand: nearest)."""
    p = k // 2
    ap = np.pad(a, p, mode="edge")
    cs = np.cumsum(np.cumsum(ap, axis=0), axis=1)
    cs = np.pad(cs, ((1, 0), (1, 0)))
    s = cs[k:, k:] - cs[:-k, k:] - cs[k:, :-k] + cs[:-k, :-k]
    return s / float(k * k)


def vrm(z, res, k):
    gy, gx = np.gradient(z, res)
    gy = -gy
    norm = np.sqrt(gx * gx + gy * gy + 1.0)
    nx, ny, nz = -gx / norm, -gy / norm, 1.0 / norm
    k = max(3, int(k) | 1)
    sx, sy, sz = box_filter(nx, k), box_filter(ny, k), box_filter(nz, k)
    return 1.0 - np.sqrt(sx * sx + sy * sy + sz * sz)


class Polyline:
    """Minimaler Ersatz fuer shapely.LineString (interpolate / substring)."""

    def __init__(self, xy):
        self.xy = np.asarray(xy, dtype=float)
        seg = np.hypot(*np.diff(self.xy, axis=0).T)
        self.cum = np.concatenate([[0.0], np.cumsum(seg)])
        self.length = float(self.cum[-1])

    def interp(self, d):
        d = np.clip(np.asarray(d, dtype=float), 0, self.length)
        x = np.interp(d, self.cum, self.xy[:, 0])
        y = np.interp(d, self.cum, self.xy[:, 1])
        return np.stack([x, y], axis=-1)

    def substring(self, s0, s1):
        inner = (self.cum > s0) & (self.cum < s1)
        pts = np.vstack([self.interp(s0)[None], self.xy[inner], self.interp(s1)[None]])
        return pts

    def split(self, chunk_len):
        if self.length <= chunk_len * 1.5:
            return [self]
        n = int(math.ceil(self.length / chunk_len))
        e = np.linspace(0, self.length, n + 1)
        return [Polyline(self.substring(e[i], e[i + 1])) for i in range(n)]


_DS = threading.local()


def open_ds(path):
    """GDAL-Dataset pro Thread nur EINMAL oeffnen und wiederverwenden.
    Ein VRT der ganzen Schweiz (40'000+ Kacheln) braucht beim Oeffnen ~2 s, das
    Lesen eines Ausschnitts nur ~0.04 s. GDAL-Handles sind nicht threadsicher,
    deshalb ein eigenes Handle pro Thread."""
    cache = getattr(_DS, "cache", None)
    if cache is None:
        cache = _DS.cache = {}
    ds = cache.get(path)
    if ds is None:
        ds = gdal.Open(path)
        if ds is None:
            raise RuntimeError("GDAL konnte %s nicht oeffnen." % path)
        cache[path] = ds
    return ds


def close_ds():
    """Handles des aufrufenden Threads freigeben."""
    if getattr(_DS, "cache", None):
        _DS.cache.clear()


def warp_grid(path, wkt, origin, size, res, alg):
    """Raster per gdal.Warp exakt auf ein Ziel-Grid (CRS wkt, Ursprung oben links,
    Groesse w x h, Zellgroesse res) in den Speicher rechnen."""
    (ox, oy), (w, h) = origin, size
    ds = gdal.Warp("", open_ds(path), format="MEM", dstSRS=wkt,
                   outputBounds=(ox, oy - h * res, ox + w * res, oy), width=w, height=h,
                   resampleAlg=alg, dstNodata=-99999, outputType=gdal.GDT_Float32)
    if ds is None:
        raise RuntimeError("GDAL konnte %s nicht lesen." % path)
    arr = ds.GetRasterBand(1).ReadAsArray().astype("float32")
    arr[arr == -99999] = np.nan
    ds = None
    return arr


_SAME_CRS = {}


def same_crs(path, wkt):
    """True, wenn das Raster bereits im Ziel-CRS liegt und nordausgerichtet ist."""
    key = (path, wkt)
    if key not in _SAME_CRS:
        ok = False
        ds = gdal.Open(path)
        if ds is not None:
            gt = ds.GetGeoTransform()
            src = osr.SpatialReference()
            dst = osr.SpatialReference()
            try:
                src.ImportFromWkt(ds.GetProjection())
                dst.ImportFromWkt(wkt)
                ok = bool(src.IsSame(dst)) and gt[2] == 0 and gt[4] == 0 and gt[5] < 0
            except Exception:
                ok = False
        ds = None
        _SAME_CRS[key] = ok
    return _SAME_CRS[key]


def read_dem(path, wkt, bounds, res):
    """DEM-Ausschnitt mit Zellgroesse res. Liegt das DEM schon im Ziel-CRS (z.B.
    swissALTI3D in EPSG:2056), wird ohne Reprojektion direkt gelesen und gemittelt
    (deutlich schneller); sonst ueber gdal.Warp.
    Rueckgabe: (array, (ox, oy), (w, h)) mit Ursprung oben links."""
    minx, miny, maxx, maxy = bounds
    if same_crs(path, wkt):
        ds = open_ds(path)
        gt = ds.GetGeoTransform()
        px, py = gt[1], -gt[5]
        f = max(1, int(round(res / px)))           # Quellzellen pro Arbeitszelle
        cell = f * px                               # tatsaechliche Arbeitszellgroesse
        # Fenster am Quellraster ausrichten -> ganzzahlige Pixeloffsets
        c0 = int(math.floor((minx - gt[0]) / px / f)) * f
        r0 = int(math.floor((gt[3] - maxy) / py / f)) * f
        w = max(2, int(math.ceil((maxx - minx) / cell)) + 1)
        h = max(2, int(math.ceil((maxy - miny) / cell)) + 1)
        arr = np.full((h, w), np.nan, dtype="float32")
        X, Y = ds.RasterXSize, ds.RasterYSize
        # Schnitt mit der Rasterausdehnung (in Arbeitszellen)
        bc0 = max(0, int(math.ceil(-c0 / f)))
        br0 = max(0, int(math.ceil(-r0 / f)))
        bc1 = min(w, (X - c0) // f)
        br1 = min(h, (Y - r0) // f)
        if bc1 > bc0 and br1 > br0:
            band = ds.GetRasterBand(1)
            sub = band.ReadAsArray(c0 + bc0 * f, r0 + br0 * f,
                                   (bc1 - bc0) * f, (br1 - br0) * f,
                                   buf_xsize=bc1 - bc0, buf_ysize=br1 - br0,
                                   resample_alg=gdal.GRIORA_Average).astype("float32")
            nod = band.GetNoDataValue()
            if nod is not None:
                sub[sub == nod] = np.nan
            arr[br0:br1, bc0:bc1] = sub
        return arr, (gt[0] + c0 * px, gt[3] - r0 * py), (w, h), cell
    w = max(2, int(math.ceil((maxx - minx) / res)))
    h = max(2, int(math.ceil((maxy - miny) / res)))
    origin = (minx, miny + h * res)
    return warp_grid(path, wkt, origin, (w, h), res, "average"), origin, (w, h), res


def hide_password(text):
    """Passwoerter aus Verbindungszeichenfolgen entfernen, bevor sie ins Protokoll gehen."""
    import re
    text = re.sub(r"(password\s*=\s*)('[^']*'|\S+)", r"\1'***'", text, flags=re.I)
    return re.sub(r"(://[^:/@\s]+:)[^@\s]+@", r"\1***@", text)


def fmt_duration(s):
    s = int(round(s))
    h, r = divmod(s, 3600)
    m, s = divmod(r, 60)
    return "%d h %02d min %02d s" % (h, m, s) if h else "%d min %02d s" % (m, s)


def cpu_name():
    """Prozessorname moeglichst lesbar (Windows, macOS, Linux)."""
    try:
        sysname = platform.system()
        if sysname == "Darwin":
            import subprocess
            out = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                 capture_output=True, text=True, timeout=5).stdout.strip()
            if out:
                return out
        elif sysname == "Windows":
            import winreg
            k = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                               r"HARDWARE\DESCRIPTION\System\CentralProcessor\0")
            return winreg.QueryValueEx(k, "ProcessorNameString")[0].strip()
        elif sysname == "Linux":
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return platform.processor() or platform.machine()


# --------------------------------------------------------------------------
# Kern
# --------------------------------------------------------------------------
def load_tile(dem_path, forest_path, wkt, bounds, p):
    """Liest DEM (und Wald) fuer ein Rechteck und berechnet die abgeleiteten Raster
    einmal. Wird von allen Wegstuecken einer Kachel gemeinsam benutzt."""
    z, (ox, oy), size, res = read_dem(dem_path, wkt, bounds, p["RES"])
    if np.isnan(z).all():
        return None
    z = np.where(np.isnan(z), np.nanmean(z), z)
    gy, gx = np.gradient(z, res)
    gy = -gy
    T = dict(z=z, gx=gx, gy=gy, ox=ox, oy=oy, res=res,
             slope=np.degrees(np.arctan(np.hypot(gx, gy))),
             fold=vrm(z, res, p["FOLD_WIN"]), fd=None)
    if forest_path:
        fd = warp_grid(forest_path, wkt, (ox, oy), size, res, "bilinear")
        T["fd"] = np.nan_to_num(fd, nan=0.0).clip(0, 100)
    return T


def line_points(line, p):
    """Startpunkte (alle STEP m), deren Laufdistanz und Normalenrichtung."""
    L = line.length
    dists = np.arange(p["STEP"] / 2.0, L, p["STEP"])
    if dists.size == 0:
        dists = np.array([L / 2.0])
    pts = line.interp(dists)
    eps = min(2.0, L / 4.0)
    t = line.interp(dists + eps) - line.interp(dists - eps)
    t /= np.linalg.norm(t, axis=1, keepdims=True).clip(min=1e-9)
    nrm = np.stack([-t[:, 1], t[:, 0]], axis=1)
    return dists, pts, nrm


def evaluate_points(T, pts, nrm, p, thr):
    """Kernrechnung fuer beliebig viele Startpunkte gleichzeitig (vektorisiert).
    Liefert pro Punkt: Klasse, risk und Kennzahlen."""
    res, ox, oy = T["res"], T["ox"], T["oy"]
    z, gx, gy, slope_r, fold_r, fd_r = T["z"], T["gx"], T["gy"], T["slope"], T["fold"], T["fd"]
    hw = p["HALFWIDTH"]

    def samp(arr, x, y):
        return bilinear(arr, (oy - np.asarray(y)) / res - 0.5, (np.asarray(x) - ox) / res - 0.5)

    # Wahrscheinlichkeitsterm
    offs = np.linspace(-hw - res / 2, hw + res / 2, 5)
    qx = pts[:, 0, None] + offs[None, :] * nrm[:, 0, None]
    qy = pts[:, 1, None] + offs[None, :] * nrm[:, 1, None]
    slope_w = samp(slope_r, qx, qy).max(axis=1)
    fold_w = samp(fold_r, qx, qy).max(axis=1)
    z0 = np.median(samp(z, qx, qy), axis=1)

    # Absturzkante direkt am Weg
    nd = np.arange(hw, p["NEAR_DIST"] + 1e-6, max(res / 2, 0.5))
    imm = np.zeros(len(pts), bool)
    for sgn in (1, -1):
        ix = pts[:, 0, None] + sgn * nd[None, :] * nrm[:, 0, None]
        iy = pts[:, 1, None] + sgn * nd[None, :] * nrm[:, 1, None]
        zz = samp(z, ix, iy)
        imm |= (samp(slope_r, ix, iy) >= p["CLIFF_DEG"]).any(axis=1) & \
               ((z0[:, None] - zz).max(axis=1) > 2.0)

    # Virtuelles Fallexperiment (beide Seiten, alle Punkte gleichzeitig)
    x = np.concatenate([pts[:, 0] + hw * nrm[:, 0], pts[:, 0] - hw * nrm[:, 0]])
    y = np.concatenate([pts[:, 1] + hw * nrm[:, 1], pts[:, 1] - hw * nrm[:, 1]])
    zprev = np.concatenate([z0, z0])
    z_start = zprev.copy()
    hk = np.zeros_like(x)
    vmax = np.zeros_like(x)
    drop_v = np.zeros_like(x)
    fd_sum = np.zeros_like(x)
    fd_n = np.zeros_like(x)
    alive = np.ones_like(x, bool)
    idx = np.arange(x.size)            # Indizes der noch laufenden Bloecke
    for i in range(int(p["MAX_PATH"] / res) + 1):
        # nur noch laufende Bloecke rechnen -> wird schnell klein
        xa, ya = x[idx], y[idx]
        if i > 0:
            ggx, ggy = samp(gx, xa, ya), samp(gy, xa, ya)
            gm = np.hypot(ggx, ggy)
            ok = gm > 1e-3
            d = np.where(ok, res / np.maximum(gm, 1e-9), 0.0)
            xa = xa - ggx * d
            ya = ya - ggy * d
            x[idx], y[idx] = xa, ya
            dl = np.where(ok, res, 0.0)
        else:
            ok = np.ones(idx.size, bool)
            dl = np.full(idx.size, hw)
        zc = samp(z, xa, ya)
        h = hk[idx]
        drag = h * dl / p["DRAG_LEN"] if p["DRAG_LEN"] > 0 else 0.0
        h = np.where(ok, h + (zprev[idx] - zc) - p["MU"] * dl - drag, h)
        stopped = (~ok) | ((h <= 0) & (i > 0))
        h = np.maximum(h, 0.0)
        hk[idx] = h
        v = np.sqrt(2 * G * h)
        better = v > vmax[idx]
        bi = idx[better]
        vmax[bi] = v[better]
        drop_v[bi] = z_start[bi] - zc[better]
        if fd_r is not None and i * res <= p["FOREST_LEN"]:
            m = ~stopped | (i == 0)
            fd_sum[idx] += np.where(m, samp(fd_r, xa, ya), 0.0)
            fd_n[idx] += m
        zprev[idx] = zc
        if i > 0:
            idx = idx[~stopped]
            if idx.size == 0:
                break

    n = len(pts)
    if fd_r is not None:
        fd = fd_sum / np.maximum(fd_n, 1)
        fmod = np.where(fd < 60, 1.0, 1.0 - (fd - 60) * 0.5 / 40)
    else:
        fd = np.zeros_like(vmax)
        fmod = np.ones_like(vmax)

    prob = np.maximum(np.maximum(slope_w, p["K1"] * fold_w), p["P_MIN"])
    rs = np.concatenate([prob, prob]) * vmax * fmod
    right = rs[n:] > rs[:n]
    pick = lambda a: np.where(right, a[n:], a[:n])
    risk, v_pt, drop_pt, fmod_pt, fd_pt = (pick(a) for a in (rs, vmax, drop_v, fmod, fd))

    cls = np.digitize(risk, thr).astype(int) + 1
    cls = np.where(slope_w > p["STEEP_DEG"], np.maximum(cls, 2), cls)
    return dict(cls=cls, risk=risk, slope=slope_w, fold=fold_w, v=v_pt, drop=drop_pt,
                fmod=fmod_pt, fd=fd_pt, imm=imm)


def make_segments(line, dists, E, p):
    """Punkte eines Wegstuecks zu SEG_LEN-Segmenten zusammenfassen."""
    L = line.length
    out = []
    nseg = max(1, int(math.ceil(L / p["SEG_LEN"])))
    edges = np.linspace(0, L, nseg + 1)
    for i in range(nseg):
        m = (dists >= edges[i]) & (dists < edges[i + 1])
        if not m.any():
            continue
        risk = E["risk"]
        o = np.argsort(risk[m])
        j = np.flatnonzero(m)[o[-2] if m.sum() >= 4 else o[-1]]
        c = np.sort(E["cls"][m])
        seg_cls = int(c[-2]) if c.size >= 4 else int(c[-1])
        out.append((line.substring(edges[i], edges[i + 1]), dict(
            exp_class=seg_cls,
            exp_label=CLASS_LABEL[seg_cls],
            risk=float(round(risk[j])),
            slope_deg=round(float(E["slope"][m].max()), 1),
            fold=round(float(E["fold"][m].max()), 3),
            v_max=round(float(E["v"][j]), 1),
            drop_m=round(float(E["drop"][j]), 1),
            forest_pct=float(round(E["fd"][j])),
            forest_mod=round(float(E["fmod"][j]), 2),
            immediate_drop=int(E["imm"][m].sum() >= (2 if m.sum() >= 4 else 1)),
        )))
    return out


def tile_reach(p):
    """Rand um eine Kachel: halbe max. Stuecklaenge + Sturzbahn + Sicherheit."""
    return 0.75 * p["CHUNK_LEN"] + p["MAX_PATH"] + 4 * p["RES"] + 20


def process_tile(chunks, bounds, dem_path, forest_path, wkt, p, thr):
    """Alle Wegstuecke einer Kachel: DEM einmal lesen, alle Punkte in einem
    Durchgang rechnen. Rueckgabe: Liste (task_info, segmente) pro Stueck."""
    T = load_tile(dem_path, forest_path, wkt, bounds, p)
    if T is None:
        return [(info, []) for info, _ in chunks]
    per = []
    for info, line in chunks:
        if line.length < 1.0:
            per.append((info, line, None, None, None))
            continue
        dists, pts, nrm = line_points(line, p)
        per.append((info, line, dists, pts, nrm))
    valid = [q for q in per if q[2] is not None]
    if not valid:
        return [(q[0], []) for q in per]
    E = evaluate_points(T, np.concatenate([q[3] for q in valid]),
                        np.concatenate([q[4] for q in valid]), p, thr)
    out, k = [], 0
    for info, line, dists, pts, nrm in per:
        if dists is None:
            out.append((info, []))
            continue
        n = len(dists)
        Ei = {key: val[k:k + n] for key, val in E.items()}
        k += n
        out.append((info, make_segments(line, dists, Ei, p)))
    return out


def process_line(line, dem_path, forest_path, wkt, p, thr):
    """Ein einzelnes Wegstueck (eigener DEM-Ausschnitt). Fuer Tests/Einzelaufrufe."""
    if line.length < 1.0:
        return []
    r = tile_reach(p)
    xy = line.xy
    bounds = (*(xy.min(axis=0) - r), *(xy.max(axis=0) + r))
    return process_tile([(None, line)], bounds, dem_path, forest_path, wkt, p, thr)[0][1]


OUT_FIELDS = [
    ("exp_class", QVariant.Int), ("exp_label", QVariant.String),
    ("risk", QVariant.Double), ("slope_deg", QVariant.Double), ("fold", QVariant.Double),
    ("v_max", QVariant.Double), ("drop_m", QVariant.Double),
    ("forest_pct", QVariant.Double), ("forest_mod", QVariant.Double),
    ("immediate_drop", QVariant.Int),
]


# --------------------------------------------------------------------------
# Processing-Algorithmus
# --------------------------------------------------------------------------
class FallExposureAlgorithm(QgsProcessingAlgorithm):
    INPUT = "INPUT"
    DEM = "DEM"
    FOREST = "FOREST"
    THR = "THR"
    THREADS = "THREADS"
    RUNLOG = "RUNLOG"
    RUNTIME = "RUNTIME_S"
    OUTPUT = "OUTPUT"

    def tr(self, s):
        return QCoreApplication.translate("Processing", s)

    def createInstance(self):
        return FallExposureAlgorithm()

    def name(self):
        return "absturzrisiko_wanderwege_v3"

    def displayName(self):
        return self.tr("Absturzrisiko Wanderwege (v3)")

    def group(self):
        return self.tr("Wanderwege")

    def groupId(self):
        return "wanderwege"

    def shortHelpString(self):
        return self.tr(
            "Berechnet das Absturzrisiko entlang von Linien (z.B. sac_scale-Wege) aus einem "
            "Hoehenmodell (z.B. swissALTI3D als VRT) und optional einem Walddichte-Raster "
            "(0-100 %, z.B. Copernicus Tree Cover Density).<br><br>"
            "risk = max(Hangneigung, k1*Faltigkeit, p_min) * v_max * Waldfaktor<br>"
            "v_max aus einem Gleitblock, der vom Wegrand der Falllinie folgt.<br><br>"
            "Ausgabe: Segmente mit exp_class 1-4 (gering/maessig/hoch/sehr hoch) im CRS des DEM. "
            "Die Klassengrenzen sind Startwerte und muessen an bekannten Wegen kalibriert werden.<br><br>"
            "Tipp: 'Nur gewaehlte Objekte' nutzen, um zuerst einen Ausschnitt zu testen. "
            "Das Ergebnis ist ein Modell und ersetzt keine Beurteilung im Gelaende.")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT, self.tr("Wege (Linien)"), [QgsProcessing.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.DEM, self.tr("Hoehenmodell (metrisches CRS, z.B. EPSG:2056)")))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.FOREST, self.tr("Walddichte 0-100 % (optional)"), optional=True))
        self.addParameter(QgsProcessingParameterString(
            self.THR, self.tr("Klassengrenzen risk (T1,T2,T3)"), "200,600,1500"))
        self.addParameter(QgsProcessingParameterNumber(
            self.THREADS, self.tr("Anzahl paralleler Threads"),
            QgsProcessingParameterNumber.Integer, max(1, os.cpu_count() or 1), False, 1, 64))
        for key, desc, default, typ, minimum in NUM_PARAMS:
            prm = QgsProcessingParameterNumber(
                key, self.tr(desc),
                QgsProcessingParameterNumber.Integer if typ == 0 else QgsProcessingParameterNumber.Double,
                default, False, minimum)
            prm.setFlags(prm.flags() | prm.FlagAdvanced)
            self.addParameter(prm)
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT, self.tr("Absturzrisiko"), QgsProcessing.TypeVectorLine))
        prm = QgsProcessingParameterFileDestination(
            self.RUNLOG, self.tr("Laufzeit-Protokoll (CSV, leer = Projektordner)"),
            "CSV (*.csv)", optional=True, createByDefault=False)
        prm.setHelp(self.tr("Leer lassen: CSV im Ordner des QGIS-Projekts; ohne gespeichertes "
                            "Projekt neben der Ausgabedatei, sonst im Benutzerordner."))
        prm.setFlags(prm.flags() | prm.FlagAdvanced)
        self.addParameter(prm)
        self.addOutput(QgsProcessingOutputNumber(self.RUNTIME, self.tr("Laufzeit [s]")))

    def processAlgorithm(self, parameters, context, feedback):
        t_start = time.time()
        source = self.parameterAsSource(parameters, self.INPUT, context)
        dem = self.parameterAsRasterLayer(parameters, self.DEM, context)
        forest = self.parameterAsRasterLayer(parameters, self.FOREST, context)
        if source is None or dem is None:
            raise QgsProcessingException("Eingabe fehlt.")
        dem_crs = dem.crs()
        if dem_crs.isGeographic():
            raise QgsProcessingException("Das DEM braucht ein metrisches CRS (z.B. EPSG:2056).")
        try:
            thr = np.array([float(v) for v in
                            self.parameterAsString(parameters, self.THR, context).split(",")])
        except ValueError:
            thr = np.array([])
        if thr.size != 3 or np.any(np.diff(thr) <= 0):
            raise QgsProcessingException("Klassengrenzen: 3 aufsteigende Zahlen, z.B. 200,600,1500")
        p = {}
        for key, _, _, typ, _ in NUM_PARAMS:
            p[key] = (self.parameterAsInt if typ == 0 else self.parameterAsDouble)(
                parameters, key, context)

        dem_path = dem.source()
        forest_path = forest.source() if forest else None
        wkt = dem_crs.toWkt(QgsCoordinateReferenceSystem.WKT_PREFERRED_GDAL)

        # Quellfelder uebernehmen. Typische Schluesselspalten (fid, id, gid, ogc_fid)
        # werden als src_<name> kopiert: Ein Weg ergibt mehrere Segmente, der Wert
        # waere in der Ausgabe nicht eindeutig und kollidiert sonst mit deren
        # Primaerschluessel (GPKG fid, PostGIS id).
        # Eigene, fortlaufende seg_id als ERSTES Feld: QGIS nimmt beim Export nach
        # PostGIS die erste Integer-Spalte als Primaerschluessel.
        key_names = {"fid", "id", "gid", "ogc_fid"}
        own_names = {"seg_id"} | {n for n, _ in OUT_FIELDS}
        fields = QgsFields()
        fields.append(QgsField("seg_id", QVariant.LongLong))
        src_idx = []      # Index des Feldes in der Quelle
        skipped = []
        for i, f in enumerate(source.fields()):
            name = f.name()
            low = name.lower()
            if low in key_names:
                name = "src_" + low
            if name.lower() in own_names:
                # Eingabe ist z.B. ein frueheres Ergebnis: alte Werte werden neu berechnet
                skipped.append(f.name())
                continue
            # Doppelte Namen (auch nur in Gross-/Kleinschreibung) eindeutig machen;
            # QgsFields.append verwirft Duplikate sonst stillschweigend und alle
            # folgenden Werte rutschen in die falsche Spalte.
            base, k = name, 1
            while fields.lookupField(name) >= 0:
                name = "%s_%d" % (base, k)
                k += 1
            nf = QgsField(f)
            # Constraints (NOT NULL/UNIQUE eines Quell-Primaerschluessels) nicht
            # uebernehmen, sonst macht der PostGIS-Export dieses Feld zum PK.
            nf.setConstraints(QgsFieldConstraints())
            nf.setName(name)
            if not fields.append(nf):
                raise QgsProcessingException("Feld %s konnte nicht angelegt werden." % name)
            src_idx.append(i)
        for name, typ in OUT_FIELDS:
            fields.append(QgsField(name, typ))
        if skipped:
            feedback.pushInfo("Hinweis: Felder %s stammen offenbar aus einem frueheren Lauf "
                              "und werden neu berechnet." % ", ".join(skipped))
        n_src = len(src_idx)
        out_idx = {name: fields.indexOf(name) for name, _ in OUT_FIELDS}
        assert fields.count() == 1 + n_src + len(OUT_FIELDS)
        sink, dest_id = self.parameterAsSink(parameters, self.OUTPUT, context, fields,
                                             QgsWkbTypes.LineString, dem_crs)
        if sink is None:
            raise QgsProcessingException("Ausgabe konnte nicht angelegt werden.")

        xform = QgsCoordinateTransform(source.sourceCrs(), dem_crs, context.transformContext())
        total = source.featureCount() or 1
        n_threads = self.parameterAsInt(parameters, self.THREADS, context)

        def tasks():
            """Liest die Wege im Hauptthread (QGIS-Objekte sind nicht threadsicher)
            und liefert reine numpy-Arbeitspakete: (Weg-Nr, Feature-ID, Attribute, Stueck)."""
            for k, feat in enumerate(source.getFeatures()):
                if feedback.isCanceled():
                    return
                geom = feat.geometry()
                if geom is None or geom.isEmpty():
                    continue
                geom = QgsGeometry(geom)
                geom.transform(xform)
                parts = geom.asMultiPolyline() if geom.isMultipart() else [geom.asPolyline()]
                all_attrs = feat.attributes()
                attrs_src = [all_attrs[i] for i in src_idx]
                for part in parts:
                    if len(part) < 2:
                        continue
                    pl = Polyline([(pt.x(), pt.y()) for pt in part])
                    for chunk in pl.split(p["CHUNK_LEN"]):
                        yield k, feat.id(), attrs_src, chunk

        # 1) Alle Wegstuecke einlesen (Hauptthread) und nach Kacheln gruppieren.
        #    Pro Kachel wird das DEM nur einmal gelesen -> viel schneller bei den
        #    vielen kurzen OSM-Wegen; die Kacheln laufen parallel in Threads.
        tile = p["TILE_SIZE"]
        tiles = {}
        n_chunks = 0
        total_m = 0.0
        for k, fid, attrs_src, chunk in tasks():
            mx, my = chunk.interp(chunk.length / 2.0)
            key = (int(math.floor(mx / tile)), int(math.floor(my / tile)))
            tiles.setdefault(key, []).append(((k, fid, attrs_src), chunk))
            n_chunks += 1
            total_m += chunk.length
        if feedback.isCanceled():
            return {self.OUTPUT: dest_id}
        t_read = time.time()
        reach = tile_reach(p)
        keys = sorted(tiles)
        feedback.pushInfo("%d Wegstuecke in %d Kacheln, %d Thread(s)."
                          % (n_chunks, len(keys), n_threads))

        def compute(key):
            # laeuft in Worker-Threads: nur numpy und GDAL, keine QGIS-Objekte
            tx, ty = key
            bounds = (tx * tile - reach, ty * tile - reach,
                      (tx + 1) * tile + reach, (ty + 1) * tile + reach)
            try:
                return process_tile(tiles[key], bounds, dem_path, forest_path, wkt, p, thr), None
            except Exception as e:  # eine Kachel darf den Lauf nicht abbrechen
                return [(info, []) for info, _ in tiles[key]], "Kachel %s: %s" % (key, e)

        nseg = 0
        done = 0

        def write(result):
            # schreibt im Hauptthread
            nonlocal nseg, done
            rows_per_chunk, err = result
            if err:
                feedback.reportError(err)
            for (k, fid, attrs_src), rows in rows_per_chunk:
                for xy, vals in rows:
                    of = QgsFeature(fields)
                    of.setGeometry(QgsGeometry.fromPolylineXY(
                        [QgsPointXY(float(a), float(b)) for a, b in xy]))
                    nseg += 1
                    attrs = [nseg] + list(attrs_src) + [None] * len(OUT_FIELDS)
                    for name, _ in OUT_FIELDS:
                        attrs[out_idx[name]] = vals[name]
                    of.setAttributes(attrs)
                    if not sink.addFeature(of, QgsFeatureSink.FastInsert):
                        raise QgsProcessingException(
                            "Schreiben in die Ausgabe fehlgeschlagen (Segment %d, Weg %s): %s"
                            % (nseg, fid, sink.lastError()))
            done += 1
            feedback.setProgress(int(100 * done / max(1, len(keys))))

        if n_threads <= 1:
            for key in keys:
                if feedback.isCanceled():
                    break
                write(compute(key))
        else:
            # Begrenzte Warteschlange: haelt den Speicher klein
            pending = deque()
            with ThreadPoolExecutor(max_workers=n_threads, initializer=close_ds) as ex:
                for key in keys:
                    if feedback.isCanceled():
                        break
                    pending.append(ex.submit(compute, key))
                    while len(pending) >= 2 * n_threads:
                        write(pending.popleft().result())
                while pending:
                    fut = pending.popleft()
                    if feedback.isCanceled():
                        fut.cancel()
                        continue
                    write(fut.result())
        close_ds()
        # Gepufferte Inserts (z.B. PostGIS) erst hier wirklich geschrieben ->
        # Fehler hier pruefen, statt sie still zu verlieren.
        if hasattr(sink, "flushBuffer") and not sink.flushBuffer():
            raise QgsProcessingException(
                "Schreiben in die Ausgabe fehlgeschlagen: %s" % sink.lastError())
        feedback.pushInfo("%d Segmente geschrieben." % nseg)

        # ---- Laufzeit protokollieren ----
        t_end = time.time()
        total_s = t_end - t_start
        km = total_m / 1000.0
        summary = {
            "zeitpunkt": datetime.datetime.now().isoformat(timespec="seconds"),
            "rechner": platform.node(),
            "system": "%s %s" % (platform.system(), platform.release()),
            "prozessor": cpu_name(),
            "kerne_logisch": os.cpu_count(),
            "threads": n_threads,
            "abgebrochen": int(feedback.isCanceled()),
            "wegstuecke": n_chunks,
            "kacheln": len(keys),
            "segmente": nseg,
            "weg_km": round(km, 1),
            "zeit_einlesen_s": round(t_read - t_start, 1),
            "zeit_rechnen_s": round(t_end - t_read, 1),
            "zeit_total_s": round(total_s, 1),
            "km_pro_s": round(km / total_s, 2) if total_s > 0 else None,
            "dem": dem_path,
            "wald": forest_path or "",
            "eingabe": source.sourceName(),
            "ausgabe": hide_password(str(dest_id)),
            "klassengrenzen": ",".join("%g" % v for v in thr),
            "parameter": " ".join("%s=%g" % (k, v) for k, v in p.items()),
        }
        feedback.pushInfo("Laufzeit: %s (Einlesen %.0f s, Rechnen %.0f s), %.0f km Weg, "
                          "%.1f km/s, %d Thread(s) auf %s"
                          % (fmt_duration(total_s), t_read - t_start, t_end - t_read,
                             km, summary["km_pro_s"] or 0, n_threads, summary["prozessor"]))
        log_path = self.parameterAsFileOutput(parameters, self.RUNLOG, context)
        if not log_path:
            # 1. Ordner des gespeicherten QGIS-Projekts
            proj = context.project() if hasattr(context, "project") else None
            proj_file = proj.fileName() if proj is not None else ""
            if proj_file and not proj_file.lower().startswith(("postgres", "geopackage:")):
                proj_dir = os.path.dirname(os.path.abspath(proj_file))
                if os.path.isdir(proj_dir):
                    log_path = os.path.join(proj_dir, "fall_exposure_laufzeiten.csv")
        if not log_path:
            # 2. Ordner der Ausgabedatei, 3. Benutzerordner
            out = parameters.get(self.OUTPUT)
            out = getattr(out, "sink", out)
            out = out.staticValue() if hasattr(out, "staticValue") else out
            out = str(out or "")
            base = out.split("|")[0]
            if base and not base.upper().startswith(("PG:", "POSTGRES", "MEMORY", "TEMPORARY")) \
                    and os.path.isdir(os.path.dirname(base) or "."):
                log_path = os.path.join(os.path.dirname(base) or ".", "fall_exposure_laufzeiten.csv")
            else:
                log_path = os.path.join(os.path.expanduser("~"), "fall_exposure_laufzeiten.csv")
        try:
            new = not os.path.exists(log_path)
            with open(log_path, "a", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(summary), delimiter=";")
                if new:
                    w.writeheader()
                w.writerow(summary)
            feedback.pushInfo("Laufzeit protokolliert in: %s" % log_path)
        except Exception as e:
            feedback.reportError("Laufzeit-Protokoll konnte nicht geschrieben werden (%s): %s"
                                 % (log_path, e))
        return {self.OUTPUT: dest_id, self.RUNTIME: round(total_s, 1)}
