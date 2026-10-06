# -*- coding: utf-8 -*-
"""
Absturzrisiko Wanderwege - QGIS Processing-Skript (v2)

Installation:
  Processing-Werkzeugkiste -> Python-Symbol -> "Skript zur Werkzeugkiste hinzufuegen..."
  diese Datei waehlen. Danach unter "Skripte > Wanderwege > Absturzrisiko Wanderwege".

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

import math

import numpy as np
from osgeo import gdal

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
    QgsProcessingParameterString,
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


def warp_grid(path, wkt, bounds, res, alg):
    """Raster auf Ziel-Grid (CRS wkt, bounds, res) in den Speicher warpen."""
    minx, miny, maxx, maxy = bounds
    w = max(2, int(math.ceil((maxx - minx) / res)))
    h = max(2, int(math.ceil((maxy - miny) / res)))
    maxx, maxy = minx + w * res, miny + h * res
    ds = gdal.Warp("", path, format="MEM", dstSRS=wkt,
                   outputBounds=(minx, miny, maxx, maxy), width=w, height=h,
                   resampleAlg=alg, dstNodata=-99999, outputType=gdal.GDT_Float32)
    if ds is None:
        raise QgsProcessingException("GDAL konnte %s nicht lesen." % path)
    arr = ds.GetRasterBand(1).ReadAsArray().astype("float32")
    arr[arr == -99999] = np.nan
    ds = None
    return arr, (minx, maxy, res)          # Ursprung oben links


# --------------------------------------------------------------------------
# Kern
# --------------------------------------------------------------------------
def process_line(line, dem_path, forest_path, wkt, p, thr):
    L = line.length
    if L < 1.0:
        return []
    res = p["RES"]
    hw = p["HALFWIDTH"]

    dists = np.arange(p["STEP"] / 2.0, L, p["STEP"])
    if dists.size == 0:
        dists = np.array([L / 2.0])
    pts = line.interp(dists)
    eps = min(2.0, L / 4.0)
    t = line.interp(dists + eps) - line.interp(dists - eps)
    t /= np.linalg.norm(t, axis=1, keepdims=True).clip(min=1e-9)
    nrm = np.stack([-t[:, 1], t[:, 0]], axis=1)

    reach = p["MAX_PATH"] + 4 * res + 10
    bounds = (*(pts.min(axis=0) - reach), *(pts.max(axis=0) + reach))
    z, (ox, oy, _) = warp_grid(dem_path, wkt, bounds, res, "average")
    if np.isnan(z).all():
        return []
    z = np.where(np.isnan(z), np.nanmean(z), z)
    gy, gx = np.gradient(z, res)
    gy = -gy
    slope_r = np.degrees(np.arctan(np.hypot(gx, gy)))
    fold_r = vrm(z, res, p["FOLD_WIN"])
    fd_r = None
    if forest_path:
        fd_r, _ = warp_grid(forest_path, wkt, bounds, res, "bilinear")
        fd_r = np.nan_to_num(fd_r, nan=0.0).clip(0, 100)

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

    # Virtuelles Fallexperiment
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
    for i in range(int(p["MAX_PATH"] / res) + 1):
        if i > 0:
            ggx, ggy = samp(gx, x, y), samp(gy, x, y)
            gm = np.hypot(ggx, ggy)
            alive &= gm > 1e-3
            d = np.where(alive, res / np.maximum(gm, 1e-9), 0.0)
            x = x - ggx * d
            y = y - ggy * d
            dl = np.where(alive, res, 0.0)
        else:
            dl = np.full_like(x, hw)
        zc = samp(z, x, y)
        drag = hk * dl / p["DRAG_LEN"] if p["DRAG_LEN"] > 0 else 0.0
        hk = np.where(alive, hk + (zprev - zc) - p["MU"] * dl - drag, hk)
        stopped = alive & (hk <= 0) & (i > 0)
        hk = np.maximum(hk, 0.0)
        alive &= ~stopped
        v = np.sqrt(2 * G * hk)
        better = v > vmax
        vmax = np.where(better, v, vmax)
        drop_v = np.where(better, z_start - zc, drop_v)
        if fd_r is not None and i * res <= p["FOREST_LEN"]:
            m = alive | (i == 0)
            fd_sum += np.where(m, samp(fd_r, x, y), 0.0)
            fd_n += m
        zprev = zc
        if i > 0 and not alive.any():
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

    out = []
    nseg = max(1, int(math.ceil(L / p["SEG_LEN"])))
    edges = np.linspace(0, L, nseg + 1)
    for i in range(nseg):
        m = (dists >= edges[i]) & (dists < edges[i + 1])
        if not m.any():
            continue
        o = np.argsort(risk[m])
        j = np.flatnonzero(m)[o[-2] if m.sum() >= 4 else o[-1]]
        c = np.sort(cls[m])
        seg_cls = int(c[-2]) if c.size >= 4 else int(c[-1])
        out.append((line.substring(edges[i], edges[i + 1]), dict(
            exp_class=seg_cls,
            exp_label=CLASS_LABEL[seg_cls],
            risk=float(round(risk[j])),
            slope_deg=round(float(slope_w[m].max()), 1),
            fold=round(float(fold_w[m].max()), 3),
            v_max=round(float(v_pt[j]), 1),
            drop_m=round(float(drop_pt[j]), 1),
            forest_pct=float(round(fd_pt[j])),
            forest_mod=round(float(fmod_pt[j]), 2),
            immediate_drop=int(imm[m].sum() >= (2 if m.sum() >= 4 else 1)),
        )))
    return out


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
    OUTPUT = "OUTPUT"

    def tr(self, s):
        return QCoreApplication.translate("Processing", s)

    def createInstance(self):
        return FallExposureAlgorithm()

    def name(self):
        return "absturzrisiko_wanderwege"

    def displayName(self):
        return self.tr("Absturzrisiko Wanderwege")

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
        for key, desc, default, typ, minimum in NUM_PARAMS:
            prm = QgsProcessingParameterNumber(
                key, self.tr(desc),
                QgsProcessingParameterNumber.Integer if typ == 0 else QgsProcessingParameterNumber.Double,
                default, False, minimum)
            prm.setFlags(prm.flags() | prm.FlagAdvanced)
            self.addParameter(prm)
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT, self.tr("Absturzrisiko"), QgsProcessing.TypeVectorLine))

    def processAlgorithm(self, parameters, context, feedback):
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
        fields = QgsFields()
        fields.append(QgsField("seg_id", QVariant.LongLong))
        src_idx = []
        for i, f in enumerate(source.fields()):
            nf = QgsField(f)
            # Constraints (NOT NULL/UNIQUE eines Quell-Primaerschluessels) nicht
            # uebernehmen, sonst macht der PostGIS-Export dieses Feld zum PK.
            nf.setConstraints(QgsFieldConstraints())
            if f.name().lower() in key_names:
                nf.setName("src_" + f.name().lower())
            fields.append(nf)
            src_idx.append(i)
        for name, typ in OUT_FIELDS:
            if fields.indexOf(name) < 0:
                fields.append(QgsField(name, typ))
        sink, dest_id = self.parameterAsSink(parameters, self.OUTPUT, context, fields,
                                             QgsWkbTypes.LineString, dem_crs)
        if sink is None:
            raise QgsProcessingException("Ausgabe konnte nicht angelegt werden.")

        xform = QgsCoordinateTransform(source.sourceCrs(), dem_crs, context.transformContext())
        total = source.featureCount() or 1
        nseg = 0
        for k, feat in enumerate(source.getFeatures()):
            if feedback.isCanceled():
                break
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
                    try:
                        rows = process_line(chunk, dem_path, forest_path, wkt, p, thr)
                    except Exception as e:  # einzelner Weg darf den Lauf nicht abbrechen
                        feedback.reportError("Weg %s: %s" % (feat.id(), e))
                        continue
                    for xy, vals in rows:
                        of = QgsFeature(fields)
                        of.setGeometry(QgsGeometry.fromPolylineXY(
                            [QgsPointXY(float(a), float(b)) for a, b in xy]))
                        nseg += 1
                        attrs = [nseg] + list(attrs_src)
                        attrs += [None] * (fields.count() - len(attrs))
                        for name, _ in OUT_FIELDS:
                            attrs[fields.indexOf(name)] = vals[name]
                        of.setAttributes(attrs)
                        if not sink.addFeature(of, QgsFeatureSink.FastInsert):
                            raise QgsProcessingException(
                                "Schreiben in die Ausgabe fehlgeschlagen (Segment %d, Weg %s): %s"
                                % (nseg, feat.id(), sink.lastError()))
            feedback.setProgress(int(100 * (k + 1) / total))
        # Gepufferte Inserts (z.B. PostGIS) erst hier wirklich geschrieben ->
        # Fehler hier pruefen, statt sie still zu verlieren.
        if hasattr(sink, "flushBuffer") and not sink.flushBuffer():
            raise QgsProcessingException(
                "Schreiben in die Ausgabe fehlgeschlagen: %s" % sink.lastError())
        feedback.pushInfo("%d Segmente geschrieben." % nseg)
        return {self.OUTPUT: dest_id}
