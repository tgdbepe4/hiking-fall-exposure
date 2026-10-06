#!/usr/bin/env python3
"""
fall_exposure.py v2 - Absturzrisiko entlang von Wanderwegen aus einem DEM ableiten.

Modell (angelehnt an Skitourenguru, Wiki-Artikel a0077, eigene Umsetzung):

    risk = max(slope, k1 * fold, p_min) * v_max * forest_mod

  slope       Hangneigung [Grad] am Weg (Maximum im Umkreis der Wegtrasse)
  fold        "Faltigkeit" des Gelaendes am Weg: Vector Ruggedness Measure
              (VRM, 0 = glatt ... 1 = extrem zerfurcht), mit k1 auf Grad skaliert
  p_min       Untergrenze der Sturzwahrscheinlichkeit (stolpern kann man ueberall)
  v_max       Maximale Geschwindigkeit [m/s] in einem virtuellen Fallexperiment:
              Ein Gleitblock startet beidseitig am Wegrand und folgt der
              Falllinie (steilster Abstieg). Energielinie mit Reibung mu:
                  h_kin += dz - mu * dl_horizontal ;  v = sqrt(2 g h_kin)
              Abbruch, sobald h_kin auf 0 faellt (Sturz gestoppt) oder max_path
              erreicht ist.
  forest_mod  Walddichte fd [%] entlang der Sturzbahn:
              fd < 60 -> 1,  sonst 1 - (fd - 60) * 0.5 / 40   (100 % -> 0.5)
              Optional, z.B. Copernicus HRL Tree Cover Density (10 m, 0-100 %).

  Zusatzregel: Hangneigung am Weg > 35 Grad -> mindestens Klasse 2.

  Klassen (risk in Grad*m/s, Startwerte - KALIBRIEREN!):
    1 gering < 200 <= 2 maessig < 600 <= 3 hoch < 1500 <= 4 sehr hoch

  immediate_drop = 1, wenn direkt neben dem Weg (near_dist) Gelaende steiler als
  cliff_deg liegt (analog "pinke" Signatur DIGIWAY).

Eingabe:
  - Wege: beliebige OGR-Quelle, z.B. "PG:service=bergiserver" --layer sac_scale_switzerland
  - DEM:  swissALTI3D (EPSG:2056) als GeoTIFF/VRT (gdalbuildvrt swissalti3d.vrt *.tif);
          wird intern auf --res gemittelt (Standard 2 m; Skitourenguru nutzt 5 m)
  - Wald: optional --forest, Raster mit Walddichte 0-100 % in beliebigem CRS

Ausgabe:
  Segmente (Standard 25 m) mit risk, Klasse und Kennzahlen nach GPKG oder PostGIS,
  optional zusaetzlich eine Zeile pro Weg (--summary-layer).

Beispiel:
  python3 fall_exposure.py --src "PG:service=bergiserver" --layer sac_scale_switzerland \\
      --dem swissalti3d.vrt --forest tcd_2018_010m_ch.tif \\
      --out "PG:service=bergiserver" --out-layer exposure_switzerland \\
      --summary-layer exposure_switzerland_ways --workers 8

Das Resultat ist ein Modell und ersetzt keine Beurteilung im Gelaende.
OSM-Lagefehler von einigen Metern koennen die Werte verfaelschen.
"""

import argparse
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import geopandas as gpd
import pyogrio
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import from_bounds
from scipy.ndimage import map_coordinates, uniform_filter
from shapely.geometry import LineString, MultiLineString
from shapely.ops import substring

G = 9.81

DEFAULTS = dict(
    seg_len=25.0,         # Laenge der Ausgabesegmente [m]
    step=5.0,             # Abstand der Startpunkte entlang des Weges [m]
    res=2.0,              # Arbeitsaufloesung DEM [m]
    trail_halfwidth=1.5,  # Wegtrasse beidseitig [m]; Start der Sturzbahn am Wegrand
    mu=0.6,               # Gleitreibung der Sturzbahn (Energielinie ~31 Grad)
    drag_len=100.0,       # Bremslaenge fuer v^2-Term [m] (0 = aus)
    max_path=400.0,       # maximale horizontale Laenge der Sturzbahn [m]
    k1=300.0,             # Skalierung VRM -> Grad
    fold_win=5.0,         # Fenster fuer VRM [Zellen]
    p_min=10.0,           # Untergrenze Wahrscheinlichkeitsterm [Grad]
    forest_len=50.0,      # Walddichte gemittelt ueber die ersten x m der Sturzbahn
    steep_deg=35.0,       # Zusatzregel: Hangneigung > steep_deg -> mind. Klasse 2
    cliff_deg=60.0,       # Absturzkante
    near_dist=5.0,        # "direkt neben dem Weg" [m]
    chunk_len=500.0,      # lange Wege werden in Stuecke dieser Laenge geteilt [m]
)
DEFAULT_THR = "200,600,1500"
CLASS_LABEL = {1: "gering", 2: "maessig", 3: "hoch", 4: "sehr hoch"}


# --------------------------------------------------------------------------
# Raster lesen
# --------------------------------------------------------------------------
def read_window(path, bounds, res, crs=None, resampling=Resampling.average):
    """Liest bounds auf einem Raster mit Zellgroesse res. crs != None -> on-the-fly
    Reprojektion (fuer Waldraster in anderem CRS). Rueckgabe (array, transform)."""
    minx, miny, maxx, maxy = bounds
    w = max(2, int(math.ceil((maxx - minx) / res)))
    h = max(2, int(math.ceil((maxy - miny) / res)))
    tr = rasterio.Affine(res, 0, minx, 0, -res, miny + h * res)
    with rasterio.open(path) as src:
        nod = src.nodata
        if crs is not None and src.crs != crs:
            # Ziel-Grid direkt als VRT definieren; Bereiche ausserhalb -> nodata
            with WarpedVRT(src, crs=crs, transform=tr, width=w, height=h,
                           resampling=resampling,
                           nodata=nod if nod is not None else -9999) as vrt:
                arr = vrt.read(1)
                if nod is None:
                    nod = -9999
        else:
            win = from_bounds(minx, miny, minx + w * res, miny + h * res, src.transform)
            arr = src.read(1, window=win, out_shape=(h, w), resampling=resampling,
                           boundless=True, fill_value=nod if nod is not None else np.nan)
    arr = arr.astype("float32")
    if nod is not None:
        arr[arr == nod] = np.nan
    return arr, tr


def vrm(z, res, win):
    """Vector Ruggedness Measure (Sappington et al. 2007)."""
    gy, gx = np.gradient(z, res)
    gy = -gy  # Zeilen laufen nach Sueden
    norm = np.sqrt(gx * gx + gy * gy + 1.0)
    nx, ny, nz = -gx / norm, -gy / norm, 1.0 / norm
    k = max(3, int(win) | 1)
    sx, sy, sz = (uniform_filter(a, k, mode="nearest") for a in (nx, ny, nz))
    return 1.0 - np.sqrt(sx * sx + sy * sy + sz * sz)


# --------------------------------------------------------------------------
# Kern: ein (Teil-)LineString in DEM-CRS
# --------------------------------------------------------------------------
def process_line(line, dem_path, forest_path, dem_crs, p, thr):
    L = line.length
    if L < 1.0:
        return []
    res = p["res"]

    dists = np.arange(p["step"] / 2.0, L, p["step"])
    if dists.size == 0:
        dists = np.array([L / 2.0])
    pts = np.array([line.interpolate(d).coords[0] for d in dists])
    eps = min(2.0, L / 4.0)
    a = np.array([line.interpolate(max(d - eps, 0)).coords[0] for d in dists])
    b = np.array([line.interpolate(min(d + eps, L)).coords[0] for d in dists])
    t = b - a
    t /= np.linalg.norm(t, axis=1, keepdims=True).clip(min=1e-9)
    nrm = np.stack([-t[:, 1], t[:, 0]], axis=1)

    reach = p["max_path"] + 4 * res + 10
    bounds = (*(pts.min(axis=0) - reach), *(pts.max(axis=0) + reach))
    z, tr = read_window(dem_path, bounds, res)
    if np.isnan(z).all():
        return []
    # Loecher (Rand/NoData) mit Mittelwert fuellen, damit Gradienten stabil bleiben
    z = np.where(np.isnan(z), np.nanmean(z), z)
    gy, gx = np.gradient(z, res)
    gy = -gy                                   # +y = Norden
    slope_r = np.degrees(np.arctan(np.hypot(gx, gy)))
    fold_r = vrm(z, res, p["fold_win"])
    fd_r = None
    if forest_path:
        fd_r, _ = read_window(forest_path, bounds, res, crs=dem_crs,
                              resampling=Resampling.bilinear)
        fd_r = np.nan_to_num(fd_r, nan=0.0).clip(0, 100)

    inv = ~tr

    def samp(arr, x, y, order=1):
        c, r = inv * (x, y)
        return map_coordinates(arr, [np.asarray(r) - 0.5, np.asarray(c) - 0.5],
                               order=order, mode="nearest")

    # --- Wahrscheinlichkeitsterm am Weg: Maximum im Umkreis der Trasse ---
    offs = np.linspace(-p["trail_halfwidth"] - res / 2, p["trail_halfwidth"] + res / 2, 5)
    qx = pts[:, 0, None] + offs[None, :] * nrm[:, 0, None]
    qy = pts[:, 1, None] + offs[None, :] * nrm[:, 1, None]
    slope_w = samp(slope_r, qx, qy).max(axis=1)
    fold_w = samp(fold_r, qx, qy).max(axis=1)
    z0 = np.median(samp(z, qx, qy), axis=1)

    # immediate_drop: Neigung > cliff_deg innerhalb near_dist beidseitig
    nd = np.arange(p["trail_halfwidth"], p["near_dist"] + 1e-6, max(res / 2, 0.5))
    imm = np.zeros(len(pts), bool)
    for sgn in (1, -1):
        ix = pts[:, 0, None] + sgn * nd[None, :] * nrm[:, 0, None]
        iy = pts[:, 1, None] + sgn * nd[None, :] * nrm[:, 1, None]
        zz = samp(z, ix, iy)
        imm |= (samp(slope_r, ix, iy) >= p["cliff_deg"]).any(axis=1) & \
               ((z0[:, None] - zz).max(axis=1) > 2.0)

    # --- Virtuelles Fallexperiment: beidseitig vom Wegrand der Falllinie folgen ---
    sx = np.concatenate([pts[:, 0] + p["trail_halfwidth"] * nrm[:, 0],
                         pts[:, 0] - p["trail_halfwidth"] * nrm[:, 0]])
    sy = np.concatenate([pts[:, 1] + p["trail_halfwidth"] * nrm[:, 1],
                         pts[:, 1] - p["trail_halfwidth"] * nrm[:, 1]])
    zprev = np.concatenate([z0, z0])
    hk = np.zeros_like(sx)          # kinetische Hoehe [m]
    vmax = np.zeros_like(sx)
    drop_at_vmax = np.zeros_like(sx)
    fd_sum = np.zeros_like(sx)
    fd_n = np.zeros_like(sx)
    alive = np.ones_like(sx, bool)
    x, y = sx.copy(), sy.copy()
    nsteps = int(p["max_path"] / res)
    z_start = zprev.copy()
    for i in range(nsteps + 1):
        if i > 0:
            ggx = samp(gx, x, y)
            ggy = samp(gy, x, y)
            gm = np.hypot(ggx, ggy)
            alive &= gm > 1e-3
            d = np.where(alive, res / np.maximum(gm, 1e-9), 0.0)
            x = x - ggx * d
            y = y - ggy * d
            dl = np.where(alive, res, 0.0)
        else:
            dl = np.full_like(x, p["trail_halfwidth"])
        zc = samp(z, x, y)
        dz = zprev - zc
        drag = hk * dl / p["drag_len"] if p["drag_len"] > 0 else 0.0
        hk = np.where(alive, hk + dz - p["mu"] * dl - drag, hk)
        stopped = alive & (hk <= 0) & (i > 0)
        hk = np.maximum(hk, 0.0)
        alive &= ~stopped
        v = np.sqrt(2 * G * hk)
        better = v > vmax
        vmax = np.where(better, v, vmax)
        drop_at_vmax = np.where(better, z_start - zc, drop_at_vmax)
        if fd_r is not None and i * res <= p["forest_len"]:
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

    prob = np.maximum.reduce([slope_w, p["k1"] * fold_w, np.full(n, p["p_min"])])
    risk_sides = np.concatenate([prob, prob]) * vmax * fmod
    side = (risk_sides[n:] > risk_sides[:n]).astype(int)
    pick = lambda arr: np.where(side == 1, arr[n:], arr[:n])
    risk, v_pt, drop_pt, fmod_pt, fd_pt = (pick(a) for a in
                                           (risk_sides, vmax, drop_at_vmax, fmod, fd))

    cls = np.digitize(risk, thr).astype(np.int8) + 1
    cls = np.where(slope_w > p["steep_deg"], np.maximum(cls, 2), cls)

    out = []
    nseg = max(1, int(math.ceil(L / p["seg_len"])))
    edges = np.linspace(0, L, nseg + 1)
    for i in range(nseg):
        m = (dists >= edges[i]) & (dists < edges[i + 1])
        if not m.any():
            continue
        # robust: zweithoechster Wert, damit ein einzelner DEM-Ausreisser nicht dominiert
        o = np.argsort(risk[m])
        j = np.flatnonzero(m)[o[-2] if m.sum() >= 4 else o[-1]]
        c = np.sort(cls[m])
        seg_cls = int(c[-2]) if c.size >= 4 else int(c[-1])
        out.append(dict(
            geometry=substring(line, edges[i], edges[i + 1]),
            exp_class=seg_cls,
            exp_label=CLASS_LABEL[seg_cls],
            risk=round(float(risk[j]), 0),
            slope_deg=round(float(slope_w[m].max()), 1),
            fold=round(float(fold_w[m].max()), 3),
            v_max=round(float(v_pt[j]), 1),
            drop_m=round(float(drop_pt[j]), 1),
            forest_pct=round(float(fd_pt[j]), 0),
            forest_mod=round(float(fmod_pt[j]), 2),
            immediate_drop=int(imm[m].sum() >= (2 if m.sum() >= 4 else 1)),
        ))
    return out


def split_long(line, chunk_len):
    L = line.length
    if L <= chunk_len * 1.5:
        return [line]
    n = int(math.ceil(L / chunk_len))
    e = np.linspace(0, L, n + 1)
    return [substring(line, e[i], e[i + 1]) for i in range(n)]


def worker(task):
    attrs, geom_wkb, dem, forest, dem_crs_wkt, p, thr = task
    from shapely import wkb
    from rasterio.crs import CRS
    dem_crs = CRS.from_wkt(dem_crs_wkt)
    geom = wkb.loads(geom_wkb)
    parts = list(geom.geoms) if isinstance(geom, MultiLineString) else [geom]
    rows = []
    for part in parts:
        if not isinstance(part, LineString) or part.is_empty:
            continue
        for chunk in split_long(part, p["chunk_len"]):
            for r in process_line(chunk, dem, forest, dem_crs, p, thr):
                r.update(attrs)
                rows.append(r)
    return rows


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Absturzrisiko entlang von Wanderwegen aus einem DEM.")
    ap.add_argument("--src", required=True, help="OGR-Quelle der Wege (z.B. PG:service=bergiserver)")
    ap.add_argument("--layer", help="Layer/Tabelle (z.B. sac_scale_switzerland)")
    ap.add_argument("--where", help="SQL-Filter, z.B. \"sac_scale <> 'hiking'\"")
    ap.add_argument("--bbox", help="minx,miny,maxx,maxy im CRS der Quelle (Testausschnitt)")
    ap.add_argument("--dem", required=True, help="DEM (GeoTIFF/VRT) in metrischem CRS")
    ap.add_argument("--forest", help="optional: Walddichte-Raster 0-100 %%")
    ap.add_argument("--out", required=True, help="Datei (.gpkg) oder OGR-Ziel (PG:...)")
    ap.add_argument("--out-layer", default="fall_exposure")
    ap.add_argument("--summary-layer", help="optional: eine Zeile pro Weg")
    ap.add_argument("--keep", default="osm_id,name,sac_scale",
                    help="Attribute aus der Quelle, die mitkopiert werden")
    ap.add_argument("--thr", default=DEFAULT_THR,
                    help="Klassengrenzen risk (3 Werte, Grad*m/s), Standard " + DEFAULT_THR)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    for k, v in DEFAULTS.items():
        ap.add_argument("--" + k.replace("_", "-"), type=float, default=v)
    args = ap.parse_args()
    p = {k: getattr(args, k) for k in DEFAULTS}
    thr = np.array([float(v) for v in args.thr.split(",")])
    if thr.size != 3 or np.any(np.diff(thr) <= 0):
        sys.exit("--thr braucht 3 aufsteigende Werte.")

    t0 = time.time()
    with rasterio.open(args.dem) as src:
        dem_crs = src.crs
    if dem_crs is None or not dem_crs.is_projected:
        sys.exit("DEM braucht ein metrisches, projiziertes CRS (z.B. EPSG:2056).")

    kw = dict(layer=args.layer, where=args.where)
    if args.bbox:
        kw["bbox"] = tuple(float(v) for v in args.bbox.split(","))
    gdf = pyogrio.read_dataframe(args.src, **{k: v for k, v in kw.items() if v is not None})
    if gdf.empty:
        sys.exit("Keine Wege gefunden.")
    if gdf.crs is None:
        sys.exit("Quelle hat kein CRS.")
    gdf = gdf.to_crs(dem_crs).reset_index(drop=True)
    gdf["_way"] = gdf.index
    keep = [c for c in args.keep.split(",") if c and c in gdf.columns]
    print(f"{len(gdf)} Wege geladen ({time.time()-t0:.1f}s), {args.workers} Worker", flush=True)

    tasks = []
    for _, row in gdf.iterrows():
        if row.geometry is None or row.geometry.is_empty:
            continue
        attrs = {c: row[c] for c in keep}
        attrs["_way"] = int(row["_way"])
        tasks.append((attrs, row.geometry.wkb, args.dem, args.forest, dem_crs.to_wkt(), p, thr))

    rows, done = [], 0
    ex = ProcessPoolExecutor(max_workers=args.workers) if args.workers > 1 else None
    it = ex.map(worker, tasks, chunksize=16) if ex else map(worker, tasks)
    for r in it:
        rows.extend(r)
        done += 1
        if done % 1000 == 0:
            print(f"  {done}/{len(tasks)} Wege, {time.time()-t0:.0f}s", flush=True)
    if ex:
        ex.shutdown()
    if not rows:
        sys.exit("Keine Segmente berechnet (DEM deckt die Wege nicht ab?).")

    seg = gpd.GeoDataFrame(rows, geometry="geometry", crs=dem_crs)
    wkw = dict(layer=args.out_layer)
    if args.out.upper().startswith("PG:"):
        wkw.update(driver="PostgreSQL",
                   layer_options={"OVERWRITE": "YES", "GEOMETRY_NAME": "geom"})
    pyogrio.write_dataframe(seg.drop(columns="_way"), args.out, **wkw)

    if args.summary_layer:
        seg["_len"] = seg.geometry.length
        g = seg.groupby("_way")
        summ = gdf.loc[g.size().index, keep + ["geometry"]].copy()
        summ["exp_max"] = g["exp_class"].max().values
        summ["risk_max"] = g["risk"].max().values
        summ["len_m"] = g["_len"].sum().round(0).values
        for c in (3, 4):
            summ[f"len_ge{c}_m"] = (seg[seg.exp_class >= c].groupby("_way")["_len"].sum()
                                    .reindex(summ.index).fillna(0).round(0).values)
        summ["immediate_drop"] = g["immediate_drop"].max().values
        wkw["layer"] = args.summary_layer
        pyogrio.write_dataframe(gpd.GeoDataFrame(summ, geometry="geometry", crs=dem_crs),
                                args.out, **wkw)

    dist = seg.groupby("exp_class").geometry.apply(lambda s: s.length.sum() / 1000)
    print("Fertig in %.0fs. Laenge je Klasse [km]:" % (time.time() - t0))
    for c, km in dist.items():
        print(f"  {c} {CLASS_LABEL[c]:<10} {km:10.1f}")


if __name__ == "__main__":
    main()
