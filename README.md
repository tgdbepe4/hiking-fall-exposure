# Absturzrisiko Wanderwege

Berechnet aus einem Höhenmodell (z.B. swissALTI3D) für jeden Abschnitt eines Wanderwegs, wie gefährlich ein Sturz an dieser Stelle wäre – als QGIS-Processing-Skript und als Kommandozeilen-Werkzeug.

> [!WARNING]
> **Die Berechnungen können Fehler enthalten und sind nicht validiert.**
>
> Die Ergebnisse stammen aus einem vereinfachten physikalischen Modell. Dessen Parameter und Klassengrenzen sind Schätzwerte und **nicht an echten Wegen kalibriert**. Lagefehler in OpenStreetMap, die Auflösung des Höhenmodells, fehlende Angaben zu Untergrund und Sicherungen (Seile, Ketten, Geländer) sowie Programmfehler können zu **falschen Einstufungen in beide Richtungen** führen. Ein als „gering“ eingestufter Weg kann gefährlich sein.
>
> **Verwende die Ergebnisse nicht als alleinige Grundlage für Tourenplanung oder Sicherheitsentscheidungen.** Sie ersetzen weder die Beurteilung vor Ort noch offizielle Wegklassifizierungen (gelb / weiss-rot-weiss / weiss-blau-weiss) oder die SAC-Wanderskala. Nutzung auf eigenes Risiko, siehe [Lizenz](LICENSE).

## Prinzip

Das Modell ist angelehnt an das Prinzip, das [Skitourenguru](https://wiki.skitourenguru.com/de/articles/a0077.html) beschreibt (eigene Umsetzung, die Originalparameter sind nicht öffentlich):

```
risk = max(Hangneigung, k1 × Faltigkeit, p_min) × v_max × Waldfaktor
```

| Term | Bedeutung |
|---|---|
| Hangneigung, Faltigkeit | Wie leicht stürzt man? Neigung des Geländes am Weg und dessen Unruhe (Vector Ruggedness Measure). |
| `v_max` | Wie schnell würde man? Ein virtueller Gleitblock folgt vom Wegrand der Falllinie, mit Reibung und geschwindigkeitsabhängiger Bremsung. |
| Waldfaktor | Dichter Wald (> 60 %) bremst, bis zur Halbierung bei 100 %. |

Ergebnis pro 25-m-Segment: Klasse 1–4 (gering / mässig / hoch / sehr hoch) und eine Markierung für Absturzkanten direkt neben dem Weg.

## Dateien

| Datei | Zweck |
|---|---|
| `fall_exposure_qgis.py` | QGIS-Processing-Skript (v3: Kacheln und Threads), braucht nur numpy und GDAL (in QGIS enthalten). |
| `fall_exposure.py` | Kommandozeile (älterer Stand, gleiches Modell ohne Kacheln). Braucht geopandas, rasterio, pyogrio, scipy. |
| `docs/fall_exposure_qgis_dokumentation.docx` | Ausführliche Dokumentation: Modell, Parameter, Datenvorbereitung, Laufzeiten, Darstellung. |
| `styles/*.qml` | QGIS-Stile für das Absturzrisiko und für `sac_scale` (siehe unten). |
| `test/make_test_terrain.py` | Erzeugt synthetisches Testgelände mit bekannten Hängen, Wald und Felswand. |

## Eingabedaten

| Layer | Pflicht | Beispiel |
|---|---|---|
| Wege (Linien) | ja | OSM-Wege mit `sac_scale`, als PostGIS-Tabelle, GeoPackage o.ä. |
| Höhenmodell | ja | [swissALTI3D](https://www.swisstopo.admin.ch/de/hoehenmodell-swissalti3d) 2 m, als **eine** GeoTIFF-Datei in einem metrischen CRS |
| Walddichte | nein | [Copernicus HRL Tree Cover Density](https://land.copernicus.eu/en/products/high-resolution-layer-tree-cover-density) (TCD, nicht Forest Type), 0–100 % |

Viele DEM-Kacheln zu einem VRT zusammenfassen:

```bash
find . -maxdepth 1 -name '*.tiff' > tiles.txt
gdalbuildvrt -input_file_list tiles.txt merged.vrt
```

Ein VRT aus 40 000 Kacheln ist langsam zu öffnen. Deshalb einmalig in eine Datei umwandeln (Windows: OSGeo4W Shell):

```bash
gdal_translate -of COG merged.vrt swissalti3d_2m.tif \
  -co COMPRESS=LERC_DEFLATE -co MAX_Z_ERROR=0.01 \
  -co BIGTIFF=YES -co NUM_THREADS=ALL_CPUS --config GDAL_CACHEMAX 4096
```

Beim Waldraster mit `gdalinfo -stats` prüfen, dass das Maximum 100 ist und NoData gesetzt ist.

## Verwendung

### QGIS

1. Processing-Werkzeugkiste → Python-Symbol → „Skript zur Werkzeugkiste hinzufügen…“ → `fall_exposure_qgis.py`
2. **Skripte → Wanderwege → Absturzrisiko Wanderwege (v3)** öffnen
3. Wege, Höhenmodell und optional Walddichte wählen. Zum Testen „Nur gewählte Objekte“ verwenden.

Die Ausgabe kann ein temporärer Layer, ein GeoPackage oder eine PostGIS-Tabelle sein. Für die ganze Schweiz besser eine Datei wählen, ein temporärer Layer mit Millionen Segmenten belegt mehrere GB Arbeitsspeicher.

Das Skript gruppiert die Wege in 1-km-Kacheln, liest das Höhenmodell pro Kachel einmal und rechnet die Kacheln parallel (Parameter „Threads“, Standard: alle Kerne). Jeder Lauf wird in `fall_exposure_laufzeiten.csv` protokolliert (im Projektordner, sonst im Ordner der Ausgabe).

Gemessene Laufzeiten, Schweiz, 2-m-DEM, mit Wald:

| Rechner | Wegnetz | Zeit | km/s |
|---|---|---|---|
| AMD Ryzen 7 5800X, 16 Threads | OSM `sac_scale`, 23 427 km | 36 min | 10.9 |
| Apple M4 (MacBook Air), 10 Threads | OSM `sac_scale`, 23 427 km | 19 min | 20.6 |
| Apple M4 (MacBook Air), 10 Threads | swissTLM3D, 67 020 km | 47 min | 24.0 |

Lüfterlose Geräte drosseln unter Dauerlast deutlich. Tipps dazu in der Dokumentation.

### Kommandozeile

```bash
pip install geopandas rasterio pyogrio scipy

python3 fall_exposure.py \
  --src "PG:service=meinserver" --layer sac_scale_switzerland \
  --dem merged.vrt --forest tcd_2018_010m.tif \
  --out exposure.gpkg --out-layer exposure \
  --summary-layer exposure_ways --workers 8
```

`python3 fall_exposure.py --help` listet alle Parameter.

## Ausgabefelder

| Feld | Bedeutung |
|---|---|
| `exp_class`, `exp_label` | Klasse 1–4 und Bezeichnung |
| `risk` | Risikowert (Grad · m/s) |
| `slope_deg`, `fold` | Hangneigung und Faltigkeit am Weg |
| `v_max`, `drop_m` | Höchstgeschwindigkeit und Fallhöhe im Fallexperiment |
| `forest_pct`, `forest_mod` | Walddichte und Waldfaktor |
| `immediate_drop` | 1 = Absturzkante direkt neben dem Weg |
| `seg_id`, `src_id` | Segmentnummer und `id` des Original-Wegs (nur QGIS-Version) |

## Stile

| Datei | Zweck |
|---|---|
| `styles/absturzrisiko_band.qml` | Absturzrisiko als breites Band unter der Weglinie (2 hellorange, 3 orange, 4 rot). |
| `styles/absturzrisiko_skitourenguru_farben.qml` | Absturzrisiko als Linie grün/gelb/orange/rot. |
| `styles/absturzrisiko_skitourenguru_farben_mit_kante.qml` | Wie oben, plus Markierung `immediate_drop`. |
| `styles/sac_scale_signalisation_lizmap.qml` | `sac_scale` nach Schweizer Signalisation mit weissem Rand und T-Beschriftung. |

Für eine kombinierte Karte liegt der Risiko-Layer mit dem Band **unter** dem `sac_scale`-Layer. Beide müssen aus demselben Wegnetz stammen.

## Bekannte Grenzen

- Parameter (Reibung, Bremslänge, Faltigkeit) und Klassengrenzen (200 / 600 / 1500) sind **Startwerte ohne Kalibrierung**.
- Untergrund (Gras, Geröll, Fels) wird nicht unterschieden.
- Sicherungen wie Seile, Ketten oder Geländer sind unbekannt.
- OSM-Lagefehler von einigen Metern können den Startpunkt auf die falsche Hangseite legen.
- Felsbänder unter etwa 2–4 m Höhe glättet ein 2-m-Höhenmodell teilweise weg.
- Kletterstellen (technische Schwierigkeit) und Gletscher (Spalten) erfasst das Modell nicht. Dafür ist `sac_scale` da.
- Mit synthetischem Gelände und mit dem ganzen Schweizer Wegnetz getestet und an bekannten Gebieten mit Skitourenguru verglichen (grosse Linie stimmt, stellenweise eine Klasse höher). Eine systematische Validierung fehlt.

Rückmeldungen zu falsch eingestuften Wegen sind als Issue willkommen.

## Quellen

- Skitourenguru Wiki, [Artikel a0077](https://wiki.skitourenguru.com/de/articles/a0077.html)
- [DIGIWAY Exposure Map](https://www.europaregion.info/en/euregio/projekte/nachhaltigkeit/digiway-geodata-for-safety-in-the-mountains/exposure-map/) (Euregio Tirol – Südtirol – Trentino)
- OpenStreetMap Wiki, [Key:hiking_exposition](https://wiki.openstreetmap.org/wiki/Key:hiking_exposition)
- [Copernicus HRL Tree Cover Density](https://land.copernicus.eu/en/products/high-resolution-layer-tree-cover-density)
- Sappington, Longshore & Thompson (2007): *Quantifying Landscape Ruggedness for Animal Habitat Analysis*, Journal of Wildlife Management 71(5)

## Lizenz

[MIT](LICENSE) © 2026 Peter Berger, BERGI IT Consulting
