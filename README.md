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
| `fall_exposure_qgis.py` | QGIS-Processing-Skript, braucht nur numpy und GDAL (in QGIS enthalten). |
| `fall_exposure.py` | Kommandozeile, mehrkernig, für grosse Gebiete. Braucht geopandas, rasterio, pyogrio, scipy. |
| `docs/fall_exposure_qgis_dokumentation.docx` | Ausführliche Dokumentation des Modells und aller Parameter. |
| `test/make_test_terrain.py` | Erzeugt synthetisches Testgelände mit bekannten Hängen, Wald und Felswand. |

## Eingabedaten

| Layer | Pflicht | Beispiel |
|---|---|---|
| Wege (Linien) | ja | OSM-Wege mit `sac_scale`, als PostGIS-Tabelle, GeoPackage o.ä. |
| Höhenmodell | ja | [swissALTI3D](https://www.swisstopo.admin.ch/de/hoehenmodell-swissalti3d) 2 m, als GeoTIFF/VRT in einem metrischen CRS |
| Walddichte | nein | [Copernicus HRL Tree Cover Density](https://land.copernicus.eu/en/products/high-resolution-layer-tree-cover-density), 0–100 % |

Viele DEM-Kacheln zu einem VRT zusammenfassen:

```bash
find . -maxdepth 1 -name '*.tiff' > tiles.txt
gdalbuildvrt -input_file_list tiles.txt merged.vrt
```

## Verwendung

### QGIS

1. Processing-Werkzeugkiste → Python-Symbol → „Skript zur Werkzeugkiste hinzufügen…“ → `fall_exposure_qgis.py`
2. **Skripte → Wanderwege → Absturzrisiko Wanderwege** öffnen
3. Wege, Höhenmodell und optional Walddichte wählen. Zum Testen „Nur gewählte Objekte“ verwenden.

Die Ausgabe kann ein temporärer Layer, ein GeoPackage oder eine PostGIS-Tabelle sein.

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

## Bekannte Grenzen

- Parameter (Reibung, Bremslänge, Faltigkeit) und Klassengrenzen (200 / 600 / 1500) sind **Startwerte ohne Kalibrierung**.
- Untergrund (Gras, Geröll, Fels) wird nicht unterschieden.
- Sicherungen wie Seile, Ketten oder Geländer sind unbekannt.
- OSM-Lagefehler von einigen Metern können den Startpunkt auf die falsche Hangseite legen.
- Felsbänder unter etwa 2–4 m Höhe glättet ein 2-m-Höhenmodell teilweise weg.
- Bisher nur mit synthetischem Gelände getestet, noch nicht systematisch mit echten Daten.

Rückmeldungen zu falsch eingestuften Wegen sind als Issue willkommen.

## Quellen

- Skitourenguru Wiki, [Artikel a0077](https://wiki.skitourenguru.com/de/articles/a0077.html)
- [DIGIWAY Exposure Map](https://www.europaregion.info/en/euregio/projekte/nachhaltigkeit/digiway-geodata-for-safety-in-the-mountains/exposure-map/) (Euregio Tirol – Südtirol – Trentino)
- OpenStreetMap Wiki, [Key:hiking_exposition](https://wiki.openstreetmap.org/wiki/Key:hiking_exposition)
- Sappington, Longshore & Thompson (2007): *Quantifying Landscape Ruggedness for Animal Habitat Analysis*, Journal of Wildlife Management 71(5)

## Lizenz

[MIT](LICENSE) © 2026 Peter Berger, BERGI IT Consulting
