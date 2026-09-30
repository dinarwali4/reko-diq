# Reko Diq copper–gold alteration map

Satellite mapping of the porphyry alteration halo at Reko Diq, Chagai District, Balochistan, Pakistan (28.88–29.06°N, 62.30–62.55°E).

**Live platform:** https://dinarwali4.github.io/reko-diq/

RegenX.eco Geospatial Computational Model · Prepared by **Dinar Wali** · dinar.wali@regenx.eco · [LinkedIn](https://www.linkedin.com/in/dinar-wali-621954102/)

## Findings
Persistent alteration was mapped from NASA EMIT hyperspectral data. A cell counts when its band depth is more than 2.5 robust standard deviations above the scene on both Oct 2023 and Oct 2024.

| Facies | Area | Largest zone (centroid, ±60 m) |
|---|---|---|
| Phyllic (sericite, 2200 nm) | 1,068 ha | RZ1: 386 ha, 29.03682 N 62.32769 E |
| Leached cap (jarosite 2265 nm + Fe³⁺) | 319 ha | RZ3: 134 ha, 28.99687 N 62.38839 E |
| Argillic (kaolinite/alunite 2165 nm) | 77 ha | RZ8: 28 ha, 28.99236 N 62.35632 E |

**Validation:**
- No alluvial control cells qualify (0 of 16,000).
- The independent Feb 2025 pass reproduces 79–87% of each facies.
- ASTER data from 2000–2003 confirms the argillic zones (z +5.4).

**Construction, 2025–26 (Sentinel-2):**
- Airstrip, about 1.9 km, at 28.918 N 62.542 E
- Camp expansion, new ponds and pads, and a new corridor

**Limits.** This maps surface alteration, not copper grade. The published reserves come from drilling: 0.48% Cu and 0.28 g/t Au (Barrick, 2025).

## Contents
| Path | What |
|---|---|
| `docs/` | Interactive platform (GitHub Pages), annotated Google Earth imagery 2003–2026, KML |
| `RekoDiq_RegenX_alteration.kml` | Facies outlines, zone labels and construction features for Google Earth |
| `data/` | Zones, 2025 validation, Sentinel-2 change patches, EMIT mineral persistence (JSON) |
| `figures/` | Facies map, spectra, change chips, quick-look maps |
| `scripts/` | Python pipeline |

## Pipeline
1. **Download:** `s2_fetch.py`, `landsat_fetch.py`, then `dl.sh`. EMIT and ASTER need a NASA Earthdata token in `~/.edl_token`.
2. **Process EMIT:** `emit_min.py`, `emit_rfl.py`, `emit_feat.py`.
3. **Map and validate:** `facies.py`, `zones.py`, `aster.py`, `validate2025.py`.
4. **Measure change:** `s2_change.py`, `landsat_change.py`.
5. **Outputs:** `figs.py`, `chips.py`, `make_kml.py`, `finish_ge.py`.

## Sources
- Barrick 2024 reserves release (Feb 2025) and the 2025 updated feasibility study.
- *Economic Geology* 109(7): Western Porphyry Cu–Au systems at Reko Diq.
- Data: NASA EMIT and ASTER (LP DAAC), Copernicus Sentinel-2 and GLO-30, USGS Landsat.
- Google Earth Pro imagery © Airbus and Maxar Technologies, used for annotation.
