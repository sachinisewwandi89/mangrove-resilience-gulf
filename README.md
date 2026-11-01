# Mangrove Resilience and Recovery in the Gulf of Carpentaria

Satellite time-series analysis of mangrove dieback and recovery following the
2015–16 mass mortality event in the Gulf of Carpentaria, northern Australia.

**ENST90007 Research Project · Semester 2, 2026 · University of Melbourne**
Supervisor: A/Prof Jagannath Aryal, Department of Infrastructure Engineering

---

## Research question

How have mangrove ecosystems in northern Australia responded to and recovered
from climate extremes since the late 1980s, and can recovery capacity be
predicted spatially from satellite and environmental data?

## Study area

Southern and western Gulf of Carpentaria coastline, Northern Territory and
Queensland: **135.0–141.5°E, 14.0–18.0°S**.

| Quantity | Value |
|---|---|
| Mangrove area (GMW v3 union) | 79,140.6 ha |
| Mangrove pixels at 30 m | 920,278 |
| Stable extent (GMW core) | 53,169.1 ha (67.2%) |
| Changed extent status | 25,971.6 ha (32.8%) |
| Time series | 1987–2025, 39 annual composites |

## Headline results

| Finding | Value |
|---|---|
| Mapped 2015–16 dieback | 2,522 ha (29,453 pixels, 952 patches) |
| Design-based dieback area | **7,335 ha** (95% CI 4,626–10,044) |
| Producer's accuracy, dieback class | 26.9 ± 10.1% |
| User's accuracy, dieback class | 78.0 ± 6.7% |
| Never regained 80% of pre-disturbance condition | 26.7% (18.2–38.6%) |
| Recovered then declined again | 9.6% (7.4–11.9%) |
| Median time to recovery | 8 years |
| Blocked cross-validated R² | 0.365 (against 0.759 under random folds) |
| Mangrove at risk, western Gulf | 46.6% (13,511 ha) |
| Mangrove at risk, southeastern Gulf | 7.3% (1,728 ha) |

---

## Repository contents

### `Scripts/` — Google Earth Engine (JavaScript)

Run in the Earth Engine Code Editor, in this order.

| File | Purpose |
|---|---|
| `00_setup_test.js` | Environment check |
| `01_annual_composites.js` | Original Phase 1 pipeline (superseded — retained for the harmonisation before/after comparison) |
| `01b_annual_composites_harmonised.js` | **Phase 1 production script.** Roy et al. harmonisation, GMW union mask, diagnostics D1–D8 |
| `01c_tidal_screening.js` | MNDWI tidal screening test. Not retained — see Appendix A of the thesis |
| `02_landtrendr.js` | Phase 2, test strip. Segmentation, three-class map, stratum weights |
| `03_change_detection.js` | Investigation of the 2011 detections and the pre-1996 mask survivorship test |
| `04_full_coast.js` | **Phase 2 production run.** Full-coast segmentation; exports the three class and history assets |
| `05_resilience_metrics.js` | Phase 3 resilience measures from fitted segments |
| `06_sample.js` | **Validation.** Draws the 400-point stratified reference sample and exports the blind interpretation sheet |

> **Earth Engine assets.** Scripts 04 and 06 write to and read from
> `projects/mangroveresilience2026/assets/mangrove-resilience/`. These assets
> are not public. To reproduce from scratch, run `04_full_coast.js` first and
> point `ASSET_DIR` in `06_sample.js` at your own asset folder.

### `notebooks/` — Python analysis

Run in order; each writes the inputs for the next.

| Notebook | Produces |
|---|---|
| `01_phase3_completion.ipynb` | The 29,453-pixel dieback census, resilience measures, Kaplan–Meier |
| `02_accuracy_assessment.ipynb` | Error matrix, design-based areas, Pontius decomposition, Bradley–Terry |
| `03_drivers.ipynb` | Covariate join, random forest, spatially blocked cross-validation, capacity projection |
| `04_maps_and_figures.ipynb` | Maps and figure panels |
| `05_statistical_robustness.ipynb` | Variogram, effective sample size, ALE, area of applicability, detection-limit tests |

### `src/` — Python scripts

| Group | Files |
|---|---|
| Figures | `fig0_study_area.py`, `fig1_trajectories.py`, `fig3_v2_detection_floor.py`, `Fig4disagreement.py`, `Fig5_Recovery.py`, `Fig6_Validation.py`, `Fig7_Capacity.py`, `fig_accuracy_plain.py`, `fig_missed_and_adjusted.py`, `fig_before_after.py`, `maps_panels.py` |
| Data retrieval | `fetch_silo_rainfall.py`, `fetch_silo_temperature.py`, `fetch_dea_tidal.py` |
| Presentation | `build_rs_deck.py`, `restyle_deck.py`, `make_motif.py`, `make_rs_backgrounds.py` |

Each figure script prints the values it produces so they can be checked
against the manuscript.

### `tables/` — Results

Accuracy summary, model performance, robustness statistics and the numbered
tables reported in the thesis, as JSON and CSV.

### `figures/` — Output figures

The seven figures reproduced in the thesis, at 300–400 dpi.

### `docs/` — Working documentation

Decisions log, draft methodology, supervisor briefings, statistical notes.

---

## Reproducing the analysis

```bash
git clone https://github.com/sachinisewwandi89/mangrove-resilience-gulf.git
cd mangrove-resilience-gulf
pip install -r requirements.txt
```

The four small data files needed by most figure scripts are included under
`data/`. The remaining intermediate files (about 0.5 GB, one of them over
GitHub's 100 MB limit) are excluded and regenerated by the notebooks — see
`data/README.md` for what each one is and which notebook writes it.

Retrieving SILO climate data requires an email address as the API username:

```bash
export SILO_EMAIL="your.address@example.com"     # Windows: set SILO_EMAIL=...
python src/fetch_silo_rainfall.py
```

---

## Method summary

**Phase 1 — Time-series construction.** Landsat Collection 2 Tier 1 Level-2
surface reflectance, 1987–2025, from Landsat 5 TM, 7 ETM+, 8 OLI and 9 OLI-2.
Cloud, shadow, cirrus and fill masked via QA_PIXEL. Landsat 8/9 harmonised to
ETM+ equivalent space using Roy et al. (2016) coefficients. NBR computed per
image. Dry-season (May–October) median composites, masked to GMW v3 union
mangrove extent.

**Phase 2 — Disturbance detection.** LandTrendr temporal segmentation
(Kennedy et al. 2010, 2018) called directly on the harmonised composites. Two
products: an unrestricted greatest-loss disturbance history, and a
year-windowed 2015–16 dieback map.

**Phase 3 — Resilience metrics.** Magnitude of loss, recovery rate, recovery
ratio and years to recovery, derived per pixel from fitted segments, with
climate adjustment by difference-in-differences against contemporaneous
undisturbed stands.

**Phase 4 — Validation.** 400-point stratified reference sample interpreted
blind. Design-based area estimation (Olofsson et al. 2014), disagreement
decomposed into quantity, allocation, exchange and shift (Pontius & Millones
2011; Pontius & Santacruz 2014), and error direction estimated with a
standardised Bradley–Terry model (Stein, Aryal & Gort 2005).

**Phase 5 — Driver modelling.** Random forest against climatic, hydrological
and geomorphic covariates under spatially blocked cross-validation, with
accumulated local effects and projection restricted to the area of
applicability.

---

## Key parameters

**Compositing**

| Parameter | Value |
|---|---|
| Compositing window | 1 May – 31 October (dry season) |
| Reducer | Median |
| Scaling | × 0.0000275, offset −0.2, applied before index |
| Roy et al. NIR | slope 0.8462, intercept 0.0412 |
| Roy et al. SWIR2 | slope 0.9071, intercept 0.0172 |

**LandTrendr segmentation**

| Parameter | Value |
|---|---|
| maxSegments | 6 |
| spikeThreshold | 0.9 |
| vertexCountOvershoot | 3 |
| preventOneYearRecovery | true |
| recoveryThreshold | 0.25 |
| pvalThreshold | 0.05 |
| bestModelProportion | 0.75 |
| minObservationsNeeded | 6 |

**Change extraction**

| Parameter | Value |
|---|---|
| delta | loss |
| sort | greatest |
| magnitude threshold | 0.15 NBR (native units, not scaled integers) |
| minimum mapping unit | 5 pixels (0.45 ha) |

**Modelling**

| Parameter | Value |
|---|---|
| Trees | 300 |
| Minimum leaf size | 5 |
| Candidate covariates per split | √p |
| Cross-validation | 5 repeats × 5 folds, 0.25° spatial blocks (~27 km) |
| Block bootstrap | 2,000 resamples at 0.25° |

---

## Reproducibility notes

Four points that are easy to get wrong and that materially affect results.

**Sensor harmonisation is required, and the LT-GEE module does not provide
it.** Module version 0.2.0 states on loading that it does not apply the Roy et
al. coefficients, so its `runLT()` wrapper rebuilds composites without
harmonisation. This project calls
`ee.Algorithms.TemporalSegmentation.LandTrendr` directly on harmonised
composites instead. Verification: mean NBR over stable mangrove shows no step
change at either the 1999 or 2013 sensor transition. An unharmonised series
overstated post-2013 decline by ~0.19 NBR units and misidentified the year of
minimum condition by one year.

**NBR must be inverted before segmentation.** The algorithm assumes increasing
values represent disturbance; NBR falls when canopy is lost. Calling the
algorithm directly makes the inversion the caller's responsibility.

**Magnitude thresholds are in native index units here.** LT-GEE defaults
assume scaled integers; a threshold of 100 against native 0–1 reflectance is
unsatisfiable and returns zero detections.

**A map pixel count is not an area estimate.** Counting the pixels the map
marks gives 2,522 ha; design-based estimation against the reference sample
gives 7,335 ha. The difference is omission of sparse and marginal canopy below
the sub-pixel detection limit, not noise.

---

## Statistical methods and reference implementations

The statistical methods were implemented in Python following their published
descriptions, guided by the reference implementations provided by the original
authors. No code from these packages was reused directly; each method was
reimplemented so that it could be applied to the pixel-level data structures
used here.

| Method | Reference implementation |
|---|---|
| Quantity / allocation / exchange / shift | diffeR (Pontius & Santacruz, 2023), R |
| Accumulated local effects | Apley & Zhu (2020); ALEPlot, R, archived from CRAN 2025 |
| Area of applicability | CAST (Meyer et al., 2025), R |
| Bradley–Terry paired comparison | BradleyTerry2 (Turner & Firth, 2012), R |
| Kaplan–Meier estimation | lifelines (Davidson-Pilon, 2019), Python |

Analysis additionally used NumPy, pandas, SciPy, scikit-learn, rasterio,
Matplotlib and Cartopy.

---

## Data sources

| Dataset | Source |
|---|---|
| Landsat Collection 2 L2 SR | USGS, via Google Earth Engine |
| Global Mangrove Watch v3 | Bunting et al. (2022), GEE community catalogue |
| SILO gridded climate | Queensland Government, Long Paddock |
| Water observations | Digital Earth Australia (Mueller et al. 2016) |
| Queensland dieback baseline | Accad et al. (2019), Queensland Herbarium |

All datasets are openly licensed. No fieldwork, no human participants, no live
organisms.

---

## References

Accad, A., Li, J., Dowling, R., Neldner, V.J. and Turpin, G. (2019) *Mangrove
dieback in the Gulf of Carpentaria: baseline for monitoring future trajectory
in Queensland.* Queensland Herbarium, Department of Environment and Science.

Bunting, P. et al. (2022) 'Global mangrove extent change 1996–2020: Global
Mangrove Watch version 3.0', *Remote Sensing*, 14(15), 3657.

Davidson-Pilon, C. (2019) 'lifelines: survival analysis in Python', *Journal
of Open Source Software*, 4(40), 1317.

Duke, N.C. et al. (2017) 'Large-scale dieback of mangroves in Australia's Gulf
of Carpentaria', *Marine and Freshwater Research*, 68(10), 1816–1829.

Duke, N.C. et al. (2022) 'ENSO-driven extreme oscillations in mean sea level
destabilise critical shoreline mangroves', *PLOS Climate*, 1(8), e0000037.

Kennedy, R.E., Yang, Z. and Cohen, W.B. (2010) 'Detecting trends in forest
disturbance and recovery using yearly Landsat time series: 1. LandTrendr',
*Remote Sensing of Environment*, 114(12), 2897–2910.

Kennedy, R.E. et al. (2018) 'Implementation of the LandTrendr algorithm on
Google Earth Engine', *Remote Sensing*, 10(5), 691.

Meyer, H. and Pebesma, E. (2021) 'Predicting into unknown space? Estimating
the area of applicability of spatial prediction models', *Methods in Ecology
and Evolution*, 12(9), 1620–1633.

Mueller, N. et al. (2016) 'Water observations from space', *Remote Sensing of
Environment*, 174, 341–352.

Olofsson, P. et al. (2014) 'Good practices for estimating area and assessing
accuracy of land change', *Remote Sensing of Environment*, 148, 42–57.

Pontius Jr, R.G. and Millones, M. (2011) 'Death to Kappa: birth of quantity
disagreement and allocation disagreement for accuracy assessment',
*International Journal of Remote Sensing*, 32(15), 4407–4429.

Pontius Jr, R.G. and Santacruz, A. (2014) 'Quantity, exchange, and shift
components of difference in a square contingency table', *International
Journal of Remote Sensing*, 35(21), 7543–7554.

Roy, D.P. et al. (2016) 'Characterization of Landsat-7 to Landsat-8 reflective
wavelength and normalized difference vegetation index continuity', *Remote
Sensing of Environment*, 185, 57–70.

Stein, A., Aryal, J. and Gort, G. (2005) 'Use of the Bradley–Terry model to
quantify association in remotely sensed images', *IEEE Transactions on
Geoscience and Remote Sensing*, 43(4), 852–856.

Turner, H. and Firth, D. (2012) 'Bradley-Terry models in R: the BradleyTerry2
package', *Journal of Statistical Software*, 48(9), 1–21.

---
