# Mangrove Resilience and Recovery in the Gulf of Carpentaria

Satellite time-series analysis of mangrove dieback and recovery following the
2015–16 mass mortality event in the Gulf of Carpentaria, northern Australia.

**ENST90007 Research Project · Semester 2, 2026 · University of Melbourne**
Supervisor: A/Prof Jagannath Aryal, Department of Infrastructure Engineering

---

## Research question

After the record 2015–16 mangrove dieback, why did some forests recover and
some not — and can recovery capacity be predicted spatially from satellite and
environmental data?

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

---

## Repository contents

### `/scripts` — Google Earth Engine (JavaScript)

| File | Purpose |
|---|---|
| `00_setup_test.js` | Environment check |
| `01_annual_composites.js` | Original Phase 1 pipeline (superseded — retained for the harmonisation before/after comparison) |
| `01b_annual_composites_harmonised.js` | **Phase 1 production script.** Roy et al. harmonisation, GMW union mask, diagnostics D1–D8 |
| `01c_tidal_screening.js` | MNDWI tidal screening test. Not retained — see decisions log #14 |
| `02_landtrendr.js` | **Phase 2 production script.** Self-contained: rebuilds composites, runs segmentation, produces the three-class map and stratum weights |

### `/docs` — Documentation

| File | Purpose |
|---|---|
| `1_decisions.docx` | Decisions log, 22 decisions with rationale and rejected alternatives |
| `2_methodology.docx` | Draft methodology, sections 3.1–3.7 |
| `3_questions.docx` | Supervisor briefing and open questions |
| `4_statistics.docx` | Statistical approach: reading of Pontius & Millones (2011) and Stein et al. (2005), and the proposed mechanism |
| `5_statistics_plain.docx` | Plain-language version of the above |

---

## Method summary

**Phase 1 — Time-series construction.** Landsat Collection 2 Tier 1 Level-2
surface reflectance, 1987–2025, from Landsat 5 TM, 7 ETM+, 8 OLI and 9 OLI-2.
Cloud, shadow, cirrus and fill masked via QA_PIXEL. Landsat 8/9 harmonised to
ETM+ equivalent space using Roy et al. (2016) coefficients. NBR computed per
image. Dry-season (May–October) median composites, masked to GMW v3 union
mangrove extent.

**Phase 2 — Disturbance detection.** LandTrendr temporal segmentation
(Kennedy et al. 2010) called directly on the harmonised composites. Two
products: an unrestricted greatest-loss disturbance history, and a
year-windowed 2015–16 dieback map.

**Phase 3 — Resilience metrics.** Magnitude of loss, recovery rate, recovery
ratio and years to recovery, derived per pixel from fitted segments. *In
progress.*

**Phase 4 — Driver modelling.** Random forest against climatic, hydrological
and geomorphic covariates, with spatially blocked cross-validation. *Not
started.*

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

---

## Reproducibility notes

Three points that are easy to get wrong and that materially affect results.

**Sensor harmonisation is required, and the LT-GEE module does not provide
it.** Module version 0.2.0 states on loading that it does not apply the Roy et
al. coefficients. Its `runLT()` wrapper therefore rebuilds composites without
harmonisation. This project calls
`ee.Algorithms.TemporalSegmentation.LandTrendr` directly on harmonised
composites instead. Verification: mean NBR over stable mangrove shows no step
change at either the 1999 or 2013 sensor transition. An unharmonised series
overstated post-2013 decline by ~0.19 NBR units and misidentified the year of
minimum condition by one year.

**NBR must be inverted before segmentation.** The algorithm assumes increasing
values represent disturbance; NBR falls when canopy is lost. Calling the
algorithm directly means the inversion is the caller's responsibility.

**Magnitude thresholds are in native index units here.** LT-GEE defaults
assume scaled integers; a threshold of 100 against native 0–1 reflectance is
unsatisfiable and returns zero detections.

---

## Data sources

| Dataset | Source |
|---|---|
| Landsat Collection 2 L2 SR | USGS, via Google Earth Engine |
| Global Mangrove Watch v3 | Bunting et al. (2022), GEE community catalogue |
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

Duke, N.C. et al. (2017) 'Large-scale dieback of mangroves in Australia's Gulf
of Carpentaria', *Marine and Freshwater Research*, 68(10), 1816–1829.

Duke, N.C. et al. (2022) 'ENSO-driven extreme oscillations in mean sea level
destabilise critical shoreline mangroves', *PLOS Climate*, 1(8), e0000037.

Kennedy, R.E., Yang, Z. and Cohen, W.B. (2010) 'Detecting trends in forest
disturbance and recovery using yearly Landsat time series: 1. LandTrendr',
*Remote Sensing of Environment*, 114(12), 2897–2910.

Lymburner, L. et al. (2020) 'Mapping the multi-decadal mangrove dynamics of
the Australian coastline', *Remote Sensing of Environment*, 238, 111185.

Pontius Jr, R.G. and Millones, M. (2011) 'Death to Kappa: birth of quantity
disagreement and allocation disagreement for accuracy assessment',
*International Journal of Remote Sensing*, 32(15), 4407–4429.

Roy, D.P. et al. (2016) 'Characterization of Landsat-7 to Landsat-8 reflective
wavelength and normalized difference vegetation index continuity', *Remote
Sensing of Environment*, 185, 57–70.

Stein, A., Aryal, J. and Gort, G. (2005) 'Use of the Bradley–Terry model to
quantify association in remotely sensed images', *IEEE Transactions on
Geoscience and Remote Sensing*, 43(4), 852–856.
