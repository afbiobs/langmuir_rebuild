# Benchmark Dataset Specification

**Status:** WP-02 deliverable. Defines the benchmark subsets used for model evaluation.
This document must not change once evaluation begins. Any revision requires a new
benchmark ID and a failure_log entry explaining why the old benchmark was retired.

**Source data:** `data/raw/observations.csv` (67 observations, spacing 37–518 m)

---

## Dataset Overview

| Field | Description |
|---|---|
| `observation_id` | Unique identifier (obs_001 … obs_067) |
| `image_date` | Date of satellite observation |
| `authoritative_lat/lng` | Location of the windrow pattern |
| `manual_spacing_m` | Manually annotated windrow spacing [m] (33 obs) |
| `wiggle_spacing_m` | Semi-automated spectral spacing estimate [m] (28 obs) |

Both spacing fields are present for 6 observations; neither is empty for all 67.
For observations where both are available, `manual_spacing_m` is treated as primary.
Where only `wiggle_spacing_m` is available, it is used as the target. No record is
excluded solely on the basis of which measurement method was used.

**Site distribution (approximate):**

| Approx. latitude | Longitude | Site identity |
|---|---|---|
| ~31°N | ~120°E | Taihu Lake, China |
| ~42°N | ~76–83°W | Lake Erie region |
| ~48°N | ~99°W | North American prairie lakes (Manitoba / North Dakota) |
| ~55°N | ~6°W | Lough Neagh, Northern Ireland (primary design site) |

---

## Benchmark Definitions

### BM-A: Full dataset (primary)

**Purpose:** Overall model performance; used for all headline metrics.
**Inclusion criteria:** All 67 observations with at least one spacing measurement.
**Target variable:** `manual_spacing_m` if available, else `wiggle_spacing_m`.
**N:** 67

This is the primary benchmark. All published model comparisons use BM-A unless
explicitly stated otherwise.

---

### BM-B: Cross-method agreement subset

**Purpose:** Assess whether `manual_spacing_m` and `wiggle_spacing_m` are
interchangeable, and whether model residuals differ systematically by method.
**Inclusion criteria:** Observations where both `manual_spacing_m` and
`wiggle_spacing_m` are non-empty.
**N:** 6

This subset is used diagnostically, not for headline scoring. If model residuals
differ significantly between the two measurement columns for the same observation,
this indicates the measurements are not measuring the same quantity and the
observation model (docs/observation_model.md §2) needs revision.

---

### BM-C: Near-onset / narrow-spacing

**Purpose:** Test model behaviour in the regime most sensitive to Ra, l_cNL, and the
linear-to-nonlinear transition. These are the cases where the CL physics matters most.
**Inclusion criteria:** Primary spacing < 65 m (below the first quartile, ~49 m).
Wait — use spacing < 75 m to give a reasonable subset size.

**Rationale:** For h = 9 m (Lough Neagh), spacing < 75 m → aspect ratio L/h < 8.3.
These are near-onset cells at l close to l_cNL. The model's ability to predict these
cases tests whether the CL instability scale (not just coarsening) is correct.

**N:** ~28 (to be confirmed by filtering the CSV at evaluation time)
**Warning:** If model performance on BM-C is substantially worse than BM-A, this
suggests the instability scale is wrong rather than the coarsening model.

---

### BM-D: Wide-spacing / coarsened cells

**Purpose:** Test the coarsening model. Wide spacings require one or more Y-junction
mergers to explain; BM-D probes whether the coarsening count and timescale are correct.
**Inclusion criteria:** Primary spacing > 150 m.
**Rationale:** For h = 9 m, spacing > 150 m → aspect ratio L/h > 16.7, which exceeds
the no-coarsening instability prediction. These observations require at least one
merger event in the model.
**N:** ~10 (estimated)
**Warning:** BM-D observations are disproportionately affected by coarsening
timescale uncertainty (AR-011) and the aspect ratio cap (AR-006). If the cap is
triggered frequently in BM-D, log it via failure_log.

---

### BM-E: Site-specific subsets

**Purpose:** Test whether the model generalises across sites. Different sites have
different depths, fetches, and wind climates.

| Subset | Lat filter | Approx. site | N (estimated) |
|---|---|---|---|
| BM-E-Taihu | 30°–33°N | Taihu Lake | ~25 |
| BM-E-NAmerica | 45°–51°N | Prairie lakes | ~18 |
| BM-E-Neagh | 53°–57°N | Lough Neagh | ~15 |

**Warning:** Site-level depth and fetch must be specified correctly for each subset.
Using Lough Neagh parameters (h = 9 m, fetch = 15 km) for Taihu Lake (h ≈ 2 m,
fetch ≈ 30 km) will silently corrupt predictions. The model must accept site-specific
depth and fetch as explicit inputs, never as hardcoded defaults.

---

## Evaluation Protocol

1. All benchmark definitions are frozen at this document's creation.
2. Model development uses BM-A for headline metrics. BM-C, BM-D, BM-E are
   diagnostic — they should not be used to tune model parameters.
3. Metrics (from `src/evaluation/metrics.py`):
   - RMSE and MAE in metres (absolute spacing error)
   - RMSE_ratio and MAE_ratio in percent (relative to observed spacing)
   - Bias (mean signed error; positive = model over-predicts)
   - Pearson r and Spearman ρ between predicted and observed spacing
   - Hit rate within ±30% of observed value
4. The old model (`lang-nonlin`) predictions on the same 67 observations form the
   baseline. Any new model must beat the old model on at least 3 of the 5 metrics
   above to be considered an improvement.
5. **Anti-saturation check:** Before reporting metrics, verify that the predicted
   spacing distribution spans at least 50% of the observed range. A model that
   predicts 80 ± 5 m for all 67 observations has a lower RMSE than a random predictor
   but is physically useless (see AGENTS.md §2).

---

## Data quality notes

- **No censored observations** in this dataset (all 67 have at least one spacing value).
  This is a selection artefact: the annotation process retained only images where a
  spacing could be measured. See docs/observation_model.md §4 on non-detections.
- **No forcing data in this file.** Forcing must be matched to observations via
  the ERA5 cache at `data/raw/era5_cache/` using the `image_date` and location.
  The matching procedure must document the temporal window used (e.g., 3-hour mean
  centred on overpass time).
- **Timing uncertainty.** The `image_date` is a date, not a timestamp. Overpass time
  must be inferred from the satellite track (MODIS ≈ 10:30 local, Landsat ≈ 10:00 local,
  Sentinel-2 ≈ 10:30 local). Propagate this ±30 min uncertainty in forcing-matched
  predictions where feasible.
