# Langmuir Rebuild

Clean-room rebuild of the Langmuir spacing workflow for shallow lakes.

## What this repo currently does

The active model workflow predicts Langmuir-spacing diagnostics from wind, depth,
fetch, wave forcing, and turbulent mixing, then compares those predictions to
observed spacing from imagery.

The important distinction is:

- `L_inst`: raw instability-selected onset width
- `coarsened_width`: raw post-merger width before the mechanical cap
- `comparison_spacing`: the final observation-facing spacing used in benchmark comparison

This repo also carries LC enhancement and visibility diagnostics, but the current
development focus is the spacing pipeline and its comparison to observations.

## Current state

- The current full benchmark reference run is in `outputs/comparison_matched_current_code_rerun/`.
- The current Neagh diagnostic reference is the grouped-point audit in `outputs/rank_audit/neagh_grouped_category_audit/`.
- The current wiggle-only cross-lake diagnostic is in `outputs/rank_audit/wiggle_crosslake_audit/`.
- The main unresolved issue is the observation/comparison layer:
  - wiggle observations can align with the current model in some cases
  - manual observations are generally much smaller than the current predictions
  - the current comparison spacing is still restoring the uncapped raw hierarchy rather than preserving the mechanically capped width

## What not to assume from older outputs

- The old "~40 m attractor" is not the current failure mode.
- Neagh analysis is no longer limited to the older benchmark-only subset; grouped point annotations are now part of the active audit workflow.
- Fetch handling for the Neagh point audit is now point-relative and uses meteorological wind direction as `from`.
- Some files in `outputs/rank_audit/` are historical or superseded. Use the canonical output directories above first.

## Current workflow

1. Rebuild grouped Neagh point observations with `outputs/rank_audit/rebuild_point_summary_enriched.py`.
2. Fill matched weather cache with `python3 -m src.data.era5_cache` when needed.
3. Run the relevant audit:
   - full benchmark comparison: `outputs/comparison_matched_current_code_rerun/`
   - Neagh grouped-point audit: `outputs/rank_audit/neagh_grouped_category_audit.py`
   - wiggle-only cross-lake audit: `outputs/rank_audit/wiggle_crosslake_audit.py`
4. Check the layer-by-layer outputs before changing model code.

## How to run the analyses

### Unit tests

Run all unit tests from the repo root:

```bash
python3 -m pytest src/ -v
```

Run only the multiscale tests:

```bash
python3 -m pytest src/hydro/tests/test_multiscale_structures.py src/prediction/tests/test_candidate_multiscale.py -v
```

Run only the CL candidate and pipeline tests:

```bash
python3 -m pytest src/prediction/tests/test_pipeline.py src/prediction/tests/test_diagnostics.py -v
```

Run forcing and hydro support tests:

```bash
python3 -m pytest src/forcing/tests/test_forcing.py src/hydro/tests/test_support_modules.py src/hydro/tests/test_literature_values.py -v
```

### Audit scripts

All audit scripts are run from the repo root with `python3`:

```bash
# 1. Rebuild grouped observation inputs (required before audits)
python3 outputs/rank_audit/rebuild_point_summary_enriched.py

# 2. Fill the ERA5 weather cache (requires network; only needed once per new observation set)
python3 -m src.data.era5_cache

# 3a. Neagh grouped-point audit (primary Neagh diagnostic)
python3 outputs/rank_audit/neagh_grouped_category_audit.py

# 3b. Wiggle-only cross-lake audit
python3 outputs/rank_audit/wiggle_crosslake_audit.py
```

To run the multiscale candidate instead of the CL candidate in audit scripts,
pass `candidate="multiscale"` to `analyse_case()` in the pipeline. The pipeline
entry point (`src/prediction/pipeline.py`) accepts `candidate="cl"`,
`candidate="scaling"`, or `candidate="multiscale"`.

### Analyses that may now be redundant

With the addition of the multiscale candidate (`candidate_multiscale.py`), the
following analyses may be partially redundant or superseded:

- **`outputs/rank_audit/neagh_fetch_convention_sensitivity.py`** — resolved the
  fetch-convention question; results are archived and not needed for ongoing work.
- **`outputs/rank_audit/rank_audit.py`** — the original single-class Neagh audit,
  superseded by the grouped-point category audit
  (`neagh_grouped_category_audit.py`).
- **`outputs/rank_audit/neagh_directional_rerun.py`** — kept only because the
  grouped audit imports helper functions from it. Should be refactored out.
- **CL-candidate-only onset-spacing diagnostics** — the CL onset scale
  (`L_inst = 2πh/l_cNL`) is now a diagnostic reference in the multiscale
  candidate, not the primary prediction. Audits that compare observations
  directly against the CL onset scale are still useful for wiggle-class
  observations but are not the main prediction path for manual-class spacings.
- **Old `outputs/comparison_matched_*` variant directories** — only
  `outputs/comparison_matched_current_code_rerun/` is canonical. Earlier variant
  runs are historical.

## Key files

- `AGENTS.md`: coding rules for this project
- `docs/assumptions_register.md`
- `docs/decision_log.md`
- `docs/failure_log.md`
- `docs/observation_model.md`
- `state-of-play-20260326.md`: current cleanup and handover note

## Scope boundary

This project is not a bloom model. It does not predict biomass, toxins, or scum
distribution. The current scientific question is how forcing maps to Langmuir
spacing and how that spacing should be compared to different observed pattern
classes.
