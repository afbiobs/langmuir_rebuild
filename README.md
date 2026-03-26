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
