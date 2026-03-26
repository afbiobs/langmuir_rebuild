# State Of Play — 2026-03-26

Purpose: prepare the repo for the next round of code changes and a cleanup pass.

This document is intentionally not a duplicate of the living documents. For the evolving scientific assumptions and modelling decisions, continue to use:

- `docs/assumptions_register.md`
- `docs/decision_log.md`
- `docs/failure_log.md`

## Current Position

- The current full-set benchmark reference is still the matched comparison rerun in `outputs/comparison_matched_current_code_rerun/`.
- The current Neagh diagnostic reference is now the grouped-point audit in `outputs/rank_audit/neagh_grouped_category_audit/`.
- The grouped-point Neagh audit uses `data/raw/point_summary_enriched.csv` directly, not the older 26-row benchmark observation subset.
- The fetch convention issue is considered resolved for current work: fetch is point-relative azimuth, and the working convention is meteorological `from`.
- The grouped-point Neagh audit now includes all 3 wiggle groups. The dated grouped set used in the latest run is:
  - `33 manual`
  - `18 stream`
  - `3 wiggle`
  - `2 manual` groups excluded because `observation_date` is missing

## What The Current Diagnostics Say

- The old “predictions clustering around ~40 m” behavior is not the active failure mode in the current code path.
- The current Neagh grouped-point audit suggests a split process:
  - `wiggle` observations are reasonably close to the raw onset scale
  - `manual` observations are much smaller than the current model outputs
- The current grouped audit supports the following working interpretation:
  - wiggle rows are being captured because they live on a larger spacing scale and all three wiggle cases have `n_events = 0`
  - manual rows are already too large at onset, then become much too large after one merger event
  - the final comparison layer currently restores the uncapped raw hierarchy width instead of preserving the mechanically capped width
- Treat the last point above as the main code-facing issue for the next step.

## Canonical Files To Keep

These are the files that should be treated as current working references.

- `data/raw/point_summary_enriched.csv`
- `outputs/rank_audit/rebuild_point_summary_enriched.py`
- `outputs/rank_audit/neagh_directional_rerun.py`
  - keep for now because the grouped audit imports helper functions from it
- `outputs/rank_audit/neagh_grouped_category_audit.py`
- `outputs/rank_audit/neagh_grouped_category_audit/neagh_grouped_category_rank_audit.csv`
- `outputs/rank_audit/neagh_grouped_category_audit/neagh_grouped_category_layer_summary.csv`
- `outputs/rank_audit/neagh_grouped_category_audit/neagh_grouped_category_wind_summary.csv`
- `outputs/rank_audit/neagh_grouped_category_audit/neagh_grouped_category_best_layer_summary.csv`
- `outputs/rank_audit/neagh_grouped_category_audit/observed_vs_layers_by_category.png`
- `outputs/rank_audit/neagh_grouped_category_audit/layer_bias_by_category.png`
- `outputs/rank_audit/neagh_grouped_category_audit/wind_duration_diagnostics.png`

## Files That Look Temporary Or Superseded

These should not be deleted blindly because the project is not currently in a git repo. Move/archive first, then delete only after the next clean rerun succeeds.

- `outputs/rank_audit/__pycache__/`
- `outputs/rank_audit/neagh_grouped_point_observations.csv`
  - temporary cache-fill input, fully regenerable
- `outputs/rank_audit/neagh_directional_rank_audit.csv`
  - superseded by grouped-point Neagh audit for class analysis
- `outputs/rank_audit/neagh_fetch_convention_sensitivity.py`
  - resolved question; probably archive after helper code is reviewed
- `outputs/rank_audit/neagh_fetch_convention_sensitivity/`
  - keep only if fetch-convention provenance is still needed
- `outputs/rank_audit/neagh_fetch_convention_sensitivity/prediction_error_by_stream_class.png`
  - stale proxy-class artifact
- `outputs/rank_audit/neagh_fetch_convention_sensitivity/predicted_vs_observed_by_fetch_convention.png`
  - superseded by category-aware outputs
- `outputs/rank_audit/neagh_fetch_convention_sensitivity/prediction_error_by_annotation_category.png`
- `outputs/rank_audit/neagh_fetch_convention_sensitivity/prediction_shift_to_minus_from.png`
- `outputs/rank_audit/neagh_fetch_convention_sensitivity/prediction_shift_to_minus_from_by_category.png`
- `outputs/rank_audit/neagh_fetch_convention_sensitivity/predicted_vs_observed_by_fetch_convention_and_category.png`
- `outputs/rank_audit/rank_audit.py`
- `outputs/rank_audit/rank_audit.csv`
  - still useful historically, but no longer the main Neagh-class diagnostic

## Recommended Cleanup Order

Because there is no `.git/` here, cleanup should be move-first, not delete-first.

1. Create an archive directory under `outputs/`, for example `outputs/archive_20260326/`.
2. Move superseded `outputs/rank_audit/` artifacts into that archive rather than deleting them.
3. Remove `__pycache__/` directories after confirming no scripts import from those paths.
4. Promote reusable audit logic out of `outputs/rank_audit/` and into a proper source location.
   - likely target: `src/evaluation/` or a new `src/tools/`
   - first candidate helpers to move: grouped-point loading, point-to-cache-row conversion, directional fetch lookup, grouped Neagh audit runner
5. Once helper code is promoted, stop importing production logic from scripts inside `outputs/`.
6. Only after the promoted code reruns cleanly, delete or archive the superseded scripts listed above.
7. After the audit tooling is stable, review the proliferation of old experiment directories under `outputs/comparison_matched_*`.
   - keep one clearly named canonical current comparison directory
   - archive one copy of historically important variant runs
   - move the rest out of the main `outputs/` surface

## Current Workflow To Use Going Forward

This is the shortest reproducible path for the current Neagh diagnostic workflow.

1. Rebuild grouped observations:
   - run `outputs/rank_audit/rebuild_point_summary_enriched.py`
   - source data come from `/home/op/PROJECTS/LN_Bathy_GEM/Data/copernicus_annotations*.csv`
2. Fill weather cache for grouped rows:
   - generate grouped observation rows as needed
   - use `python3 -m src.data.era5_cache` against those rows
3. Run grouped Neagh diagnostic:
   - run `outputs/rank_audit/neagh_grouped_category_audit.py`
4. Inspect only these outputs first:
   - grouped case table
   - layer summary
   - wind summary
   - observed-vs-layers plot
   - layer-bias plot
5. When prediction code changes, rerun the full benchmark comparison separately.
   - do not use the grouped Neagh audit as a substitute for the full-set regression check

## Immediate Technical Priorities

- Refactor the grouped audit so it no longer depends on helper imports from `outputs/rank_audit/neagh_directional_rerun.py`.
- Re-express the comparison layer in code so it is obvious whether the final prediction should use:
  - onset width
  - raw coarsened hierarchy width
  - mechanically capped width
  - some explicit observation-model transform
- Test the hypothesis that the old ~40 m attractor may have been closer to the manual-spacing process, while the current code is closer to the wiggle-spacing process.
  - treat this as a hypothesis to test, not a settled conclusion
- Keep manual and wiggle classes separated in all future diagnostics.
  - do not collapse them into one target spacing during debugging

## The One Sentence Handover

The repo should now be cleaned so that the grouped-point Neagh audit becomes the canonical diagnostic path, the old fetch-sensitivity detours are archived, and the next code change can focus on the observation/comparison layer plus the manual-vs-wiggle process split without working through stale outputs.
