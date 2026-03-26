# Langmuir Analysis: Current State and Tuning Plan

**Project:** Clean-room rebuild of Langmuir circulation spacing prediction for shallow lakes.
**Phase:** Tuning and diagnostics — the physics engine is built, the task is now to understand and fix the observation mapping.
**Primary references:** `AGENTS.md` (governance), `docs/observation_model.md`, `docs/assumptions_register.md`, `docs/decision_log.md`, `docs/failure_log.md`

---

## 0. What Was Built

The system predicts satellite-visible windrow spacing from wind, depth, fetch, wave forcing, and turbulent mixing. It was built through six work packages (WP-01 to WP-06), all now complete.

### Architecture

```
Wind/waves/depth/fetch  →  ForcingState (Ra, La_t, u*, ν_T, profiles)
                                    ↓
                        CL nonlinear solver → l_cNL → L_inst = 2πh/l_cNL
                                    ↓
                        Coarsening (Y-junction mergers) → coarsened_width
                                    ↓
                        Mechanical cap (12×depth) → mechanical_cell_width
                                    ↓
                        Observation model → comparison_spacing
```

Two physics candidates (CL nonlinear solver, La_t scaling laws) plus three baselines (constant, linear-in-wind, depth-scaled).

### Repository Structure

```
langmuir_rebuild/
├── AGENTS.md                              # Coding governance (read first)
├── src/
│   ├── forcing/                           # Wind, waves, currents, eddy viscosity
│   │   ├── wind.py, waves.py, currents.py, eddy_viscosity.py
│   │   └── tests/
│   ├── hydro/                             # CL solver, coarsening, profiles, BCs
│   │   ├── rayleigh.py, profiles.py, robin_bc.py
│   │   ├── linear_solver.py, nonlinear_solver.py, galerkin.py
│   │   ├── scaling_laws.py, coarsening.py
│   │   └── tests/
│   ├── prediction/                        # Candidate pipelines, visibility, enhancement
│   │   ├── candidate_cl.py, candidate_scaling.py, baseline.py
│   │   ├── common.py, pipeline.py, visibility.py, lc_enhancement.py
│   │   └── tests/
│   └── evaluation/                        # Metrics, comparison, plots
│       ├── metrics.py, comparison.py, plots.py
├── tests/                                 # Integration, saturation, dimensional
├── data/raw/                              # Observations, ERA5 cache (67/67 populated)
├── outputs/                               # Benchmark comparisons, rank audits
└── docs/                                  # Assumptions, decisions, failures, references
```

---

## 1. Current Benchmark Results

Test suite: **95 passed**, 149 warnings. Matched weather cache: **67/67** cases.

### BM-A_full metrics (matched rerun)

| Model | RMSE [m] | Range coverage | Tail coverage | Spearman ρ |
|---|---:|---:|---:|---:|
| baseline_depth_scaled | 83.81 | 0.334 | 0.395 | — |
| baseline_linear_wind | 84.25 | 0.328 | 0.279 | — |
| baseline_constant | 91.43 | 0.000 | 0.163 | — |
| scaling | 126.02 | 0.180 | 0.047 | — |
| **cl** | **159.88** | **0.973** | **0.070** | **-0.464** |

**Key finding:** CL has excellent dynamic range (0.973 coverage) but **wrong ordering** — predictions are anti-correlated with observations. The model assigns large spacings to cases that should be small and vice versa.

---

## 2. The Diagnostic Picture

### Manual vs Wiggle Process Split

The grouped Neagh point audit (`outputs/rank_audit/neagh_grouped_category_audit/`) reveals that observations split into distinct classes with very different model behaviour:

| Category | n | Median obs [m] | Median L_inst [m] | Ratio | Median coarsened [m] | Median mech_cap [m] |
|----------|--:|---------------:|-------------------:|------:|---------------------:|--------------------:|
| wiggle | 3 | 214 | 228 | 0.94 | 228 | 188 |
| stream | 18 | 69 | 195 | 2.4 | 225 | 161 |
| manual | 33 | 43 | 196 | 3.9 | 388 | 166 |

**Wiggle** observations are well-matched by the raw onset scale (L_inst). No coarsening needed — all three wiggle cases have n_events = 0. The CL solver's onset width ≈ 228 m vs observed ≈ 214 m is a good match.

**Manual** observations are much smaller than any model layer. Even L_inst (the smallest model scale, before coarsening) is 3.9× the observed spacing. Coarsening makes it worse (7.8× at coarsened_width). The mechanical cap helps (3.3× at mechanical_cell_width) but not enough.

**Stream** observations fall between wiggle and manual.

### The comparison_spacing Bug

The comparison_spacing layer currently equals coarsened_width, not mechanical_cell_width. This means the 12×depth aspect-ratio cap is bypassed in the final comparison. Fixing this is the first code task.

### Layer Summary (Spearman ρ by category)

| Category | L_inst | coarsened | mech_cap | comparison |
|----------|-------:|----------:|---------:|-----------:|
| all | 0.18 | -0.22 | 0.10 | -0.22 |
| manual | 0.19 | 0.01 | 0.13 | 0.01 |
| stream | 0.13 | 0.07 | 0.14 | 0.07 |
| wiggle | 0.50 | 0.50 | 0.50 | 0.50 |

Coarsening degrades the rank ordering. The L_inst layer has weakly positive correlation that coarsening destroys.

---

## 3. Working Hypotheses

### H1: Coarsening initialization is too large for manual-class observations

The onset width L_inst ≈ 196 m for manual cases (depth ≈ 9 m) implies an aspect ratio of ~22 before any coarsening. This is already well above the 12× cap. The CL solver may be producing l_cNL values that are too small (cells too wide) for the conditions that produce manual-class observations.

### H2: Wiggle captures a different physical process

Wiggle observations (~214 m) match the raw onset scale with no coarsening. This may reflect newly-formed LC patterns at their instability wavelength. Manual observations (~43 m) may reflect a more evolved or different-scale process.

### H3: The fix direction is smaller initialization, not more coarsening

The user's insight: wiggle does better with no coarsening steps. Manual currently initializes with too large a value and then coarsens, pushing predictions further away. The productive direction is figuring out how to initialize at smaller sizes and not so quickly coarsen, rather than adding more coarsening machinery.

**Background:** It looks like there is a better link between the wiggle spacing process and the manual spacing process than previously thought, but the apparent difference may be partly a coarsening artefact. Wiggle does better with no coarsening steps. Manual initializes too large and then coarsens. If we can figure out how to initialize at smaller sizes and avoid premature coarsening, the two classes may converge.

### H4: comparison_spacing should use mechanical_cell_width

Fixing the comparison layer to use the capped width should improve results for manual cases, where the uncapped coarsened_width is 7.8× observed but mechanical_cell_width is 3.3× (still too large, but closer).

---

## 4. Investigation Tasks

These replace the old work package structure. They are not strictly sequential — some can be explored in parallel.

### Task A: Fix comparison layer (code fix)

**What:** Make comparison_spacing use mechanical_cell_width instead of uncapped coarsened_width.
**Where:** `src/prediction/common.py` or `src/evaluation/comparison.py` — wherever the final comparison spacing is assembled.
**Validation:** Rerun benchmark, check Spearman ρ by class. The fix should improve manual-class correlation.

### Task B: Onset width audit (diagnostic)

**What:** Understand why L_inst ≈ 196 m for typical Neagh conditions. At depth 9 m this is aspect ratio 22.
**Questions:**
- What l_cNL does the solver produce? What Ra, profiles, and BCs drive it?
- How does l_cNL vary across the 67 cases? Is there meaningful variation?
- Which forcing/profile conditions produce smaller l_cNL (tighter cells)?
- Are the Robin BC closures (derive_robin_bc_from_forcing) producing reasonable γ_s, γ_b?
**Method:** Dump l_cNL, Ra, γ_s, γ_b, profile coefficients for all 67 cases. Correlate with observed spacing.

### Task C: Layer-by-layer rank audit (diagnostic)

**What:** Identify which layer first breaks correct ordering.
**Method:** For each of the 67 cases, record the model quantity at each layer. Compute Spearman ρ(model, observed) at each layer, separately for manual/wiggle/stream.
**Key question:** Is forcing (Ra) already mis-ordered, or does the ordering break at the solver or coarsening stage?

### Task D: Coarsening sensitivity (diagnostic)

**What:** Determine whether the problem is onset width, merger speed, or both.
**Experiments:**
- Run with n_events forced to 0 (no coarsening) for all cases
- Run with onset width artificially halved
- Run with merger timescale doubled
- Compare Spearman ρ at each variant
**Key question:** If we skip coarsening entirely, does rank ordering improve?

### Task E: Profile and BC sensitivity (investigation)

**What:** The onset width l_cNL depends on the affine profiles D'(z), U'(z) and the Robin BCs γ_s, γ_b. Small changes in these can change l_cNL (and hence L_inst) significantly.
**Questions:**
- κ ranges from 1.425 to 1.944 depending on profiles — which profiles does the current code use for Neagh?
- Are the forcing-derived Robin BCs (AR-028) producing γ values in the literature-supported range?
- What profile/BC combination would produce L_inst ≈ 50 m at depth 9 m?

---

## 5. Decision Criteria for This Phase

### Primary metric: Spearman ρ

Rank correlation is more informative than RMSE at this stage. A model with correct ordering but a scale offset is much closer to success than one with low RMSE but wrong ordering.

### Separate evaluation by observation class

Every diagnostic and metric must be reported separately for manual, wiggle, and stream. Pooled metrics can mask class-specific problems.

### Parameter changes require justification

Per AGENTS.md §6.3 and §8.4: no parameter is changed solely because it improves a metric. Physical justification is required and must be documented in `docs/decision_log.md`.

### Success criteria

| Criterion | Target |
|-----------|--------|
| Spearman ρ (manual) | > 0 (positive correlation) |
| Spearman ρ (wiggle) | > 0 (positive correlation) |
| Spearman ρ (all) | > 0 (positive correlation) |
| Dynamic range preserved | Range coverage > 0.5 |
| All existing tests pass | 95+ tests |
| No new silent saturation | Saturation audit clean |

Reaching positive Spearman ρ across all classes with preserved dynamic range is the gate condition for promoting a production candidate.

---

## 6. Canonical Files and Workflow

### Current working references

- `data/raw/point_summary_enriched.csv` — grouped observation data
- `outputs/rank_audit/neagh_grouped_category_audit.py` — main diagnostic script
- `outputs/rank_audit/neagh_grouped_category_audit/` — current diagnostic outputs
- `outputs/comparison_matched_current_code_rerun/` — current benchmark reference

### Diagnostic workflow

1. Rebuild grouped observations: `outputs/rank_audit/rebuild_point_summary_enriched.py`
2. Fill weather cache: `python3 -m src.data.era5_cache`
3. Run grouped diagnostic: `outputs/rank_audit/neagh_grouped_category_audit.py`
4. Inspect: layer summary, wind summary, observed-vs-layers plot, bias plot
5. Rerun full benchmark comparison separately for regression check

### Key living documents

- `docs/assumptions_register.md` — current through AR-035
- `docs/decision_log.md` — current through 2026-03-24
- `docs/failure_log.md` — current through FL-021
- `docs/observation_model.md` — observation interpretation
- `docs/wp06_handover_2026-03-24.md` — detailed build-phase handover (archival reference)
