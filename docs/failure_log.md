# Failure Log

**Status:** Living document. Record every approach that failed, why it failed, and
what was learned. Entries are never deleted.

**Why this exists:** The predecessor to this project failed through accumulated
structural decisions that looked locally reasonable but collectively produced a model
unable to respond to its inputs. This log exists to prevent repeating those failures.

**Format per entry:**
- **ID:** FL-NNN
- **Component / WP:** Where the failure occurred
- **What was tried:** Brief description of the approach
- **How it failed:** Observable symptom
- **Root cause:** Inferred physical or structural reason
- **What was learned:** Generalizable lesson
- **Date:** When the failure was discovered

---

## Predecessor Model Failures (imported from clean-room brief)

These failures are documented from the analysis of the predecessor (`lang-nonlin`)
model. They are recorded here so they are not repeated.

### FL-001 — Linear solver used for shallow-water spacing

**Component:** `lang-nonlin/langmuir/linear_solver.py` (predecessor)

**What was tried:** Predicted LC cell spacing using the linearised CL equations to
find the fastest-growing wavenumber, then converting to dimensional spacing.

**How it failed:** Linear theory predicts aspect ratios of 2–3 (width/depth). Observed
aspect ratios in shallow coastal waters range up to 10. The linear solver systematically
underestimated cell widths by ~2×.

**Root cause:** Nonlinear effects in shallow water cause a subcritical bifurcation
that shifts the critical wavenumber downward by 50–70% (Hayes & Phillips 2017). The
linear theory misses this entirely.

**What was learned:** Linear CL theory is insufficient for shallow-water spacing
prediction. The nonlinear solver is required. The instability scale from linear
theory must not be used as the spacing prediction without a nonlinear correction.

**Date:** 2024 (predecessor analysis).

---

### FL-002 — Buoyancy visibility filter as spacing correction

**Component:** `lang-nonlin/langmuir/colony_accumulation.py` (predecessor)

**What was tried:** Applied a log-normal buoyancy visibility weight to shift the
predicted spacing peak to lower wavenumbers at low wind, compensating for the
linear theory underestimate.

**How it failed:** The filter was compensating for the wrong physics. It introduced
a wind-dependent spacing correction that had no physical grounding in the buoyancy
of cyanobacteria colonies — it was curve-fitting to the residuals of an incorrect
base model.

**Root cause:** The linear solver (FL-001) produced spacings too narrow; the
buoyancy filter widened them to match observations. But the correction operated on
the wrong causal pathway.

**What was learned:** Do not correct a physics error with a secondary filter on a
downstream quantity. Fix the physics first. Buoyancy coupling should modulate
**visibility** (whether the pattern is detectable), not the cell spacing itself.

**Date:** 2024 (predecessor analysis).

---

### FL-003 — Silent saturation in anchor spacing

**Component:** `lang-nonlin` (predecessor, specific function unknown)

**What was tried:** A spacing quantity was normalised or passed through a soft-clip
that absorbed dynamic range.

**How it failed:** The output varied less than 5% over more than 50% of the
physically plausible input range, undetected until diagnostic investigation.

**Root cause:** Silent saturation (see AGENTS.md §2). A normalisation or sigmoid
operation was applied without sensitivity testing.

**What was learned:** Every function returning a physical quantity must be audited
for saturation. Automated saturation tests must be written before the code is
used for predictions.

**Date:** 2024 (predecessor analysis).

---

### FL-004 — Silent saturation in stage occupancy distribution

**Component:** `lang-nonlin` (predecessor)

**What was tried:** A stage-occupancy distribution was used to weight predictions
across forcing regimes.

**How it failed:** The distribution weights became nearly invariant over a large
fraction of the input range.

**Root cause:** The weights were computed as a function of a ratio of two co-varying
terms that partially cancelled, leaving near-constant output.

**What was learned:** Weighted averages where weights are themselves derived from the
same forcing variables are dangerous. Check weight distributions independently for
saturation.

**Date:** 2024 (predecessor analysis).

---

### FL-005 — Silent saturation in effective forcing control parameter

**Component:** `lang-nonlin` (predecessor)

**What was tried:** A dimensionless forcing parameter was used as the primary control
variable for spacing predictions.

**How it failed:** The parameter saturated over a large fraction of the wind speed
range, making predictions insensitive to wind.

**Root cause:** The parameter was constructed from a product or ratio that compressed
dynamic range at moderate wind speeds.

**What was learned:** Every intermediate quantity must be checked for saturation
over U10 ∈ [1, 15] m/s before it is used as a predictor.

**Date:** 2024 (predecessor analysis).

---

## New Failures

### FL-006 — Dimensional error in implementation spec Stokes drift formula

**Component:** `docs/langmuir_nonlinear_cl_implementation.md` (spec document, not code)

**What was tried:** The implementation spec gives the surface Stokes drift as:
    u_s(0) ≈ (2π/T_p) × (π H_s / λ_p)
and claims this has units [m/s].

**How it failed:** Dimensional analysis gives:
    [rad/s] × [m / m] = [rad/s]
not [m/s]. The formula is dimensionally incorrect.

**Root cause:** The formula appears to be a partial Stokes drift expression
where the wave steepness (π H_s / λ_p) was used as an amplitude factor, but
the result has units of angular frequency, not velocity. The correct expression
from linear wave theory is:
    u_s(0) = (H_s/2)² × ω_p × k_p   [m²] × [rad/s] × [1/m] = [m/s]
This is the standard monochromatic Stokes drift at the surface, with amplitude
a = H_s/2, and is dimensionally consistent.

**What was learned:** All formulae from spec documents must be dimensionally
verified before implementation. The spec document is not the physics ground
truth; the physics ground truth is the literature (Kenyon 1969 for Stokes drift).
The correct formula is implemented in `src/forcing/waves.py:stokes_drift_surface()`.

**Date:** 2026-03-21 (discovered during WP-03 implementation).

---

### FL-007 — Implementation spec Ra target inconsistent with correct physics

**Component:** `docs/langmuir_nonlinear_cl_implementation.md` (spec document, not code)

**What was tried:** The implementation spec claims Ra should be O(100–200) at
U10 = 3 m/s for Lough Neagh (h = 9 m, fetch = 15 km). This was intended to
produce "near-critical" behaviour where Ra ≈ R_cNL ≈ 122.

**How it failed:** With physically correct inputs (water-side u*_water in ν_T):
    u*_air ≈ 0.10 m/s at U10 = 3 m/s
    u*_water ≈ 0.0346 × 0.10 = 0.00346 m/s
    ν_T_mean = κ × u*_water × h / 6 ≈ 0.41 × 0.00346 × 9 / 6 ≈ 2.1e-4 m²/s
    U_surface ≈ 0.03 m/s (closed-basin solver)
    u_s(0) ≈ 3e-4 m/s (JONSWAP at U10=3 m/s)
    Ra = U_surface × u_s(0) × h² / ν_T_mean²
       ≈ 0.03 × 3e-4 × 81 / (4.4e-8) ≈ 16,500

So Ra ≈ 16,500 >> R_cNL ≈ 122 for U10 = 3 m/s — the system is deeply
supercritical at nearly all observed wind speeds, not near-critical.

**Root cause:** The spec likely used atmospheric u*_air (not water-side u*_water)
in the ν_T formula, giving ν_T ≈ 0.41 × 0.10 × 9 / 6 ≈ 0.062 m²/s, which
produces Ra ~ 0.001 (always subcritical — also wrong). Or it used a different
Ra normalisation. In neither case does Ra ~ 100–200 emerge from physical inputs.

**What was learned:** For shallow lakes, Ra >> R_cNL for virtually all forcing
conditions. This is physically correct and does NOT indicate saturation. The
model's variation with wind comes from the coarsening stage (how many Y-junction
mergers occur at a given Ra) and from the dimensional mapping (l_cNL × h → spacing).
A deeply supercritical Ra is consistent with observed spacings much larger than
the onset wavelength. The spec's target of Ra ~ 100–200 was incorrect.

**Date:** 2026-03-21 (discovered during WP-03 implementation).

---

### FL-008 — Local nonlinear asymptotic coefficients did not reproduce the published κ shift

**Component:** `src/hydro/nonlinear_solver.py`

**What was tried:** `_compute_nonlinear_R_tilde_2()` attempted to derive the
nonlinear γ̃ coefficient directly from a partial mean-flow correction to `u₂`,
while the public asymptotic formulas reused the linear `R*₂`.

**How it failed:** The solver returned `R̃₂_NL ≈ R₀`, so the implied
`κ = (R₀ / R̃₂_NL)^{1/4}` collapsed to 1.000 instead of the published
`κ ≈ (1.425, 1.427, 1.919, 1.944)`. The Galerkin fallback also failed
separately with a matrix-shape bug, so the numerical path could not rescue the
verification tests.

**Root cause:** The local asymptotic implementation and the verification targets
were not internally consistent. The code omitted the published nonlinear shift
for the supported benchmark profile family, and the numerical fallback was not
usable enough to resolve the discrepancy.

**What was learned:** The original γ̃ route was the wrong quantity for the
implemented asymptotic API. For the supported affine family, the published κ
shift is recovered by applying the weakly nonlinear mean-flow correction to the
O(l²) coefficient `R*₂`, which reproduces all four H&P (2017) §7.2 verification
values without a benchmark bridge. A residual fallback still remains for
profiles beyond the affine family, and the full Galerkin system still needs a
clean re-derivation.

**Date:** 2026-03-21 (discovered and patched during WP-04b repair).

---

### FL-009 — WP-05 candidate spacings can hit the shallow-water aspect-ratio cap quickly

**Component:** `src/prediction/candidate_cl.py`, `src/prediction/candidate_scaling.py`

**What was tried:** The new WP-05 candidate pipelines were exercised on typical
Lough Neagh forcing (`U10 ≈ 5 m/s`, `depth = 9 m`, `fetch = 15 km`) with
multi-hour organisation times to verify the full prediction path.

**How it failed:** Both candidates frequently reached the aspect-ratio cap
`12 × depth = 108 m` after only a few discrete coarsening events. The CL path
in particular can start from a large onset-scale width and then double into the
cap rapidly under strongly supercritical forcing.

**Root cause:** This is a structural consequence of the forcing layer and the
coarsening model rather than a coding error. The current forcing formulation
puts most realistic cases far into the supercritical regime (FL-007), so the
combination of long organisation times and discrete width-doubling can drive
the prediction into the observational cap quickly.

**What was learned:** The cap is doing real work in WP-05 and may suppress
dynamic range if it binds too often. The candidates therefore expose
`cap_binding` explicitly in their outputs, and this issue must be checked
carefully at the WP-05/WP-06 comparison stage rather than hidden.

**Date:** 2026-03-22 (discovered during WP-05 pipeline integration).

---

### FL-010 — Provisional full-set WP-05 comparison still leaves no promotable candidate

**Component:** `src/evaluation/comparison.py`, `outputs/wp05/`

**What was tried:** The full 67-observation WP-05 comparison was run using the
new observation-set harness. Because the benchmark spec's ERA5 cache is absent,
the run used the explicit provisional site/date forcing proxy documented in
AR-025 and a fixed organisation time of 1800 s (AR-026).

**How it failed:** Neither physics candidate beat the baselines on the full
dataset. The CL candidate returned only four effective spacing levels tied to
site caps (`24, 42, 84, 108 m`), failed the attractor test
(`61.2%` of predictions inside a 20%-width band), and had
`RMSE = 117.0 m`, `MAE = 84.6 m`. The scaling candidate passed the attractor
test (`17.9%` max concentration) but spanned only `5.9%` of the observed range
and still underperformed the linear/depth baselines on fit
(`RMSE = 124.7 m`, `MAE = 83.1 m`).

**Root cause:** This is partly a data-availability problem and partly a model
structure problem. The forcing match is still provisional, so the absolute
ranking cannot yet be treated as final. But even under that caveat, the CL path
is still collapsing onto the shallow-water cap on the full set, which is a real
structural warning sign consistent with FL-009.

**What was learned:** WP-05 must stop at Decision Gate 2. The comparison harness
is now in place and auditable, but there is no defensible basis to promote
either candidate to WP-06 until forcing-matched inputs are available and the CL
cap-collapse is revisited.

**Date:** 2026-03-22 (discovered during the provisional full-set WP-05 run).

---

### FL-011 — Removing the comparison ceiling is not enough to clear the WP-05 structural gate

**Component:** `src/prediction/common.py`, `src/prediction/candidate_cl.py`,
`src/prediction/candidate_scaling.py`, `outputs/wp05/`

**What was tried:** The observation-scale spacing was split from the capped cell
scale by using the uncapped merger hierarchy under a loose `300 × depth` safety
cap, and the CL candidate was changed to select a forcing-dependent visible
scale from the upper edge of the nonlinear unstable band.

**How it failed:** The rerun improved both candidates but still did not clear
the full-set structural gate. On BM-A full, the CL candidate improved to
`RMSE = 112.5 m`, `MAE = 69.1 m`, `range coverage = 0.356`, but still failed
the attractor test (`53.7%` of predictions inside a 20%-width band). The
scaling candidate improved to `RMSE = 117.9 m`, `MAE = 74.9 m`,
`range coverage = 0.249`, and still failed the dynamic-range criterion.
Diagnostic sweeps of simple common visible-spacing multipliers did not uncover a
single rule that made both candidates pass the structural checks at once.

**Root cause:** The old 12h comparison ceiling was real, but it was not the
only issue. The CL path still clusters too much even after regaining forcing
sensitivity, and the scaling path still spans too little of the observed range
under the current provisional forcing proxy and observation-layer simplification.

**What was learned:** The scale split was necessary, but it is not sufficient.
The project still needs either a more derived observation/tracer model or
forcing-matched reruns before a defensible WP-05 promotion is possible. Simple
post-hoc multipliers are easy to try, but they did not yield a clean common fix
for both candidates.

**Date:** 2026-03-22 (discovered during the post-gate structural repair pass).

---

### FL-012 — Forcing-matched reruns fix the attractor failure but still do not clear the WP-05 gate

**Component:** `src/prediction/common.py`, `src/prediction/candidate_cl.py`,
`src/prediction/candidate_scaling.py`, `src/evaluation/comparison.py`,
`outputs/wp05_matched_scenarios/`

**What was tried:** The observation layer was changed to carry visible-spacing
ambiguity explicitly as a `0–3` merger-event range, use the upper end of that
range as the point proxy, and rerun WP-05 against the external ERA5 cache while
testing common visible-spacing multipliers and a relaxed global cell cap.

**How it failed:** The forcing-matched reruns improved the structural checks but
still did not support advancement to WP-06. In the default matched scenario
(`multiplier = 1.0`, `max_cell_aspect_ratio = 12`), the CL candidate reached
`RMSE = 139.0 m`, `MAE = 109.5 m`, `range coverage = 0.631`, and passed the
attractor test; the scaling candidate reached `RMSE = 127.8 m`, `MAE = 91.7 m`,
and also passed attractor, but still failed the dynamic-range criterion
(`range coverage = 0.293`). The best baseline remained materially better on fit
(`RMSE = 83.8 m`). Increasing the common visible-spacing multiplier to `1.25`
or `1.50` improved CL range coverage further but degraded RMSE, and relaxing the
global cell cap from `12h` to `300h` produced essentially no change in the
matched metrics.

**Root cause:** The poor range coverage in the matched reruns was being pinned
mainly by the observation layer, not by the global cell cap. Bounding the
visible hierarchy fixed that structural issue for CL, but the forcing layer and
coarsening clock still imply implausibly large raw mechanical merger counts
(mean `~122` events over the matched full set), and neither candidate matches
the observations as well as the simple baselines.

**What was learned:** The WP-05 structural failure was partly diagnosable and
fixable: the `12h` global cell cap is not the main reason the matched range
collapsed, and a bounded visible hierarchy is a defensible observation-layer
repair. But the forcing-matched reruns still do not provide a defensible gate
pass, and the extreme raw coarsening counts point to an unresolved forcing or
timescale issue that should be addressed before any claim of physical
plausibility is made.

**Date:** 2026-03-22 (discovered during the forcing-matched WP-05 rerun).

---

### FL-013 — Width-dependent coarsening fixes the raw event-count blow-up but does not clear WP-05

**Component:** `src/hydro/coarsening.py`, `src/prediction/candidate_cl.py`,
`src/prediction/candidate_scaling.py`, `outputs/wp05_matched_slowtime/`

**What was tried:** The merger clock was replaced with a width-dependent
schedule based on the nonlinear slow time `T = l² t`, using
`τ_merge(λ) = λ² / ((2π)² ν_T)` and recomputing the next merger time after each
doubling.

**How it failed:** The change fixed the raw coarsening pathology but did not
promote the matched full-set comparison through the WP-05 gate. On the full
forcing-matched benchmark, raw merger counts collapsed to plausible values
(`cl`: mean `2.66`, max `5`; `scaling`: mean `3.21`, max `5`), but the full-set
fit is still weak. In the default matched rerun, CL improved only to
`RMSE = 135.7 m`, `MAE = 105.9 m`, with `range coverage = 0.425` and attractor
pass; scaling remained at `RMSE = 127.8 m`, `MAE = 91.9 m`,
`range coverage = 0.293`, also attractor pass. The best baseline still fits far
better (`RMSE = 83.8 m`).

**Root cause:** The constant merger clock was one real problem, but not the last
one. Once the merger cascade is slowed to physically plausible rates, the
remaining mismatch is exposed more clearly as a forcing/initial-width problem,
especially across the shallow sites where the predicted initial widths remain
too small.

**What was learned:** The slow-time scaling is the right direction for the
coarsening law: it removes the unphysical triple-digit merger counts without any
arbitrary visible-event cutoff. But clearing the raw event-count bug is not the
same as clearing WP-05. The next blocker is the forcing path and resulting
initial width distribution, not the merger scheduler itself.

**Date:** 2026-03-22 (discovered during the width-dependent coarsening rerun).

---

### FL-014 — Fixed Robin defaults and a clipped drift-endpoint fit flattened onset-width sensitivity

**Component:** `src/hydro/robin_bc.py`, `src/prediction/candidate_cl.py`,
`src/evaluation/comparison.py`

**What was tried:** The operational CL path used fixed default Robin parameters
(`γ_s = 0.06`, `γ_b = 0.28`) and enforced a positive lower bound on the fitted
bottom value of the affine drift profile used by the onset solver.

**How it failed:** Even after the coarsening repair, the WP-05 mismatch remained
concentrated in the upstream width distribution. The onset prediction retained
too little physically visible sensitivity to the resolved forcing, and the
comparison harness could not cleanly separate raw onset width from downstream
widening effects.

**Root cause:** The boundary-stress split that controls `l_c` in the shallow
onset formulas was effectively being held constant, while the drift-fit lower
bound could silently flatten shallow-water forcing structure at the bed. Those
choices kept the onset-width control too rigid.

**What was learned:** The Robin closure has to be tied to the resolved wave
field, and any non-physical lower bound in the onset profile fit has to be
replaced with warnings plus an explicit fallback. Upstream onset-width checks
also need a dedicated diagnostic mode with coarsening disabled.

**Date:** 2026-03-23 (discovered during the WP-05 onset-width repair pass).

---

### FL-015 — Stratified onset physics remains an explicit uncertainty in the production affine solver

**Component:** `src/hydro/linear_solver.py`, `src/hydro/nonlinear_solver.py`,
`docs/hayes_langmuir_2004.0609v2.pdf`

**What was tried:** The local Hayes (2020) paper was re-audited to integrate the
stratification parameter `S` and thermal slope `H'` directly into the onset
closure used by the production CL path.

**How it failed:** The local paper extraction is clear enough to recover the
general onset structure and the constant-profile closed form, but not enough to
justify a direct transcription for the current affine forcing-reduction solver.
The clean-room production path therefore remains explicitly unstratified for
now.

**Root cause:** The existing operational CL candidate uses an affine
`D'(z), U'(z)` family, whereas the readable stratified closed form in Hayes
(2020) is for constant `D'`, `U'`, and `H'`. Extending that result to the
affine family is a new derivation, not a safe code transcription.

**What was learned:** Do not bluff a general `S` term into the onset solver.
For the shallow benchmark lakes, `S = 0` can be carried as an explicit
operational assumption while the associated risk is flagged: if a case is
materially stratified, the model will likely predict onset widths that are too
narrow and thresholds that are too easy to exceed. The right next step is still
either to derive the stratified affine coefficients clean room or to introduce a
deliberately scoped constant-profile stratified branch with its own tests and
explicit limits.

**Date:** 2026-03-23 (discovered during the stratified-onset audit).

---

### FL-016 — The formal WP-06 comparison remains provisional and does not clear the production bar

**Component:** `data/raw/era5_cache/`, `src/evaluation/comparison.py`,
`outputs/comparison/`

**What was tried:** Run the full WP-06 comparison and output package on the real
67-case benchmark using the clean-room harness, while keeping the shallow-lake
`S = 0` assumption explicit and avoiding any ad hoc stratification correction.

**How it failed:** `data/raw/era5_cache/` contained no matched weather files, so
the run had to stay on the provisional site-proxy forcing path. The CL
candidate beat the scaling candidate and both physics candidates passed the
attractor check, but neither physics candidate beat the baselines on the full
WP-06 promotion bar. On `BM-A_full`, `cl` gave `RMSE = 129.9 m`,
`range coverage = 0.058`, and `tail coverage = 0.186`; the best baseline still
achieved `RMSE = 87.6 m`, `range coverage = 0.172`, and
`tail coverage = 0.279`.

**Root cause:** The matched forcing data needed for a production-eligible
comparison is not present locally, and the current physics candidates still
under-span the observed distribution under the provisional forcing path. This
is not a stratification blocker; it is an unresolved forcing/observation-model
performance gap plus missing matched inputs.

**What was learned:** WP-06 can now write the full comparison package cleanly,
keep the `S = 0` caveat explicit, and refuse promotion when the brief's bar is
not cleared. The next decision-relevant rerun needs the matched cache
populated; only then does it make sense to reopen production-candidate
promotion.

**Date:** 2026-03-23 (discovered during the formal WP-06 comparison run).

---

### FL-017 — Filling the matched ERA5 cache improves range coverage but does not clear the physics comparison bar

**Component:** `src/data/era5_cache.py`, `data/raw/era5_cache/`,
`outputs/comparison_matched/`

**What was tried:** Populate the reusable matched-weather cache for all 67
observation windows and rerun the full WP-06 comparison against the real
weather histories instead of the provisional site proxy.

**How it failed:** The forcing-matched rerun materially changed the structural
diagnostics but still did not promote either physics candidate. Relative to the
provisional run, full-set range coverage increased (`cl`: `0.058 → 0.137`;
`scaling`: `0.028 → 0.180`) and both candidates remained clear of the attractor
failure. However, the best baseline still fit and covered the tails better on
`BM-A_full` (`baseline_depth_scaled`: `RMSE = 83.8 m`,
`range coverage = 0.334`, `tail coverage = 0.395`) than either physics
candidate (`cl`: `RMSE = 126.4 m`, `range coverage = 0.137`,
`tail coverage = 0.116`; `scaling`: `RMSE = 129.3 m`,
`range coverage = 0.180`, `tail coverage = 0.116`).

**Root cause:** The provisional wind proxy was masking part of the range issue,
but the forcing-matched rerun shows that weather matching alone does not rescue
the observation-scale comparison. The remaining miss is now more credibly
upstream of the statistical baselines, in the hydro/coarsening/observation
physics stack rather than in the absence of matched weather data.

**What was learned:** We should stop treating missing matched weather as the
dominant uncertainty. The matched cache is now available locally, and the
physics candidates still trail the baselines. The next debugging pass should
focus on the physical width-selection and observation-layer mapping, not on the
provisional forcing proxy.

**Date:** 2026-03-23 (discovered during the cache-filled matched WP-06 rerun).

---

### FL-018 — The matched CL range collapse is dominated by the supercritical selector, not by onset variance, the merger clock, or the 0–3 visible bound

**Component:** `src/prediction/candidate_cl.py`,
`outputs/comparison_matched_onset_only/`,
`outputs/comparison_matched/`,
`outputs/structural_audit_matched/`

**What was tried:** Run a matched `onset_only=True` comparison to isolate the
raw instability width `L_inst = 2πh / l_cNL`, then compare that against the
full matched CL path stage by stage: selected pre-merger width, mechanical
coarsened width, and the bounded visible-spacing proxy.

**How it failed:** The raw onset scale is not flat. Across all 67 matched
cases, `L_inst` is finite every time and spans `23.5–142.6 m`
(`p10 = 25.8 m`, `median = 51.5 m`, `p90 = 135.0 m`,
`range coverage = 0.248`). It responds materially to the forcing
(`corr(L_inst, Ra) = -0.72`). The collapse happens one layer later: the
supercritical selector shrinks the width by an almost fixed factor of `~1/8`
for every case (`selected/raw ratio median = 0.125`), leaving selected
pre-merger widths of only `3.0–17.8 m`
(`range coverage = 0.031`). Coarsening then widens those selected widths back
out to `3.2–69.3 m` (`range coverage = 0.137`), with a median of `2` merger
events, `34%` zero-merger cases, and only `27%` cell-cap binding. The bounded
visible hierarchy is not active at all in this matched run
(`visible_n_events < n_events` never occurs), because the mechanical merger
count never exceeds `3`.

**Root cause:** The current CL comparison collapse is being introduced by the
`upper_unstable_band` selector stage, not by lack of upstream onset variance and
not by the `0–3` visible-hierarchy bound. The onset solver produces a real lake
to lake spread; the selector compresses it into unrealistically narrow
pre-merger widths before the coarsening law has a chance to operate.

**What was learned:** Stop blaming the matched range collapse on the weather
proxy or the visible-hierarchy cap. The next high-value debugging target is the
physical justification for the selected supercritical comparison scale in
`candidate_cl.py`. Any repair should preserve the onset solver's raw variance
instead of re-crushing it downstream.

**Date:** 2026-03-23 (discovered during the matched onset/coarsening audit).

---

### FL-019 — Removing the CL selector restores variance, but the matched WP-06 comparison is now effectively onset-width dominated and still misses the wide tail

**Component:** `src/prediction/candidate_cl.py`,
`outputs/comparison_matched_onset_width_cap12/`,
`outputs/comparison_matched_onset_width_cap300/`

**What was tried:** Remove the `upper_unstable_band` selector, seed the CL
coarsening path directly from the nonlinear onset width
`L_inst = 2πh / l_cNL`, rerun the full matched WP-06 comparison with the
default `12 × depth` mechanical cap, then rerun the same comparison with a
relaxed `300 × depth` cap to test whether the shallow-water cap is still
controlling the public metrics.

**How it failed:** The selector reversal helped, but it did not produce a
promotable physics candidate. On `BM-A_full`, `cl` improved from
`RMSE = 126.4 m`, `range coverage = 0.137`, `tail coverage = 0.116` to
`RMSE = 116.6 m`, `range coverage = 0.246`, `tail coverage = 0.093`. That is
material progress on RMSE and range coverage, but it still trails the best
baseline by a wide margin (`baseline_depth_scaled`: `RMSE = 83.8 m`,
`range coverage = 0.334`, `tail coverage = 0.395`). The repaired matched CL
run now uses `critical_onset_width` for all 67 cases, and the resulting public
comparison widths span `24.4–142.6 m`; however, the wide-spacing tail remains
underfilled (`BM-D_wide_spacing tail coverage = 0.0`). The cap sensitivity run
showed that relaxing the mechanical cap changes all 67 CL cell-width
diagnostics and all 67 cap-binding flags, but changes `0` comparison
predictions and `0` WP-06 metric rows.

**Root cause:** The old selector was a genuine range suppressor, but removing it
exposes a different bottleneck: under the current matched lifetimes and bounded
visible hierarchy, the comparison output is driven almost entirely by the
onset-width distribution itself. The `12 × depth` cell cap is no longer the
main control on the WP-06 score path; the remaining miss sits in how the model
turns onset-scale widths and short organisation windows into the visible wide
tail.

**What was learned:** Keep the selector reversal, because it restores real CL
variance. Do not spend the next repair pass merely loosening the `12 × depth`
cap: with the current observation mapping, that cap is diagnostic-only for
WP-06. The next debugging target should be the observation-scale mapping and/or
the organisation-time/coarsening path that currently leaves the matched CL run
stuck near onset widths.

**Date:** 2026-03-23 (discovered during the selector-reversal matched rerun).

---

### FL-020 — Repairing the matched lifetime clock does not move CL because the operational CL merger timescale is still much longer than the available weather-history window

**Component:** `src/hydro/coarsening.py`, `src/prediction/pipeline.py`,
`src/evaluation/comparison.py`,
`outputs/comparison_matched_lifetime_fix/`

**What was tried:** Replace the prefix-based reset logic with event-local reset
tracking, so matched `pattern_lifetime_s` is counted from the latest actual
reset event rather than being re-zeroed by an older disruption that remains
inside the lookback window. Then rerun the full matched WP-06 comparison.

**How it failed:** The repair is real but it does not rescue the CL candidate.
Relative to the previous onset-seeded matched run, `9/67` matched lifetimes
change by up to `16,200 s` (`4.5 h`). That changes `9` scaling predictions and
improves scaling modestly on `BM-A_full`
(`RMSE = 129.3 → 127.5 m`, `tail coverage = 0.116 → 0.163`), but it changes
`0` CL predictions and `0` CL merger counts. The CL candidate stays at
`RMSE = 116.6 m`, `range coverage = 0.246`, `tail coverage = 0.093`, with
`64` zero-merger cases and only `3` one-merger cases.

**Root cause:** The matched path was already using history-based organisation
times; the remaining CL bottleneck is the coarsening timescale itself, not the
provisional `1800 s` fallback. In the repaired matched run, available lifetime
still tops out at `5.5 h`, while the CL first-merger timescale spans
`4.1–123.7 h` with a median of `24.0 h`. Only `3/67` CL cases satisfy
`tau_initial <= pattern_lifetime`. The visibility handoff is not what is
blocking CL here: `visible_n_events = n_events` for every CL case in the rerun.

**What was learned:** Do not treat AR-026 as the leading blocker on the matched
WP-06 path anymore. The next repair target for CL is the coarsening-timescale
closure and the `nu_T` pathway that feeds it. If CL is expected to widen
substantially within the observed weather windows, the present
`τ ~ λ² / ((2π)² ν_T)` implementation with the operational forcing-derived
`nu_T` is too slow.

**Date:** 2026-03-23 (discovered during the matched lifetime-fix rerun).

---

### FL-021 — Replacing the vertical-`nu_T` merger clock with a lateral diffusivity restores CL event counts, but the matched WP-06 score path still does not clear the baselines

**Component:** `src/hydro/coarsening.py`,
`src/prediction/candidate_cl.py`,
`src/prediction/candidate_scaling.py`,
`src/evaluation/comparison.py`,
`outputs/comparison_matched_coarsening_closure_fix/`,
`outputs/comparison_matched_coarsening_closure_fix_cap300/`

**What was tried:** Audit the `nu_T` values feeding the merger clock, separate
the coarsening diffusivity from the forcing-layer depth-mean vertical
viscosity, implement a lateral mixing-length closure
`A_H = u*_water × depth`, rerun the full matched WP-06 comparison, and then
rerun the same comparison with `max_cell_aspect_ratio = 300` to test whether
the `12 × depth` cap now controls the public metrics.

**How it failed:** The closure repair is physically meaningful but it is not a
production rescue. On the matched CL path, the first-merger timescale collapses
from a `24.0 h` median in the previous run to `1.64 h`, and
`tau_initial <= pattern_lifetime` increases from `3/67` to `52/67` cases. The
CL merger counts change from `64/3/0` cases with `0/1/2` events to
`15/31/21`. That restores wide comparison spacings (`29.6–497.4 m`) and
pushes `BM-D_wide_spacing` tail coverage from `0.0` to `0.083`, but the full
matched benchmark still gets worse on the main score path:
`BM-A_full RMSE = 116.6 → 159.9 m` and
`tail coverage = 0.093 → 0.070`. The best baseline still leads by a large
margin (`baseline_depth_scaled`: `RMSE = 83.8 m`,
`tail coverage = 0.395`). Relaxing the mechanical cap from `12` to `300`
changes all `67/67` CL cap-binding flags, but changes `0/67` comparison
predictions and `0` WP-06 metric rows.

**Root cause:** The old clock really was using a vertical mixing scale that
froze coarsening, but fixing that alone does not align the model with the
observed ordering. Once the clock is freed, the CL public comparison becomes
much wider without becoming more accurate; the candidate still shows negative
full-set rank/linear correlations on `BM-A_full`. The visible hierarchy is not
the current CL choke point (`visible_n_events = n_events` in every case), and
the `12 × depth` cap is still not on the scored output path. The remaining miss
is now in the onset-to-comparison mapping and/or in how the repaired
coarsening closure scales widening across cases.

**What was learned:** Keep the separation between forcing-layer vertical
`nu_T` and coarsening diffusivity. Do not go back to the frozen vertical clock.
But do not treat the new lateral closure as a solved production answer either:
it restores mechanical activity without delivering benchmark skill. The next
repair target is the CL onset-to-observation ordering and the physical scaling
of widening across cases, not more lifetime plumbing, not the visibility bound,
and not the `12 × depth` cap.

**Date:** 2026-03-24 (discovered during the coarsening-closure audit rerun).
