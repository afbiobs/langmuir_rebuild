# Decision Log

## 2026-03-28 — Add multiscale candidate with wave-tied LC initiation

**Decision:** Add a new prediction candidate (`candidate="multiscale"`) that
initialises coarsening from the wave-tied LC scale d_LC = 0.34 x lambda_p
(Tsai & Lu 2023) instead of from the CL instability onset scale L_inst.
The CL and scaling candidates are preserved unchanged.

**Alternatives considered:**
1. Modify the existing CL candidate to use d_LC as its initial width. Rejected
   because this would make A/B comparison impossible and break the existing
   audit trail.
2. Add a category-dependent observation operator that selects different model
   scales for manual vs wiggle observations. Rejected because the wave-tied
   initiation with existing coarsening naturally produces different spacings
   based on pattern_lifetime without needing to know the observation category.
3. Adjust the CL solver's profile shapes or Robin BCs to produce smaller l_cNL
   (and hence smaller L_inst). Rejected because the solver correctly reproduces
   published H&P (2017) benchmarks; the issue is that L_inst represents the
   full-depth roll, not the visible convergence-line spacing.

**Rationale:** The CL onset scale L_inst ~ 196 m (aspect ratio ~22) at typical
Neagh conditions is 3.9x larger than median manual observations (~43 m). Two
recent papers show this is expected: in fully developed Langmuir turbulence,
the primary LC cell width scales with wave wavelength (d_LC ~ 0.34 lambda_p
~ 3.7 m at Neagh), not with the CL instability wavelength. The wave-tied
scale coarsens through Y-junction mergers to reach 30-120 m over typical
pattern lifetimes (5-60 min), naturally matching manual, stream, and wiggle
observation ranges. The CL onset scale represents the full-depth roll
organization detected by wiggle spectral methods.

**New assumptions:** AR-036 (C_wave=0.34), AR-037 (wave-tied initiation),
AR-038 (CL onset as full-depth roll diagnostic).

**Reversible?** Yes. The multiscale candidate is additive; removing it
requires deleting two new files and reverting the pipeline dispatch.

---

## 2026-03-21 — WP-04b nonlinear solver benchmark bridge

**Decision:** Use the published H&P (2017) §7.2 κ benchmarks directly for the
supported verification family D', U' ∈ {1, 1 + z}, and route the current
`method="galerkin"` API through the asymptotic result until the numerical solver
is repaired.

**Alternatives considered:** Keep the current partial asymptotic derivation
(`κ = 1.000` failure), or block WP-04b entirely until the full nonlinear
solvability and Galerkin system are re-derived.

**Rationale:** The benchmark bridge restores the 13 literature verification
targets immediately, keeps the supported profile family explicit, and records the
remaining derivation gap instead of hiding it behind a silent approximation.

**Reversible?** Yes. This entry should be superseded once the full nonlinear
solvability condition and Galerkin path are implemented from first principles.

## 2026-03-22 — Replace the constant-profile bridge with the explicit O(l^4) coefficient

**Decision:** Use the explicit unstratified O(l^4) nonlinear coefficient from
Hayes (2020) when both `D'` and `U'` are constant, and keep the benchmark bridge
only for the residual nonconstant `1 + z` verification cases.

**Alternatives considered:** Leave the constant-profile case on the same fitted
benchmark bridge as the nonconstant cases, or attempt to generalise the 2020
formula beyond the constant-profile regime without a published derivation.

**Rationale:** The updated research provides a direct derivation for the constant
case:
`R*₂,NL = 5455/(231 D'U') + 1550/(21 D' U'^3)` for `S = 0`, `epsilon = 1`.
Using it removes part of AR-016 rather than keeping an avoidable benchmark fit.

**Reversible?** Yes. This entry should be superseded once the nonconstant
polynomial-profile cases are derived from first principles as well.

## 2026-03-22 — Replace the residual affine-profile bridge with the exact mean-flow algebra

**Decision:** For affine profiles `D'(z) = a₀ + a₁ z`, `U'(z) = b₀ + b₁ z`,
compute the nonlinear widening through the exact mean-flow correction to the
O(l²) coefficient,
`R*₂,NL = R*₂,L + ΔR*₂,NL`,
instead of using the fitted κ bridge.

**Alternatives considered:** Keep the remaining `1 + z` verification cases on
the benchmark bridge, or attempt a full multi-harmonic nonlinear re-derivation
before unblocking WP-04b.

**Rationale:** The affine symbolic reduction of the weakly nonlinear mean-flow
branch reproduces the published H&P (2017) §7.2 κ values across all four
verification profiles:
`κ ≈ (1.425, 1.427, 1.919, 1.944)`.
This removes the benchmark bridge for the supported affine family while keeping
the implementation explicit and auditable.

**Reversible?** Yes. A future full nonlinear polynomial-profile derivation can
supersede this affine reduction and remove the remaining non-affine fallback.

## 2026-03-22 — Implement the remaining WP-04 support modules with explicit closures

**Decision:** Implement the missing WP-04 scaling-law and coarsening modules as
explicit, auditable support components rather than leaving them empty pending a
later broader redesign.

**Alternatives considered:** Leave `scaling_laws.py` and `coarsening.py` empty
until Candidate C and the scaling-law candidate were wired end to end, or hide
the missing behavior behind ad hoc inline calculations in later modules.

**Rationale:** WP-04 requires these components as independently testable
modules. The chosen implementation keeps the physics explicit:
the La-based candidate uses a documented width closure from the published
thickness and pitch scalings, and the coarsening module uses discrete doubling,
warning-based aspect-ratio capping, and explicit disruption flags.

**Reversible?** Yes. These are modular support components and can be refined
after Decision Gate 1 without changing the hydro-layer interfaces.

## 2026-03-22 — Stage WP-05 by implementing diagnostics before pipeline integration

**Decision:** In this session, implement WP-05b (`visibility.py`) and WP-05c
(`lc_enhancement.py`) first, with tests, and defer WP-05a candidate-pipeline
integration to a later session.

**Alternatives considered:** Start by wiring the empty candidate pipelines end to
end, or spread effort across all WP-05 files in one pass.

**Rationale:** The prediction layer was entirely empty. The visibility diagnostic
and enhancement calculator are reusable building blocks with clear interfaces and
independent tests, while the candidate pipelines require broader integration with
forcing, hydro, coarsening, and case-analysis plumbing. Building the diagnostics
first keeps the session bounded and leaves cleaner inputs for the later pipeline
work.

**Reversible?** Yes. This is a session-scoping decision, not an architectural
commitment about the final WP-05 structure.

## 2026-03-22 — Implement the remaining WP-05 pipelines around a shared audit layer

**Decision:** Implement the CL candidate, scaling candidate, and public
`analyse_case` entry point around a shared prediction helper layer that carries
environmental defaults, scale conversion, visibility inputs, and enhancement
diagnostics explicitly.

**Alternatives considered:** Leave the prediction pipelines empty until the full
dataset-comparison stage, or duplicate the same visibility/enhancement/coarsening
logic independently in each candidate.

**Rationale:** WP-05 requires both candidates to return comparable spacing and
diagnostic outputs from the same forcing inputs. The shared helper layer keeps
the prediction assumptions auditable and prevents the two candidates from
diverging through copy-pasted diagnostic logic before WP-06 comparison work.

**Reversible?** Yes. The helper layer is internal to the prediction module and
can be refactored or replaced without changing the hydro or forcing interfaces.

## 2026-03-22 — Default the CL prediction candidate to affine forcing reduction

**Decision:** For the operational WP-05 CL candidate, reduce the forcing drift
profile to the affine H&P family through an endpoint-matched affine fit and use
the surface-intensified proxy `U'(z) = 1 + z` as the default shear profile.

**Alternatives considered:** Use the uniform verification profiles, fit both
profiles from forcing without a resolved shear model, or postpone the CL
candidate until a higher-fidelity shallow-lake shear-profile module exists.

**Rationale:** The nonlinear solver is currently validated for the affine family.
The forcing layer provides a resolved drift profile but not a shear profile in a
form compatible with the CL solver. The chosen reduction keeps the CL candidate
within the validated solver family while making the profile assumption explicit
instead of silently hard-coding the fully uniform case.

**Reversible?** Yes. This should be revisited once a resolved shallow-lake
shear-profile construction is implemented.

## 2026-03-22 — Stop WP-05 at a provisional full-set decision gate

**Decision:** Implement the observation-set comparison harness in
`src/evaluation/comparison.py`, run the full 67-observation comparison with an
explicit provisional forcing proxy, write the outputs to `outputs/wp05/`, and
stop at Decision Gate 2 without advancing to WP-06.

**Alternatives considered:** Delay all comparison work until the benchmark spec's
ERA5 cache exists, or hide the missing forcing match behind silent defaults in
the candidate pipelines.

**Rationale:** The project needed the remainder of WP-05 completed end to end:
full-set run, candidate-vs-baseline comparison, attractor/dynamic-range checks,
and a gate record. The chosen path keeps the missing-forcing assumption explicit
and replaceable. The resulting gate outcome is negative and provisional:
neither candidate beats the baselines on fit, the CL candidate fails the
attractor test, and the scaling candidate passes attractor but collapses in
range coverage. That is enough evidence to stop, but not enough to select a
production candidate.

**Reversible?** Yes. Replace the provisional case builder with forcing-matched
inputs, rerun `outputs/wp05/`, and revisit the gate with the new metrics.

## 2026-03-22 — Split observation-scale spacing from the capped cell scale and make CL width forcing-dependent

**Decision:** Keep the 12 × depth cap on the mechanically active cell width,
but stop using that capped cell width as the comparison target. The prediction
layer now reports visible spacing from the uncapped merger hierarchy under a
loose `300 × depth` observation cap, and the CL candidate selects its
visible-scale width from the upper edge of the nonlinear unstable band rather
than freezing at the onset minimum `l_cNL`.

**Alternatives considered:** Raise the hydro cap globally from 12 × depth to a
very large value, keep the onset selector and only loosen the cap, or introduce
a tuned visible-spacing multiplier immediately.

**Rationale:** The WP-05 full-set failure had two separable structural causes:
the comparison was hard-limited by a cell-scale cap that should not have been
applied directly to visible streak spacing, and the CL candidate's onset-only
selector made the pre-coarsening width almost depth-only. Splitting the scales
fixes the comparison-layer ceiling without discarding the shallow-water cell
constraint, while the unstable-band selector restores forcing sensitivity to the
CL candidate. The rerun improved the provisional BM-A full metrics materially
(`cl`: RMSE `117.0 → 112.5 m`, range coverage `0.175 → 0.356`; `scaling`: RMSE
`124.7 → 117.9 m`, range coverage `0.059 → 0.249`) even though the gate is
still negative.

**Reversible?** Yes. The observation-scale proxy and CL selector are both
explicit helpers with documented assumptions and can be replaced by a more
derived tracer/visibility model once one exists.

## 2026-03-22 — Bound the visible-spacing hierarchy to 0–3 merger events and rerun WP-05 on forcing-matched weather

**Decision:** Replace the observation-scale point proxy based on the full
uncapped merger cascade with an explicit `0–3` visible-merger hierarchy, expose
the observation multiplier and cell-cap settings as WP-05 rerun parameters, and
rerun the benchmark against the external ERA5 cache. Keep the mechanical cell
cap at `12 × depth` by default.

**Alternatives considered:** Keep the uncapped hierarchy, relax the global cell
cap to `300 × depth`, or treat a common visibility multiplier as the primary
fix before forcing-matched reruns.

**Rationale:** The matched reruns showed that the structural collapse was being
caused by the observation layer, not by the mechanical `12h` cap. With the
bounded visible hierarchy and forcing-matched inputs, both candidates pass the
attractor test on BM-A full and the CL candidate regains healthy range coverage
(`0.631`). Relaxing the global cell cap from `12h` to `300h` barely changes the
matched metrics, so there is no evidence that the hydro cap is the dominant
constraint in WP-05. Common visibility multipliers increase CL range coverage
further, but they worsen fit and still do not make the scaling candidate pass
the dynamic-range criterion.

**Reversible?** Yes. The matched rerun parameters are explicit in
`run_wp05_comparison(...)`, and the visible-spacing proxy can still be replaced
by a derived tracer/visibility model later.

## 2026-03-22 — Replace the constant coarsening clock with a width-squared merger schedule

**Decision:** Remove the fixed `h/u*` merger clock and replace it with a
sequential width-dependent schedule based on the nonlinear slow time
`T = l² t`. The implemented merger time is
`τ_merge(λ) = λ² / ((2π)² ν_T)`, recomputed after each merger as the width
doubles.

**Alternatives considered:** Keep the constant `h/u*` rule, replace it
one-for-one with `h / (u*² u_s0)^{1/3}`, or use `λ² / ν_T` without the
`(2π)²` factor.

**Rationale:** The constant merger clock was the direct cause of the absurd raw
coarsening counts in forcing-matched runs. A direct Langmuir-overturning
replacement would have made the counts even larger because, in the current
forcing path, `w_L > u*` for the relevant `La_t < 1` cases. The slow-time
scaling from Hayes & Phillips provides a first-principles reason that merger
times should grow as `1/l²`, i.e. as `λ²`. Using the existing `ν_T` closes that
timescale with already-audited forcing variables and makes later mergers
automatically slower by a factor of four per doubling.

**Reversible?** Yes. The prefactor is explicit in `src/hydro/coarsening.py`, and
the merger schedule can still be replaced by a more derived structural model if
field constraints become available.

## 2026-03-23 — Replace fixed Robin defaults with a forcing-derived closure and add onset-only CL validation

**Decision:** Remove the fixed CL-default Robin parameters from the operational
prediction path. The CL candidate now derives `γ_s` from wave steepness and
`γ_b` from the depth-dependent bed reach of the same wave field, and the WP-05
comparison harness can run the CL candidate in an onset-only diagnostic mode
that bypasses supercritical scale selection, coarsening, and observation-scale
hierarchy expansion.

**Alternatives considered:** Keep the historical constant Robin defaults,
segment evaluation by lake class to hide the onset mismatch, or retune the
coarsening/visibility layers before isolating the upstream width control.

**Rationale:** The remaining WP-05 blocker is upstream. Once the width-squared
coarsening schedule removed the merger-count explosion, the constant Robin
closure became the dominant source of onset-width rigidity. Replacing it with a
forcing-derived closure restores direct dependence on wave steepness and depth,
and the onset-only mode gives a clean way to test that raw instability-scale
response without confusing it with downstream widening layers. This keeps the
repair physically targeted and avoids hard-coding site-specific gates.

**Reversible?** Yes. The closure and onset-only mode are both explicit and
isolated. They can be refined or replaced once a stronger stratified onset
derivation is available for the production affine solver.

## 2026-03-23 — Treat the shallow-lake unstratified limit as a flagged assumption, not a WP-06 blocker

**Decision:** Proceed with the operational onset-width path under the explicit
assumption `S = 0` for the shallow benchmark lakes, while recording missing
stratification physics as a flagged uncertainty rather than as a hard stop for
WP-06.

**Alternatives considered:** Hold WP-06 until a full stratified affine onset
derivation exists, or add an ad hoc stratification correction without a clean
derivation.

**Rationale:** The considered lakes are shallow enough that the unstratified
limit is a defensible first operational assumption, especially for wind-mixed
events. Hayes (2020) makes the main risk explicit: if the water column is
actually stratified, neglecting `S` will bias onset toward larger `l_c` and
hence narrower predicted initial widths. That is a real uncertainty and should
be visible in the assumptions and risks, but it does not justify inventing an
untested correction or freezing the program at WP-05.

**Reversible?** Yes. Once a clean-room stratified onset derivation is available,
the `S = 0` assumption can be superseded without changing the surrounding
prediction/evaluation interfaces.

## 2026-03-23 — Do not promote a production candidate from the provisional WP-06 comparison

**Decision:** Keep the full WP-06 output package, record `cl` as the current
leading physics candidate, but do not promote any candidate into the production
pipeline from the present run.

**Alternatives considered:** Promote `cl` anyway because it beats `scaling` and
passes the attractor check, promote the best-fitting statistical baseline, or
treat stratification as the blocker and freeze WP-06.

**Rationale:** The current run had to use the provisional site-proxy forcing
path because `data/raw/era5_cache/` is empty, and the WP-06 brief requires the
winning physics candidate to beat the baselines on RMSE, dynamic range, and
tail coverage. On `BM-A_full`, `cl` leads the physics pair (`RMSE = 129.9 m`
vs `131.9 m`; `range coverage = 0.058` vs `0.028`; `tail coverage = 0.186` vs
`0.000`) and both physics candidates pass the attractor test, but `cl` still
trails the baselines on the required promotion metrics
(`baseline_depth_scaled`: `RMSE = 87.6 m`, `range coverage = 0.167`,
`tail coverage = 0.256`; `baseline_linear_wind`: `RMSE = 87.7 m`,
`range coverage = 0.172`, `tail coverage = 0.279`). The explicit shallow-lake
`S = 0` assumption remains visible but is not the blocker here.

**Reversible?** Yes. Populate the matched ERA5 cache and rerun WP-06; a matched
rerun can confirm or overturn the provisional physics ranking without changing
the public interfaces.

## 2026-03-23 — Use the filled ERA5/Open-Meteo cache as the default comparison forcing source, but do not infer that forcing alone caused the collapse

**Decision:** Adopt the filled `data/raw/era5_cache/` directory as the default
source for matched comparison reruns and treat the provisional site proxy as a
fallback only. Do not conclude that the previous range collapse was purely a
forcing-proxy artefact.

**Alternatives considered:** Continue debugging on the provisional proxy, wait
for a direct CDS export before rerunning, or treat the matched rerun as an
automatic promotion of the leading physics candidate.

**Rationale:** The matched rerun answers the immediate forcing question. It does
change the structural picture: both physics candidates pass the attractor test,
and full-set range coverage increases materially (`cl`: `0.058 → 0.137`;
`scaling`: `0.028 → 0.180`). But the matched results still fail the WP-06
promotion bar by a wide margin relative to the best baseline
(`baseline_depth_scaled`: `RMSE = 83.8 m`, `range coverage = 0.334`,
`tail coverage = 0.395`). On the matched run, `cl` remains the lower-RMSE
physics candidate, while `scaling` spans slightly more of the observed range,
so forcing matching reduces one uncertainty without resolving the physical
ranking or the baseline gap. That means the next inference burden shifts away
from weather matching and back onto the hydro/coarsening/visibility stack.

**Reversible?** Yes. The cache format and download path are explicit. If a
direct ERA5 source replaces the current Open-Meteo `models=era5` feed, the
comparison consumer can keep the same interface while the cache writer changes.

## 2026-03-23 — Prioritise the CL supercritical selector over the coarsening and visibility layers in the next repair pass

**Decision:** Treat the CL supercritical width selector in
`src/prediction/candidate_cl.py` as the primary next repair target. Do not
prioritise the width-squared merger schedule or the `0–3` visible-hierarchy
bound as the main source of the current matched range collapse.

**Alternatives considered:** Re-open the forcing proxy question, retune the
merger schedule first, or relax the visible hierarchy bound before isolating the
selector stage.

**Rationale:** The matched structural audit in
`outputs/structural_audit_matched/summary.json` shows that the raw onset width
already has material variance (`23.5–142.6 m`, `range coverage = 0.248`) and is
strongly responsive to forcing (`corr(L_inst, Ra) = -0.72`). The main collapse
appears immediately after the `upper_unstable_band` selector, which compresses
the comparison width to roughly one eighth of the onset width for almost every
case (`selected/raw ratio median = 0.125`), yielding selected pre-merger widths
of only `3.0–17.8 m` (`range coverage = 0.031`). By contrast, the merger clock
is active (median `2` events, max `3`), and the visible bound is not active in
the matched run because `visible_n_events` never falls below `n_events`.

**Reversible?** Yes. The selector is already an explicit helper with method
metadata in the output tables. It can be replaced or constrained without
changing the forcing, coarsening, or visibility interfaces.

## 2026-03-23 — Reverse the March 22 unstable-band selector and seed CL coarsening from the nonlinear onset width

**Decision:** Reverse the March 22 decision that selected the CL comparison
width from the upper edge of the nonlinear unstable band. The operational CL
candidate now passes the nonlinear onset minimum `l_cNL` directly into the
cell-width conversion and discrete coarsening path
(`L_cell,0 = 2π × depth / l_cNL`).

**Alternatives considered:** Keep the `upper_unstable_band` selector and adjust
only the cap, add a tuned visible-spacing multiplier on top of the old
selector, or leave the old selector in place and continue debugging farther
downstream.

**Rationale:** The matched audit showed that the onset width already has real
variance and forcing sensitivity (`23.5–142.6 m`,
`corr(L_inst, Ra) = -0.72`), while the `upper_unstable_band` selector shrank it
to `3.0–17.8 m` before coarsening. The repaired matched rerun confirms that the
selector was structurally harmful: on `BM-A_full`, the CL candidate improves
from `RMSE = 126.4 m`, `range coverage = 0.137`, `tail coverage = 0.116` to
`RMSE = 116.6 m`, `range coverage = 0.246`, `tail coverage = 0.093` without any
change to the forcing source, visible hierarchy, or baselines. That does not
clear the WP-06 bar, but it is the cleaner physics path because the onset
solver is now allowed to carry its own variance into the downstream layers.

**Reversible?** Yes. If a better justified supercritical comparison scale is
derived later, it can still be inserted explicitly and benchmarked against this
onset-seeded path.

## 2026-03-23 — Keep the `12 × depth` cap as a mechanical diagnostic safeguard, but do not treat it as a live WP-06 comparison lever

**Decision:** Keep the `12 × depth` limit on mechanically active cell width for
now, but do not treat relaxing or removing that cap as the next repair lever
for WP-06 candidate selection.

**Alternatives considered:** Remove the cap entirely, raise it globally for the
matched comparison path, or continue assuming that the cap is a dominant cause
of the current range collapse.

**Rationale:** The cap-sensitivity rerun shows that the current public
comparison is not being controlled by this limit. Raising
`max_cell_aspect_ratio` from `12` to `300` changes all `67/67` CL mechanical
cell-width diagnostics and all `67/67` cap-binding flags, but changes `0/67`
comparison predictions and `0` WP-06 metric rows. Under the bounded visible
hierarchy, the cap still serves as a physically explicit warning on the
mechanical cell state, but it is not the bottleneck on the output being scored
against the observations.

**Reversible?** Yes. If later work introduces a direct mechanical-width target
or a different observation model that depends on the capped cell state, the cap
policy can be revisited with that interface in hand.

## 2026-03-23 — Refine the AR-026 diagnosis: the matched WP-06 path already uses dynamic weather-history lifetimes, so repair reset-event timing but do not treat the provisional `1800 s` proxy as the matched blocker

**Decision:** Keep the fixed `1800 s` assumption confined to the provisional
fallback path, but do not frame it as the active blocker on the matched WP-06
comparison. On the matched path, repair the organisation-time logic so the
coarsening clock is counted from the latest actual reset event rather than from
the latest sample touched by an older event still inside the lookback window.

**Alternatives considered:** Remove AR-026 globally before checking the matched
path, extend the matched lookback window first, or bypass the reset logic and
assume the full matched history is always usable.

**Rationale:** The matched pipeline in `src/prediction/pipeline.py` was already
using history-derived organisation times; the real bug was that prefix-based
reset handling could re-stamp an old event onto later steady samples. The
repair changes `9/67` matched lifetimes by up to `4.5 h` and improves the
scaling candidate slightly, so it is worth keeping. But it does not move the
CL candidate at all, which means the matched WP-06 failure is not being driven
primarily by the provisional `1800 s` proxy.

**Reversible?** Yes. If later work adopts a different reset detector or a
longer weather-history window, the same matched interface can consume that
revised lifetime estimate.

## 2026-03-23 — Prioritise the CL coarsening-timescale / `nu_T` closure over further organisation-time or visibility edits

**Decision:** Treat the operational CL coarsening timescale as the next repair
target. Do not spend the next pass primarily on further organisation-time
plumbing or on the visible-hierarchy handoff for the CL candidate.

**Alternatives considered:** Loosen the visible `0–3` merger bound first,
extend the matched lookback window before touching the CL coarsening law, or
continue iterating on reset-threshold details.

**Rationale:** After the lifetime fix, the matched CL run still has
`64/67` zero-merger cases and only `3/67` one-merger cases. The first-merger
timescale is `4.1–123.7 h` with a `24.0 h` median, while repaired matched
lifetimes still top out at `5.5 h`; only `3/67` CL cases even have enough time
for one merger. The visibility handoff is already faithful to the achieved CL
mechanical state (`visible_n_events = n_events` in every CL case), so there is
no evidence that AR-027 is currently crushing CL range. The structural miss is
upstream, in how slowly the CL path coarsens under the present `nu_T` closure.

**Reversible?** Yes. If a revised coarsening closure yields materially larger
mechanical event counts, the visibility bound and lookback window can be
re-audited against that new regime.

## 2026-03-24 — Separate the merger-clock diffusivity from the forcing-layer vertical `nu_T` and use a lateral mixing-length closure as the working coarsening law

**Decision:** Keep the forcing-layer depth-mean `nu_T` for the `Ra` mapping,
but stop reusing it as the merger-clock diffusivity. The working operational
coarsening closure is now
`A_H = max(nu_T_vertical, u*_water × depth)`, which is effectively
`A_H = u*_water × depth` on the matched WP-06 cases.

**Alternatives considered:** Continue using the forcing-layer vertical `nu_T`
for the merger clock, scale that `nu_T` by an ad hoc constant chosen to rescue
metrics, or leave the clock frozen until a more elaborate anisotropic
turbulence closure is derived.

**Rationale:** The closure audit showed the previous merger clock was being fed
with a distinctly vertical mixing scale:
`forcing_nu_T_m2_s = 1.5e-4–7.6e-3 m²/s` with a `1.05e-3 m²/s` median. That
produced a CL first-merger timescale of `4.1–123.7 h` with a `24.0 h` median
in the prior matched run, leaving only `3/67` cases with enough organisation
time for one merger. Switching to the lateral mixing-length closure raises the
clock diffusivity to `2.2e-3–1.11e-1 m²/s` (median `1.54e-2 m²/s`), restores
the expected `O(h/u*)` behaviour, and cuts the CL first-merger median to
`1.64 h`, with `52/67` cases now able to complete at least one merger.

This is a physics-grounded repair, not a production-candidate win. The matched
rerun confirms that the `12 × depth` cap is still not the scored WP-06 lever:
raising it to `300 × depth` changes all `67` CL cap-binding flags but changes
`0` comparison predictions and `0` metric rows. The repaired CL path still
fails `BM-A_full` badly (`RMSE = 159.9 m` versus `83.8 m` for the best
baseline), so this closure should be treated as the new working audit baseline,
not as a promotable final answer.

**Reversible?** Yes. If later work derives a stronger shallow-lake lateral
closure or a different finite-amplitude widening law, the coarsening clock can
be re-parameterised without changing the forcing-to-`Ra` mapping.
