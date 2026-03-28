# Assumptions Register

**Status:** Living document. Every assumption introduced by any work package must be
recorded here. Entries are never deleted — only superseded (mark old entry SUPERSEDED
and add new one).

**Format per entry:**
- **ID:** AR-NNN
- **Assumption:** Statement of what is assumed
- **Source:** Literature reference, physical argument, or calibration
- **Used in:** Which module(s) depend on this assumption
- **Testable?** Whether the assumption can be falsified with available data
- **Risk:** What goes wrong if the assumption is violated
- **Status:** Active / Superseded / Under review

---

## Boundary Conditions

### AR-001 — Robin BC: free-surface parameter γ_s

**Assumption:** The Robin boundary condition parameter at the free surface is
γ_s ≈ 0.06.

**Source:** Cox & Leibovich (1993), physical argument based on wave slope and surface
momentum exchange. Hayes & Phillips (2016, 2017) adopt these values without revision
for the s = 1 (shallow coastal) case.

**Used in:** `src/hydro/robin_bc.py`, `src/hydro/linear_solver.py`,
`src/hydro/nonlinear_solver.py`

**Testable?** Partially. The dependence of l_cNL on γ_s is smooth and can be
bracketed by sensitivity runs. Phillips & Dai (2014) derive γ as a function of
wavenumber and wave parameters; the constant approximation is testable by comparison
with their full expression.

**Risk:** γ_s sets the free-surface restoring force for the perturbed flow. If the
true γ_s is substantially larger (e.g. 0.2–0.5 as suggested by Phillips & Dai 2014
for longer waves), the onset wavenumber increases and predicted cell widths decrease.
This would narrow the predicted spacing distribution.

**Status:** Superseded by AR-028.

---

### AR-002 — Robin BC: bottom parameter γ_b

**Assumption:** The Robin boundary condition parameter at the bottom boundary is
γ_b ≈ 0.28.

**Source:** Cox & Leibovich (1993), physical argument. Represents partial stress
reflection from a rough bottom in the presence of turbulent bottom boundary layer.

**Used in:** `src/hydro/robin_bc.py`, `src/hydro/linear_solver.py`,
`src/hydro/nonlinear_solver.py`

**Testable?** Weakly. The γ_b value affects l_cNL through its interaction with the
bottom boundary layer. Sensitivity analysis across γ_b ∈ [0, 1] should be run.

**Risk:** A rigid no-slip bottom (Neumann, γ_b → ∞) would produce l_c → 0 (infinite
cell width) without the surface Robin condition providing a lower bound. The physical
bottom condition in a shallow lake is partially rough; γ_b ≈ 0.28 is a reasonable
mid-range estimate but uncertain by factor ~2.

**Status:** Superseded by AR-028.

---

## Instability Physics

### AR-003 — CL2 instability assumed dominant

**Assumption:** Langmuir circulation in the study lakes arises through the CL2
(Craik–Leibovich type 2) instability, driven by the interaction of wind-driven shear
with the depth-averaged Stokes drift gradient D'(z). The CL1 instability (requiring
cross-wind drift variation) is neglected.

**Source:** Craik & Leibovich (1976): cross-wind drift variations phase-average to
zero in a random wave spectrum, leaving CL2 as the generic mechanism. Confirmed for
shallow coastal layers by Phillips & Dai (2014).

**Used in:** All hydro module solvers.

**Testable?** The CL1 contribution requires coherent cross-wind wave structure
(swell from a single direction). Its presence could be inferred from strong directional
wave data, not available here.

**Risk:** Low for random-wave conditions. In cases of strong, directional swell the
CL1 mechanism may augment CL2 and produce narrower, more intense cells.

**Status:** Active.

---

### AR-004 — Shear index s = 1 (shallow coastal)

**Assumption:** The shear index s = 1, appropriate for shallow coastal layers where
wave amplitude/depth is not negligible, as opposed to s = 2 for the open ocean.

**Source:** Phillips (1998a), Phillips & Dai (2014), Hayes & Phillips (2017) §2.1.
The s = 1 limit means the LC timescale scales as O(ε^{3/2}), slightly faster than
the open-ocean O(ε²).

**Used in:** `src/hydro/profiles.py`, Rayleigh number definition.

**Testable?** Indirectly: if the deep-water s = 2 equations were used instead,
the predicted onset threshold Ra would differ. No direct observational test available.

**Risk:** Low for h < λ_p/2 (shallow water wave condition). If the peak wave
wavelength is short enough that the site is effectively deep (λ_p << 2h), s = 2
would be more appropriate.

**Status:** Active.

---

## Coarsening

### AR-005 — Y-junction merger as coarsening mechanism

**Assumption:** The dominant mechanism by which LC cell width increases with time
is the Y-junction merger event described by Thorpe (2004): adjacent cell pairs of
slightly different width merge at a Y-junction defect, doubling the local spacing.
Each merger event approximately doubles the cell width.

**Source:** Thorpe (2004), field and laboratory observations.

**Used in:** `src/hydro/coarsening.py`

**Testable?** Yes, in principle: the coarsened width after n events is 2^n × L_inst.
If the time between mergers is O(h/u*), then wind episode duration constrains n.
This is testable against observations if LC age is approximately known.

**Risk:** If the dominant coarsening mechanism is not Y-junction merger but rather
direct widening through turbulent diffusion, the doubling-per-event relationship
fails. The model would then overestimate the discreteness of the width distribution.

**Status:** Active.

---

### AR-006 — Aspect ratio cap at 12

**Assumption:** LC cell spacing is capped at 12 × depth. Cells wider than this are
not observed in shallow water.

**Source:** Marmorino et al. (2005), from infrared and visible imagery of shallow
coastal LC.

**Used in:** `src/hydro/coarsening.py`

**Testable?** Weakly: no observation in the current dataset shows spacing/depth > 12.
But the dataset may not sample extreme forcing regimes.

**Risk:** The cap must emit a warning if triggered rather than silently truncating
the spacing. A silent cap is a saturation mechanism. See AGENTS.md §2.

**Status:** Active.

---

## Observation Model

### AR-007 — Visible windrow spacing ≈ mature cell spacing (with uncertainty)

**Assumption:** The satellite-visible windrow spacing recorded in the observations
is approximately equal to the LC cell width after coarsening (L_cell), subject to
a tracer-accumulation factor and visibility threshold.

**Source:** Observation model document (`docs/observation_model.md`), §5.

**Used in:** `src/prediction/visibility.py`, evaluation metrics.

**Testable?** Yes: if the visibility factor is set to 1 (direct comparison), residuals
should be uncorrelated with tracer indicators. If they are correlated, the tracer
model needs adjustment.

**Risk:** If the visible spacing consistently differs from cell spacing by a fixed
factor, this will appear as a systematic bias in the metrics. The bias should be
diagnosed as a model deficiency, not tuned away.

**Status:** Active.

---

### AR-008 — Non-detections treated as censored, not subcritical

**Assumption:** Observation records with no spacing measurement are censored
(LC may or may not be present; the satellite image was uninformative) rather than
confirmed subcritical events.

**Source:** Physical reasoning: see `docs/observation_model.md` §4.

**Used in:** `src/evaluation/metrics.py`, benchmark dataset construction.

**Testable?** Weakly: if forcing data exists for non-detection records, the Ra
distribution of non-detections can be compared to that of detected events.

**Risk:** If non-detections are actually dominated by subcritical forcing, treating
them as censored inflates the precision of near-critical predictions. The sensitivity
to this assumption should be tested during evaluation.

**Status:** Active.

---

## Forcing Estimation

### AR-009 — Closed-basin surface current (no 3% rule)

**Assumption:** The surface current in a closed basin (no net throughflow) is
estimated by solving the depth-integrated momentum equation with wind stress at the
surface and a closed-basin constraint (zero net transport), rather than applying
the oceanic 3% heuristic (U_surface = 0.03 × U10).

**Source:** Physical reasoning for enclosed lakes. The 3% rule applies to open-ocean
Ekman drift and is inappropriate here.

**Used in:** `src/forcing/currents.py`

**Testable?** Yes: the closed-basin solution must satisfy ∫_{-h}^{0} u(z) dz = 0
by construction (machine precision test).

**Risk:** Low for the physics. The uncertainty in eddy viscosity propagates into the
surface current estimate and thence into Ra.

**Status:** Active.

---

### AR-010 — Parabolic eddy viscosity profile

**Assumption:** The eddy viscosity profile is parabolic:
ν_T(z) = κ u* (z + h)(1 − (z + h)/h), where κ = 0.41 is the von Kármán constant.

**Source:** Standard log-layer turbulence model (Prandtl mixing length). This is a
leading-order approximation; KPP-type profiles are more accurate but require
additional parameters.

**Used in:** `src/forcing/eddy_viscosity.py`

**Testable?** The depth-averaged ν_T should scale as u* × h to within a factor of 2,
consistent with turbulence scaling. The functional form is not directly testable
without in-situ turbulence measurements.

**Risk:** The parabolic profile underestimates ν_T near the surface in the presence
of strong Langmuir circulation (which enhances vertical mixing). This creates a
potential circular dependency: LC modifies ν_T, which modifies Ra, which modifies LC.
This feedback is neglected at present and flagged for future work.

**Status:** Active.

---

## Scaling Parameters

### AR-011 — Coarsening follows the nonlinear slow-time scaling τ ~ λ² / ((2π)² ν_T)

**Assumption:** The time for one Y-junction merger event is controlled by the
current horizontal cell width λ and the effective turbulent eddy viscosity ν_T,
not by a fixed vertical turnover time. The implemented merger schedule uses

τ_merge(λ) = λ² / ((2π)² ν_T)

and updates τ after each merger as the cell width doubles.

**Source:** Hayes & Phillips (2017), Eq. (27), introduces the nonlinear slow
time `T = l² t`; with `λ = 2πh/l`, this implies a physical timescale that grows
as `λ²`. The use of `ν_T` supplies the dimensional diffusion scale already
present in the forcing layer.

**Used in:** `src/hydro/coarsening.py`, `src/prediction/candidate_cl.py`,
`src/prediction/candidate_scaling.py`

**Testable?** Yes. If matched weather histories are available, the resulting
merger counts can be compared against the ratio of observed to initial spacing.
The key falsification check is whether the law keeps raw merger counts in the
same order as the `0–3` visible hierarchy expected from the imagery.

**Risk:** If ν_T is too small, the law will still over-predict coarsening; if
ν_T is too large, coarsening will stall too early. The `1/(2π)²` factor is a
first-principles mapping from the slow-time scaling, but it is still a closure
assumption rather than a direct field measurement.

**Status:** Active.

---

### AR-012 — Fetch-limited wave spectrum

**Assumption:** The surface wave spectrum is fully determined by U10 and fetch, using
a fetch-limited empirical spectrum (JONSWAP or equivalent). Swell contributions from
outside the basin are neglected.

**Source:** Standard practice for enclosed lakes of limited fetch (< 100 km).

**Used in:** `src/forcing/waves.py`

**Testable?** Yes: wave period and height predictions from fetch-limited formulae can
be compared to any available wave buoy data for the site.

**Risk:** If swell penetrates from outside the basin (unlikely for Lough Neagh but
possible for other sites), the Stokes drift profile will be underestimated, and Ra
will be underestimated accordingly.

**Status:** Active.

---

---

## Growth Conditions Diagnostics

### AR-013 — Three-pathway growth enhancement model

**Assumption:** The effect of LC on cyanobacterial growth conditions is decomposed
into three independent pathways: light exposure, nutrient availability, and temperature
distribution. Each pathway is computed as an enhancement triple (value_LC,
value_static, ratio). A composite favourability index is a weighted combination of
the three ratios.

**Source:** Governing plan (WP-05, `lc_enhancement.py` specification). Pathway
decomposition is justified by the mechanistic independence of the three processes:
photosynthesis is light-limited, nutrient flux is diffusion-limited, and thermal
optimum is a separate growth-rate function.

**Used in:** `src/prediction/lc_enhancement.py`

**Testable?** Each pathway individually: light enhancement is verifiable from
published photic-zone depth estimates and surface irradiance measurements. Nutrient
enhancement is order-of-magnitude verifiable from eddy diffusivity ratios. Temperature
enhancement is verifiable if surface and bottom temperature data are available.

**Risk:** The pathways are not strictly independent. LC-enhanced nutrient flux affects
growth rate, which affects buoyancy, which affects how cells respond to LC. These
feedbacks are neglected in this first-order model. Stated clearly in the composite
index output.

**Status:** Active.

---

### AR-014 — Static case: buoyant cells pool at surface (photoinhibited)

**Assumption:** In the absence of LC, positively buoyant cyanobacteria pool at the
water surface where they experience surface irradiance reduced by photoinhibition
factor f_photoinhibition ≈ 0.3–0.7.

**Source:** Literature consensus on cyanobacterial surface pooling (see
literature_review.md §4.1). f_photoinhibition range from published Microcystis
irradiance-response curves.

**Used in:** `src/prediction/lc_enhancement.py` (`light_enhancement` function)

**Testable?** The photoinhibition factor is observable in growth-irradiance curves.
The surface-pooling assumption is testable if near-surface chlorophyll profiles are
available for calm conditions.

**Risk:** If buoyant cells are dispersed by wind mixing even in the absence of LC
(e.g., by surface Langmuir-scale waves below the LC threshold), the static irradiance
is not the photoinhibited surface value. This would reduce the enhancement ratio.

**Status:** Active.

---

### AR-015 — LC downwelling velocity as key coupling variable

**Assumption:** The characteristic downwelling velocity w_d is the primary coupling
variable between the hydrodynamic prediction and the growth conditions diagnostics.
w_d determines: (1) whether cells circulate (w_d > w_rise), (2) the circulation
period T = 2h/w_d, (3) the effective eddy viscosity ν_T_lc.

**Source:** Governing plan (WP-05). Derived from the nonlinear eigenfunction structure:
w_d is proportional to u* × (nonlinear eigenfunction amplitude) × l_cNL.

**Used in:** `src/prediction/lc_enhancement.py`, `src/prediction/visibility.py`

**Testable?** w_d can be estimated from field measurements of LC-induced tracer
displacement rates. Published values for shallow lakes: w_d ~ 1–10 mm/s (O(u*)).

**Risk:** The nonlinear eigenfunction amplitude requires the full nonlinear solver
(WP-03 code). Before the solver is available, w_d is estimated as w_d ~ α × u*
with α ~ 0.1–0.5 (order-of-magnitude). Using the linear eigenfunction amplitude
(FL-001) underestimates w_d by a factor related to the k=0 base-flow correction
(Hayes & Phillips 2017, Figure 6a, curve i).

**Status:** Active.

---

### AR-016 — Residual κ bridge outside the derived affine-profile family

**Assumption:** For nonlinear runs where either `D'` or `U'` is higher than
affine in `z`, the solver still falls back to a benchmark-style κ bridge based
on the profile endpoint contrasts ΔD = D'(0) - D'(-1) and ΔU = U'(0) - U'(-1).
This fallback is no longer used for the supported WP-04b verification family:
all affine profiles, including `D', U' ∈ {1, 1 + z}`, now use the explicit
mean-flow weakly nonlinear correction derived from the H&P solvability system.

**Source:** H&P (2017) §7.2 benchmark κ values for the legacy non-affine
fallback; affine weakly nonlinear solvability reduction implemented in
`src/hydro/nonlinear_solver.py`.

**Used in:** `src/hydro/nonlinear_solver.py`

**Testable?** Partially. The affine derivation is directly falsified by the
literature verification tests in `src/hydro/tests/test_literature_values.py`.
The residual non-affine bridge remains an explicit fallback until the full
higher-order polynomial solvability is implemented.

**Risk:** The supported verification family is now algebraic rather than
benchmark-fitted, but predictions for profiles beyond the affine family are
still not physically validated by a complete nonlinear derivation.

**Status:** Active.

---

### AR-017 — La-based width closure from thickness and pitch

**Assumption:** The non-CL scaling-law candidate closes the cell-pair width as
`estimated_cell_width = 2 × downwelling_thickness × pitch`, where
`downwelling_thickness ~ h × La_t^{1/2}` and `pitch ~ La_t^{1/6}`. The LES
literature gives the proportional scalings for thickness, velocity and pitch
but not a unique shallow-water width formula, so this closure keeps the width
explicitly tied to those published ingredients.

**Source:** McWilliams, Sullivan & Moeng (1997) turbulent Langmuir number
framework; scaling summary in `docs/literature_review.md` §1.3.

**Used in:** `src/hydro/scaling_laws.py`

**Testable?** Yes. The resulting width must vary monotonically and avoid
saturation across plausible `La_t`, `depth` and `u*` ranges. It can also be
compared against the CL-based width after the WP-04 Decision Gate.

**Risk:** If the real shallow-water width depends on a different combination of
the LES geometry scalings, this candidate will mis-rank the scaling-law
alternative even if the underlying `La_t` trends are correct.

**Status:** Active.

---

### AR-018 — Rapid wind-increase disruption trigger

**Assumption:** The coarsening reset logic flags a rapid wind-speed increase when
the latest `U10` is at least 50% and 3 m/s above the recent median over the
lookback window. Direction changes >45° and drops below 2 m/s use the literature
thresholds directly; the rapid-increase trigger is a practical detection rule
for the otherwise qualitative “rapid increase” mechanism.

**Source:** Thorpe (2004) and Marmorino et al. (2005) disruption mechanisms,
summarised in `docs/literature_review.md` §2.2.

**Used in:** `src/hydro/coarsening.py`

**Testable?** Indirectly. If wind histories with known pattern resets become
available, the trigger can be checked against observed decay-and-reformation
episodes.

**Risk:** If real structures tolerate faster wind increases without resetting,
the current rule will reset the coarsening clock too aggressively.

**Status:** Active.

---

### AR-019 — Visibility diagnostic thresholds

**Assumption:** A Langmuir pattern is treated as satellite-visible only when
surface convergence is at least O(10⁻⁴ m/s), the 10 m wind speed is not so high
that surface roughness obscures streaks (`U10 \lesssim 10 m/s`), and the pattern
lifetime exceeds the tracer-accumulation time. The resulting visibility
confidence is the product of the convergence, wind and persistence scores.

**Source:** WP-05b specification structure in `docs/clean_room_rebuild_brief.md`;
AR-007 states that visibility is diagnostic, not a spacing correction.

**Used in:** `src/prediction/visibility.py`

**Testable?** Yes, in principle. If image metadata and wind conditions are
available, detectability failures can be compared against these diagnosed
limiting factors.

**Risk:** The thresholds are first-order heuristics. If the true visibility
transition occurs at different convergence or wind levels, the confidence score
will be miscalibrated even if the spacing prediction itself is correct.

**Status:** Active.

---

### AR-020 — Temperature enhancement uses inverse distance to optimum

**Assumption:** Thermal favourability in the LC enhancement calculator is
represented by inverse distance to the growth optimum,
`proximity = 1 / (1 + |T - T_opt|)`, and the enhancement ratio is the LC/static
ratio of that proximity measure. This avoids introducing an additional free
thermal-width parameter while preserving the sign of “closer to optimum” versus
“further from optimum”.

**Source:** WP-05c specification in `docs/clean_room_rebuild_brief.md`, which
defines enhancement through proximity to the optimum rather than absolute
temperature change.

**Used in:** `src/prediction/lc_enhancement.py`

**Testable?** Yes. If surface and mixed-layer temperatures are observed together
with growth-rate response curves, the proxy can be compared against a fuller
thermal performance function.

**Risk:** Real cyanobacterial growth responds nonlinearly and asymmetrically to
temperature, so this symmetric inverse-distance proxy may understate penalties
far above the optimum or benefits near the peak.

**Status:** Active.

---

### AR-021 — CL candidate reduces forcing profiles to the affine H&P family

**Assumption:** The WP-05 CL candidate projects the resolved forcing drift
profile onto the supported affine solver family by normalising `D'(z)` with its
surface magnitude and matching the bottom and surface endpoints on the
nondimensional interval `z ∈ [-1, 0]`. Because the forcing layer does not yet
provide a resolved mean-shear profile compatible with the CL solver, the
candidate uses the surface-intensified affine proxy `U'(z) = 1 + z` by default.

**Source:** H&P (2017) verification family for affine profiles; forcing drift
profile from `src/forcing/__init__.py`; unresolved shallow-lake shear profile
noted in AGENTS.md §13.

**Used in:** `src/prediction/candidate_cl.py`

**Testable?** Yes. The projection is falsified if a later resolved shallow-lake
shear/drift profile module materially changes the predicted spacing or κ values
for the same forcing histories.

**Risk:** The default `U'(z) = 1 + z` proxy may bias the CL candidate toward the
surface-intensified verification family. If the real lake shear profile is more
uniform or more strongly surface-confined, the predicted spacing may shift by a
nontrivial factor.

**Status:** Active.

---

### AR-022 — Instability wavenumber maps to cell-pair width as `2πh / l_c`

**Assumption:** The nondimensional CL instability wavenumber is converted to the
pre-coarsening cell-pair width through
`L_cell = 2π × depth / l_c`.

**Source:** H&P (2017) spanwise-wavenumber definition; AGENTS.md §4 and §5.3 on
explicit scale conversion.

**Used in:** `src/prediction/common.py`

**Testable?** Yes. The resulting widths can be checked against literature
aspect-ratio ranges and against any systematic factor-of-two bias in the
prediction residuals.

**Risk:** If the instability mode corresponds to half a cell pair rather than a
full pair, the visible-spacing prediction would be biased by a factor of two.

**Status:** Active.

---

### AR-023 — Visibility convergence estimated from continuity scaling

**Assumption:** The visibility diagnostic uses the continuity-based surface
convergence proxy
`convergence_velocity = |U_surface| × depth / predicted_spacing`,
so stronger shear and narrower cells give faster tracer accumulation.

**Source:** Mass-continuity scaling between horizontal convergence and vertical
circulation; WP-05b visibility diagnostic requires a convergence-velocity input.

**Used in:** `src/prediction/common.py`

**Testable?** Yes. If observed detectability metadata or resolved LC velocity
fields become available, this proxy can be compared directly to diagnosed
surface-convergence rates.

**Risk:** The proxy is order-of-magnitude only. If actual surface convergence is
systematically weaker or stronger than this continuity estimate, the visibility
confidence will be miscalibrated.

**Status:** Active.

---

### AR-024 — Default environmental context for single-case prediction runs

**Assumption:** When a case is analysed without explicit environmental inputs,
the prediction layer uses the following defaults:
`surface_irradiance = 400 W/m²`,
`K_d = 0.5 1/m`,
`photoinhibition_factor = 0.5`,
`surface_temperature = 20°C`,
`bottom_temperature = 16°C`,
`T_optimum = 25°C`,
`nutrient_gradient = 10 µg/L per m`,
`tracer_accumulation_time = 1800 s`,
`static_nu_T = 1e-4 m²/s`.

**Source:** WP-05b/c specification structure and the static-column description
in `docs/clean_room_rebuild_brief.md`; the chosen `static_nu_T` sits at the top
of the stated weak-mixing range `10⁻⁵–10⁻⁴ m²/s`.

**Used in:** `src/prediction/common.py`, `src/prediction/pipeline.py`

**Testable?** Partially. The defaults should be superseded by case-specific
observations whenever irradiance, attenuation, temperature, nutrient, or tracer
timescale data are available.

**Risk:** These defaults make the diagnostic outputs reproducible but not
site-specific. If the true environmental state differs substantially, the
enhancement and visibility diagnostics will be biased even when the spacing
prediction is reasonable.

**Status:** Active.

---

### AR-025 — WP-05 full-set comparison uses a provisional site/date forcing proxy

**Assumption:** Because the benchmark spec's ERA5 cache is not present in this
workspace, the full observation-set WP-05 comparison classifies each observation
into one of four site buckets (`taihu`, `erie`, `prairie`, `neagh`), assigns an
explicit representative depth/fetch pair, assumes a 10:30 local overpass, and
uses a clipped seasonal cosine proxy for representative wind speed:
`U10 = mean_site + amplitude_site × cos(2π (doy - 30) / 365.25)`,
clipped to `[2.5, 10.5] m/s`.

**Source:** `data/benchmarks/benchmark_spec.md` requires site-specific depth and
fetch and notes morning overpass times near 10:30 local. Taihu (`h ≈ 2 m`,
`fetch ≈ 30 km`) and Lough Neagh (`h ≈ 9 m`, `fetch ≈ 15 km`) come from the
benchmark spec; Erie-region and prairie-lake values are explicit provisional
engineering placeholders until matched site metadata and ERA5 forcing are added.

**Used in:** `src/evaluation/comparison.py`

**Testable?** Yes. This assumption should be superseded immediately once the
ERA5 cache and site-matched depth/fetch metadata exist, at which point the proxy
run can be compared directly to the forcing-matched run.

**Risk:** All full-set comparison metrics produced under this proxy are
provisional. A forcing-matched rerun may change both the candidate ranking and
the structural diagnostics.

**Status:** Uncertain.

---

### AR-026 — Provisional WP-05 comparison fixes organisation time at 1800 s and scores no-LC as 0 m

**Assumption:** The provisional full-set comparison evaluates both candidates at
a fixed organisation time `pattern_lifetime = 1800 s`, and when a candidate
returns `predicted_spacing_m = NaN` (no coherent LC), the comparison layer keeps
the raw NaN output but scores the case as `comparison_spacing_m = 0 m`.

**Source:** `docs/observation_model.md` notes overpass-timing uncertainty of
order tens of minutes to hours; `1800 s` is the same order as the default tracer
accumulation time in WP-05b/c. The `0 m` comparison score is an explicit penalty
for a "no coherent pattern" prediction against a detected spacing observation.

**Used in:** `src/evaluation/comparison.py`

**Testable?** Yes. With matched weather histories, the fixed lifetime can be
replaced by an explicit organisation-time estimate and the no-LC scoring choice
can be checked against any future censored/non-detection records.

**Risk:** If the true organisation times are systematically longer or shorter
than 1800 s, the current comparison will respectively understate or overstate
coarsening. The `0 m` scoring rule will also increase error metrics for any case
where the model is subcritical but the observation is a true detection.

**Status:** Uncertain.

---

### AR-027 — Observation-scale spacing uses a bounded visible merger hierarchy plus a loose safety cap

**Assumption:** The 12 × depth aspect-ratio limit applies to the mechanically
active LC cell width, but the satellite-visible streak spacing should be
compared against an explicitly bounded observation-scale hierarchy rather than
against the full uncapped `2^n` merger cascade. The prediction layer therefore
reports a plausible visible-spacing range corresponding to `0–3` merger events,
uses the upper edge of that range as the point comparison proxy, and applies
only a loose safety cap of `300 × depth` at the observation scale.

**Source:** `docs/observation_model.md` §2, §5, and §6 distinguish mature cell
width from expected visible spacing and state that downstream comparisons should
carry coarsening ambiguity explicitly as a `0–3` merger-event range. The
`300 × depth` limit remains a numerical safeguard, not a physical cell-scale
constraint.

**Used in:** `src/prediction/common.py`, `src/prediction/candidate_cl.py`,
`src/prediction/candidate_scaling.py`, `src/evaluation/comparison.py`

**Testable?** Yes. The forcing-matched reruns can compare the bounded visible
hierarchy against alternative point rules and common visible-spacing
multipliers. If a derived tracer model supersedes the hierarchy proxy, this
assumption should be replaced.

**Risk:** This fixes the observation-layer attractor failure without curing the
remaining fit problem. The matched reruns show that bounding the visible
hierarchy materially improves CL range coverage, but the mechanical coarsening
clock still produces implausibly large raw event counts under some matched wind
histories. The visibility proxy should therefore remain explicit and uncertain
until the forcing/coarsening timescale issue is resolved.

**Status:** Uncertain.

---

### AR-028 — Robin BC closure is derived from wave steepness and bottom orbital reach

**Assumption:** The operational CL candidate no longer uses fixed Robin
parameters. Instead, the surface Robin parameter is derived from the resolved
wave steepness,
`γ_s = a k_p = (H_s / 2) k_p`,
and the bottom Robin parameter is derived from the bed-reaching fraction of that
wave field,
`γ_b = γ_s / sinh(k_p h)`.
The raw values are preserved without clipping; unusually large shallow-water
bottom coupling is exposed by warnings rather than saturation.

**Source:** Cox & Leibovich (1993) and Hayes & Phillips (2017) require Robin
stress terms to keep onset at finite `l_c > 0`. The clean-room closure maps
those terms onto the resolved wave field using linear-wave steepness and the
standard vertical decay of orbital motion with depth.

**Used in:** `src/hydro/robin_bc.py`, `src/prediction/candidate_cl.py`,
`src/evaluation/comparison.py`

**Testable?** Yes. The closure should show explicit sensitivity to `U10`,
`k_p`, and `h`, and the resulting onset-width distribution can be checked
without coarsening via the onset-only comparison mode.

**Risk:** If the bottom-coupling factor is physically incomplete, the shallow
sites will still carry a width bias. This closure is mechanistic and auditable,
but it is still only a first-pass mapping from forcing to Robin stress.

**Status:** Uncertain.

---

### AR-029 — WP-05 onset-width validation can bypass coarsening and visibility

**Assumption:** For diagnostic reruns only, the CL candidate may be evaluated in
an onset-only mode that reports the raw instability width
`L_inst = 2πh / l_cNL`
with supercritical selector, coarsening, and visibility-hierarchy expansion all
disabled.

**Source:** WP-05 repair work showed that the remaining blocker sits upstream in
the onset-width closure. A diagnostic mode is needed to test whether the raw
instability scale regains dynamic range before any downstream widening layer is
re-engaged.

**Used in:** `src/prediction/candidate_cl.py`, `src/prediction/pipeline.py`,
`src/evaluation/comparison.py`

**Testable?** Yes. The onset-only rerun should move the raw width distribution
with forcing and depth while remaining invariant to pattern lifetime.

**Risk:** This is not a production observation model. If used outside diagnostic
comparison, it would understate mature visible spacing by construction.

**Status:** Uncertain.

---

### AR-030 — Operational shallow-lake onset runs assume negligible stratification (`S = 0`)

**Assumption:** For the present study lakes, the operational onset-width path is
run in the unstratified limit `S = 0`. This is treated as a practical
shallow-lake approximation rather than as a claim that the water column is never
thermally stratified.

**Source:** The benchmark lakes are shallow enough that well-mixed conditions
are often plausible over the observed wind events, while the current clean-room
production solver is derived and verified only for the unstratified affine
family. Hayes (2020) shows that stratification would act mainly as an
onset-side stabiliser by lowering `l_c` and raising `R_c`.

**Used in:** `src/hydro/linear_solver.py`, `src/hydro/nonlinear_solver.py`,
`src/prediction/candidate_cl.py`

**Testable?** Partially. If independent temperature profiles or mixed-layer
depth estimates become available, the onset-width residuals can be checked for a
systematic warm-season or weak-wind widening bias consistent with missing
stratification.

**Risk:** If a site is materially stratified, the current `S = 0` assumption
will tend to over-predict `l_c` and therefore under-predict onset width. In
practice that means the model may produce cells that are too narrow and may
understate the critical forcing threshold during stratified episodes.

**Status:** Uncertain.

---

### AR-031 — WP-06 production promotion requires forcing-matched comparison

**Assumption:** A WP-06 run using the provisional site-proxy forcing path may
generate the full comparison output package and a provisional physics ranking,
but it may not promote any candidate into production. Production promotion
requires the ERA5 cache-matched path so the public weather-to-prediction
pipeline and the comparison numbers are evaluated on the same forcing history.

**Source:** `docs/clean_room_rebuild_brief.md` WP-06 requires production cleanup
only after the formal comparison is locked, and
`docs/wp06_handover_2026-03-23.md` explicitly instructs the next session to
prefer the forcing-matched workflow and to keep outputs provisional if the cache
is unavailable.

**Used in:** `src/evaluation/comparison.py`,
`outputs/comparison/production_candidate_assessment.json`

**Testable?** Yes. Once `data/raw/era5_cache/` is populated, the matched rerun
can either confirm the provisional leading candidate or overturn it. The
production-reproduction check should only move from `skipped` or `blocked` to
`passed` on that matched rerun.

**Risk:** If this boundary is ignored, a candidate could be promoted on the
basis of hand-built site proxies rather than the public weather-history
pipeline that production actually uses. That would hide forcing uncertainty
behind a false sense of finality.

**Status:** Active.

---

### AR-032 — The reusable matched-weather cache is filled from the Open-Meteo archive API with `models=era5`

**Assumption:** The reusable weather cache consumed by the matched comparison
path is populated from the Open-Meteo archive endpoint with `timezone = GMT`
and `models = era5`, using the hourly fields required by the clean-room
pipeline (`wind_speed_10m`, `wind_direction_10m`, `wind_gusts_10m`,
`surface_pressure`, `temperature_2m`, `shortwave_radiation`, `cloud_cover`,
`precipitation`).

**Source:** Open-Meteo archive API response schema, exercised through
`src/data/era5_cache.py`, and the existing cache consumer in
`src/evaluation/comparison.py`.

**Used in:** `src/data/era5_cache.py`, `src/evaluation/comparison.py`,
`data/raw/era5_cache/`

**Testable?** Yes. The cache filler has direct tests for request construction
and on-disk schema, and the matched comparison can be rerun deterministically
from the saved JSON files.

**Risk:** This is still a processed API layer rather than a direct CDS export.
If future work requires native ERA5 metadata or additional variables, the cache
writer may need to be replaced while preserving the consumer-facing schema.

**Status:** Active.

---

### AR-033 — Operational CL coarsening is seeded from the nonlinear onset width `L_inst = 2πh / l_cNL`

**Assumption:** In the operational CL candidate, the nonlinear critical
wavenumber `l_cNL` defines the initial cell-pair width passed into the discrete
coarsening clock:
`L_cell,0 = 2π × depth / l_cNL`.
The model no longer replaces that onset width with a separate supercritical
selection from the upper edge of the unstable band before coarsening.

**Source:** The matched structural audit
(`outputs/structural_audit_matched/summary.json`) showed that the raw onset
width already carries material lake-to-lake and forcing-driven variance, while
the previous `upper_unstable_band` selector crushed that variance by an almost
fixed factor of `~1/8`. The matched rerun with onset-seeded coarsening improved
the CL full-set metrics materially (`RMSE: 126.4 → 116.6 m`,
`range coverage: 0.137 → 0.246`).

**Used in:** `src/prediction/candidate_cl.py`,
`src/evaluation/comparison.py`

**Testable?** Yes. The full CL path should now share its initial width with the
`onset_only` diagnostic mode and differ only through coarsening and
visibility-layer transformations. Comparison reruns can test whether this
restored onset variance improves coverage without introducing a new attractor.

**Risk:** If the physically relevant supercritical comparison scale should sit
systematically above the onset minimum, this choice may still under-predict the
widest observed spacings even though it restores lost variance.

**Status:** Active.

---

### AR-034 — Matched organisation time is counted from the latest actual reset event, not repeatedly re-stamped by older events still inside the lookback window

**Assumption:** On the matched weather-history path, the organisation time used
for coarsening is the elapsed time since the latest actual reset event in the
selected history window. A reset event is either the most recent subcritical
forcing sample or the latest weather-disruption sample
(`direction_change`, `low_wind_shutdown`, or `rapid_speed_increase`) that
occurs on the current post-reset episode. Once a reset happens, later steady
samples are evaluated only against the new episode rather than against stale
pre-reset history.

**Source:** WP-06 lifetime audit in
`outputs/comparison_matched_lifetime_fix/organisation_time_audit.json` and the
event-local regression tests added in
`src/hydro/tests/test_support_modules.py` and
`src/prediction/tests/test_pipeline.py`.

**Used in:** `src/hydro/coarsening.py`, `src/prediction/pipeline.py`,
`src/evaluation/comparison.py`

**Testable?** Yes. Synthetic histories with one early reset followed by steady
winds should produce organisation times that grow from the actual event time
instead of collapsing back to zero at later samples.

**Risk:** If real structures require additional decay or re-formation lag after
the detected reset sample, this rule can still overestimate the usable
organisation time by up to roughly one sampling interval.

**Status:** Active.

---

### AR-035 — The coarsening clock uses a lateral mixing-length diffusivity `A_H = u*_water × depth` rather than the forcing-layer depth-mean vertical `nu_T`

**Assumption:** The discrete merger clock models horizontal cell-pair
reorganisation, so it should be dimensionalised with a lateral turbulent
diffusivity rather than with the depth-mean vertical eddy viscosity used in the
forcing-to-`Ra` mapping. The operational coarsening closure therefore uses
`A_H = max(nu_T_vertical, u*_water × depth)` in
`τ_merge = λ² / ((2π)² A_H)`. Across the matched WP-06 cases, the lateral term
dominates in every case, so the practical closure is `A_H = u*_water × depth`.

**Source:** WP-06 coarsening-closure audit in
`outputs/comparison_matched_coarsening_closure_fix/coarsening_closure_audit.json`.
That audit shows the forcing-layer `nu_T` distribution is
`1.5e-4–7.6e-3 m²/s` with a `1.05e-3 m²/s` median, while the lateral
coarsening diffusivity is `2.2e-3–1.11e-1 m²/s` with a `1.54e-2 m²/s`
median. The anisotropy ratio is a constant `6/κ ≈ 14.63`, and the CL
first-merger timescale falls from a `24 h` median in the previous matched run
to `1.64 h`, recovering the `O(h/u*)` order anticipated in the rebuild brief.

**Used in:** `src/hydro/coarsening.py`,
`src/prediction/candidate_cl.py`,
`src/prediction/candidate_scaling.py`,
`src/evaluation/comparison.py`

**Testable?** Yes. Regression tests check that the coarsening diffusivity is
reported separately from the forcing-layer `nu_T` and that it equals
`u*_water × depth` for the operational closure. Matched comparison audits can
verify that `A_H >= nu_T_vertical` in every case and that the resulting
first-merger timescale is in the `O(1–4 h)` range rather than `O(1 day)`.

**Risk:** This closure is physically more defensible than reusing the vertical
`nu_T`, but it is not yet a production-clearing fix. In the matched WP-06
rerun it frees the CL clock enough to produce `0/1/2` mergers in `15/31/21`
cases, but the public comparison still over-expands the CL range and worsens
`BM-A_full` RMSE (`116.6 → 159.9 m`). If real lateral merger rates depend more
strongly on wave state, anisotropy, or finite-amplitude saturation than this
mixing-length closure captures, the present `A_H` can still over-coarsen.

**Status:** Active.

---

## Multi-scale Structure Predictions

### AR-036 — Wave-tied LC spacing coefficient C_wave = 0.34

**Assumption:** The primary LC cell-pair width in fully developed Langmuir
turbulence scales as d_LC = C_wave x lambda_p, where C_wave = 0.34 and
lambda_p is the peak surface wave wavelength. This is the dominant spacing of
counter-rotating LC vortex pairs driven by the CL-II instability mechanism.

**Source:** Tsai & Lu (2023), *J. Fluid Mech.* 969, A30. DNS at Re_tau ~ 530
gives d_s / lambda_wave ~ 0.34 (nondimensional wavenumber l_s ~ 2.9),
consistent across wave steepness ak = 0.135 to 0.22.

**Used in:** `src/hydro/multiscale_structures.py`,
`src/prediction/candidate_multiscale.py`

**Testable?** Yes. C_wave can be directly compared against the DNS result.
The scaling d_LC ~ lambda_p can be tested via the dependence of predicted
spacing on wind speed and fetch (which control lambda_p). If manual and
stream observations vary as lambda_p, this assumption is supported.

**Risk:** The coefficient was determined at laboratory-scale Re and for
specific wave steepness values. At field-scale Re (O(10^6)) or in
shallow-water conditions (k_p H < pi), the coefficient may differ. Treating
C_wave as a constant may miss La_t-dependent or depth-dependent corrections.

**Status:** Active.

### AR-037 — Coarsening initialised from wave-tied scale d_LC, not CL onset L_inst

**Assumption:** The multiscale candidate initialises the Y-junction merger
clock from d_LC = C_wave x lambda_p rather than from L_inst = 2*pi*h / l_cNL.
The physical rationale is that in the supercritical regime (Ra >> R_c), the
CL instability sets the large-scale full-depth roll organization, but the
primary visible LC convergence lines form at the wave-tied scale and coarsen
through discrete mergers from that smaller initial width.

**Source:** Tsai & Lu (2023) for the wave-tied scale; Xuan & Shen (2025),
*J. Fluid Mech.* 1023, A4, resolvent analysis showing that in fully developed
turbulence the CL instability scale represents the full-depth roll, while
smaller structures coexist near the surface and set the visible convergence
pattern.

**Used in:** `src/prediction/candidate_multiscale.py`

**Testable?** Yes. If coarsening from d_LC produces predictions in the right
range (30-120 m for typical Neagh conditions) while coarsening from L_inst
overshoots, this assumption is supported. The CL candidate (using L_inst) is
preserved for direct A/B comparison.

**Risk:** If d_LC is too small, the coarsening model must produce many
mergers (5+) to reach observed scales, making predictions very sensitive to
pattern_lifetime and A_H. If the merger timescale is wrong, predictions
will cluster at the wrong scale.

**Status:** Active.

### AR-038 — CL onset scale represents full-depth roll organization

**Assumption:** The CL instability onset width L_inst = 2*pi*h / l_cNL
represents the wavelength of the full-depth counter-rotating roll structure,
not the spacing of satellite-visible surface convergence lines. In the
supercritical regime, L_inst is a real physical scale but is not directly
comparable to manual or stream observation categories. It may correspond to
the coarser spacing detected by wiggle spectral methods.

**Source:** Xuan & Shen (2025), *J. Fluid Mech.* 1023, A4. Resolvent analysis
shows that streamwise-invariant full-depth rolls (lambda_z ~ H to 2H) are the
most amplified large-scale structures, while smaller near-surface structures
have distinct (smaller) spanwise wavelengths.

**Used in:** `src/prediction/candidate_multiscale.py` (as a diagnostic
reference, not as the coarsening initial condition)

**Testable?** Yes. If wiggle observations (median ~214 m) correlate better
with L_inst than with the coarsened wave-tied scale, this interpretation is
supported. Layer-by-layer Spearman rho analysis can test this per category.

**Risk:** If the CL onset scale is not physically realized in the flow (i.e.
the instability saturates at a different wavelength), retaining it as a
diagnostic may be misleading. The nonlinear solver's l_cNL was validated
against H&P (2017) benchmarks but not against field observations of full-depth
roll spacing.

**Status:** Active.

---

## Superseded Entries

- **AR-001** — Superseded by AR-028. The free-surface Robin parameter is no
  longer a fixed constant in the operational CL path.
- **AR-002** — Superseded by AR-028. The bottom Robin parameter is no longer a
  fixed constant in the operational CL path.
