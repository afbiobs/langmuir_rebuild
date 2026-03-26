# AGENTS.md — Coding Governance for the Langmuir Clean-Room Rebuild

**Scope:** This document governs all code written for the `langmuir_rebuild/` project. Every agentic coder working on this project must read this file before writing any code and must treat its rules as hard constraints, not suggestions.

**Why this document exists:** The predecessor to this project failed not because of any single bug but because of accumulated structural decisions that looked locally reasonable but collectively produced a model that could not respond to its own inputs. This document encodes the lessons from that failure as enforceable rules.

---

## Context Efficiency Rules (for agentic coders)

These rules apply to any LLM agent working in this codebase. Violating them wastes context budget and can cause output truncation mid-implementation.

1. **Never read a whole file when you only need a section.** Use `offset` and `limit` parameters on the Read tool. Large reference documents (`hayes_phillips_2017.md`, `clean_room_rebuild_brief.md`, AGENTS.md itself) cost 5–27k tokens each if read whole.

2. **Read docs once per session, then work from memory.** After reading a reference document, do not re-read it. Take notes in your working memory during the first read.

3. **Read source files before editing.** But use targeted reads: if you only need to check a function signature, read 20–30 lines around it, not the whole file.

4. **Prefer Grep over Read for locating content.** Use `Grep` to find which file/line contains what you need, then `Read` only that section.

5. **Write implementations in one pass per file.** Plan the full file content before calling Write. Do not write a stub and then immediately re-read and re-edit.

6. **Break large implementations into separate files.** Each file is a separate Write call. Do not attempt to write the entire hydro module in one response.

---

## Part I — Foundational Principles

### 1. The Hierarchy of Concerns

When making any decision, apply these priorities in order. A lower-priority concern never overrides a higher one.

```
1. Physical correctness   — Does this respect conservation laws, dimensional
                            consistency, and known constraints?
2. Falsifiability         — Can this choice be shown wrong with available data?
3. Transparency           — Can a reviewer see exactly why a prediction was made?
4. Robustness             — Does this work across the observed input range
                            without hidden saturation or clipping?
5. Accuracy               — Does this reduce prediction error?
6. Elegance               — Is this simple and clean?
```

**Accuracy is fifth, not first.** The old model repeatedly improved aggregate RMSE while becoming physically less informative. A model that is transparent, falsifiable, and robust but slightly less accurate is strictly preferred over one that is accurate but opaque and fragile.

### 2. The Core Failure Mode: Silent Saturation

The single most dangerous pattern in this project is a quantity that stops varying in response to its inputs without anyone noticing. This killed the old model at least three separate times: in the anchor spacing, in the stage occupancy distribution, and in the effective forcing control parameter.

Silent saturation occurs when:

- A quantity is normalised and the normalisation absorbs the dynamic range.
- A quantity is passed through a sigmoid, tanh, soft clip, or min/max that compresses it into a narrow band.
- A quantity is computed as a ratio of two co-varying terms that partially cancel.
- A quantity depends on a weighted average where the weights are themselves nearly invariant.
- A quantity is correct at two extremes but passes through a flat region in between.

**The saturation rule:** Every function that returns a physical or diagnostic quantity must be auditable for input-output sensitivity. If the output varies less than 5% over more than 50% of the physically plausible input range, the function is considered saturated and must be either fixed or explicitly justified in the assumptions register.

This is not a guideline. It is a hard test that runs automatically.

### 3. The Separation Principle

This project models a system with at least five conceptually distinct layers:

```
External forcing  →  Hydrodynamic response  →  Coherent structure state
                                                        ↓
                                            Visibility / tracer expression
                                                        ↓
                                            Cyanobacteria accumulation
```

**The separation rule:** These layers must be implemented as independent modules that can be called, tested, and validated in isolation. No layer may import from a layer that depends on it (no circular dependencies). No layer may assume the internal state representation of another layer.

The interface between layers is a data structure with named, typed, unit-documented fields — not a tuple, not a dict with ad-hoc keys, not a global variable.

Why: the old model blended these layers into a small number of state variables. When the output collapsed, it was impossible to determine which layer was responsible. Independent testability is not a software nicety; it is a scientific requirement.

### 4. The Scale Distinction

This project involves at least four distinct spatial scales that must never be silently conflated:

| Scale | What it is | Where it comes from |
|---|---|---|
| Instability scale | The wavelength selected by the CL instability mechanism | Linear or nonlinear CL solver |
| Cell scale | The physical width of a mature LC cell pair | Instability scale + coarsening + finite-amplitude effects |
| Visible spacing | The spacing between surface features visible from satellite | Cell scale + tracer accumulation + observability filter |
| Accumulation scale | The spacing between cyanobacterial surface scum bands | Visible spacing + biological buoyancy + residence time |

Any function that converts between these scales must:

1. Be a named, documented function (not an inline multiplication).
2. State which two scales it converts between.
3. State what assumptions underlie the conversion.
4. Be registered in the assumptions register.

**The conflation test:** If you find yourself writing `spacing = 2 * pi * depth / l_c` without specifying which of the four scales `spacing` represents, you are conflating scales. Stop, name the output, and document the conversion.

---

## Part II — Code Patterns and Anti-Patterns

### 5. Mandatory Code Patterns

#### 5.1 Every function documents its units

```python
# CORRECT
def friction_velocity(U10: float) -> float:
    """
    Compute friction velocity from 10-metre wind speed.

    Parameters:
        U10: Wind speed at 10 m height [m/s]

    Returns:
        u_star: Friction velocity [m/s]

    Source: COARE 3.5 (Fairall et al. 2003)
    Assumptions: Neutral stability, open water fetch
    """

# WRONG — no units, no source
def friction_velocity(U10):
    """Compute u_star from wind speed."""
```

This is enforced by the automated test `test_dimensional_consistency.py`, which parses docstrings for unit annotations.

#### 5.2 Dataclasses for layer interfaces

```python
from dataclasses import dataclass
from datetime import datetime
import numpy as np

@dataclass(frozen=True)
class ForcingState:
    """
    Complete description of external forcing at a single time.
    All fields dimensional with documented units.
    Immutable once created — no downstream mutation.
    """
    U10: float              # [m/s] 10-metre wind speed
    u_star: float           # [m/s] friction velocity
    U_surface: float        # [m/s] surface current velocity
    # ... etc

    # PROVENANCE — not physics, just audit trail
    drag_method: str        # which drag model was used
    drift_method: str       # which Stokes drift model was used
```

Use `frozen=True` to prevent any downstream code from mutating a forcing state. If a downstream layer needs to modify something, it creates a new object.

#### 5.3 Explicit conversion functions between scales

```python
def instability_scale_to_cell_width(
    l_c: float,
    depth: float,
    n_coarsening_events: int = 0,
) -> float:
    """
    Convert nondimensional critical wavenumber to dimensional cell width.

    Input scale:  INSTABILITY SCALE (from CL solver, nondimensional)
    Output scale: CELL SCALE (physical width of one cell pair) [m]

    Conversion:
        base_width = 2 * pi * depth / l_c        [m]
        cell_width = base_width * 2^n_coarsening  [m]

    Assumptions:
        - Coarsening doubles width per event (Y-junction merging)
        - n_coarsening_events is estimated externally

    Registered as assumption A07 in assumptions_register.md
    """
    base_width = 2.0 * np.pi * depth / l_c
    return base_width * (2.0 ** n_coarsening_events)
```

#### 5.4 Raw values always preserved alongside any normalisation

```python
# CORRECT — raw value accessible
result = {
    "growth_rate_raw": sigma_r,           # [1/s]
    "growth_rate_max": sigma_r_max,       # [1/s]
    "growth_rate_normalised": sigma_r / sigma_r_max,  # [dimensionless]
}

# WRONG — raw value discarded
result = {
    "growth_rate": sigma_r / sigma_r_max,  # where did the amplitude go?
}
```

If you normalise, the raw numerator and denominator must both be in the output. If a downstream consumer only uses the normalised form, that is the downstream consumer's choice — the information must not be destroyed at source.

#### 5.5 Test-first for published verification targets

```python
# Write this BEFORE implementing the solver:

def test_kappa_uniform_profiles():
    """
    Hayes & Phillips (2017), Section 7.2:
    For D' = U' = 1, kappa ≈ 1.425.
    """
    result = nonlinear_solver.compute_kappa(
        profile=ShearDriftProfile(a_coeffs=[1], b_coeffs=[1]),
        gamma_s=0.0001,
        gamma_b=0.0,
    )
    assert abs(result - 1.425) < 0.005, f"kappa = {result}, expected ~1.425"
```

The test exists and fails before the solver is written. The solver is complete when the test passes. Do not adjust the test target to match the solver output.

### 6. Forbidden Anti-Patterns

#### 6.1 NEVER: Hidden bounding

```python
# FORBIDDEN — sigmoid compresses the output without physical justification
def effective_spacing(raw_spacing):
    return 40 + 120 / (1 + np.exp(-0.05 * (raw_spacing - 100)))

# FORBIDDEN — clip masks that the physics produced an out-of-range value
def effective_spacing(raw_spacing):
    return np.clip(raw_spacing, 40, 160)
```

If the physics produces an out-of-range value, that is diagnostic information. Clipping it hides a problem. Instead:

```python
# CORRECT — flag the issue, return the raw value, let the caller decide
def effective_spacing(raw_spacing):
    if raw_spacing < 20 or raw_spacing > 300:
        warnings.warn(
            f"Spacing {raw_spacing:.1f} m outside expected range [20, 300]. "
            "Check forcing inputs and solver convergence."
        )
    return raw_spacing  # unmodified
```

#### 6.2 NEVER: Compensatory memory wrappers

```python
# FORBIDDEN — adding exponential smoothing to fix a dynamical problem
spacing = alpha * previous_spacing + (1 - alpha) * instantaneous_spacing
```

If the instantaneous prediction is too volatile, the problem is in the instantaneous prediction, not in the smoothing. Memory is only acceptable as an explicit state variable with a physically justified timescale:

```python
# ACCEPTABLE — explicit state with justified timescale
@dataclass
class StructureState:
    cell_width: float       # [m] current cell width
    time_since_onset: float # [s] time since LC onset
    coarsening_tau: float   # [s] timescale for next merger = O(h / u*)

def evolve_structure(state: StructureState, dt: float,
                     forcing: ForcingState) -> StructureState:
    """
    Advance the structure state by dt seconds.
    Coarsening occurs discretely when time_since_onset exceeds coarsening_tau.
    """
```

#### 6.3 NEVER: Chaining untested components

```python
# FORBIDDEN — three transformations chained without intermediate validation
spacing = observation_model(
    coarsened_width(
        cl_solver(
            forcing_to_rayleigh(weather_data)
        )
    )
)
```

Each link in the chain must have been independently tested:

```python
# CORRECT — each step validated separately, intermediates preserved
ra = forcing_to_rayleigh(weather_data)
assert ra > 0, "Rayleigh number must be positive"

cl_result = cl_solver(ra, profile, bcs)
assert cl_result.lcNL > 0, "Critical wavenumber must be positive"

cell_width = instability_scale_to_cell_width(cl_result.lcNL, depth)
assert 10 < cell_width < 500, f"Cell width {cell_width} m outside plausible range"

visible = observation_model(cell_width, coarsening_events)
```

#### 6.4 NEVER: Tuning to aggregate error

```python
# FORBIDDEN — adjusting a parameter because it improves RMSE
gamma_s = 0.08  # tuned from 0.06 to improve RMSE by 2 m

# REQUIRED — adjusting a parameter only with physical justification
gamma_s = 0.06  # Cox & Leibovich (1993), physical estimate
# Sensitivity: lc ~ gamma^{1/4}, so 0.06 → 0.08 changes lc by ~7%
# Decision: retain literature value. See parameter_register.md row P03.
```

If a parameter change improves RMSE but has no physical justification, it is overfitting. Document it in the failure log and do not apply it.

#### 6.5 NEVER: Implicit resolution of ambiguity

```python
# FORBIDDEN — silently assuming visible spacing equals cell spacing
predicted_spacing = 2 * np.pi * depth / l_c  # which scale is this?

# REQUIRED — explicit about the assumption
cell_spacing = instability_scale_to_cell_width(l_c, depth)
visible_spacing = cell_spacing_to_visible(
    cell_spacing,
    model_type="direct",  # assumes 1:1 mapping — see observation_model.md
)
```

---

## Part III — Workflow Discipline

### 7. The Work Package Protocol

This project proceeds through numbered work packages (WPs) with decision gates. The protocol for each WP is:

```
1. READ the WP specification in the engineering plan.
2. CHECK entry conditions. If not met, stop and report.
3. WRITE tests first (for code WPs) or outline first (for documentation WPs).
4. IMPLEMENT the deliverable.
5. RUN all tests, including:
   - Unit tests for the current WP
   - Saturation audit for any new functions
   - Dimensional consistency check
   - Integration tests if this WP connects to previous WPs
6. UPDATE living documents:
   - assumptions_register.md (any new assumptions)
   - parameter_register.md (any new parameters)
   - failure_log.md (anything that didn't work)
   - decision_log.md (any significant choices)
7. CHECK exit conditions. If not met, iterate.
8. REPORT completion and present outputs.
9. STOP at decision gates. Do not proceed without human approval.
```

### 8. Decision Gates Are Hard Stops

There are four decision gates in the engineering plan. At each gate:

1. **Stop coding.**
2. Present the deliverables produced so far.
3. Summarise what worked, what didn't, and what is uncertain.
4. State the specific questions that require human judgement.
5. **Wait for human response before proceeding.**

Do not:

- Assume the answer to a gate question and continue.
- Skip a gate because "the answer seems obvious."
- Batch multiple gates together.
- Proceed with one candidate while waiting for the gate decision on another.

### 9. The Living Documents

Four documents are updated continuously throughout the project. They are not afterthoughts — they are primary deliverables.

#### 9.1 `assumptions_register.md`

Every modelling assumption gets a row. Format:

```markdown
| ID   | Assumption | Source | Falsification criterion | Status |
|------|-----------|--------|------------------------|--------|
| A01  | Robin BC surface parameter γ_s ≈ 0.06 | Cox & Leibovich 1993 | If spacing is insensitive to γ_s across the observed range, this value is not constrained by data | Active |
| A02  | Visible spacing = cell spacing (1:1) | WP-02 default | If model systematically under- or over-predicts by a consistent factor ≠ 1 | Active |
| A03  | Coarsening doubles width per event | Thorpe 2004 Y-junction merging | If observed spacing shows no doubling structure | Active |
```

**When to add a row:** Every time you write code that depends on a choice that could reasonably have been made differently.

**Status values:** `Active`, `Uncertain`, `Falsified`, `Superseded`

#### 9.2 `parameter_register.md`

Every parameter that is not derived from first principles. Format:

```markdown
| ID  | Parameter | Value | Units | Source | How it enters | Sensitivity | Status |
|-----|-----------|-------|-------|--------|--------------|-------------|--------|
| P01 | γ_s | 0.06 | — | Cox & Leibovich 1993 | lc = (γ R̃₂/R*₂)^{1/4} | lc ~ γ^{1/4}: ×2 in γ → ×1.19 in lc | Fixed |
| P02 | γ_b | 0.28 | — | Cox & Leibovich 1993 | Same as P01 | Same as P01 | Fixed |
| P03 | C_d (low wind) | 1.0e-3 | — | Sheltered lake estimate | u* = C_d^{1/2} × U10 | u* ~ C_d^{1/2}: ×2 in C_d → ×1.41 in u* | Estimated |
```

**Status values:** `Fixed` (from literature, not adjustable), `Estimated` (defensible but uncertain), `Calibrated` (fit to data — requires strong justification), `Legacy` (carried from old model — must be re-derived or removed)

**No parameter may have status `Legacy`.** If a value comes from the old model, it must be re-derived from literature or flagged as `Estimated` with an independent source.

#### 9.3 `failure_log.md`

Everything that did not work. Format:

```markdown
| Date | WP | What was attempted | What happened | Root cause | Resolution |
|------|----|--------------------|---------------|------------|------------|
| ... | WP-07 | Galerkin solver with J=8 Legendre modes | Newton iteration diverged for l > 0.2 | Insufficient resolution for higher wavenumbers | Increased to J=13 per paper recommendation |
```

**This log has no shame.** Failures are information. An empty failure log is suspicious, not virtuous.

#### 9.4 `decision_log.md`

Every architectural or modelling choice. Format:

```markdown
| Date | WP | Decision | Alternatives considered | Rationale | Reversible? |
|------|----|---------|-----------------------|-----------|-------------|
| ... | WP-06 | Use COARE 3.5 for drag | Smith 1980, constant C_d | COARE 3.5 is continuous (no step functions), well-validated for moderate winds | Yes — drag model is a parameter |
```

### 10. The Saturation Audit Protocol

After completing any module that contains functions returning physical or diagnostic quantities, run the saturation audit:

```python
def saturation_audit(func, input_ranges: dict, n_samples: int = 100,
                     threshold_variation: float = 0.05,
                     threshold_fraction: float = 0.50) -> dict:
    """
    Sweep func over its input ranges and check for saturation.

    For each input dimension:
    1. Hold all other inputs at their midpoint.
    2. Vary this input linearly from min to max (n_samples points).
    3. Compute the output at each point.
    4. Compute output_range = (max - min) / mean.
    5. Find the longest contiguous sub-interval where
       output varies < threshold_variation * output_range.
    6. If that sub-interval spans > threshold_fraction of the input range,
       FLAG this input-output pair as saturated.

    Returns:
        dict mapping (input_name, output_name) → {
            "saturated": bool,
            "flat_fraction": float,  # fraction of input range that is flat
            "output_range": float,   # total output variation
            "input_range": tuple,    # (min, max) of input
        }
    """
```

**If any output is flagged as saturated:**

1. Stop and investigate.
2. Is the saturation physically correct? (e.g., spacing should not respond to wind speed in the subcritical regime — that flat region is real.)
3. If yes: document in assumptions register with justification.
4. If no: fix the function before proceeding.
5. If uncertain: flag as `Uncertain` in assumptions register and continue, but add a note to revisit at the next decision gate.

---

## Part IV — Physics-Specific Guidance

### 11. The CL Solver Is Not The Model

The nonlinear CL solver (Hayes & Phillips 2017) is a component of Candidate C. It is not the model. The model includes the forcing layer, the hydro layer (of which the CL solver is part), the structure layer, the visibility layer, and the accumulation layer. The CL solver's output is a nondimensional critical wavenumber $l_{c_{NL}}$. Converting that to a dimensional spacing prediction requires the forcing-to-Ra mapping, the scale conversion functions, the coarsening model, and the observation model. Each of those components has its own assumptions and failure modes.

**Practical consequence:** Do not spend disproportionate effort on the CL solver at the expense of the other layers. A perfectly accurate CL solver with a broken observation model will produce wrong predictions. Budget effort according to uncertainty, not mathematical complexity.

### 12. The Robin Boundary Conditions Are Not Optional

The paper proves that Neumann boundary conditions ($\gamma = 0$) predict onset at $l = 0$ regardless of nonlinearities (Chapman & Proctor 1980). This means infinite cell spacing — physically meaningless. Robin boundary conditions are what ensure a finite, nonzero preferred wavenumber at onset.

**Practical consequence:** If the Robin BC implementation is wrong, every downstream prediction is wrong. Test this component with extreme care. Specifically verify that:

- With $\gamma = 0$ (Neumann): the neutral curve $R(l)$ is monotonically decreasing. There is no minimum except at $l = 0$.
- With $\gamma > 0$ (Robin): the neutral curve has a well-defined minimum at finite $l > 0$.
- The location of the minimum scales as $\gamma^{1/4}$.
- The value $\kappa = l_{c_L} / l_{c_{NL}}$ is independent of $\gamma$.

### 13. The Shear and Drift Profiles Control The Answer

The ratio $\kappa$ ranges from 1.425 to 1.944 depending on the vertical profiles of $U'(z)$ and $D'(z)$. For the same boundary conditions, changing the profiles can change the predicted cell width by nearly 40%. Furthermore, the paper shows that the shear profile $U'$ is more influential than the drift profile $D'$ in suppressing instability (curves d vs e in their Figure 2).

**Practical consequence:** The profiles chosen for Lough Neagh must be physically justified, not defaulted to uniform ($D' = U' = 1$). Uniform profiles are the simplest case and useful for verification, but they do not represent a wind-driven shallow lake. The `shallow_lake_profile` factory function in `profiles.py` must be implemented with care, and its sensitivity to input parameters must be characterised.

### 14. The Rayleigh Number Is The Master Control

Every downstream quantity ultimately depends on $Ra = U D h^2 / \nu_T^2$. Each of the four terms in this expression carries uncertainty:

| Term | Source of uncertainty | Order of magnitude |
|---|---|---|
| $U$ (surface velocity) | Depends on drag model, closed-basin assumption | Factor of 2–3 |
| $D$ (drift magnitude) | Depends on wave model, shallow-water corrections | Factor of 2–3 |
| $h$ (depth) | Usually well known | <10% |
| $\nu_T$ (eddy viscosity) | Least constrained. Parabolic vs KPP can differ by factor of 3 | Factor of 2–5 |

Because $Ra \propto \nu_T^{-2}$, the eddy viscosity enters quadratically. A factor-of-2 error in $\nu_T$ produces a factor-of-4 error in $Ra$. This means the regime classification (subcritical/near-onset/supercritical) can shift dramatically with modest changes to the turbulence closure.

**Practical consequence:**

1. Always report $Ra$ and its constituent terms alongside every prediction.
2. Characterise the sensitivity of the final prediction to $\nu_T$ explicitly.
3. Consider running predictions at $\nu_T$, $0.5 \nu_T$, and $2 \nu_T$ to bracket the uncertainty.
4. If the answer is "subcritical at $\nu_T$ but supercritical at $0.5 \nu_T$", the regime is genuinely uncertain and the model should say so.

### 15. The Bio Module Computes Enhancement, Not Accumulation

The biological component of this model does not predict blooms, scum formation, colony dynamics, or biomass distribution. It computes how much Langmuir circulation enhances growing conditions relative to a static water column.

**The Triple Rule:** Every biological metric is a named triple: `(value_LC, value_static, enhancement_ratio)`. A function that returns only the LC value is incomplete. A function that returns only the ratio has discarded the absolute values. Both are forbidden. The automated test suite enforces this: any scalar return from `lc_enhancement.py` that represents a biological quantity is a test failure.

The three enhancement dimensions are:

1. **Light exposure.** LC cycles cells through the photic zone, avoiding photoinhibition at the surface and light starvation at depth. The static column has buoyant cells pooling at the surface under photoinhibiting irradiance. Enhancement ratio = mean irradiance experienced by a cycling cell / effective irradiance of a surface-pooled cell.

2. **Nutrient upwelling.** LC drives vertical circulation that replenishes surface nutrients from depth. The static column relies on molecular diffusion (orders of magnitude slower). Enhancement ratio = nutrient flux with LC / nutrient flux without LC.

3. **Temperature distribution.** LC homogenises the mixed layer. The static column stratifies (hot surface, cool bottom). Enhancement ratio = proximity of mixed-layer temperature to growth optimum / proximity of surface temperature to growth optimum. This can be <1 if the surface temperature is already near optimum and mixing cools it away.

**The Comparison Principle:** The model answers "how much better are conditions with LC than without?" — not "will there be a bloom?" Whether enhanced growing conditions produce a bloom depends on initial biomass, grazing, nutrient loading, and other factors outside this model's scope.

**The Scope Boundary:** The model's biological relevance ends at enhancement ratios. It does not extend to predicting biomass concentration, bloom timing, toxin production, or spatial scum distribution. These are downstream questions. The LC enhancement calculator provides one input to those questions, not the answer.

### 16. Coarsening Is Discrete, Not Continuous

LC cells merge through Y-junction events (Thorpe 2004). Each merger approximately doubles the cell width. This is a discrete process, not a continuous exponential growth. Modelling it as continuous (e.g., `spacing *= exp(rate * t)`) obscures the mechanism and makes it harder to constrain.

**Preferred implementation:**

```python
def count_coarsening_events(time_available: float,
                            tau_coarsening: float) -> int:
    """
    Number of discrete merger events that have occurred.
    Each event takes approximately tau_coarsening seconds.

    tau_coarsening = O(h / u*) for turbulent adjustment.
    """
    if tau_coarsening <= 0:
        return 0
    return int(time_available / tau_coarsening)

def coarsened_width(initial_width: float, n_events: int,
                    max_aspect_ratio: float = 12.0,
                    depth: float = 9.0) -> float:
    """
    Width after n coarsening events.
    Capped at max_aspect_ratio × depth.

    The cap is an ASSUMPTION (A-coarsening-cap) that must be registered.
    Justification: aspect ratios > 12 are not observed (Marmorino et al. 2005).
    """
    width = initial_width * (2.0 ** n_events)
    cap = max_aspect_ratio * depth
    if width > cap:
        warnings.warn(
            f"Coarsened width {width:.0f} m exceeds cap {cap:.0f} m. "
            f"Capping at aspect ratio {max_aspect_ratio}."
        )
        width = cap
    return width
```

Note that the cap is a `warnings.warn`, not a silent clip. The diagnostic information (how far past the cap) is preserved.

---

## Part V — Testing Philosophy

### 17. The Test Taxonomy

Tests in this project fall into five categories. All are mandatory.

#### 16.1 Verification tests (against published results)

These test that the solver reproduces known results from the literature. They are written before the solver exists. They define success.

```python
class TestHayesPhillips2017:
    """Verification against Hayes & Phillips (2017) published values."""

    def test_R0_uniform(self): ...
    def test_critical_linear_values(self): ...
    def test_critical_nonlinear_values(self): ...
    def test_kappa_all_profiles(self): ...
    def test_supercritical_stability(self): ...
    def test_aspect_ratio_range(self): ...
    def test_asymptotic_numeric_agreement(self): ...
```

**Rule:** Verification test targets are never adjusted to match the implementation. If the implementation does not match the published value, the implementation is wrong.

#### 16.2 Dimensional consistency tests

Every function that accepts or returns physical quantities is tested for correct dimensions.

```python
def test_rayleigh_number_dimensions():
    """Ra = U * D * h^2 / nu_T^2 should be dimensionless."""
    # Compute Ra with known inputs
    # Verify the result is dimensionless (no residual units)
    # Verify it scales correctly: doubling h quadruples Ra
```

#### 16.3 Monotonicity and constraint tests

Physical monotonicity relationships that must hold regardless of parameter values.

```python
def test_wider_cells_at_lower_wind():
    """For fixed depth, lower U10 → wider cells (at least in the near-onset regime)."""

def test_supercritical_above_subcritical():
    """R_bar(l) >= R(l) for all l > 0 with Robin BCs."""

def test_zero_mass_flux():
    """Surface current profile integrates to zero in closed basin."""

def test_instability_requires_same_sign():
    """D'U' > 0 is necessary for instability."""
```

#### 16.4 Saturation tests

Automated sweep tests that detect hidden saturation (see section 10).

```python
def test_forcing_layer_no_saturation():
    """Sweep U10 from 1 to 15 m/s. No output saturates."""

def test_hydro_layer_no_saturation():
    """Sweep Ra from 50 to 2000. Predicted spacing varies throughout."""

def test_candidate_c_no_attractor():
    """Run Candidate C on full observation set.
    Predictions must not cluster in <20% of observed range."""
```

#### 16.5 Integration tests

End-to-end tests that verify the full pipeline produces sensible results.

```python
def test_pipeline_subcritical():
    """U10 = 1 m/s, depth = 9 m → subcritical regime → no LC prediction."""

def test_pipeline_moderate_wind():
    """U10 = 5 m/s, depth = 9 m → moderate regime → spacing in [60, 120] m."""

def test_pipeline_strong_wind():
    """U10 = 10 m/s, depth = 9 m → supercritical → spacing in [30, 80] m."""

def test_pipeline_sensitivity():
    """±10% U10 perturbation → detectable change in output."""
```

### 18. The "Explain Your Prediction" Test

Every candidate model must be able to produce a one-sentence explanation for every prediction. This is not a cosmetic requirement — it is a diagnostic test.

```python
def explain_prediction(result: dict) -> str:
    """
    Generate a human-readable explanation of why this prediction was made.

    Example outputs:
        "Ra = 185 (near-onset). Nonlinear CL solver predicts lcNL = 0.72,
         giving cell width = 78 m at depth 9 m. No coarsening (insufficient
         time since onset). Visible spacing = 78 m (direct mapping)."

        "Ra = 45 (subcritical). No Langmuir circulation expected.
         Prediction: no coherent surface pattern."

        "Ra = 1200 (strongly supercritical). CL solver predicts lcNL = 1.1,
         giving cell width = 51 m. 2 coarsening events in 45 min since onset.
         Coarsened width = 204 m. Capped at 108 m (aspect ratio limit).
         WARNING: aspect ratio cap was binding — prediction may be unreliable."
    """
```

If a candidate model cannot produce such explanations, it fails the transparency requirement (Principle 1, priority 3).

---

## Part VI — Numerical Discipline

### 19. Floating-Point Hygiene

The CL solver involves iterated polynomial integrals and high-order expansions ($O(l^{16})$). Numerical issues are likely.

**Rules:**

1. **Use SymPy for the asymptotic expansion.** The polynomial coefficients at high orders involve large intermediate values that cancel. Symbolic computation avoids catastrophic cancellation. Convert to numerical only at the final evaluation step.

2. **Use double precision (float64) everywhere.** Do not use float32 for any physical computation.

3. **Monitor condition numbers.** The Galerkin system produces a matrix that may become ill-conditioned. Compute and report the condition number. If it exceeds $10^{10}$, flag a warning.

4. **Validate by cross-checking.** The asymptotic expansion and the Galerkin numerical solver must agree in their overlap region. If they disagree by more than 1%, something is wrong.

5. **Beware of $l = 0$ and $l \to 0$.** The Robin BC parameters $\tilde{\gamma}_s = \gamma_s / l^4$ diverge as $l \to 0$. The expansion is designed for this limit, but numerical evaluation at very small $l$ may still be problematic. Test at $l = 10^{-4}$ and $l = 10^{-2}$ as well as at $l = O(l_c)$.

### 20. Solver Convergence

The Galerkin Newton solver may fail to converge, especially far from the initial guess.

**Rules:**

1. **Always initialise from the linear solution.** The linear eigenfunctions provide the best starting point for the Newton iteration.

2. **Use continuation in $Ra$.** To trace the nonlinear neutral curve, start at $Ra$ slightly above the linear neutral curve and increase gradually. Use the converged solution at $Ra_n$ as the initial guess for $Ra_{n+1}$.

3. **Report convergence.** Every solver call must return a convergence flag and the number of iterations. Log any non-convergence in the failure log.

4. **Do not catch and silence convergence failures.** If the solver fails, the prediction is "solver did not converge" — not a default value.

```python
# CORRECT
try:
    result = nonlinear_solver.find_steady_state(l, Ra, profile, bcs)
except ConvergenceError as e:
    return PredictionResult(
        spacing=None,
        regime="solver_failure",
        explanation=f"Newton solver did not converge: {e}",
    )

# WRONG — silent fallback hides the failure
try:
    result = nonlinear_solver.find_steady_state(l, Ra, profile, bcs)
except ConvergenceError:
    result = linear_solver.find_neutral(l, Ra, profile, bcs)  # silent downgrade
```

---

## Part VII — When To Stop, When To Redirect

### 21. Recognising Structural Problems Early

The old model exhibited several warning signs that were noticed too late. Watch for these in the new model:

| Warning sign | What it means | What to do |
|---|---|---|
| Predictions cluster in a narrow band despite varied inputs | Hidden saturation or weak forcing sensitivity | Run saturation audit. Stop and diagnose before proceeding. |
| RMSE improves but dynamic range shrinks | Model is converging to the mean, not tracking variation | Reject the change. RMSE improvement without dynamic range improvement is not improvement. |
| A parameter change improves one subset but worsens another | The parameter is compensating for a structural error | Do not tune. Investigate the structural error. |
| Two candidates produce nearly identical outputs despite different physics | Either both are dominated by the same forcing feature, or both are insensitive to their own physics | Check whether forcing layer outputs vary enough. Check saturation in both candidates. |
| The observation model assumption changes the answer substantially | The model is not well constrained by the physics — it is constrained by the mapping assumption | This is important information. Report it at the decision gate. Do not hide it by choosing one assumption. |

### 22. When To Invoke Human Judgement

Beyond the formal decision gates, invoke human judgement when:

1. **A verification test fails and you cannot identify the bug.** Do not adjust the test or bypass it.
2. **The saturation audit flags a function and you are unsure whether the saturation is physical.** The human may have domain knowledge you lack.
3. **Two implementation approaches are both physically defensible and produce materially different results.** Document both, present the difference, and let the human choose.
4. **The literature is contradictory.** Present both sources and their implications. Do not silently prefer one.
5. **You are about to introduce a parameter that is not in the parameter register.** Every new parameter needs a source and a justification. If you cannot provide one, ask.

### 23. When To Abandon A Candidate

A candidate model should be abandoned (moved to `failure_log.md`) if:

1. It fails the attractor test (predictions cluster in <20% of observed range) and the saturation audit cannot identify a fixable cause.
2. It fails more than three verification tests after reasonable debugging effort.
3. It produces physically nonsensical intermediate states (e.g., negative eddy viscosity, negative cell width, infinite growth rate) that cannot be traced to a coding error.
4. Its explanation strings are uninformative or circular (e.g., "spacing is 80 m because the model predicted 80 m").

Abandoning a candidate is not failure — it is the comparison framework working as designed. Document what was learned and move on.

---

## Part VIII — Quick Reference Card

Print this section and keep it visible during coding sessions.

```
BEFORE WRITING ANY FUNCTION:
  □ What are the units of every input and output?
  □ What physical principle governs this computation?
  □ What assumption am I making? (Add to register)
  □ What test would show this function is wrong?
  □ Does this function normalise, clip, or bound anything? (Justify or remove)

BEFORE COMMITTING ANY CODE:
  □ All existing tests still pass
  □ New tests written for new functions
  □ Saturation audit run on new functions
  □ Dimensional consistency check passes
  □ Living documents updated (assumptions, parameters, failures, decisions)

BEFORE PROCEEDING TO NEXT WP:
  □ Entry conditions for next WP are met
  □ Exit conditions for current WP are met
  □ If a decision gate: STOP and present to human

WHEN SOMETHING GOES WRONG:
  □ Log it in failure_log.md FIRST
  □ Do not adjust test targets to match broken code
  □ Do not add compensatory wrappers
  □ Do not silently fall back to simpler methods
  □ Ask for help if stuck for more than 30 minutes on the same issue
```

---

## Appendix: File Header Template

Every Python file in `src/` must begin with this header:

```python
"""
Module: [module name]
Layer: [forcing | hydro | structure | visibility | bio | candidates | evaluation]
Purpose: [one sentence]

Assumptions made in this module:
    [A-XX]: [brief description] — see assumptions_register.md

Parameters used in this module:
    [P-XX]: [brief description] — see parameter_register.md

Key references:
    [Author Year] — [what is used from this reference]

Last reviewed: [date]
"""
```

This header exists so that anyone reading the file can immediately understand where it sits in the architecture, what assumptions it carries, and where to look for justification.
