# Code Engineering Plan: Clean-Room Langmuir Analysis Rebuild

**For execution by:** Frontier-level agentic LLM coder (Claude Code or equivalent)
**Governing document:** `clean_room_rebuild_brief.md`
**Physics references:** `hayes_phillips_2017.md`, background formulations document
**Constraint:** No code from the existing model may be copied. Every equation must be re-derived from literature or dimensional reasoning.

---

## 0. How To Read This Document

This plan is structured as a sequence of **work packages** (WPs), each with explicit entry conditions, deliverables, test criteria, and exit conditions. An agentic coder should execute them in order and **stop at any decision gate** that requires human judgement.

The plan follows the clean-room brief's recommended workflow but translates each phase into concrete engineering tasks. It also respects the brief's core warnings:

- Do not assume CL theory is the correct organising physics — build it as one candidate among several.
- Do not conflate hydrodynamic instability scale with observed spacing.
- Do not allow hidden saturation.
- Validate layers independently before coupling them.

The background physics formulations are referenced as **candidate inputs**, not mandatory components. Each must earn its place through the validation process.

---

## 1. Repository Structure

Create this structure at project initialisation. All new code lives here. The old model is frozen /home/op/PROJECTS/lang-nonlin as a reference-only dependency - weather api code and other non-biasing information can be recycled from old code.

```
langmuir_rebuild/
├── AGENTS.md                          # Coding governance (read before writing code)
├── README.md                          # Problem statement (WP-01 output)
├── docs/
│   ├── observation_model.md           # What "spacing" means (WP-01)
│   ├── assumptions_register.md        # Assumptions, parameters, decisions (living)
│   └── failure_log.md                 # What didn't work and why (living)
│
├── data/
│   ├── raw/                           # Observation CSVs, weather cache (reused)
│   └── benchmarks/                    # Defined subsets (WP-02)
│
├── src/
│   ├── __init__.py
│   │
│   ├── forcing/                       # Module 1
│   │   ├── __init__.py
│   │   ├── wind.py                    # Drag, friction velocity
│   │   ├── waves.py                   # Spectrum, Stokes drift profiles
│   │   ├── currents.py                # Closed-basin surface current
│   │   ├── eddy_viscosity.py          # ν_T models
│   │   └── tests/
│   │
│   ├── hydro/                         # Module 2
│   │   ├── __init__.py
│   │   ├── rayleigh.py                # Ra computation, regime classification
│   │   ├── profiles.py                # U'(z), D'(z) polynomial profiles
│   │   ├── robin_bc.py                # Robin boundary conditions
│   │   ├── linear_solver.py           # Linear CL perturbation solver
│   │   ├── nonlinear_solver.py        # Nonlinear CL steady-state solver
│   │   ├── galerkin.py                # Spectral method infrastructure
│   │   ├── scaling_laws.py            # La-dependent geometry (non-CL candidate)
│   │   ├── coarsening.py              # Y-junction merging, temporal evolution
│   │   └── tests/
│   │
│   ├── prediction/                    # Module 3 + calculators
│   │   ├── __init__.py
│   │   ├── candidate_cl.py            # CL-based prediction pipeline
│   │   ├── candidate_scaling.py       # Scaling-law prediction pipeline
│   │   ├── baseline.py                # Simple baselines (constant, linear-in-wind)
│   │   ├── visibility.py              # Visibility diagnostic (calculator)
│   │   ├── lc_enhancement.py          # LC enhancement calculator
│   │   ├── pipeline.py                # Entry point
│   │   └── tests/
│   │
│   └── evaluation/                    # Metrics, comparison, plots
│       ├── __init__.py
│       ├── metrics.py
│       ├── comparison.py
│       └── plots.py
│
└── tests/
    ├── test_saturation_audit.py
    ├── test_dimensional_consistency.py
    └── test_integration.py
```

---

### WP-01: Foundation Documents (No Code)

Write three short documents before any code is written.

**1. `README.md` — Problem statement (1 page)**

Must answer:

- What exactly is being predicted? (LC cell spacing under specified forcing)
- What is being diagnosed? (Whether LC enhances cyanobacterial growing conditions)
- What is out of scope? (Bloom prediction, biomass quantification, toxin production, scum distribution)
- What are the primary physical variables? (Ra, La_t, cell width, downwelling velocity)
- What are secondary observable proxies? (Satellite-visible windrow spacing)

Must not reference any variable names, function names, or state definitions from the old model.

**2. `docs/observation_model.md` — What "spacing" means (1 page)**

Must explicitly choose between (or acknowledge ambiguity among):

| Interpretation | Modelling implication |
|---|---|
| Dominant coherent LC cell spacing | Model predicts cell width directly |
| Visible streak/windrow spacing | Model must account for tracer accumulation |
| Mixed-scale visual composite | Model should predict a distribution |

Must also address: measurement uncertainty, directional biases (low-wind observations may favour wider visible patterns), and what a non-detection means.

**3. Initial `docs/assumptions_register.md`**

Seed with known assumptions:

- Robin BC parameters (γ_s ≈ 0.06, γ_b ≈ 0.28 from Cox & Leibovich 1993)
- Observation model interpretation choice
- Coarsening mechanism (Y-junction doubling, Thorpe 2004)
- Aspect ratio cap (~12, Marmorino et al. 2005)

**Exit condition:** Human reviews and approves. This gates all code.

---

### WP-02: Benchmarks and Metrics

Implement the evaluation framework before any physics code.

**Benchmark subsets from existing data:**

```python
BENCHMARK_SUBSETS = {
    "full": "All observations",
    "low_spacing": "Observed spacing < 60 m",
    "high_spacing": "Observed spacing > 120 m",
    "low_wind": "Representative U10 < 4 m/s",
    "high_wind": "Representative U10 > 8 m/s",
}
```

**Evaluation metrics (`src/evaluation/metrics.py`):**

```python
def spacing_rmse(predicted, observed): ...
def spacing_mae(predicted, observed): ...

def tail_coverage(predicted, observed, low_threshold=60, high_threshold=120):
    """Fraction of tail observations predicted within ±30%."""

def dynamic_range(predicted):
    """Ratio of 90th to 10th percentile of predictions."""

def attractor_test(predicted, observed_range, band_fraction=0.20):
    """FAIL if >50% of predictions fall in a band spanning
    <20% of the observed range. This is the primary structural
    health check — it catches the failure mode that killed
    the old model."""

def saturation_audit(func, input_ranges, threshold_variation=0.05,
                     threshold_fraction=0.50):
    """Sweep func over input ranges. Flag any output that varies
    <5% over >50% of the input range."""
```

**Baselines (`src/prediction/baseline.py`):**

```python
def baseline_constant(observations):
    """Predict the mean observed spacing for all cases."""

def baseline_linear_wind(wind_speeds, observations):
    """OLS regression of spacing on wind speed."""

def baseline_depth_scaled(wind_speeds, depths, observations):
    """OLS regression of spacing on wind speed and depth."""
```

Run baselines on all subsets. Record numbers. These are the floor to beat.

**Exit condition:** Metrics tested on synthetic data. Baselines produce numbers on real data.

---

### WP-03: Forcing Module

Implement `src/forcing/`. This module is shared by both candidates.

**Wind stress (`wind.py`):**

```python
def friction_velocity(U10: float, method: str = "coare35") -> float:
    """
    Compute u* from 10-metre wind speed.

    Parameters:
        U10: Wind speed at 10 m height [m/s]

    Returns:
        u_star: Friction velocity [m/s]

    Methods:
        "coare35": COARE 3.5 continuous formulation.
            Neutral drag: C_DN = [κ / ln(z/z0)]^2
            Charnock roughness: α = m × U10N + b
            where m = 0.017 m⁻¹s, b = -0.005
        "lake_low": Sheltered lake regime.
            C_d = 1.0e-3 for U10 < 5 m/s, Charnock above.

    Source: Fairall et al. (2003). See assumptions_register.md.
    """
```

**Wave spectrum and Stokes drift (`waves.py`):**

```python
def jonswap_parameters(U10: float, fetch: float, depth: float) -> dict:
    """Fetch-limited JONSWAP parameters with shallow-water correction.
    Returns: H_s [m], T_p [s], λ_p [m], f_p [Hz], γ_jonswap [-]."""

def stokes_drift_profile(z: np.ndarray, U10: float, fetch: float,
                         depth: float, method: str = "webb_fox_kemper") -> np.ndarray:
    """
    Vertical profile of Stokes drift [m/s].

    Methods:
        "monochromatic": u_s(z) = u_s(0) exp(2kz). Included as comparator only.
            WARNING: overestimates surface drift, underestimates at depth.
        "webb_fox_kemper": Exponential integral approximation for JONSWAP
            spectra (Webb & Fox-Kemper 2011).
            u_s(z) = u_s0 × exp(2 k_e z) / (1 - 8 k_e z)
    """

def differential_drift(z: np.ndarray, U10: float, fetch: float,
                       depth: float) -> np.ndarray:
    """D'(z) = du_s/dz [1/s]. Used by CL instability calculations."""
```

**Surface current (`currents.py`):**

```python
def surface_velocity_1d(U10: float, depth: float,
                        nu_T_profile: Callable) -> float:
    """
    Surface velocity from 1D vertically-resolved momentum balance [m/s].

    Solves:
        d/dz [ν_T(z) du/dz] = 0  (interior)
        ρ ν_T du/dz|_{z=0} = ρ_a C_D U10²  (surface BC)
        ∫_{-h}^{0} u(z) dz = 0  (closed-basin return flow)

    NOT the oceanic 3% heuristic.
    """
```

**Eddy viscosity (`eddy_viscosity.py`):**

```python
def parabolic_nu_T(z: np.ndarray, u_star: float, depth: float,
                   kappa: float = 0.41) -> np.ndarray:
    """Parabolic eddy viscosity profile [m²/s].
    ν_T(z) = κ u* (z+h)(1 - (z+h)/h)"""

def representative_nu_T(profile: np.ndarray, z_grid: np.ndarray) -> float:
    """Depth-averaged ν_T for use in Ra computation [m²/s]."""
```

**Composite forcing output:**

```python
@dataclass(frozen=True)
class ForcingState:
    """Complete forcing description at a single time. All dimensional."""
    U10: float                          # [m/s]
    u_star: float                       # [m/s]
    U_surface: float                    # [m/s]
    stokes_drift_surface: float         # [m/s]
    stokes_drift_profile: np.ndarray    # [m/s] on z_grid
    differential_drift_profile: np.ndarray  # [1/s] on z_grid
    z_grid: np.ndarray                  # [m] from -depth to 0
    depth: float                        # [m]
    fetch: float                        # [m]
    H_s: float                          # [m]
    T_p: float                          # [s]
    La_t: float                         # [-] turbulent Langmuir number
    nu_T: float                         # [m²/s] representative
    nu_T_profile: np.ndarray            # [m²/s] on z_grid
    Ra: float                           # [-] Rayleigh number
    timestamp: datetime
    drag_method: str                    # provenance
    drift_method: str                   # provenance
```

**Tests (`src/forcing/tests/`):**

- u* increases monotonically with U10
- Closed-basin integral equals zero to machine precision
- Shallow-water drift differs from deep-water drift when h < λ_p/2
- No output saturates over U10 ∈ [1, 15] m/s (saturation audit)

**Exit condition:** All tests pass. Every output varies meaningfully across the input range.

---

### WP-04: Hydrodynamics Module

The core work package. Implement `src/hydro/`.

**4a. Rayleigh number and regime classification (`rayleigh.py`):**

```python
def compute_rayleigh(U_surface: float, D_max: float,
                     depth: float, nu_T: float) -> float:
    """Ra = U D h² / ν_T² [-]. Fully dimensional inputs."""

def classify_regime(Ra: float, R0: float, RcNL: float) -> str:
    """
    Returns:
        "subcritical"    — Ra < R0, no LC
        "near_onset"     — R0 ≤ Ra < 1.5 × RcNL
        "moderate"       — 1.5 × RcNL ≤ Ra < 5 × RcNL
        "supercritical"  — Ra ≥ 5 × RcNL
    """

def unstable_band(Ra: float, neutral_curve: Callable,
                  l_array: np.ndarray) -> tuple[float, float]:
    """[l_min, l_max] where Ra > R̄(l). Returns (nan, nan) if subcritical."""
```

**4b. CL solver suite:**

Follow `langmuir_nonlinear_cl_implementation.md` exactly. This includes:

- `profiles.py` — Polynomial shear/drift profiles (equation 7a,b)
- `robin_bc.py` — Robin boundary conditions (equations 2–3)
- `linear_solver.py` — Small-l expansion to O(l¹⁶) (section 3)
- `nonlinear_solver.py` — Nonlinear expansion (sections 4–5) + Galerkin numerical (section 6)
- `galerkin.py` — Shifted Legendre basis, inner products, trigonometric product rules

Verification tests (write before implementing):

```python
def test_R0_uniform():
    """R0 ≈ 120 for D' = U' = 1."""

def test_critical_values():
    """RcL ≈ 121.068, lcL ≈ 0.150, RcNL ≈ 122.194, lcNL ≈ 0.105
    for D' = U' = 1, γ_s = 0.0001, γ_b = 0."""

def test_kappa_values():
    """κ ≈ (1.425, 1.427, 1.919, 1.944) for four profile pairs."""

def test_supercritical_stability():
    """R̄(l) > R(l) for all l > 0 with Robin BCs."""

def test_aspect_ratio_range():
    """L = 2π/lcNL ∈ [5, 11] with γ_s = 0.06, γ_b = 0.28."""

def test_asymptotic_numeric_agreement():
    """Asymptotic expansion and Galerkin agree to 4 significant figures."""
```

**4c. Scaling-law alternative (`scaling_laws.py`):**

```python
def la_dependent_geometry(La_t: float, depth: float, u_star: float) -> dict:
    """
    Non-CL cell geometry from La-dependent scaling laws.

    Returns:
        downwelling_thickness: ~ depth × La^{1/2}  [m]
        downwelling_velocity_max: ~ u* × La^{-1/3}  [m/s]
        pitch: ~ La^{1/6}  [-]
        estimated_cell_width: derived from above  [m]
    """

def empirical_spacing_wind(U10: float, depth: float) -> dict:
    """Published empirical spacing-wind relationships. Cites sources."""
```

**4d. Coarsening (`coarsening.py`):**

```python
def coarsening_timescale(depth: float, u_star: float) -> float:
    """Time for one Y-junction merger event [s]. O(h/u*)."""

def count_coarsening_events(time_available: float, tau: float) -> int:
    """Discrete count of mergers since onset."""

def coarsened_width(initial_width: float, n_events: int,
                    depth: float, max_aspect_ratio: float = 12.0) -> float:
    """Width after n mergers [m]. Capped at max_aspect_ratio × depth.
    Cap emits a warning, does not silently clip."""

def disruption_check(forcing_history: list, lookback_hours: float = 3.0) -> dict:
    """Check for wind direction change >45°, wind drop below 2 m/s,
    or rapid speed increase. Returns disruption flags and reset time."""
```

**Tests:**

- All Hayes & Phillips verification targets pass
- Asymptotic and Galerkin solutions agree
- Scaling laws produce monotonic, unsaturated outputs
- Coarsening produces discrete doublings
- No function saturates across the relevant input range

**Exit condition:** All verification tests pass. Both the CL solver and scaling laws produce spacing predictions spanning a physically plausible range.

---

### ===== DECISION GATE 1 =====

**Stop.** Present WP-01 through WP-04 results to the human. Questions requiring judgement:

1. Does the CL solver reproduce the published results?
2. Does the forcing-to-Ra mapping produce sensible regimes for Lough Neagh (h ≈ 9 m, fetch ≈ 15 km)?
3. Is there a credible path from the hydro outputs to the observed spacing range (40–160 m)?
4. Does the scaling-law candidate look competitive or clearly inferior to CL?

**Do not proceed until this gate is passed.**

---

### WP-05: Prediction Module + LC Enhancement Calculator

**5a. Two candidate pipelines:**

`candidate_cl.py`: forcing → Ra → nonlinear CL solver → coarsening → dimensional spacing → diagnostics

`candidate_scaling.py`: forcing → La → scaling laws → coarsening → dimensional spacing → diagnostics

Both share the forcing module, coarsening logic, and evaluation metrics. They differ only in how they determine initial cell width.

**5b. Visibility diagnostic (`visibility.py`):**

A function, not a module. Takes convergence velocity and wind speed, returns a diagnostic flag and confidence estimate for whether the pattern is satellite-detectable. This is attached to predictions — it does not modify the spacing.

```python
def is_pattern_visible(convergence_v: float, U10: float,
                       tracer_accumulation_time: float,
                       pattern_lifetime: float) -> dict:
    """
    Diagnostic: is this LC pattern detectable from satellite?

    Returns:
        visible: bool
        confidence: float (0-1)
        limiting_factor: str ("convergence_too_weak" / "wind_obscuring" /
                              "insufficient_accumulation_time" / "visible")
    """
```

**5c. LC Enhancement Calculator (`lc_enhancement.py`):**

This module computes how much Langmuir circulation enhances cyanobacterial growing conditions relative to a static water column. It does not model bloom dynamics, scum aggregation, or colony-level processes. It computes physical enhancement ratios.

**The triple rule:** Every biological metric is a named triple `(value_LC, value_static, enhancement_ratio)`. No function in this module may return a single scalar representing a biological quantity.

```python
@dataclass
class EnhancementTriple:
    """One enhancement metric. Always a comparison."""
    lc_value: float
    static_value: float
    ratio: float            # lc_value / static_value; >1 means LC enhances
    units: str
    name: str
    interpretation: str     # one-sentence explanation
```

**Light exposure enhancement:**

```python
def light_enhancement(
    cell_depth: float,           # [m]
    w_circulation: float,        # [m/s] characteristic vertical velocity
    surface_irradiance: float,   # [W/m²] PAR
    K_d: float,                  # [1/m] attenuation coefficient
    photoinhibition_factor: float = 0.5,  # [-] reduction at surface (literature: 0.3–0.7)
) -> EnhancementTriple:
    """
    LC case: Cells cycle through the photic zone. Time-averaged irradiance:
        I_lc = I_0 × (1 - exp(-K_d × h)) / (K_d × h)
    This avoids photoinhibition at the surface and light starvation at depth.
    Cycling period: T = 2h / w_circulation.

    Static case: Buoyant cells pool at the surface under photoinhibiting irradiance:
        I_static = I_0 × photoinhibition_factor

    Enhancement ratio = I_lc / I_static.
    """
```

**Nutrient upwelling enhancement:**

```python
def nutrient_enhancement(
    nu_T_lc: float,              # [m²/s] eddy viscosity with LC
    nu_T_static: float,          # [m²/s] background mixing without LC
    nutrient_gradient: float,    # [µg/L per m] vertical gradient
) -> EnhancementTriple:
    """
    LC case: Vertical flux = ν_T_lc × nutrient_gradient [µg/(m²·s)]
        (eddy diffusion + advective upwelling from organised circulation)

    Static case: Vertical flux = ν_T_static × nutrient_gradient [µg/(m²·s)]
        (molecular diffusion + weak wind mixing; ν_T_static ≈ 1e-5 to 1e-4 m²/s)

    Enhancement ratio = ν_T_lc / ν_T_static.
    Can be very large (10×–1000×) because LC mixing >> molecular diffusion.
    """
```

**Temperature distribution enhancement:**

```python
def temperature_enhancement(
    nu_T_lc: float,              # [m²/s]
    nu_T_static: float,          # [m²/s]
    surface_temperature: float,  # [°C]
    bottom_temperature: float,   # [°C]
    T_optimum: float = 25.0,     # [°C] Microcystis growth optimum
) -> EnhancementTriple:
    """
    LC case: Mixing homogenises temperature.
        T_lc ≈ (T_surface + T_bottom) / 2

    Static case: Stratified. Buoyant cells at surface experience:
        T_static = T_surface

    Enhancement = proximity of T_lc to T_optimum / proximity of T_static to T_optimum.
    Context-dependent: if surface is already near optimum, mixing may cool it away (ratio < 1).
    """
```

**Composite development index:**

```python
def lc_development_index(
    light: EnhancementTriple,
    nutrients: EnhancementTriple,
    temperature: EnhancementTriple,
    regime: str,
    pattern_lifetime: float,     # [s]
    reference_timescale: float = 3600.0,  # [s] 1 hour reference
) -> dict:
    """
    Composite index of how favourable LC conditions are for
    cyanobacterial growth relative to a static water column.

    Components:
        light.ratio, nutrients.ratio, temperature.ratio
        persistence_factor = min(pattern_lifetime / reference_timescale, 1.0)
        regime_factor: subcritical→0, near_onset→0.3, moderate→0.7, supercritical→1.0

    Composite:
        index = regime_factor × persistence_factor ×
                geometric_mean(light.ratio, nutrients.ratio, temperature.ratio)

    Geometric mean because all three must be favourable for growth
    (avoids one very large ratio dominating).

    Returns dict with:
        - All three EnhancementTriples
        - persistence_factor, regime_factor
        - development_index (composite)
        - components_consistent: bool (all ratios > 1 or all < 1)
        - interpretation: str
    """
```

**What "static column" means:**

The static reference is not "no wind." It is the same depth and surface forcing but without the organised overturning circulation. Physically: weak turbulent mixing (ν_T ~ 10⁻⁵ to 10⁻⁴ m²/s), buoyant cells at the surface, nutrients depleting in the surface layer, thermal stratification developing.

**5d. Pipeline entry point (`pipeline.py`):**

```python
def analyse_case(
    weather_data: pd.DataFrame,
    observation_time: datetime,
    depth: float,
    fetch: float,
    candidate: str = "cl",      # "cl" or "scaling"
    environmental: dict = None,  # irradiance, K_d, temperatures, nutrients
) -> dict:
    """
    Full analysis for a single observation case.

    Returns:
        predicted_spacing [m]
        regime classification
        Ra and constituent terms
        forcing summary
        coarsening details
        visibility diagnostic
        lc_enhancement (all three triples + composite index)
        explanation (human-readable sentence)
        all intermediate quantities (for auditing)
    """
```

**Tests:**

- Attractor test: neither candidate clusters in <20% of observed range
- Dynamic range: predictions span ≥60% of observed range
- Sensitivity: ±10% U10 perturbation produces detectable output change
- Enhancement triple rule: every bio output is a 3-tuple (automated check)
- Enhancement ratios > 1 when LC is present, = 1 when subcritical
- Subcritical forcing → no-LC prediction
- Strong wind → moderate spacing; low wind in shallow water → wide spacing

**Exit condition:** Both candidates produce spacing predictions and enhancement indices on the full observation set.

---

### ===== DECISION GATE 2 =====

**Stop.** Present WP-05 results to the human. Questions requiring judgement:

1. Does the CL candidate outperform the scaling candidate? (If not, CL is demoted to comparator.)
2. Do both candidates beat the baselines?
3. Does the attractor test pass for both candidates?
4. Do the LC enhancement ratios look physically sensible?
5. Does the development index correlate with observed bloom conditions?

**Do not proceed until this gate is passed.**

---

### WP-06: Formal Comparison and Production

**6a. Run comparison:**

Both candidates + all baselines on all benchmark subsets. Compute all metrics from WP-02.

**6b. Produce comparison outputs:**

```
outputs/comparison/
├── metrics_table.csv                              # All metrics × candidates × subsets
├── predicted_vs_observed_{candidate}.png           # Scatter plots
├── spacing_vs_wind_{candidate}.png                 # Key diagnostic
├── dynamic_range_comparison.png
├── tail_coverage_comparison.png
├── enhancement_index_timeseries_{case_id}.png      # Representative cases
├── case_diagnostics/
│   └── case_{id}.json                             # Full diagnostic per case
└── attractor_diagnostic_{candidate}.png
```

**6c. Select production candidate.** Document rationale in `assumptions_register.md`.

**6d. Clean up winning candidate into production pipeline.**

The final system takes weather data and returns:

- Spacing prediction with confidence
- Regime classification
- Full forcing and hydrodynamic diagnostics
- LC enhancement triples (light, nutrients, temperature)
- Composite development index
- Human-readable explanation

**Exit condition:** Production pipeline runs on the full dataset. Results match the comparison numbers exactly (sanity check).

---

## 4. Cross-Cutting Requirements

### Saturation Audit

Runs after every module is completed. Sweeps every function over its physically plausible input range. Flags any output varying <5% over >50% of the input range. This is the primary structural health check. See AGENTS.md §2 and §10 for details.

### Dimensional Consistency

Every function documents units for all inputs and outputs. Automated test parses docstrings and fails if any return value lacks units. See AGENTS.md §5.1.

### Living Documents

**`assumptions_register.md`**: Every modelling assumption, parameter value, and architectural decision with source, falsification criterion, and status. Updated every work package.

**`failure_log.md`**: Everything that didn't work, with root cause and resolution. An empty failure log is suspicious.

### Anti-Patterns Enforced by AGENTS.md

- No hidden bounding (sigmoid, tanh, clip) without physical justification
- No compensatory memory wrappers (exponential smoothing to fix dynamics)
- No chaining untested components
- No tuning to aggregate error without physical rationale
- No implicit resolution of ambiguity between spatial scales
- No normalisation that discards amplitude information
- Raw values always preserved alongside any normalisation

---

## 5. Candidate Comparison Logic

Two physics candidates plus three baselines compete:

| Candidate | How it determines cell width | Advantages | Risks |
|---|---|---|---|
| **CL (nonlinear)** | Hayes & Phillips 2017 nonlinear CL solver with Robin BCs | Strongest physical basis; predicts aspect ratios 5–11 | Mapping from instability scale to observed spacing may not be direct |
| **Scaling laws** | La-dependent geometry (downwelling ~ La^{1/2}, velocity ~ La^{-1/3}) | Simpler; no PDE solver; published LES support | Limited explanatory depth; no boundary condition physics |
| Baseline: constant | Mean observed spacing | Floor | No physics |
| Baseline: linear-in-wind | OLS on U10 | Floor | No physics |
| Baseline: depth-scaled | OLS on U10 + depth | Floor | Minimal physics |

The CL candidate earns its place by outperforming the scaling candidate. If it does not, it is demoted to comparator status. Either way, the winning candidate must beat all three baselines on RMSE, dynamic range, and tail coverage.

---

## 6. Agent Instruction Summary

When handing this plan to the agentic coder:

> You are implementing a clean-room rebuild of a Langmuir circulation analysis system. You have access to:
>
> 1. This engineering plan (follow work packages in order)
> 2. `AGENTS.md` — coding governance (read before writing any code)
> 3. `hayes_phillips_2017.md` — primary physics reference for the CL solver
> 4. `langmuir_nonlinear_cl_implementation.md` — detailed CL solver implementation spec
> 5. Background physics formulations — candidate inputs, not mandatory
> 6. Observation datasets and weather cache — reusable data infrastructure
> 7. The old model codebase — FOR BENCHMARK COMPARISON ONLY, never copy code
>
> **Critical rules:**
> - Read AGENTS.md before writing any code.
> - Execute work packages in order. Stop at every decision gate.
> - Write documentation before code for WP-01.
> - Write tests before implementations for WP-02 through WP-06.
> - Every function must have units in its docstring.
> - Run the saturation audit after every module is complete.
> - Log every assumption, parameter, and failure in the living documents.
> - Every biological output must be a triple: (LC_value, static_value, ratio).
> - If a verification test fails, stop and report. Do not adjust the test.
> - If you are uncertain about a modelling choice, flag it in the assumptions register with status "Uncertain" and continue.
