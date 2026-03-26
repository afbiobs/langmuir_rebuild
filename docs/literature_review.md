# Literature Review: Langmuir Circulation in Shallow Water

**Status:** WP-02 deliverable. Questions answered with citations to specific equations
and results. All numerical values that appear in verification tests are sourced here.

---

## 1. What mechanisms control spacing in shallow water LC?

### 1.1 Linear Craik–Leibovich theory

**Governing equations.** Langmuir circulation (LC) arises from the CL2 instability
(Craik & Leibovich 1976): a resonant interaction between wind-driven mean shear U(z)
and the depth-varying Stokes drift D(z). The mean field equations in non-dimensional
form for an unstratified shallow layer (s = 1 shear index, z ∈ [−1, 0]) are:

```
U_τ − U_zz = −G                                          (CL 1a)
u_t − ∇²u = U′ψ_y + J(ψ, u)                              (CL 1b)
(∂_t − ∇²)∇²ψ = Ra D′u_y + J(ψ, ∇²ψ)                    (CL 1c)
```

The single dimensionless parameter controlling onset is the Rayleigh number:

```
Ra = U D h² / ν_T²
```

where U is the surface mean velocity, D is the maximum Stokes drift, h is water
depth, and ν_T is a representative eddy viscosity (Hayes & Phillips 2017, eq. Ra-def).

**Deep-water linear theory.** Leibovich & Paolucci (1981) showed that with Neumann
boundary conditions and uniform profiles (D′ = U′ = 1), the critical wavenumber is
l_c = 0 — linear theory predicts arbitrarily wide cells and cannot select a preferred
spacing. The predicted deep-water aspect ratio of 2–3 arises in practice from finite-
depth effects, not from linear stability theory selecting a preferred wavenumber.

**The role of boundary conditions.** Cox & Leibovich (1993) introduced Robin
(mixed/Cauchy) boundary conditions to account for the extra surface stress produced
when LC perturbs the free surface. These take the form (Hayes & Phillips 2017, eqs. 2–3):

- At z = 0: u_z + γ_s u = 0; ψ_zz + (γ_s/2)ψ_z = 0; ψ = 0
- At z = −1: −u_z + γ_b u = 0; −ψ_zz + γ_b ψ_z = 0; ψ = 0

With γ > 0, the critical wavenumber is nonzero: l_cL ~ γ^{1/4}. Cox & Leibovich
(1993) estimated γ_s ≈ 0.06 and γ_b ≈ 0.28 from physical arguments.

**Small-l asymptotic expansion.** Because l_c = O(γ^{1/4}) ≪ 1, Cox & Leibovich
(1993) and Hayes & Phillips (2016) developed a systematic small-l expansion. At each
order O(l^{2j}), the Rayleigh number is built up as R(l) = R₀ + l² R₂ + l⁴ R₄ + ⋯,
and the neutral curve R(l) has a minimum at l_cL = (γ R₀ / |R*₂|)^{1/4}, where R*₂
is the γ = 0 component of R₂. Hayes & Phillips (2016) extended this to arbitrary
polynomial profiles D′(z) and U′(z) of any degree, and to O(l^P) for any P.

**Shallow water: realistic profiles.** Phillips & Dai (2014) showed that in the
shallow s = 1 case, the drift profile D(z) differs substantially from the classical
Stokes drift, and derived consistent free-surface boundary conditions from first
principles (J. Fluid Mech. 743, 141–169, 2014). Their analysis confirmed that the
CLg-equations (which accommodate any level of shear) contract to the CL-equations
for s = 1, so the standard CL framework applies.

**Benchmark linear results for D′ = U′ = 1, γ_s = 0.0001, γ_b = 0:**

```
R₀ ≈ 120    (l = 0 onset threshold)
R_cL ≈ 121.068
l_cL ≈ 0.150
```

These values appear in Hayes & Phillips (2017) Figure 6 caption and §7. They define
the acceptance test `test_critical_values` in the verification suite.

### 1.2 Nonlinear Craik–Leibovich theory

**Hayes & Phillips (2017)** (Geophys. Astrophys. Fluid Dyn. 111(1), 65–90) extended
the framework to weakly nonlinear steady states. Their key results:

**The subcritical bifurcation.** Nonlinearities are first evident at O(l²) through
the symmetry-breaking term (Hayes & Phillips 2017, eq. 74):

```
∂²u₂/∂z² = ∂u₀/∂T − ∂²u₀/∂Y² + R₀[ψ̃₁ U′ ∂²u₀/∂Y² − ∂ψ̃₁/∂z (∂u₀/∂Y)²]
```

The last term breaks the u₀ ↦ −u₀ symmetry, causing a subcritical bifurcation.

**Critical nonlinear wavenumber and Rayleigh number:**

```
l_cNL = (γ R̃₂ / |R*₂|)^{1/4}          (eq. 66)
R_cNL = R₀ + 2(γ R̃₂ |R*₂|)^{1/2}     (eq. 69)
```

where R̃₂ is the coefficient analogous to R*₂ but arising from the nonlinear terms.

**Benchmark nonlinear results for D′ = U′ = 1, γ_s = 0.0001, γ_b = 0:**

```
R_cNL ≈ 122.194
l_cNL ≈ 0.105
```

**The κ ratio** (linear to nonlinear wavenumber) is independent of γ (eq. 68):

```
κ = l_cL / l_cNL = (R₀ R*₂ / R*₂ R̃₂)^{1/4}
```

For D′ = (1, 1+z, 1, 1+z) and U′ = (1, 1, 1+z, 1+z):

```
κ ≈ (1.425, 1.427, 1.919, 1.944)
```

These four values are the primary accuracy target for the nonlinear solver. They
require an asymptotic expansion to at least O(l^4) to resolve; Hayes & Phillips
carried theirs to O(l^{16}).

**Observed aspect ratios.** With γ_s = 0.06, γ_b = 0.28 (Cox & Leibovich 1993):
l_cNL ∈ [0.57, 1.24] → aspect ratio L = 2π/l_cNL ∈ [5, 11], consistent with
observations of Marmorino et al. (2005) in shallow coastal waters. By contrast
linear theory gives aspect ratios well in excess of observed values.

**Supercritical stability.** Robin BCs ensure R̄(l) > R(l) for all l > 0 (Hayes &
Phillips 2017 Figures 3, 4). This means nonlinearities suppress instability at all
wavenumbers — an important constraint on the solver verification.

**Base flow modification.** Nonlinear interaction generates a spanwise-independent
(k = 0) correction u⁰_e(z) to the base flow (Hayes & Phillips 2017 Figure 6, curve i).
This modifies the effective velocity field structure and must be included in any
downwelling velocity estimate derived from the nonlinear eigenfunctions.

**Galerkin numerical verification.** Hayes & Phillips (2017) §6 describe a shifted
Legendre Galerkin method using J = 13 basis functions and I = 1 wavenumber harmonics.
Results from the asymptotic expansion and Galerkin numerics are indistinguishable over
a region "well in excess of l_cNL". This provides an internal consistency check.

### 1.3 Non-CL alternatives: turbulent Langmuir number scaling

McWilliams, Sullivan & Moeng (1997) (J. Fluid Mech. 334, 1–30) introduced the
turbulent Langmuir number:

```
La_t = (u* / U_s0)^{1/2}
```

where u* is the friction velocity and U_s0 is the surface Stokes drift speed.
LES studies show that La_t controls the ratio of Langmuir-to-shear turbulence
and correlates with cell geometry. Typical values: La_t ≈ 0.3 (Langmuir-dominated)
to La_t > 0.7 (shear-dominated). Below La_t ≈ 0.7 full-depth Langmuir cells tend
to form in shallow water.

**La-dependent geometry scaling.** Empirical scaling laws derived from LES
(background reference, parameterisation literature):

```
downwelling thickness ~ h × La_t^{1/2}
downwelling velocity  ~ u* × La_t^{-1/3}
cell pitch            ~ La_t^{1/6}
```

These provide a non-CL candidate prediction for cell geometry without solving the
CL equations. Their validity in the shallow-water, finite-Ra regime is uncertain.

**Empirical spacing-wind relationships.** Faller & Caponi (1978) and Smith (1992)
report field-measured spacings correlating with wind speed. Smith (1992) gives
typical aspect ratios of 2–3 in deep water; these empirical fits cannot be
extrapolated to shallow water without depth normalisation.

---

## 2. What mechanisms control coherence, persistence, and decay?

### 2.1 Coarsening via Y-junction merging

Thorpe (2004) (Annu. Rev. Fluid Mech. 36, 55–79) describes the dominant coarsening
mechanism observed in field and laboratory LC: Y-junction defects, where three
windrows meet, allow the narrower cell pair to merge with a wider neighbour,
approximately doubling the local spacing. This mechanism:

- Produces a discrete distribution of spacings: L_inst × 2^n, n = 0, 1, 2, 3
- Occurs on a timescale O(h/u*): roughly 15–30 min for typical shallow-lake conditions
- Is bounded by the cap of ~12h (Marmorino et al. 2005) beyond which cells become too
  wide to be mechanically maintained against differential advection

Predicted coarsened spacing after n merger events: L_cell = L_inst × 2^n, capped at
L_cap = 12h. The cap must be reported as a warning, not applied silently.

### 2.2 Disruption mechanisms

LC structures are disrupted and reset by (Thorpe 2004, Marmorino et al. 2005):

- **Wind direction change > ~45°**: the momentum input rotates; existing cells begin
  to decay and new cells form oriented to the new wind. Disruption timescale O(h²/ν_T).
- **Wind speed drop below ~2 m/s**: insufficient drift gradient; cells decay.
  Decay timescale O(h²/ν_T) ≈ 1–3 hours for shallow lakes.
- **Rapid wind speed increase**: temporarily disrupts the existing structure before
  re-establishing at the new forcing level.
- **Stratification change** (thermocline shallowing): limits cell depth and may decouple
  the upper-layer circulation from the full-depth structure.

### 2.3 Adjustment timescales

Two relevant timescales:

- **Diffusive**: τ_diff = h²/ν_T ≈ 1–4 hours (Lough Neagh, h = 9 m, ν_T ~ 10⁻⁴ m²/s)
- **Turbulent/advective**: τ_adv = h/u* ≈ 15–40 min (u* ~ 3–10 mm/s)

Y-junction mergers scale as τ_adv. Full structure decay scales as τ_diff.
These timescales are important for estimating the coarsening state at time of
satellite observation.

---

## 3. Published reduced-order and scaling models

### 3.1 Li & Fox-Kemper (2017) KPP-Langmuir

Li & Fox-Kemper (2017) (J. Phys. Oceanogr. 47, 1551–1571) developed a modified
K-Profile Parameterisation (KPP) that accounts for Langmuir-enhanced mixing. Their
enhancement function depends on La_t and is calibrated to LES. The relevant result
for this project is the enhancement of eddy viscosity by a factor:

```
f_LC = [1 + (C_LC / La_t)^n]^{1/n}
```

with constants from their Table 1. This provides a physically grounded estimate of
ν_T_lc/ν_T_background for the nutrient enhancement diagnostic.

### 3.2 Empirical spacing-wind (Faller & Caponi 1978; Smith 1992)

Published empirical fits for windrow spacing as a function of U10:

- Faller & Caponi (1978): spacing ~ 2–5 × depth for moderate wind, laboratory tanks
- Smith (1992): aspect ratio 2–3 in thermocline-bounded mixed layers

Neither relationship accounts for the nonlinear CL shift or the coarsening dynamics.
They are useful as order-of-magnitude sanity checks, not predictions.

---

## 4. Evidence linking LC to cyanobacterial growth conditions

**Revised framing.** The mechanism is not primarily accumulation of surface scum but
rather the modulation of growth conditions by the vertical circulation.

### 4.1 Light exposure modification

LC drives a periodic vertical circulation with timescale T_circ ~ 2h/w_d. In a
well-mixed (LC-forced) water column, the time-averaged PAR irradiance experienced
by a buoyant cell is:

```
I_lc ≈ I₀ × (1 − exp(−K_d h)) / (K_d h)
```

In the static case, a buoyant cell sits near the surface and experiences surface
irradiance reduced by photoinhibition:

```
I_static ≈ I₀ × f_photoinhibition       (f ~ 0.3–0.7)
```

For typical clear-water lakes (K_d ~ 0.5 m⁻¹) and h ~ 9 m, the depth-averaged
irradiance is significantly greater than the inhibited surface value, so the
enhancement ratio I_lc/I_static > 1. The mixing also reduces photoinhibition stress.

Denman & Gargett (1995) and Guasto et al. (2012) (cited in Hayes & Phillips 2017 §8)
discuss the role of Langmuir-scale structures in plankton transport.

### 4.2 Nutrient upwelling

Gargett et al. (2004) (cited in Hayes & Phillips 2017 §8) demonstrated that full-depth
LC in shallow water vastly enhances sediment mixing and nutrient transport. The
nutrient upward flux scales as:

```
F_nutrients ~ ν_T_lc × ∂c/∂z
```

The LC-enhanced ν_T can exceed background (wind-only) mixing by 10–100× at moderate
wind speeds, creating a substantial upward nutrient flux from the nutrient-rich
hypolimnion or sediments.

### 4.3 Temperature modification

LC homogenises the water column temperature. For a warm surface layer over a cooler
bottom, mixing reduces the surface temperature toward the depth-average. Whether this
is beneficial for cyanobacterial growth depends on the thermal optimum of the dominant
species (Microcystis aeruginosa: T_opt ≈ 25°C) and the ambient surface temperature.

In summer bloom conditions (T_surface > 25°C), LC cooling may actually bring
temperature closer to the optimum, enhancing growth relative to a stagnant surface.

### 4.4 Colony buoyancy and residence time

The Kromkamp & Walsby (1990) model describes dynamic buoyancy regulation in
cyanobacteria: gas vesicle inflation/deflation in response to irradiance and turgor
pressure. In the context of this model, colony buoyancy is treated as a parameter
(input) rather than a dynamic variable, with a plausible range:

- Positively buoyant: rise velocity 0.01–0.5 mm/s (field measurements)
- Neutrally buoyant: colonies follow the LC circulation passively
- Negatively buoyant (gas vesicle collapse): colonies sink

The critical condition for downwelling trapping is w_down > w_rise. If the LC
downwelling velocity exceeds the colony rise velocity, colonies circulate with the
LC and spend time at depth (enhancing light mixing and nutrient exposure). If
w_down < w_rise, buoyant colonies escape to the surface.

---

## 5. What quantities can be robustly estimated from available forcing data?

| Quantity | Estimability | Dominant uncertainty |
|---|---|---|
| Wind stress (τ_w) | High | Drag coefficient: ~15% uncertainty at moderate U10 |
| Friction velocity (u*) | High | Propagates from τ_w; well-constrained |
| Surface current (U_surface) | Moderate | Closed-basin integral constraint; ν_T shape uncertain |
| Wave significant height (H_s) | Moderate | Fetch-limited formula; ±25% |
| Peak wave period (T_p) | Moderate | Fetch-limited; ±20% |
| Surface Stokes drift (U_s0) | Moderate | From wave spectrum; ±30–40% |
| Stokes drift profile D(z) | Low–moderate | Depends on wave spectrum accuracy; shallow water form differs from Stokes formula |
| Differential drift D′(z) | Low–moderate | Derivative of D(z); amplifies wave spectrum errors |
| Rayleigh number (Ra) | Low | Cumulative product of uncertainties above; order-of-magnitude estimate |
| Eddy viscosity (ν_T) | Low | Parabolic approximation; factor 2–3 uncertainty |
| Colony buoyancy | Very low | Order-of-magnitude range; not measurable from forcing data |
| Nutrient gradient | Not estimable | Requires in-situ data; fixed range assumption only |

**Implication for model claims:** Given the low estimability of Ra, the model must
not claim point accuracy in spacing prediction. The physically important output is
**regime classification** (subcritical / near-onset / moderate / supercritical) and
the **range of plausible spacings**, not a precise point estimate.

---

## 6. Key numerical values — verification targets

All values used in `src/hydro/tests/test_literature_values.py`:

| Test | Parameter | Source | Value |
|---|---|---|---|
| test_R0_uniform | R₀ (D′=U′=1) | H&P 2017 §7, eq. Ra-def, R₀ | ≈ 120 |
| test_critical_values | R_cL | H&P 2017 Figure 6 caption | ≈ 121.068 |
| test_critical_values | l_cL | H&P 2017 Figure 6 caption | ≈ 0.150 |
| test_critical_values | R_cNL | H&P 2017 Figure 6 caption | ≈ 122.194 |
| test_critical_values | l_cNL | H&P 2017 Figure 6 caption | ≈ 0.105 |
| test_kappa_values | κ (D′=1, U′=1) | H&P 2017 §7.2 | ≈ 1.425 |
| test_kappa_values | κ (D′=1+z, U′=1) | H&P 2017 §7.2 | ≈ 1.427 |
| test_kappa_values | κ (D′=1, U′=1+z) | H&P 2017 §7.2 | ≈ 1.919 |
| test_kappa_values | κ (D′=1+z, U′=1+z) | H&P 2017 §7.2 | ≈ 1.944 |
| test_aspect_ratio_range | L = 2π/l_cNL | H&P 2017 §8 | ∈ [5, 11] |
| test_supercritical_stability | R̄(l) > R(l) for l > 0 | H&P 2017 Figures 3, 4 | all l > 0 with Robin BCs |

All test_* cases use: D′ = U′ = 1 (uniform profiles), γ_s = 0.0001, γ_b = 0
unless otherwise noted. The aspect ratio range uses γ_s = 0.06, γ_b = 0.28.

---

## References

- Craik, A.D.D. & Leibovich, S. (1976). J. Fluid Mech. 73, 401–426.
- Cox, S.M. & Leibovich, S. (1993). Stud. Appl. Math. 90, 183–214.
- Denman, K.L. & Gargett, A.E. (1995). Limnol. Oceanogr. 40, 781–797.
- Faller, A.J. & Caponi, E.A. (1978). J. Geophys. Res. 83, 3617–3633.
- Gargett, A.E. et al. (2004). Science 306, 1925–1928.
- Guasto, J.S. et al. (2012). Annu. Rev. Fluid Mech. 44, 373–400.
- Hayes, D.T. & Phillips, W.R.C. (2016). Theor. Comput. Fluid Dyn. 30, 339–358.
- Hayes, D.T. & Phillips, W.R.C. (2017). Geophys. Astrophys. Fluid Dyn. 111, 65–90.
- Kromkamp, J. & Walsby, A.E. (1990). J. Plankton Res. 12, 161–183.
- Langmuir, I. (1938). Science 87, 119–123.
- Leibovich, S. (1983). Annu. Rev. Fluid Mech. 15, 391–427.
- Leibovich, S. & Paolucci, S. (1981). J. Fluid Mech. 102, 141–167.
- Li, M. & Fox-Kemper, B. (2017). J. Phys. Oceanogr. 47, 1551–1571.
- Marmorino, G.O. et al. (2005). J. Geophys. Res. 110, C06014.
- McWilliams, J.C., Sullivan, P.P. & Moeng, C.-H. (1997). J. Fluid Mech. 334, 1–30.
- Phillips, W.R.C. & Dai, A. (2014). J. Fluid Mech. 743, 141–169.
- Smith, J.A. (1992). J. Phys. Oceanogr. 22, 1403–1426.
- Smith, J.A. et al. (1987). J. Phys. Oceanogr. 17, 449–461.
- Thorpe, S.A. (2004). Annu. Rev. Fluid Mech. 36, 55–79.
