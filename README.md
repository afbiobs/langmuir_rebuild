# Langmuir Circulation Analysis — Clean-Room Rebuild

## What this model predicts

This model predicts the **hydrodynamic cell spacing** of Langmuir circulation (LC)
in shallow, wind-forced lakes under specified forcing conditions. From that spacing
it derives secondary diagnostics about whether LC is likely to enhance the accumulation
of positively buoyant cyanobacteria at the water surface.

The primary predicted quantity is the **LC cell width** (cross-wind span of a single
circulation cell, in metres), computed from wind speed, water depth, fetch, and a
wave spectrum. The model also predicts regime (subcritical / near-onset / moderate /
supercritical) and the associated downwelling velocity.

## What this model diagnoses

Given the predicted cell geometry and forcing, the model diagnoses whether conditions
are **favourable for the growth of positively buoyant cyanobacteria** through the
effects of LC-driven vertical circulation dynamics on three pathways:

- **Light exposure**: LC cycles buoyant cells through the photic zone on a timescale
  of O(h/w_d), reducing photoinhibition at the surface while preventing light
  starvation at depth. The time-averaged irradiance experienced by a circulating
  cell differs systematically from that of a cell pooled at the surface.
- **Nutrient availability**: Organised vertical circulation enhances upward transport
  of bottom-sourced nutrients (phosphorus, silica) relative to the static or
  weakly mixed water column. The enhancement factor is proportional to the ratio of
  LC-driven eddy viscosity to background mixing.
- **Temperature distribution**: LC mixes the water column toward a more uniform
  temperature profile. Whether this is favourable depends on the thermal optimum of
  the organism and the ambient stratification; the model computes both the LC and
  static temperature experienced and reports the proximity of each to the optimum.

Each pathway is reported as an **enhancement triple**: (value with LC, value without
LC, ratio). A composite favourability index is derived from the three triples weighted
by the regime (subcritical / near-onset / moderate / supercritical).

A secondary diagnostic for **surface-pattern visibility** is also computed: whether
Langmuir cells are likely to produce satellite-visible windrow patterns. This is
relevant to the spacing-prediction validation workflow but is not the primary purpose.

This is a **physical favourability index**, not a bloom prediction. The model does not
estimate biomass, chlorophyll, toxin concentration, or the spatial distribution of
surface scum.

## What is out of scope

- Bloom initiation, growth, or collapse (biological population dynamics)
- Biomass or chlorophyll estimation
- Toxin or microcystin production
- Spatial distribution or drift of surface scum
- Colony-level buoyancy regulation dynamics (Kromkamp–Walsby model; treated as a
  parameter, not modelled dynamically)
- Prediction beyond the forcing window (no biological memory)
- Hydrodynamics beyond the mixed layer (stratification dynamics, deep circulation)
- Surface accumulation (scum formation) — this is a downstream consequence of the
  growth conditions, not the diagnostic target

## Primary physical variables

| Variable | Symbol | Units | Physical meaning |
|---|---|---|---|
| Rayleigh number | Ra | — | Ratio of LC-driving to viscous-damping forces |
| Turbulent Langmuir number | La_t | — | Ratio of friction velocity to Stokes drift speed |
| LC cell width (instability scale) | L_inst | m | Wavelength selected by the CL instability |
| LC cell width (coarsened) | L_cell | m | Physical cell width after merger events |
| Downwelling velocity | w_d | m/s | Vertical velocity at the convergence zone |

## Secondary observable proxies

| Quantity | Relationship to primary variables |
|---|---|
| Satellite-visible windrow spacing | L_cell + tracer accumulation + visibility filter |
| Windrow spacing measured from imagery | Proxy for L_cell subject to observability bias |

The distinction between **hydrodynamic cell spacing** and **satellite-visible windrow
spacing** is fundamental to this model. They are not the same quantity. The observation
model document (`docs/observation_model.md`) defines this relationship explicitly.

## Governing physics

Langmuir circulation arises from the CL2 (Craik–Leibovich type 2) instability: an
interaction between wind-driven mean shear and the Stokes drift of surface waves that
produces counter-rotating longitudinal vortices. In shallow water, nonlinear effects
shift the critical wavenumber downward by 50–70% relative to linear theory (Hayes &
Phillips 2017), producing aspect ratios (width/depth) of 5–11 consistent with
observation rather than the 2–3 predicted by linear theory.

The model implements both a nonlinear CL solver and a scaling-law alternative as
competing candidates. Their predictions are compared against observations to determine
which physics is dominant.

## Site context

The primary validation site is **Lough Neagh, Northern Ireland**: depth h ≈ 9 m,
dominant fetch ≈ 15 km (SW wind direction), latitude ≈ 54.6°N. The observation dataset
spans multiple shallow lakes; the model is designed to generalise to any shallow lake
with known depth and fetch.

## Governing documents

- `AGENTS.md` — coding rules and anti-patterns (read before writing any code)
- `docs/observation_model.md` — defines what the measurements represent
- `docs/assumptions_register.md` — all assumptions with sources and justification
- `docs/failure_log.md` — what has been tried and failed
