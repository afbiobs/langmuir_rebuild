# Observation Model: What the Spacing Measurements Represent

## 1. What was measured

The dataset contains 67 observations drawn from satellite imagery of shallow lakes.
Each observation records the **spacing between visible wind-aligned surface streaks**
(windrows) as seen from above. Two measurement variants exist:

- `manual_spacing_m` — spacing measured by human annotation of the dominant visible
  streak pattern. Present in ~60% of records.
- `wiggle_spacing_m` — spacing estimated from a semi-automated spectral/correlation
  method applied to the same imagery. Present in ~50% of records; occasionally
  co-present with manual measurements.

Both variants measure **the distance between adjacent visible surface features** in
the cross-wind direction. They do not directly measure the hydrodynamic cell boundary.

The observed range is approximately **37–518 m** (median ~79 m), with the large majority
of observations falling in 50–200 m.

## 2. The four interpretations and their modelling implications

| Interpretation | What it means | Modelling implication |
|---|---|---|
| **A. Hydrodynamic instability scale** | The wavelength selected by the CL instability mechanism at onset | Predicted directly by the CL solver as `2π/l_c × h` (dimensional) |
| **B. Mature cell width** | Physical cross-wind width of a mature LC cell pair after nonlinear saturation and coarsening | Predicted by CL solver output × coarsening model; typically 1–4× the instability scale depending on age |
| **C. Visible streak spacing** | Spacing between surface accumulation lines detectable in satellite imagery | Requires: mature cell width + tracer accumulation dynamics + observability threshold |
| **D. Mixed-scale composite** | Observer captures the most prominent feature, which may be a coarsened cell, a secondary structure, or a tracer band at a different scale | Requires a distribution over scales; point prediction is inappropriate |

### Current working interpretation: **C with acknowledged ambiguity toward D**

The manual and wiggle measurements capture **whatever streak is most visually prominent
in the image**. In most cases this is likely to be interpretation C (a visible convergence
line marking an LC cell boundary). However, several factors introduce ambiguity:

1. **Coarsening state is unknown.** If the image was acquired hours after LC onset, the
   cells may have undergone one or more Y-junction merger events (Thorpe 2004), doubling
   the apparent spacing relative to the instability scale. Without knowledge of LC age,
   the relationship to the instability scale is uncertain by a factor of 1–4.

2. **Wide-spacing bias at low wind.** Narrow windrows fade first as wind decreases;
   surviving features at low wind tend to be wider, older cells. Annotations therefore
   likely over-represent wide spacings in the low-wind regime.

3. **Tracer dependence.** Windrows are only visible if a tracer (bubbles, foam, surface
   scum, thermal contrast) has accumulated at the convergence zone. A cell may exist
   without being observable. The non-detection problem (see §4) is not solved.

4. **Scale mixing.** In some images, streaks at two or more spacings are simultaneously
   present. The annotation captures one; the other is lost. The `wiggle_spacing_m`
   spectral method may capture a different mode than `manual_spacing_m`.

**Consequence for modelling:** The model must not assume that its predicted hydrodynamic
cell spacing is directly comparable to the observed windrow spacing. The comparison
requires an explicit intermediate step that converts cell spacing to expected visible
spacing. That conversion carries uncertainty and must be documented.

## 3. Measurement uncertainty

### Annotation precision
Manual measurements of streak spacing from satellite imagery carry typical uncertainties
of ±10–30% depending on image resolution, streak clarity, and measurement method. A
reported spacing of 80 m may represent anything from 55–105 m.

### Timing uncertainty
The exact time offset between the satellite overpass and the peak of the LC episode is
generally unknown. LC onset, growth, and coarsening occur over timescales of O(h/u*)
to O(h²/ν_T) — roughly 15 minutes to several hours at Lough Neagh. An observation taken
at the end of a 3-hour wind event may represent a more coarsened state than one taken
at the beginning.

### Scale mixture uncertainty
When `manual_spacing_m` and `wiggle_spacing_m` are both present for the same observation
and agree within 20%, they likely capture the same dominant mode. When they disagree by
more than 50%, scale mixing or method divergence is probable. Both values should be
retained for comparison; neither should be silently preferred.

### Directional measurement bias
Spacing measurements from satellite imagery are inherently two-dimensional projections.
If the satellite viewing geometry and the LC orientation are not perfectly perpendicular,
the measured spacing is `L_true / cos(θ)`, where θ is the angle between the satellite
track and the cross-wind direction. For typical polar-orbit satellites and moderate
wind directions this error is small (<15%) but non-negligible.

## 4. Non-detection

An observation record without any spacing value (`manual_spacing_m` and `wiggle_spacing_m`
both empty) may represent any of the following:

- **Subcritical forcing:** Ra below threshold; no LC present.
- **LC present but invisible:** tracer insufficient for visual detection.
- **LC present and visible but not annotated:** image quality, resolution, or analyst
  choice.
- **Wrong image timing:** LC present at some hours but not at the time of the satellite
  overpass.

**The model must not treat non-detections as evidence that Ra < R_c.** They are
censored observations. If a non-detection record has forcing data, it should be retained
in the dataset and flagged as censored rather than silently excluded.

## 5. Mapping to model outputs

The model predicts the following distinct spacing quantities, which must not be silently
conflated:

| Model quantity | Symbol | What it is |
|---|---|---|
| Instability wavelength | L_inst | `2π/l_cNL × h` [m], from CL solver at onset |
| Coarsened cell width | L_cell | L_inst × 2^n_mergers [m], capped at 12h |
| Expected visible spacing | L_vis | L_cell × visibility_factor(tracer, wind) [m] |

The comparison to observations uses **L_vis**, not L_inst or L_cell directly.
The difference between these quantities is a scientific question, not a calibration
parameter. Converting L_cell to L_vis requires explicit physical assumptions about
tracer accumulation rates and observability thresholds, documented in the assumptions
register.

## 6. Ambiguity carriage policy

Where the observation model is genuinely ambiguous (as it is for the coarsening state),
the downstream models must **carry that ambiguity explicitly**:

- Predictions should report L_cell with an estimated range reflecting 0–3 merger events.
- Comparisons should show sensitivity to the coarsening assumption.
- No single number should be presented as "the predicted spacing" without stating which
  of the quantities in §5 it represents.
