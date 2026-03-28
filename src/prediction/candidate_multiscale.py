"""
Module: candidate_multiscale
Layer: prediction
Purpose: Multi-scale prediction pipeline using wave-tied LC initiation.

This candidate replaces the CL onset scale as the coarsening initial condition
with the wave-tied LC scale d_LC = C_wave x lambda_p (Tsai & Lu 2023).  The
existing coarsening, visibility, and enhancement machinery is reused without
modification.

The key physical difference from the CL candidate:
    CL candidate:        L_inst (CL onset, ~196 m at Neagh) -> coarsening
    Multiscale candidate: d_LC   (wave-tied, ~3.7 m at Neagh) -> coarsening

The wave-tied scale coarsens through Y-junction mergers to reach 30-120 m over
typical pattern lifetimes (5-60 min), naturally matching the observed manual and
stream spacing ranges.  The CL onset scale is preserved as a diagnostic for the
full-depth roll organization detected by wiggle spectral methods.

Assumptions: AR-031, AR-032, AR-033 -- see assumptions_register.md
Key references: Tsai & Lu (2023, JFM 969 A30), Xuan & Shen (2025, JFM 1023 A4)
"""

from __future__ import annotations

import math
import warnings

from src.forcing import ForcingState
from src.hydro.coarsening import (
    coarsened_width,
    coarsening_diffusivity,
    coarsening_schedule,
)
from src.hydro.multiscale_structures import (
    compute_multiscale_result,
    wave_tied_lc_spacing,
)
from src.hydro.nonlinear_solver import (
    compute_kappa,
    critical_nonlinear,
)
from src.hydro.profiles import PolynomialProfile, polynomial_profile
from src.hydro.rayleigh import classify_regime
from src.hydro.robin_bc import RobinBC, derive_robin_bc_from_forcing
from src.prediction.common import (
    EnvironmentalContext,
    build_lc_enhancement,
    build_visibility_diagnostic,
    cell_width_to_visible_spacing,
    initial_cell_width_to_visible_hierarchy_width,
    instability_scale_to_cell_width,
)


def _fit_affine_drift_profile(forcing: ForcingState) -> tuple[PolynomialProfile, dict]:
    """
    Fit the resolved forcing drift profile with an affine polynomial on [-1, 0].

    Reuses the same logic as the CL candidate for consistency.

    Parameters:
        forcing: Current forcing state with dimensional drift profile

    Returns:
        (PolynomialProfile, diagnostics dict)
    """
    import numpy as np

    zeta = forcing.z_grid / forcing.depth
    drift = np.asarray(forcing.differential_drift_profile, dtype=float)
    scale = float(np.max(np.abs(drift)))
    if scale <= 0.0 or not np.isfinite(scale):
        warnings.warn(
            "Differential-drift profile is degenerate; falling back to D'(z)=1+z.",
            stacklevel=2,
        )
        return polynomial_profile([1.0, 1.0], name="affine_drift_fallback"), {
            "normalisation_scale": scale,
            "rmse": float("nan"),
            "source": "fallback",
        }

    normalised = drift / scale
    bottom_ratio = float(normalised[0])
    if bottom_ratio <= 0.0 or not np.isfinite(bottom_ratio):
        warnings.warn(
            "Affine drift fit encountered a non-positive bottom ratio; falling "
            "back to D'(z)=1+z instead of clipping the forcing shape.",
            stacklevel=2,
        )
        return polynomial_profile([1.0, 1.0], name="affine_drift_fallback"), {
            "normalisation_scale": scale,
            "rmse": float("nan"),
            "bottom_ratio": bottom_ratio,
            "source": "fallback_nonpositive_bottom_ratio",
        }
    intercept = 1.0
    slope = 1.0 - bottom_ratio
    fit = intercept + slope * zeta

    return polynomial_profile([intercept, slope], name="affine_drift_fit"), {
        "normalisation_scale": scale,
        "rmse": float(np.sqrt(np.mean((fit - normalised) ** 2))),
        "bottom_ratio": bottom_ratio,
        "source": "forcing_endpoint_fit",
    }


def _default_shear_profile() -> PolynomialProfile:
    """Surface-intensified affine shear profile: U'(z) = 1 + z on [-1, 0]."""
    return polynomial_profile([1.0, 1.0], name="surface_intensified_shear")


def analyse_candidate_multiscale(
    forcing: ForcingState,
    pattern_lifetime: float,
    environmental: EnvironmentalContext,
    bc: RobinBC | None = None,
    D_profile: PolynomialProfile | None = None,
    U_profile: PolynomialProfile | None = None,
    n_pts: int = 2000,
    visible_spacing_multiplier: float = 1.0,
    max_visible_mergers: int = 3,
    max_visible_aspect_ratio: float = 300.0,
    max_cell_aspect_ratio: float = 12.0,
) -> dict:
    """
    Analyse a single case with the multi-scale wave-tied candidate.

    This candidate initialises coarsening from the wave-tied LC scale
    d_LC = C_wave x lambda_p (Tsai & Lu 2023) instead of the CL onset scale.
    The CL onset scale is still computed and reported as a diagnostic for the
    full-depth roll pattern.

    Parameters:
        forcing:           Current forcing state
        pattern_lifetime:  Time available for organised LC [s]
        environmental:     Environmental context dataclass
        bc:                Robin boundary conditions [-]
        D_profile:         Optional override for D'(z) [-]
        U_profile:         Optional override for U'(z) [-]
        n_pts:             Integration resolution for the CL solver
        visible_spacing_multiplier: Observation-scale spacing multiplier [-]
        max_visible_mergers:        Maximum visible merger levels carried [-]
        max_visible_aspect_ratio:   Observation-scale safety cap [width/depth]
        max_cell_aspect_ratio:      Mechanical cell-scale cap [width/depth]

    Returns:
        Audit-friendly dict containing spacing, diagnostics, and intermediates.
        The predicted_spacing is based on the coarsened wave-tied scale.
    """
    if pattern_lifetime < 0.0:
        raise ValueError("pattern_lifetime must be non-negative.")
    if visible_spacing_multiplier <= 0.0:
        raise ValueError("visible_spacing_multiplier must be positive.")
    if max_visible_mergers < 0:
        raise ValueError("max_visible_mergers must be non-negative.")
    if max_visible_aspect_ratio <= 0.0:
        raise ValueError("max_visible_aspect_ratio must be positive.")
    if max_cell_aspect_ratio <= 0.0:
        raise ValueError("max_cell_aspect_ratio must be positive.")

    # --- Robin BCs and profiles (same as CL candidate) ---
    robin_diagnostics = None
    if bc is None:
        bc, robin_diagnostics = derive_robin_bc_from_forcing(forcing)
    else:
        robin_diagnostics = {
            "source": "user_override",
            "gamma_s_raw": float(bc.gamma_s),
            "gamma_b_raw": float(bc.gamma_b),
            "gamma_total_raw": float(bc.gamma),
        }

    coarsening_closure = coarsening_diffusivity(
        nu_T_vertical=forcing.nu_T,
        u_star_water=forcing.u_star_water,
        depth=forcing.depth,
    )
    regime = classify_regime(forcing.Ra)

    drift_fit = None
    if D_profile is None:
        D_profile, drift_fit = _fit_affine_drift_profile(forcing)
    if U_profile is None:
        U_profile = _default_shear_profile()

    # --- Forcing summary (shared across all outcomes) ---
    forcing_summary = {
        "U10": float(forcing.U10),
        "u_star_water": float(forcing.u_star_water),
        "U_surface": float(forcing.U_surface),
        "stokes_drift_surface": float(forcing.stokes_drift_surface),
        "nu_T": float(forcing.nu_T),
        "depth": float(forcing.depth),
        "fetch": float(forcing.fetch),
        "k_p": float(forcing.k_p),
        "lambda_p": float(forcing.lambda_p),
        "wave_steepness": float(forcing.wave_steepness),
        "drag_method": str(forcing.drag_method),
        "drift_method": str(forcing.drift_method),
    }

    # --- Subcritical: no LC ---
    if regime == "subcritical":
        enhancement = build_lc_enhancement(
            forcing=forcing,
            regime=regime,
            pattern_lifetime=0.0,
            environmental=environmental,
            circulation_velocity=0.0,
        )
        visibility = build_visibility_diagnostic(
            forcing=forcing,
            regime=regime,
            predicted_spacing=float("nan"),
            pattern_lifetime=0.0,
            environmental=environmental,
        )
        return {
            "candidate": "multiscale",
            "has_lc": False,
            "predicted_spacing": float("nan"),
            "predicted_spacing_m": float("nan"),
            "regime": regime,
            "Ra": float(forcing.Ra),
            "La_t": float(forcing.La_t),
            "forcing_summary": forcing_summary,
            "multiscale": {
                "wave_tied_initial_m": float("nan"),
                "cl_onset_spacing_m": float("nan"),
                "C_wave": 0.34,
                "lambda_p_m": float(forcing.lambda_p),
            },
            "coarsening": {
                "initial_cell_width_m": float("nan"),
                "tau_coarsening_s": float("nan"),
                "tau_next_coarsening_s": float("nan"),
                "time_used_for_coarsening_s": 0.0,
                "pattern_lifetime_s": 0.0,
                "forcing_nu_T_m2_s": float(coarsening_closure["vertical_nu_T_m2_s"]),
                "coarsening_diffusivity_m2_s": float(
                    coarsening_closure["diffusivity_m2_s"]
                ),
                "lateral_mixing_length_diffusivity_m2_s": float(
                    coarsening_closure["lateral_mixing_length_diffusivity_m2_s"]
                ),
                "coarsening_anisotropy_ratio": float(
                    coarsening_closure["anisotropy_ratio"]
                ),
                "coarsening_diffusivity_method": str(coarsening_closure["method"]),
                "n_events": 0,
                "visible_n_events": 0,
                "max_visible_mergers": int(max_visible_mergers),
                "raw_coarsened_width_m": float("nan"),
                "coarsened_width_m": float("nan"),
                "visible_spacing_lower_m": float("nan"),
                "visible_spacing_upper_m": float("nan"),
                "cap_binding": False,
                "max_cell_aspect_ratio": float(max_cell_aspect_ratio),
                "max_visible_aspect_ratio": float(max_visible_aspect_ratio),
            },
            "visibility": visibility,
            "lc_enhancement": enhancement,
            "intermediate": {
                "D_profile": D_profile,
                "U_profile": U_profile,
                "drift_fit": drift_fit,
                "robin_bc": bc,
                "robin_closure": robin_diagnostics,
            },
            "explanation": (
                f"Ra = {forcing.Ra:.0f} (subcritical). No organised Langmuir "
                "circulation is predicted, so spacing is undefined and the "
                "enhancement ratios reduce to the static reference."
            ),
        }

    # --- Supercritical: compute both scales ---

    # 1. CL onset scale (diagnostic reference for full-depth roll)
    critical = critical_nonlinear(
        D_prime=D_profile,
        U_prime=U_profile,
        bc=bc,
        method="asymptotic",
        n_pts=n_pts,
    )
    cl_onset_width = instability_scale_to_cell_width(critical.l_c, forcing.depth)

    # 2. Wave-tied LC initial scale (Tsai & Lu 2023: d_LC = 0.34 x lambda_p)
    wave_tied_initial = wave_tied_lc_spacing(forcing.lambda_p)

    # 3. Assemble multi-scale result
    multiscale = compute_multiscale_result(
        lambda_p=forcing.lambda_p,
        depth=forcing.depth,
        cl_onset_spacing=cl_onset_width,
    )

    # 4. Coarsening from the wave-tied initial scale (NOT from CL onset)
    schedule = coarsening_schedule(
        time_available=pattern_lifetime,
        initial_width=wave_tied_initial,
        horizontal_diffusivity=float(coarsening_closure["diffusivity_m2_s"]),
    )
    n_events = int(schedule["n_events"])
    raw_coarsened = float(wave_tied_initial * (2.0 ** n_events))
    cell_spacing = coarsened_width(
        initial_width=wave_tied_initial,
        n_events=n_events,
        depth=forcing.depth,
        max_aspect_ratio=max_cell_aspect_ratio,
    )
    cap_binding = raw_coarsened > cell_spacing + 1.0e-12

    # 5. Visible spacing hierarchy (same machinery as CL candidate)
    kappa = compute_kappa(D_profile, U_profile, bc=bc, n_pts=n_pts)
    visible_hierarchy = initial_cell_width_to_visible_hierarchy_width(
        initial_cell_width=wave_tied_initial,
        n_events=n_events,
        max_visible_mergers=max_visible_mergers,
    )

    visibility = build_visibility_diagnostic(
        forcing=forcing,
        regime=regime,
        predicted_spacing=cell_spacing,
        pattern_lifetime=pattern_lifetime,
        environmental=environmental,
    )

    visible_spacing_lower = cell_width_to_visible_spacing(
        cell_width=cell_spacing,
        raw_hierarchy_width=wave_tied_initial,
        depth=forcing.depth,
        visible_spacing_multiplier=visible_spacing_multiplier,
        max_visible_aspect_ratio=max_visible_aspect_ratio,
    )
    visible_spacing = cell_width_to_visible_spacing(
        cell_width=cell_spacing,
        raw_hierarchy_width=visible_hierarchy["visible_hierarchy_width_m"],
        depth=forcing.depth,
        visible_spacing_multiplier=visible_spacing_multiplier,
        max_visible_aspect_ratio=max_visible_aspect_ratio,
    )
    predicted_spacing = float(visible_spacing["visible_spacing_m"])

    enhancement = build_lc_enhancement(
        forcing=forcing,
        regime=regime,
        pattern_lifetime=pattern_lifetime,
        environmental=environmental,
        circulation_velocity=visibility["convergence_velocity"],
    )

    # --- Explanation ---
    explanation = (
        f"Ra = {forcing.Ra:.0f} ({regime.replace('_', ' ')}). The multiscale "
        f"candidate uses d_LC = 0.34 x lambda_p = {wave_tied_initial:.2f} m as "
        f"the wave-tied LC initial scale (Tsai & Lu 2023), instead of the CL "
        f"onset scale L_inst = {cl_onset_width:.1f} m. The merger clock uses "
        f"lateral diffusivity A_H = "
        f"{coarsening_closure['diffusivity_m2_s']:.4f} m^2/s, yielding "
        f"{n_events} merger event(s) in {pattern_lifetime / 3600.0:.2f} h, "
        f"capped mature cell width = {cell_spacing:.1f} m, "
        f"and plausible visible spacing range = "
        f"{visible_spacing_lower['visible_spacing_m']:.1f}"
        f"--{predicted_spacing:.1f} m."
    )
    if cap_binding:
        explanation += " The mechanical cell-scale cap is binding."
    if visible_hierarchy["visible_n_events"] < n_events:
        explanation += (
            f" The observation layer limits the visible hierarchy to the first "
            f"{visible_hierarchy['max_visible_mergers']} merger levels."
        )

    return {
        "candidate": "multiscale",
        "has_lc": True,
        "predicted_spacing": float(predicted_spacing),
        "predicted_spacing_m": float(predicted_spacing),
        "regime": regime,
        "Ra": float(forcing.Ra),
        "La_t": float(forcing.La_t),
        "forcing_summary": forcing_summary,
        "multiscale": {
            "wave_tied_initial_m": float(wave_tied_initial),
            "cl_onset_spacing_m": float(cl_onset_width),
            "C_wave": float(multiscale.C_wave),
            "lambda_p_m": float(forcing.lambda_p),
            "aspect_ratio_wave": float(multiscale.aspect_ratio_wave),
            "aspect_ratio_onset": float(multiscale.aspect_ratio_onset),
        },
        "coarsening": {
            "initial_cell_width_m": float(wave_tied_initial),
            "tau_coarsening_s": float(schedule["tau_initial_s"]),
            "tau_next_coarsening_s": float(schedule["tau_next_s"]),
            "time_used_for_coarsening_s": float(schedule["time_used_s"]),
            "pattern_lifetime_s": float(pattern_lifetime),
            "forcing_nu_T_m2_s": float(coarsening_closure["vertical_nu_T_m2_s"]),
            "coarsening_diffusivity_m2_s": float(
                coarsening_closure["diffusivity_m2_s"]
            ),
            "lateral_mixing_length_diffusivity_m2_s": float(
                coarsening_closure["lateral_mixing_length_diffusivity_m2_s"]
            ),
            "coarsening_anisotropy_ratio": float(
                coarsening_closure["anisotropy_ratio"]
            ),
            "coarsening_diffusivity_method": str(coarsening_closure["method"]),
            "n_events": int(n_events),
            "visible_n_events": int(visible_hierarchy["visible_n_events"]),
            "max_visible_mergers": int(visible_hierarchy["max_visible_mergers"]),
            "raw_coarsened_width_m": float(raw_coarsened),
            "coarsened_width_m": float(cell_spacing),
            "cell_width_m": float(cell_spacing),
            "cap_binding": cap_binding,
            "visible_hierarchy_width_m": float(
                visible_hierarchy["visible_hierarchy_width_m"]
            ),
            "visible_spacing_lower_m": float(
                visible_spacing_lower["visible_spacing_m"]
            ),
            "visible_spacing_upper_m": float(predicted_spacing),
            "visible_spacing_m": float(predicted_spacing),
            "visible_spacing_multiplier": float(
                visible_spacing["visible_spacing_multiplier"]
            ),
            "visible_cap_binding": bool(visible_spacing["visible_cap_binding"]),
            "max_cell_aspect_ratio": float(max_cell_aspect_ratio),
            "max_visible_aspect_ratio": float(max_visible_aspect_ratio),
        },
        "visibility": visibility,
        "lc_enhancement": enhancement,
        "intermediate": {
            "critical_result": critical,
            "kappa": float(kappa),
            "cl_onset_wavenumber": float(critical.l_c),
            "multiscale_result": multiscale,
            "D_profile": D_profile,
            "U_profile": U_profile,
            "drift_fit": drift_fit,
            "robin_bc": bc,
            "robin_closure": robin_diagnostics,
        },
        "explanation": explanation,
    }
