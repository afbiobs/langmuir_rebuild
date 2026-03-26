"""
Scaling-law prediction pipeline for Langmuir-cell spacing.

This candidate uses the forcing-derived turbulent Langmuir number, applies the
explicit La-based geometry closure, then passes the resulting cell width through
the same coarsening and diagnostic layers as the CL candidate.
"""

from __future__ import annotations

from src.forcing import ForcingState
from src.hydro.coarsening import (
    coarsened_width,
    coarsening_diffusivity,
    coarsening_schedule,
)
from src.hydro.rayleigh import classify_regime
from src.hydro.scaling_laws import empirical_spacing_wind, la_dependent_geometry
from src.prediction.common import (
    EnvironmentalContext,
    build_lc_enhancement,
    build_visibility_diagnostic,
    cell_width_to_visible_spacing,
    initial_cell_width_to_visible_hierarchy_width,
)


def analyse_candidate_scaling(
    forcing: ForcingState,
    pattern_lifetime: float,
    environmental: EnvironmentalContext,
    visible_spacing_multiplier: float = 1.0,
    max_visible_mergers: int = 3,
    max_visible_aspect_ratio: float = 300.0,
    max_cell_aspect_ratio: float = 12.0,
) -> dict:
    """
    Analyse a single case with the La-based scaling candidate.

    Parameters:
        forcing:           Current forcing state
        pattern_lifetime:  Time available for organised LC [s]
        environmental:     Environmental context dataclass
        visible_spacing_multiplier: Common observation-scale spacing multiplier [-]
        max_visible_mergers:        Maximum visible merger levels carried [-]
        max_visible_aspect_ratio:   Observation-scale safety cap [width/depth]
        max_cell_aspect_ratio:      Mechanical cell-scale cap [width/depth]

    Returns:
        Audit-friendly dict containing spacing, diagnostics, and intermediates
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

    coarsening_closure = coarsening_diffusivity(
        nu_T_vertical=forcing.nu_T,
        u_star_water=forcing.u_star_water,
        depth=forcing.depth,
    )
    regime = classify_regime(forcing.Ra)
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
            "candidate": "scaling",
            "has_lc": False,
            "predicted_spacing": float("nan"),
            "predicted_spacing_m": float("nan"),
            "regime": regime,
            "Ra": float(forcing.Ra),
            "La_t": float(forcing.La_t),
            "forcing_summary": {
                "U10": float(forcing.U10),
                "u_star_water": float(forcing.u_star_water),
                "U_surface": float(forcing.U_surface),
                "stokes_drift_surface": float(forcing.stokes_drift_surface),
                "nu_T": float(forcing.nu_T),
                "depth": float(forcing.depth),
                "fetch": float(forcing.fetch),
                "drag_method": str(forcing.drag_method),
                "drift_method": str(forcing.drift_method),
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
                "geometry": None,
                "empirical_references": empirical_spacing_wind(
                    U10=forcing.U10,
                    depth=forcing.depth,
                ),
            },
            "explanation": (
                f"Ra = {forcing.Ra:.0f} (subcritical). No organised Langmuir "
                "circulation is predicted, so the scaling-law candidate does "
                "not return a coherent spacing."
            ),
        }

    geometry = la_dependent_geometry(
        La_t=forcing.La_t,
        depth=forcing.depth,
        u_star=forcing.u_star_water,
    )
    base_width = float(geometry["estimated_cell_width"])
    schedule = coarsening_schedule(
        time_available=pattern_lifetime,
        initial_width=base_width,
        horizontal_diffusivity=float(coarsening_closure["diffusivity_m2_s"]),
    )
    n_events = int(schedule["n_events"])
    raw_coarsened_width = float(base_width * (2.0 ** n_events))
    cell_spacing = coarsened_width(
        initial_width=base_width,
        n_events=n_events,
        depth=forcing.depth,
        max_aspect_ratio=max_cell_aspect_ratio,
    )
    cap_binding = raw_coarsened_width > cell_spacing + 1.0e-12
    visible_hierarchy = initial_cell_width_to_visible_hierarchy_width(
        initial_cell_width=base_width,
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
        raw_hierarchy_width=base_width,
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
        circulation_velocity=geometry["downwelling_velocity_max"],
    )

    explanation = (
        f"Ra = {forcing.Ra:.0f} ({regime.replace('_', ' ')}). The scaling-law "
        f"candidate uses La_t = {forcing.La_t:.3f}, giving initial cell width "
        f"{base_width:.1f} m. The merger clock uses lateral diffusivity "
        f"A_H = {coarsening_closure['diffusivity_m2_s']:.4f} m²/s, yielding "
        f"capped mature cell width {cell_spacing:.1f} m, "
        f"and plausible visible spacing range "
        f"{visible_spacing_lower['visible_spacing_m']:.1f}–{predicted_spacing:.1f} m. "
        f"The point proxy uses {visible_hierarchy['visible_n_events']} visible "
        f"merger event(s) out of {n_events} mechanical event(s) over "
        f"{pattern_lifetime / 3600.0:.2f} h."
    )
    if cap_binding:
        explanation += " The mechanical cell-scale cap is binding."
    if visible_hierarchy["visible_n_events"] < n_events:
        explanation += (
            f" The observation layer limits the visible hierarchy to the first "
            f"{visible_hierarchy['max_visible_mergers']} merger levels."
        )
    if visible_spacing["visible_cap_binding"]:
        explanation += " The loose observation-scale cap is also binding."

    return {
        "candidate": "scaling",
        "has_lc": True,
        "predicted_spacing": float(predicted_spacing),
        "predicted_spacing_m": float(predicted_spacing),
        "regime": regime,
        "Ra": float(forcing.Ra),
        "La_t": float(forcing.La_t),
        "forcing_summary": {
            "U10": float(forcing.U10),
            "u_star_water": float(forcing.u_star_water),
            "U_surface": float(forcing.U_surface),
            "stokes_drift_surface": float(forcing.stokes_drift_surface),
            "nu_T": float(forcing.nu_T),
            "depth": float(forcing.depth),
            "fetch": float(forcing.fetch),
            "drag_method": str(forcing.drag_method),
            "drift_method": str(forcing.drift_method),
        },
        "coarsening": {
            "initial_cell_width_m": float(base_width),
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
            "raw_coarsened_width_m": float(raw_coarsened_width),
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
            "geometry": geometry,
            "empirical_references": empirical_spacing_wind(
                U10=forcing.U10,
                depth=forcing.depth,
            ),
        },
        "explanation": explanation,
    }
