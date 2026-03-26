"""
Shared helpers for the WP-05 prediction pipelines.

This module keeps the prediction-layer assumptions explicit:
    - environmental inputs are carried in a frozen dataclass
    - instability-scale to cell-scale conversion is a named function
    - visibility and LC-enhancement diagnostics preserve their raw inputs
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import warnings
from typing import Any, Mapping

from src.forcing import ForcingState
from src.prediction.lc_enhancement import (
    EnhancementTriple,
    lc_development_index,
    light_enhancement,
    nutrient_enhancement,
    temperature_enhancement,
)
from src.prediction.visibility import is_pattern_visible


@dataclass(frozen=True)
class EnvironmentalContext:
    """
    Environmental inputs for the WP-05 diagnostics.

    Fields:
        surface_irradiance:       Surface PAR proxy [W/m²]
        attenuation_coefficient:  Light attenuation coefficient K_d [1/m]
        photoinhibition_factor:   Static-column surface light reduction [-]
        surface_temperature:      Surface temperature [°C]
        bottom_temperature:       Bottom temperature [°C]
        temperature_optimum:      Cyanobacterial growth optimum [°C]
        nutrient_gradient:        Vertical nutrient gradient [µg/L per m]
        tracer_accumulation_time: Time needed to build a visible surface tracer [s]
        static_nu_T:              Background diffusivity without organised LC [m²/s]
    """

    surface_irradiance: float = 400.0
    attenuation_coefficient: float = 0.5
    photoinhibition_factor: float = 0.5
    surface_temperature: float = 20.0
    bottom_temperature: float = 16.0
    temperature_optimum: float = 25.0
    nutrient_gradient: float = 10.0
    tracer_accumulation_time: float = 1800.0
    static_nu_T: float = 1.0e-4

    def __post_init__(self) -> None:
        if self.surface_irradiance < 0.0:
            raise ValueError("surface_irradiance must be non-negative.")
        if self.attenuation_coefficient < 0.0:
            raise ValueError("attenuation_coefficient must be non-negative.")
        if not (0.0 < self.photoinhibition_factor <= 1.0):
            raise ValueError("photoinhibition_factor must lie in (0, 1].")
        if self.tracer_accumulation_time <= 0.0:
            raise ValueError("tracer_accumulation_time must be positive.")
        if self.static_nu_T < 0.0:
            raise ValueError("static_nu_T must be non-negative.")


def build_environmental_context(
    environmental: Mapping[str, Any] | None = None,
) -> EnvironmentalContext:
    """
    Create an EnvironmentalContext from optional overrides.

    Parameters:
        environmental: Optional mapping with any EnvironmentalContext field names

    Returns:
        EnvironmentalContext with defaults filled explicitly
    """
    if environmental is None:
        return EnvironmentalContext()

    defaults = EnvironmentalContext()
    values = {
        "surface_irradiance": float(
            environmental.get("surface_irradiance", defaults.surface_irradiance)
        ),
        "attenuation_coefficient": float(
            environmental.get(
                "attenuation_coefficient",
                environmental.get("K_d", defaults.attenuation_coefficient),
            )
        ),
        "photoinhibition_factor": float(
            environmental.get(
                "photoinhibition_factor", defaults.photoinhibition_factor
            )
        ),
        "surface_temperature": float(
            environmental.get("surface_temperature", defaults.surface_temperature)
        ),
        "bottom_temperature": float(
            environmental.get("bottom_temperature", defaults.bottom_temperature)
        ),
        "temperature_optimum": float(
            environmental.get("temperature_optimum", defaults.temperature_optimum)
        ),
        "nutrient_gradient": float(
            environmental.get("nutrient_gradient", defaults.nutrient_gradient)
        ),
        "tracer_accumulation_time": float(
            environmental.get(
                "tracer_accumulation_time", defaults.tracer_accumulation_time
            )
        ),
        "static_nu_T": float(
            environmental.get("static_nu_T", defaults.static_nu_T)
        ),
    }
    return EnvironmentalContext(**values)


def instability_scale_to_cell_width(l_c: float, depth: float) -> float:
    """
    Convert nondimensional instability scale to dimensional cell scale [m].

    Input scale:  instability scale from the CL solver (`l_c`) [-]
    Output scale: cell-pair width before coarsening [m]

    Conversion:
        width = 2πh / l_c

    Parameters:
        l_c:   Critical spanwise wavenumber from the CL solver [-]
        depth: Water depth h [m]

    Returns:
        Cell-pair width [m]
    """
    if l_c <= 0.0 or not math.isfinite(l_c):
        raise ValueError(f"l_c must be positive and finite, got {l_c}.")
    if depth <= 0.0:
        raise ValueError(f"depth must be positive, got {depth:.6g} m")
    return float(2.0 * math.pi * depth / l_c)


def initial_cell_width_to_visible_hierarchy_width(
    initial_cell_width: float,
    n_events: int,
    max_visible_mergers: int = 3,
) -> dict:
    """
    Convert initial cell width to a bounded visible merger-hierarchy width [m].

    Input scale:
        initial_cell_width: Instability-selected cell width before late coarsening [m]

    Output scale:
        visible_hierarchy_width: Observation-scale hierarchy width [m]

    Conversion:
        visible_n_events = min(n_events, max_visible_mergers)
        visible_hierarchy_width = initial_cell_width × 2^visible_n_events

    Assumptions:
        - Y-junction mergers double the characteristic spacing.
        - The observation layer carries only the first 0–3 merger levels as a
          plausible visible hierarchy unless stated otherwise.
    """
    if initial_cell_width <= 0.0 or not math.isfinite(initial_cell_width):
        raise ValueError(
            "initial_cell_width must be positive and finite, "
            f"got {initial_cell_width:.6g} m"
        )
    if n_events < 0:
        raise ValueError(f"n_events must be non-negative, got {n_events}")
    if max_visible_mergers < 0:
        raise ValueError(
            f"max_visible_mergers must be non-negative, got {max_visible_mergers}"
        )

    visible_n_events = min(int(n_events), int(max_visible_mergers))
    visible_hierarchy_width = float(initial_cell_width * (2.0 ** visible_n_events))
    return {
        "visible_hierarchy_width_m": visible_hierarchy_width,
        "visible_n_events": int(visible_n_events),
        "mechanical_n_events": int(n_events),
        "max_visible_mergers": int(max_visible_mergers),
    }


def cell_width_to_visible_spacing(
    cell_width: float,
    raw_hierarchy_width: float,
    depth: float,
    visible_spacing_multiplier: float = 1.0,
    max_visible_aspect_ratio: float = 300.0,
) -> dict:
    """
    Convert the mature cell scale to expected visible spacing [m].

    Input scales:
        cell_width:            Mechanically active mature cell width [m]
        raw_hierarchy_width:   Selected observation-scale hierarchy width [m]
        depth:                 Water depth h [m]

    Output scale:
        visible_spacing: Expected spacing between satellite-visible streaks [m]

    Conversion:
        visible_base = raw_hierarchy_width
        visible_spacing = visible_base × visible_spacing_multiplier

    Notes:
        The 12 × h cap applies to the active LC cell width, not necessarily to
        the observed tracer-band spacing. The visible-spacing proxy is therefore
        driven by an explicitly selected observation-scale hierarchy rather than
        by the mechanically active width itself. The observation-scale safeguard
        remains intentionally much looser than the cell-scale cap.
    """
    if cell_width <= 0.0 or not math.isfinite(cell_width):
        raise ValueError(
            f"cell_width must be positive and finite, got {cell_width:.6g} m"
        )
    if raw_hierarchy_width <= 0.0 or not math.isfinite(raw_hierarchy_width):
        raise ValueError(
            "raw_hierarchy_width must be positive and finite, "
            f"got {raw_hierarchy_width:.6g} m"
        )
    if depth <= 0.0:
        raise ValueError(f"depth must be positive, got {depth:.6g} m")
    if visible_spacing_multiplier <= 0.0 or not math.isfinite(visible_spacing_multiplier):
        raise ValueError(
            "visible_spacing_multiplier must be positive and finite, "
            f"got {visible_spacing_multiplier:.6g}"
        )
    if max_visible_aspect_ratio <= 0.0 or not math.isfinite(max_visible_aspect_ratio):
        raise ValueError(
            "max_visible_aspect_ratio must be positive and finite, "
            f"got {max_visible_aspect_ratio:.6g}"
        )

    visible_base = float(raw_hierarchy_width)
    visible_spacing = float(visible_base * visible_spacing_multiplier)
    visible_cap = float(max_visible_aspect_ratio * depth)
    visible_cap_binding = visible_spacing > visible_cap + 1.0e-12
    if visible_cap_binding:
        warnings.warn(
            f"Visible spacing {visible_spacing:.1f} m exceeds observation cap "
            f"{visible_cap:.1f} m. Capping at aspect ratio "
            f"{max_visible_aspect_ratio:.1f}.",
            stacklevel=2,
        )
        visible_spacing = visible_cap

    return {
        "visible_spacing_m": visible_spacing,
        "visible_base_m": visible_base,
        "cell_width_m": float(cell_width),
        "visible_spacing_multiplier": float(visible_spacing_multiplier),
        "visible_cap_m": visible_cap,
        "visible_cap_binding": bool(visible_cap_binding),
    }


def surface_convergence_velocity(
    U_surface: float,
    depth: float,
    cell_width: float,
) -> float:
    """
    Estimate surface convergence velocity from horizontal shear continuity [m/s].

    Parameters:
        U_surface: Wind-driven surface current scale [m/s]
        depth:     Water depth h [m]
        cell_width: Cell-pair width [m]

    Returns:
        Estimated surface convergence velocity [m/s]
    """
    if depth <= 0.0:
        raise ValueError(f"depth must be positive, got {depth:.6g} m")
    if cell_width <= 0.0 or not math.isfinite(cell_width):
        return 0.0
    return float(abs(U_surface) * depth / cell_width)


def _neutral_triple(
    value: float,
    *,
    units: str,
    name: str,
    interpretation: str,
) -> EnhancementTriple:
    """Construct a 1:1 LC/static comparison when no LC is present."""
    return EnhancementTriple(
        lc_value=float(value),
        static_value=float(value),
        ratio=1.0,
        units=units,
        name=name,
        interpretation=interpretation,
    )


def build_lc_enhancement(
    forcing: ForcingState,
    regime: str,
    pattern_lifetime: float,
    environmental: EnvironmentalContext,
    circulation_velocity: float,
) -> dict:
    """
    Build all LC-enhancement diagnostics for a prediction result.

    Parameters:
        forcing:              Current forcing state
        regime:               Hydrodynamic regime label [-]
        pattern_lifetime:     Time available for an organised pattern [s]
        environmental:        Environmental context dataclass
        circulation_velocity: Characteristic overturning velocity [m/s]

    Returns:
        Dict from `lc_development_index`, plus the raw environmental inputs.
    """
    if regime == "subcritical":
        light = _neutral_triple(
            value=environmental.surface_irradiance * environmental.photoinhibition_factor,
            units="W/m²",
            name="light_exposure",
            interpretation="No organised LC is present, so light exposure equals the static reference.",
        )
        nutrients = _neutral_triple(
            value=environmental.static_nu_T * environmental.nutrient_gradient,
            units="µg/(L·s)",
            name="nutrient_flux",
            interpretation="No organised LC is present, so nutrient transport equals the static reference.",
        )
        static_temperature_proximity = 1.0 / (
            1.0
            + abs(environmental.surface_temperature - environmental.temperature_optimum)
        )
        temperature = _neutral_triple(
            value=static_temperature_proximity,
            units="proximity_to_optimum",
            name="temperature_condition",
            interpretation="No organised LC is present, so the temperature condition equals the static reference.",
        )
    else:
        light = light_enhancement(
            cell_depth=forcing.depth,
            w_circulation=max(circulation_velocity, 1.0e-9),
            surface_irradiance=environmental.surface_irradiance,
            K_d=environmental.attenuation_coefficient,
            photoinhibition_factor=environmental.photoinhibition_factor,
        )
        nutrients = nutrient_enhancement(
            nu_T_lc=forcing.nu_T,
            nu_T_static=environmental.static_nu_T,
            nutrient_gradient=environmental.nutrient_gradient,
        )
        temperature = temperature_enhancement(
            nu_T_lc=forcing.nu_T,
            nu_T_static=environmental.static_nu_T,
            surface_temperature=environmental.surface_temperature,
            bottom_temperature=environmental.bottom_temperature,
            T_optimum=environmental.temperature_optimum,
        )

    summary = lc_development_index(
        light=light,
        nutrients=nutrients,
        temperature=temperature,
        regime=regime,
        pattern_lifetime=pattern_lifetime,
    )
    summary["environmental_context"] = environmental
    return summary


def build_visibility_diagnostic(
    forcing: ForcingState,
    regime: str,
    predicted_spacing: float,
    pattern_lifetime: float,
    environmental: EnvironmentalContext,
) -> dict:
    """
    Build the satellite-visibility diagnostic with preserved raw inputs.

    Parameters:
        forcing:            Current forcing state
        regime:             Hydrodynamic regime label [-]
        predicted_spacing:  Predicted cell-pair spacing [m]
        pattern_lifetime:   Time available for an organised pattern [s]
        environmental:      Environmental context dataclass

    Returns:
        Visibility diagnostic dict with raw convergence input included
    """
    if regime == "subcritical" or not math.isfinite(predicted_spacing):
        return {
            "visible": False,
            "confidence": 0.0,
            "limiting_factor": "convergence_too_weak",
            "convergence_velocity": 0.0,
            "tracer_accumulation_time": float(environmental.tracer_accumulation_time),
        }

    convergence_velocity = surface_convergence_velocity(
        U_surface=forcing.U_surface,
        depth=forcing.depth,
        cell_width=predicted_spacing,
    )
    diagnostic = is_pattern_visible(
        convergence_v=convergence_velocity,
        U10=forcing.U10,
        tracer_accumulation_time=environmental.tracer_accumulation_time,
        pattern_lifetime=pattern_lifetime,
    )
    diagnostic["convergence_velocity"] = float(convergence_velocity)
    diagnostic["tracer_accumulation_time"] = float(
        environmental.tracer_accumulation_time
    )
    return diagnostic
