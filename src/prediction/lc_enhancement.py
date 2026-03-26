"""
LC enhancement calculator for cyanobacterial growth conditions.

This module quantifies how Langmuir circulation changes the physical growth
environment relative to a static water column. It does not predict biomass or
bloom state. Every biological output is returned as an explicit LC/static
comparison triple.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class EnhancementTriple:
    """
    One enhancement metric comparing LC and static conditions.

    Fields:
        lc_value:         Value in the LC case [units given by `units`]
        static_value:     Value in the static-column case [same units]
        ratio:            lc_value / static_value [-]
        units:            Units string for both values
        name:             Metric name
        interpretation:   One-sentence explanation
    """

    lc_value: float
    static_value: float
    ratio: float
    units: str
    name: str
    interpretation: str


def _build_triple(
    *,
    lc_value: float,
    static_value: float,
    units: str,
    name: str,
    interpretation: str,
) -> EnhancementTriple:
    """Construct a validated enhancement triple."""
    if static_value == 0.0:
        ratio = 1.0 if lc_value == 0.0 else float("inf")
    else:
        ratio = float(lc_value / static_value)
    return EnhancementTriple(
        lc_value=float(lc_value),
        static_value=float(static_value),
        ratio=float(ratio),
        units=units,
        name=name,
        interpretation=interpretation,
    )


def light_enhancement(
    cell_depth: float,
    w_circulation: float,
    surface_irradiance: float,
    K_d: float,
    photoinhibition_factor: float = 0.5,
) -> EnhancementTriple:
    """
    Light-exposure enhancement from LC-driven vertical circulation.

    Parameters:
        cell_depth:             Mixed-layer or cell depth h [m]
        w_circulation:          Characteristic vertical circulation speed [m/s]
        surface_irradiance:     Surface irradiance I0 [W/m²]
        K_d:                    Light attenuation coefficient [1/m]
        photoinhibition_factor: Static-case surface reduction factor [-]

    Returns:
        EnhancementTriple with irradiance values [W/m²].

    Source:
        `docs/clean_room_rebuild_brief.md` WP-05c; supporting discussion in
        `docs/literature_review.md` §4.1.
    """
    if cell_depth <= 0.0:
        raise ValueError(f"cell_depth must be positive, got {cell_depth:.6g} m")
    if w_circulation <= 0.0:
        raise ValueError(
            f"w_circulation must be positive, got {w_circulation:.6g} m/s"
        )
    if surface_irradiance < 0.0:
        raise ValueError(
            "surface_irradiance must be non-negative, "
            f"got {surface_irradiance:.6g} W/m²"
        )
    if K_d < 0.0:
        raise ValueError(f"K_d must be non-negative, got {K_d:.6g} 1/m")
    if not (0.0 < photoinhibition_factor <= 1.0):
        raise ValueError(
            "photoinhibition_factor must lie in (0, 1], "
            f"got {photoinhibition_factor:.6g}"
        )

    if K_d == 0.0:
        lc_value = surface_irradiance
    else:
        lc_value = surface_irradiance * (1.0 - math.exp(-K_d * cell_depth)) / (
            K_d * cell_depth
        )
    static_value = surface_irradiance * photoinhibition_factor
    circulation_period = 2.0 * cell_depth / w_circulation

    return _build_triple(
        lc_value=lc_value,
        static_value=static_value,
        units="W/m²",
        name="light_exposure",
        interpretation=(
            "LC circulation period "
            f"{circulation_period:.0f} s mixes cells through the photic zone, "
            "increasing time-averaged irradiance relative to a photoinhibited "
            "surface-pooled static column."
        ),
    )


def nutrient_enhancement(
    nu_T_lc: float,
    nu_T_static: float,
    nutrient_gradient: float,
) -> EnhancementTriple:
    """
    Nutrient-flux enhancement from LC-enhanced mixing.

    Parameters:
        nu_T_lc:           Eddy viscosity with LC [m²/s]
        nu_T_static:       Background eddy viscosity without LC [m²/s]
        nutrient_gradient: Vertical nutrient gradient [µg/L per m]

    Returns:
        EnhancementTriple with vertical nutrient-flux values [µg/(L·s)].

    Source:
        `docs/clean_room_rebuild_brief.md` WP-05c; supporting discussion in
        `docs/literature_review.md` §4.2.
    """
    if nu_T_lc < 0.0:
        raise ValueError(f"nu_T_lc must be non-negative, got {nu_T_lc:.6g} m²/s")
    if nu_T_static < 0.0:
        raise ValueError(
            f"nu_T_static must be non-negative, got {nu_T_static:.6g} m²/s"
        )

    lc_value = nu_T_lc * nutrient_gradient
    static_value = nu_T_static * nutrient_gradient
    return _build_triple(
        lc_value=lc_value,
        static_value=static_value,
        units="µg/(L·s)",
        name="nutrient_flux",
        interpretation=(
            "LC enhances vertical nutrient transport in direct proportion to the "
            "increase in effective eddy diffusivity."
        ),
    )


def temperature_enhancement(
    nu_T_lc: float,
    nu_T_static: float,
    surface_temperature: float,
    bottom_temperature: float,
    T_optimum: float = 25.0,
) -> EnhancementTriple:
    """
    Temperature-condition enhancement from LC homogenisation.

    Parameters:
        nu_T_lc:              Eddy viscosity with LC [m²/s]
        nu_T_static:          Background eddy viscosity without LC [m²/s]
        surface_temperature:  Surface temperature [°C]
        bottom_temperature:   Bottom temperature [°C]
        T_optimum:            Growth optimum [°C]

    Returns:
        EnhancementTriple with temperature-proximity values [arbitrary units].

    Notes:
        The returned `lc_value` and `static_value` are proximity-to-optimum
        diagnostics, not temperatures. Temperatures still enter transparently
        through the homogenised and surface cases.

    Source:
        `docs/clean_room_rebuild_brief.md` WP-05c; supporting discussion in
        `docs/literature_review.md` §4.3.
    """
    if nu_T_lc < 0.0:
        raise ValueError(f"nu_T_lc must be non-negative, got {nu_T_lc:.6g} m²/s")
    if nu_T_static < 0.0:
        raise ValueError(
            f"nu_T_static must be non-negative, got {nu_T_static:.6g} m²/s"
        )

    lc_temperature = 0.5 * (surface_temperature + bottom_temperature)
    static_temperature = surface_temperature

    # Assumption: thermal favourability is represented by inverse distance to
    # the growth optimum, avoiding an arbitrary thermal width parameter.
    lc_value = 1.0 / (1.0 + abs(lc_temperature - T_optimum))
    static_value = 1.0 / (1.0 + abs(static_temperature - T_optimum))

    return _build_triple(
        lc_value=lc_value,
        static_value=static_value,
        units="proximity_to_optimum",
        name="temperature_condition",
        interpretation=(
            "LC mixes the surface and bottom temperatures; enhancement is positive "
            "only when the mixed temperature lies closer to the growth optimum "
            "than the stagnant surface temperature."
        ),
    )


def lc_development_index(
    light: EnhancementTriple,
    nutrients: EnhancementTriple,
    temperature: EnhancementTriple,
    regime: str,
    pattern_lifetime: float,
    reference_timescale: float = 3600.0,
) -> dict:
    """
    Composite favourability index for LC-modified growth conditions.

    Parameters:
        light:                Light-exposure enhancement triple [-]
        nutrients:            Nutrient-flux enhancement triple [-]
        temperature:          Temperature-condition enhancement triple [-]
        regime:               Hydrodynamic regime label
        pattern_lifetime:     Expected LC persistence [s]
        reference_timescale:  Lifetime for full persistence credit [s]

    Returns:
        dict containing the three triples, the persistence and regime factors,
        the composite index, a consistency flag, and a short interpretation.

    Source:
        WP-05c specification in `docs/clean_room_rebuild_brief.md`.
    """
    if pattern_lifetime < 0.0:
        raise ValueError(
            f"pattern_lifetime must be non-negative, got {pattern_lifetime:.6g} s"
        )
    if reference_timescale <= 0.0:
        raise ValueError(
            "reference_timescale must be positive, "
            f"got {reference_timescale:.6g} s"
        )

    regime_factor_map = {
        "subcritical": 0.0,
        "near_onset": 0.3,
        "moderate": 0.7,
        "supercritical": 1.0,
    }
    if regime not in regime_factor_map:
        raise ValueError(f"Unknown regime '{regime}'.")

    ratios = [light.ratio, nutrients.ratio, temperature.ratio]
    if any(ratio <= 0.0 for ratio in ratios):
        raise ValueError("Enhancement ratios must be positive to form a geometric mean.")

    persistence_factor = min(pattern_lifetime / reference_timescale, 1.0)
    regime_factor = regime_factor_map[regime]
    geometric_mean = (ratios[0] * ratios[1] * ratios[2]) ** (1.0 / 3.0)
    development_index = regime_factor * persistence_factor * geometric_mean
    components_consistent = (
        all(ratio > 1.0 for ratio in ratios)
        or all(ratio < 1.0 for ratio in ratios)
    )

    if regime == "subcritical":
        interpretation = (
            "No organised LC is predicted, so the composite development index is zero."
        )
    else:
        interpretation = (
            f"{regime.replace('_', ' ')} LC with persistence factor "
            f"{persistence_factor:.2f} gives development index {development_index:.2f} "
            "relative to the static column."
        )

    return {
        "light": light,
        "nutrients": nutrients,
        "temperature": temperature,
        "persistence_factor": float(persistence_factor),
        "regime_factor": float(regime_factor),
        "development_index": float(development_index),
        "components_consistent": components_consistent,
        "interpretation": interpretation,
    }
