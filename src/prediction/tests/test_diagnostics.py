"""
Tests for the bounded WP-05 visibility and enhancement slice.
"""

from __future__ import annotations

from src.prediction.lc_enhancement import (
    EnhancementTriple,
    lc_development_index,
    light_enhancement,
    nutrient_enhancement,
    temperature_enhancement,
)
from src.prediction.visibility import is_pattern_visible


class TestVisibilityDiagnostic:
    def test_visible_case(self):
        """Moderate wind, adequate convergence and persistence should be visible."""
        result = is_pattern_visible(
            convergence_v=2.0e-4,
            U10=5.0,
            tracer_accumulation_time=1800.0,
            pattern_lifetime=3600.0,
        )
        assert result["visible"] is True
        assert result["limiting_factor"] == "visible"
        assert 0.0 <= result["confidence"] <= 1.0

    def test_convergence_limited_case(self):
        """Weak convergence should be identified as the limiting factor."""
        result = is_pattern_visible(
            convergence_v=1.0e-5,
            U10=5.0,
            tracer_accumulation_time=1800.0,
            pattern_lifetime=3600.0,
        )
        assert result["visible"] is False
        assert result["limiting_factor"] == "convergence_too_weak"

    def test_wind_and_persistence_limits(self):
        """High wind and short lifetime should each be reported explicitly."""
        windy = is_pattern_visible(
            convergence_v=2.0e-4,
            U10=12.0,
            tracer_accumulation_time=1800.0,
            pattern_lifetime=3600.0,
        )
        short = is_pattern_visible(
            convergence_v=2.0e-4,
            U10=5.0,
            tracer_accumulation_time=3600.0,
            pattern_lifetime=900.0,
        )
        assert windy["limiting_factor"] == "wind_obscuring"
        assert short["limiting_factor"] == "insufficient_accumulation_time"


class TestEnhancementTriples:
    def test_light_enhancement_returns_triple_and_ratio_gt_one(self):
        """LC light exposure should exceed the photoinhibited static case."""
        triple = light_enhancement(
            cell_depth=9.0,
            w_circulation=0.003,
            surface_irradiance=400.0,
            K_d=0.05,
            photoinhibition_factor=0.5,
        )
        assert isinstance(triple, EnhancementTriple)
        assert triple.ratio > 1.0

    def test_nutrient_enhancement_tracks_diffusivity_ratio(self):
        """The nutrient enhancement ratio reduces to ν_T,lc / ν_T,static."""
        triple = nutrient_enhancement(
            nu_T_lc=1.0e-3,
            nu_T_static=2.0e-4,
            nutrient_gradient=10.0,
        )
        assert isinstance(triple, EnhancementTriple)
        assert triple.ratio == 5.0

    def test_temperature_enhancement_can_help_or_hurt(self):
        """Mixing can move temperature toward or away from the growth optimum."""
        helpful = temperature_enhancement(
            nu_T_lc=1.0e-3,
            nu_T_static=1.0e-4,
            surface_temperature=30.0,
            bottom_temperature=22.0,
            T_optimum=25.0,
        )
        harmful = temperature_enhancement(
            nu_T_lc=1.0e-3,
            nu_T_static=1.0e-4,
            surface_temperature=25.0,
            bottom_temperature=15.0,
            T_optimum=25.0,
        )
        assert helpful.ratio > 1.0
        assert harmful.ratio < 1.0

    def test_development_index_respects_regime_and_triple_rule(self):
        """The composite index uses triples and collapses to zero when subcritical."""
        light = light_enhancement(9.0, 0.003, 400.0, 0.5)
        nutrients = nutrient_enhancement(1.0e-3, 2.0e-4, 10.0)
        temperature = temperature_enhancement(1.0e-3, 1.0e-4, 30.0, 22.0)

        subcritical = lc_development_index(
            light, nutrients, temperature, regime="subcritical", pattern_lifetime=7200.0
        )
        supercritical = lc_development_index(
            light, nutrients, temperature, regime="supercritical", pattern_lifetime=7200.0
        )

        assert isinstance(subcritical["light"], EnhancementTriple)
        assert isinstance(subcritical["nutrients"], EnhancementTriple)
        assert isinstance(subcritical["temperature"], EnhancementTriple)
        assert subcritical["development_index"] == 0.0
        assert supercritical["development_index"] > 0.0
