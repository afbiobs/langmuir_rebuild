"""Tests for the multiscale wave-tied prediction candidate."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import math

import numpy as np
import pandas as pd
import pytest

from src.forcing import ForcingState
from src.hydro.multiscale_structures import wave_tied_lc_spacing
from src.prediction.candidate_multiscale import analyse_candidate_multiscale
from src.prediction.common import EnvironmentalContext
from src.prediction.pipeline import analyse_case


def _synthetic_forcing_state(Ra: float, lambda_p: float = 15.7) -> ForcingState:
    """Create a minimal forcing state for branch tests."""
    k_p = 2.0 * math.pi / lambda_p
    z_grid = np.linspace(-9.0, 0.0, 5)
    stokes_profile = np.linspace(5.0e-4, 1.0e-3, 5)
    drift_profile = np.full(5, 1.0e-4)
    nu_profile = np.full(5, 1.0e-4)
    return ForcingState(
        U10=5.0,
        u_star_air=0.02,
        u_star_water=7.0e-4,
        U_surface=0.01,
        stokes_drift_surface=1.0e-3,
        stokes_drift_profile=stokes_profile,
        differential_drift_profile=drift_profile,
        z_grid=z_grid,
        depth=9.0,
        fetch=15000.0,
        H_s=0.3,
        T_p=2.5,
        f_p=1.0 / 2.5,
        omega_p=2.0 * math.pi / 2.5,
        k_p=k_p,
        lambda_p=lambda_p,
        X_tilde=1000.0,
        wave_steepness=0.5 * 0.3 * k_p,
        La_t=0.44,
        nu_T=4.0e-3,
        nu_T_profile=nu_profile,
        Ra=Ra,
        timestamp=datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc),
        drag_method="coare35",
        drift_method="webb_fox_kemper",
    )


class TestSubcritical:
    """Subcritical regime returns no LC."""

    def test_subcritical_returns_no_lc(self):
        forcing = _synthetic_forcing_state(Ra=50.0)
        env = EnvironmentalContext()
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=3600.0,
            environmental=env,
        )
        assert result["candidate"] == "multiscale"
        assert result["has_lc"] is False
        assert math.isnan(result["predicted_spacing"])
        assert result["regime"] == "subcritical"

    def test_subcritical_has_nan_coarsening(self):
        forcing = _synthetic_forcing_state(Ra=50.0)
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=3600.0,
            environmental=EnvironmentalContext(),
        )
        assert math.isnan(result["coarsening"]["initial_cell_width_m"])
        assert result["coarsening"]["n_events"] == 0


class TestSupercritical:
    """Supercritical regime with wave-tied initiation."""

    def test_returns_finite_spacing(self):
        forcing = _synthetic_forcing_state(Ra=20000.0)
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=3600.0,
            environmental=EnvironmentalContext(),
        )
        assert result["has_lc"] is True
        assert math.isfinite(result["predicted_spacing"])
        assert result["predicted_spacing"] > 0.0

    def test_initial_width_is_wave_tied(self):
        """Coarsening initial width must be d_LC = 0.34 * lambda_p, not L_inst."""
        forcing = _synthetic_forcing_state(Ra=20000.0, lambda_p=11.0)
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=600.0,
            environmental=EnvironmentalContext(),
        )
        expected_d_lc = wave_tied_lc_spacing(11.0)
        assert abs(result["coarsening"]["initial_cell_width_m"] - expected_d_lc) < 1e-10
        # Wave-tied initial must be much smaller than CL onset
        cl_onset = result["multiscale"]["cl_onset_spacing_m"]
        assert expected_d_lc < cl_onset * 0.5

    def test_both_scales_in_output(self):
        forcing = _synthetic_forcing_state(Ra=20000.0)
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=3600.0,
            environmental=EnvironmentalContext(),
        )
        ms = result["multiscale"]
        assert math.isfinite(ms["wave_tied_initial_m"])
        assert math.isfinite(ms["cl_onset_spacing_m"])
        assert ms["C_wave"] == 0.34
        assert ms["wave_tied_initial_m"] < ms["cl_onset_spacing_m"]

    def test_coarsening_from_wave_tied_produces_mergers(self):
        """Small wave-tied initial + reasonable lifetime should produce multiple mergers."""
        forcing = _synthetic_forcing_state(Ra=20000.0, lambda_p=11.0)
        # 10 minutes of pattern lifetime
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=600.0,
            environmental=EnvironmentalContext(),
        )
        n_events = result["coarsening"]["n_events"]
        # d_LC ~ 3.7m, A_H = max(nu_T, u*×h) = max(0.004, 0.0063)
        # tau_initial = 3.7^2 / (39.5 * 0.0063) ~ 55s
        # In 600s we should get multiple mergers
        assert n_events >= 1, f"Expected mergers, got n_events={n_events}"

    def test_zero_pattern_lifetime_gives_initial_scale(self):
        """With no time for coarsening, result should reflect initial wave-tied scale."""
        forcing = _synthetic_forcing_state(Ra=20000.0, lambda_p=11.0)
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=0.0,
            environmental=EnvironmentalContext(),
        )
        assert result["coarsening"]["n_events"] == 0
        # Predicted spacing should be close to wave-tied initial
        expected_d_lc = wave_tied_lc_spacing(11.0)
        lower = result["coarsening"]["visible_spacing_lower_m"]
        assert abs(lower - expected_d_lc) < 1e-6

    def test_lambda_p_sensitivity(self):
        """Different wave wavelengths produce different initial scales and predictions."""
        results = []
        for lp in [5.0, 11.0, 20.0]:
            f = _synthetic_forcing_state(Ra=20000.0, lambda_p=lp)
            r = analyse_candidate_multiscale(
                forcing=f,
                pattern_lifetime=600.0,
                environmental=EnvironmentalContext(),
            )
            results.append(r["multiscale"]["wave_tied_initial_m"])
        # Must be strictly increasing with lambda_p
        assert results[0] < results[1] < results[2]


class TestOutputCompatibility:
    """Output dict has all required keys for the evaluation module."""

    def test_required_keys_present(self):
        forcing = _synthetic_forcing_state(Ra=20000.0)
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=3600.0,
            environmental=EnvironmentalContext(),
        )
        required_top = [
            "candidate", "has_lc", "predicted_spacing", "predicted_spacing_m",
            "regime", "Ra", "La_t", "forcing_summary", "coarsening",
            "visibility", "lc_enhancement", "explanation",
        ]
        for key in required_top:
            assert key in result, f"Missing top-level key: {key}"

    def test_coarsening_keys_present(self):
        forcing = _synthetic_forcing_state(Ra=20000.0)
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=3600.0,
            environmental=EnvironmentalContext(),
        )
        coarsening_keys = [
            "initial_cell_width_m", "n_events", "coarsened_width_m",
            "cap_binding", "max_cell_aspect_ratio", "visible_spacing_m",
        ]
        for key in coarsening_keys:
            assert key in result["coarsening"], f"Missing coarsening key: {key}"

    def test_candidate_label(self):
        forcing = _synthetic_forcing_state(Ra=20000.0)
        result = analyse_candidate_multiscale(
            forcing=forcing,
            pattern_lifetime=3600.0,
            environmental=EnvironmentalContext(),
        )
        assert result["candidate"] == "multiscale"


class TestPipelineIntegration:
    """The multiscale candidate works through the pipeline entry point."""

    def test_analyse_case_multiscale(self):
        ts = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
        weather = pd.DataFrame({
            "timestamp": [ts - timedelta(hours=1), ts],
            "U10": [5.0, 5.0],
        })
        result = analyse_case(
            weather_data=weather,
            observation_time=ts,
            depth=9.0,
            fetch=15000.0,
            candidate="multiscale",
        )
        assert result["candidate_selected"] == "multiscale"
        assert result["candidate"] == "multiscale"
        assert "multiscale" in result

    def test_invalid_candidate_raises(self):
        ts = datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc)
        weather = pd.DataFrame({
            "timestamp": [ts],
            "U10": [5.0],
        })
        with pytest.raises(ValueError, match="Unknown candidate"):
            analyse_case(
                weather_data=weather,
                observation_time=ts,
                depth=9.0,
                fetch=15000.0,
                candidate="invalid",
            )
