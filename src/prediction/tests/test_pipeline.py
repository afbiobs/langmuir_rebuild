"""
Integration tests for the remaining WP-05 prediction pipeline slice.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import math

import numpy as np
import pandas as pd
import pytest

from src.forcing import ForcingState, compute_forcing
from src.hydro.robin_bc import RobinBC
from src.prediction.candidate_cl import analyse_candidate_cl
from src.prediction.candidate_scaling import analyse_candidate_scaling
from src.prediction.common import (
    build_environmental_context,
    cell_width_to_visible_spacing,
    instability_scale_to_cell_width,
)
from src.prediction.pipeline import analyse_case


def _synthetic_forcing_state(Ra: float) -> ForcingState:
    """Create a minimal forcing state for branch tests."""
    z_grid = np.linspace(-9.0, 0.0, 5)
    stokes_profile = np.linspace(5.0e-4, 1.0e-3, 5)
    drift_profile = np.full(5, 1.0e-4)
    nu_profile = np.full(5, 1.0e-4)
    return ForcingState(
        U10=1.0,
        u_star_air=0.02,
        u_star_water=7.0e-4,
        U_surface=0.01,
        stokes_drift_surface=1.0e-3,
        stokes_drift_profile=stokes_profile,
        differential_drift_profile=drift_profile,
        z_grid=z_grid,
        depth=9.0,
        fetch=15000.0,
        H_s=0.05,
        T_p=1.5,
        f_p=1.0 / 1.5,
        omega_p=2.0 * math.pi / 1.5,
        k_p=0.4,
        lambda_p=2.0 * math.pi / 0.4,
        X_tilde=1000.0,
        wave_steepness=0.5 * 0.05 * 0.4,
        La_t=0.8,
        nu_T=1.0e-4,
        nu_T_profile=nu_profile,
        Ra=Ra,
        timestamp=datetime(2026, 3, 22, 12, 0, tzinfo=timezone.utc),
        drag_method="coare35",
        drift_method="webb_fox_kemper",
    )


def _weather_history(final_u10: float = 5.0) -> pd.DataFrame:
    """Small weather history for public pipeline tests."""
    start = datetime(2026, 3, 22, 9, 0, tzinfo=timezone.utc)
    timestamps = [start + timedelta(minutes=30 * i) for i in range(7)]
    winds = [4.0, 4.2, 4.5, 4.8, 5.0, 5.1, final_u10]
    directions = [180.0] * len(winds)
    return pd.DataFrame(
        {
            "timestamp": timestamps,
            "U10": winds,
            "wind_direction_deg": directions,
        }
    )


def _direction_reset_then_steady_history() -> pd.DataFrame:
    """History with one early direction reset followed by steady supercritical winds."""
    start = datetime(2026, 3, 22, 9, 0, tzinfo=timezone.utc)
    timestamps = [start + timedelta(hours=i) for i in range(5)]
    return pd.DataFrame(
        {
            "timestamp": timestamps,
            "U10": [6.0, 6.1, 6.1, 6.0, 6.0],
            "wind_direction_deg": [180.0, 240.0, 242.0, 243.0, 245.0],
        }
    )


def _rapid_increase_then_steady_history() -> pd.DataFrame:
    """History with one early rapid-increase reset followed by steady winds."""
    start = datetime(2026, 3, 22, 9, 0, tzinfo=timezone.utc)
    timestamps = [start + timedelta(hours=i) for i in range(5)]
    return pd.DataFrame(
        {
            "timestamp": timestamps,
            "U10": [3.0, 3.2, 6.8, 6.9, 7.0],
            "wind_direction_deg": [180.0, 180.0, 180.0, 180.0, 180.0],
        }
    )


class TestCandidatePipelines:
    def test_cl_candidate_subcritical_returns_no_lc_and_unity_ratios(self):
        forcing = _synthetic_forcing_state(Ra=50.0)
        result = analyse_candidate_cl(
            forcing=forcing,
            pattern_lifetime=3600.0,
            environmental=build_environmental_context(),
            bc=RobinBC(gamma_s=0.06, gamma_b=0.28),
        )

        assert result["regime"] == "subcritical"
        assert not result["has_lc"]
        assert math.isnan(result["predicted_spacing_m"])
        assert not result["visibility"]["visible"]
        assert result["lc_enhancement"]["development_index"] == 0.0
        assert result["lc_enhancement"]["light"].ratio == 1.0
        assert result["lc_enhancement"]["nutrients"].ratio == 1.0
        assert result["lc_enhancement"]["temperature"].ratio == 1.0

    def test_cl_candidate_supercritical_produces_spacing(self):
        forcing = compute_forcing(U10=5.0, depth=9.0, fetch=15000.0)
        result = analyse_candidate_cl(
            forcing=forcing,
            pattern_lifetime=7200.0,
            environmental=build_environmental_context(),
        )

        assert result["candidate"] == "cl"
        assert result["has_lc"]
        assert result["predicted_spacing_m"] > 0.0
        assert (
            result["coarsening"]["visible_spacing_upper_m"]
            >= result["coarsening"]["visible_spacing_lower_m"]
        )
        assert result["predicted_spacing_m"] == pytest.approx(
            result["coarsening"]["visible_spacing_upper_m"]
        )
        assert result["intermediate"]["kappa"] > 1.0
        assert "l_cNL" in result["explanation"] or "predicted spacing" in result["explanation"]

    def test_cl_candidate_defaults_to_forcing_derived_robin_closure(self):
        forcing = compute_forcing(U10=5.0, depth=9.0, fetch=15_000.0)
        result = analyse_candidate_cl(
            forcing=forcing,
            pattern_lifetime=1800.0,
            environmental=build_environmental_context(),
        )

        closure = result["intermediate"]["robin_closure"]
        assert closure["source"] == "forcing_wave_steepness_bottom_reach"
        assert closure["gamma_s_raw"] == pytest.approx(forcing.wave_steepness)
        assert closure["gamma_b_raw"] > 0.0
        assert closure["gamma_total_raw"] == pytest.approx(
            closure["gamma_s_raw"] + closure["gamma_b_raw"]
        )

    def test_cl_candidate_uses_critical_onset_width_before_coarsening(self):
        forcing = compute_forcing(U10=5.0, depth=9.0, fetch=15000.0)
        onset = analyse_candidate_cl(
            forcing=forcing,
            pattern_lifetime=1800.0,
            environmental=build_environmental_context(),
            onset_only=True,
        )
        result = analyse_candidate_cl(
            forcing=forcing,
            pattern_lifetime=1800.0,
            environmental=build_environmental_context(),
        )

        critical_l = result["intermediate"]["critical_result"].l_c
        selected_l = result["intermediate"]["selected_wavenumber"]
        onset_width = instability_scale_to_cell_width(critical_l, forcing.depth)

        assert result["intermediate"]["selected_wavenumber_method"] == "critical_onset_width"
        assert selected_l == pytest.approx(critical_l)
        assert result["coarsening"]["initial_cell_width_m"] == pytest.approx(onset_width)
        assert result["coarsening"]["initial_cell_width_m"] == pytest.approx(
            onset["predicted_spacing_m"]
        )

    def test_cl_candidate_onset_only_bypasses_selector_and_coarsening(self):
        forcing = compute_forcing(U10=5.0, depth=9.0, fetch=15_000.0)
        short = analyse_candidate_cl(
            forcing=forcing,
            pattern_lifetime=600.0,
            environmental=build_environmental_context(),
            onset_only=True,
        )
        long = analyse_candidate_cl(
            forcing=forcing,
            pattern_lifetime=7200.0,
            environmental=build_environmental_context(),
            onset_only=True,
        )

        assert short["onset_only"] is True
        assert short["intermediate"]["selected_wavenumber_method"] == "critical_onset_only"
        assert short["coarsening"]["n_events"] == 0
        assert short["predicted_spacing_m"] == pytest.approx(
            short["coarsening"]["initial_cell_width_m"]
        )
        assert long["predicted_spacing_m"] == pytest.approx(short["predicted_spacing_m"])
        assert "Onset-only validation mode is active" in short["explanation"]

    def test_scaling_candidate_spacing_grows_with_coarsening_time(self):
        forcing = compute_forcing(U10=5.0, depth=9.0, fetch=15000.0)
        short = analyse_candidate_scaling(
            forcing=forcing,
            pattern_lifetime=600.0,
            environmental=build_environmental_context(),
        )
        long = analyse_candidate_scaling(
            forcing=forcing,
            pattern_lifetime=7200.0,
            environmental=build_environmental_context(),
        )

        assert short["predicted_spacing_m"] > 0.0
        assert long["predicted_spacing_m"] >= short["predicted_spacing_m"]

    def test_visible_spacing_proxy_can_exceed_capped_cell_width(self):
        forcing = compute_forcing(U10=5.5, depth=2.0, fetch=30000.0)
        result = analyse_candidate_scaling(
            forcing=forcing,
            pattern_lifetime=4.0 * 3600.0,
            environmental=build_environmental_context(),
        )

        assert result["coarsening"]["cap_binding"] is True
        assert result["coarsening"]["raw_coarsened_width_m"] > result["coarsening"]["coarsened_width_m"]
        assert result["predicted_spacing_m"] == pytest.approx(
            result["coarsening"]["visible_spacing_upper_m"]
        )
        assert result["coarsening"]["visible_cap_binding"] is False
        assert result["coarsening"]["visible_n_events"] <= 3

    def test_visible_spacing_proxy_limits_visible_merger_levels(self):
        forcing = compute_forcing(U10=11.0, depth=9.0, fetch=15000.0)
        result = analyse_candidate_scaling(
            forcing=forcing,
            pattern_lifetime=12.0 * 3600.0,
            environmental=build_environmental_context(),
        )

        assert result["coarsening"]["n_events"] > result["coarsening"]["visible_n_events"]
        assert result["coarsening"]["visible_n_events"] == 3
        assert result["predicted_spacing_m"] == pytest.approx(
            result["coarsening"]["visible_spacing_upper_m"]
        )
        assert result["coarsening"]["visible_spacing_upper_m"] < result["coarsening"]["raw_coarsened_width_m"]

    def test_cl_candidate_detects_a_ten_percent_wind_change(self):
        base = analyse_candidate_cl(
            forcing=compute_forcing(U10=5.0, depth=9.0, fetch=15000.0),
            pattern_lifetime=1300.0,
            environmental=build_environmental_context(),
        )
        perturbed = analyse_candidate_cl(
            forcing=compute_forcing(U10=5.5, depth=9.0, fetch=15000.0),
            pattern_lifetime=1300.0,
            environmental=build_environmental_context(),
        )

        assert base["predicted_spacing_m"] != perturbed["predicted_spacing_m"]

    def test_cell_width_to_visible_spacing_uses_loose_observation_cap(self):
        result = cell_width_to_visible_spacing(
            cell_width=24.0,
            raw_hierarchy_width=650.0,
            depth=2.0,
        )

        assert result["visible_spacing_m"] == pytest.approx(600.0)
        assert result["visible_cap_binding"] is True


class TestPublicPipeline:
    def test_analyse_case_runs_for_cl_candidate(self):
        history = _weather_history(final_u10=5.0)
        result = analyse_case(
            weather_data=history,
            observation_time=history.iloc[-1]["timestamp"].to_pydatetime(),
            depth=9.0,
            fetch=15000.0,
            candidate="cl",
        )

        assert result["candidate"] == "cl"
        assert result["candidate_selected"] == "cl"
        assert result["history_points"] >= 2
        assert result["pattern_lifetime_s"] >= 0.0
        assert "visibility" in result
        assert "lc_enhancement" in result

    def test_analyse_case_passes_onset_only_to_cl_candidate(self):
        history = _weather_history(final_u10=5.0)
        result = analyse_case(
            weather_data=history,
            observation_time=history.iloc[-1]["timestamp"].to_pydatetime(),
            depth=9.0,
            fetch=15_000.0,
            candidate="cl",
            onset_only=True,
        )

        assert result["onset_only"] is True
        assert result["coarsening"]["n_events"] == 0
        assert result["intermediate"]["selected_wavenumber_method"] == "critical_onset_only"

    def test_analyse_case_runs_for_scaling_candidate(self):
        history = _weather_history(final_u10=5.8)
        result = analyse_case(
            weather_data=history,
            observation_time=history.iloc[-1]["timestamp"].to_pydatetime(),
            depth=9.0,
            fetch=15000.0,
            candidate="scaling",
            environmental={"surface_irradiance": 350.0, "K_d": 0.4},
        )

        assert result["candidate"] == "scaling"
        assert result["candidate_selected"] == "scaling"
        assert result["predicted_spacing_m"] > 0.0
        assert result["lc_enhancement"]["development_index"] >= 0.0

    def test_analyse_case_uses_time_since_last_direction_reset(self):
        history = _direction_reset_then_steady_history()
        result = analyse_case(
            weather_data=history,
            observation_time=history.iloc[-1]["timestamp"].to_pydatetime(),
            depth=9.0,
            fetch=15_000.0,
            candidate="cl",
            onset_only=True,
            lookback_hours=6.0,
        )

        assert result["pattern_lifetime_s"] == pytest.approx(3.0 * 3600.0)
        assert result["disruption"]["direction_change"] is True
        assert result["disruption"]["reset_time"] == history.iloc[1]["timestamp"].to_pydatetime()
        assert result["disruption"]["event_count"] == 1
        assert result["disruption"]["reset_causes"] == ["direction_change"]

    def test_analyse_case_uses_time_since_last_rapid_increase_reset(self):
        history = _rapid_increase_then_steady_history()
        result = analyse_case(
            weather_data=history,
            observation_time=history.iloc[-1]["timestamp"].to_pydatetime(),
            depth=9.0,
            fetch=15_000.0,
            candidate="cl",
            onset_only=True,
            lookback_hours=6.0,
        )

        assert result["pattern_lifetime_s"] == pytest.approx(2.0 * 3600.0)
        assert result["disruption"]["rapid_speed_increase"] is True
        assert result["disruption"]["reset_time"] == history.iloc[2]["timestamp"].to_pydatetime()
        assert result["disruption"]["event_count"] == 1
        assert result["disruption"]["reset_causes"] == ["rapid_speed_increase"]
