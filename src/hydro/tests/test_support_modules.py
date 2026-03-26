"""
Support-module tests for WP-04 hydro deliverables outside the H&P verification file.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import numpy as np
import pytest

from src.evaluation.metrics import saturation_audit
from src.forcing import compute_forcing
from src.hydro.coarsening import (
    coarsened_width,
    coarsening_diffusivity,
    coarsening_schedule,
    coarsening_timescale,
    count_coarsening_events,
    disruption_check,
)
from src.hydro.rayleigh import classify_regime, compute_rayleigh, unstable_band
from src.hydro.robin_bc import derive_robin_bc_from_forcing
from src.hydro.scaling_laws import empirical_spacing_wind, la_dependent_geometry


class TestRayleighSupport:
    def test_compute_rayleigh_dimensional_scaling(self):
        """Ra scales linearly in U and D, quadratically in h and ν_T^{-2}."""
        base = compute_rayleigh(U_surface=0.05, D_max=0.02, depth=9.0, nu_T=1.0e-3)
        assert compute_rayleigh(0.10, 0.02, 9.0, 1.0e-3) == pytest.approx(2.0 * base)
        assert compute_rayleigh(0.05, 0.02, 18.0, 1.0e-3) == pytest.approx(4.0 * base)
        assert compute_rayleigh(0.05, 0.02, 9.0, 2.0e-3) == pytest.approx(base / 4.0)

    def test_classify_regime_thresholds(self):
        """Regime labels follow the WP-04 thresholds exactly."""
        assert classify_regime(100.0, R0=120.0, RcNL=122.194) == "subcritical"
        assert classify_regime(120.0, R0=120.0, RcNL=122.194) == "near_onset"
        assert classify_regime(190.0, R0=120.0, RcNL=122.194) == "moderate"
        assert classify_regime(700.0, R0=120.0, RcNL=122.194) == "supercritical"

    def test_unstable_band_returns_nan_when_subcritical(self):
        """No unstable band is reported when Ra stays below the neutral curve."""
        band = unstable_band(
            Ra=90.0,
            neutral_curve=lambda l: 100.0 + 0.0 * l,
            l_array=[0.1, 0.2, 0.3],
        )
        assert band[0] != band[0]
        assert band[1] != band[1]

    def test_unstable_band_returns_extent_when_supercritical(self):
        """Unstable-band bounds span the portion of l where Ra exceeds R̄(l)."""
        l_array = [0.1, 0.2, 0.3, 0.4, 0.5]
        band = unstable_band(
            Ra=120.0,
            neutral_curve=lambda l: 100.0 + 400.0 * (l - 0.3) ** 2,
            l_array=l_array,
        )
        assert band == pytest.approx((0.1, 0.5))


class TestScalingLaws:
    def test_la_dependent_geometry_is_positive_and_monotonic(self):
        """The La-based geometry follows the documented power laws."""
        low = la_dependent_geometry(La_t=0.3, depth=9.0, u_star=0.004)
        high = la_dependent_geometry(La_t=0.9, depth=9.0, u_star=0.004)

        assert high["downwelling_thickness"] > low["downwelling_thickness"]
        assert high["downwelling_velocity_max"] < low["downwelling_velocity_max"]
        assert high["pitch"] > low["pitch"]
        assert high["estimated_cell_width"] > low["estimated_cell_width"]
        assert all(value > 0.0 for value in low.values())

    def test_la_dependent_geometry_no_saturation_over_relevant_ranges(self):
        """The scalar scaling-law outputs vary meaningfully over plausible inputs."""
        width_audit = saturation_audit(
            lambda La_t, depth: {
                "estimated_cell_width": la_dependent_geometry(
                    La_t=La_t, depth=depth, u_star=0.004
                )["estimated_cell_width"]
            },
            input_ranges={"La_t": (0.2, 1.2), "depth": (5.0, 15.0)},
        )
        vel_audit = saturation_audit(
            lambda La_t, u_star: {
                "downwelling_velocity_max": la_dependent_geometry(
                    La_t=La_t, depth=9.0, u_star=u_star
                )["downwelling_velocity_max"]
            },
            input_ranges={"La_t": (0.2, 1.2), "u_star": (0.001, 0.02)},
        )

        assert width_audit["passed"], width_audit["findings"]
        assert vel_audit["passed"], vel_audit["findings"]

    def test_empirical_spacing_ranges_are_ordered(self):
        """The empirical references return ordered, depth-scaled spacing bounds."""
        refs = empirical_spacing_wind(U10=6.0, depth=9.0)
        assert refs["faller_caponi_min_spacing"] == pytest.approx(18.0)
        assert refs["faller_caponi_max_spacing"] == pytest.approx(45.0)
        assert refs["smith_min_spacing"] == pytest.approx(18.0)
        assert refs["smith_max_spacing"] == pytest.approx(27.0)
        assert refs["faller_caponi_min_spacing"] < refs["faller_caponi_max_spacing"]
        assert refs["smith_min_spacing"] < refs["smith_max_spacing"]


class TestRobinClosure:
    def test_forcing_derived_robin_bc_tracks_wave_steepness_and_bed_reach(self):
        """The Robin closure should expose the raw wave-steepness split without clipping."""
        forcing = compute_forcing(U10=5.0, depth=9.0, fetch=15_000.0)
        bc, diagnostics = derive_robin_bc_from_forcing(forcing)

        assert bc.gamma_s == pytest.approx(forcing.wave_steepness)
        assert diagnostics["gamma_s_raw"] == pytest.approx(forcing.wave_steepness)
        assert diagnostics["gamma_b_raw"] == pytest.approx(
            forcing.wave_steepness / np.sinh(forcing.k_p * forcing.depth)
        )
        assert diagnostics["gamma_total_raw"] == pytest.approx(bc.gamma)
        assert diagnostics["bottom_coupling_factor"] > 0.0
        assert diagnostics["source"] == "forcing_wave_steepness_bottom_reach"

    def test_forcing_derived_robin_bc_has_no_hidden_saturation(self):
        """The Robin closure should retain dynamic range over plausible wind/depth changes."""
        audit = saturation_audit(
            lambda U10, depth: {
                key: value
                for key, value in derive_robin_bc_from_forcing(
                    compute_forcing(U10=U10, depth=depth, fetch=15_000.0)
                )[1].items()
                if key in {"gamma_s_raw", "gamma_b_raw", "gamma_total_raw"}
            },
            input_ranges={"U10": (1.0, 15.0), "depth": (2.0, 15.0)},
        )

        assert audit["passed"], audit["findings"]


class TestCoarsening:
    def test_coarsening_diffusivity_uses_lateral_mixing_length_scale(self):
        """The merger clock should use a lateral diffusivity scale larger than ν_T,z."""
        u_star_water = 4.0e-3
        depth = 9.0
        nu_T_vertical = 0.41 * u_star_water * depth / 6.0

        closure = coarsening_diffusivity(
            nu_T_vertical=nu_T_vertical,
            u_star_water=u_star_water,
            depth=depth,
        )

        assert closure["method"] == "mixing_length_lateral"
        assert closure["vertical_nu_T_m2_s"] == pytest.approx(nu_T_vertical)
        assert closure["lateral_mixing_length_diffusivity_m2_s"] == pytest.approx(
            u_star_water * depth
        )
        assert closure["diffusivity_m2_s"] == pytest.approx(u_star_water * depth)
        assert closure["anisotropy_ratio"] == pytest.approx(6.0 / 0.41)

    def test_coarsening_timescale_and_event_count(self):
        """Merger time follows the slow nonlinear λ²/((2π)²ν_T) scaling."""
        tau = coarsening_timescale(cell_width=20.0, nu_T=1.0e-2)
        expected = 20.0 ** 2 / (((2.0 * 3.141592653589793) ** 2) * 1.0e-2)
        assert tau == pytest.approx(expected)

        schedule = coarsening_schedule(
            time_available=4.9 * tau,
            initial_width=20.0,
            nu_T=1.0e-2,
        )
        assert schedule["n_events"] == 1
        assert schedule["tau_next_s"] == pytest.approx(4.0 * tau)
        assert count_coarsening_events(
            time_available=5.1 * tau,
            initial_width=20.0,
            nu_T=1.0e-2,
        ) == 2

    def test_coarsened_width_doubles_and_warns_at_cap(self):
        """Each merger doubles the width and the aspect-ratio cap warns when binding."""
        assert coarsened_width(initial_width=20.0, n_events=2, depth=20.0) == pytest.approx(80.0)

        with pytest.warns(UserWarning):
            capped = coarsened_width(
                initial_width=20.0,
                n_events=2,
                depth=9.0,
                max_aspect_ratio=6.0,
            )
        assert capped == pytest.approx(54.0)

    def test_disruption_check_flags_direction_change(self):
        """Direction changes >45° reset the structure."""
        t0 = datetime(2026, 3, 22, 9, 0, tzinfo=timezone.utc)
        history = [
            {"timestamp": t0, "U10": 5.0, "wind_direction_deg": 10.0},
            {"timestamp": t0 + timedelta(hours=1), "U10": 5.2, "wind_direction_deg": 20.0},
            {"timestamp": t0 + timedelta(hours=2), "U10": 5.1, "wind_direction_deg": 80.0},
        ]
        result = disruption_check(history, lookback_hours=3.0)
        assert result["direction_change"] is True
        assert result["disrupted"] is True
        assert result["reset_time"] == history[-1]["timestamp"]

    def test_disruption_check_reports_latest_reset_time_after_earlier_direction_change(self):
        """An earlier direction reset should not be re-timed to the latest stable sample."""
        t0 = datetime(2026, 3, 22, 9, 0, tzinfo=timezone.utc)
        history = [
            {"timestamp": t0, "U10": 5.0, "wind_direction_deg": 180.0},
            {"timestamp": t0 + timedelta(hours=1), "U10": 5.1, "wind_direction_deg": 240.0},
            {"timestamp": t0 + timedelta(hours=2), "U10": 5.2, "wind_direction_deg": 242.0},
            {"timestamp": t0 + timedelta(hours=3), "U10": 5.3, "wind_direction_deg": 245.0},
        ]

        result = disruption_check(history, lookback_hours=6.0)

        assert result["direction_change"] is True
        assert result["disrupted"] is True
        assert result["event_count"] == 1
        assert result["reset_time"] == history[1]["timestamp"]
        assert result["reset_causes"] == ["direction_change"]
        assert result["current_direction_deg"] == pytest.approx(245.0)

    def test_disruption_check_flags_low_wind_and_rapid_increase(self):
        """Low-wind shutdowns and sharp wind increases are both detected."""
        t0 = datetime(2026, 3, 22, 9, 0, tzinfo=timezone.utc)
        low_wind = [
            {"timestamp": t0, "U10": 4.5},
            {"timestamp": t0 + timedelta(hours=1), "U10": 1.8},
        ]
        rapid = [
            {"timestamp": t0, "U10": 3.0},
            {"timestamp": t0 + timedelta(hours=1), "U10": 3.2},
            {"timestamp": t0 + timedelta(hours=2), "U10": 6.8},
        ]

        low_result = disruption_check(low_wind, lookback_hours=3.0)
        rapid_result = disruption_check(rapid, lookback_hours=3.0)

        assert low_result["low_wind_shutdown"] is True
        assert low_result["disrupted"] is True
        assert rapid_result["rapid_speed_increase"] is True
        assert rapid_result["disrupted"] is True

    def test_disruption_check_reports_latest_reset_time_after_earlier_rapid_increase(self):
        """A rapid-increase reset should preserve its own event time after later steady winds."""
        t0 = datetime(2026, 3, 22, 9, 0, tzinfo=timezone.utc)
        history = [
            {"timestamp": t0, "U10": 3.0, "wind_direction_deg": 180.0},
            {"timestamp": t0 + timedelta(hours=1), "U10": 3.2, "wind_direction_deg": 180.0},
            {"timestamp": t0 + timedelta(hours=2), "U10": 6.8, "wind_direction_deg": 180.0},
            {"timestamp": t0 + timedelta(hours=3), "U10": 6.9, "wind_direction_deg": 180.0},
            {"timestamp": t0 + timedelta(hours=4), "U10": 7.0, "wind_direction_deg": 180.0},
        ]

        result = disruption_check(history, lookback_hours=6.0)

        assert result["rapid_speed_increase"] is True
        assert result["disrupted"] is True
        assert result["event_count"] == 1
        assert result["reset_time"] == history[2]["timestamp"]
        assert result["reset_causes"] == ["rapid_speed_increase"]
