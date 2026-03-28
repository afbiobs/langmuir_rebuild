"""Tests for the multiscale coherent structure module."""

import math
import pytest

from src.hydro.multiscale_structures import (
    MultiscaleResult,
    compute_multiscale_result,
    wave_tied_lc_spacing,
)


class TestWaveTiedLCSpacing:
    """Tests for wave_tied_lc_spacing()."""

    def test_scales_linearly_with_lambda_p(self):
        """Doubling lambda_p doubles d_LC."""
        d1 = wave_tied_lc_spacing(10.0)
        d2 = wave_tied_lc_spacing(20.0)
        assert abs(d2 / d1 - 2.0) < 1.0e-12

    def test_default_coefficient_is_034(self):
        """Default C_wave = 0.34 from Tsai & Lu (2023)."""
        d = wave_tied_lc_spacing(100.0)
        assert abs(d - 34.0) < 1.0e-12

    def test_tsai_lu_nondimensional_wavenumber(self):
        """C_wave = 0.34 corresponds to l_s = 2*pi/0.34 ~ 18.48 -> l_s/2pi ~ 2.94."""
        # Tsai & Lu report l_s ~ 2.9.  Check: lambda/d_s = 1/C_wave = 2.94
        assert abs(1.0 / 0.34 - 2.941) < 0.01

    def test_custom_coefficient(self):
        d = wave_tied_lc_spacing(10.0, C_wave=0.5)
        assert abs(d - 5.0) < 1.0e-12

    def test_negative_lambda_p_raises(self):
        with pytest.raises(ValueError, match="lambda_p"):
            wave_tied_lc_spacing(-1.0)

    def test_zero_lambda_p_raises(self):
        with pytest.raises(ValueError, match="lambda_p"):
            wave_tied_lc_spacing(0.0)

    def test_nan_lambda_p_raises(self):
        with pytest.raises(ValueError, match="lambda_p"):
            wave_tied_lc_spacing(float("nan"))

    def test_negative_C_wave_raises(self):
        with pytest.raises(ValueError, match="C_wave"):
            wave_tied_lc_spacing(10.0, C_wave=-0.1)

    def test_typical_neagh_value(self):
        """Sanity check: lambda_p ~ 11 m at Neagh -> d_LC ~ 3.7 m."""
        d = wave_tied_lc_spacing(10.8)
        assert 3.0 < d < 5.0


class TestComputeMultiscaleResult:
    """Tests for compute_multiscale_result()."""

    def test_returns_frozen_dataclass(self):
        r = compute_multiscale_result(lambda_p=10.0, depth=9.0, cl_onset_spacing=196.0)
        assert isinstance(r, MultiscaleResult)
        with pytest.raises(AttributeError):
            r.wave_lc_spacing_m = 999.0  # type: ignore[misc]

    def test_wave_tied_spacing_matches_standalone(self):
        r = compute_multiscale_result(lambda_p=10.0, depth=9.0, cl_onset_spacing=196.0)
        expected = wave_tied_lc_spacing(10.0)
        assert abs(r.wave_lc_spacing_m - expected) < 1.0e-12

    def test_cl_onset_preserved(self):
        r = compute_multiscale_result(lambda_p=10.0, depth=9.0, cl_onset_spacing=196.0)
        assert abs(r.cl_onset_spacing_m - 196.0) < 1.0e-12

    def test_aspect_ratios_computed(self):
        r = compute_multiscale_result(lambda_p=10.0, depth=5.0, cl_onset_spacing=100.0)
        assert abs(r.aspect_ratio_wave - r.wave_lc_spacing_m / 5.0) < 1.0e-12
        assert abs(r.aspect_ratio_onset - 20.0) < 1.0e-12

    def test_wave_tied_smaller_than_onset(self):
        """At typical conditions, wave-tied scale << CL onset scale."""
        r = compute_multiscale_result(lambda_p=11.0, depth=9.0, cl_onset_spacing=196.0)
        assert r.wave_lc_spacing_m < r.cl_onset_spacing_m

    def test_negative_depth_raises(self):
        with pytest.raises(ValueError, match="depth"):
            compute_multiscale_result(lambda_p=10.0, depth=-1.0, cl_onset_spacing=100.0)

    def test_negative_cl_onset_raises(self):
        with pytest.raises(ValueError, match="cl_onset_spacing"):
            compute_multiscale_result(lambda_p=10.0, depth=9.0, cl_onset_spacing=-1.0)

    def test_no_saturation_across_wind_range(self):
        """Wave-tied scale varies meaningfully across realistic lambda_p range."""
        # lambda_p ~ 4-27 m for U10 = 3-10 m/s, fetch = 5-30 km
        spacings = [wave_tied_lc_spacing(lp) for lp in [4.0, 8.0, 12.0, 20.0, 27.0]]
        # Check that spread is > 5% of mean
        spread = max(spacings) - min(spacings)
        mean = sum(spacings) / len(spacings)
        assert spread / mean > 0.5, (
            f"Insufficient dynamic range: spread={spread:.2f}, mean={mean:.2f}"
        )
