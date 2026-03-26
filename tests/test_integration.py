"""
Integration tests for the WP-05 observation-set comparison harness.
"""

from __future__ import annotations

from datetime import timedelta
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.evaluation.comparison import (
    _era5_cache_key,
    load_observations,
    run_wp05_comparison,
    verify_production_candidate_reproduction,
)


def _mini_observations() -> pd.DataFrame:
    """Minimal four-site observation table for comparison-harness tests."""
    return pd.DataFrame(
        [
            {
                "observation_id": "obs_taihu",
                "image_date": "2021-08-17",
                "authoritative_lat": 31.29812,
                "authoritative_lng": 120.01688,
                "manual_spacing_m": float("nan"),
                "wiggle_spacing_m": 320.0,
            },
            {
                "observation_id": "obs_neagh",
                "image_date": "2024-02-12",
                "authoritative_lat": 54.62757,
                "authoritative_lng": -6.47327,
                "manual_spacing_m": 55.0,
                "wiggle_spacing_m": float("nan"),
            },
            {
                "observation_id": "obs_prairie",
                "image_date": "2023-07-04",
                "authoritative_lat": 48.01353,
                "authoritative_lng": -98.95995,
                "manual_spacing_m": 78.0,
                "wiggle_spacing_m": 82.0,
            },
            {
                "observation_id": "obs_erie",
                "image_date": "2023-09-01",
                "authoritative_lat": 41.85671,
                "authoritative_lng": -83.12900,
                "manual_spacing_m": 410.0,
                "wiggle_spacing_m": float("nan"),
            },
        ]
    )


def _write_spinup_cache(cache_dir: Path, observation_row: pd.Series) -> None:
    """Write a minimal ERA5-style cache file for one observation."""
    image_date = pd.Timestamp(observation_row["image_date"]).date()
    start = pd.Timestamp(image_date - timedelta(days=10), tz="UTC")
    end = pd.Timestamp(image_date, tz="UTC") + pd.Timedelta(hours=23)
    times = pd.date_range(start=start, end=end, freq="1h", tz="UTC")
    hours = np.arange(len(times), dtype=float)
    base_speed = 4.0 + 0.01 * hours
    payload = {
        "latitude": float(observation_row["authoritative_lat"]),
        "longitude": float(observation_row["authoritative_lng"]),
        "generationtime_ms": 0.0,
        "utc_offset_seconds": 0,
        "timezone": "GMT",
        "timezone_abbreviation": "GMT",
        "elevation": 0.0,
        "hourly_units": {
            "time": "iso8601",
            "wind_speed_10m": "m/s",
            "wind_direction_10m": "deg",
            "wind_gusts_10m": "m/s",
            "surface_pressure": "hPa",
            "temperature_2m": "degC",
            "shortwave_radiation": "W/m2",
            "cloud_cover": "%",
            "precipitation": "mm",
        },
        "hourly": {
            "time": [stamp.strftime("%Y-%m-%dT%H:%M") for stamp in times],
            "wind_speed_10m": base_speed.tolist(),
            "wind_direction_10m": (180.0 + 0.1 * hours).tolist(),
            "wind_gusts_10m": (base_speed + 1.0).tolist(),
            "surface_pressure": (1015.0 + 0.0 * hours).tolist(),
            "temperature_2m": (18.0 + 0.0 * hours).tolist(),
            "shortwave_radiation": (150.0 + 0.0 * hours).tolist(),
            "cloud_cover": (40.0 + 0.0 * hours).tolist(),
            "precipitation": (0.0 * hours).tolist(),
        },
    }
    key = _era5_cache_key(
        float(observation_row["authoritative_lat"]),
        float(observation_row["authoritative_lng"]),
        (image_date - timedelta(days=10)).isoformat(),
        image_date.isoformat(),
    )
    cache_file = cache_dir / f"{key}.json"
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(json.dumps(payload), encoding="utf-8")


def test_load_observations_prefers_manual_spacing_and_tracks_method():
    loaded = load_observations(_mini_observations())
    loaded = loaded.set_index("observation_id")

    assert loaded.loc["obs_taihu", "target_spacing_m"] == 320.0
    assert loaded.loc["obs_taihu", "measurement_method"] == "wiggle"
    assert loaded.loc["obs_prairie", "target_spacing_m"] == 78.0
    assert loaded.loc["obs_prairie", "measurement_method"] == "manual"
    assert bool(loaded.loc["obs_prairie", "has_manual"])
    assert bool(loaded.loc["obs_prairie", "has_wiggle"])


def test_run_wp05_comparison_returns_all_models_and_writes_outputs(tmp_path: Path):
    result = run_wp05_comparison(_mini_observations(), output_dir=tmp_path)

    assert result["forcing_mode"] == "provisional_site_proxy"
    assert len(result["case_table"]) == 4
    assert len(result["prediction_table"]) == 20
    assert set(result["prediction_table"]["model"]) == {
        "cl",
        "scaling",
        "baseline_constant",
        "baseline_linear_wind",
        "baseline_depth_scaled",
    }

    full_rows = result["metrics_table"].loc[
        result["metrics_table"]["subset_id"] == "BM-A_full"
    ]
    assert len(full_rows) == 5
    assert result["decision_gate"]["provisional"] is True
    assert "question_1_cl_vs_scaling" in result["decision_gate"]

    assert (tmp_path / "case_inputs.csv").exists()
    assert (tmp_path / "observation_predictions.csv").exists()
    assert (tmp_path / "metrics_table.csv").exists()
    assert (tmp_path / "decision_gate_summary.json").exists()


def test_run_wp05_comparison_accepts_matched_cache_dir(tmp_path: Path):
    observations = _mini_observations()
    cache_dir = tmp_path / "era5_cache"
    for _, row in observations.iterrows():
        _write_spinup_cache(cache_dir, row)

    result = run_wp05_comparison(
        observations,
        output_dir=tmp_path,
        era5_cache_dir=cache_dir,
        lookback_hours=6.0,
    )

    assert result["forcing_mode"] == "era5_cache_matched"
    assert result["era5_cache_dir"] == str(cache_dir)
    assert result["decision_gate"]["provisional"] is False
    assert result["scenario"]["max_visible_mergers"] == 3
    assert set(result["case_table"]["forcing_mode"]) == {"era5_cache_matched"}
    candidate_rows = result["prediction_table"].loc[
        result["prediction_table"]["model"].isin(["cl", "scaling"])
    ]
    assert set(candidate_rows["forcing_mode"]) == {"era5_cache_matched"}
    assert candidate_rows["matched_cache_file"].astype(str).str.len().gt(0).all()
    assert candidate_rows["pattern_lifetime_s"].notna().all()
    assert candidate_rows["visible_spacing_upper_m"].notna().all()
    assert (candidate_rows["n_visible_events"] <= candidate_rows["n_coarsening_events"]).all()
    assert (candidate_rows["coarsening_diffusivity_m2_s"] >= candidate_rows["forcing_nu_T_m2_s"]).all()
    assert set(candidate_rows["coarsening_diffusivity_method"]) == {"mixing_length_lateral"}


def test_run_wp05_comparison_supports_onset_only_cl_validation(tmp_path: Path):
    result = run_wp05_comparison(
        _mini_observations(),
        output_dir=tmp_path,
        onset_only=True,
    )

    assert result["scenario"]["onset_only"] is True
    cl_rows = result["prediction_table"].loc[result["prediction_table"]["model"] == "cl"]
    assert cl_rows["onset_only"].all()
    assert (cl_rows["n_coarsening_events"] == 0).all()
    assert set(cl_rows["selected_wavenumber_method"]) == {"critical_onset_only"}
    assert (cl_rows["robin_gamma_total"] > 0.0).all()


def test_run_wp05_comparison_writes_wp06_output_package(tmp_path: Path):
    result = run_wp05_comparison(_mini_observations(), output_dir=tmp_path)

    expected_files = [
        "case_inputs.csv",
        "observation_predictions.csv",
        "metrics_table.csv",
        "decision_gate_summary.json",
        "coarsening_closure_audit.json",
        "predicted_vs_observed_cl.png",
        "predicted_vs_observed_scaling.png",
        "spacing_vs_wind_cl.png",
        "spacing_vs_wind_scaling.png",
        "dynamic_range_comparison.png",
        "tail_coverage_comparison.png",
        "attractor_diagnostic_cl.png",
        "attractor_diagnostic_scaling.png",
        "production_candidate_assessment.json",
    ]
    for relative_path in expected_files:
        assert (tmp_path / relative_path).exists(), relative_path

    representative_case_ids = result["representative_case_ids"]
    assert len(representative_case_ids) == 3
    for observation_id in representative_case_ids:
        assert (tmp_path / f"enhancement_index_timeseries_{observation_id}.png").exists()

    case_diagnostics_dir = tmp_path / "case_diagnostics"
    case_files = sorted(case_diagnostics_dir.glob("case_*.json"))
    assert len(case_files) == len(_mini_observations())

    case_payload = json.loads(
        (case_diagnostics_dir / "case_obs_neagh.json").read_text(encoding="utf-8")
    )
    assert case_payload["comparison_context"]["forcing_mode"] == "provisional_site_proxy"
    assert case_payload["assumptions"]["shallow_lake_stratification"]["model_parameter"] == "S = 0"
    assert case_payload["models"]["cl"]["prediction_row"]["model"] == "cl"
    assert case_payload["models"]["scaling"]["prediction_row"]["model"] == "scaling"


def test_run_wp05_comparison_marks_provisional_assessment_as_non_promotable(tmp_path: Path):
    result = run_wp05_comparison(_mini_observations(), output_dir=tmp_path)

    assessment = result["production_candidate_assessment"]
    assert assessment["forcing_mode"] == "provisional_site_proxy"
    assert assessment["promotable"] is False
    assert assessment["selected_candidate"] is None
    assert any("S = 0" in note for note in assessment["notes"])


def test_verify_production_candidate_reproduction_matches_matched_comparison(tmp_path: Path):
    observations = _mini_observations()
    cache_dir = tmp_path / "era5_cache"
    for _, row in observations.iterrows():
        _write_spinup_cache(cache_dir, row)

    result = run_wp05_comparison(
        observations,
        output_dir=tmp_path,
        era5_cache_dir=cache_dir,
        lookback_hours=6.0,
    )

    verification = verify_production_candidate_reproduction(result, candidate="cl")

    assert verification["status"] == "passed"
    assert verification["candidate"] == "cl"
    assert verification["n_cases"] == len(observations)
    assert verification["max_abs_spacing_diff_m"] == pytest.approx(0.0)
    assert verification["max_abs_ra_diff"] == pytest.approx(0.0)


def test_real_observation_file_loads_expected_case_count():
    loaded = load_observations(Path("data/raw/observations.csv"))
    assert len(loaded) == 67
    assert loaded["target_spacing_m"].min() > 0.0
