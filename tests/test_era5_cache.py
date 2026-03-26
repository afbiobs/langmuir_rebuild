"""
Tests for the reusable ERA5/Open-Meteo cache filler.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from src.data.era5_cache import (
    ERA5_ARCHIVE_HOURLY_FIELDS,
    build_era5_cache_requests,
    fill_era5_cache,
)


def _mini_observations() -> pd.DataFrame:
    """Minimal observation table for cache-filler tests."""
    return pd.DataFrame(
        [
            {
                "observation_id": "obs_a",
                "image_date": "2021-08-17",
                "authoritative_lat": 31.29812,
                "authoritative_lng": 120.01688,
            },
            {
                "observation_id": "obs_b",
                "image_date": "2024-02-12",
                "authoritative_lat": 54.62757,
                "authoritative_lng": -6.47327,
            },
        ]
    )


class _FakeResponse:
    """Minimal requests-style response for downloader tests."""

    def __init__(self, payload: dict):
        self._payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict:
        return self._payload


class _FakeSession:
    """Session double that records requests and returns deterministic payloads."""

    def __init__(self) -> None:
        self.calls: list[dict] = []

    def get(self, url: str, params: dict, timeout: float) -> _FakeResponse:
        self.calls.append({"url": url, "params": dict(params), "timeout": timeout})
        start_date = pd.Timestamp(params["start_date"], tz="UTC")
        end_date = pd.Timestamp(params["end_date"], tz="UTC") + pd.Timedelta(hours=23)
        hourly_times = pd.date_range(start=start_date, end=end_date, freq="1h", tz="UTC")
        payload = {
            "latitude": params["latitude"],
            "longitude": params["longitude"],
            "generationtime_ms": 1.0,
            "utc_offset_seconds": 0,
            "timezone": "GMT",
            "timezone_abbreviation": "GMT",
            "elevation": 0.0,
            "hourly_units": {
                "time": "iso8601",
                **{field: "-" for field in ERA5_ARCHIVE_HOURLY_FIELDS},
            },
            "hourly": {
                "time": [stamp.strftime("%Y-%m-%dT%H:%M") for stamp in hourly_times],
                **{
                    field: [float(index) for index, _ in enumerate(hourly_times)]
                    for field in ERA5_ARCHIVE_HOURLY_FIELDS
                },
            },
        }
        return _FakeResponse(payload)


def test_build_era5_cache_requests_uses_expected_spinup_window():
    requests = build_era5_cache_requests(_mini_observations())

    assert len(requests) == 2
    assert requests[0].observation_id == "obs_a"
    assert requests[0].start_date == "2021-08-07"
    assert requests[0].end_date == "2021-08-17"
    assert requests[0].cache_file_name == "31p298_120p017_20210807_20210817.json"
    assert requests[1].cache_file_name == "54p628_6p473_20240202_20240212.json"


def test_fill_era5_cache_writes_expected_json_files(tmp_path: Path):
    session = _FakeSession()

    summary = fill_era5_cache(
        observations=_mini_observations(),
        cache_dir=tmp_path,
        session=session,
        sleep_seconds=0.0,
    )

    assert summary["n_total_requests"] == 2
    assert summary["n_downloaded"] == 2
    assert summary["n_skipped_existing"] == 0
    assert len(session.calls) == 2

    cache_file = tmp_path / "31p298_120p017_20210807_20210817.json"
    assert cache_file.exists()
    payload = json.loads(cache_file.read_text(encoding="utf-8"))
    assert payload["timezone"] == "GMT"
    assert set(payload["hourly"]).issuperset({"time", *ERA5_ARCHIVE_HOURLY_FIELDS})
    assert len(payload["hourly"]["time"]) == 264


def test_fill_era5_cache_skips_existing_files(tmp_path: Path):
    session = _FakeSession()
    fill_era5_cache(
        observations=_mini_observations(),
        cache_dir=tmp_path,
        session=session,
        sleep_seconds=0.0,
    )

    second_session = _FakeSession()
    summary = fill_era5_cache(
        observations=_mini_observations(),
        cache_dir=tmp_path,
        session=second_session,
        sleep_seconds=0.0,
    )

    assert summary["n_downloaded"] == 0
    assert summary["n_skipped_existing"] == 2
    assert second_session.calls == []
