"""
Reusable ERA5/Open-Meteo cache builder for observation-matched weather history.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import timedelta
import argparse
import json
from pathlib import Path
import time
from typing import Any

import pandas as pd
import requests


ERA5_ARCHIVE_URL = "https://archive-api.open-meteo.com/v1/archive"

ERA5_ARCHIVE_HOURLY_FIELDS = (
    "wind_speed_10m",
    "wind_direction_10m",
    "wind_gusts_10m",
    "surface_pressure",
    "temperature_2m",
    "shortwave_radiation",
    "cloud_cover",
    "precipitation",
)


@dataclass(frozen=True)
class Era5CacheRequest:
    """
    One hourly weather-history request for the reusable cache.

    Parameters:
        observation_id: Observation identifier [-]
        latitude_deg:   Observation latitude [deg]
        longitude_deg:  Observation longitude [deg]
        start_date:     Inclusive history-window start date [YYYY-MM-DD]
        end_date:       Inclusive history-window end date [YYYY-MM-DD]
        cache_key:      Stable cache key shared with the comparison harness [-]
        cache_file_name: Output JSON file name [-]
    """

    observation_id: str
    latitude_deg: float
    longitude_deg: float
    start_date: str
    end_date: str
    cache_key: str
    cache_file_name: str


def era5_cache_key(
    latitude_deg: float,
    longitude_deg: float,
    start_date: str,
    end_date: str,
) -> str:
    """
    Build the stable weather-cache key used by the comparison harness.

    Parameters:
        latitude_deg:  Latitude [deg]
        longitude_deg: Longitude [deg]
        start_date:    Inclusive history-window start [YYYY-MM-DD]
        end_date:      Inclusive history-window end [YYYY-MM-DD]

    Returns:
        Cache key string [-]
    """
    return (
        f"{latitude_deg:.3f}_{longitude_deg:.3f}_{start_date}_{end_date}"
        .replace("-", "")
        .replace(".", "p")
    )


def _observation_frame(
    observations: str | Path | pd.DataFrame,
) -> pd.DataFrame:
    """Load the raw observation table required to build cache requests."""
    if isinstance(observations, pd.DataFrame):
        frame = observations.copy()
    else:
        frame = pd.read_csv(observations)

    required = {"observation_id", "image_date", "authoritative_lat", "authoritative_lng"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(
            "observations are missing required columns: "
            + ", ".join(sorted(missing))
        )
    return frame.reset_index(drop=True)


def build_era5_cache_requests(
    observations: str | Path | pd.DataFrame = "data/raw/observations.csv",
    spinup_days: int = 10,
) -> list[Era5CacheRequest]:
    """
    Build the unique cache-request set implied by the observation table.

    Parameters:
        observations: Raw observation table or CSV path [-]
        spinup_days:  Number of days of weather history to include before the
                      observation date [days]

    Returns:
        Ordered list of unique cache requests [-]
    """
    if spinup_days <= 0:
        raise ValueError(f"spinup_days must be positive, got {spinup_days}.")

    frame = _observation_frame(observations)
    requests_by_key: dict[str, Era5CacheRequest] = {}
    for _, row in frame.iterrows():
        image_date = pd.Timestamp(row["image_date"]).date()
        start_date = (image_date - timedelta(days=spinup_days)).isoformat()
        end_date = image_date.isoformat()
        cache_key = era5_cache_key(
            float(row["authoritative_lat"]),
            float(row["authoritative_lng"]),
            start_date,
            end_date,
        )
        requests_by_key.setdefault(
            cache_key,
            Era5CacheRequest(
                observation_id=str(row["observation_id"]),
                latitude_deg=float(row["authoritative_lat"]),
                longitude_deg=float(row["authoritative_lng"]),
                start_date=start_date,
                end_date=end_date,
                cache_key=cache_key,
                cache_file_name=f"{cache_key}.json",
            ),
        )
    return list(requests_by_key.values())


def _archive_request_params(request: Era5CacheRequest) -> dict[str, Any]:
    """Build the Open-Meteo archive query for one request."""
    return {
        "latitude": request.latitude_deg,
        "longitude": request.longitude_deg,
        "start_date": request.start_date,
        "end_date": request.end_date,
        "hourly": ",".join(ERA5_ARCHIVE_HOURLY_FIELDS),
        "timezone": "GMT",
        "models": "era5",
        "wind_speed_unit": "ms",
        "precipitation_unit": "mm",
    }


def _validate_archive_payload(payload: dict[str, Any], request: Era5CacheRequest) -> None:
    """Validate the minimum JSON structure expected by the cache consumer."""
    if "hourly" not in payload:
        raise ValueError(f"{request.cache_file_name}: payload is missing 'hourly'.")
    hourly = payload["hourly"]
    required_fields = {"time", *ERA5_ARCHIVE_HOURLY_FIELDS}
    missing = required_fields.difference(hourly.keys())
    if missing:
        raise ValueError(
            f"{request.cache_file_name}: payload is missing hourly fields "
            + ", ".join(sorted(missing))
        )
    n_times = len(hourly["time"])
    if n_times == 0:
        raise ValueError(f"{request.cache_file_name}: payload contains no hourly rows.")
    bad_lengths = [
        field for field in required_fields if len(hourly[field]) != n_times
    ]
    if bad_lengths:
        raise ValueError(
            f"{request.cache_file_name}: hourly field lengths do not match time for "
            + ", ".join(sorted(bad_lengths))
        )


def fetch_era5_archive_payload(
    request: Era5CacheRequest,
    *,
    session: requests.Session | Any | None = None,
    timeout_seconds: float = 90.0,
    retries: int = 3,
    backoff_seconds: float = 1.0,
) -> dict[str, Any]:
    """
    Fetch one weather-history payload from the Open-Meteo archive API.

    Parameters:
        request:         Cache request metadata [-]
        session:         Optional requests-compatible session [-]
        timeout_seconds: Per-request timeout [s]
        retries:         Number of request attempts [-]
        backoff_seconds: Base retry backoff [s]

    Returns:
        JSON payload compatible with the comparison cache loader [-]
    """
    if retries <= 0:
        raise ValueError(f"retries must be positive, got {retries}.")
    if timeout_seconds <= 0.0:
        raise ValueError(f"timeout_seconds must be positive, got {timeout_seconds}.")

    http = session if session is not None else requests.Session()
    if hasattr(http, "headers"):
        http.headers.setdefault(
            "User-Agent",
            "langmuir-rebuild-era5-cache/1.0",
        )

    last_error: Exception | None = None
    for attempt_index in range(retries):
        try:
            response = http.get(
                ERA5_ARCHIVE_URL,
                params=_archive_request_params(request),
                timeout=timeout_seconds,
            )
            response.raise_for_status()
            payload = response.json()
            _validate_archive_payload(payload, request)
            return payload
        except Exception as exc:  # pragma: no cover - exercised via retry path
            last_error = exc
            if attempt_index == retries - 1:
                break
            time.sleep(backoff_seconds * float(2**attempt_index))

    assert last_error is not None
    raise RuntimeError(
        f"Failed to fetch {request.cache_file_name} after {retries} attempt(s)."
    ) from last_error


def fill_era5_cache(
    observations: str | Path | pd.DataFrame = "data/raw/observations.csv",
    cache_dir: str | Path = "data/raw/era5_cache",
    *,
    spinup_days: int = 10,
    overwrite: bool = False,
    session: requests.Session | Any | None = None,
    sleep_seconds: float = 0.25,
    timeout_seconds: float = 90.0,
    retries: int = 3,
    backoff_seconds: float = 1.0,
) -> dict[str, Any]:
    """
    Populate the reusable observation-matched weather cache on disk.

    Parameters:
        observations:     Raw observation table or CSV path [-]
        cache_dir:        Destination directory for cache files [-]
        spinup_days:      Number of days to fetch before each image date [days]
        overwrite:        If True, re-download files even if they exist [-]
        session:          Optional requests-compatible session [-]
        sleep_seconds:    Delay between successful downloads [s]
        timeout_seconds:  Per-request timeout [s]
        retries:          Number of request attempts per file [-]
        backoff_seconds:  Base retry backoff [s]

    Returns:
        Summary dict with counts, file paths, and any failures [-]
    """
    if sleep_seconds < 0.0:
        raise ValueError(f"sleep_seconds must be non-negative, got {sleep_seconds}.")

    requests_to_fill = build_era5_cache_requests(observations, spinup_days=spinup_days)
    output_dir = Path(cache_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    downloaded_files: list[str] = []
    skipped_existing: list[str] = []
    failed_requests: list[dict[str, str]] = []

    for request in requests_to_fill:
        cache_file = output_dir / request.cache_file_name
        if cache_file.exists() and not overwrite:
            skipped_existing.append(str(cache_file))
            continue

        try:
            payload = fetch_era5_archive_payload(
                request,
                session=session,
                timeout_seconds=timeout_seconds,
                retries=retries,
                backoff_seconds=backoff_seconds,
            )
            cache_file.write_text(json.dumps(payload), encoding="utf-8")
            downloaded_files.append(str(cache_file))
            if sleep_seconds > 0.0:
                time.sleep(sleep_seconds)
        except Exception as exc:
            failed_requests.append(
                {
                    "cache_file": str(cache_file),
                    "observation_id": request.observation_id,
                    "error": str(exc),
                }
            )

    return {
        "cache_dir": str(output_dir),
        "n_total_requests": len(requests_to_fill),
        "n_downloaded": len(downloaded_files),
        "n_skipped_existing": len(skipped_existing),
        "n_failed": len(failed_requests),
        "downloaded_files": downloaded_files,
        "skipped_existing": skipped_existing,
        "failed_requests": failed_requests,
    }


def _build_arg_parser() -> argparse.ArgumentParser:
    """Create the CLI argument parser for the cache filler."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--observations",
        default="data/raw/observations.csv",
        help="Observation CSV path.",
    )
    parser.add_argument(
        "--cache-dir",
        default="data/raw/era5_cache",
        help="Directory where cache JSON files will be written.",
    )
    parser.add_argument(
        "--spinup-days",
        type=int,
        default=10,
        help="Number of days of history to include before each image date.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download existing cache files.",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=0.25,
        help="Delay between successful downloads in seconds.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=90.0,
        help="Per-request timeout in seconds.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of request attempts per file.",
    )
    parser.add_argument(
        "--backoff-seconds",
        type=float,
        default=1.0,
        help="Base exponential backoff in seconds.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the cache filler as a CLI."""
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    summary = fill_era5_cache(
        observations=args.observations,
        cache_dir=args.cache_dir,
        spinup_days=args.spinup_days,
        overwrite=bool(args.overwrite),
        sleep_seconds=float(args.sleep_seconds),
        timeout_seconds=float(args.timeout_seconds),
        retries=int(args.retries),
        backoff_seconds=float(args.backoff_seconds),
    )
    print(json.dumps(summary, indent=2))
    return 0 if summary["n_failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
