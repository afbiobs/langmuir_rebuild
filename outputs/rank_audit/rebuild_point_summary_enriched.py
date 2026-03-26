#!/usr/bin/env python3
"""
Rebuild data/raw/point_summary_enriched.csv from Copernicus annotation CSVs.

Workflow:
- read all `copernicus_annotations*.csv` files from LN_Bathy_GEM/Data
- keep only point annotations
- group by `group_id`
- aggregate point-level geometry and spacing statistics
- retain the point `category` so groups remain distinguishable as
  `manual`, `stream`, or `wiggle`
- join bathymetry and fetch distances from the existing
  `data/raw/point_summary_enriched.csv` by exact `group_id`
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _first_non_null(series: pd.Series):
    """Return the first non-null value from `series`, else NA."""
    non_null = series.dropna()
    if non_null.empty:
        return pd.NA
    return non_null.iloc[0]


def _unique_join(series: pd.Series) -> str:
    """Join unique non-null values into a stable `|`-separated string."""
    values = sorted({str(value) for value in series if pd.notna(value)})
    return "|".join(values)


def load_annotation_points(annotation_dir: Path) -> pd.DataFrame:
    """
    Load and concatenate all point annotations from the source CSV files.

    Parameters:
        annotation_dir: Directory containing `copernicus_annotations*.csv` [path]

    Returns:
        Concatenated point-annotation table [mixed]
    """
    files = sorted(annotation_dir.glob("copernicus_annotations*.csv"))
    if not files:
        raise FileNotFoundError(
            f"No copernicus annotation CSVs found under {annotation_dir}"
        )
    frames = []
    for path in files:
        frame = pd.read_csv(path)
        frame["source_file"] = path.name
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True, sort=False)
    points = combined.loc[combined["type"] == "point"].copy()
    if points.empty:
        raise ValueError("No point annotations found in source files.")
    return points


def aggregate_points(points: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate point annotations to one row per `group_id`.

    Parameters:
        points: Point annotation table [mixed]

    Returns:
        Grouped point summary with retained category [mixed]
    """
    category_counts = points.groupby("group_id")["category"].nunique(dropna=True)
    mixed = category_counts.loc[category_counts > 1]
    if not mixed.empty:
        raise ValueError(
            "Point groups contain mixed categories: "
            f"{mixed.index.tolist()}"
        )

    grouped = (
        points.groupby("group_id", dropna=False)
        .agg(
            category=("category", _first_non_null),
            label=("label", _first_non_null),
            observation_date=("observation_date", _unique_join),
            n_points=("id", "size"),
            lat=("lat", "mean"),
            lng=("lng", "mean"),
            lat1=("lat1", "mean"),
            lng1=("lng1", "mean"),
            lat2=("lat2", "mean"),
            lng2=("lng2", "mean"),
            length_m=("length_m", "mean"),
            bearing_deg=("bearing_deg", "mean"),
            group_index=("group_index", "mean"),
            dist_from_prev_m=("dist_from_prev_m", "mean"),
            meters_per_px=("meters_per_px", "mean"),
            zoom=("zoom", "mean"),
            map_center_lat=("map_center_lat", "mean"),
            map_center_lng=("map_center_lng", "mean"),
            first_timestamp=("timestamp", "min"),
            last_timestamp=("timestamp", "max"),
            source_files=("source_file", _unique_join),
        )
        .reset_index()
    )
    grouped["category"] = grouped["category"].astype("string")
    return grouped


def merge_fetch_fields(
    grouped: pd.DataFrame,
    current_enriched_path: Path,
) -> pd.DataFrame:
    """
    Join bathymetry/fetch fields from the existing enriched point summary.

    Parameters:
        grouped: Grouped annotation summary [mixed]
        current_enriched_path: Existing enriched point summary [path]

    Returns:
        Rebuilt enriched point summary [mixed]
    """
    current = pd.read_csv(current_enriched_path)
    fetch_columns = [
        "bathymetry",
        "fetch_0",
        "fetch_20",
        "fetch_40",
        "fetch_60",
        "fetch_80",
        "fetch_100",
        "fetch_120",
        "fetch_140",
        "fetch_160",
        "fetch_180",
        "fetch_200",
        "fetch_220",
        "fetch_240",
        "fetch_260",
        "fetch_280",
        "fetch_300",
        "fetch_320",
        "fetch_340",
        "x",
        "y",
        "dist_to_grid",
    ]
    join_frame = current[["group_id", *fetch_columns]].copy()

    grouped_ids = set(grouped["group_id"])
    current_ids = set(join_frame["group_id"])
    if grouped_ids != current_ids:
        missing_in_current = sorted(grouped_ids - current_ids)
        missing_in_grouped = sorted(current_ids - grouped_ids)
        raise ValueError(
            "Group-id mismatch between annotations and current enriched table. "
            f"Missing in current: {missing_in_current[:10]}; "
            f"missing in grouped: {missing_in_grouped[:10]}"
        )

    rebuilt = grouped.merge(join_frame, on="group_id", how="left", validate="one_to_one")
    return rebuilt


def rebuild_point_summary(annotation_dir: Path, current_enriched_path: Path) -> pd.DataFrame:
    """Run the full point-summary rebuild pipeline."""
    points = load_annotation_points(annotation_dir)
    grouped = aggregate_points(points)
    rebuilt = merge_fetch_fields(grouped, current_enriched_path)
    column_order = [
        "group_id",
        "category",
        "label",
        "observation_date",
        "n_points",
        "lat",
        "lng",
        "lat1",
        "lng1",
        "lat2",
        "lng2",
        "length_m",
        "bearing_deg",
        "group_index",
        "dist_from_prev_m",
        "meters_per_px",
        "zoom",
        "map_center_lat",
        "map_center_lng",
        "first_timestamp",
        "last_timestamp",
        "source_files",
        "bathymetry",
        "fetch_0",
        "fetch_20",
        "fetch_40",
        "fetch_60",
        "fetch_80",
        "fetch_100",
        "fetch_120",
        "fetch_140",
        "fetch_160",
        "fetch_180",
        "fetch_200",
        "fetch_220",
        "fetch_240",
        "fetch_260",
        "fetch_280",
        "fetch_300",
        "fetch_320",
        "fetch_340",
        "x",
        "y",
        "dist_to_grid",
    ]
    return rebuilt[column_order].sort_values("group_id").reset_index(drop=True)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for the rebuild script."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--annotation-dir",
        type=Path,
        default=Path("/home/op/PROJECTS/LN_Bathy_GEM/Data"),
        help="Directory containing copernicus annotation CSVs [path]",
    )
    parser.add_argument(
        "--current-enriched",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "point_summary_enriched.csv",
        help="Existing enriched point summary used as the fetch/bathymetry source [path]",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=REPO_ROOT / "data" / "raw" / "point_summary_enriched.csv",
        help="Destination for the rebuilt enriched point summary [path]",
    )
    return parser.parse_args()


def main() -> None:
    """Rebuild and write the enriched point summary."""
    args = parse_args()
    rebuilt = rebuild_point_summary(
        annotation_dir=args.annotation_dir,
        current_enriched_path=args.current_enriched,
    )
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    rebuilt.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(rebuilt)} rows to {args.output_csv}")
    print(rebuilt[["category", "n_points"]].groupby("category").agg(n_groups=("n_points", "size"), median_n_points=("n_points", "median")).to_string())


if __name__ == "__main__":
    main()
