"""
Tests for the structural evaluation metrics used by the WP-05 comparison harness.
"""

from __future__ import annotations

import pytest

from src.evaluation.metrics import (
    hit_rate_within_fraction,
    pearson_correlation,
    range_coverage,
    spacing_bias,
    spacing_mae_ratio,
    spacing_rmse_ratio,
    spearman_correlation,
)


def test_range_coverage_matches_fraction_of_observed_span():
    predicted = [40.0, 70.0, 100.0]
    observed_range = (20.0, 120.0)
    assert range_coverage(predicted, observed_range) == pytest.approx(0.60)


def test_relative_error_metrics_and_bias_are_consistent():
    observed = [50.0, 100.0]
    predicted = [60.0, 80.0]

    assert spacing_bias(predicted, observed) == pytest.approx(-5.0)
    assert spacing_mae_ratio(predicted, observed) == pytest.approx(20.0)
    assert spacing_rmse_ratio(predicted, observed) == pytest.approx(20.0)


def test_hit_rate_and_correlations_reach_unity_for_perfect_ordered_fit():
    observed = [40.0, 80.0, 120.0, 160.0]
    predicted = [40.0, 80.0, 120.0, 160.0]

    assert hit_rate_within_fraction(predicted, observed) == 1.0
    assert pearson_correlation(predicted, observed) == pytest.approx(1.0)
    assert spearman_correlation(predicted, observed) == pytest.approx(1.0)
