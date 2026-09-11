import numpy as np
import pytest

from penindex.stats import fit_distribution, linear_fit


def test_linear_fit_perfect_line():
    x = np.array([0, 1, 2, 3, 4])
    y = 2 * x + 1  # slope=2, intercept=1

    result = linear_fit(x, y)

    assert result['slope'] == pytest.approx(2.0)
    assert result['intercept'] == pytest.approx(1.0)
    assert result['r2'] == pytest.approx(1.0)
    assert result['n'] == 5


def test_linear_fit_matches_pearson_r_squared():
    x = np.array([0, 1, 2, 3, 4])
    y = np.array([3, 1, 4, 1, 5])  # not perfectly correlated with x

    result = linear_fit(x, y)
    expected_r2 = np.corrcoef(x, y)[0, 1] ** 2

    assert result['r2'] == pytest.approx(expected_r2, abs=1e-9)


def test_linear_fit_requires_at_least_two_points():
    with pytest.raises(ValueError):
        linear_fit([1.0], [2.0])


def test_linear_fit_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        linear_fit([1.0, 2.0], [1.0, 2.0, 3.0])


def test_fit_distribution_recovers_gaussian_parameters():
    rng = np.random.default_rng(42)
    true_center, true_width = 50.0, 8.0
    values = rng.normal(loc=true_center, scale=true_width, size=5000)

    result = fit_distribution(values, distribution="gaussian", bins=40)

    assert result['center'] == pytest.approx(true_center, abs=1.0)
    assert result['width'] == pytest.approx(true_width, abs=1.5)


def test_fit_distribution_rejects_unknown_distribution():
    with pytest.raises(ValueError):
        fit_distribution([1, 2, 3, 4, 5], distribution="unknown")
