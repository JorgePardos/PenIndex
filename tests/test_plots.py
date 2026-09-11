import numpy as np
import pytest

from penindex.viz.plots import (
    plot_descriptor_correlation,
    plot_pi_histogram,
    plot_pi_timeseries,
    plot_r2_summary,
)


def test_plot_descriptor_correlation_produces_file_and_fits(tmp_path):
    distance = [1.0, 1.5, 2.0, 2.5, 3.0]
    pi = [90.0, 70.0, 50.0, 30.0, 10.0]
    descriptor = [10.0, 8.0, 6.0, 4.0, 2.0]  # perfectly linear in both distance and pi here

    out_file = tmp_path / "correlation.png"
    fits = plot_descriptor_correlation(distance, pi, descriptor, "Descriptor (a.u.)", str(out_file))

    assert out_file.exists()
    assert fits['pi']['r2'] > 0.99
    assert fits['distance']['r2'] > 0.99


def test_plot_r2_summary_produces_file(tmp_path):
    out_file = tmp_path / "summary.png"
    plot_r2_summary(["X2", "HX"], [0.5, 0.6], [0.95, 0.98], str(out_file))
    assert out_file.exists()


def test_plot_pi_histogram_without_fit(tmp_path):
    values = [10, 20, 30, 40, 50, 20, 30, 40, 30]
    out_file = tmp_path / "hist.png"

    result = plot_pi_histogram(values, str(out_file))

    assert out_file.exists()
    assert result is None


def test_plot_pi_histogram_with_gaussian_fit(tmp_path):
    rng = np.random.default_rng(0)
    values = rng.normal(loc=60.0, scale=5.0, size=2000)
    out_file = tmp_path / "hist_fit.png"

    result = plot_pi_histogram(values, str(out_file), distribution_fit="gaussian")

    assert out_file.exists()
    assert result is not None
    assert result['center'] == pytest.approx(60.0, abs=2.0)


def test_plot_pi_timeseries_produces_file(tmp_path):
    frames = list(range(10))
    pi_values = [100 - i * 5 for i in frames]
    out_file = tmp_path / "timeseries.png"

    plot_pi_timeseries(frames, pi_values, str(out_file), label="Test bond")

    assert out_file.exists()
