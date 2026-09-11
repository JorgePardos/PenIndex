"""Publication figures.

plot_descriptor_correlation reproduces the flagship figure type of Pardos,
Gonzalo, Merino & Echeverria, Dalton Trans. 2026 (their Fig. 2-7): a
two-panel scatter, distance-vs-descriptor next to PI-vs-descriptor, each
with an ordinary-least-squares fit line and R^2 annotation - the exact
comparison used throughout that paper to show PI is a better covalency
predictor than raw distance. plot_r2_summary reproduces their Fig. 8,
a bar-chart comparison of R^2(distance) vs R^2(PI) across datasets.

plot_pi_histogram and plot_pi_timeseries reproduce the other figure type
used throughout Echeverria & Alvarez, Chem. Sci. 2023: distribution
histograms with a peak/width fit overlay (e.g. Fig. 14c, 30), and PI
evolving along a trajectory/optimization path (relevant to reactive_mode).
"""

from typing import Any, Dict, Optional, Sequence

import numpy as np
import matplotlib
matplotlib.use("Agg")  # figures are written to file; no display needed
import matplotlib.pyplot as plt

from penindex.stats import fit_distribution, linear_fit

# Colorblind-safe palette (Okabe-Ito).
_COLOR_DISTANCE = "#0072B2"
_COLOR_PI = "#D55E00"


def plot_descriptor_correlation(
    distance: Sequence[float],
    pi: Sequence[float],
    descriptor: Sequence[float],
    descriptor_label: str,
    out_file: str,
    distance_label: str = "Distance (Å)",
    pi_label: str = "Penetration Index (%)",
) -> Dict[str, Dict[str, float]]:
    """Two-panel scatter: distance-vs-descriptor (left) and PI-vs-descriptor
    (right), each with a fitted line and R^2 annotation.

    Saves a figure to out_file (format chosen by its extension - use .pdf
    or .svg for a vector file suitable for a manuscript, .png for a raster
    preview) and returns the two linear_fit() results as
    {'distance': {...}, 'pi': {...}}.
    """
    fit_distance = linear_fit(distance, descriptor)
    fit_pi = linear_fit(pi, descriptor)

    fig, (ax_dist, ax_pi) = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)

    for ax, x, fit, color, xlabel in (
        (ax_dist, distance, fit_distance, _COLOR_DISTANCE, distance_label),
        (ax_pi, pi, fit_pi, _COLOR_PI, pi_label),
    ):
        x = np.asarray(x, dtype=float)
        ax.scatter(x, descriptor, color=color, edgecolor="black", linewidth=0.5, s=40, zorder=3)

        x_line = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_line, fit['slope'] * x_line + fit['intercept'], color=color, linewidth=1.5, zorder=2)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(descriptor_label)
        ax.text(0.05, 0.95, f"$R^2$ = {fit['r2']:.2f}", transform=ax.transAxes,
                va='top', ha='left', fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.savefig(out_file, dpi=300)
    plt.close(fig)

    return {'distance': fit_distance, 'pi': fit_pi}


def plot_r2_summary(
    dataset_labels: Sequence[str],
    r2_distance: Sequence[float],
    r2_pi: Sequence[float],
    out_file: str,
) -> None:
    """Bar chart comparing R^2(distance) vs R^2(PI) across several datasets -
    the same comparison as Fig. 8 of Pardos et al., Dalton Trans. 2026."""
    n = len(dataset_labels)
    x = np.arange(n)
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(6, n * 0.9), 4), constrained_layout=True)
    ax.bar(x - width / 2, r2_distance, width, label="Distance", color=_COLOR_DISTANCE)
    ax.bar(x + width / 2, r2_pi, width, label="Penetration Index", color=_COLOR_PI)
    ax.set_xticks(x)
    ax.set_xticklabels(dataset_labels, rotation=30, ha='right')
    ax.set_ylabel("$R^2$")
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.savefig(out_file, dpi=300)
    plt.close(fig)


def plot_pi_histogram(
    values: Sequence[float],
    out_file: str,
    bins: int = 30,
    distribution_fit: Optional[str] = None,
    pi_label: str = "Penetration Index (%)",
) -> Optional[Dict[str, Any]]:
    """Histogram of PI values, with an optional Gaussian/Lorentzian fit
    overlay (distribution_fit='gaussian'|'lorentzian') to extract a peak
    position and width - reproduces the distribution-characterization
    figures of Echeverria & Alvarez, Chem. Sci. 2023 (e.g. Fig. 14c, 30).

    Returns the penindex.stats.fit_distribution() result if a fit was
    requested, else None.
    """
    values = np.asarray(values, dtype=float)
    fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)

    fit_result = None
    if distribution_fit is not None:
        fit_result = fit_distribution(values, distribution=distribution_fit, bins=bins)
        bin_width = fit_result['bin_centers'][1] - fit_result['bin_centers'][0]
        ax.bar(fit_result['bin_centers'], fit_result['bin_counts'], width=bin_width,
               color=_COLOR_PI, alpha=0.6, edgecolor='black', linewidth=0.3, zorder=2)

        x_fit = np.linspace(values.min(), values.max(), 300)
        center, width_param = fit_result['center'], fit_result['width']
        if distribution_fit == 'gaussian':
            y_fit = fit_result['amplitude'] * np.exp(-0.5 * ((x_fit - center) / width_param) ** 2)
        else:
            y_fit = fit_result['amplitude'] * (width_param ** 2) / ((x_fit - center) ** 2 + width_param ** 2)

        ax.plot(x_fit, y_fit, color='black', linewidth=1.5, zorder=3,
                label=f"{distribution_fit.capitalize()} fit\ncenter={center:.1f}%, width={width_param:.1f}%")
        ax.legend(frameon=False, fontsize=8)
    else:
        ax.hist(values, bins=bins, color=_COLOR_PI, alpha=0.8, edgecolor='black', linewidth=0.3)

    ax.set_xlabel(pi_label)
    ax.set_ylabel("Count")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.savefig(out_file, dpi=300)
    plt.close(fig)

    return fit_result


def plot_pi_timeseries(
    frames: Sequence[int],
    pi_values: Sequence[float],
    out_file: str,
    label: Optional[str] = None,
    pi_label: str = "Penetration Index (%)",
) -> None:
    """Plots PI vs. frame/step for one tracked interaction across a
    trajectory or QM optimization/IRC path - useful to visualize a bond
    breaking/forming event (reactive_mode) alongside the discrete event
    table from penindex.stats.aggregate_*_reactive."""
    frames = np.asarray(frames)
    pi_values = np.asarray(pi_values, dtype=float)

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.plot(frames, pi_values, color=_COLOR_PI, linewidth=1.5, marker='o', markersize=3)
    ax.set_xlabel("Frame / step")
    ax.set_ylabel(pi_label)
    if label:
        ax.set_title(label, fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.savefig(out_file, dpi=300)
    plt.close(fig)
