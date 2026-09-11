"""Aggregation of per-frame/per-interaction PI data into summary statistics.

Every value reported in the source literature (Echeverria & Alvarez, Chem.
Sci. 2023) is a mean with its standard deviation, e.g. "47(7)% for the
H...N pair" - a CSD survey average with its spread. Std_PI below is the
per-pair equivalent of that spread across the frames/structures where the
interaction was detected: pandas' default sample standard deviation (ddof=1),
which is NaN for pairs detected in only one frame (undefined for n=1), same
convention the literature itself uses for a single measurement.
"""

from typing import Any, Dict, Sequence

import numpy as np
import pandas as pd

from penindex.exceptions import PenIndexError


def aggregate_hbond_standard(df: pd.DataFrame, total_frames: int) -> pd.DataFrame:
    """Aggregates hydrogen-bond PI data for stable interactions in classical MD."""
    return df.groupby(['Donor_Idx', 'Hydrogen_Idx', 'Acceptor_Idx']).agg(
        Occupancy=('Frame', lambda x: (len(x) / total_frames) * 100),
        Avg_Dist_HA=('Distance_HA', 'mean'),
        Avg_PI=('PI(%)', 'mean'),
        Std_PI=('PI(%)', 'std'),
    ).reset_index().sort_values(by='Occupancy', ascending=False)


def aggregate_hbond_reactive(df: pd.DataFrame, step_size: int) -> pd.DataFrame:
    """Detects and records independent hydrogen-bond breaking/forming events."""
    df = df.sort_values(by=['Donor_Idx', 'Hydrogen_Idx', 'Acceptor_Idx', 'Frame']).reset_index(drop=True)

    df['Prev_Frame'] = df.groupby(['Donor_Idx', 'Hydrogen_Idx', 'Acceptor_Idx'])['Frame'].shift(1)
    df['New_Event'] = (df['Frame'] - df['Prev_Frame'] > step_size) | (df['Prev_Frame'].isna())
    df['Event_ID'] = df.groupby(['Donor_Idx', 'Hydrogen_Idx', 'Acceptor_Idx'])['New_Event'].cumsum()

    return df.groupby(['Donor_Idx', 'Hydrogen_Idx', 'Acceptor_Idx', 'Event_ID']).agg(
        Formed_At_Frame=('Frame', 'min'),
        Broken_At_Frame=('Frame', 'max'),
        Duration_Frames=('Frame', lambda x: len(x) * step_size),
        Max_Penetration_Index=('PI(%)', 'max'),
        Avg_Distance_HA=('Distance_HA', 'mean')
    ).reset_index().drop(columns=['Event_ID'])


def aggregate_generic_standard(df: pd.DataFrame, total_frames: int) -> pd.DataFrame:
    """Aggregates generic atom-pair PI data for stable/repeated contacts.

    If df carries the optional P_Atom1(%)/P_Atom2(%)/Crust_Width_Ratio
    columns (core.detect_generic_interactions with
    include_individual_penetrations=True), their per-pair averages are
    carried through too; Crust_Width_Ratio is constant per pair (it only
    depends on the two elements' radii, not the per-frame distance), so
    averaging it is equivalent to just picking it out.
    """
    agg_spec = {
        'Occupancy': ('Frame', lambda x: (len(x) / total_frames) * 100),
        'Avg_Distance_AB': ('Distance_AB', 'mean'),
        'Avg_PI': ('PI(%)', 'mean'),
        'Std_PI': ('PI(%)', 'std'),
    }
    if 'P_Atom1(%)' in df.columns:
        agg_spec['Avg_P_Atom1'] = ('P_Atom1(%)', 'mean')
        agg_spec['Avg_P_Atom2'] = ('P_Atom2(%)', 'mean')
        agg_spec['Crust_Width_Ratio'] = ('Crust_Width_Ratio', 'mean')

    return (df.groupby(['Atom1_Idx', 'Atom2_Idx']).agg(**agg_spec)
            .reset_index().sort_values(by='Avg_PI', ascending=False))


def aggregate_generic_reactive(df: pd.DataFrame, step_size: int) -> pd.DataFrame:
    """Detects and records independent breaking/forming events for generic
    atom-pair interactions (e.g. covalent bonds forming/breaking in a
    reactive QM/MM trajectory)."""
    df = df.sort_values(by=['Atom1_Idx', 'Atom2_Idx', 'Frame']).reset_index(drop=True)

    df['Prev_Frame'] = df.groupby(['Atom1_Idx', 'Atom2_Idx'])['Frame'].shift(1)
    df['New_Event'] = (df['Frame'] - df['Prev_Frame'] > step_size) | (df['Prev_Frame'].isna())
    df['Event_ID'] = df.groupby(['Atom1_Idx', 'Atom2_Idx'])['New_Event'].cumsum()

    return df.groupby(['Atom1_Idx', 'Atom2_Idx', 'Event_ID']).agg(
        Formed_At_Frame=('Frame', 'min'),
        Broken_At_Frame=('Frame', 'max'),
        Duration_Frames=('Frame', lambda x: len(x) * step_size),
        Max_Penetration_Index=('PI(%)', 'max'),
        Avg_Distance_AB=('Distance_AB', 'mean')
    ).reset_index().drop(columns=['Event_ID'])


def linear_fit(x: Sequence[float], y: Sequence[float]) -> Dict[str, float]:
    """Ordinary least-squares linear regression y = slope*x + intercept.

    Returns slope, intercept, r2 (coefficient of determination, the square
    of the Pearson correlation coefficient) and n (number of points). r2 is
    the exact quantity used throughout the penetration-index literature
    (e.g. Pardos, Gonzalo, Merino & Echeverria, Dalton Trans. 2026) to
    compare how well raw distance vs. PI predicts an external covalency
    descriptor - the core "does PI beat distance" comparison this project
    exists to reproduce.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) != len(y):
        raise ValueError(f"x and y must have the same length (got {len(x)} and {len(y)})")
    if len(x) < 2:
        raise ValueError("linear_fit needs at least 2 points")

    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')

    return {'slope': float(slope), 'intercept': float(intercept), 'r2': float(r2), 'n': len(x)}


def fit_distribution(values: Sequence[float], distribution: str = "gaussian", bins: int = 30) -> Dict[str, Any]:
    """Fits a Gaussian or Lorentzian curve to the histogram of `values`,
    extracting a peak position and width - the same kind of distribution
    characterization used in the source literature to summarize a class of
    interactions (Echeverria & Alvarez, Chem. Sci. 2023, e.g. Fig. 14c/30:
    "Fitting of the M/M penetration index distribution").

    Returns a dict with the fitted 'center' and 'width' (plus 'amplitude',
    'distribution', and the histogram's 'bin_centers'/'bin_counts' used for
    the fit, so a caller can overlay both on a plot without recomputing the
    histogram).
    """
    from scipy.optimize import curve_fit

    if distribution == "gaussian":
        def model(x, amplitude, center, width):
            return amplitude * np.exp(-0.5 * ((x - center) / width) ** 2)
    elif distribution == "lorentzian":
        def model(x, amplitude, center, width):
            return amplitude * (width ** 2) / ((x - center) ** 2 + width ** 2)
    else:
        raise ValueError(f"Unknown distribution '{distribution}'; use 'gaussian' or 'lorentzian'")

    values = np.asarray(values, dtype=float)
    counts, edges = np.histogram(values, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2

    # Reasonable initial guesses: amplitude = tallest bin, center = its
    # location, width = the sample's own standard deviation.
    p0 = [float(counts.max()), float(centers[np.argmax(counts)]), float(np.std(values)) or 1.0]
    try:
        popt, _ = curve_fit(model, centers, counts, p0=p0, maxfev=10000)
    except RuntimeError as e:
        raise PenIndexError(f"Distribution fit did not converge: {e}") from e

    amplitude, center, width = popt
    return {
        'distribution': distribution,
        'amplitude': float(amplitude),
        'center': float(center),
        'width': float(abs(width)),
        'bin_centers': centers,
        'bin_counts': counts,
    }
