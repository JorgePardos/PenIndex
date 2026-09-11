"""Batch single-pair mode: computes the Penetration Index for one specific,
user-identified bond in each of many small systems (diatomics, phosphine
chalcogenides, hydrogen-bonded complexes, halonium ions, ...).

This is the workflow behind Pardos, Gonzalo, Merino & Echeverria, Dalton
Trans. 2026: a curated dataset of small molecules/ions, each reduced to a
single PI value for the bond of interest, then correlated (see
penindex.stats.linear_fit and penindex.viz.plots) against an externally
computed covalency descriptor (penindex.io.descriptors) for the same
systems. This is a different shape of problem from penindex.core's
structural scanning (one specific bond per small system, rather than every
contact within one large structure), so it gets its own entry points here
rather than being folded into detect_generic_interactions.
"""

import warnings
from dataclasses import dataclass
from typing import Any, Dict, List

import numpy as np
import pandas as pd

from penindex import core
from penindex.exceptions import PenIndexError
from penindex.io.loaders import build_universe


@dataclass
class BatchSystem:
    """One entry in a batch run: a structure file plus two MDAnalysis atom
    selectors identifying the single bond of interest in it. Each selector
    must resolve to exactly one atom."""
    label: str
    structure: str
    atom_1: str
    atom_2: str
    dynamic_carbon_radii: bool = False


def load_batch_systems_csv(path: str) -> List[BatchSystem]:
    """Loads a batch systems list from a CSV with columns: label, structure,
    atom_1, atom_2, and an optional dynamic_carbon_radii (true/false)."""
    try:
        df = pd.read_csv(path)
    except Exception as e:
        raise PenIndexError(f"Could not read batch systems CSV '{path}': {e}") from e

    required = {'label', 'structure', 'atom_1', 'atom_2'}
    missing = required - set(df.columns)
    if missing:
        raise PenIndexError(f"Batch systems CSV '{path}' is missing columns: {sorted(missing)}")

    has_dynamic_col = 'dynamic_carbon_radii' in df.columns
    systems = []
    for _, row in df.iterrows():
        dynamic = bool(row['dynamic_carbon_radii']) if has_dynamic_col else False
        systems.append(BatchSystem(row['label'], row['structure'], row['atom_1'], row['atom_2'], dynamic))
    return systems


def run_batch(systems: List[BatchSystem]) -> pd.DataFrame:
    """Computes distance and PI for the specified bond in each system.

    Returns one row per system (label, structure, Distance_AB, PI(%)).
    Systems that fail to load, or whose atom_1/atom_2 selectors don't
    resolve to exactly one atom each, are skipped with a warning rather
    than aborting the whole batch - a single malformed entry in a dataset
    of dozens of small systems shouldn't lose the rest of the results.
    """
    rows = []
    for system in systems:
        try:
            rows.append(_process_one_system(system))
        except PenIndexError as e:
            warnings.warn(f"Skipping '{system.label}': {e}")

    return pd.DataFrame(rows, columns=['label', 'structure', 'Distance_AB', 'PI(%)'])


def _process_one_system(system: BatchSystem) -> Dict[str, Any]:
    config = {'system': {'topology': system.structure,
                          'dynamic_carbon_radii': system.dynamic_carbon_radii}}
    universe = build_universe(config)

    cov_radii, vdw_radii, _ = core.precompute_radii_arrays(
        universe,
        dynamic_carbon_radii=system.dynamic_carbon_radii,
        topology_path=system.structure,
    )

    sel1 = universe.select_atoms(system.atom_1)
    sel2 = universe.select_atoms(system.atom_2)
    if len(sel1) != 1 or len(sel2) != 1:
        raise PenIndexError(
            f"atom_1/atom_2 selectors must each resolve to exactly one atom "
            f"(got {len(sel1)} and {len(sel2)} for '{system.label}')"
        )

    idx1, idx2 = sel1.indices[0], sel2.indices[0]
    distance = float(np.linalg.norm(sel1.positions[0] - sel2.positions[0]))
    pi = core.compute_pi(distance, cov_radii[idx1], vdw_radii[idx1], cov_radii[idx2], vdw_radii[idx2])

    return {'label': system.label, 'structure': system.structure, 'Distance_AB': distance, 'PI(%)': pi}
