"""The Penetration Index engine: radii setup and the two interaction
detectors (hydrogen-bond triads and generic atom-pair scanning). Pure
computation - every function here returns numpy arrays or DataFrames, does
no printing of results, and writes no files. Diagnostic/progress messages
about what was detected during setup (dynamic carbon radii) are kept as
prints for interactive-CLI parity with the previous single-file version.
"""

import warnings
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from scipy.spatial import cKDTree

from penindex.constants import RADII, D_H_BONDS, CARBON_HYBRID_COV_RADII


def compute_pi(distance: float, cov_radius_a: float, vdw_radius_a: float,
                cov_radius_b: float, vdw_radius_b: float) -> float:
    """Computes the Penetration Index for a single atom pair A-B (eqn 1 in
    Echeverria & Alvarez, Chem. Sci. 2023, 14, 11647-11688):

        PI = 100 * (v_A + v_B - d_AB) / (v_A + v_B - r_A - r_B)

    where r is the covalent radius, v the van der Waals radius, and d_AB the
    interatomic distance. This is the same formula used internally by
    post_process_hbonds/detect_generic_interactions, exposed standalone for
    single-pair calculations (validation tests, the batch module).
    """
    return (100 * (vdw_radius_a + vdw_radius_b - distance)
            / (vdw_radius_a + vdw_radius_b - cov_radius_a - cov_radius_b))


def compute_pi_for_elements(element_1: str, element_2: str, distance: float) -> float:
    """Convenience wrapper around compute_pi() that looks up both atoms in
    the built-in RADII table by element symbol (case-insensitive)."""
    r1 = RADII[element_1.upper()]
    r2 = RADII[element_2.upper()]
    return compute_pi(distance, r1['r'], r1['v'], r2['r'], r2['v'])


def compute_individual_penetrations(
    distance: float, cov_radius_a: float, vdw_radius_a: float,
    cov_radius_b: float, vdw_radius_b: float,
) -> Dict[str, float]:
    """Computes the joint Penetration Index p_AB (identical to compute_pi())
    plus the two individual, atom-specific penetration indices introduced
    by Echeverria & Alvarez (Chem. Sci. 2024, 15, 12166-12168, eqn 1-5) to
    give extra information for atom pairs with very different ("misfit")
    van der Waals crust widths - e.g. a light atom in contact with a heavy
    metal, where p_AB alone can be harder to interpret.

    p_a is atom A's penetration scaled by atom B's own crust width (and
    vice versa for p_b) - each atom's individual index reflects how much of
    the OTHER atom's crust has been used up, not its own. When the two
    crust widths are similar, p_a ~ p_b ~ p_ab and this extra detail adds
    little; the wider the gap, the more p_a and p_b diverge from p_ab and
    from each other (their ratio equals the crust-width ratio, eqn 9).

    Returns p_ab, p_a, p_b, and crust_width_ratio (>= 1, a quick "how
    misfit is this pair" gauge; the source paper's own worked examples are
    ~1.01 for a near-identical-width Ta-O pair and ~1.92 for the more
    mismatched Zn-Te pair).
    """
    w_a = vdw_radius_a - cov_radius_a
    w_b = vdw_radius_b - cov_radius_b
    i_ab = vdw_radius_a + vdw_radius_b - distance

    p_ab = 100 * i_ab / (w_a + w_b)
    p_a = 100 * i_ab / (2 * w_b)
    p_b = 100 * i_ab / (2 * w_a)
    ratio = max(w_a, w_b) / min(w_a, w_b)

    return {'p_ab': p_ab, 'p_a': p_a, 'p_b': p_b, 'crust_width_ratio': ratio}


def precompute_radii_arrays(
    universe: mda.Universe,
    dynamic_carbon_radii: bool = False,
    topology_path: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precomputes per-atom covalent radius, van der Waals radius, and
    donor-hydrogen bond length arrays for O(1) vectorized lookups.

    If dynamic_carbon_radii is True, refines carbon's covalent radius with
    the per-atom hybridization (sp/sp2/sp3) detected from topology_path
    (see apply_dynamic_carbon_radii).
    """
    n_atoms = len(universe.atoms)
    cov_radii = np.zeros(n_atoms)
    vdw_radii = np.zeros(n_atoms)
    dh_lengths = np.ones(n_atoms)

    for atom in universe.atoms:
        elem = (atom.element.upper() or atom.name[0].upper())
        idx = atom.index
        cov_radii[idx] = RADII.get(elem, {}).get('r', np.nan)
        vdw_radii[idx] = RADII.get(elem, {}).get('v', np.nan)
        dh_lengths[idx] = D_H_BONDS.get(elem, 1.0)

    if dynamic_carbon_radii:
        apply_dynamic_carbon_radii(cov_radii, topology_path)

    return cov_radii, vdw_radii, dh_lengths


def apply_dynamic_carbon_radii(cov_radii: np.ndarray, topology_path: str) -> None:
    """Opt-in refinement: overrides the flat sp2 default for carbon covalent
    radii (in-place, on the array precompute_radii_arrays already built)
    with the per-atom hybridization (sp/sp2/sp3) detected from topology_path
    via penindex.io.hybridization.HybridizationAnalyzer.

    Only overrides the covalent radius (used for the 100%-penetration
    reference point); the van der Waals radius is hybridization-independent
    in this radii set. Silently keeps the flat default if the topology
    format isn't supported by HybridizationAnalyzer, or if the optional
    RDKit/ParmEd/cclib dependencies aren't installed.
    """
    try:
        from penindex.io.hybridization import HybridizationAnalyzer
    except ImportError as e:
        warnings.warn(f"WARNING: dynamic_carbon_radii requires penindex.io.hybridization "
                      f"and its dependencies (rdkit, parmed, cclib); skipping ({e}).")
        return

    try:
        hyb_map = HybridizationAnalyzer(topology_path).get_carbon_hybridizations()
    except Exception as e:
        warnings.warn(f"WARNING: could not determine carbon hybridization ({e}); "
                      f"keeping the flat default covalent radius for carbon.")
        return

    if not hyb_map:
        warnings.warn("WARNING: no carbon hybridization data returned; "
                      "keeping the flat default covalent radius for carbon.")
        return

    n_atoms = len(cov_radii)
    counts = {'SP3': 0, 'SP2': 0, 'SP': 0, 'UNRESOLVED': 0}
    for idx, hyb in hyb_map.items():
        if idx >= n_atoms:
            continue  # index misalignment guard: source file doesn't match this Universe
        r = CARBON_HYBRID_COV_RADII.get(hyb)
        if r is not None:
            cov_radii[idx] = r
            counts[hyb] += 1
        else:
            counts['UNRESOLVED'] += 1

    print(f"Dynamic carbon radii applied: sp3={counts['SP3']}, sp2={counts['SP2']}, "
          f"sp={counts['SP']}, unresolved={counts['UNRESOLVED']}")


def detect_hbonds(
    universe: mda.Universe,
    donors_sel: str,
    acceptors_sel: str,
    d_a_cutoff: float,
    d_h_a_angle_cutoff: float,
    start: int,
    stop: Optional[int],
    step: int,
    reactive_mode: bool = False,
) -> Tuple[np.ndarray, int]:
    """Detects D-H...A hydrogen-bond triads with MDAnalysis's native
    algorithm. Returns (raw_hbonds_array, n_frames_analyzed)."""
    hbonds = HydrogenBondAnalysis(
        universe=universe,
        donors_sel=donors_sel,
        hydrogens_sel="element H or name H*",
        acceptors_sel=acceptors_sel,
        d_a_cutoff=d_a_cutoff,
        d_h_a_angle_cutoff=d_h_a_angle_cutoff,
        update_selections=reactive_mode,
    )
    hbonds.run(start=start, stop=stop, step=step)
    return hbonds.results.hbonds, hbonds.n_frames


def post_process_hbonds(
    raw_data: np.ndarray,
    cov_radii: np.ndarray,
    vdw_radii: np.ndarray,
    dh_lengths: np.ndarray,
) -> pd.DataFrame:
    """Applies C-level NumPy vectorization to reconstruct the exact
    Hydrogen...Acceptor distance (law of cosines) and compute the PI for
    every detected hydrogen-bond triad."""
    frames = raw_data[:, 0].astype(int)
    donor_idx = raw_data[:, 1].astype(int)
    hydro_idx = raw_data[:, 2].astype(int)
    acc_idx = raw_data[:, 3].astype(int)
    dist_da = raw_data[:, 4]
    rad_angles = np.radians(raw_data[:, 5])

    # O(1) array lookups
    d_dh = dh_lengths[donor_idx]
    v_a = vdw_radii[acc_idx]
    r_a = cov_radii[acc_idx]
    v_h, r_h = RADII['H']['v'], RADII['H']['r']

    # Vectorized law of cosines
    dist_ha = d_dh * np.cos(rad_angles) + np.sqrt(dist_da**2 - (d_dh * np.sin(rad_angles))**2)

    # Vectorized PI equation
    pi_values = 100 * (v_h + v_a - dist_ha) / (v_h + v_a - r_h - r_a)

    df = pd.DataFrame({
        'Frame': frames, 'Donor_Idx': donor_idx, 'Hydrogen_Idx': hydro_idx,
        'Acceptor_Idx': acc_idx, 'Distance_DA': dist_da, 'Angle_DHA': raw_data[:, 5],
        'Distance_HA': dist_ha, 'PI(%)': pi_values
    })

    return df.dropna(subset=['PI(%)'])


def detect_generic_interactions(
    universe: mda.Universe,
    selection_1: str,
    selection_2: str,
    cutoff: float,
    start: int,
    stop: Optional[int],
    step: int,
    cov_radii: np.ndarray,
    vdw_radii: np.ndarray,
    include_individual_penetrations: bool = False,
) -> Tuple[Optional[pd.DataFrame], int]:
    """Detects the Penetration Index between any two atom selections.

    Unlike the H-bond triad detector, this scans raw atom-atom distances
    directly (no D-H-A geometry), so the same PI equation captures the
    whole continuum: weak van der Waals contacts (PI ~ 0%), hydrogen bonds
    and other secondary interactions (PI ~ 20-70%), and covalent or
    metal-metal bonds (PI >= ~90%).

    If include_individual_penetrations is True, also adds the P_Atom1(%),
    P_Atom2(%) and Crust_Width_Ratio columns (see
    compute_individual_penetrations) and warns once if any detected pair
    has a notably mismatched ("misfit") crust-width ratio, for which the
    plain PI(%) column alone can be harder to interpret.

    Returns (DataFrame_or_None, n_frames_analyzed).
    """
    sel1 = universe.select_atoms(selection_1)
    sel2 = universe.select_atoms(selection_2)

    if len(sel1) == 0 or len(sel2) == 0:
        warnings.warn("WARNING: One or both generic-interaction selections are empty.")
        return None, 0

    idx1_global = sel1.indices
    idx2_global = sel2.indices

    frame_chunks, a_chunks, b_chunks, dist_chunks = [], [], [], []
    n_frames = 0

    for ts in universe.trajectory[start:stop:step]:
        n_frames += 1
        # Sparse distance matrix does the O(N log N) neighbor search AND
        # the distance calculation in a single vectorized C call.
        tree1 = cKDTree(sel1.positions)
        tree2 = cKDTree(sel2.positions)
        sdm = tree1.sparse_distance_matrix(tree2, max_distance=cutoff)
        if len(sdm) == 0:
            continue

        pairs = np.array(list(sdm.keys()), dtype=np.int64)
        dists = np.fromiter(sdm.values(), dtype=float, count=len(sdm))

        gi = idx1_global[pairs[:, 0]]
        gj = idx2_global[pairs[:, 1]]

        mask = gi != gj
        gi, gj, dists = gi[mask], gj[mask], dists[mask]

        frame_chunks.append(np.full(len(gi), ts.frame))
        a_chunks.append(np.minimum(gi, gj))
        b_chunks.append(np.maximum(gi, gj))
        dist_chunks.append(dists)

    if not a_chunks:
        warnings.warn("WARNING: No generic interactions found within the given cutoff.")
        return None, n_frames

    df = pd.DataFrame({
        'Frame': np.concatenate(frame_chunks),
        'Atom1_Idx': np.concatenate(a_chunks),
        'Atom2_Idx': np.concatenate(b_chunks),
        'Distance_AB': np.concatenate(dist_chunks),
    })
    # If the two selections overlap (e.g. both "all"), the same contact is
    # found from both directions; canonicalizing to (min, max) above and
    # dropping duplicates here collapses it to a single row per frame.
    df = df.drop_duplicates(subset=['Frame', 'Atom1_Idx', 'Atom2_Idx']).reset_index(drop=True)

    atom1_idx = df['Atom1_Idx'].to_numpy()
    atom2_idx = df['Atom2_Idx'].to_numpy()
    v_a, v_b = vdw_radii[atom1_idx], vdw_radii[atom2_idx]
    r_a, r_b = cov_radii[atom1_idx], cov_radii[atom2_idx]

    df['PI(%)'] = 100 * (v_a + v_b - df['Distance_AB']) / (v_a + v_b - r_a - r_b)
    df = df.dropna(subset=['PI(%)'])

    if include_individual_penetrations:
        atom1_idx = df['Atom1_Idx'].to_numpy()
        atom2_idx = df['Atom2_Idx'].to_numpy()
        v_a, v_b = vdw_radii[atom1_idx], vdw_radii[atom2_idx]
        r_a, r_b = cov_radii[atom1_idx], cov_radii[atom2_idx]

        # See compute_individual_penetrations: each atom's individual index
        # is scaled by the OTHER atom's crust width.
        w_a, w_b = v_a - r_a, v_b - r_b
        i_ab = v_a + v_b - df['Distance_AB']
        df['P_Atom1(%)'] = 100 * i_ab / (2 * w_b)
        df['P_Atom2(%)'] = 100 * i_ab / (2 * w_a)
        df['Crust_Width_Ratio'] = np.maximum(w_a, w_b) / np.minimum(w_a, w_b)

        n_misfit = int((df['Crust_Width_Ratio'] > 3.0).sum())
        if n_misfit > 0:
            warnings.warn(f"WARNING: {n_misfit} pair(s) have a van der Waals crust-width "
                          f"ratio > 3 (a notably 'misfit' pair per Echeverria & Alvarez, "
                          f"Chem. Sci. 2024); PI(%) alone may be less informative for "
                          f"these - see the P_Atom1(%)/P_Atom2(%) columns.")

    return df, n_frames
