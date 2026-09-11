"""Tests for detect_generic_interactions(include_individual_penetrations=True):
the vectorized P_Atom1(%)/P_Atom2(%)/Crust_Width_Ratio columns must agree
with the scalar core.compute_individual_penetrations() reference
implementation, and the misfit-crust warning must fire only when warranted.
"""

import numpy as np
import pytest
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader

from penindex import core
from penindex.constants import RADII


def _two_atom_universe(distance: float) -> mda.Universe:
    u = mda.Universe.empty(2, n_residues=1, atom_resindex=np.zeros(2, dtype=int), trajectory=True)
    u.add_TopologyAttr('name', ['A1', 'A2'])
    u.add_TopologyAttr('resname', ['X'])
    u.add_TopologyAttr('resid', [1])
    u.add_TopologyAttr('segid', ['S'])
    coords = np.array([[[0.0, 0.0, 0.0], [0.0, 0.0, distance]]], dtype=np.float32)
    u.load_new(coords, format=MemoryReader)
    return u


def test_vectorized_matches_scalar_reference_for_zn_te():
    distance = 2.6
    u = _two_atom_universe(distance)
    r_zn, r_te = RADII['ZN'], RADII['TE']
    cov_radii = np.array([r_zn['r'], r_te['r']])
    vdw_radii = np.array([r_zn['v'], r_te['v']])

    df, n_frames = core.detect_generic_interactions(
        u, selection_1="all", selection_2="all", cutoff=4.0,
        start=0, stop=None, step=1,
        cov_radii=cov_radii, vdw_radii=vdw_radii,
        include_individual_penetrations=True,
    )

    expected = core.compute_individual_penetrations(
        distance, r_zn['r'], r_zn['v'], r_te['r'], r_te['v'])

    assert 'P_Atom1(%)' in df.columns
    assert 'Crust_Width_Ratio' in df.columns
    row = df.iloc[0]
    assert row['P_Atom1(%)'] == pytest.approx(expected['p_a'])
    assert row['P_Atom2(%)'] == pytest.approx(expected['p_b'])
    assert row['Crust_Width_Ratio'] == pytest.approx(expected['crust_width_ratio'])
    assert row['Crust_Width_Ratio'] == pytest.approx(1.92, abs=0.02)


def test_no_columns_added_when_flag_is_off():
    u = _two_atom_universe(2.6)
    r_zn, r_te = RADII['ZN'], RADII['TE']
    cov_radii = np.array([r_zn['r'], r_te['r']])
    vdw_radii = np.array([r_zn['v'], r_te['v']])

    df, _ = core.detect_generic_interactions(
        u, selection_1="all", selection_2="all", cutoff=4.0,
        start=0, stop=None, step=1,
        cov_radii=cov_radii, vdw_radii=vdw_radii,
    )

    assert 'P_Atom1(%)' not in df.columns


def test_warns_on_notably_mismatched_crust_widths():
    # Synthetic, deliberately extreme radii (not tied to a real element) to
    # exercise the crust_width_ratio > 3 warning path directly.
    u = _two_atom_universe(1.0)
    cov_radii = np.array([1.0, 0.5])
    vdw_radii = np.array([1.1, 2.0])  # widths: 0.1 and 1.5 -> ratio = 15

    with pytest.warns(UserWarning, match="misfit"):
        df, _ = core.detect_generic_interactions(
            u, selection_1="all", selection_2="all", cutoff=4.0,
            start=0, stop=None, step=1,
            cov_radii=cov_radii, vdw_radii=vdw_radii,
            include_individual_penetrations=True,
        )

    assert df['Crust_Width_Ratio'].iloc[0] > 3.0
