"""Tests for the QM ingestion path. cclib.io.ccread is mocked with a
synthetic ccData-like object so these run without a real Gaussian log file -
the same technique already used manually to validate this path earlier in
the project (cclib itself is a mature, independently-tested library; what's
under test here is penindex's own glue code around it).
"""

from unittest.mock import patch

import numpy as np
import pytest
import cclib

from penindex.exceptions import StructureLoadError
from penindex.io.loaders import load_qm_universe, write_frame_as_pdb


class _FakeWaterCCData:
    atomnos = np.array([8, 1, 1])
    atomcoords = np.array([
        [[0.0, 0.0, 0.1173], [0.0, 0.7572, -0.4692], [0.0, -0.7572, -0.4692]],
        [[0.0, 0.0, 0.1150], [0.0, 0.7600, -0.4650], [0.0, -0.7600, -0.4650]],
    ])


def test_load_qm_universe_builds_expected_universe(tmp_path):
    fake_path = str(tmp_path / "fake.log")

    with patch.object(cclib.io, 'ccread', return_value=_FakeWaterCCData()):
        universe = load_qm_universe(fake_path)

    assert len(universe.atoms) == 3
    assert len(universe.trajectory) == 2
    assert list(universe.atoms.elements) == ['O', 'H', 'H']
    np.testing.assert_allclose(universe.atoms.positions, _FakeWaterCCData.atomcoords[0], atol=1e-4)


def test_load_qm_universe_raises_on_unparseable_file(tmp_path):
    fake_path = str(tmp_path / "empty.log")

    class _EmptyData:
        pass

    with patch.object(cclib.io, 'ccread', return_value=_EmptyData()):
        with pytest.raises(StructureLoadError):
            load_qm_universe(fake_path)


class _FakeWaterCCDataWithEnergies(_FakeWaterCCData):
    scfenergies = np.array([-2079.1, -2079.3])  # eV, one per geometry


def test_load_qm_universe_attaches_matching_energies(tmp_path):
    fake_path = str(tmp_path / "fake_with_energy.log")

    with patch.object(cclib.io, 'ccread', return_value=_FakeWaterCCDataWithEnergies()):
        universe = load_qm_universe(fake_path)

    np.testing.assert_allclose(universe.trajectory.qm_energies, [-2079.1, -2079.3])


class _FakeWaterCCDataMismatchedEnergies(_FakeWaterCCData):
    scfenergies = np.array([-2079.1, -2079.2, -2079.3])  # 3 energies, 2 geometries


def test_load_qm_universe_drops_mismatched_energies_with_warning(tmp_path):
    fake_path = str(tmp_path / "fake_mismatched.log")

    with patch.object(cclib.io, 'ccread', return_value=_FakeWaterCCDataMismatchedEnergies()):
        with pytest.warns(UserWarning, match="SCF energies"):
            universe = load_qm_universe(fake_path)

    assert universe.trajectory.qm_energies is None


def test_load_qm_universe_energies_none_when_absent(tmp_path):
    fake_path = str(tmp_path / "fake_no_energy.log")

    with patch.object(cclib.io, 'ccread', return_value=_FakeWaterCCData()):
        universe = load_qm_universe(fake_path)

    assert universe.trajectory.qm_energies is None


def test_write_frame_as_pdb_writes_correct_geometry(tmp_path):
    fake_path = str(tmp_path / "fake.log")

    with patch.object(cclib.io, 'ccread', return_value=_FakeWaterCCData()):
        universe = load_qm_universe(fake_path)

    pdb_path = str(tmp_path / "frame1.pdb")
    write_frame_as_pdb(universe, frame_index=1, out_file=pdb_path)

    # A fresh Universe loaded from the written PDB must reproduce frame 1's
    # geometry (this is exactly what PyMOL will load in the generated .pml).
    import MDAnalysis as mda
    reloaded = mda.Universe(pdb_path)
    assert len(reloaded.atoms) == 3
    np.testing.assert_allclose(
        reloaded.atoms.positions, _FakeWaterCCData.atomcoords[1], atol=1e-2)
