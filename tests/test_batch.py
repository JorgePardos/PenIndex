import numpy as np
import pytest

from penindex import core
from penindex.batch import BatchSystem, run_batch


def _write_xyz(path, symbol_1, symbol_2, distance):
    path.write_text(f"2\ntest\n{symbol_1} 0.0 0.0 0.0\n{symbol_2} 0.0 0.0 {distance:.5f}\n")


def test_run_batch_matches_direct_pi_calculation(tmp_path):
    h2_path = tmp_path / "h2.xyz"
    hf_path = tmp_path / "hf.xyz"
    _write_xyz(h2_path, "H", "H", 0.741)
    _write_xyz(hf_path, "H", "F", 0.917)

    systems = [
        BatchSystem(label="H2", structure=str(h2_path), atom_1="index 0", atom_2="index 1"),
        BatchSystem(label="HF", structure=str(hf_path), atom_1="index 0", atom_2="index 1"),
    ]

    result = run_batch(systems)

    assert list(result['label']) == ["H2", "HF"]
    expected_h2_pi = core.compute_pi_for_elements("H", "H", 0.741)
    expected_hf_pi = core.compute_pi_for_elements("H", "F", 0.917)
    np.testing.assert_allclose(result.loc[result['label'] == 'H2', 'PI(%)'].iloc[0], expected_h2_pi)
    np.testing.assert_allclose(result.loc[result['label'] == 'HF', 'PI(%)'].iloc[0], expected_hf_pi)


def test_run_batch_skips_bad_selectors_with_warning(tmp_path):
    h2_path = tmp_path / "h2.xyz"
    _write_xyz(h2_path, "H", "H", 0.741)

    systems = [
        # "index 0 1" resolves to 2 atoms, not the required 1.
        BatchSystem(label="bad_selector", structure=str(h2_path), atom_1="index 0 1", atom_2="index 1"),
    ]

    with pytest.warns(UserWarning, match="bad_selector"):
        result = run_batch(systems)

    assert len(result) == 0
