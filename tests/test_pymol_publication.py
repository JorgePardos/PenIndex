import pandas as pd

from penindex.viz.pymol_publication import (
    GENERIC_COLOR_THRESHOLDS,
    HBOND_COLOR_THRESHOLDS,
    generate_pi_color_legend,
    write_generic_pymol_script,
    write_hbond_pymol_script,
)


def _hbond_agg_df():
    return pd.DataFrame({
        'Donor_Idx': [0, 1],
        'Hydrogen_Idx': [1, 2],
        'Acceptor_Idx': [5, 6],
        'Occupancy': [80.0, 10.0],
        'Avg_Dist_HA': [1.9, 2.5],
        'Avg_PI': [55.0, 10.0],
    })


def test_write_hbond_pymol_script_basic_structure(tmp_path):
    out_file = tmp_path / "hbonds.pml"
    write_hbond_pymol_script(_hbond_agg_df(), "system.pdb", str(out_file), total_frames=10)

    content = out_file.read_text()
    assert "bg_color white" in content
    assert "load system.pdb, system" in content
    assert "show cartoon, polymer" in content
    # Only the row with Occupancy >= 30 should produce a distance object.
    assert "distance hb_2_6" in content
    assert "hb_3_7" not in content
    assert "show sticks, byres (id 2+6)" in content


def test_write_hbond_pymol_script_with_png_and_session(tmp_path):
    out_file = tmp_path / "hbonds.pml"
    write_hbond_pymol_script(
        _hbond_agg_df(), "system.pdb", str(out_file), total_frames=1,
        png_out="fig.png", session_out="fig.pse", width_inches=4.0, dpi=150,
    )

    content = out_file.read_text()
    assert "ray 600" in content  # 4.0 in * 150 dpi
    assert "png fig.png, dpi=150" in content
    assert "save fig.pse" in content


def test_write_generic_pymol_script_basic_structure(tmp_path):
    agg_df = pd.DataFrame({
        'Atom1_Idx': [0, 3],
        'Atom2_Idx': [4, 7],
        'Occupancy': [100.0, 5.0],
        'Avg_Distance_AB': [1.5, 3.9],
        'Avg_PI': [95.0, 2.0],
    })
    out_file = tmp_path / "generic.pml"
    write_generic_pymol_script(agg_df, "system.pdb", str(out_file), total_frames=20)

    content = out_file.read_text()
    assert "distance gi_1_5" in content
    assert "color firebrick, gi_1_5" in content
    assert "gi_4_8" not in content  # filtered out by occupancy < 30


def test_generate_pi_color_legend_produces_file(tmp_path):
    out_file = tmp_path / "legend.png"
    generate_pi_color_legend(HBOND_COLOR_THRESHOLDS, str(out_file))
    assert out_file.exists()

    out_file_2 = tmp_path / "legend_generic.png"
    generate_pi_color_legend(GENERIC_COLOR_THRESHOLDS, str(out_file_2))
    assert out_file_2.exists()
