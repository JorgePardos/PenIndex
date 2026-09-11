"""PyMOL script generation for publication-quality interaction figures.

Produces a .pml script that, when run in PyMOL (interactively, or headless
via `pymol -cq script.pml`), renders a clean, print-sized figure following
the current community-standard recipe for PyMOL publication figures: white
background, a stripped cartoon+sticks representation (rather than the
default all-lines view), dashed contacts colored by Penetration Index, and
optionally a ray-traced PNG export and an editable .pse session.

This project only *writes* the script - it does not invoke PyMOL itself, so
no PyMOL installation is required to use PenIndex; you run the generated
.pml file in your own PyMOL. generate_pi_color_legend produces a small,
separate matplotlib legend mapping color -> PI range (PyMOL has no good
native way to draw one); compose it with the rendered structure in a figure
layout tool.
"""

from typing import List, Optional, Sequence, Tuple

import pandas as pd

# (PyMOL color name, upper-bound-exclusive, human-readable label), low to
# high PI; the last entry's upper bound is None (catches everything above
# the previous threshold).
ColorThresholds = List[Tuple[str, Optional[float], str]]

HBOND_COLOR_THRESHOLDS: ColorThresholds = [
    ("skyblue", 30.0, "PI < 30% (weak / van der Waals)"),
    ("orange", 50.0, "30-50% (moderate)"),
    ("firebrick", None, "PI >= 50% (strong)"),
]

GENERIC_COLOR_THRESHOLDS: ColorThresholds = [
    ("skyblue", 20.0, "PI < 20% (van der Waals)"),
    ("orange", 70.0, "20-70% (non-covalent / secondary)"),
    ("firebrick", None, "PI >= 70% (covalent-like)"),
]

# PyMOL color name -> hex, so the matplotlib legend matches the render.
_PYMOL_COLOR_HEX = {"skyblue": "#87CEEB", "orange": "#FFA500", "firebrick": "#B22222"}


def _color_for_pi(pi: float, thresholds: ColorThresholds) -> str:
    for color, upper, _ in thresholds:
        if upper is None or pi < upper:
            return color
    return thresholds[-1][0]


def write_hbond_pymol_script(
    agg_df: pd.DataFrame,
    topology_path: str,
    out_file: str,
    total_frames: int,
    png_out: Optional[str] = None,
    session_out: Optional[str] = None,
    width_inches: float = 8.0,
    dpi: int = 300,
    label_mode: str = 'distance',
) -> None:
    """Publication-quality PyMOL script for hydrogen-bond contacts.

    label_mode controls what text appears on each contact: 'distance' (the
    default - PyMOL's own auto-generated H...A distance in Angstrom), 'pi'
    (the Penetration Index, as a percentage), or 'both'.
    """
    _write_script(
        agg_df, topology_path, out_file, total_frames,
        id_columns=('Hydrogen_Idx', 'Acceptor_Idx'),
        distance_column='Avg_Dist_HA',
        name_prefix='hb',
        color_thresholds=HBOND_COLOR_THRESHOLDS,
        static_pi_threshold=30.0,
        dynamic_occupancy_threshold=30.0,
        png_out=png_out, session_out=session_out,
        width_inches=width_inches, dpi=dpi,
        label_mode=label_mode,
    )
    print(f"PyMOL script generated: '{out_file}'")


def write_generic_pymol_script(
    agg_df: pd.DataFrame,
    topology_path: str,
    out_file: str,
    total_frames: int,
    png_out: Optional[str] = None,
    session_out: Optional[str] = None,
    width_inches: float = 8.0,
    dpi: int = 300,
    label_mode: str = 'distance',
) -> None:
    """Publication-quality PyMOL script for generic atom-pair contacts,
    coloring by Penetration Index across the whole covalent/non-covalent
    continuum (not just H-bonds).

    label_mode controls what text appears on each contact: 'distance' (the
    default - PyMOL's own auto-generated A-B distance in Angstrom), 'pi'
    (the Penetration Index, as a percentage), or 'both'.
    """
    _write_script(
        agg_df, topology_path, out_file, total_frames,
        id_columns=('Atom1_Idx', 'Atom2_Idx'),
        distance_column='Avg_Distance_AB',
        name_prefix='gi',
        color_thresholds=GENERIC_COLOR_THRESHOLDS,
        static_pi_threshold=20.0,
        dynamic_occupancy_threshold=30.0,
        png_out=png_out, session_out=session_out,
        width_inches=width_inches, dpi=dpi,
        label_mode=label_mode,
    )
    print(f"PyMOL script generated: '{out_file}'")


def _custom_label_commands(name: str, id1: int, id2: int, distance: float, pi: float, label_mode: str) -> List[str]:
    """PyMOL doesn't let a `distance` object's auto-label be overridden with
    arbitrary text directly, so this hides that auto-label and places a
    separate labeled pseudoatom at the contact's midpoint instead. The
    midpoint is computed by PyMOL itself at render time (via an inline
    `python`/`python end` block) rather than passed in from here, since this
    module never has the atoms' 3D coordinates - only the aggregated
    DataFrame."""
    if label_mode == 'pi':
        text = f"{pi:.0f}%"
    else:  # 'both'
        text = f"{distance:.2f} / {pi:.0f}%"

    label_name = f"lbl_{name}"
    return [
        f"hide labels, {name}",
        "python",
        f"_c1 = cmd.get_atom_coords('id {id1}')".format(id1=id1),
        f"_c2 = cmd.get_atom_coords('id {id2}')".format(id2=id2),
        "_mid = [(_a + _b) / 2.0 for _a, _b in zip(_c1, _c2)]",
        f"cmd.pseudoatom('{label_name}', pos=_mid, label='{text}')",
        "python end",
        f"hide nonbonded, {label_name}",
        f"set label_color, black, {label_name}",
        f"set label_size, 14, {label_name}",
    ]


def _write_script(
    agg_df: pd.DataFrame,
    topology_path: str,
    out_file: str,
    total_frames: int,
    id_columns: Tuple[str, str],
    distance_column: str,
    name_prefix: str,
    color_thresholds: ColorThresholds,
    static_pi_threshold: float,
    dynamic_occupancy_threshold: float,
    png_out: Optional[str],
    session_out: Optional[str],
    width_inches: float,
    dpi: int,
    label_mode: str = 'distance',
) -> None:
    if label_mode not in ('distance', 'pi', 'both'):
        raise ValueError(f"Unknown label_mode '{label_mode}'; use 'distance', 'pi', or 'both'")

    col_a, col_b = id_columns

    lines = [
        "# Auto-generated PyMOL publication script (PenIndex)",
        f"load {topology_path}, system",
        "bg_color white",
        "set ray_opaque_background, 1",
        "hide everything",
        "show cartoon, polymer",
        "color grey80, polymer",
        "set cartoon_transparency, 0.2",
        "set dash_radius, 0.08",
        "set dash_gap, 0.2",
        "",
    ]

    involved_ids = set()
    dash_lines = []
    for _, row in agg_df.iterrows():
        # CONTEXT-AWARE FILTERING: static single-structure runs filter by PI
        # itself (no time dimension to judge stability); trajectory runs
        # filter by occupancy (removes transient, thermally-driven contacts).
        if total_frames == 1:
            if row['Avg_PI'] < static_pi_threshold:
                continue
        else:
            if row['Occupancy'] < dynamic_occupancy_threshold:
                continue

        # PyMOL is 1-indexed, MDAnalysis is 0-indexed
        id1, id2 = int(row[col_a]) + 1, int(row[col_b]) + 1
        involved_ids.update((id1, id2))

        name = f"{name_prefix}_{id1}_{id2}"
        color = _color_for_pi(row['Avg_PI'], color_thresholds)
        dash_lines.append(f"distance {name}, id {id1}, id {id2}")
        dash_lines.append(f"color {color}, {name}")

        if label_mode != 'distance':
            dash_lines.extend(_custom_label_commands(
                name, id1, id2, row[distance_column], row['Avg_PI'], label_mode))

    if involved_ids:
        # Show the residues actually involved in a highlighted contact as
        # sticks - under cartoon-only representation their sidechain atoms
        # (the dash endpoints) would otherwise not be drawn at all.
        id_list = "+".join(str(i) for i in sorted(involved_ids))
        lines.append(f"show sticks, byres (id {id_list})")
        lines.append(f"util.cnc byres (id {id_list})")
        lines.append("")

    lines.extend(dash_lines)
    lines.append("")
    lines.append("orient")
    lines.append("zoom visible, 3")

    if png_out:
        width_px = int(round(width_inches * dpi))
        lines.append(f"ray {width_px}")
        lines.append(f"png {png_out}, dpi={dpi}")

    if session_out:
        lines.append(f"save {session_out}")

    with open(out_file, 'w') as f:
        f.write("\n".join(lines) + "\n")


def generate_pi_color_legend(
    color_thresholds: ColorThresholds,
    out_file: str,
    title: str = "Penetration Index",
) -> None:
    """Small standalone legend matching the color coding used in the .pml
    script. PyMOL has no good native way to draw a figure legend; compose
    this image with the rendered structure in a figure layout tool."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    handles = [
        Patch(facecolor=_PYMOL_COLOR_HEX.get(color, color), edgecolor='black', label=label)
        for color, _, label in color_thresholds
    ]

    fig, ax = plt.subplots(figsize=(3.5, 0.35 * len(handles) + 0.3))
    ax.axis('off')
    ax.legend(handles=handles, loc='center', frameon=False, title=title,
              fontsize=9, title_fontsize=10)
    fig.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
