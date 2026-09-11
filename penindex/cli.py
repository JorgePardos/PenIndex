"""Command-line entry point. This is the only module allowed to print
progress messages, write output files, or exit the process - everything
else (penindex.io, penindex.core, penindex.stats) is a plain library that
raises exceptions and returns data.
"""

import argparse
import os
import sys
import warnings
from typing import Any, Dict, List, Optional, Tuple

import yaml

from penindex import core, provenance, stats
from penindex.exceptions import ConfigError, PenIndexError
from penindex.io.loaders import QM_OUTPUT_EXTENSIONS, build_universe, write_frame_as_pdb
from penindex.viz import pymol_publication


def load_config(path: str) -> Dict[str, Any]:
    """Loads and parses the YAML configuration file."""
    try:
        with open(path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        raise ConfigError(f"Error reading YAML input '{path}': {e}") from e


def run(config_path: str) -> None:
    """Main execution pipeline. Dispatches the H-bond triad detector and/or
    the generic atom-pair detector according to the config switches."""
    config = load_config(config_path)
    universe = build_universe(config)

    sys_cfg = config['system']
    cov_radii, vdw_radii, dh_lengths = core.precompute_radii_arrays(
        universe,
        dynamic_carbon_radii=sys_cfg.get('dynamic_carbon_radii', False),
        topology_path=sys_cfg['topology'],
    )

    analysis_cfg = config['analysis']
    is_reactive = analysis_cfg.get('reactive_mode', False)

    # Every _output_path() call below records here, so a provenance
    # manifest (if enabled) can list exactly what this run produced.
    output_files: List[str] = []

    if analysis_cfg.get('hbond_enabled', True):
        _run_hbond_pipeline(config, universe, cov_radii, vdw_radii, dh_lengths, is_reactive, output_files)

    generic_cfg = analysis_cfg.get('generic_interactions', {})
    if generic_cfg.get('enabled', False):
        _run_generic_pipeline(config, universe, cov_radii, vdw_radii, generic_cfg, is_reactive, output_files)

    if config.get('provenance', {}).get('enabled', False) and output_files:
        manifest = provenance.build_manifest(config, config_path, output_files)
        manifest_path = _output_path(config, "run_manifest.json")
        provenance.write_manifest(manifest, manifest_path)
        print(f"Run manifest saved to '{manifest_path}'")


def _frame_slice(sys_cfg: Dict[str, Any]) -> Tuple[int, Optional[int], int]:
    start = sys_cfg.get('start_frame', 0)
    stop = sys_cfg.get('stop_frame', None)
    stop = None if stop == -1 else stop
    step = sys_cfg.get('step_frame', 1)
    return start, stop, step


def _output_path(config: Dict[str, Any], filename: str, track: Optional[List[str]] = None) -> str:
    """Resolves an output filename against the optional 'output' config
    block (directory + filename prefix), creating the directory if needed.
    Defaults (directory='.', prefix='') reproduce the previous hardcoded
    cwd-only behavior exactly. If `track` is given, the resolved path is
    appended to it (see run()'s output_files list)."""
    output_cfg = config.get('output') or {}
    directory = output_cfg.get('directory', '.')
    prefix = output_cfg.get('prefix', '')
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, f"{prefix}{filename}")
    if track is not None:
        track.append(path)
    return path


def _pymol_structure_path(config: Dict[str, Any], universe, start_frame: int, output_files: List[str]) -> str:
    """Returns the path a generated .pml script should `load`.

    For formats PyMOL understands natively (PDB, Amber, Gromacs, ...) this
    is just the original topology_path, unchanged. For a QM output
    (.log/.out, loaded via cclib), PyMOL cannot parse the file itself - its
    own format plugins are independent of cclib's and don't cover Gaussian
    logs, so `load <the .log path>` silently loads zero atoms. In that case
    a companion PDB of the exact analyzed frame is written instead (frame
    = start_frame, i.e. the first frame in the analyzed range - the single
    frame for a static/single-point analysis, or a defensible representative
    snapshot for a multi-frame QM reactive_mode run) and its path is
    returned so the render actually shows the same geometry that was analyzed.
    """
    topology_path = config['system']['topology']
    ext = os.path.splitext(topology_path)[1].lower()
    if ext not in QM_OUTPUT_EXTENSIONS:
        return topology_path

    pdb_path = _output_path(config, "_qm_structure.pdb", output_files)
    write_frame_as_pdb(universe, start_frame, pdb_path)
    return pdb_path


def _pymol_render_options(config: Dict[str, Any], base_stem: str, output_files: List[str]) -> Dict[str, Any]:
    """Resolves the optional analysis.pymol_render config block into kwargs
    for the write_*_pymol_script() functions, deriving PNG/session
    filenames from base_stem (e.g. 'visualize_hbonds'). Ray-tracing/session
    export are off by default (they require the .pml script to actually be
    run in a working PyMOL install, which this project doesn't assume), but
    label_mode defaults to 'pi' here - showing the Penetration Index on
    each contact, rather than PyMOL's native raw-distance label, is this
    program's normal behavior; pass label_mode='distance' explicitly if you
    want PyMOL's default instead."""
    render_cfg = config.get('analysis', {}).get('pymol_render') or {}
    kwargs: Dict[str, Any] = {
        'width_inches': render_cfg.get('width_inches', 8.0),
        'dpi': render_cfg.get('dpi', 300),
        'label_mode': render_cfg.get('label_mode', 'pi'),
    }
    if render_cfg.get('png', False):
        kwargs['png_out'] = _output_path(config, f"{base_stem}.png", output_files)
    if render_cfg.get('session', False):
        kwargs['session_out'] = _output_path(config, f"{base_stem}.pse", output_files)
    return kwargs


def _run_hbond_pipeline(config, universe, cov_radii, vdw_radii, dh_lengths, is_reactive: bool,
                         output_files: List[str]) -> None:
    """Detects D-H...A hydrogen bonds and their Penetration Index."""
    sys_cfg = config['system']
    start, stop, step = _frame_slice(sys_cfg)

    print("\nStep 1: Extracting Hydrogen-Bond Geometric Time-Series (MDAnalysis)...")
    raw_data, n_frames = core.detect_hbonds(
        universe,
        donors_sel=config['analysis']['donors'],
        acceptors_sel=config['analysis']['acceptors'],
        d_a_cutoff=config['hbond_params']['max_distance'],
        d_h_a_angle_cutoff=config['hbond_params']['min_angle'],
        start=start, stop=stop, step=step,
        reactive_mode=is_reactive,
    )

    if len(raw_data) == 0:
        warnings.warn("WARNING: No hydrogen bonds detected in the given selections.")
        return

    print("Step 2: Vectorized Penetration Index Calculation (H-bonds)...")
    df = core.post_process_hbonds(raw_data, cov_radii, vdw_radii, dh_lengths)

    if is_reactive:
        print("Step 3: Mode -> REACTIVE (Tracking H-Bond Breaking/Formation)...")
        event_df = stats.aggregate_hbond_reactive(df, step)
        out_filename = _output_path(config, "hbonds_reactive_events.csv", output_files)
        event_df.to_csv(out_filename, index=False)
        print(f"Reactive events tracking saved to '{out_filename}'")
        print(f"   Total chemical interaction events recorded: {len(event_df)}")
    else:
        print("Step 3: Mode -> STANDARD (Equilibrium Aggregation, H-bonds)...")
        agg_df = stats.aggregate_hbond_standard(df, n_frames)
        csv_filename = _output_path(config, "hbonds_standard_stats.csv", output_files)
        agg_df.to_csv(csv_filename, index=False)
        print(f"Statistics saved to '{csv_filename}'")

        pml_filename = _output_path(config, "visualize_hbonds.pml", output_files)
        structure_path = _pymol_structure_path(config, universe, start, output_files)
        render_kwargs = _pymol_render_options(config, "visualize_hbonds", output_files)
        pymol_publication.write_hbond_pymol_script(
            agg_df, structure_path, pml_filename, n_frames, **render_kwargs)
        if 'png_out' in render_kwargs or 'session_out' in render_kwargs:
            legend_filename = _output_path(config, "visualize_hbonds_legend.png", output_files)
            pymol_publication.generate_pi_color_legend(pymol_publication.HBOND_COLOR_THRESHOLDS, legend_filename)
            print(f"Legend saved to '{legend_filename}'")


def _run_generic_pipeline(config, universe, cov_radii, vdw_radii, generic_cfg: Dict[str, Any], is_reactive: bool,
                           output_files: List[str]) -> None:
    """Detects the Penetration Index between any pair of atoms across two
    selections, covering the full covalent/non-covalent continuum."""
    sys_cfg = config['system']
    start, stop, step = _frame_slice(sys_cfg)

    print(f"\nStep 4: Scanning generic atom-pair interactions "
          f"('{generic_cfg.get('selection_1', 'all')}' <-> '{generic_cfg.get('selection_2', 'all')}')...")
    df, n_frames = core.detect_generic_interactions(
        universe,
        selection_1=generic_cfg.get('selection_1', 'all'),
        selection_2=generic_cfg.get('selection_2', 'all'),
        cutoff=generic_cfg.get('max_distance', 4.0),
        start=start, stop=stop, step=step,
        cov_radii=cov_radii, vdw_radii=vdw_radii,
        include_individual_penetrations=generic_cfg.get('include_individual_penetrations', False),
    )
    if df is None or len(df) == 0:
        return

    if is_reactive:
        print("Step 5: Mode -> REACTIVE (Tracking Generic Interaction Breaking/Formation)...")
        event_df = stats.aggregate_generic_reactive(df, step)
        out_filename = _output_path(config, "generic_interactions_reactive_events.csv", output_files)
        event_df.to_csv(out_filename, index=False)
        print(f"Generic reactive events tracking saved to '{out_filename}'")
        print(f"   Total generic interaction events recorded: {len(event_df)}")
    else:
        print("Step 5: Mode -> STANDARD (Equilibrium Aggregation, Generic Interactions)...")
        agg_df = stats.aggregate_generic_standard(df, n_frames)
        csv_filename = _output_path(config, "generic_interactions_standard_stats.csv", output_files)
        agg_df.to_csv(csv_filename, index=False)
        print(f"Generic interaction statistics saved to '{csv_filename}'")

        pml_filename = _output_path(config, "visualize_generic_interactions.pml", output_files)
        structure_path = _pymol_structure_path(config, universe, start, output_files)
        render_kwargs = _pymol_render_options(config, "visualize_generic_interactions", output_files)
        pymol_publication.write_generic_pymol_script(
            agg_df, structure_path, pml_filename, n_frames, **render_kwargs)
        if 'png_out' in render_kwargs or 'session_out' in render_kwargs:
            legend_filename = _output_path(config, "visualize_generic_interactions_legend.png", output_files)
            pymol_publication.generate_pi_color_legend(pymol_publication.GENERIC_COLOR_THRESHOLDS, legend_filename)
            print(f"Legend saved to '{legend_filename}'")


def main() -> None:
    parser = argparse.ArgumentParser(description="Unified Penetration Index H-Bond Analyzer")
    parser.add_argument("-c", "--config", required=True, help="Path to the YAML configuration file")
    args = parser.parse_args()

    try:
        run(args.config)
    except PenIndexError as e:
        sys.exit(f"CRITICAL ERROR: {e}")


if __name__ == "__main__":
    main()
