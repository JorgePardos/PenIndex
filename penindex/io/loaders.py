"""Builds an MDAnalysis Universe from any of the input formats PenIndex
supports: MDAnalysis-native topologies (PDB, Amber, Gromacs, ...) with an
optional trajectory, or a QM output (Gaussian .log/.out, and by extension
anything cclib understands) loaded entirely in memory.
"""

import os
import warnings
from typing import Any, Dict

import numpy as np
import MDAnalysis as mda

from penindex.exceptions import StructureLoadError

# Extensions handled via cclib instead of MDAnalysis's native readers.
QM_OUTPUT_EXTENSIONS = ('.log', '.out')


def build_universe(config: Dict[str, Any]) -> mda.Universe:
    """Instantiates the MDAnalysis Universe for a run and ensures topology
    integrity (charges, bonds). Raises StructureLoadError on failure instead
    of exiting the process, so this is safe to call as a library."""
    sys_cfg = config['system']
    topology_path = sys_cfg['topology']
    ext = os.path.splitext(topology_path)[1].lower()

    try:
        if ext in QM_OUTPUT_EXTENSIONS:
            if sys_cfg.get('trajectory'):
                warnings.warn("WARNING: 'trajectory' is ignored for QM outputs; every "
                              "geometry found in the file (opt/IRC/scan steps) is already "
                              "used as a trajectory frame.")
            u = load_qm_universe(topology_path)
        elif sys_cfg.get('trajectory'):
            u = mda.Universe(topology_path, sys_cfg['trajectory'])
            print(f"Loaded trajectory with {len(u.trajectory)} frames.")
        else:
            u = mda.Universe(topology_path)
            print(f"Loaded static topology: {topology_path}")

        # Topology patching for minimal files (e.g., standard PDBs or QM/MM outputs)
        if not hasattr(u.atoms, 'charges'):
            u.add_TopologyAttr('charges', np.zeros(len(u.atoms)))

        if not hasattr(u, 'bonds') or len(u.bonds) == 0:
            from MDAnalysis.topology.guessers import guess_bonds
            bonds = guess_bonds(u.atoms, u.atoms.positions)
            u.add_TopologyAttr('bonds', bonds)
            print(f"Auto-generated {len(u.bonds)} covalent bonds based on distances.")

        return u

    except StructureLoadError:
        raise
    except Exception as e:
        raise StructureLoadError(f"Could not load the molecular system: {e}") from e


def load_qm_universe(path: str) -> mda.Universe:
    """Builds an in-memory MDAnalysis Universe from any quantum-chemistry
    output cclib understands (Gaussian .log/.out, and by extension ORCA,
    Psi4, etc. via the same cclib.io.ccread entry point).

    Every geometry stored in the file - each optimization step, IRC point,
    or relaxed-scan point - becomes a trajectory frame. This means
    start_frame/stop_frame/step_frame in the config select QM steps
    exactly like they would select MD frames, and reactive_mode can track
    PI changes along an optimization or IRC path the same way it tracks
    bond breaking/forming in a QM/MM MD trajectory.
    """
    import cclib
    from cclib.parser.utils import PeriodicTable

    data = cclib.io.ccread(path)
    if data is None or not hasattr(data, 'atomcoords') or not hasattr(data, 'atomnos'):
        raise StructureLoadError(f"cclib could not parse atomic coordinates from '{path}'. "
                                  f"Check that it's a supported, complete QM output file.")

    periodic_table = PeriodicTable()
    elements = [periodic_table.element[int(z)] for z in data.atomnos]
    coords = np.asarray(data.atomcoords, dtype=np.float32)
    n_atoms = len(elements)

    u = mda.Universe.empty(
        n_atoms,
        n_residues=1,
        atom_resindex=np.zeros(n_atoms, dtype=int),
        trajectory=True,
    )
    u.add_TopologyAttr('name', elements)
    u.add_TopologyAttr('type', elements)
    u.add_TopologyAttr('elements', elements)
    u.add_TopologyAttr('resname', ['QM'])
    u.add_TopologyAttr('resid', [1])
    u.add_TopologyAttr('segid', ['SYSTEM'])

    from MDAnalysis.coordinates.memory import MemoryReader
    u.load_new(coords, format=MemoryReader)

    # Best-effort: cclib parses a per-geometry energy series from the same
    # file it just gave us coordinates from. Attaching it here (rather than
    # discarding it) makes a PI-vs-energy analysis of an optimization/IRC
    # path a one-liner: df['QM_Energy'] = u.trajectory.qm_energies[df['Frame']],
    # mirroring the energy-vs-PI fits used for weak interactions in
    # Echeverria & Alvarez, Chem. Sci. 2023 (eqn 3-4).
    u.trajectory.qm_energies = _extract_qm_energies(data, n_frames=coords.shape[0])

    print(f"Loaded QM output '{path}' via cclib: {n_atoms} atoms, "
          f"{len(u.trajectory)} geometries (opt/IRC/scan steps).")
    return u


def write_frame_as_pdb(universe: mda.Universe, frame_index: int, out_file: str) -> None:
    """Writes a single trajectory frame out as a standalone PDB file.

    Used to give PyMOL something it can actually load when the original
    topology_path is a format PyMOL doesn't understand (a Gaussian
    .log/.out): PyMOL's own file-format plugins are independent of, and
    generally narrower than, cclib's, so a generated .pml script that does
    `load <original .log path>` will silently load zero atoms for such
    inputs (confirmed: PyMOL's GAMESS plugin rejects a Gaussian log and the
    load simply produces an empty scene, no error).
    """
    universe.trajectory[frame_index]
    universe.atoms.write(out_file)


def _extract_qm_energies(data, n_frames: int):
    """Best-effort extraction of a per-geometry energy series (cclib's
    scfenergies, normalized to eV) from a cclib ccData object. Returns None
    (with a warning) if cclib didn't parse one, or if its length doesn't
    match the number of geometries loaded - a QM output can log a
    different number of energy entries than optimized geometries depending
    on the calculation type, and silently misaligning the two would corrupt
    any PI-vs-energy correlation."""
    energies = getattr(data, 'scfenergies', None)
    if energies is None:
        return None

    energies = np.asarray(energies, dtype=float)
    if len(energies) != n_frames:
        warnings.warn(f"WARNING: found {len(energies)} SCF energies but {n_frames} "
                      f"geometries in this QM output; cannot reliably align them "
                      f"frame-by-frame, so trajectory.qm_energies will be None.")
        return None

    return energies
