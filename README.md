# PenIndex

Analyzer of covalent and non-covalent interactions (hydrogen bonds, van der Waals contacts, covalent bonds, metallophilic contacts, etc.) for computational-chemistry and biochemistry structures and trajectories, based on the **Penetration Index (PI)**.

> Echeverría, J., & Alvarez, S. (2023). *The borderless world of chemical bonding across the van der Waals crust and the valence region*. Chemical Science, 14, 11647–11688.
> Echeverría, J., & Alvarez, S. (2024). *A further focus on penetration indices of misfit van der Waals crusts*. Chemical Science, 15, 12166–12168.
> Pardos, J., Gonzalo, J., Merino, P., & Echeverría, J. (2026). *Beyond bond distances: a purely geometrical descriptor correlating with covalency-related bonding terms*. Dalton Transactions, 55, 3967–3974.

---

## 1. What is the Penetration Index?

The PI measures, for a pair of atoms A–B, how much their "van der Waals crusts" (the region between the covalent radius and the van der Waals radius) interpenetrate:

```
PI(%) = 100 * (v_A + v_B − d_AB) / (v_A + v_B − r_A − r_B)
```

where `v` is the van der Waals radius, `r` the covalent radius, and `d_AB` the actual interatomic distance. It places on the same scale:

| PI range | Interpretation |
|---|---|
| `< 0%` | No real contact / repulsive |
| `~ 0%` | Pure van der Waals contact |
| `~ 20–70%` | Non-covalent/secondary interactions: hydrogen, halogen, metallophilic, agostic bonds |
| `~ 100%` | Typical single covalent bond |
| `> 100%` | Multiple or very short covalent bond |

Beyond describing bonding regimes, PI **correlates strongly with independently computed covalency descriptors** (ALMO-EDA charge transfer, IQA exchange-correlation energy, DDEC6 bond orders, 3c–4e bonding percentages) — consistently better than raw interatomic distance (Pardos et al., Dalton Trans. 2026). This is the primary evidence behind PenIndex's design: it isn't just a structural-scanning tool, it's built to support that correlation workflow directly (see section 7).

---

## 2. Installation

Base environment (required):

```bash
conda env create -f enviroment.yml
conda activate hbond_pi_env
```

This installs Python, MDAnalysis, pandas, numpy, and pyyaml. Also install `pytest`, `scipy`, and `matplotlib` if not already present (scipy ships as an MDAnalysis dependency; matplotlib is needed for the plotting/PyMOL-legend features):

```bash
pip install pytest matplotlib
```

**Optional** dependencies, only needed for certain features:

| Feature | Requires | Install |
|---|---|---|
| `dynamic_carbon_radii` on `.pdb`/`.xyz` | `rdkit` | `conda install -c conda-forge rdkit` |
| `dynamic_carbon_radii` on `.prmtop`/`.top` | `parmed` | `conda install -c conda-forge parmed` (needs a C compiler; on Windows requires Visual C++ Build Tools) |
| Loading a Gaussian output (`.log`/`.out`), or `dynamic_carbon_radii` on it | `cclib` | `conda install -c conda-forge cclib` |

If an optional dependency is missing, the program **doesn't crash**: it warns via the console and continues with the default behavior.

---

## 3. Package layout

```
penindex/
    constants.py            RADII table, D-H bond lengths, carbon hybridization radii
    exceptions.py            PenIndexError and subclasses (raised by the library, caught by cli.py)
    core.py                  the PI engine: radii setup, hbond/generic detectors, compute_pi(),
                              compute_individual_penetrations(); pure computation, no I/O
    stats.py                 aggregation (occupancy/mean/std), linear_fit(), fit_distribution()
    batch.py                 batch single-pair mode (one bond per small system, e.g. a curated
                              diatomics/H-bond-complex dataset)
    provenance.py            run manifest (config snapshot, package versions, timestamp)
    cli.py                   entry point - the only module that prints, writes files, or exits
    io/
        loaders.py            builds an MDAnalysis Universe from any supported input
        hybridization.py       carbon sp/sp2/sp3 detection (RDKit/cclib/ParmEd backends)
        descriptors.py         loads/joins externally computed covalency descriptors
    viz/
        plots.py               matplotlib publication figures
        pymol_publication.py   PyMOL script generation (+ a matching legend image)
tests/                        pytest suite, including regression tests against literature values
config.yaml                  configuration for a structural-scan run (see section 5)
hbond_analyzer.py            thin backward-compatible shim: `python hbond_analyzer.py -c config.yaml`
```

`hbond_analyzer.py` and `hybridization_analyzer.py` at the project root are unchanged entry points that just import from the `penindex` package — old commands and scripts still work. `distances.py` and `penIndex_jorge.py` are earlier prototypes, kept locally for reference but not tracked in this repository. Likewise, `earp.pdb` (the example structure used during development) and its example outputs (`hbonds_standard_stats.csv`, `visualize_hbonds.pml`) are local-only — they're real, unpublished research data, not generic example data, so run the quickstart below against your own structure file.

---

## 4. Quick start (structural scan)

1. Edit `config.yaml` (see the full reference below).
2. Run:

```bash
python hbond_analyzer.py -c config.yaml
```

3. Check the generated CSV files and, to visually inspect the contacts, open the `.pml` in PyMOL:

```bash
pymol visualize_hbonds.pml
```

---

## 5. `config.yaml` reference

### `system`

| Key | Required | Description |
|---|---|---|
| `topology` | Yes | PDB/Amber/Gromacs/etc. (via MDAnalysis), or a Gaussian output (`.log`/`.out`, via `cclib`). |
| `trajectory` | No | MD trajectory (MDAnalysis-native topologies only). |
| `start_frame` / `stop_frame` / `step_frame` | No | Trajectory slicing (0-indexed; `-1` = through the end). For a QM output, each geometry (opt/IRC/scan step) counts as a frame. |
| `dynamic_carbon_radii` | No (`false`) | Refine carbon's covalent radius by detected hybridization (sp/sp2/sp3) instead of the flat sp2 default. |

### `hbond_params`

`max_distance` / `min_angle`: first-pass D-H···A detection cutoffs (hydrogen-bond detector only).

### `analysis`

| Key | Description |
|---|---|
| `donors` / `acceptors` | MDAnalysis selections for the hydrogen-bond detector. |
| `reactive_mode` | `false` = aggregated averages (with `Std_PI`) over the whole run. `true` = discrete formation/breaking events. Applies to both detectors. |
| `hbond_enabled` | Turn off the hydrogen-bond detector if you only want the generic scan. |
| `pymol_render` | Optional: `png`/`session`/`width_inches`/`dpi` — have the generated `.pml` also ray-trace a print-sized PNG and/or save a `.pse` session when it's run in PyMOL (this project only writes the script; you run it). A matching legend PNG is generated alongside. |
| `generic_interactions.enabled` | Turn on the generic PI detector (any atom pair, not just H-bonds). |
| `generic_interactions.selection_1` / `selection_2` | MDAnalysis selections between which pairs are searched. |
| `generic_interactions.max_distance` | First-pass distance cutoff (Å) for the KD-tree search. |
| `generic_interactions.include_individual_penetrations` | Optional: adds the individual `P_Atom1(%)`/`P_Atom2(%)`/`Crust_Width_Ratio` columns (Echeverría & Alvarez, 2024) — useful for atom pairs with very mismatched crust widths (e.g. a light atom against a heavy metal). Warns if any detected pair is notably mismatched. |

### `output` (optional)

`directory` / `prefix` — where result files land, so multiple runs/systems don't overwrite each other. Omit to keep writing to the current directory with no prefix.

### `provenance` (optional)

`enabled: true` writes a `run_manifest.json` alongside the other outputs: exact resolved config, timestamp, and key package versions — a reproducibility record for a manuscript's methods section.

---

## 6. The two structural detectors

- **Hydrogen bonds** (`hbond_enabled`): D-H···A triads via MDAnalysis's `HydrogenBondAnalysis`, exact H···A distance via the law of cosines, PI vectorized in NumPy.
- **Generic scan** (`generic_interactions`): PI between **any** pair of atoms, no D-H-A geometry required — the full covalent/non-covalent continuum, using `scipy.spatial.cKDTree.sparse_distance_matrix` for the neighbor search + distance calculation in one vectorized call.

Both write a standard-mode CSV (occupancy, mean, `Std_PI`) or a reactive-mode event table (see caveat on event definition in section 10), plus a `.pml` script.

---

## 7. Batch mode + external descriptor correlation

This is the workflow behind Pardos, Gonzalo, Merino & Echeverría (Dalton Trans. 2026): a curated set of small systems, each reduced to a single PI value for the bond of interest, correlated against an externally computed covalency descriptor (from AIMAll/IQA, Q-Chem ALMO-EDA, DDEC6, etc.).

**1. A batch systems CSV** (one bond per structure file):

```csv
label,structure,atom_1,atom_2
F2,F2.xyz,index 0,index 1
Cl2,Cl2.xyz,index 0,index 1
Br2,Br2.xyz,index 0,index 1
```

**2. A descriptor CSV** (same `label`s, any numeric column name):

```csv
label,delta_E_SC
F2,-42.0
Cl2,-38.0
Br2,-33.0
```

**3. Run it:**

```python
from penindex.batch import load_batch_systems_csv, run_batch
from penindex.io.descriptors import load_descriptor_table, join_batch_with_descriptors
from penindex.viz.plots import plot_descriptor_correlation

systems = load_batch_systems_csv("systems.csv")
batch_df = run_batch(systems)                    # computes Distance_AB and PI(%) per system

descriptors = load_descriptor_table("descriptors.csv")
merged = join_batch_with_descriptors(batch_df, descriptors)

fits = plot_descriptor_correlation(
    merged['Distance_AB'], merged['PI(%)'], merged['delta_E_SC'],
    descriptor_label=r"$\Delta E_{SC}$ (kcal/mol)",
    out_file="correlation.png",   # or .pdf/.svg for a vector figure
)
print(fits['distance']['r2'], fits['pi']['r2'])   # the core "does PI beat distance" comparison
```

This produces the same two-panel figure (distance-vs-descriptor next to PI-vs-descriptor, each with a fit line and R² annotation) used throughout that paper. `plot_r2_summary()` reproduces its Fig. 8-style R²-by-dataset bar chart across several such runs.

Systems that fail to load, or whose `atom_1`/`atom_2` selectors don't resolve to exactly one atom each, are skipped with a warning rather than aborting the whole batch.

---

## 8. Publication figures

`penindex.viz.plots` (matplotlib, `Agg` backend, colorblind-safe palette):

- `plot_descriptor_correlation` — the two-panel correlation figure above.
- `plot_r2_summary` — R²(distance) vs. R²(PI) bar chart across datasets.
- `plot_pi_histogram` — PI distribution, with an optional Gaussian/Lorentzian fit overlay (`distribution_fit='gaussian'|'lorentzian'`) to extract a peak position and width, matching how Echeverría & Alvarez characterize a class of interactions. **Fit a single distribution only to a chemically coherent population** (e.g. one bond type) — an unfiltered "all pairs" scan is naturally bimodal (van der Waals hump + covalent-bond peak) and a single Gaussian/Lorentzian won't describe it well.
- `plot_pi_timeseries` — PI vs. frame/step for one tracked interaction, useful alongside a `reactive_mode` event table.

`penindex.viz.pymol_publication` writes a `.pml` script following the current community-standard recipe for PyMOL publication figures: white background, a stripped cartoon+sticks representation (not the default all-lines view), and — if `pymol_render` is enabled in the config — a ray-traced, print-sized PNG export and/or an editable `.pse` session. A separate legend PNG (`generate_pi_color_legend`) maps color → PI range, since PyMOL has no good native way to draw one; compose it with the rendered structure in a figure layout tool. This project only writes the script — running it (`pymol -cq script.pml` or interactively) requires your own PyMOL install.

---

## 9. QM energies and individual penetrations

- Loading a Gaussian/QM output also attaches `universe.trajectory.qm_energies` (a NumPy array aligned with the frame index, from cclib's `scfenergies`, in eV) when cclib parsed one matching the number of geometries. Combine with a per-frame PI series (`df['QM_Energy'] = universe.trajectory.qm_energies[df['Frame']]`) for a PI-vs-energy fit along an optimization/IRC path, mirroring the modified Lennard-Jones fitting used for weak interactions in Echeverría & Alvarez (2023, eqn 3–4).
- `core.compute_individual_penetrations()` (and the `include_individual_penetrations` config flag for the generic scan) exposes the individual, atom-specific penetration indices `p_A`/`p_B` introduced in Echeverría & Alvarez (2024) for pairs with mismatched van der Waals crust widths, alongside the standard, symmetric `p_AB`.

---

## 10. Known limitations

- The generic scan uses plain Euclidean distances (no periodic boundary conditions).
- `RADII` uses fixed values per element (except carbon, with `dynamic_carbon_radii`); Mn/Fe/Co default to the high-spin radius.
- ParmEd (`.prmtop`/`.top`) couldn't be tested in this environment (needs a C++ compiler not present here); the code path fails gracefully (warning + fallback) rather than crashing if something doesn't fit, but hasn't been exercised against a real file.
- The "event" in `reactive_mode` is defined by continuity of detection within the first-pass cutoff, not by a PI-value threshold crossing.

**Validated against real data, not just synthetic tests:** the Gaussian ingestion path and the PyMOL rendering have both been run against two real DNA base-triple QM outputs (T·A·T Hoogsteen adducts, 57 atoms, 17–21 optimization steps each) with a real PyMOL install. This surfaced and fixed a real bug: a generated `.pml` script always did `load <original topology_path>`, which silently loaded **zero atoms** for a QM source (PyMOL's own format plugins don't understand Gaussian logs — confirmed via its GAMESS plugin rejecting the file). The fix (`io.loaders.write_frame_as_pdb` + `cli._pymol_structure_path`): when the analyzed structure came from a QM output, a companion PDB of the exact analyzed frame is written and loaded instead, guaranteeing the render matches what was analyzed. Both files rendered correctly afterward (labeled H-bond distances, correct element coloring, ray-traced PNG + `.pse` session).

---

## 11. Testing

```bash
pytest tests/
```

Includes regression tests against exact literature values (He₂ dimer → PI = −4.3%; H₂⁺ → PI = 75%; Zn–Te/Ta–O crust-width ratios) from Echeverría & Alvarez (2023, 2024), so a change to the radii table or formula that breaks agreement with the published numbers fails the suite.

---

## 12. Roadmap

- Packaging as an installable Python package (deferred by choice).
- ML ideas from a local brainstorming doc (not tracked in this repo): clustering reactive events, ML-based dynamic van der Waals radii — phase 2, doesn't block current use.
