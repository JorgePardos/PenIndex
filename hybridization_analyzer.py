"""Backward-compatible entry point.

The implementation has moved to penindex/io/hybridization.py. This file is
kept so `from hybridization_analyzer import HybridizationAnalyzer` and
direct script usage keep working unchanged.
"""

from penindex.io.hybridization import HybridizationAnalyzer

if __name__ == "__main__":
    # To try it out, just create an instance with the path to your file.
    # Examples (uncomment once you have the files):

    # analyzer_xyz = HybridizationAnalyzer("molecule.xyz")
    # analyzer_xyz.analyze()

    # analyzer_out = HybridizationAnalyzer("optimization.log")
    # analyzer_out.analyze()

    # analyzer_pdb = HybridizationAnalyzer("protein.pdb")
    # analyzer_pdb.analyze()

    # analyzer_prmtop = HybridizationAnalyzer("complex.prmtop")
    # analyzer_prmtop.analyze()
    pass
