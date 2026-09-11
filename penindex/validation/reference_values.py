"""Literature reference values for regression-testing the PI engine.

Sources:
- Echeverria, J., & Alvarez, S. (2023). Chem. Sci., 14, 11647-11688.
- Echeverria, J., & Alvarez, S. (2024). Chem. Sci., 15, 12166-12168.
"""

from typing import List, NamedTuple


class ExactPIReference(NamedTuple):
    """A single (element_1, element_2, distance) -> PI value that is exactly
    reproducible from the formula and radii table alone (a specific
    structure/optimized geometry reported in the literature, not a survey
    average)."""
    description: str
    element_1: str
    element_2: str
    distance_angstrom: float
    expected_pi_percent: float
    tolerance_percent: float
    source: str


# Both values below were hand-verified against penindex.core.compute_pi_for_elements
# while writing this module: He2 -> -4.35% (paper: -4.3%), H2+ -> 75.4% (paper: 75%).
EXACT_VALUES: List[ExactPIReference] = [
    ExactPIReference(
        description="He2 dimer at the DFT-optimized interatomic distance",
        element_1="HE", element_2="HE",
        distance_angstrom=2.96,
        expected_pi_percent=-4.3,
        tolerance_percent=0.2,
        source="Echeverria & Alvarez, Chem. Sci. 2023, p. 11650: 'High-level calculations "
               "place the minimal energy for the He2 dimer at 2.96 A, i.e., a penetration "
               "index of -4.3%.'",
    ),
    ExactPIReference(
        description="H2+ dihydrogen cation, gas-phase interatomic distance",
        element_1="H", element_2="H",
        distance_angstrom=1.058,
        expected_pi_percent=75.0,
        tolerance_percent=0.5,
        source="Echeverria & Alvarez, Chem. Sci. 2023, p. 11664: 'Its interatomic distance "
               "determined in the gas phase is 1.058 A, that corresponds to a penetration "
               "index of 75%.'",
    ),
]


class RangeReference(NamedTuple):
    """A literature mean(std) range from a CSD structural survey - not
    reproducible from a single distance, but useful as a sanity-range check
    on aggregated occupancy-weighted PI statistics for the same bond type."""
    description: str
    mean_percent: float
    std_percent: float
    source: str


CSD_HBOND_RANGES: List[RangeReference] = [
    RangeReference(
        description="O-H...N hydrogen bonds, H...N penetration (CSD survey average)",
        mean_percent=47, std_percent=7,
        source="Echeverria & Alvarez, Chem. Sci. 2023, p. 11652.",
    ),
    RangeReference(
        description="O-H...N hydrogen bonds, O...N penetration (CSD survey average)",
        mean_percent=18, std_percent=6,
        source="Echeverria & Alvarez, Chem. Sci. 2023, p. 11652.",
    ),
    RangeReference(
        description="O-H...O hydrogen bonds, H...O penetration (CSD survey average)",
        mean_percent=42, std_percent=8,
        source="Echeverria & Alvarez, Chem. Sci. 2023, p. 11652.",
    ),
    RangeReference(
        description="O-H...O hydrogen bonds, O...O penetration (CSD survey average)",
        mean_percent=13, std_percent=6,
        source="Echeverria & Alvarez, Chem. Sci. 2023, p. 11652.",
    ),
]


class CrustWidthRatioReference(NamedTuple):
    """A literature (element_1, element_2) -> van der Waals crust-width
    ratio, exactly reproducible from the radii table alone (w = v - r for
    each element; ratio = max(w1, w2) / min(w1, w2))."""
    description: str
    element_1: str
    element_2: str
    expected_ratio: float
    tolerance: float
    source: str


CRUST_WIDTH_RATIOS: List[CrustWidthRatioReference] = [
    CrustWidthRatioReference(
        description="Ta-O pair, near-identical crust widths",
        element_1="TA", element_2="O",
        expected_ratio=1.01,
        tolerance=0.02,
        source="Echeverria & Alvarez, Chem. Sci. 2024, p. 12168: 'wTa = 0.83, wO = 0.84 A' "
               "(Fig. 3 caption), crust-width ratio 1.01.",
    ),
    CrustWidthRatioReference(
        description="Zn-Te pair, notably mismatched crust widths",
        element_1="ZN", element_2="TE",
        expected_ratio=1.92,
        tolerance=0.02,
        source="Echeverria & Alvarez, Chem. Sci. 2024, p. 12168: 'wZn = 1.17 Å, wTe = 0.61 A' "
               "(Fig. 3 caption), crust-width ratio 1.92.",
    ),
]
