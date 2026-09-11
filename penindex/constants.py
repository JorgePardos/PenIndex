"""Scientific constants for the Penetration Index engine.

Covalent (r) and van der Waals (v) radii in Angstroms, from Alvarez et al.
(Cordero et al., Dalton Trans. 2008, 2832 for r; Alvarez, Dalton Trans. 2013,
8617 for v), the same empirical sets referenced by the Penetration Index
papers (Echeverria & Alvarez, Chem. Sci. 2023, 2024; Pardos et al., Dalton
Trans. 2026).
"""

from typing import Dict

# 'C' uses the sp2 value as a generic default; overridden per-atom when
# dynamic_carbon_radii is enabled (see io.hybridization / core.apply_dynamic_carbon_radii).
# 'MN'/'FE'/'CO' default to the high-spin radii; override per-system if the
# metal centre is known to be low-spin.
RADII: Dict[str, Dict[str, float]] = {
    'H':  {'r': 0.31, 'v': 1.20}, 'HE': {'r': 0.28, 'v': 1.43},
    'LI': {'r': 1.28, 'v': 2.12}, 'BE': {'r': 0.96, 'v': 1.98},
    'B':  {'r': 0.84, 'v': 1.91}, 'C':  {'r': 0.73, 'v': 1.77},
    'N':  {'r': 0.71, 'v': 1.66}, 'O':  {'r': 0.66, 'v': 1.50},
    'F':  {'r': 0.57, 'v': 1.46}, 'NE': {'r': 0.58, 'v': 1.58},
    'NA': {'r': 1.66, 'v': 2.50}, 'MG': {'r': 1.41, 'v': 2.51},
    'AL': {'r': 1.21, 'v': 2.25}, 'SI': {'r': 1.11, 'v': 2.19},
    'P':  {'r': 1.07, 'v': 1.90}, 'S':  {'r': 1.05, 'v': 1.89},
    'CL': {'r': 1.02, 'v': 1.82}, 'AR': {'r': 1.06, 'v': 1.94},
    'K':  {'r': 2.03, 'v': 2.73}, 'CA': {'r': 1.76, 'v': 2.62},
    'SC': {'r': 1.70, 'v': 2.58}, 'TI': {'r': 1.60, 'v': 2.46},
    'V':  {'r': 1.53, 'v': 2.42}, 'CR': {'r': 1.39, 'v': 2.45},
    'MN': {'r': 1.61, 'v': 2.45}, 'FE': {'r': 1.52, 'v': 2.44},
    'CO': {'r': 1.50, 'v': 2.40}, 'NI': {'r': 1.24, 'v': 2.40},
    'CU': {'r': 1.32, 'v': 2.38}, 'ZN': {'r': 1.22, 'v': 2.39},
    'GA': {'r': 1.22, 'v': 2.32}, 'GE': {'r': 1.20, 'v': 2.29},
    'AS': {'r': 1.19, 'v': 1.88}, 'SE': {'r': 1.20, 'v': 1.82},
    'BR': {'r': 1.20, 'v': 1.86}, 'KR': {'r': 1.16, 'v': 2.07},
    'RB': {'r': 2.20, 'v': 3.21}, 'SR': {'r': 1.95, 'v': 2.84},
    'Y':  {'r': 1.90, 'v': 2.75}, 'ZR': {'r': 1.75, 'v': 2.52},
    'NB': {'r': 1.64, 'v': 2.56}, 'MO': {'r': 1.54, 'v': 2.45},
    'TC': {'r': 1.47, 'v': 2.44}, 'RU': {'r': 1.46, 'v': 2.46},
    'RH': {'r': 1.42, 'v': 2.44}, 'PD': {'r': 1.39, 'v': 2.15},
    'AG': {'r': 1.45, 'v': 2.53}, 'CD': {'r': 1.44, 'v': 2.49},
    'IN': {'r': 1.42, 'v': 2.43}, 'SN': {'r': 1.39, 'v': 2.42},
    'SB': {'r': 1.39, 'v': 2.47}, 'TE': {'r': 1.38, 'v': 1.99},
    'I':  {'r': 1.39, 'v': 2.04}, 'XE': {'r': 1.40, 'v': 2.28},
    'CS': {'r': 2.44, 'v': 3.48}, 'BA': {'r': 2.15, 'v': 3.03},
    'LA': {'r': 2.07, 'v': 2.98}, 'CE': {'r': 2.04, 'v': 2.88},
    'PR': {'r': 2.03, 'v': 2.92}, 'ND': {'r': 2.01, 'v': 2.95},
    'SM': {'r': 1.98, 'v': 2.90}, 'EU': {'r': 1.98, 'v': 2.87},
    'GD': {'r': 1.96, 'v': 2.83}, 'TB': {'r': 1.94, 'v': 2.79},
    'DY': {'r': 1.92, 'v': 2.87}, 'HO': {'r': 1.92, 'v': 2.81},
    'ER': {'r': 1.89, 'v': 2.83}, 'TM': {'r': 1.90, 'v': 2.79},
    'YB': {'r': 1.87, 'v': 2.80}, 'LU': {'r': 1.87, 'v': 2.74},
    'HF': {'r': 1.75, 'v': 2.63}, 'TA': {'r': 1.70, 'v': 2.53},
    'W':  {'r': 1.62, 'v': 2.57}, 'RE': {'r': 1.51, 'v': 2.49},
    'OS': {'r': 1.44, 'v': 2.48}, 'IR': {'r': 1.41, 'v': 2.41},
    'PT': {'r': 1.36, 'v': 2.29}, 'AU': {'r': 1.36, 'v': 2.32},
    'HG': {'r': 1.32, 'v': 2.45}, 'TL': {'r': 1.45, 'v': 2.47},
    'PB': {'r': 1.46, 'v': 2.60}, 'BI': {'r': 1.48, 'v': 2.54},
    'RN': {'r': 1.50, 'v': 2.40}, 'AC': {'r': 2.15, 'v': 2.80},
    'TH': {'r': 2.06, 'v': 2.93}, 'PA': {'r': 2.00, 'v': 2.88},
    'U':  {'r': 1.96, 'v': 2.71}, 'NP': {'r': 1.90, 'v': 2.82},
    'PU': {'r': 1.97, 'v': 2.81}, 'AM': {'r': 1.80, 'v': 2.83},
    'CM': {'r': 1.69, 'v': 3.05}, 'CF': {'r': 1.81, 'v': 3.05},
    'ES': {'r': 1.88, 'v': 2.70},
}

# Typical donor-hydrogen bond lengths (Angstroms), used to reconstruct the
# exact H...Acceptor distance from the reported Donor-Acceptor distance and
# D-H-A angle via the law of cosines (see core.post_process_hbonds).
D_H_BONDS: Dict[str, float] = {
    'O': 0.96, 'N': 1.01, 'S': 1.34, 'C': 1.09
}

# Carbon covalent radius depends on hybridization (Cordero et al., Dalton
# Trans. 2008); the van der Waals radius is treated as constant across
# hybridizations, so only 'r' changes when dynamic_carbon_radii overrides the
# flat RADII['C'] default.
CARBON_HYBRID_COV_RADII: Dict[str, float] = {
    'SP': 0.69, 'SP2': 0.73, 'SP3': 0.76,
}
