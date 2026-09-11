"""Loads externally computed per-bond covalency descriptors (from AIMAll/IQA,
Q-Chem ALMO-EDA, DDEC6, generalized-Badger-rule bond-strength orders, or any
other source) so they can be correlated against the Penetration Index
computed by penindex.batch for the same systems - the core workflow of
Pardos, Gonzalo, Merino & Echeverria, Dalton Trans. 2026.
"""

import pandas as pd

from penindex.exceptions import PenIndexError


def load_descriptor_table(path: str, label_column: str = "label") -> pd.DataFrame:
    """Loads a CSV of external descriptors.

    Must contain `label_column` (matching the `label` used in the batch
    systems list from penindex.batch) plus one or more numeric descriptor
    columns with whatever names you choose (e.g. delta_E_SC, V_XC,
    delta_E_CT, DDEC6_BO, pct_3c4e).
    """
    try:
        df = pd.read_csv(path)
    except Exception as e:
        raise PenIndexError(f"Could not read descriptor table '{path}': {e}") from e

    if label_column not in df.columns:
        raise PenIndexError(
            f"Descriptor table '{path}' has no '{label_column}' column to join on "
            f"(found columns: {list(df.columns)})"
        )
    return df


def join_batch_with_descriptors(
    batch_df: pd.DataFrame,
    descriptor_df: pd.DataFrame,
    label_column: str = "label",
) -> pd.DataFrame:
    """Inner-joins a penindex.batch.run_batch() result with a descriptor
    table on label_column. Systems present in only one of the two tables
    are silently dropped (mismatched labels are the most common data-entry
    error here, so this returns the overlap rather than raising)."""
    merged = batch_df.merge(descriptor_df, on=label_column, how='inner', suffixes=('', '_descriptor'))
    if len(merged) == 0:
        raise PenIndexError(
            f"No matching '{label_column}' values between the batch results and the "
            f"descriptor table - check that both use the exact same system labels."
        )
    return merged
