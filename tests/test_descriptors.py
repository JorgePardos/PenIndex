import pandas as pd
import pytest

from penindex.exceptions import PenIndexError
from penindex.io.descriptors import join_batch_with_descriptors, load_descriptor_table


def test_load_descriptor_table_requires_label_column(tmp_path):
    csv_path = tmp_path / "descriptors.csv"
    csv_path.write_text("system,value\nA,1.0\nB,2.0\n")

    with pytest.raises(PenIndexError):
        load_descriptor_table(str(csv_path))


def test_load_descriptor_table_reads_valid_csv(tmp_path):
    csv_path = tmp_path / "descriptors.csv"
    csv_path.write_text("label,delta_E_CT\nA,1.5\nB,2.5\n")

    df = load_descriptor_table(str(csv_path))

    assert list(df['label']) == ['A', 'B']
    assert list(df['delta_E_CT']) == [1.5, 2.5]


def test_join_batch_with_descriptors_inner_join():
    batch_df = pd.DataFrame({'label': ['A', 'B', 'C'], 'PI(%)': [10.0, 20.0, 30.0]})
    descriptor_df = pd.DataFrame({'label': ['A', 'B'], 'delta_E_CT': [1.5, 2.5]})

    merged = join_batch_with_descriptors(batch_df, descriptor_df)

    assert list(merged['label']) == ['A', 'B']
    assert 'delta_E_CT' in merged.columns


def test_join_batch_with_descriptors_raises_on_no_overlap():
    batch_df = pd.DataFrame({'label': ['A'], 'PI(%)': [10.0]})
    descriptor_df = pd.DataFrame({'label': ['Z'], 'delta_E_CT': [1.5]})

    with pytest.raises(PenIndexError):
        join_batch_with_descriptors(batch_df, descriptor_df)
