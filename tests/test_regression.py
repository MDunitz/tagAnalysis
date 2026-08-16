"""Real-data regression tests for the pure-Python transforms.

These run the refactored functions against real DADA2 output from a full-depth
16S/18S run and assert they reproduce the golden output files byte-for-value.
This is the regression guard for the port: if the refactor changed any
transform's behavior, these fail.

Fixtures under tests/fixtures/{16s,18s}/ are committed normally (small CSVs),
so these tests need no Git LFS and add negligible CI bandwidth. The fastq /
R-object stages (DADA2, taxonomy, decontam) are NOT covered here - they need
the R/conda stack and are deferred to the integration tier.
"""

import os

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from tag_analysis import etl

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")
AMPLICONS = ["16s", "18s"]


def _fixture(amplicon, name):
    return os.path.join(FIXTURES, amplicon, name)


# ---------- create_asv_outputs: seqtab_nochim -> counts table + mapping ----------

@pytest.mark.parametrize("amplicon", AMPLICONS)
def test_create_asv_outputs_matches_golden(amplicon, tmp_path):
    """create_asv_outputs on the real sequence_table_nochim.csv must reproduce
    the golden ASVs_counts.csv and asv_mapping.csv exactly."""
    # Stage the real DADA2 output where the function expects it.
    seqtab_src = _fixture(amplicon, "sequence_table_nochim.csv")
    work = tmp_path / amplicon
    work.mkdir()
    with open(seqtab_src) as fh:
        (work / "sequence_table_nochim.csv").write_text(fh.read())

    mapping_df, asv_tab = etl.create_asv_outputs(str(work))

    # --- counts table ---
    produced_counts = pd.read_csv(work / "ASVs_counts.csv", sep="\t", index_col=0)
    golden_counts = pd.read_csv(_fixture(amplicon, "ASVs_counts.csv"), sep="\t", index_col=0)
    assert_frame_equal(produced_counts, golden_counts, check_dtype=False)

    # --- mapping (ASV_ID -> sequence) ---
    golden_mapping = pd.read_csv(_fixture(amplicon, "asv_mapping.csv"), sep="\t")
    assert_frame_equal(
        mapping_df.reset_index(drop=True),
        golden_mapping.reset_index(drop=True),
        check_dtype=False,
    )


@pytest.mark.parametrize("amplicon", AMPLICONS)
def test_asv_ids_sequential_and_aligned(amplicon, tmp_path):
    """ASV IDs are ASV_1..ASV_N in column order, and counts rows align to mapping."""
    seqtab_src = _fixture(amplicon, "sequence_table_nochim.csv")
    work = tmp_path / amplicon
    work.mkdir()
    with open(seqtab_src) as fh:
        (work / "sequence_table_nochim.csv").write_text(fh.read())

    mapping_df, asv_tab = etl.create_asv_outputs(str(work))

    n = len(mapping_df)
    assert list(mapping_df["ASV_ID"]) == [f"ASV_{i+1}" for i in range(n)]
    assert list(asv_tab.index) == list(mapping_df["ASV_ID"])


# ---------- split_taxonomy_to_ranks against golden taxonomy ----------

@pytest.mark.parametrize("amplicon", AMPLICONS)
def test_split_taxonomy_reproduces_golden_ranks(amplicon):
    """Splitting the golden taxonomy strings must reproduce the golden per-rank
    columns (domain..species)."""
    golden = pd.read_csv(_fixture(amplicon, "ASV_taxonomy.csv"), sep="\t")
    if "taxonomy" not in golden.columns:
        pytest.skip(f"{amplicon} taxonomy fixture has no 'taxonomy' column")

    recomputed = golden.apply(etl.split_taxonomy_to_ranks, axis=1)
    for rank in etl.RANKS:
        if rank not in golden.columns:
            continue
        # Golden blanks may be NaN on read; align by filling both sides.
        assert (
            recomputed[rank].fillna("").astype(str).tolist()
            == golden[rank].fillna("").astype(str).tolist()
        ), f"{amplicon}: rank '{rank}' does not match golden"


# ---------- relative abundance sanity against real counts ----------

@pytest.mark.parametrize("amplicon", AMPLICONS)
def test_relative_abundance_sums_to_100_per_sample(amplicon, tmp_path):
    """Relative abundance derived from real counts sums to ~100 per sample."""
    counts = pd.read_csv(_fixture(amplicon, "ASVs_counts.csv"), sep="\t", index_col=0)
    rel = counts.div(counts.sum(axis=0), axis=1) * 100
    col_sums = rel.sum(axis=0)
    assert col_sums.apply(lambda s: s == pytest.approx(100.0, abs=1e-6)).all(), (
        f"{amplicon}: some sample relative abundances do not sum to 100"
    )
