"""Unit tests for the pure-Python ETL functions (no R / DADA2 required).

These are the transforms that the refactor is most likely to have broken,
so they get real fixtures and value assertions rather than smoke tests.
"""

import os
import pandas as pd
import pytest

from tag_analysis import etl


# ---------- create_asv_outputs ----------

def test_create_asv_outputs_writes_counts_fasta_mapping(tmp_path):
    # seqtab_nochim: rows = samples, cols = ASV sequences (as DADA2 emits)
    seqtab = pd.DataFrame(
        {"ACGT": [10, 0], "TTTT": [5, 7]},
        index=["sampleA", "sampleB"],
    )
    seqtab.to_csv(os.path.join(tmp_path, "sequence_table_nochim.csv"))

    mapping_df, asv_tab = etl.create_asv_outputs(str(tmp_path))

    # ASV IDs assigned in column order
    assert list(mapping_df["ASV_ID"]) == ["ASV_1", "ASV_2"]
    assert list(mapping_df["sequence"]) == ["ACGT", "TTTT"]

    # counts table is transposed: rows = ASVs, cols = samples
    assert list(asv_tab.index) == ["ASV_1", "ASV_2"]
    assert list(asv_tab.columns) == ["sampleA", "sampleB"]
    assert asv_tab.loc["ASV_1", "sampleA"] == 10
    assert asv_tab.loc["ASV_2", "sampleB"] == 7

    # FASTA written with headers matching mapping
    fasta = open(os.path.join(tmp_path, "ASVs.fa")).read()
    assert ">ASV_1\nACGT\n" in fasta
    assert ">ASV_2\nTTTT\n" in fasta

    # counts file on disk is tab-separated
    written = pd.read_csv(os.path.join(tmp_path, "ASVs_counts.csv"), sep="\t", index_col=0)
    assert written.loc["ASV_1", "sampleA"] == 10


# ---------- split_taxonomy_to_ranks ----------

def test_split_taxonomy_pads_with_last_taxon():
    row = pd.Series({"taxonomy": "Bacteria;Proteobacteria;Gammaproteobacteria"})
    result = etl.split_taxonomy_to_ranks(row)
    assert result["domain"] == "Bacteria"
    assert result["phylum"] == "Proteobacteria"
    assert result["class"] == "Gammaproteobacteria"
    # padding repeats the last observed taxon down through species
    assert result["order"] == "Gammaproteobacteria"
    assert result["species"] == "Gammaproteobacteria"


def test_split_taxonomy_empty_string_yields_all_blank():
    row = pd.Series({"taxonomy": ""})
    result = etl.split_taxonomy_to_ranks(row)
    assert all(result[r] == "" for r in etl.RANKS)


def test_split_taxonomy_truncates_over_length():
    long_tax = ";".join(["a", "b", "c", "d", "e", "f", "g", "h", "i"])
    row = pd.Series({"taxonomy": long_tax})
    result = etl.split_taxonomy_to_ranks(row)
    assert len(result) == len(etl.RANKS)
    assert result["species"] == "g"  # 7th element, extras dropped


# ---------- prepare_relative_abundance_data ----------

def test_prepare_relative_abundance_data_normalizes_and_melts(tmp_path):
    # Real create_asv_outputs writes an UNNAMED index; prepare_relative_abundance_data
    # relies on reset_index() producing a column literally named 'index'. Match that
    # contract here rather than naming the index.
    counts = pd.DataFrame(
        {"s1": [30, 10], "s2": [0, 20]},
        index=["ASV_1", "ASV_2"],
    )
    counts_path = os.path.join(tmp_path, "counts.csv")
    counts.to_csv(counts_path, sep="\t")

    tax = pd.DataFrame(
        {"domain": ["Bacteria", "Archaea"], "phylum": ["Firmicutes", "Euryarchaeota"]},
        index=["ASV_1", "ASV_2"],
    )
    tax.index.name = "ASV_ID"
    tax_path = os.path.join(tmp_path, "tax.csv")
    tax.to_csv(tax_path, sep="\t")

    long_df = etl.prepare_relative_abundance_data(counts_path, tax_path)

    # long format: one row per ASV x sample
    assert len(long_df) == 4
    # relative abundance: s1 column sums to 40 -> ASV_1 = 75%, ASV_2 = 25%
    s1_asv1 = long_df[(long_df["sample"] == "s1") & (long_df["ASV"] == "ASV_1")]["relabund"].iloc[0]
    assert s1_asv1 == pytest.approx(75.0)
    # taxonomy merged in
    assert long_df[long_df["ASV"] == "ASV_1"]["domain"].iloc[0] == "Bacteria"


# ---------- prepare_data_for_contamination_plot ----------

def test_prepare_contamination_plot_aggregates_by_status():
    relative_df = pd.DataFrame(
        {"realSample": [60.0, 40.0], "controlBlank": [90.0, 10.0]},
        index=["ASV_1", "ASV_2"],
    )
    contam_asvs = ["ASV_2"]
    # predicted_controls aligned to columns: realSample=False, controlBlank=True
    predicted_controls = [False, True]

    agg = etl.prepare_data_for_contamination_plot(relative_df, contam_asvs, predicted_controls)

    # For realSample: contaminant (ASV_2)=40, non-contaminant (ASV_1)=60
    real_contam = agg[(agg["Sample"] == "realSample") & (agg["decontam"])]["relabund"].iloc[0]
    assert real_contam == pytest.approx(40.0)
    real_clean = agg[(agg["Sample"] == "realSample") & (~agg["decontam"])]["relabund"].iloc[0]
    assert real_clean == pytest.approx(60.0)
    # library_type derived from predicted_controls
    assert agg[agg["Sample"] == "controlBlank"]["library_type"].unique().tolist() == ["Control"]
