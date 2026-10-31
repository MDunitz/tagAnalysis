"""Batch-aware decontam: extraction sheet -> library batches -> per-batch calls."""

from unittest import mock

import numpy as np
import pandas as pd
import pytest

from tag_analysis import decontaminate
from tag_analysis.extraction_metadata import (
    assign_extraction_batches,
    google_sheet_csv_url,
    library_tag_id,
    parse_tag_ids,
    read_extraction_sheet,
    tag_id_to_batch,
)

SHEET_ID = "1P0FmwyJ3e6OaDkcs_5RxYZqTm7p2SfeNE3L1koMkqbo"

LIBRARIES = [
    "GP_1_16S_PCR_Neg_S97_R1_filtered.fastq.gz",
    "GP_2_Neg_1_S98_R1_filtered.fastq.gz",
    "GP_3_Exp03_B01_T0_S99_R1_filtered.fastq.gz",
    "GP_12_Neg_2_S108_R1_filtered.fastq.gz",
    "GP_13_Exp04_B01_T0_S109_R1_filtered.fastq.gz",
]
IS_CONTROL = [True, True, False, True, False]


@pytest.fixture
def sheet():
    return pd.DataFrame(
        {
            "Extraction Batch ID": [1, 1, 2, 2, 2],
            "Sample ID": [None, "Exp03_B01", None, "Exp04_B01", "Exp04_B02"],
            "TAG ID": ["2, 51", "3;52", "12 61", "13", None],
        }
    )


def test_google_sheet_csv_url_keeps_gid():
    url = f"https://docs.google.com/spreadsheets/d/{SHEET_ID}/edit?gid=123#gid=123"
    assert google_sheet_csv_url(url) == (
        f"https://docs.google.com/spreadsheets/d/{SHEET_ID}/export?format=csv&gid=123"
    )


def test_google_sheet_csv_url_defaults_to_first_tab():
    url = f"https://docs.google.com/spreadsheets/d/{SHEET_ID}/edit"
    assert google_sheet_csv_url(url).endswith("export?format=csv&gid=0")


def test_read_extraction_sheet_rewrites_google_url():
    url = f"https://docs.google.com/spreadsheets/d/{SHEET_ID}/edit?gid=0#gid=0"
    with mock.patch.object(pd, "read_csv", return_value=pd.DataFrame()) as m:
        read_extraction_sheet(url)
    assert m.call_args.args[0] == google_sheet_csv_url(url)


def test_parse_tag_ids_accepts_mixed_separators():
    assert parse_tag_ids("3, 52") == [3, 52]
    assert parse_tag_ids("12 61;7") == [12, 61, 7]
    assert parse_tag_ids("13") == [13]


def test_parse_tag_ids_numeric_cell():
    """A TAG ID column with blanks reads back as float (13.0)."""
    assert parse_tag_ids(13.0) == [13]


def test_tag_id_to_batch_maps_every_listed_tag(sheet):
    assert tag_id_to_batch(sheet) == {2: 1, 51: 1, 3: 1, 52: 1, 12: 2, 61: 2, 13: 2}


def test_tag_id_to_batch_rejects_tag_in_two_batches(sheet):
    sheet.loc[4, "TAG ID"] = "3"
    with pytest.raises(ValueError, match="TAG ID 3"):
        tag_id_to_batch(sheet)


def test_library_tag_id():
    assert library_tag_id(LIBRARIES[2]) == 3
    assert library_tag_id("algaebrick_B06") is None


def test_assign_batches_pcr_negative_is_shared(sheet):
    batches = assign_extraction_batches(LIBRARIES, IS_CONTROL, tag_id_to_batch(sheet))
    assert batches == [None, 1, 1, 2, 2]


def test_assign_batches_rejects_unplaced_sample(sheet):
    libs = LIBRARIES + ["GP_14_Exp04_B02_T0_S110_R1_filtered.fastq.gz"]
    with pytest.raises(ValueError, match="GP_14_Exp04_B02"):
        assign_extraction_batches(libs, IS_CONTROL + [False], tag_id_to_batch(sheet))


def test_combine_batch_p_values_takes_minimum_and_skips_nan():
    results = {
        1: pd.DataFrame(
            {"ASV": ["ASV_1", "ASV_2", "ASV_3"], "p_prev": [0.9, 0.05, np.nan]}
        ),
        2: pd.DataFrame(
            {"ASV": ["ASV_1", "ASV_2", "ASV_3"], "p_prev": [0.3, 0.8, 0.6]}
        ),
    }
    out = decontaminate.combine_batch_p_values(results, threshold=0.1).set_index("ASV")
    assert out["p_prev"].tolist() == [0.3, 0.05, 0.6]
    assert out["contaminant"].tolist() == [False, True, False]
    assert out.loc["ASV_3", "p_prev_batch_1"] != out.loc["ASV_3", "p_prev_batch_1"]


def test_by_batch_runs_each_batch_against_own_and_shared_controls(tmp_path):
    counts = pd.DataFrame(
        np.ones((2, len(LIBRARIES)), dtype=int),
        index=["ASV_1", "ASV_2"],
        columns=LIBRARIES,
    )
    seen = {}

    def fake_prevalence(path, results_path, controls, threshold):
        seen[path] = (list(pd.read_csv(path, sep="\t", index_col=0).columns), controls)
        return pd.DataFrame({"ASV": ["ASV_1", "ASV_2"], "p_prev": [0.5, 0.5]})

    with mock.patch.object(
        decontaminate, "run_decontam_prevalence", side_effect=fake_prevalence
    ):
        decontaminate.run_decontam_prevalence_by_batch(
            counts, str(tmp_path), IS_CONTROL, [None, 1, 1, 2, 2], 0.1
        )
    batch1 = seen[str(tmp_path / "ASVs_counts_batch_1.csv")]
    batch2 = seen[str(tmp_path / "ASVs_counts_batch_2.csv")]
    assert batch1 == ([LIBRARIES[1], LIBRARIES[2], LIBRARIES[0]], [True, False, True])
    assert batch2 == ([LIBRARIES[3], LIBRARIES[4], LIBRARIES[0]], [True, False, True])


def test_by_batch_rejects_batch_without_any_control(tmp_path):
    counts = pd.DataFrame([[1, 1]], index=["ASV_1"], columns=["GP_3_a", "GP_13_b"])
    with pytest.raises(ValueError, match="batch 1 has no negative control"):
        decontaminate.run_decontam_prevalence_by_batch(
            counts, str(tmp_path), [False, False], [1, 2], 0.1
        )


def test_remove_contaminants_rejects_sheet_with_conc(tmp_path):
    with pytest.raises(ValueError, match="mutually exclusive"):
        decontaminate.remove_contaminants(
            "c.csv", "t.csv", str(tmp_path), conc=[1.0], extraction_sheet="s.csv"
        )


def test_remove_contaminants_uses_extraction_batches(tmp_path, sheet):
    counts = pd.DataFrame(
        [[5, 1, 9, 1, 8], [0, 3, 0, 4, 0]], index=["ASV_1", "ASV_2"], columns=LIBRARIES
    )
    counts_path, tax_path, sheet_path = (
        tmp_path / n for n in ("c.csv", "t.csv", "s.csv")
    )
    counts.to_csv(counts_path, sep="\t")
    pd.DataFrame({"ASV_ID": ["ASV_1", "ASV_2"], "phylum": ["A", "B"]}).to_csv(
        tax_path, sep="\t", index=False
    )
    sheet.to_csv(sheet_path, index=False)

    def fake_prevalence(path, results_path, controls, threshold):
        return pd.DataFrame({"ASV": ["ASV_1", "ASV_2"], "p_prev": [0.9, 0.01]})

    with mock.patch.object(
        decontaminate, "run_decontam_prevalence", side_effect=fake_prevalence
    ):
        clean, _, contam, _ = decontaminate.remove_contaminants(
            str(counts_path),
            str(tax_path),
            str(tmp_path),
            extraction_sheet=str(sheet_path),
        )
    assert contam == ["ASV_2"]
    assert list(clean.index) == ["ASV_1"]
    assert (tmp_path / "contam_results.csv").exists()
