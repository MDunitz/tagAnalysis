"""Batch-aware decontam: library-batch table -> library batches -> per-batch calls."""

from unittest import mock

import numpy as np
import pandas as pd
import pytest

from tag_analysis import decontaminate
from tag_analysis.library_batches import (
    assign_extraction_batches,
    library_id,
    read_library_batches,
)

LIBRARIES = [
    "GP_1_16S_PCR_Neg_S97_R1_filtered.fastq.gz",
    "GP_2_Neg_1_S98_R1_filtered.fastq.gz",
    "GP_3_Exp03_B01_T0_S99_R1_filtered.fastq.gz",
    "GP_12_Neg_2_S108_R1_filtered.fastq.gz",
    "GP_13_Exp04_B01_T0_S109_R1_filtered.fastq.gz",
]
IS_CONTROL = [True, True, False, True, False]


@pytest.fixture
def table_path(tmp_path):
    path = tmp_path / "library_batches.csv"
    pd.DataFrame(
        {
            "library": [
                "GP_1_16S_PCR_Neg_S97",
                "GP_2_Neg_1_S98",
                "GP_3_Exp03_B01_T0_S99",
                "GP_12_Neg_2_S108",
                "GP_13_Exp04_B01_T0_S109",
            ],
            "extraction_batch": [None, 1, 1, 2, 2],
        }
    ).to_csv(path, index=False)
    return path


def test_library_id_strips_raw_and_filtered_suffixes():
    assert library_id("GP_3_Exp03_B01_T0_S99_R1_001.fastq.gz") == "GP_3_Exp03_B01_T0_S99"
    assert library_id(LIBRARIES[2]) == "GP_3_Exp03_B01_T0_S99"


def test_read_library_batches_blank_batch_is_shared(table_path):
    batches = read_library_batches(table_path)
    assert batches["GP_1_16S_PCR_Neg_S97"] is None
    assert batches["GP_13_Exp04_B01_T0_S109"] == 2


def test_read_library_batches_rejects_duplicate_library(tmp_path):
    path = tmp_path / "dup.csv"
    pd.DataFrame({"library": ["a", "a"], "extraction_batch": [1, 2]}).to_csv(path, index=False)
    with pytest.raises(ValueError, match="more than once"):
        read_library_batches(path)


def test_assign_batches_matches_columns_by_library_id(table_path):
    batches = assign_extraction_batches(LIBRARIES, IS_CONTROL, read_library_batches(table_path))
    assert batches == [None, 1, 1, 2, 2]


def test_assign_batches_control_missing_from_table_is_shared(table_path):
    libs = LIBRARIES + ["GP_49_18S_PCR_Neg_S144_R1_filtered.fastq.gz"]
    batches = assign_extraction_batches(libs, IS_CONTROL + [True], read_library_batches(table_path))
    assert batches[-1] is None


def test_assign_batches_rejects_unplaced_sample(table_path):
    libs = LIBRARIES + ["GP_14_Exp04_B02_T0_S110_R1_filtered.fastq.gz"]
    with pytest.raises(ValueError, match="GP_14_Exp04_B02"):
        assign_extraction_batches(libs, IS_CONTROL + [False], read_library_batches(table_path))


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


def test_remove_contaminants_rejects_table_with_conc(tmp_path):
    with pytest.raises(ValueError, match="mutually exclusive"):
        decontaminate.remove_contaminants(
            "c.csv", "t.csv", str(tmp_path), conc=[1.0], library_batches="b.csv"
        )


def test_remove_contaminants_uses_library_batches(tmp_path, table_path):
    counts = pd.DataFrame(
        [[5, 1, 9, 1, 8], [0, 3, 0, 4, 0]], index=["ASV_1", "ASV_2"], columns=LIBRARIES
    )
    counts_path, tax_path = (tmp_path / n for n in ("c.csv", "t.csv"))
    counts.to_csv(counts_path, sep="\t")
    pd.DataFrame({"ASV_ID": ["ASV_1", "ASV_2"], "phylum": ["A", "B"]}).to_csv(
        tax_path, sep="\t", index=False
    )

    def fake_prevalence(path, results_path, controls, threshold):
        return pd.DataFrame({"ASV": ["ASV_1", "ASV_2"], "p_prev": [0.9, 0.01]})

    with mock.patch.object(
        decontaminate, "run_decontam_prevalence", side_effect=fake_prevalence
    ):
        clean, _, contam, _ = decontaminate.remove_contaminants(
            str(counts_path),
            str(tax_path),
            str(tmp_path),
            library_batches=str(table_path),
        )
    assert contam == ["ASV_2"]
    assert list(clean.index) == ["ASV_1"]
    assert (tmp_path / "contam_results.csv").exists()
