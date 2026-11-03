"""redecontam: rerun decontam on persisted outputs without the upstream stages."""

import json
from unittest import mock

import pandas as pd
import pytest

from tag_analysis import decontaminate, redecontam
from tag_analysis.manifest import MANIFEST_FILE

LIBS = [
    "GP_1_16S_PCR_Neg_S97_R1_filtered.fastq.gz",
    "GP_2_Neg_1_S98_R1_filtered.fastq.gz",
    "GP_3_Exp03_B01_T0_S99_R1_filtered.fastq.gz",
]


@pytest.fixture
def run_dir(tmp_path):
    pd.DataFrame([[5, 1, 9], [0, 3, 0]], index=["ASV_1", "ASV_2"], columns=LIBS).to_csv(
        tmp_path / "ASVs_counts.csv", sep="\t"
    )
    pd.DataFrame(
        {"ASV_ID": ["ASV_1", "ASV_2"], "domain": ["Bacteria"] * 2, "phylum": ["A", "B"]}
    ).to_csv(tmp_path / "ASV_taxonomy.csv", sep="\t", index=False)
    pd.DataFrame(
        {"library": ["GP_2_Neg_1_S98", "GP_3_Exp03_B01_T0_S99"], "extraction_batch": [1, 1]}
    ).to_csv(tmp_path / "library_batches.csv", index=False)
    (tmp_path / MANIFEST_FILE).write_text(json.dumps({"reference": {"file": "SILVA.RData"}}))
    return tmp_path


def _fake_prevalence(path, results_path, controls, threshold):
    """Shape of run_decontam_prevalence's R output (contaminant at threshold 0.1)."""
    return pd.DataFrame(
        {"ASV": ["ASV_1", "ASV_2"], "contaminant": [False, True], "p_prev": [0.9, 0.01]}
    )


def test_redecontam_rewrites_clean_counts_and_keeps_upstream_manifest(run_dir):
    with mock.patch.object(decontaminate, "run_decontam_prevalence", side_effect=_fake_prevalence):
        contam = redecontam.redecontaminate(
            str(run_dir), str(run_dir / "library_batches.csv"), levels=["phylum"]
        )
    assert contam == ["ASV_2"]
    clean = pd.read_csv(run_dir / "ASVs_counts_clean.csv", sep="\t", index_col=0)
    assert list(clean.index) == ["ASV_1"]
    assert (run_dir / "imgs" / "Contamination_plot.html").exists()

    manifest = json.loads((run_dir / MANIFEST_FILE).read_text())
    assert manifest["reference"] == {"file": "SILVA.RData"}
    decon = manifest["decontamination"]
    assert decon["contaminant_asvs"] == ["ASV_2"]
    assert decon["library_batches"]["file"] == "library_batches.csv"
    assert "rerun_at" in decon and "rerun_tag_analysis_version" in decon


def test_redecontam_without_manifest_writes_decontamination_only(run_dir):
    (run_dir / MANIFEST_FILE).unlink()
    with mock.patch.object(decontaminate, "run_decontam_prevalence", side_effect=_fake_prevalence):
        redecontam.redecontaminate(str(run_dir), levels=["phylum"])
    manifest = json.loads((run_dir / MANIFEST_FILE).read_text())
    assert list(manifest) == ["decontamination"]
    assert manifest["decontamination"]["library_batches"] is None


def test_redecontam_cli_passes_batch_table(run_dir):
    with mock.patch.object(redecontam, "redecontaminate", return_value=[]) as m:
        redecontam.main([str(run_dir), "--library-batches", "b.csv", "--levels", "phylum"])
    assert m.call_args.args[:3] == (str(run_dir), "b.csv", ["phylum"])
