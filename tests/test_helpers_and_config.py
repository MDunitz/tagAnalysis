"""Tests for helper_functions path generation, RunConfig derived paths,
and the R-execution boundary (mocked - no Rscript in test env)."""

import os
from unittest import mock

import pytest

from tag_analysis import helper_functions as hf
from tag_analysis.config import RunConfig, PRIMERS_16S, PRIMERS_18S


# ---------- generate_file_paths_and_samples ----------

def test_generate_file_paths_and_samples(tmp_path):
    # Create fake paired fastqs matching the Illumina naming the regex expects
    names = [
        "SampleX_S1_L001_R1_001.fastq.gz",
        "SampleX_S1_L001_R2_001.fastq.gz",
        "SampleY_S2_L001_R1_001.fastq.gz",
        "SampleY_S2_L001_R2_001.fastq.gz",
    ]
    for n in names:
        (tmp_path / n).write_text("")

    info = hf.generate_file_paths_and_samples(str(tmp_path))

    assert len(info["forward_reads"]) == 2
    assert len(info["reverse_reads"]) == 2
    # sample names strip _S#_L... suffix
    assert info["samples"] == ["SampleX", "SampleY"]
    # filtered paths derived by substituting _001.fastq -> _filtered.fastq
    assert all("_filtered.fastq" in f for f in info["filtered_forward_reads"])
    # R-formatted strings are quoted, comma-joined
    assert info["r_samples"] == '"SampleX", "SampleY"'


def test_generate_file_paths_empty_dir(tmp_path):
    info = hf.generate_file_paths_and_samples(str(tmp_path))
    assert info["forward_reads"] == []
    assert info["samples"] == []
    assert info["r_samples"] == ""


# ---------- _execute_r_script (mock the subprocess boundary) ----------

def test_execute_r_script_calls_rscript_and_cleans_up():
    with mock.patch("tag_analysis.helper_functions.subprocess.run") as m_run, \
         mock.patch("tag_analysis.helper_functions.os.unlink") as m_unlink:
        hf._execute_r_script("cat('hi')", success_message="ok")
        # Rscript invoked with check=True
        args, kwargs = m_run.call_args
        assert args[0][0] == "Rscript"
        assert kwargs.get("check") is True
        # temp file cleaned up
        assert m_unlink.called


def test_execute_r_script_cleans_up_on_failure():
    import subprocess
    with mock.patch("tag_analysis.helper_functions.subprocess.run",
                    side_effect=subprocess.CalledProcessError(1, "Rscript")), \
         mock.patch("tag_analysis.helper_functions.os.unlink") as m_unlink:
        with pytest.raises(subprocess.CalledProcessError):
            hf._execute_r_script("stop('x')")
        # finally-block still unlinks the temp file
        assert m_unlink.called


# ---------- RunConfig ----------

def test_runconfig_derived_paths():
    cfg = RunConfig(
        data_path="/proj/HelenThesis/16S",
        output_path="/proj/HelenThesis/output/16s",
        dataset_name="HelenThesis",
        reference_db_path="/refs/SILVA.RData",
    )
    assert cfg.deprimered_path == "/proj/HelenThesis/16S/fastq_cutadapt"
    assert cfg.counts_file_path == "/proj/HelenThesis/output/16s/ASVs_counts.csv"
    assert cfg.taxonomy_file_path == "/proj/HelenThesis/output/16s/ASV_taxonomy.csv"
    assert cfg.img_dir == "/proj/HelenThesis/output/16s/imgs"


def test_runconfig_ensure_dirs(tmp_path):
    cfg = RunConfig(
        data_path=str(tmp_path / "data" / "16S"),
        output_path=str(tmp_path / "out"),
        dataset_name="t",
        reference_db_path="/refs/x.RData",
    )
    cfg.ensure_dirs()
    assert os.path.isdir(cfg.output_path)
    assert os.path.isdir(cfg.img_dir)
    assert os.path.isdir(cfg.deprimered_path)


def test_primer_sets_carry_expected_sequences():
    assert PRIMERS_16S.name == "16S"
    assert PRIMERS_16S.fwd.startswith("GTGYCAGC")
    assert PRIMERS_18S.name == "18S"
    assert PRIMERS_18S.fwd.startswith("CCAGCAGC")
