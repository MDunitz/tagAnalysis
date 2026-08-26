"""Tests for helper_functions path generation, RunConfig derived paths,
and the R-execution boundary (mocked - no Rscript in test env)."""

import os
import subprocess
from unittest import mock

import pandas as pd
import pytest

from tag_analysis import helper_functions as hf, pipelines, process
from tag_analysis.config import RunConfig, PrimerSet, PRIMERS_16S, PRIMERS_18S


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


# ---------- configurable primers (RunConfig.primers) ----------

def test_runconfig_primers_defaults_none():
    cfg = RunConfig(data_path="/d", output_path="/o", dataset_name="t",
                    reference_db_path="/r.RData")
    assert cfg.primers is None


def test_process_16s_injects_standard_pair_when_unset():
    cfg = RunConfig(data_path="/d", output_path="/o", dataset_name="t",
                    reference_db_path="/r.RData")
    with mock.patch.object(pipelines, "_run") as m:
        pipelines.process_16s(cfg)
    assert cfg.primers is PRIMERS_16S


def test_process_18s_injects_standard_pair_when_unset():
    cfg = RunConfig(data_path="/d", output_path="/o", dataset_name="t",
                    reference_db_path="/r.RData")
    with mock.patch.object(pipelines, "_run") as m:
        pipelines.process_18s(cfg)
    assert cfg.primers is PRIMERS_18S


def test_explicit_primers_are_respected_by_process_16s():
    """An explicitly set primer pair is NOT overwritten by process_16s."""
    custom = PrimerSet(name="custom", fwd="AAAA", rev="TTTT", fwd_rc="TTTT", rev_rc="AAAA")
    cfg = RunConfig(data_path="/d", output_path="/o", dataset_name="t",
                    reference_db_path="/r.RData", primers=custom)
    with mock.patch.object(pipelines, "_run") as m:
        pipelines.process_16s(cfg)
    assert cfg.primers is custom


def test_process_requires_primers_set():
    """The generic process() raises if no primers are configured."""
    cfg = RunConfig(data_path="/d", output_path="/o", dataset_name="t",
                    reference_db_path="/r.RData")
    with pytest.raises(ValueError, match="config.primers"):
        process(cfg)


def test_custom_primers_flow_to_cutadapt(tmp_path):
    """A custom primer pair on the config reaches remove_primers_cutadapt."""
    custom = PrimerSet(name="v3v4", fwd="CCTACGGG", rev="GACTACHV",
                       fwd_rc="CCCGTAGG", rev_rc="DBGTAGTC")
    cfg = RunConfig(data_path=str(tmp_path / "d"), output_path=str(tmp_path / "o"),
                    dataset_name="t", reference_db_path=str(tmp_path / "r.RData"),
                    primers=custom)
    with mock.patch.object(pipelines, "remove_primers_cutadapt") as m_cut, \
         mock.patch.object(pipelines, "run_dada2_pipeline", return_value={}), \
         mock.patch.object(pipelines, "create_asv_outputs",
                           return_value=(pd.DataFrame(
                               {"ASV_ID": ["ASV_1"], "sequence": ["ACGT"]}), None)), \
         mock.patch.object(pipelines, "assign_taxonomy",
                           return_value=pd.DataFrame()), \
         mock.patch.object(pipelines, "remove_contaminants",
                           return_value=(pd.DataFrame(),
                                         pd.DataFrame(), [], [False])), \
         mock.patch.object(pipelines, "prepare_data_for_contamination_plot",
                           return_value=pd.DataFrame()), \
         mock.patch.object(pipelines, "create_contamination_plot"), \
         mock.patch.object(pipelines, "prepare_relative_abundance_data",
                           return_value=pd.DataFrame()), \
         mock.patch.object(pipelines, "create_relative_abundance_stackbars"), \
         mock.patch.object(pipelines.pd, "read_csv",
                           return_value=pd.DataFrame(
                               {"s1": [1.0]}, index=["ASV_1"])):
        process(cfg)
    args, _ = m_cut.call_args
    assert args[1] == "CCTACGGG" and args[2] == "GACTACHV"
