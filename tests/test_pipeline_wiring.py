"""Mock-boundary test for pipeline orchestration.

Verifies process_16s / process_18s call the pipeline stages with the
correct config-derived arguments and correct primer set, without invoking
any external tool. This catches refactor breakage in the wiring layer
(the part most changed by parameterization). It does NOT verify scientific
correctness of the R/DADA2 steps - that requires a real run.
"""

import os
from unittest import mock

import pandas as pd
import pytest

from tag_analysis import pipelines
from tag_analysis.config import RunConfig


@pytest.fixture
def cfg(tmp_path):
    return RunConfig(
        data_path=str(tmp_path / "16S"),
        output_path=str(tmp_path / "out"),
        dataset_name="TestSet",
        reference_db_path="/refs/SILVA.RData",
    )


def _patch_all():
    """Patch every external-touching stage inside pipelines with mocks."""
    patchers = {
        "remove_primers_cutadapt": mock.patch.object(pipelines, "remove_primers_cutadapt"),
        "run_dada2_pipeline": mock.patch.object(pipelines, "run_dada2_pipeline", return_value={}),
        "create_asv_outputs": mock.patch.object(
            pipelines, "create_asv_outputs",
            return_value=(pd.DataFrame({"ASV_ID": ["ASV_1"], "sequence": ["ACGT"]}), None),
        ),
        "assign_taxonomy": mock.patch.object(pipelines, "assign_taxonomy", return_value=pd.DataFrame()),
        "remove_contaminants": mock.patch.object(
            pipelines, "remove_contaminants",
            return_value=(pd.DataFrame(), pd.DataFrame(), [], [False]),
        ),
        "prepare_data_for_contamination_plot": mock.patch.object(
            pipelines, "prepare_data_for_contamination_plot", return_value=pd.DataFrame()),
        "create_contamination_plot": mock.patch.object(pipelines, "create_contamination_plot"),
        "prepare_relative_abundance_data": mock.patch.object(
            pipelines, "prepare_relative_abundance_data", return_value=pd.DataFrame()),
        "create_relative_abundance_stackbars": mock.patch.object(
            pipelines, "create_relative_abundance_stackbars"),
    }
    return patchers


def _run_with_mocks(fn, cfg):
    patchers = _patch_all()
    started = {k: p.start() for k, p in patchers.items()}
    # counts file is read via pd.read_csv between stages; stub it
    read_csv = mock.patch.object(
        pipelines.pd, "read_csv",
        return_value=pd.DataFrame({"s1": [1.0]}, index=["ASV_1"]),
    )
    started["read_csv"] = read_csv.start()
    try:
        result = fn(cfg)
    finally:
        for p in patchers.values():
            p.stop()
        read_csv.stop()
    return started, result


def test_process_16s_uses_16s_primers_and_config_paths(cfg):
    mocks, result = _run_with_mocks(pipelines.process_16s, cfg)

    # primer removal called with 16S forward primer and config's deprimered path
    args, kwargs = mocks["remove_primers_cutadapt"].call_args
    assert args[0] == cfg.data_path
    assert args[1].startswith("GTGYCAGC")  # 16S FWD
    assert kwargs["output_directory"] == cfg.deprimered_path

    # dada2 called with config's output + dataset name
    _, d_kwargs = mocks["run_dada2_pipeline"].call_args
    assert d_kwargs["path_to_output_dir"] == cfg.output_path
    assert d_kwargs["dataset_name"] == "TestSet"

    # taxonomy called with the config's reference DB (not a hardcoded global)
    _, t_kwargs = mocks["assign_taxonomy"].call_args
    assert t_kwargs["reference_db_path"] == "/refs/SILVA.RData"

    # output dirs actually created
    assert os.path.isdir(cfg.output_path)


def test_process_18s_uses_18s_primers(cfg):
    mocks, result = _run_with_mocks(pipelines.process_18s, cfg)
    args, _ = mocks["remove_primers_cutadapt"].call_args
    assert args[1].startswith("CCAGCAGC")  # 18S FWD, differs from 16S


def test_full_stage_sequence_invoked(cfg):
    """Every stage runs exactly once (no stage silently dropped in refactor)."""
    mocks, _ = _run_with_mocks(pipelines.process_16s, cfg)
    for stage in ["remove_primers_cutadapt", "run_dada2_pipeline", "create_asv_outputs",
                  "assign_taxonomy", "remove_contaminants", "create_contamination_plot",
                  "create_relative_abundance_stackbars"]:
        assert mocks[stage].call_count == 1, f"{stage} called {mocks[stage].call_count}x"
