"""Run manifest: records what produced a set of pipeline outputs."""

import hashlib
import json

import pandas as pd
import pytest

from tag_analysis import __version__, constants, decontaminate
from tag_analysis.config import PRIMERS_16S, RunConfig
from tag_analysis.manifest import MANIFEST_FILE, build_manifest, file_sha256, write_manifest

REFERENCE_BYTES = b"not a real training set"
CLEAN_COUNT_FILE = "ASVs_counts_clean.csv"


@pytest.fixture
def cfg(tmp_path):
    reference = tmp_path / "SILVA_SSU_r138_2_2024.RData"
    reference.write_bytes(REFERENCE_BYTES)
    config = RunConfig(
        data_path=str(tmp_path / "16S"),
        output_path=str(tmp_path / "out"),
        dataset_name="TestSet",
        reference_db_path=str(reference),
        primers=PRIMERS_16S,
    )
    config.ensure_dirs()
    return config


@pytest.fixture
def counts_df():
    return pd.DataFrame(
        {"sampleA": [10, 0], "PCRBlank": [1, 3], "sampleB": [5, 7]},
        index=["ASV_1", "ASV_2"],
    )


def _manifest(cfg, counts_df, contam_asvs=("ASV_2",)):
    return build_manifest(
        cfg, {"truncLen": (230, 200)}, counts_df, list(contam_asvs),
        [False, True, False], CLEAN_COUNT_FILE,
    )


def test_file_sha256_matches_hashlib(cfg):
    assert file_sha256(cfg.reference_db_path) == hashlib.sha256(REFERENCE_BYTES).hexdigest()


def test_manifest_records_version_gene_and_reference(cfg, counts_df):
    manifest = _manifest(cfg, counts_df)
    assert manifest["tag_analysis_version"] == __version__
    assert manifest["gene"] == "16S"
    assert manifest["primers"] == {"forward": PRIMERS_16S.fwd, "reverse": PRIMERS_16S.rev}
    assert manifest["reference"] == {
        "file": "SILVA_SSU_r138_2_2024.RData",
        "sha256": hashlib.sha256(REFERENCE_BYTES).hexdigest(),
    }
    assert manifest["taxonomy_confidence_threshold"] == constants.TAXONOMY_CONFIDENCE_THRESHOLD


def test_manifest_records_decontamination(cfg, counts_df):
    decontamination = _manifest(cfg, counts_df)["decontamination"]
    assert decontamination == {
        "prevalence_threshold": constants.DECONTAM_PREVALENCE_THRESHOLD,
        "control_libraries": ["PCRBlank"],
        "contaminant_asvs": ["ASV_2"],
    }


def test_manifest_names_output_files(cfg, counts_df):
    assert _manifest(cfg, counts_df)["outputs"] == {
        "counts": "ASVs_counts.csv",
        "clean_counts": CLEAN_COUNT_FILE,
        "taxonomy": "ASV_taxonomy.csv",
    }


def test_write_manifest_round_trips_as_json(cfg, counts_df):
    path = write_manifest(cfg, _manifest(cfg, counts_df))
    assert path.endswith(MANIFEST_FILE)
    with open(path) as handle:
        written = json.load(handle)
    assert written["dada2"] == {"truncLen": [230, 200]}
    assert written["dataset_name"] == "TestSet"


def test_clean_counts_written_when_no_contaminants(cfg, counts_df):
    """Without controls nothing is removed, but the clean file must still exist:
    the stackbar step and the manifest both point at it."""
    no_controls = counts_df.rename(columns={"PCRBlank": "sampleC"})
    no_controls.to_csv(cfg.counts_file_path, sep="\t")
    pd.DataFrame({"taxonomy": ["a", "b"]}, index=no_controls.index).to_csv(cfg.taxonomy_file_path, sep="\t")

    _, _, contam_asvs, _ = decontaminate.remove_contaminants(
        cfg.counts_file_path, cfg.taxonomy_file_path, cfg.output_path, clean_count_file=CLEAN_COUNT_FILE,
    )

    assert contam_asvs == []
    clean = pd.read_csv(f"{cfg.output_path}/{CLEAN_COUNT_FILE}", index_col=0, sep="\t")
    pd.testing.assert_frame_equal(clean, no_controls)
