"""Run manifest: what produced a set of pipeline outputs.

Written next to the outputs as run_manifest.json so anything that reads them
later (a database loader, a re-analysis) knows the package version, the
reference training set and the settings behind the counts and taxonomy,
without inferring them from file names.
"""

import hashlib
import json
import os
from datetime import datetime, timezone
from importlib.metadata import version

from .config import RunConfig
from .constants import DECONTAM_PREVALENCE_THRESHOLD, TAXONOMY_CONFIDENCE_THRESHOLD

MANIFEST_FILE = "run_manifest.json"
PACKAGE_NAME = "tag-analysis"
_HASH_CHUNK_BYTES = 1 << 20


def file_sha256(path):
    """SHA-256 of a file, read in chunks (reference training sets are large)."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(_HASH_CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decontamination_record(library_batches, counts_df, contam_asvs, predicted_controls):
    """Manifest section describing one decontam pass.

    counts_df is the raw ASV-by-library counts table; predicted_controls is
    aligned with its columns; library_batches is the batch-table path or None.
    """
    controls = [library for library, is_control in zip(counts_df.columns, predicted_controls) if is_control]
    return {
        "prevalence_threshold": DECONTAM_PREVALENCE_THRESHOLD,
        "library_batches": None if library_batches is None else {
            "file": os.path.basename(library_batches),
            "sha256": file_sha256(library_batches),
        },
        "control_libraries": controls,
        "contaminant_asvs": list(contam_asvs),
    }


def build_manifest(config: RunConfig, dada2_kwargs, counts_df, contam_asvs,
                   predicted_controls, clean_count_file):
    """
    The manifest for one finished run.

    counts_df is the raw ASV-by-library counts table; predicted_controls is
    aligned with its columns.
    """
    return {
        "tag_analysis_version": version(PACKAGE_NAME),
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "dataset_name": config.dataset_name,
        "gene": config.primers.name,
        "primers": {"forward": config.primers.fwd, "reverse": config.primers.rev},
        "reference": {
            "file": os.path.basename(config.reference_db_path),
            "sha256": file_sha256(config.reference_db_path),
        },
        "dada2": dada2_kwargs,
        "taxonomy_confidence_threshold": TAXONOMY_CONFIDENCE_THRESHOLD,
        "decontamination": decontamination_record(
            config.library_batches, counts_df, contam_asvs, predicted_controls
        ),
        "outputs": {
            "counts": os.path.basename(config.counts_file_path),
            "clean_counts": clean_count_file,
            "taxonomy": os.path.basename(config.taxonomy_file_path),
        },
    }


def write_manifest(config: RunConfig, manifest):
    """Write the manifest into the run's output directory; returns its path."""
    path = os.path.join(config.output_path, MANIFEST_FILE)
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2)
    return path


def record_redecontamination(output_dir, record):
    """Replace the decontamination section of an existing run's manifest.

    Used when decontam is rerun on persisted outputs (see redecontam.py): the
    upstream sections (reference, DADA2 settings) still describe how the counts
    were made, so only decontamination is replaced, stamped with when and by
    which package version. A run from before manifests existed gets a manifest
    holding only that section.
    """
    path = os.path.join(output_dir, MANIFEST_FILE)
    manifest = {}
    if os.path.exists(path):
        with open(path) as handle:
            manifest = json.load(handle)
    manifest["decontamination"] = {
        **record,
        "rerun_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "rerun_tag_analysis_version": version(PACKAGE_NAME),
    }
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2)
    return path
