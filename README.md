# tag_analysis

Importable 16S/18S amplicon (DADA2) processing pipeline for Orphan lab tag
data. Extracted and refactored from `MDunitz/algaeBricks`
(`experiments/processing/genome_analysis/v2`) into a standalone, installable
package so it can be used as a dependency across projects.

## What changed from the algaeBricks version

- Intra-repo absolute imports (`experiments.processing.genome_analysis.v2.*`)
  rewritten as a proper package (`tag_analysis`) with a unique top-level
  import name (parallel to `labdata`), so installing it into another repo
  never collides with that repo's own modules.
- Per-run configuration (data path, output path, dataset name, reference DB)
  moved out of mutated module globals into a `RunConfig` dataclass. No project
  needs to edit package source to configure a run.
- The taxonomy reference DB is now an explicit per-run parameter
  (`reference_db_path`) instead of a single hardcoded `PATH_TO_SILVA_DB`, so
  18S runs can point at PR2 rather than being forced onto SILVA.
- `process_16s` / `process_18s` are thin wrappers over a shared core; the two
  amplicons differ only in primer set (verified identical downstream).
- Primer reverse-complement constants are unit-tested against computed RCs.
- Dead/broken `process_reads.py` scratch file dropped; the full pipeline
  sequence was reconstructed from its (complete) logic into `pipelines.py`.

## Install

    pip install git+https://github.com/MDunitz/tagAnalysis.git

Or editable, for development:

    pip install -e ".[test]"

## External tools (NOT pip-installable)

This package orchestrates external programs that must be present and
configured separately:

- **R** with packages `dada2`, `DECIPHER`, `decontam` (`Rscript` on `PATH`)
- **cutadapt**
- A **taxonomy training set**: SILVA `.RData` for 16S, PR2 for 18S

## Usage

    from tag_analysis import RunConfig, process_16s, process_18s

    cfg16 = RunConfig(
        data_path="data/local_data/<project>/16S",
        output_path="data/local_data/<project>/output/16s",
        dataset_name="<project>",
        reference_db_path="data/local_data/SILVA_SSU_r138_2_2024.RData",
    )
    process_16s(cfg16)

    cfg18 = RunConfig(
        data_path="data/local_data/<project>/18S",
        output_path="data/local_data/<project>/output/18s",
        dataset_name="<project>",
        reference_db_path="data/local_data/PR2_version_5.RData",  # PR2 for 18S
    )
    process_18s(cfg18)

DADA2 parameters can be overridden per run. Note the original defaults
(`truncLen=(230,200)`, etc.) were tuned for a specific 16S run and should be
re-checked against your own read quality profiles per amplicon:

    process_16s(cfg16, dada2_kwargs={"truncLen": (240, 180), "maxEE": (2, 2)})

## Testing

    pytest

The suite covers the pure-Python transforms (ASV table construction,
taxonomy rank splitting, relative-abundance melting, contamination
aggregation, file/sample path generation) and mock-boundary tests for the
R/subprocess calls and pipeline wiring. It does **not** exercise R, DADA2,
cutadapt, or a real reference DB — scientific correctness of those stages
requires a real run against fastq data.
