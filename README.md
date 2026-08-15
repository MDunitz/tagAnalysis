# tag_analysis

Importable 16S/18S amplicon (DADA2) processing pipeline for Orphan lab tag
data. A standalone, installable package usable as a dependency across projects.

## Install

    pip install git+https://github.com/MDunitz/tagAnalysis.git

Or set up a dev environment (editable install + test tooling):

    pip install -r requirements.txt

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

DADA2 parameters can be overridden per run. The default `truncLen=(230, 200)`
was tuned for a specific 16S run and should be re-checked against your own
read-quality profiles per amplicon:

    process_16s(cfg16, dada2_kwargs={"truncLen": (240, 180), "maxEE": (2, 2)})

## Testing

    pytest

The suite covers the pure-Python transforms (ASV table construction,
taxonomy rank splitting, relative-abundance melting, contamination
aggregation, file/sample path generation) and mock-boundary tests for the
R/subprocess calls and pipeline wiring. It does **not** exercise R, DADA2,
cutadapt, or a real reference DB — scientific correctness of those stages
requires a real run against fastq data.

## Provenance

See `docs/PROVENANCE.md` for origin and refactor history.
