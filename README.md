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

## Batch-aware decontamination

By default every negative control in the run (libraries matching
`identify_control_samples`) is pooled against every sample. When DNA was
extracted in several batches, each with its own kit negative, pass the
extraction sheet so each library is judged against its own batch's negatives:

    cfg16 = RunConfig(..., extraction_sheet="https://docs.google.com/spreadsheets/d/<id>/edit")

`extraction_sheet` is a CSV path or a Google Sheets link (shared readable by
link; it is fetched as CSV). The sheet needs a `TAG ID` column and an
`Extraction Batch ID` column. `TAG ID` holds the sequencing tag number(s) of
the libraries made from that extraction, i.e. the `<n>` in `GP_<n>_...` library
names. List several numbers in one cell when one extraction went into more
than one library (e.g. `3, 52` for its 16S and 18S libraries).

- Prevalence decontam runs once per batch, on that batch's libraries plus any
  control with no extraction row (e.g. a PCR negative, which applies to every
  batch).
- An ASV's p-value is its minimum across batches (decontam's
  `batch.combine="minimum"`). Per-batch p-values are kept as
  `p_prev_batch_<b>` columns in `contam_results.csv`.
- A non-control library with no `TAG ID` row raises, naming the library.
  Libraries are never silently pooled.

## Run manifest

Every run writes `run_manifest.json` next to its outputs, recording what
produced them: the package version, gene and primers, the reference training
set (file name and SHA-256), the DADA2 settings, the taxonomy confidence
threshold, the control libraries and contaminant ASVs from decontamination,
and the names of the counts, cleaned counts and taxonomy files. Anything that
reads the outputs later (for example the orphanDB loader) takes the reference
version from here instead of guessing it from a path.

`ASVs_counts_clean.csv` is always written; when no contaminant is found it
equals `ASVs_counts.csv`.

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
