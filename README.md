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
extracted in several batches, each with its own kit negative, pass a
library -> extraction batch table so each library is judged against its own
batch's negatives:

    cfg16 = RunConfig(..., library_batches="path/to/library_batches.csv")

The table is a CSV with two columns:

| column | meaning |
|---|---|
| `library` | library ID: the forward-read fastq name without `_R1_001.fastq.gz` (e.g. `GP_3_Exp03_B01_T0_S99`) |
| `extraction_batch` | the DNA extraction batch; blank for a shared control (e.g. a PCR negative) |

The package does not read any project's extraction records. Each project builds
this table from its own records, so a change in how a lab keeps its sheet is a
change to that project's code, not to this package.

- Prevalence decontam runs once per batch, on that batch's libraries plus the
  shared controls (blank batch, or a control absent from the table).
- An ASV's p-value is its minimum across batches (decontam's
  `batch.combine="minimum"`). Per-batch p-values are kept as
  `p_prev_batch_<b>` columns in `contam_results.csv`.
- A non-control library missing from the table raises, naming it; a library
  listed twice raises.
- `run_manifest.json` records the table's file name and SHA-256, so a run can
  be traced to the exact table it used.

### Applying it to a finished run

Decontam reads only `ASVs_counts.csv` and `ASV_taxonomy.csv`, so a finished
run can be re-decontaminated without primer removal, DADA2 or taxonomy:

    make redecontam INPUT_DIR=path/to/output LIBRARY_BATCHES=path/to/library_batches.csv
    # or: python -m tag_analysis.redecontam path/to/output --library-batches path/to/library_batches.csv

This rewrites `ASVs_counts_clean.csv`, `contam_results.csv`, the
contamination plot and the stackbars. It replaces only the
`decontamination` section of `run_manifest.json`, stamped with `rerun_at`
and `rerun_tag_analysis_version`; the upstream sections still describe how the
counts were made. Needs R with `decontam`.

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
