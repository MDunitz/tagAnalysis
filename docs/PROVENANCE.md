# Provenance and refactor notes

`tag_analysis` was extracted and refactored from the private `MDunitz/algaeBricks`
repo, from `experiments/processing/genome_analysis/v2`.

## What changed in the port

- Intra-repo absolute imports (`experiments.processing.genome_analysis.v2.*`)
  rewritten as a proper package (`tag_analysis`) with a unique top-level import
  name, so installing it into another repo never collides with that repo's own
  modules (the same convention as `labdata` from DataExtractionHelpers).
- Per-run configuration (data path, output path, dataset name, reference DB)
  moved out of mutated module globals into a `RunConfig` dataclass. No project
  edits package source to configure a run.
- The taxonomy reference DB is an explicit per-run parameter
  (`reference_db_path`) instead of a single hardcoded SILVA path, so 18S runs
  can point at PR2 rather than being forced onto SILVA.
- `process_16s` / `process_18s` are thin wrappers over a shared `_run()`; the
  two amplicons differ only in primer set (verified identical downstream).
- Primer reverse-complement constants are unit-tested against computed RCs.
- A dead, un-importable scratch entry script was dropped; the full pipeline
  sequence was reconstructed from its (complete) logic into `pipelines.py`.

## Open verification items

- The R / DADA2 / cutadapt / reference-DB stages are mock-tested only.
  Acceptance test before trusting results: re-run a known dataset and diff the
  ASV counts table against the original v2 output.
- The carried-over DADA2 `truncLen=(230, 200)` default was tuned for one 16S
  run; it should be re-checked per amplicon against read-quality profiles
  (the 18S V4 amplicon length differs from 16S V4).
