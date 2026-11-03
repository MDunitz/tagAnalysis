# Changelog

## 0.5.0

Batch-aware decontamination, decontam reruns on finished runs, and PR2 rank
names for 18S.

- Batch-aware decontam (#21, #22). `RunConfig.library_batches` takes a CSV of
  `library,extraction_batch`. Prevalence decontam then runs within each DNA
  extraction batch against that batch's own kit negatives. PCR negatives and
  other blank-batch controls are shared across batches, and per-batch p-values
  combine by minimum. Without a table, the pooled behavior is unchanged.
  `run_manifest.json` records the table's SHA-256.
- `redecontam` (#23): `make redecontam INPUT_DIR=... LIBRARY_BATCHES=...`, or
  `python -m tag_analysis.redecontam`, reruns decontam and its plots on a
  finished run's `ASVs_counts.csv` and `ASV_taxonomy.csv`, without primer
  removal, DADA2 or taxonomy. It replaces only the manifest's decontamination
  section.
- PR2 rank names for 18S (#20). 18S runs against PR2 split taxonomy into PR2's
  nine levels (domain, supergroup, division, subdivision, class, order,
  family, genus, species) instead of SILVA's seven names. Override with
  `RunConfig.ranks`. 18S taxonomy tables written by 0.5.0 have different
  column names than 0.4.0 ones: there is no `phylum` column.

## 0.4.0

- Run manifest (`run_manifest.json`) written next to every run's outputs (#19).
