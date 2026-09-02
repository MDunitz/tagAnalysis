# analysis/

Project-side analysis scripts that consume persisted pipeline outputs
(sequence tables, ASV_taxonomy.csv). Package code that runs *inside* the
pipeline stays in `tag_analysis/`; anything here reads its outputs.

## bottle_timecourse.py

Per-batch community timecourse grid: rows = incubation batches ordered by
descending water activity, columns = amplicons (16S left, 18S right),
each cell a stacked relative abundance bar chart across timepoints.

Batch metadata (a_w, salt composition) and gas production timeseries both
come from the saltyBiomass transformed data (`ISQ_DATA_*.ecsv`), so
nothing is hand-entered. Keys are canonicalized to ints: `Exp03` (fastq
names), `Exp003` (ECSV Experiment), and `Exp_03` (Sample IDs) all merge.

```python
from analysis.bottle_timecourse import (
    load_amplicon_run,
    load_batch_metadata,
    load_gas_timeseries,
    batch_order_by_water_activity,
    create_batch_timecourse_grid,
)

df16 = load_amplicon_run("out16/16s_sequence_table_nochim.csv", "out16/16s_ASV_taxonomy.csv")
df18 = load_amplicon_run("out18/18s_sequence_table_nochim.csv", "out18/18s_ASV_taxonomy.csv")

meta = load_batch_metadata("path/to/ISQ_DATA_<latest>.ecsv")
order = batch_order_by_water_activity(meta)

create_batch_timecourse_grid(
    {"16S": df16, "18S": df18},
    order,
    "outputs/batch_timecourse_phylum.html",
    batch_metadata=meta,
    taxonomic_level="phylum",
)

gas, units = load_gas_timeseries("path/to/ISQ_DATA_<latest>.ecsv")
```

CH4 production uses the FID channel (`CH4_FID`); the MS m/z 15 channel is
diagnostic-only in saltyBiomass and is excluded here for the same reason.
Data files are untracked (too large) -- point the loaders at local copies.
