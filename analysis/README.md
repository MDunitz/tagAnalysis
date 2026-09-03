# analysis/

Project-side analysis scripts that consume persisted pipeline outputs
(sequence tables, ASV_taxonomy.csv) and the saltyBiomass incubation data
(`ISQ_DATA_*.ecsv`). Package code that runs *inside* the amplicon pipeline
stays in `tag_analysis/`; anything here reads its outputs.

## amplicon_data.py — loaders and key joins

Three sources meet on canonical integer `(experiment, batch)` keys, so the
`Exp03` / `Exp003` / `Exp_03_B01_R04` spellings all merge:

- `load_amplicon_run(seqtab, taxonomy)` — long relative-abundance frame,
  negatives / algae bricks / `_wash` samples dropped by the name parser.
- `load_batch_metadata(ecsvs)` — water activity, salt composition per batch.
- `load_gas_timeseries(ecsvs)` — cumulative moles vs days, per replicate.
- `load_sampling_events(ecsvs)` — the pressure sheet's `Sampled` flag marks
  a **destructive** sample (bottle opened, sampled, N2-flushed). These rows
  carry a Replicate ID and are the only record tying a sequencing timepoint
  to a specific bottle.
- `map_timepoints_to_days(events)` — T0 at day 0 (setup), T1..Tn onto the
  batch's sampling events in date order.
- `unmapped_timepoints(amplicon_df, tmap)` — sequencing timepoints with no
  sampling record; a non-empty result means those samples cannot be placed
  on the incubation day axis.

The `.ecsv` files are per-experiment, so every loader takes a path or a
list of paths.

## bottle_timecourse.py — the grid

Rows = batches ordered by descending a_w, columns = amplicons (16S left,
18S right), each cell a stacked relative-abundance bar chart across
timepoints. Passing `gas_df` adds a rightmost gas-production column:
cumulative moles vs days, one line per gas, the DNA-source bottle solid
and the other replicates faded, with dashed vertical lines at the
destructive-sampling days.

```python
import glob
from analysis.amplicon_data import (
    load_amplicon_run, load_batch_metadata, load_gas_timeseries,
    load_sampling_events, map_timepoints_to_days, unmapped_timepoints,
)
from analysis.bottle_timecourse import (
    batch_order_by_water_activity, create_batch_timecourse_grid,
)

ecsvs = sorted(glob.glob("data/transformed/ISQ/ISQ_DATA_*_<commit>.ecsv"))

df16 = load_amplicon_run("out16/16s_sequence_table_nochim.csv", "out16/16s_ASV_taxonomy.csv")
df18 = load_amplicon_run("out18/18s_sequence_table_nochim.csv", "out18/18s_ASV_taxonomy.csv")

meta = load_batch_metadata(ecsvs)
gas, units = load_gas_timeseries(ecsvs)
tmap = map_timepoints_to_days(load_sampling_events(ecsvs))
print(unmapped_timepoints(df18, tmap))  # check before trusting the day axis

create_batch_timecourse_grid(
    {"16S": df16, "18S": df18},
    batch_order_by_water_activity(meta),
    "outputs/batch_timecourse_phylum.html",
    batch_metadata=meta,
    gas_df=gas,
    timepoint_map=tmap,
    moles_unit=units["Cumulative Moles"],
    taxonomic_level="phylum",
)
```

CH4 production uses the FID channel (`CH4_FID`); the MS m/z 15 channel is
diagnostic-only in saltyBiomass and is excluded for the same reason.
Sequencing and gas data files are untracked — point the loaders at local
copies.
