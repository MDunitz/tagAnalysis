# analysis/

Project-side analysis scripts that consume persisted pipeline outputs
(sequence tables, ASV_taxonomy.csv). Package code that runs *inside* the
pipeline stays in `tag_analysis/`; anything here reads its outputs.

## bottle_timecourse.py

Per-bottle community timecourse grid: rows = bottles (ordered by
descending water activity via `bottle_order_by_water_activity`),
columns = amplicons (16S left, 18S right), each cell a stacked relative
abundance bar chart across timepoints.

```python
from analysis.bottle_timecourse import (
    load_amplicon_run,
    bottle_order_by_water_activity,
    create_bottle_timecourse_grid,
)
import pandas as pd

df16 = load_amplicon_run("out16/sequence_table_nochim.csv", "out16/ASV_taxonomy.csv")
df18 = load_amplicon_run("out18/sequence_table_nochim.csv", "out18/ASV_taxonomy.csv")

meta = pd.read_csv("analysis/bottle_metadata.csv")  # fill from template
order = bottle_order_by_water_activity(meta)

create_bottle_timecourse_grid(
    {"16S": df16, "18S": df18},
    order,
    "outputs/bottle_timecourse_phylum.html",
    taxonomic_level="phylum",
)
```

`bottle_metadata_template.csv` lists the bottles present in the current
sequencing run; copy to `bottle_metadata.csv` and fill `water_activity`
with the WP4C meter-read a_w per bottle (dimensionless). Data files are
untracked (too large) — point the loaders at local copies.
