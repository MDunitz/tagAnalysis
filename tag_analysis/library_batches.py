"""Map sequenced libraries to their DNA extraction batch.

Each extraction batch carries its own kit negative, so contamination is
batch-specific: a library should be judged against the negatives extracted
alongside it, not against every negative in the run.

The package does not know how any project records its extractions. A project
supplies a library -> batch table (CSV, columns LIBRARY_BATCHES_LIBRARY_COL and
LIBRARY_BATCHES_BATCH_COL), built however that project's records require; this
module only reads it and places the counts-table columns.
"""

import pandas as pd

from .constants import (
    LIBRARY_BATCHES_BATCH_COL,
    LIBRARY_BATCHES_LIBRARY_COL,
    LIBRARY_ID_SUFFIX_RE,
)


def library_id(fastq_name):
    """Library ID of a forward-read fastq name (raw or filtered); see LIBRARY_ID_SUFFIX_RE."""
    return LIBRARY_ID_SUFFIX_RE.sub("", fastq_name)


def read_library_batches(path):
    """{library ID: extraction batch or None} from a library-batch CSV.

    A blank batch is a shared control (never extracted, e.g. a PCR negative).
    A library listed twice is a table error and raises.
    """
    table = pd.read_csv(path, dtype={LIBRARY_BATCHES_LIBRARY_COL: str})
    duplicated = table[LIBRARY_BATCHES_LIBRARY_COL][table[LIBRARY_BATCHES_LIBRARY_COL].duplicated()]
    if not duplicated.empty:
        raise ValueError(f"Libraries listed more than once in {path}: {sorted(duplicated)}")
    return {
        lib: (None if pd.isna(batch) else batch)
        for lib, batch in zip(table[LIBRARY_BATCHES_LIBRARY_COL], table[LIBRARY_BATCHES_BATCH_COL])
    }


def assign_extraction_batches(libraries, is_control, batch_by_library):
    """Extraction batch for each counts-table column, aligned with `libraries`.

    Columns are matched to the table by library ID. A control absent from the
    table, or listed with a blank batch, is a shared control (None). A
    non-control library with no batch cannot be placed and raises, listing
    every such library.
    """
    batches, unplaced = [], []
    for library, control in zip(libraries, is_control):
        batch = batch_by_library.get(library_id(library))
        if batch is None and not control:
            unplaced.append(library)
        batches.append(batch)
    if unplaced:
        raise ValueError(f"No extraction batch in the library-batch table for: {unplaced}")
    return batches
