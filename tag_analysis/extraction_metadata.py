"""Map sequenced libraries to their DNA extraction batch.

Each extraction batch carries its own kit negative, so contamination is
batch-specific: a library should be judged against the negatives extracted
alongside it, not against every negative in the sequencing run. The mapping
comes from the lab's extraction sheet (one row per extraction), joined to the
counts-table columns by the sequencing-facility tag number.

The sheet is a CSV path or a Google Sheets link (any share/edit URL; it is
rewritten to the sheet's CSV export, so the sheet must be readable by link).
"""

import re

import pandas as pd

from .constants import (
    EXTRACTION_BATCH_COL,
    EXTRACTION_TAG_ID_COL,
    LIBRARY_TAG_ID_RE,
)

# Google Sheets share/edit URL: .../spreadsheets/d/<sheet id>/... with an
# optional gid=<tab id> anywhere after it (query string or #fragment).
_GOOGLE_SHEET_RE = re.compile(
    r"docs\.google\.com/spreadsheets/d/(?P<sheet_id>[A-Za-z0-9_-]+)"
    r"(?:.*?[?#&]gid=(?P<gid>\d+))?"
)

# A TAG ID cell lists one or more tag numbers -- one per amplicon library made
# from that extraction (e.g. "3, 52" for its 16S and 18S libraries) --
# separated by commas, semicolons and/or whitespace.
_TAG_ID_SPLIT_RE = re.compile(r"[,;\s]+")


def google_sheet_csv_url(url):
    """CSV export URL for a Google Sheets link; first tab when no gid is given."""
    m = _GOOGLE_SHEET_RE.search(url)
    if m is None:
        raise ValueError(f"Not a Google Sheets URL: {url}")
    gid = m.group("gid") or "0"
    return (
        f"https://docs.google.com/spreadsheets/d/{m.group('sheet_id')}"
        f"/export?format=csv&gid={gid}"
    )


def read_extraction_sheet(source):
    """Extraction sheet as a DataFrame from a CSV path or a Google Sheets URL."""
    if _GOOGLE_SHEET_RE.search(source):
        source = google_sheet_csv_url(source)
    return pd.read_csv(source)


def parse_tag_ids(cell):
    """Tag numbers listed in one TAG ID cell (see _TAG_ID_SPLIT_RE).

    A column holding only single numbers and blanks reads back as float, so a
    numeric cell is one tag number.
    """
    if isinstance(cell, (int, float)):
        return [int(cell)]
    return [int(t) for t in _TAG_ID_SPLIT_RE.split(str(cell).strip()) if t]


def tag_id_to_batch(
    sheet, tag_col=EXTRACTION_TAG_ID_COL, batch_col=EXTRACTION_BATCH_COL
):
    """{tag number: extraction batch} from every sheet row with a TAG ID.

    A tag number listed on rows of two different batches is a data-entry error
    and raises rather than silently picking one.
    """
    mapping = {}
    for cell, batch in sheet.loc[
        sheet[tag_col].notna(), [tag_col, batch_col]
    ].itertuples(index=False):
        for tag in parse_tag_ids(cell):
            if tag in mapping and mapping[tag] != batch:
                raise ValueError(
                    f"TAG ID {tag} is listed in extraction batches {mapping[tag]} and {batch}"
                )
            mapping[tag] = batch
    return mapping


def library_tag_id(library, tag_re=LIBRARY_TAG_ID_RE):
    """Tag number in a library (counts-table column) name, or None if absent."""
    m = tag_re.search(library)
    return None if m is None else int(m.group("tag_id"))


def assign_extraction_batches(
    libraries, is_control, tag_to_batch, tag_re=LIBRARY_TAG_ID_RE
):
    """Extraction batch for each library, aligned with `libraries`.

    A control with no extraction row (e.g. a PCR negative) is never extracted,
    so it belongs to no batch and is returned as None: it is a shared control,
    applied to every batch. A non-control library with no extraction row
    cannot be placed and raises, listing every such library.
    """
    batches, unplaced = [], []
    for library, control in zip(libraries, is_control):
        batch = tag_to_batch.get(library_tag_id(library, tag_re))
        if batch is None and not control:
            unplaced.append(library)
        batches.append(batch)
    if unplaced:
        raise ValueError(
            "No extraction batch for these libraries (fill their TAG ID in the "
            f"extraction sheet): {unplaced}"
        )
    return batches
