"""Per-bottle community timecourse grid, keyed to saltyBiomass batch data.

Lays out one row per incubation batch and one column per amplicon
(16S on the left, 18S on the right). Each cell is a stacked relative
abundance bar chart across timepoints (T0, T1, T2, ...), so community
shifts within a batch read left-to-right and the same batch's
prokaryotic vs. eukaryotic response reads across a row.

Row order comes from measured water activity in the saltyBiomass
transformed data (`ISQ_DATA_*.ecsv`): highest a_w on top. The same
file supplies gas-production timeseries for community-vs-gas
comparison (`load_gas_timeseries`).

Keys are canonicalized to integers so the three ID spellings in play
("Exp03" in fastq names, "Exp003" in the ECSV Experiment column,
"Exp_03_B01_R01" in Sample IDs) all merge:
  experiment: int (3), batch: int (1)

Amplicon inputs are the persisted DADA2 pipeline outputs:
  - sequence table (samples x sequences, chimera-removed), comma-separated
  - ASV_taxonomy.csv (tab-separated, as written by etl.assign_taxonomy)
"""

import re

import pandas as pd
from astropy.table import Table
from bokeh.io import output_file
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool, Title
from bokeh.plotting import figure, save
import seaborn as sns

from tag_analysis.constants import COLORS

SAMPLE_NAME_RE = re.compile(
    r"GP_\d+_(?P<experiment>Exp_?\d+)_(?P<batch>B\d+)_(?P<timepoint>T\d+)"
    r"(?P<wash>_wash)?_S\d+"
)

# Gas channels treated as production series. CH4 uses the FID channel;
# the MS m/z 15 CH4 channel is diagnostic-only in saltyBiomass
# (constants.PRODUCTION_GASES) and is excluded here for the same reason.
PRODUCTION_GASES = ["CO2", "N2O", "CH4_FID"]


def canonical_experiment(value):
    """'Exp03' / 'Exp003' / 'Exp_03' / 3 -> 3."""
    if isinstance(value, str):
        return int(value.replace("Exp", "").replace("_", ""))
    return int(value)


def canonical_batch(value):
    """'B01' / 'B1' / 1 -> 1."""
    if isinstance(value, str):
        return int(value.lstrip("B"))
    return int(value)


def parse_sample_name(name):
    """Parse experiment/batch/timepoint from a filtered-fastq sample name.

    Returns None for anything that is not an experiment batch sample
    (PCR negatives, extraction negatives, algae brick reads, washes) so
    callers can drop non-timecourse samples with a single filter.
    """
    m = SAMPLE_NAME_RE.search(name)
    if m is None or m.group("wash"):
        return None
    return {
        "experiment": canonical_experiment(m.group("experiment")),
        "batch": canonical_batch(m.group("batch")),
        "timepoint": m.group("timepoint"),
    }


def load_amplicon_run(seqtab_path, taxonomy_path):
    """Load a sequence table + taxonomy into a long relative-abundance frame.

    Relative abundance is per-sample: relabund_i = 100 * n_i / sum_j(n_j),
    where n_i is the read count of ASV i in that sample (dimensionless %).

    Returns long df with columns:
      sample, experiment, batch, timepoint, ASV, relabund, <RANKS...>
    Non-batch samples (negatives, algae bricks, washes) are excluded.
    """
    seqtab = pd.read_csv(seqtab_path, index_col=0)
    taxonomy = pd.read_csv(taxonomy_path, sep="\t")

    seq_to_asv = dict(zip(taxonomy["sequence"], taxonomy["ASV_ID"]))
    counts = seqtab.rename(columns=seq_to_asv)

    relabund = counts.div(counts.sum(axis=1), axis=0) * 100

    long_df = (
        relabund.stack().rename("relabund").rename_axis(["sample", "ASV"]).reset_index()
    )

    parsed = long_df["sample"].map(parse_sample_name)
    long_df = long_df[parsed.notna()].copy()
    meta = pd.DataFrame(list(parsed.dropna()), index=parsed.dropna().index)
    long_df = long_df.join(meta)

    rank_cols = [
        c
        for c in taxonomy.columns
        if c not in ("ASV_ID", "sequence", "taxonomy", "confidence")
    ]
    long_df = long_df.merge(
        taxonomy[["ASV_ID"] + rank_cols], left_on="ASV", right_on="ASV_ID"
    ).drop(columns="ASV_ID")
    return long_df


def load_isq_ecsv(ecsv_path):
    """Read a saltyBiomass ISQ_DATA_*.ecsv into pandas, keeping unit info.

    astropy's to_pandas() drops Quantity units, so they are captured from
    the Table first and returned alongside.

    Returns (df, units) where units maps column name -> astropy unit string
    for every column that carries one, and df gains integer `experiment`
    and `batch` key columns.
    """
    table = Table.read(ecsv_path)
    units = {
        name: str(col.unit)
        for name, col in table.columns.items()
        if getattr(col, "unit", None) is not None
    }
    df = table.to_pandas()
    df["experiment"] = df["Experiment"].map(canonical_experiment)
    df["batch"] = df["Batch ID"].map(canonical_batch)
    return df, units


def load_batch_metadata(ecsv_path):
    """Per-batch metadata (one row per experiment x batch) from the ECSV.

    Water Activity is batch-level in saltyBiomass (broadcast from the
    batch sheet), so duplicates across measurement rows collapse to one.
    A batch with more than one distinct a_w in the file is a data error
    upstream and fails loudly here.
    """
    df, _ = load_isq_ecsv(ecsv_path)
    keep = ["experiment", "batch", "Water Activity"]
    optional = ["Salt Composition", "Salt Makeup", "Contains Sulfate"]
    keep += [c for c in optional if c in df.columns]
    meta = df[keep].drop_duplicates()

    counts = meta.groupby(["experiment", "batch"])["Water Activity"].nunique()
    conflicted = counts[counts > 1]
    if not conflicted.empty:
        raise ValueError(
            f"Multiple Water Activity values per batch: {conflicted.to_dict()}"
        )
    return meta.rename(columns={"Water Activity": "water_activity"})


def load_gas_timeseries(ecsv_path, gases=None):
    """Per-batch gas production timeseries from the ECSV.

    Returns (df, units): rows are individual measurements with columns
      experiment, batch, Replicate ID, Molecule, Days since start,
      Cumulative Moles
    averaging nothing — replicate bottles within a batch stay separate
    so the caller decides how to aggregate against pooled DNA samples.
    """
    df, units = load_isq_ecsv(ecsv_path)
    gases = PRODUCTION_GASES if gases is None else gases
    cols = [
        "experiment",
        "batch",
        "Replicate ID",
        "Molecule",
        "Days since start",
        "Cumulative Moles",
    ]
    gas_df = df[df["Molecule"].isin(gases)][cols].copy()
    return gas_df, units


def batch_order_by_water_activity(batch_metadata):
    """Order (experiment, batch) keys by descending water activity.

    batch_metadata: df from load_batch_metadata. water_activity is the
    meter-read a_w (dimensionless, 0-1).
    """
    ordered = batch_metadata.sort_values("water_activity", ascending=False)
    return list(zip(ordered["experiment"], ordered["batch"]))


def _aggregate_for_batch(long_df, experiment, batch, taxonomic_level, min_abundance):
    """Sum relabund to taxonomic_level per timepoint; lump minors to Other."""
    sub = long_df[(long_df["experiment"] == experiment) & (long_df["batch"] == batch)]
    agg = sub.groupby(["timepoint", taxonomic_level])["relabund"].sum().reset_index()
    agg[taxonomic_level] = agg[taxonomic_level].fillna("").replace("", "unclassified")
    peak = agg.groupby(taxonomic_level)["relabund"].max()
    minor = peak[peak < min_abundance].index
    agg.loc[agg[taxonomic_level].isin(minor), taxonomic_level] = "Other"
    return agg.groupby(["timepoint", taxonomic_level])["relabund"].sum().reset_index()


def _taxon_color_map(long_df, taxonomic_level, min_abundance):
    """One consistent taxon->color map per amplicon across all batches."""
    level = long_df[taxonomic_level].fillna("").replace("", "unclassified")
    per_sample = (
        long_df.assign(**{taxonomic_level: level})
        .groupby(["sample", taxonomic_level])["relabund"]
        .sum()
    )
    peak = per_sample.groupby(taxonomic_level).max()
    taxa = sorted(peak[peak >= min_abundance].index)
    palette = (
        COLORS[: len(taxa)]
        if len(taxa) <= len(COLORS)
        else sns.color_palette("husl", len(taxa)).as_hex()
    )
    color_map = dict(zip(taxa, palette))
    color_map["Other"] = "#808080"
    return color_map


def _batch_panel(
    agg, timepoints, taxonomic_level, color_map, title, width, height, show_x_labels
):
    """One stacked-bar cell: x = timepoints, stacks = taxa."""
    taxa = [t for t in color_map if t in set(agg[taxonomic_level])]
    wide = (
        agg.pivot_table(
            index="timepoint",
            columns=taxonomic_level,
            values="relabund",
            fill_value=0,
        )
        .reindex(timepoints, fill_value=0)
        .reindex(columns=taxa, fill_value=0)
        .reset_index()
    )
    source = ColumnDataSource(wide)

    p = figure(
        x_range=timepoints,
        width=width,
        height=height,
        toolbar_location=None,
        tools="",
    )
    p.add_tools(
        HoverTool(
            tooltips=[
                ("Timepoint", "@timepoint"),
                (taxonomic_level.title(), "$name"),
                ("Abundance", "@$name{0.0}%"),
            ]
        )
    )
    p.vbar_stack(
        taxa,
        x="timepoint",
        width=0.8,
        source=source,
        color=[color_map[t] for t in taxa],
        alpha=0.85,
    )
    p.y_range.start = 0
    p.y_range.end = 100
    p.add_layout(Title(text=title, text_font_size="9pt"), "above")
    p.xaxis.visible = show_x_labels
    p.yaxis.axis_label = "%"
    p.yaxis.axis_label_text_font_size = "7pt"
    p.yaxis.major_label_text_font_size = "7pt"
    return p


def _legend_panel(color_map, taxonomic_level, height):
    """Standalone legend column for one amplicon's taxon colors."""
    p = figure(
        width=260,
        height=height,
        toolbar_location=None,
        tools="",
        x_range=(0, 1),
        y_range=(0, 1),
    )
    p.axis.visible = False
    p.grid.visible = False
    p.outline_line_color = None
    for taxon, color in color_map.items():
        p.scatter(x=[-1], y=[-1], color=color, size=10, legend_label=str(taxon)[:40])
    p.legend.location = "top_left"
    p.legend.label_text_font_size = "8pt"
    p.legend.spacing = 0
    p.add_layout(
        Title(text=f"{taxonomic_level.title()} legend", text_font_size="9pt"),
        "above",
    )
    return p


def create_batch_timecourse_grid(
    runs,
    batch_order,
    output_path,
    batch_metadata=None,
    taxonomic_level="phylum",
    min_abundance=2.0,
    panel_width=240,
    panel_height=170,
    title=None,
):
    """Grid of stacked-bar timecourses: rows = batches, columns = amplicons.

    runs: dict of column label -> long df from load_amplicon_run, in
        left-to-right column order (e.g. {"16S": df16, "18S": df18}).
    batch_order: list of (experiment, batch) int tuples, top row first --
        typically batch_order_by_water_activity(load_batch_metadata(ecsv)).
    batch_metadata: optional df from load_batch_metadata; when given, row
        titles carry the batch's a_w.
    min_abundance: peak per-sample relabund (%) below which a taxon is
        lumped into "Other".

    Batches absent from a run render as an empty cell in that column.
    """
    color_maps = {
        label: _taxon_color_map(df, taxonomic_level, min_abundance)
        for label, df in runs.items()
    }
    timepoints = sorted(
        {tp for df in runs.values() for tp in df["timepoint"].unique()},
        key=lambda t: int(t[1:]),
    )

    aw_lookup = {}
    if batch_metadata is not None:
        aw_lookup = {
            (row.experiment, row.batch): row.water_activity
            for row in batch_metadata.itertuples()
        }

    grid = []
    for i, (experiment, batch) in enumerate(batch_order):
        last_row = i == len(batch_order) - 1
        aw = aw_lookup.get((experiment, batch))
        aw_txt = f"  a_w={aw:.3f}" if aw is not None else ""
        row = []
        for label, df in runs.items():
            agg = _aggregate_for_batch(
                df, experiment, batch, taxonomic_level, min_abundance
            )
            if agg.empty:
                row.append(None)
                continue
            row.append(
                _batch_panel(
                    agg,
                    timepoints,
                    taxonomic_level,
                    color_maps[label],
                    title=f"Exp{experiment:03d} B{batch:02d} — {label}{aw_txt}",
                    width=panel_width,
                    height=panel_height,
                    show_x_labels=last_row,
                )
            )
        grid.append(row)

    legend_row = [
        _legend_panel(
            color_maps[label], taxonomic_level, height=max(300, panel_height * 2)
        )
        for label in runs
    ]
    grid.append(legend_row)

    output_file(output_path, title=title or f"Batch timecourses ({taxonomic_level})")
    save(gridplot(grid, merge_tools=True))
    return output_path
