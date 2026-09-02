"""Per-bottle community timecourse grid.

Lays out one row per incubation bottle and one column per amplicon
(16S on the left, 18S on the right). Each cell is a stacked relative
abundance bar chart across timepoints (T0, T1, T2, ...), so community
shifts within a bottle read left-to-right and the same bottle's
prokaryotic vs. eukaryotic response reads across a row.

Row order is caller-controlled; `bottle_order_by_water_activity` derives
it from a bottle metadata table (highest a_w on top).

Inputs are the persisted DADA2 pipeline outputs:
  - sequence table (samples x sequences, chimera-removed), comma-separated
  - ASV_taxonomy.csv (tab-separated, as written by etl.assign_taxonomy)

Sample names are the filtered fastq filenames, e.g.
  GP_52_Exp03_B01_T0_S147_R1_filtered.fastq.gz
"""

import re

import pandas as pd
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool, Title
from bokeh.plotting import figure, save
from bokeh.io import output_file
import seaborn as sns

from tag_analysis.constants import COLORS

SAMPLE_NAME_RE = re.compile(
    r"GP_\d+_(?P<experiment>Exp\d+)_(?P<bottle>B\d+)_(?P<timepoint>T\d+)"
    r"(?P<wash>_wash)?_S\d+"
)


def parse_sample_name(name):
    """Parse experiment/bottle/timepoint from a filtered-fastq sample name.

    Returns None for anything that is not an experiment bottle sample
    (PCR negatives, extraction negatives, algae brick reads, washes) so
    callers can drop non-timecourse samples with a single filter.
    """
    m = SAMPLE_NAME_RE.search(name)
    if m is None or m.group("wash"):
        return None
    return {
        "experiment": m.group("experiment"),
        "bottle": m.group("bottle"),
        "timepoint": m.group("timepoint"),
    }


def load_amplicon_run(seqtab_path, taxonomy_path):
    """Load a sequence table + taxonomy into a long relative-abundance frame.

    Relative abundance is per-sample: relabund_i = 100 * n_i / sum_j(n_j),
    where n_i is the read count of ASV i in that sample (dimensionless %).

    Returns long df with columns:
      sample, experiment, bottle, timepoint, ASV, relabund, <RANKS...>
    Non-bottle samples (negatives, algae bricks, washes) are excluded.
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


def bottle_order_by_water_activity(bottle_metadata):
    """Order (experiment, bottle) keys by descending water activity.

    bottle_metadata: df with columns experiment, bottle, water_activity.
    water_activity is the WP4C meter-read a_w (dimensionless, 0-1).
    Every bottle to be plotted must appear here; missing bottles fail
    at lookup time in the grid builder, not silently at the bottom.
    """
    ordered = bottle_metadata.sort_values("water_activity", ascending=False)
    return list(zip(ordered["experiment"], ordered["bottle"]))


def _aggregate_for_bottle(long_df, experiment, bottle, taxonomic_level, min_abundance):
    """Sum relabund to taxonomic_level per timepoint; lump minors to Other."""
    sub = long_df[(long_df["experiment"] == experiment) & (long_df["bottle"] == bottle)]
    agg = sub.groupby(["timepoint", taxonomic_level])["relabund"].sum().reset_index()
    agg[taxonomic_level] = agg[taxonomic_level].fillna("").replace("", "unclassified")
    peak = agg.groupby(taxonomic_level)["relabund"].max()
    minor = peak[peak < min_abundance].index
    agg.loc[agg[taxonomic_level].isin(minor), taxonomic_level] = "Other"
    return agg.groupby(["timepoint", taxonomic_level])["relabund"].sum().reset_index()


def _taxon_color_map(long_df, taxonomic_level, min_abundance):
    """One consistent taxon->color map per amplicon across all bottles."""
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


def _bottle_panel(
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


def create_bottle_timecourse_grid(
    runs,
    bottle_order,
    output_path,
    taxonomic_level="phylum",
    min_abundance=2.0,
    panel_width=240,
    panel_height=170,
    title=None,
):
    """Grid of stacked-bar timecourses: rows = bottles, columns = amplicons.

    runs: dict of column label -> long df from load_amplicon_run, in
        left-to-right column order (e.g. {"16S": df16, "18S": df18}).
    bottle_order: list of (experiment, bottle) tuples, top row first --
        typically bottle_order_by_water_activity(metadata).
    min_abundance: peak per-sample relabund (%) below which a taxon is
        lumped into "Other".

    Bottles absent from a run render as an empty cell in that column
    (e.g. a bottle sequenced for 16S but not 18S).
    """
    color_maps = {
        label: _taxon_color_map(df, taxonomic_level, min_abundance)
        for label, df in runs.items()
    }
    timepoints = sorted(
        {tp for df in runs.values() for tp in df["timepoint"].unique()},
        key=lambda t: int(t[1:]),
    )

    grid = []
    for i, (experiment, bottle) in enumerate(bottle_order):
        last_row = i == len(bottle_order) - 1
        row = []
        for label, df in runs.items():
            agg = _aggregate_for_bottle(
                df, experiment, bottle, taxonomic_level, min_abundance
            )
            if agg.empty:
                row.append(None)
                continue
            row.append(
                _bottle_panel(
                    agg,
                    timepoints,
                    taxonomic_level,
                    color_maps[label],
                    title=f"{experiment} {bottle} — {label}",
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

    output_file(output_path, title=title or f"Bottle timecourses ({taxonomic_level})")
    save(gridplot(grid, merge_tools=True))
    return output_path
