"""Per-batch community timecourse grid, optionally beside gas production.

Rows are incubation batches, ordered by descending water activity
(`batch_order_by_water_activity`); columns are amplicons, 16S left and 18S
right, each cell a stacked relative-abundance bar chart across sequencing
timepoints. Passing a gas frame adds a rightmost column of cumulative-moles
curves with the DNA sampling days marked, so a community shift and the
bottle's gas record line up on one row.

Loaders live in `analysis.amplicon_data`.
"""

import pandas as pd
from bokeh.io import output_file
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool, Span, Title
from bokeh.plotting import figure, save
import seaborn as sns

from tag_analysis.constants import COLORS

# Matches saltyBiomass incubations.constants.GAS_COLORS so a curve keeps its
# color between this grid and the saltyBiomass dashboards.
GAS_COLORS = {
    "CO2": "#1f77b4",
    "N2O": "#ff7f0e",
    "CH4": "#2ca02c",
    "CH4_FID": "#8c564b",
}


def batch_order_by_water_activity(batch_metadata, within_experiment=True):
    """Order (experiment, batch) keys by descending water activity.

    batch_metadata: df from amplicon_data.load_batch_metadata.
    water_activity is the meter-read a_w (dimensionless, 0-1).

    The matched-a_w contrasts (sulfate vs. none, Na vs. Mg) are designed
    *within* an experiment, so `within_experiment=True` blocks rows by
    experiment first and keeps each designed pair on adjacent rows. A
    global a_w sort interleaves experiments and splits pairs apart.
    """
    keys = ["experiment", "water_activity"] if within_experiment else ["water_activity"]
    ordered = batch_metadata.sort_values(
        keys, ascending=[True, False] if within_experiment else False
    )
    return list(zip(ordered["experiment"], ordered["batch"]))


def _aggregate_for_batch(long_df, experiment, batch, taxonomic_level, min_abundance):
    """Sum relabund to taxonomic_level per timepoint; lump minors to Other."""
    sub = long_df[(long_df["experiment"] == experiment) & (long_df["batch"] == batch)]
    agg = sub.groupby(["timepoint", taxonomic_level])["relabund"].sum().reset_index()
    agg[taxonomic_level] = agg[taxonomic_level].fillna("").replace("", "unclassified")
    peak = agg.groupby(taxonomic_level)["relabund"].max()
    agg.loc[
        agg[taxonomic_level].isin(peak[peak < min_abundance].index), taxonomic_level
    ] = "Other"
    return agg.groupby(["timepoint", taxonomic_level])["relabund"].sum().reset_index()


def _taxon_color_map(long_df, taxonomic_level, min_abundance):
    """One consistent taxon->color map per amplicon across all batches."""
    level = long_df[taxonomic_level].fillna("").replace("", "unclassified")
    peak = (
        long_df.assign(**{taxonomic_level: level})
        .groupby(["sample", taxonomic_level])["relabund"]
        .sum()
        .groupby(taxonomic_level)
        .max()
    )
    taxa = sorted(peak[peak >= min_abundance].index)
    palette = (
        COLORS[: len(taxa)]
        if len(taxa) <= len(COLORS)
        else sns.color_palette("husl", len(taxa)).as_hex()
    )
    color_map = dict(zip(taxa, palette))
    color_map["Other"] = "#808080"
    return color_map


def _community_panel(
    agg, timepoints, taxonomic_level, color_map, title, width, height, show_x_labels
):
    """One stacked-bar cell: x = timepoints, stacks = taxa, y = relabund %."""
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

    p = figure(
        x_range=timepoints, width=width, height=height, toolbar_location=None, tools=""
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
        source=ColumnDataSource(wide),
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


def _gas_panel(
    gas_df, sampling_days, dna_replicate, moles_unit, width, height, show_x_labels
):
    """Cumulative production for one batch: x = days, y = moles, line per gas.

    Every replicate bottle is drawn; the bottle the DNA came from is drawn
    solid and heavier, the rest faded, because only that bottle's gas record
    is the same physical vessel as the community measurement. Vertical spans
    mark destructive-sampling days.
    """
    p = figure(width=width, height=height, toolbar_location=None, tools="")
    for (molecule, replicate), series in gas_df.groupby(["Molecule", "replicate"]):
        series = series.sort_values("Days since start")
        is_dna_bottle = replicate == dna_replicate
        p.line(
            series["Days since start"],
            series["Cumulative Moles"],
            color=GAS_COLORS[molecule],
            line_width=2 if is_dna_bottle else 1,
            alpha=1.0 if is_dna_bottle else 0.25,
            legend_label=molecule,
        )
    for day in sampling_days:
        p.add_layout(
            Span(
                location=day,
                dimension="height",
                line_color="#444444",
                line_dash="dashed",
                line_width=1,
            )
        )
    p.add_layout(Title(text="Cumulative production", text_font_size="9pt"), "above")
    p.xaxis.visible = show_x_labels
    p.xaxis.axis_label = "Days since start"
    p.yaxis.axis_label = moles_unit
    p.axis.axis_label_text_font_size = "7pt"
    p.axis.major_label_text_font_size = "7pt"
    p.legend.visible = False
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
        Title(text=f"{taxonomic_level.title()} legend", text_font_size="9pt"), "above"
    )
    return p


def _gas_legend_panel(gases, height):
    """Standalone legend for the gas column."""
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
    for gas in gases:
        p.line(
            [-1, -1], [-1, -1], color=GAS_COLORS[gas], line_width=2, legend_label=gas
        )
    p.legend.location = "top_left"
    p.legend.label_text_font_size = "8pt"
    p.add_layout(Title(text="Gas (solid = DNA bottle)", text_font_size="9pt"), "above")
    return p


def create_batch_timecourse_grid(
    runs,
    batch_order,
    output_path,
    batch_metadata=None,
    gas_df=None,
    timepoint_map=None,
    moles_unit="mol",
    taxonomic_level="phylum",
    min_abundance=2.0,
    panel_width=240,
    panel_height=170,
    title=None,
):
    """Grid of community timecourses, rows = batches, columns = amplicons.

    runs: dict of column label -> long df from load_amplicon_run, in
        left-to-right column order, e.g. {"16S": df16, "18S": df18}.
    batch_order: list of (experiment, batch) int tuples, top row first.
    batch_metadata: df from load_batch_metadata; adds a_w to row titles.
    gas_df: df from load_gas_timeseries; adds a gas production column.
    timepoint_map: df from map_timepoints_to_days; supplies the DNA sampling
        days marked on the gas panels and the bottle each came from.
    min_abundance: peak per-sample relabund (%) below which a taxon lumps
        into "Other".

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

    label_lookup = {}
    if batch_metadata is not None:
        for record in batch_metadata.to_dict("records"):
            label_lookup[(record["experiment"], record["batch"])] = (
                record["water_activity"],
                record.get("Salt Makeup"),
            )

    grid = []
    for i, (experiment, batch) in enumerate(batch_order):
        last_row = i == len(batch_order) - 1
        label = label_lookup.get((experiment, batch))
        if label is None:
            aw_txt = ""
        else:
            water_activity, salt_makeup = label
            aw_txt = f"  a_w={water_activity:.3f}"
            if salt_makeup:
                aw_txt += f"  {salt_makeup}"
        row = []
        for label, df in runs.items():
            agg = _aggregate_for_batch(
                df, experiment, batch, taxonomic_level, min_abundance
            )
            row.append(
                None
                if agg.empty
                else _community_panel(
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

        if gas_df is not None:
            batch_gas = gas_df[
                (gas_df["experiment"] == experiment) & (gas_df["batch"] == batch)
            ]
            sampling = (
                pd.DataFrame(columns=["days_since_start", "replicate"])
                if timepoint_map is None
                else timepoint_map[
                    (timepoint_map["experiment"] == experiment)
                    & (timepoint_map["batch"] == batch)
                    & timepoint_map["replicate"].notna()
                ]
            )
            row.append(
                None
                if batch_gas.empty
                else _gas_panel(
                    batch_gas,
                    sampling["days_since_start"].tolist(),
                    dna_replicate=(
                        None if sampling.empty else sampling["replicate"].iloc[0]
                    ),
                    moles_unit=moles_unit,
                    width=panel_width,
                    height=panel_height,
                    show_x_labels=last_row,
                )
            )
        grid.append(row)

    legend_height = max(300, panel_height * 2)
    legend_row = [
        _legend_panel(color_maps[label], taxonomic_level, legend_height)
        for label in runs
    ]
    if gas_df is not None:
        legend_row.append(
            _gas_legend_panel(sorted(gas_df["Molecule"].unique()), legend_height)
        )
    grid.append(legend_row)

    output_file(output_path, title=title or f"Batch timecourses ({taxonomic_level})")
    save(gridplot(grid, merge_tools=True))
    return output_path
