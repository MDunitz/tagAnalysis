"""Loaders joining amplicon runs to saltyBiomass incubation data.

Three sources meet here, keyed on canonical integer (experiment, batch):

1. Amplicon runs   -- DADA2 sequence table + ASV taxonomy
2. Batch metadata  -- water activity, salt composition, from ISQ_DATA_*.ecsv
3. Gas production  -- cumulative moles vs days since start, same .ecsv

The .ecsv files are per-experiment (one per pipeline run per experiment),
so every loader takes a path or an iterable of paths and concatenates.

ID spellings differ across sources and are canonicalized on load:
  fastq sample names  Exp03_B01_T0    -> experiment 3, batch 1
  .ecsv Experiment    "Exp004"        -> 4
  .ecsv Batch ID      int64 1         -> 1
  GC Sample IDs       Exp_04_B01_R04  -> experiment 4, batch 1, replicate 4
"""

import re

import pandas as pd
from astropy.table import Table

SAMPLE_NAME_RE = re.compile(
    r"GP_\d+_(?P<experiment>Exp_?\d+)_(?P<batch>B\d+)_(?P<timepoint>T\d+)"
    r"(?P<wash>_wash)?_S\d+"
)

# Gas channels treated as production series. CH4 uses the FID channel; the
# MS m/z 15 channel is diagnostic-only in saltyBiomass
# (incubations.constants.PRODUCTION_GASES) and is excluded for the same reason.
PRODUCTION_GASES = ["CO2", "N2O", "CH4_FID"]

# Sequencing timepoint taken at batch setup, before the incubation clock
# starts. It has no destructive-sampling row in the pressure sheet because
# no headspace existed yet.
T0_DAYS_SINCE_START = 0


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

    Returns None for anything that is not an experiment batch sample (PCR
    negatives, extraction negatives, algae brick reads, washes) so callers
    drop non-timecourse samples with a single filter.
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

    Returns a long df with columns:
      sample, ASV, relabund, experiment, batch, timepoint, <RANKS...>
    """
    seqtab = pd.read_csv(seqtab_path, index_col=0)
    taxonomy = pd.read_csv(taxonomy_path, sep="\t")

    counts = seqtab.rename(columns=dict(zip(taxonomy["sequence"], taxonomy["ASV_ID"])))
    relabund = counts.div(counts.sum(axis=1), axis=0) * 100

    long_df = (
        relabund.stack().rename("relabund").rename_axis(["sample", "ASV"]).reset_index()
    )

    parsed = long_df["sample"].map(parse_sample_name)
    long_df = long_df[parsed.notna()].copy()
    long_df = long_df.join(
        pd.DataFrame(list(parsed.dropna()), index=parsed.dropna().index)
    )

    rank_cols = [
        c
        for c in taxonomy.columns
        if c not in ("ASV_ID", "sequence", "taxonomy", "confidence")
    ]
    return long_df.merge(
        taxonomy[["ASV_ID"] + rank_cols], left_on="ASV", right_on="ASV_ID"
    ).drop(columns="ASV_ID")


def load_isq_ecsv(paths):
    """Read one or more saltyBiomass ISQ_DATA_*.ecsv into pandas.

    astropy's to_pandas() drops Quantity units, so they are captured from
    the Table metadata first.

    Returns (df, units): units maps column name -> unit string for every
    column carrying one; df gains integer `experiment` / `batch` /
    `replicate` key columns.
    """
    if isinstance(paths, str):
        paths = [paths]

    frames, units = [], {}
    for path in paths:
        table = Table.read(path)
        units.update(
            {
                name: str(col.unit)
                for name, col in table.columns.items()
                if getattr(col, "unit", None) is not None
            }
        )
        frames.append(table.to_pandas())

    df = pd.concat(frames, ignore_index=True)
    df["experiment"] = df["Experiment"].map(canonical_experiment)
    df["batch"] = df["Batch ID"].map(canonical_batch)
    df["replicate"] = df["Replicate ID"].astype(int)
    return df, units


def load_batch_metadata(paths):
    """Per-batch metadata (one row per experiment x batch) from the .ecsv.

    Water Activity is batch-level in saltyBiomass (broadcast from the batch
    sheet), so it collapses to one row per batch. A batch carrying more than
    one distinct a_w is an upstream data error and fails here rather than
    silently picking one.
    """
    df, _ = load_isq_ecsv(paths)
    keep = ["experiment", "batch", "Water Activity"]
    keep += [
        c
        for c in ("Salt Composition", "Salt Makeup", "Contains Sulfate")
        if c in df.columns
    ]
    meta = df[keep].drop_duplicates()

    conflicted = meta.groupby(["experiment", "batch"])["Water Activity"].nunique()
    conflicted = conflicted[conflicted > 1]
    if not conflicted.empty:
        raise ValueError(
            f"Multiple Water Activity values per batch: {conflicted.to_dict()}"
        )
    return meta.rename(columns={"Water Activity": "water_activity"})


def load_gas_timeseries(paths, gases=None):
    """Per-replicate gas production timeseries from the .ecsv.

    Rows are individual measurements; replicates within a batch stay
    separate so the caller decides how to aggregate against pooled DNA.

    Returns (df, units) with columns experiment, batch, replicate, Molecule,
    Days since start, Cumulative Moles [mol], Date.
    """
    df, units = load_isq_ecsv(paths)
    gases = PRODUCTION_GASES if gases is None else gases
    cols = [
        "experiment",
        "batch",
        "replicate",
        "Molecule",
        "Days since start",
        "Cumulative Moles",
        "Date",
    ]
    return df[df["Molecule"].isin(gases)][cols].copy(), units


def load_sampling_events(paths):
    """Destructive-sampling events per replicate bottle, from the .ecsv.

    The pressure sheet's `Sampled` flag marks a *destructive* sample: the
    bottle was opened, sampled, and flushed with N2 (saltyBiomass
    `moles.calculate_post_degassing_moles`, case 2). These are the events
    that produced DNA aliquots, and they carry a Replicate ID -- the only
    record tying a sequencing timepoint to a specific bottle.

    Returns one row per (experiment, batch, replicate, date), ordered, with
    a 1-based `event_index` per batch.
    """
    df, _ = load_isq_ecsv(paths)
    events = (
        df[df["Sampled"]][
            ["experiment", "batch", "replicate", "Date", "Days since start"]
        ]
        .drop_duplicates()
        .sort_values(["experiment", "batch", "Date"])
        .reset_index(drop=True)
    )
    events["event_index"] = (
        events.groupby(["experiment", "batch"])["Date"].rank(method="dense").astype(int)
    )
    return events


def map_timepoints_to_days(sampling_events):
    """Map sequencing timepoint labels onto incubation days.

    T0 is the setup sample at day `T0_DAYS_SINCE_START`; T1..Tn map onto the
    batch's destructive-sampling events in date order (event_index n -> Tn).

    Returns a df with experiment, batch, timepoint, days_since_start, date,
    replicate -- `replicate` being the bottle the DNA came from (NaN for T0,
    which predates any bottle-specific sampling).
    """
    mapped = sampling_events.assign(
        timepoint="T" + sampling_events["event_index"].astype(str)
    ).rename(columns={"Days since start": "days_since_start", "Date": "date"})

    t0 = (
        sampling_events[["experiment", "batch"]]
        .drop_duplicates()
        .assign(
            timepoint="T0",
            days_since_start=T0_DAYS_SINCE_START,
            date=pd.NaT,
            replicate=pd.NA,
        )
    )
    cols = [
        "experiment",
        "batch",
        "timepoint",
        "days_since_start",
        "date",
        "replicate",
    ]
    timepoint_map = (
        pd.concat([t0, mapped[cols]], ignore_index=True)
        .sort_values(["experiment", "batch", "days_since_start"])
        .reset_index(drop=True)
    )

    # More than one row per batch-timepoint means two replicate bottles were
    # destructively sampled on the same date, so the timepoint label alone
    # cannot say which bottle the DNA came from. Sequencing names carry no
    # replicate, so this is unresolvable here rather than merely ambiguous.
    duplicated = timepoint_map[
        timepoint_map.duplicated(["experiment", "batch", "timepoint"], keep=False)
    ]
    if not duplicated.empty:
        raise ValueError(
            "Multiple replicate bottles sampled on the same date; timepoint "
            f"cannot be attributed to one bottle:\n{duplicated.to_string(index=False)}"
        )
    return timepoint_map


def unmapped_timepoints(amplicon_df, timepoint_map):
    """Amplicon (experiment, batch, timepoint) keys absent from the map.

    A non-empty result means sequencing exists for a batch-timepoint with no
    corresponding destructive-sampling record in the pressure sheet, so that
    sample cannot be placed on the incubation day axis.
    """
    seq_keys = amplicon_df[["experiment", "batch", "timepoint"]].drop_duplicates()
    return seq_keys.merge(
        timepoint_map[["experiment", "batch", "timepoint", "days_since_start"]],
        on=["experiment", "batch", "timepoint"],
        how="left",
    ).query("days_since_start.isna()")[["experiment", "batch", "timepoint"]]
