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

import numpy as np
import pandas as pd
from astropy import constants as const
from astropy import units as u
from astropy.table import Table

SAMPLE_NAME_RE = re.compile(
    r"GP_\d+_(?P<experiment>Exp_?\d+)_(?P<batch>B\d+)_(?P<timepoint>T\d+)"
    r"(?P<wash>_wash)?_S\d+"
)

# Gas channels treated as production series. CH4 uses the FID channel; the
# MS m/z 15 channel is diagnostic-only in saltyBiomass
# (incubations.constants.PRODUCTION_GASES) and is excluded for the same reason.
PRODUCTION_GASES = ["CO2", "N2O", "CH4_FID"]

# Partial molar volume of water at 25 C, used to convert the WP4C water
# potential reading to water activity. Pure-water value: the brine's true
# partial molar volume differs, which is a known approximation in the
# meter-to-a_w conversion (saltyBiomass #539).
WATER_MOLAR_VOLUME = 1.8068e-5 * u.m**3 / u.mol
WP4C_REFERENCE_TEMPERATURE = 298.15 * u.K

# WP4C stated accuracy: +/-0.05 MPa from 0 to -5 MPa, +/-1% of reading
# beyond that. Used to decide whether two batches are matched in a_w by
# design or genuinely separated.
WP4C_LOW_RANGE_ACCURACY = 0.05 * u.MPa
WP4C_LOW_RANGE_LIMIT = 5.0 * u.MPa
WP4C_HIGH_RANGE_FRACTION = 0.01

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


def water_activity_to_potential(water_activity, temperature=WP4C_REFERENCE_TEMPERATURE):
    """Kelvin equation, inverted: water activity -> water potential.

        a_w = exp(psi * V_m / (R * T))   =>   psi = ln(a_w) * R * T / V_m

    a_w         : water activity (dimensionless, 0-1)
    psi         : water potential [MPa], negative
    V_m         : partial molar volume of water [m^3/mol]
    R           : molar gas constant
    T           : sample temperature [K]

    The WP4C measures psi directly; a_w in the batch sheet is the forward
    conversion. Working in psi is what makes meter accuracy comparable
    across the a_w range, since a fixed psi error maps to an a_w error
    that grows ~7x from a_w=1 to a_w=0.7.
    """
    return (np.log(water_activity) * const.R * temperature / WATER_MOLAR_VOLUME).to(
        u.MPa
    )


def meter_accuracy(potential):
    """WP4C 1-sigma accuracy at a given water potential [MPa]."""
    magnitude = np.abs(potential)
    return (
        np.where(
            magnitude <= WP4C_LOW_RANGE_LIMIT,
            WP4C_LOW_RANGE_ACCURACY.value,
            WP4C_HIGH_RANGE_FRACTION * magnitude.value,
        )
        * u.MPa
    )


def matched_aw_groups(batch_metadata, temperature=WP4C_REFERENCE_TEMPERATURE):
    """Label batches that share a water activity within meter accuracy.

    Two batches are treated as matched when their water potentials differ
    by less than the combined WP4C accuracy, sigma_combined =
    sqrt(sigma_1^2 + sigma_2^2). Grouping is within an experiment, since
    the matched-a_w contrasts (sulfate vs. none, Na vs. Mg) are designed
    within an experiment, not across.

    Returns batch_metadata with added columns:
      psi_MPa, meter_sigma_MPa, aw_group (int, 1-based within experiment)
    Batches in the same aw_group are a designed matched pair; a pair the
    design intends but that lands in different groups is not matched at
    the resolution the meter provides.
    """
    meta = batch_metadata.copy()
    meta["psi_MPa"] = water_activity_to_potential(
        meta["water_activity"].to_numpy(), temperature
    ).value
    meta["meter_sigma_MPa"] = meter_accuracy(meta["psi_MPa"].to_numpy() * u.MPa).value

    groups = []
    for experiment, block in meta.groupby("experiment", sort=True):
        block = block.sort_values("psi_MPa", ascending=False)
        group_id, previous = 1, None
        labels = []
        for row in block.itertuples():
            if previous is not None:
                combined = np.hypot(previous.meter_sigma_MPa, row.meter_sigma_MPa)
                if abs(row.psi_MPa - previous.psi_MPa) >= combined:
                    group_id += 1
            labels.append(group_id)
            previous = row
        groups.append(block.assign(aw_group=labels))
    return pd.concat(groups, ignore_index=True)


def aw_match_report(batch_metadata, group_column="Salt Makeup"):
    """Per matched-a_w group, the contrast it supports and its psi spread.

    A group with one member supports no within-a_w contrast. A group whose
    members all share the same `group_column` value is a replicate, not a
    contrast.
    """
    grouped = matched_aw_groups(batch_metadata)
    return (
        grouped.groupby(["experiment", "aw_group"])
        .agg(
            batches=("batch", lambda s: sorted(s)),
            contrast=(group_column, lambda s: sorted(set(s))),
            psi_spread_MPa=("psi_MPa", lambda s: s.max() - s.min()),
            aw_range=("water_activity", lambda s: (s.min(), s.max())),
        )
        .reset_index()
    )
