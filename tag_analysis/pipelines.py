"""Parameterized 16S/18S pipeline entry points.

Runs are driven entirely by RunConfig + PrimerSet, so downstream projects
call process_16s(config) / process_18s(config) without touching source.

16S and 18S are identical except for the primer set and (optionally) the
reference DB, so both delegate to _run().
"""

import os
import pandas as pd

from .config import RunConfig, PrimerSet, PRIMERS_16S, PRIMERS_18S
from .constants import COLORS, RANKS
from .etl import (
    remove_primers_cutadapt,
    create_asv_outputs,
    assign_taxonomy,
    prepare_data_for_contamination_plot,
    prepare_relative_abundance_data,
)
from .dada2_pipeline import run_dada2_pipeline
from .decontaminate import remove_contaminants
from .plotting import (
    create_contamination_plot,
    create_relative_abundance_stackbars,
)

# DADA2 defaults carried over verbatim from the original entry scripts.
_DEFAULT_DADA2_KWARGS = dict(
    quality_plot_count=3,
    maxEE=(2, 2),
    truncLen=(230, 200),
    minLen=150,
    multithread=8,
    plot_output="read_count_tracking.html",
)


def _run(config: RunConfig, clean_count_file="ASVs_counts_clean.csv",
         dada2_kwargs=None):
    """Full amplicon pipeline for one run. Behavior-equivalent to the
    original process_reads.py, parameterized by config (primers come from
    config.primers)."""
    if config.primers is None:
        raise ValueError(
            "config.primers is not set. Use process_16s/process_18s, or set "
            "config.primers to a PrimerSet and call process()."
        )
    primers = config.primers
    config.ensure_dirs()
    dada2_kwargs = {**_DEFAULT_DADA2_KWARGS, **(dada2_kwargs or {})}
    # RunConfig.binned_quality_bins sets the dada2 default for the dataset; an
    # explicit dada2_kwargs["binned_quality_bins"] (if passed) still wins.
    dada2_kwargs.setdefault("binned_quality_bins", config.binned_quality_bins)

    # 1. Primer removal (cutadapt)
    remove_primers_cutadapt(
        config.data_path,
        primers.fwd, primers.rev, primers.rev_rc, primers.fwd_rc,
        output_directory=config.deprimered_path,
        min_length=100,
    )

    # 2. DADA2: filter/trim -> denoise -> merge -> chimera removal -> summary
    results = run_dada2_pipeline(
        path_to_fastq_files=config.deprimered_path,
        path_to_output_dir=config.output_path,
        dataset_name=config.dataset_name,
        img_dir=config.img_dir,
        **dada2_kwargs,
    )

    # 3. ASV outputs (counts table, fasta, mapping)
    mapping_df, clean_asvs = create_asv_outputs(config.output_path)

    # 4. Taxonomy assignment against the run's reference DB
    taxonomy_df = assign_taxonomy(
        config.output_path, mapping_df, reference_db_path=config.reference_db_path
    )

    # 5. Decontamination (decontam via R)
    clean_counts_df, clean_relative_df, contam_asvs, predicted_controls = remove_contaminants(
        config.counts_file_path, config.taxonomy_file_path, config.output_path,
        clean_count_file=clean_count_file,
    )

    # 6. Contamination plot
    counts_df = pd.read_csv(config.counts_file_path, index_col=0, sep="\t")
    relative_df = counts_df.div(counts_df.sum(axis=0), axis=1) * 100
    contamination_df = prepare_data_for_contamination_plot(
        relative_df, contam_asvs, predicted_controls
    )
    create_contamination_plot(
        contamination_df, os.path.join(config.img_dir, "Contamination_plot.html")
    )

    # 7. Relative-abundance stackbars
    relative_abundance_df = prepare_relative_abundance_data(
        os.path.join(config.output_path, clean_count_file), config.taxonomy_file_path
    )
    create_relative_abundance_stackbars(
        relative_abundance_df, config.output_path,
        taxonomic_levels=RANKS, colors=COLORS,
        img_dir=config.img_dir,
    )

    return {
        "results": results,
        "mapping_df": mapping_df,
        "taxonomy_df": taxonomy_df,
        "clean_counts_df": clean_counts_df,
        "contam_asvs": contam_asvs,
        "predicted_controls": predicted_controls,
    }


def process(config: RunConfig, **kwargs):
    """Run the full pipeline for any primer pair. Requires config.primers to be
    set to a PrimerSet. Use this for non-standard primers; process_16s/process_18s
    are convenience wrappers that inject the standard pairs."""
    return _run(config, **kwargs)


def process_16s(config: RunConfig, **kwargs):
    """Run the full 16S pipeline. Injects the standard 515F-Y/926R pair when
    config.primers is unset; an explicitly set config.primers is respected.
    config.reference_db_path should point at a 16S training set (e.g. SILVA)."""
    if config.primers is None:
        config.primers = PRIMERS_16S
    return _run(config, **kwargs)


def process_18s(config: RunConfig, **kwargs):
    """Run the full 18S pipeline. Injects the standard V4 pair when config.primers
    is unset; an explicitly set config.primers is respected. config.reference_db_path
    should point at an 18S training set (PR2 recommended; SILVA SSU also contains 18S)."""
    if config.primers is None:
        config.primers = PRIMERS_18S
    return _run(config, **kwargs)
