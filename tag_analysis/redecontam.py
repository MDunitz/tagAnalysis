"""CLI: rerun decontamination on a finished run's persisted outputs.

Usage:
    python -m tag_analysis.redecontam OUTPUT_DIR [--library-batches CSV] [--levels ...]

Decontam reads only the pre-decontam counts (ASVs_counts.csv) and the
taxonomy (ASV_taxonomy.csv), so a change to how contaminants are called -- e.g.
adding a library -> extraction batch table -- does not need primer removal,
DADA2 or taxonomy rerun. This rewrites, in OUTPUT_DIR:

- ASVs_counts_clean.csv and contam_results.csv (plus per-batch files when
  batch-aware),
- the contamination plot and relative-abundance stackbars (in --img-dir,
  default OUTPUT_DIR/imgs),
- the decontamination section of run_manifest.json.
"""

import argparse
import os

import pandas as pd

from .constants import COLORS, RANKS
from .decontaminate import remove_contaminants
from .etl import prepare_data_for_contamination_plot
from .manifest import decontamination_record, record_redecontamination
from .plotting import create_contamination_plot, replot_relative_abundance

RAW_COUNTS = "ASVs_counts.csv"
CLEAN_COUNTS = "ASVs_counts_clean.csv"
TAXONOMY = "ASV_taxonomy.csv"


def redecontaminate(output_dir, library_batches=None, levels=RANKS, img_dir=None, dataset_name=None):
    """Rerun decontam and its downstream products on persisted outputs; returns contaminant ASVs."""
    img_dir = img_dir or os.path.join(output_dir, "imgs")
    os.makedirs(img_dir, exist_ok=True)
    counts_path = os.path.join(output_dir, RAW_COUNTS)
    taxonomy_path = os.path.join(output_dir, TAXONOMY)

    _, _, contam_asvs, predicted_controls = remove_contaminants(
        counts_path, taxonomy_path, output_dir,
        clean_count_file=CLEAN_COUNTS, library_batches=library_batches,
    )

    counts_df = pd.read_csv(counts_path, index_col=0, sep="\t")
    relative_df = counts_df.div(counts_df.sum(axis=0), axis=1) * 100
    create_contamination_plot(
        prepare_data_for_contamination_plot(relative_df, contam_asvs, predicted_controls),
        os.path.join(img_dir, "Contamination_plot.html"),
    )
    replot_relative_abundance(
        counts_file_path=os.path.join(output_dir, CLEAN_COUNTS),
        taxonomy_file_path=taxonomy_path,
        output_dir=output_dir,
        taxonomic_levels=levels,
        colors=COLORS,
        dataset_name=dataset_name,
        img_dir=img_dir,
    )
    record_redecontamination(
        output_dir,
        decontamination_record(library_batches, counts_df, contam_asvs, predicted_controls),
    )
    return contam_asvs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_dir", help="Finished run's output dir (ASVs_counts.csv + ASV_taxonomy.csv)")
    parser.add_argument("--library-batches", default=None,
                        help="Library -> extraction batch CSV for batch-aware decontam (omit for pooled)")
    parser.add_argument("--levels", nargs="+", default=RANKS,
                        help=f"Taxonomic levels for the stackbars (default: {RANKS})")
    parser.add_argument("--dataset-name", default=None, help="Title prefix for plots")
    parser.add_argument("--img-dir", default=None, help="Plot directory (default: OUTPUT_DIR/imgs)")
    args = parser.parse_args(argv)
    contam = redecontaminate(
        args.output_dir, args.library_batches, args.levels, args.img_dir, args.dataset_name
    )
    print(f"Removed {len(contam)} contaminant ASVs; outputs rewritten in {args.output_dir}")


if __name__ == "__main__":
    main()
