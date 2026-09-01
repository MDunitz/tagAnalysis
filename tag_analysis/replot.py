"""CLI for regenerating relative-abundance plots from persisted pipeline
outputs, without rerunning the pipeline.

Usage:
    python -m tag_analysis.replot INPUT_DIR [--counts-file NAME] [--levels ...]

INPUT_DIR is a pipeline output directory containing the counts and
taxonomy CSVs. Plots are written to INPUT_DIR/imgs/, matching the
pipeline's layout. Prefers the decontaminated counts file
(ASVs_counts_clean.csv) when present, falling back to ASVs_counts.csv.
"""
import argparse
import os

from .constants import RANKS, COLORS
from .plotting import replot_relative_abundance

CLEAN_COUNTS = "ASVs_counts_clean.csv"
RAW_COUNTS = "ASVs_counts.csv"
TAXONOMY = "ASV_taxonomy.csv"


def default_counts_file(input_dir):
    """Prefer decontaminated counts when present, matching pipeline behavior."""
    if os.path.exists(os.path.join(input_dir, CLEAN_COUNTS)):
        return CLEAN_COUNTS
    return RAW_COUNTS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", help="Pipeline output dir with counts + taxonomy CSVs")
    parser.add_argument("--counts-file", default=None,
                        help=f"Counts filename within INPUT_DIR (default: {CLEAN_COUNTS} if present, else {RAW_COUNTS})")
    parser.add_argument("--levels", nargs="+", default=RANKS,
                        help=f"Taxonomic levels to plot (default: {RANKS})")
    parser.add_argument("--dataset-name", default=None, help="Title prefix for plots")
    args = parser.parse_args()

    counts_file = args.counts_file or default_counts_file(args.input_dir)
    plots = replot_relative_abundance(
        counts_file_path=os.path.join(args.input_dir, counts_file),
        taxonomy_file_path=os.path.join(args.input_dir, TAXONOMY),
        output_dir=args.input_dir,
        taxonomic_levels=args.levels,
        colors=COLORS,
        dataset_name=args.dataset_name,
    )
    print(f"Regenerated {len(plots)} plots from {counts_file}")


if __name__ == "__main__":
    main()
