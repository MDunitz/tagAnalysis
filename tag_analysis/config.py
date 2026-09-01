"""Per-run configuration for the amplicon pipeline.

Per-run values (data path, output path, dataset name, reference DB) are
carried on RunConfig rather than module-level globals, so no project needs
to edit package source to configure a run.
"""

import os
from dataclasses import dataclass, field

from . import constants


@dataclass(frozen=True)
class PrimerSet:
    """A forward/reverse primer pair with precomputed reverse complements.

    fwd / rev:        primer sequences (5'->3'), may contain IUPAC degenerate bases
    fwd_rc / rev_rc:  reverse complements, used by cutadapt for read-through trimming
    """
    name: str
    fwd: str
    rev: str
    fwd_rc: str
    rev_rc: str


PRIMERS_16S = PrimerSet(
    name="16S",
    fwd=constants.FWD_SEQUENCE__16s,
    rev=constants.REV_SEQUENCE__16s,
    fwd_rc=constants.FWD_RC__16s,
    rev_rc=constants.REV_RC__16s,
)

PRIMERS_18S = PrimerSet(
    name="18S",
    fwd=constants.FWD_SEQUENCE__18s,
    rev=constants.REV_SEQUENCE__18s,
    fwd_rc=constants.FWD_RC__18s,
    rev_rc=constants.REV_RC__18s,
)


@dataclass
class RunConfig:
    """Everything a single amplicon run needs, with no global state.

    data_path:      directory of raw fastq.gz for this amplicon
    output_path:    directory for pipeline outputs (created if absent)
    img_dir:        directory for generated plots/images. Defaults to
                    output_path/imgs; set it elsewhere (e.g. an untracked
                    outputs/ dir) to keep small analysis products and large
                    regenerable images in different locations.
    dataset_name:   label used in plot titles / filenames
    reference_db_path: path to the taxonomy training set (e.g. SILVA .RData, PR2).
                    Made explicit so each run chooses its own reference.
    primers:        the PrimerSet used to generate these reads. Optional: if left
                    None, process_16s / process_18s inject their standard pair, so
                    existing calls keep working. Set it explicitly (or use the
                    generic process()) to run any other primer pair without
                    editing package source.
    binned_quality_bins: for NextSeq/NovaSeq data whose quality scores are binned,
                    the bin values Illumina collapses Q-scores onto (e.g.
                    [2, 12, 24, 40]). When set, DADA2 error learning uses a binned
                    error function. Leave None for unbinned MiSeq data. Confirm the
                    exact bins against the quality-profile heatmap per run.
    """
    data_path: str
    output_path: str
    dataset_name: str
    reference_db_path: str
    primers: "PrimerSet | None" = None
    binned_quality_bins: "list | None" = None

    # Derived paths (populated in __post_init__ so callers never build them).
    deprimered_path: str = field(init=False)
    counts_file_path: str = field(init=False)
    taxonomy_file_path: str = field(init=False)

    img_dir: "str | None" = None

    def __post_init__(self):
        self.deprimered_path = os.path.join(self.data_path, "fastq_cutadapt")
        self.counts_file_path = os.path.join(self.output_path, "ASVs_counts.csv")
        self.taxonomy_file_path = os.path.join(self.output_path, "ASV_taxonomy.csv")
        if self.img_dir is None:
            self.img_dir = os.path.join(self.output_path, "imgs")

    def ensure_dirs(self):
        """Create output directories if they do not yet exist."""
        os.makedirs(self.output_path, exist_ok=True)
        os.makedirs(self.img_dir, exist_ok=True)
        os.makedirs(self.deprimered_path, exist_ok=True)
