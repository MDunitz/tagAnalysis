"""Reference constants for 16S/18S amplicon processing.

Only genuinely-constant reference data lives here (primer sequences,
taxonomic ranks, plot palette). Per-run configuration (data paths,
dataset name, reference DB path) lives in `config.RunConfig`.
"""

import re

# --- 16S rRNA primers (515F/926R, Parada et al. 2016; Orphan lab tag set) ---
# Forward 515F, reverse 926R (Y/M/R = IUPAC degenerate bases).
# Amplifies the V4-V5 region (~373 bp between primers). NB: the reverse primer
# is 926R (CCGYCAATTYMTTTRAGTTT), NOT 806R (GGACTACNVGGGTWTCTAAT); the two give
# different amplicon lengths, which the truncLen/minOverlap defaults depend on.
FWD_SEQUENCE__16s = "GTGYCAGCMGCCGCGGTAA"
REV_SEQUENCE__16s = "CCGYCAATTYMTTTRAGTTT"
# Reverse complements (verified against sequences by test_primers.py).
FWD_RC__16s = "TTACCGCGGCKGCTGRCAC"
REV_RC__16s = "AAACTYAAAKRAATTGRCGG"

# --- 18S rRNA V4 primers (Stoeck et al. 2010; Nature Sci Rep 10:6519) ---
# Forward V4 (5'-CCAGCAGCCGCGGTAATTCC-3'), reverse V4 (5'-ACTTTCGTTCTTGATTAA-3').
FWD_SEQUENCE__18s = "CCAGCAGCCGCGGTAATTCC"
REV_SEQUENCE__18s = "ACTTTCGTTCTTGATTAA"
# Reverse complements (verified against sequences by test_primers.py).
FWD_RC__18s = "GGAATTACCGCGGCTGCTGG"
REV_RC__18s = "TTAATCAAGAACGAAAGT"

# Taxonomic ranks emitted by IdTaxa classification, in order.
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# IdTaxa confidence (percent) below which a rank is left unclassified. 40 keeps
# lineages the stricter default of 60 would drop.
TAXONOMY_CONFIDENCE_THRESHOLD = 40

# decontam prevalence-method probability threshold for calling a contaminant.
DECONTAM_PREVALENCE_THRESHOLD = 0.5

# Library -> DNA extraction batch table (see library_batches.py): one row per
# library, keyed by library ID; a blank batch marks a shared control (e.g. a
# PCR negative) that applies to every batch.
LIBRARY_BATCHES_LIBRARY_COL = "library"
LIBRARY_BATCHES_BATCH_COL = "extraction_batch"

# A library ID is its forward-read fastq name without the read/extension
# suffix: the raw "<id>_R1_001.fastq.gz" and the filtered
# "<id>_R1_filtered.fastq.gz" (the counts-table column) share the same <id>,
# e.g. "GP_3_Exp03_B01_T0_S99". Keep this comment in sync with the pattern.
LIBRARY_ID_SUFFIX_RE = re.compile(r"_R1_(?:001|filtered)\.fastq\.gz$")

# Stackbar palette (teal -> green -> gold gradient family).
COLORS = [
    "#33CCCC", "#009999", "#006666", "#669999", "#76DDDA",
    "#0F793D", "#064C26", "#013333", "#354C3E", "#0F2B1B",
    "#99CC33", "#669933", "#CDDE60", "#669966", "#6ABD45",
    "#FFCC33", "#FF9900", "#CC9933", "#996600", "#724419",
]
