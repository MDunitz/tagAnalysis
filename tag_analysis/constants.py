"""Reference constants for 16S/18S amplicon processing.

Only genuinely-constant reference data lives here (primer sequences,
taxonomic ranks, plot palette). Per-run configuration (data paths,
dataset name, reference DB path) lives in `config.RunConfig`.
"""

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

# Stackbar palette (teal -> green -> gold gradient family).
COLORS = [
    "#33CCCC", "#009999", "#006666", "#669999", "#76DDDA",
    "#0F793D", "#064C26", "#013333", "#354C3E", "#0F2B1B",
    "#99CC33", "#669933", "#CDDE60", "#669966", "#6ABD45",
    "#FFCC33", "#FF9900", "#CC9933", "#996600", "#724419",
]
