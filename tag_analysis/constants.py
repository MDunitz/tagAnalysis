"""Reference constants for 16S/18S amplicon processing.

Only genuinely-constant reference data lives here (primer sequences,
taxonomic ranks, plot palette). Per-run configuration (data paths,
dataset name, reference DB path) lives in `config.RunConfig`.
"""

# --- 16S rRNA primers (515F-Y / 926R; Parada et al. 2016, Environ. Microbiol.) ---
# 515F-Y forward, 926R reverse -> V4-V5 (~374 bp amplicon). 926R differs from the
# V4-only 806R (GGACTACNVGGGTWTCTAAT); do not confuse the two. (Y/M/R = IUPAC.)
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
