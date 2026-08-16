"""Integration tests that exercise the external-tool stages on a small,
committed subsample of real 18S reads (tests/fixtures/integration/).

Tiers:
  - cutadapt: runs in the standard CI (cutadapt is pip-installable). Verifies
    remove_primers_cutadapt actually strips primers and produces valid output.
  - DADA2 / taxonomy: require R + Bioconductor (dada2/DECIPHER) and, for
    taxonomy, a reference DB. These are skipped unless that environment is
    present, so they run only under the conda CI tier (or locally). They are
    scaffolds: correct-by-construction but not yet executed in this repo's CI.

The subsample is intentionally small (2000 read pairs) - these test the
plumbing of each stage, NOT reproduction of the full-depth golden ASV table
(that exact regression lives in test_regression.py against the golden CSVs).
"""

import gzip
import os
import shutil
import subprocess

import pytest

from tag_analysis import etl
from tag_analysis.config import PRIMERS_18S

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "integration", "18S")

_HAS_CUTADAPT = shutil.which("cutadapt") is not None
_HAS_RSCRIPT = shutil.which("Rscript") is not None


def _has_r_package(pkg):
    """True if an R package is installed (requires Rscript)."""
    if not _HAS_RSCRIPT:
        return False
    r = subprocess.run(
        ["Rscript", "-e", f'quit(status = !requireNamespace("{pkg}", quietly = TRUE))'],
        capture_output=True,
    )
    return r.returncode == 0


def _read_count(path):
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def _first_seq(path):
    with gzip.open(path, "rt") as fh:
        next(fh)              # header
        return next(fh).strip()


def _stage_inputs(tmp_path):
    """Copy the subsample fixture into a writable input dir; return its path."""
    dst = tmp_path / "input"
    dst.mkdir()
    for f in os.listdir(FIXTURES):
        shutil.copy(os.path.join(FIXTURES, f), dst / f)
    return str(dst)


# ============================ cutadapt (runs in CI) ============================

@pytest.mark.skipif(not _HAS_CUTADAPT, reason="cutadapt not installed")
class TestCutadaptPrimerRemoval:

    def test_produces_output_pair(self, tmp_path):
        """remove_primers_cutadapt writes a trimmed R1/R2 pair (guards against the
        function's silent-failure path, which only prints on cutadapt error)."""
        data_path = _stage_inputs(tmp_path)
        out = str(tmp_path / "cutadapt")
        etl.remove_primers_cutadapt(
            data_path, PRIMERS_18S.fwd, PRIMERS_18S.rev,
            PRIMERS_18S.rev_rc, PRIMERS_18S.fwd_rc,
            output_directory=out, min_length=100,
        )
        outs = sorted(os.listdir(out))
        assert any("_R1_" in f for f in outs), f"no R1 output produced: {outs}"
        assert any("_R2_" in f for f in outs), f"no R2 output produced: {outs}"

    def test_primer_is_stripped_from_reads(self, tmp_path):
        """After trimming, R1 reads no longer begin with the forward primer stem."""
        data_path = _stage_inputs(tmp_path)
        out = str(tmp_path / "cutadapt")
        etl.remove_primers_cutadapt(
            data_path, PRIMERS_18S.fwd, PRIMERS_18S.rev,
            PRIMERS_18S.rev_rc, PRIMERS_18S.fwd_rc,
            output_directory=out, min_length=100,
        )
        r1_out = os.path.join(out, "MD-sub-20250619-18s_S236_L001_R1_001.fastq.gz")
        seq = _first_seq(r1_out)
        # Forward primer core (ignoring the leading degenerate/N base of the primer).
        assert not seq.startswith(PRIMERS_18S.fwd[1:]), (
            "forward primer still present at read start after trimming"
        )

    def test_read_count_not_inflated(self, tmp_path):
        """Trimmed output has <= input reads (--discard-untrimmed only removes)."""
        data_path = _stage_inputs(tmp_path)
        out = str(tmp_path / "cutadapt")
        etl.remove_primers_cutadapt(
            data_path, PRIMERS_18S.fwd, PRIMERS_18S.rev,
            PRIMERS_18S.rev_rc, PRIMERS_18S.fwd_rc,
            output_directory=out, min_length=100,
        )
        r1_in = os.path.join(data_path, "MD-sub-20250619-18s_S236_L001_R1_001.fastq.gz")
        r1_out = os.path.join(out, "MD-sub-20250619-18s_S236_L001_R1_001.fastq.gz")
        n_in, n_out = _read_count(r1_in), _read_count(r1_out)
        assert 0 < n_out <= n_in, f"unexpected read counts in={n_in} out={n_out}"


# ===================== DADA2 (conda tier - skipped in pip CI) =====================

@pytest.mark.skipif(
    not _has_r_package("dada2"),
    reason="R package dada2 not available (conda integration tier)",
)
class TestDada2Smoke:

    def test_dada2_emits_sequence_table(self, tmp_path):
        """Smoke test: cutadapt -> DADA2 on the subsample yields a non-empty
        sequence table whose columns are DNA. NOT an exact-count regression.

        Scaffold: correct-by-construction; not yet executed in CI (needs R/dada2).
        """
        from tag_analysis.dada2_pipeline import run_dada2_pipeline  # noqa: F401

        data_path = _stage_inputs(tmp_path)
        deprimered = str(tmp_path / "cutadapt")
        etl.remove_primers_cutadapt(
            data_path, PRIMERS_18S.fwd, PRIMERS_18S.rev,
            PRIMERS_18S.rev_rc, PRIMERS_18S.fwd_rc,
            output_directory=deprimered, min_length=100,
        )
        out = str(tmp_path / "dada2_out")
        os.makedirs(out, exist_ok=True)
        run_dada2_pipeline(
            path_to_fastq_files=deprimered,
            path_to_output_dir=out,
            dataset_name="subsample_smoke",
            truncLen=(150, 120), maxEE=(2, 2), minLen=100, multithread=2,
        )
        seqtab = os.path.join(out, "sequence_table_nochim.csv")
        assert os.path.exists(seqtab), "DADA2 did not produce sequence_table_nochim.csv"

        import pandas as pd
        df = pd.read_csv(seqtab, index_col=0)
        assert df.shape[1] > 0, "sequence table has no ASV columns"
        assert all(set(c.upper()) <= set("ACGTN") for c in df.columns), (
            "ASV column headers are not DNA sequences"
        )


# ================= taxonomy (needs R + reference DB - skipped) =================

_REF_DB = os.environ.get("TAG_ANALYSIS_REFERENCE_DB")


@pytest.mark.skipif(
    not (_has_r_package("DECIPHER") and _REF_DB and os.path.exists(_REF_DB)),
    reason="DECIPHER + TAG_ANALYSIS_REFERENCE_DB not available (conda integration tier)",
)
def test_assign_taxonomy_runs(tmp_path):
    """Scaffold: taxonomy assignment against a real reference DB.

    Blocked until a SILVA/PR2 training set is provided via the
    TAG_ANALYSIS_REFERENCE_DB env var. Not yet executed in this repo.
    """
    import pandas as pd

    mapping = pd.DataFrame(
        {"ASV_ID": ["ASV_1"], "sequence": [_first_seq(
            os.path.join(FIXTURES, "MD-sub-20250619-18s_S236_L001_R1_001.fastq.gz")
        )]}
    )
    out = str(tmp_path / "tax_out")
    os.makedirs(out, exist_ok=True)
    tax_df = etl.assign_taxonomy(out, mapping, reference_db_path=_REF_DB)
    assert "domain" in tax_df.columns
