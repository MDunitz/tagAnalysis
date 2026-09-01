"""Regression tests for stackbar plot scaling.

The original implementation created one GlyphRenderer + one-row
ColumnDataSource per (sample x taxon) bar segment: O(samples * taxa)
Bokeh objects. At 48 samples x 188 taxa that serialized 9,024 renderers
and a 20 MB HTML file that hung the browser on load. The fixed
implementation stacks from a single wide-format source: O(taxa)
renderers. These tests pin that scaling so it can't silently regress.
"""
import os

import numpy as np
import pandas as pd
import pytest
from bokeh.models import GlyphRenderer

from tag_analysis.plotting import create_stackbar_plot

N_SAMPLES = 50
N_TAXA = 200

# A per-bar implementation at this scale produces a ~20 MB file; the
# single-source implementation produces ~1 MB. 5 MB splits the two with
# margin on both sides.
MAX_HTML_BYTES = 5_000_000


@pytest.fixture
def relabund_long_df():
    """Synthetic long-format relative abundance df: every sample sums to 100%."""
    rng = np.random.default_rng(42)
    samples = [f"sample_{i:03d}" for i in range(N_SAMPLES)]
    taxa = [f"taxon_{i:03d}" for i in range(N_TAXA)]
    rows = []
    for s in samples:
        abund = rng.dirichlet(np.ones(N_TAXA)) * 100
        rows.extend((s, t, a) for t, a in zip(taxa, abund))
    return pd.DataFrame(rows, columns=["sample", "species", "relabund"])


def test_stackbar_renderer_count_scales_with_taxa_not_bars(tmp_path, relabund_long_df):
    out = str(tmp_path / "species.html")
    p = create_stackbar_plot(relabund_long_df, "species", out, min_abundance=0)

    glyph_renderers = [r for r in p.renderers if isinstance(r, GlyphRenderer)]
    assert len(glyph_renderers) == N_TAXA, (
        f"expected one renderer per taxon ({N_TAXA}), got {len(glyph_renderers)}; "
        f"per-bar regression would give {N_SAMPLES * N_TAXA}"
    )

    # All stackers must share a single ColumnDataSource.
    sources = {id(r.data_source) for r in glyph_renderers}
    assert len(sources) == 1, f"expected 1 shared ColumnDataSource, got {len(sources)}"


def test_stackbar_html_output_stays_small(tmp_path, relabund_long_df):
    out = str(tmp_path / "species.html")
    create_stackbar_plot(relabund_long_df, "species", out, min_abundance=0)

    size = os.path.getsize(out)
    assert size < MAX_HTML_BYTES, (
        f"stackbar HTML is {size / 1e6:.1f} MB (limit {MAX_HTML_BYTES / 1e6:.0f} MB); "
        "check for per-bar renderer regression"
    )


def test_stackbar_preserves_totals(tmp_path, relabund_long_df):
    out = str(tmp_path / "species.html")
    p = create_stackbar_plot(relabund_long_df, "species", out, min_abundance=0)

    source = next(
        r.data_source for r in p.renderers if isinstance(r, GlyphRenderer)
    )
    # ColumnDataSource adds the dataframe index as an "index" column;
    # exclude it and the sample labels before summing taxa columns.
    wide = pd.DataFrame(dict(source.data)).drop(columns=["sample", "index"])
    row_sums = wide.sum(axis=1).to_numpy()
    assert np.allclose(row_sums, 100.0), "stacked segments must sum to 100% per sample"


def test_stackbars_write_to_custom_img_dir(tmp_path, relabund_long_df):
    """img_dir separates regenerable images from the analysis-product dir."""
    out_dir = tmp_path / "transformed"
    img_dir = tmp_path / "outputs"
    out_dir.mkdir()

    from tag_analysis.plotting import create_relative_abundance_stackbars

    plots = create_relative_abundance_stackbars(
        relabund_long_df, str(out_dir),
        taxonomic_levels=["species"],
        img_dir=str(img_dir),
    )
    assert plots == [str(img_dir / "relative_abundance_species.html")]
    assert (img_dir / "relative_abundance_species.html").exists()
    # nothing image-shaped leaked into the output dir
    assert list(out_dir.iterdir()) == []


def test_runconfig_img_dir_default_and_override(tmp_path):
    from tag_analysis import RunConfig

    default_cfg = RunConfig(
        data_path=str(tmp_path / "d"), output_path=str(tmp_path / "o"),
        dataset_name="x", reference_db_path="ref.RData",
    )
    assert default_cfg.img_dir == str(tmp_path / "o" / "imgs")

    custom = RunConfig(
        data_path=str(tmp_path / "d"), output_path=str(tmp_path / "o"),
        dataset_name="x", reference_db_path="ref.RData",
        img_dir=str(tmp_path / "outputs"),
    )
    assert custom.img_dir == str(tmp_path / "outputs")
