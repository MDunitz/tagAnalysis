"""
Bokeh plotting utilities and project conventions.

Standard figure creation helpers that enforce consistent styling
across the project. Key conventions:
    - Water column / depth profiles: x-axis on top
    - Calibration curves: include R² annotation
    - All figures: tooltips enabled by default
"""

import numpy as np
from bokeh.plotting import figure
from bokeh.models import (
    HoverTool,
    Label,
    ColumnDataSource,
)
from bokeh.palettes import Category10, Spectral

# Default figure dimensions
DEFAULT_WIDTH = 700
DEFAULT_HEIGHT = 450


def styled_figure(
    title="",
    x_label="",
    y_label="",
    width=DEFAULT_WIDTH,
    height=DEFAULT_HEIGHT,
    tooltips=None,
):
    """Create a Bokeh figure with standard project styling.

    Parameters
    ----------
    title : str
    x_label : str
    y_label : str
    width : int
    height : int
    tooltips : list of tuples, optional
        HoverTool tooltip specifications.

    Returns
    -------
    bokeh.plotting.Figure
    """
    p = figure(
        title=title,
        x_axis_label=x_label,
        y_axis_label=y_label,
        width=width,
        height=height,
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )
    p.title.text_font_size = "14pt"
    p.xaxis.axis_label_text_font_size = "12pt"
    p.yaxis.axis_label_text_font_size = "12pt"

    if tooltips:
        p.add_tools(HoverTool(tooltips=tooltips))

    return p


def depth_profile_figure(
    title="",
    x_label="",
    depth_label="Depth",
    width=DEFAULT_WIDTH,
    height=DEFAULT_HEIGHT,
    tooltips=None,
):
    """Create a figure for water column / depth profiles with x-axis on top.

    Convention: depth increases downward (y-axis inverted), measurement
    axis on top.

    Parameters
    ----------
    title : str
    x_label : str
        Label for the measurement axis (top).
    depth_label : str
        Label for the depth axis (left, inverted).
    width : int
    height : int
    tooltips : list of tuples, optional

    Returns
    -------
    bokeh.plotting.Figure
    """
    p = figure(
        title=title,
        x_axis_label=x_label,
        y_axis_label=depth_label,
        width=width,
        height=height,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        y_axis_location="left",
    )

    # Invert y-axis so depth increases downward
    p.y_range.flipped = True

    # Move x-axis to top
    p.xaxis.axis_label_text_font_size = "12pt"
    p.yaxis.axis_label_text_font_size = "12pt"
    p.above = p.below
    p.below = []

    if tooltips:
        p.add_tools(HoverTool(tooltips=tooltips))

    return p


def plot_calibration_curve(
    cal_params,
    title="Calibration Curve",
    x_label="Concentration",
    y_label="Response",
    width=DEFAULT_WIDTH,
    height=DEFAULT_HEIGHT,
):
    """Plot a calibration curve with fit line and R² annotation.

    Parameters
    ----------
    cal_params : dict
        Output of calibration.fit_linear_calibration().
    title : str
    x_label : str
    y_label : str
    width : int
    height : int

    Returns
    -------
    bokeh.plotting.Figure
    """
    conc = cal_params["concentrations"]
    resp = cal_params["responses"]

    source = ColumnDataSource(data=dict(x=conc, y=resp))

    p = styled_figure(
        title=title,
        x_label=x_label,
        y_label=y_label,
        width=width,
        height=height,
        tooltips=[("Conc", "@x{0.2f}"), ("Response", "@y{0.1f}")],
    )

    # Data points
    p.circle(
        "x",
        "y",
        source=source,
        size=8,
        color="steelblue",
        alpha=0.8,
        legend_label="Standards",
    )

    # Fit line
    x_fit = np.linspace(min(conc), max(conc), 100)
    y_fit = cal_params["slope"] * x_fit + cal_params["intercept"]
    p.line(x_fit, y_fit, color="tomato", line_width=2, legend_label="Fit")

    # R² annotation
    r2_label = Label(
        x=0.05,
        y=0.92,
        x_units="screen",
        y_units="screen",
        text=f"R² = {cal_params['r_squared']:.4f}",
        text_font_size="12pt",
        background_fill_color="white",
        background_fill_alpha=0.7,
    )
    p.add_layout(r2_label)

    # Equation annotation
    eq_label = Label(
        x=0.05,
        y=0.85,
        x_units="screen",
        y_units="screen",
        text=f"y = {cal_params['slope']:.2e}x + {cal_params['intercept']:.2e}",
        text_font_size="10pt",
        background_fill_color="white",
        background_fill_alpha=0.7,
    )
    p.add_layout(eq_label)

    p.legend.location = "bottom_right"
    return p


def get_palette(n):
    """Return a color palette with n colors.

    Parameters
    ----------
    n : int
        Number of colors needed.

    Returns
    -------
    list of str
        Hex color strings.
    """
    if n <= 10:
        return Category10[max(n, 3)][:n]
    return Spectral[min(n, 11)]
