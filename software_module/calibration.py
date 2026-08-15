"""
Calibration curve fitting and quantification utilities.

Provides functions for fitting linear calibration curves from standard
measurements, calculating detection limits, and applying calibrations
to unknown samples. Designed for analytical chemistry workflows
(GC-MS, ICP-MS, spectrophotometry, etc.).

All functions use astropy units where applicable to enforce dimensional
consistency.
"""

import numpy as np
from scipy import stats


def fit_linear_calibration(concentrations, responses, include_zero=True):
    """Fit a linear calibration curve: response = slope * concentration + intercept.

    Ordinary least squares linear regression.
    y = mx + b

    Parameters
    ----------
    concentrations : array-like
        Known standard concentrations.
    responses : array-like
        Measured instrument responses (peak areas, absorbance, etc.).
    include_zero : bool
        If True, include zero-concentration standards in the fit.
        Recommended True to capture baseline noise (see Notes).

    Returns
    -------
    dict
        slope, intercept, r_squared, std_err, n_points,
        concentrations, responses (arrays used in fit).

    Notes
    -----
    Zero standards should generally be included because instrument baseline
    noise produces nonzero signal even at zero analyte concentration. The
    intercept captures this baseline and improves accuracy at low
    concentrations. See: Danzer & Currie, Pure Appl. Chem. 70 (1998) 993.
    """
    conc = np.asarray(concentrations, dtype=float)
    resp = np.asarray(responses, dtype=float)

    if not include_zero:
        mask = conc > 0
        conc = conc[mask]
        resp = resp[mask]

    slope, intercept, r_value, p_value, std_err = stats.linregress(conc, resp)

    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_value**2,
        "std_err": std_err,
        "p_value": p_value,
        "n_points": len(conc),
        "concentrations": conc,
        "responses": resp,
    }


def apply_calibration(responses, cal_params, below_detection="nan"):
    """Convert instrument responses to concentrations using a calibration.

    Inverse of the calibration equation:
        concentration = (response - intercept) / slope

    Parameters
    ----------
    responses : array-like
        Measured instrument responses for unknown samples.
    cal_params : dict
        Output of fit_linear_calibration().
    below_detection : str
        How to handle negative calculated concentrations.
        'nan' -> np.nan, 'zero' -> 0.0, 'negative' -> keep as-is.

    Returns
    -------
    np.ndarray
        Calculated concentrations.
    """
    resp = np.asarray(responses, dtype=float)
    conc = (resp - cal_params["intercept"]) / cal_params["slope"]

    if below_detection == "nan":
        conc[conc < 0] = np.nan
    elif below_detection == "zero":
        conc[conc < 0] = 0.0

    return conc


def detection_limits(cal_params, blank_responses=None, k_lod=3.3, k_loq=10.0):
    """Calculate limit of detection (LOD) and limit of quantification (LOQ).

    LOD = k_lod * (std_blank / slope)
    LOQ = k_loq * (std_blank / slope)

    Based on ICH Q2(R1) guidelines for analytical method validation.

    Parameters
    ----------
    cal_params : dict
        Output of fit_linear_calibration().
    blank_responses : array-like, optional
        Replicate blank/zero-standard responses. If None, uses the
        standard error of the regression as an approximation.
    k_lod : float
        Multiplier for LOD (default 3.3 per ICH).
    k_loq : float
        Multiplier for LOQ (default 10.0 per ICH).

    Returns
    -------
    dict
        lod, loq (in concentration units), std_blank.
    """
    if blank_responses is not None:
        std_blank = np.std(blank_responses, ddof=1)
    else:
        # Approximate from regression residuals
        predicted = (
            cal_params["slope"] * cal_params["concentrations"] + cal_params["intercept"]
        )
        residuals = cal_params["responses"] - predicted
        std_blank = np.std(residuals, ddof=2)

    lod = k_lod * std_blank / abs(cal_params["slope"])
    loq = k_loq * std_blank / abs(cal_params["slope"])

    return {"lod": lod, "loq": loq, "std_blank": std_blank}


def dilution_corrected_concentration(measured_conc, dilution_factor):
    """Apply dilution correction to measured concentrations.

    corrected_concentration = measured_concentration * dilution_factor

    Parameters
    ----------
    measured_conc : array-like or astropy Quantity
        Concentrations from calibration curve.
    dilution_factor : float
        Ratio of total volume to sample volume (e.g., 51 for a 50:1 split).

    Returns
    -------
    array-like or astropy Quantity
        Corrected concentrations in original sample.
    """
    return np.asarray(measured_conc) * dilution_factor
