"""
Tests for the software_module calibration utilities.
"""

import numpy as np
import pytest
from software_module.calibration import (
    fit_linear_calibration,
    apply_calibration,
    detection_limits,
    dilution_corrected_concentration,
)


class TestFitLinearCalibration:
    """Tests for fit_linear_calibration."""

    def test_perfect_linear_data(self):
        conc = [0, 10, 20, 30, 40, 50]
        resp = [100, 200, 300, 400, 500, 600]  # slope=10, intercept=100
        result = fit_linear_calibration(conc, resp)

        assert result["r_squared"] == pytest.approx(1.0, abs=1e-10)
        assert result["slope"] == pytest.approx(10.0, abs=1e-10)
        assert result["intercept"] == pytest.approx(100.0, abs=1e-10)
        assert result["n_points"] == 6

    def test_exclude_zeros(self):
        conc = [0, 0, 10, 20, 30]
        resp = [5, 3, 104, 203, 298]
        result = fit_linear_calibration(conc, resp, include_zero=False)

        assert result["n_points"] == 3

    def test_include_zeros(self):
        conc = [0, 0, 10, 20, 30]
        resp = [5, 3, 104, 203, 298]
        result = fit_linear_calibration(conc, resp, include_zero=True)

        assert result["n_points"] == 5


class TestApplyCalibration:
    """Tests for apply_calibration."""

    def setup_method(self):
        self.cal = fit_linear_calibration([0, 10, 20, 30], [100, 200, 300, 400])

    def test_known_values(self):
        result = apply_calibration([200, 300, 400], self.cal)
        np.testing.assert_allclose(result, [10, 20, 30], atol=1e-10)

    def test_below_detection_nan(self):
        result = apply_calibration([50], self.cal, below_detection="nan")
        assert np.isnan(result[0])

    def test_below_detection_zero(self):
        result = apply_calibration([50], self.cal, below_detection="zero")
        assert result[0] == 0.0

    def test_below_detection_negative(self):
        result = apply_calibration([50], self.cal, below_detection="negative")
        assert result[0] < 0


class TestDetectionLimits:
    """Tests for detection_limits."""

    def test_with_blank_responses(self):
        cal = fit_linear_calibration([0, 10, 20], [0, 100, 200])
        blanks = [2, 3, 1, 4, 2, 3]
        result = detection_limits(cal, blank_responses=blanks)

        assert result["lod"] > 0
        assert result["loq"] > result["lod"]
        assert result["loq"] / result["lod"] == pytest.approx(10.0 / 3.3, rel=0.01)

    def test_without_blanks(self):
        cal = fit_linear_calibration([0, 10, 20, 30], [1, 102, 199, 305])
        result = detection_limits(cal)

        assert result["lod"] > 0
        assert result["loq"] > result["lod"]


class TestDilutionCorrection:
    """Tests for dilution_corrected_concentration."""

    def test_simple_dilution(self):
        result = dilution_corrected_concentration([10, 20], 51)
        np.testing.assert_allclose(result, [510, 1020])

    def test_no_dilution(self):
        result = dilution_corrected_concentration([10], 1)
        np.testing.assert_allclose(result, [10])
