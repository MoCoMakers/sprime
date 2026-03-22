"""
Tests for hill_fitting module using validated examples from 4PL.py.

These tests use examples with known expected results to verify
that the hill_fitting implementation produces correct outputs.
"""

import pytest

from sprime import HillCurveParams, hill_fitting


class TestExample1_MedCalc:
    """
    Test Example 1 from 4PL.py

    Source: https://www.medcalc.org/manual/nonlinearregression-example.php
    Expected output:
        A = 0.1536 (low-dose limit / zero_asymptote in sprime)
        B = 1.7718 (hill coefficient)
        C = 19.3494 (EC50)
        D = 28.4479 (upper asymptote)
    """

    @pytest.fixture
    def data(self):
        """Test data from Example 1"""
        concentrations = [0, 1.3, 2.8, 5, 10.2, 16.5, 21.3, 31.8, 52.2]
        responses = [0.1, 0.5, 0.9, 2.6, 7.1, 12.3, 15.3, 20.4, 24.2]
        return concentrations, responses

    @pytest.fixture
    def expected_params(self):
        """Expected parameters from Example 1"""
        return {
            "zero_asymptote": 0.1536,
            "steepness_coefficient": 1.7718,
            "ec50": 19.3494,
            "inf_asymptote": 28.4479,
        }

    def test_fit_hill_curve_example1(self, data, expected_params):
        """Test that fitting produces parameters close to expected values"""
        concentrations, responses = data

        result = hill_fitting.fit_hill_curve(concentrations, responses)

        # Verify result is HillCurveParams
        assert isinstance(result, HillCurveParams)

        # Check parameters with reasonable tolerance (0.1% relative or 0.01 absolute)
        assert abs(result.zero_asymptote - expected_params["zero_asymptote"]) < max(
            0.01, abs(expected_params["zero_asymptote"]) * 0.01
        )
        assert abs(result.steepness_coefficient - expected_params["steepness_coefficient"]) < max(
            0.01, abs(expected_params["steepness_coefficient"]) * 0.01
        )
        assert abs(result.ec50 - expected_params["ec50"]) < max(
            0.1, abs(expected_params["ec50"]) * 0.01
        )
        assert abs(result.inf_asymptote - expected_params["inf_asymptote"]) < max(
            0.01, abs(expected_params["inf_asymptote"]) * 0.01
        )

        # Verify r-squared is reasonable (should be high for good fit)
        assert result.r_squared is not None
        assert result.r_squared > 0.9  # Should have good fit

    def test_fit_hill_curve_example1_curve_direction_up(self, data, expected_params):
        """Test fitting with explicit 'up' curve direction"""
        concentrations, responses = data

        result = hill_fitting.fit_hill_curve(concentrations, responses, curve_direction="up")

        assert isinstance(result, HillCurveParams)
        # Should still produce reasonable fit
        assert result.r_squared is not None
        assert result.r_squared > 0.8

    def test_fit_hill_curve_example1_with_initial_guesses(self, data, expected_params):
        """Test fitting with initial guesses close to expected values"""
        concentrations, responses = data

        result = hill_fitting.fit_hill_curve(
            concentrations,
            responses,
            initial_zero_asymptote=0.15,
            initial_inf_asymptote=28.0,
            initial_ec50=20.0,
            initial_steepness_coefficient=1.8,
        )

        assert isinstance(result, HillCurveParams)
        # With good initial guesses, should converge to expected values
        assert abs(result.ec50 - expected_params["ec50"]) < 1.0


class TestExample2_Cardillo:
    """
    Test Example 2 from 4PL.py

    Source: Cardillo G. (2012) Four parameters logistic regression - There and back again
    https://it.mathworks.com/matlabcentral/fileexchange/38122

    Expected output:
        A = 0.001002 (low-dose limit / zero_asymptote in sprime)
        B = 1.515 (hill coefficient)
        C = 108 (EC50)
        D = 3.784 (high-dose limit / inf_asymptote in sprime)
        r-squared = 0.9998
        RMSE = 0.0200
    """

    @pytest.fixture
    def data(self):
        """Test data from Example 2"""
        concentrations = [0, 4.5, 10.6, 19.7, 40, 84, 210]
        responses = [0.0089, 0.0419, 0.0873, 0.2599, 0.7074, 1.528, 2.7739]
        return concentrations, responses

    @pytest.fixture
    def expected_params(self):
        """Expected parameters from Example 2"""
        return {
            "zero_asymptote": 0.001002,
            "steepness_coefficient": 1.515,
            "ec50": 108.0,
            "inf_asymptote": 3.784,
            "r_squared": 0.9998,
            "rmse": 0.0200,
        }

    def test_fit_hill_curve_example2(self, data, expected_params):
        """Test that fitting produces parameters close to expected values"""
        concentrations, responses = data

        # Use explicit 'up' direction since data is increasing
        result = hill_fitting.fit_hill_curve(concentrations, responses, curve_direction="up")

        # Verify result is HillCurveParams
        assert isinstance(result, HillCurveParams)

        # Check parameters with reasonable tolerance
        # Zero asymptote is very small, so use absolute tolerance
        assert abs(result.zero_asymptote - expected_params["zero_asymptote"]) < 0.01

        # Hill coefficient (allow for sign difference if direction was swapped)
        assert (
            abs(abs(result.steepness_coefficient) - expected_params["steepness_coefficient"]) < 0.1
        )

        # EC50 (larger value, so allow more tolerance)
        assert abs(result.ec50 - expected_params["ec50"]) < 5.0

        # Inf asymptote
        assert abs(result.inf_asymptote - expected_params["inf_asymptote"]) < 0.1

        # r-squared should be very high (close to 0.9998)
        assert result.r_squared is not None
        assert result.r_squared > 0.99  # Should be very close to 0.9998

    def test_fit_hill_curve_example2_high_quality_fit(self, data, expected_params):
        """Verify that Example 2 produces a high-quality fit"""
        concentrations, responses = data

        result = hill_fitting.fit_hill_curve(concentrations, responses)

        # This example is known to have excellent fit
        assert result.r_squared is not None
        assert result.r_squared > 0.999  # Should be very close to 0.9998


class TestHillFittingGeneral:
    """General tests for hill_fitting module"""

    def test_fit_hill_curve_basic(self):
        """Test basic functionality with simple data"""
        concentrations = [0.1, 1, 10, 100, 1000]
        responses = [5, 10, 50, 90, 95]

        result = hill_fitting.fit_hill_curve(concentrations, responses)

        assert isinstance(result, HillCurveParams)
        assert result.ec50 > 0
        assert result.inf_asymptote > result.zero_asymptote
        assert result.r_squared is not None
        assert 0 <= result.r_squared <= 1

    def test_fit_hill_curve_minimum_data_points(self):
        """Test that at least 4 data points are required"""
        concentrations = [1, 10, 100]
        responses = [10, 50, 90]

        with pytest.raises(ValueError, match="at least 4 data points"):
            hill_fitting.fit_hill_curve(concentrations, responses)

    def test_fit_hill_curve_mismatched_lengths(self):
        """Test that mismatched array lengths raise error"""
        concentrations = [1, 10, 100, 1000]
        responses = [10, 50, 90]  # Wrong length

        with pytest.raises(ValueError, match="same length"):
            hill_fitting.fit_hill_curve(concentrations, responses)

    def test_fit_hill_curve_with_bounds(self):
        """Test fitting with parameter bounds"""
        concentrations = [0.1, 1, 10, 100, 1000]
        responses = [5, 10, 50, 90, 95]

        # Constrain EC50 to be between 1 and 100
        # Use wider bounds that accommodate initial guesses (inf_asymptote can be up to 110)
        bounds = (
            [0, 0.5, 1, 0],  # Lower bounds: [zero_asymptote, hill, ec50, inf_asymptote]
            [100, 5.0, 100, 110],  # Upper bounds (allow inf_asymptote up to 110 for initial guess)
        )

        result = hill_fitting.fit_hill_curve(
            concentrations,
            responses,
            curve_direction="up",  # Specify direction to avoid auto-detection issues with bounds
            bounds=bounds,
        )

        assert isinstance(result, HillCurveParams)
        assert 1 <= result.ec50 <= 100

    def test_fit_hill_curve_with_custom_maxfev(self):
        """Test fitting with custom maxfev parameter"""
        concentrations = [0.1, 1, 10, 100, 1000]
        responses = [5, 10, 50, 90, 95]

        result = hill_fitting.fit_hill_curve(
            concentrations,
            responses,
            maxfev=10000,  # Lower than default
        )

        assert isinstance(result, HillCurveParams)
        assert result.r_squared is not None

    def test_fit_hill_curve_reverse_order(self):
        """Test that data in reverse order (high to low) is handled correctly"""
        concentrations = [1000, 100, 10, 1, 0.1]  # Reverse order
        responses = [95, 90, 50, 10, 5]

        result = hill_fitting.fit_hill_curve(concentrations, responses)

        assert isinstance(result, HillCurveParams)
        assert result.r_squared is not None

    def test_fit_hill_curve_zero_concentration(self):
        """Test handling of zero concentrations"""
        concentrations = [0, 1, 10, 100, 1000]
        responses = [5, 10, 50, 90, 95]

        result = hill_fitting.fit_hill_curve(concentrations, responses)

        assert isinstance(result, HillCurveParams)
        # Should handle zero without error
        assert result.ec50 > 0

    def test_fit_hill_curve_auto_direction_detection(self):
        """Test that auto-detection selects appropriate curve direction"""
        # Increasing curve
        concentrations = [0.1, 1, 10, 100, 1000]
        responses = [5, 10, 50, 90, 95]

        result = hill_fitting.fit_hill_curve(concentrations, responses)

        assert isinstance(result, HillCurveParams)
        assert result.inf_asymptote > result.zero_asymptote  # Should be increasing

        # Decreasing curve
        responses_decreasing = [95, 90, 50, 10, 5]

        result2 = hill_fitting.fit_hill_curve(concentrations, responses_decreasing)

        assert isinstance(result2, HillCurveParams)
        assert result2.inf_asymptote < result2.zero_asymptote  # Should be decreasing
