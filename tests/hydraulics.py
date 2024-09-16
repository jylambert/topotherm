"""Test the hydraulic precalculation module."""
from topotherm import precalculation_hydraulic as precalc
import topotherm as tt

settings = tt.settings.Settings()
settings.temperatures.ambient = -20
settings.temperatures.supply = 90
settings.temperatures.return_ = 55


def test_feed_line_temp():
    """Test feed line temperature calculation"""
    # Optionally: calculate the supply temperatures for a given max supply t and
    # return temperature according to ambient temperature
    variable_feed_temp = precalc.determine_feed_line_temp(
        0,
        temp_sup_high=90,
        temp_sup_low=80,
        temp_turn_high=-14,
        temp_turn_low=6
        )
    assert variable_feed_temp == 83.0


def test_regression():
    """Test thermal capacity and thermal losses regression"""
    r_thermal_cap = precalc.regression_thermal_capacity(settings)
    assert abs(r_thermal_cap['a'] - 0.018377) < 0.01 * 0.018377
    assert abs(r_thermal_cap['b'] - 567.335) < 0.01 * 567.335
    assert abs(r_thermal_cap['r2'] - 0.869778) < 0.01 * 0.869778
    # heat loss regression
    r_heat_loss = precalc.regression_heat_losses(
        settings, thermal_capacity=r_thermal_cap)
    # assert that each params item is within 1% of the expected value
    assert abs(r_heat_loss['a'] - 4.348000e-07) < 0.01 * 4.348000e-07
    assert abs(r_heat_loss['b'] - 0.02189) < 0.01 * 0.02189
    assert abs(r_heat_loss['r2'] - 0.658745) < 0.01 * 0.658745
