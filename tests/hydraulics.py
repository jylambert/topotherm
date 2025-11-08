"""Test the hydraulic precalculation module."""

from pytest import approx

from topotherm import hydraulic
import topotherm as tt


settings = tt.settings.Settings()
settings.temperatures.ambient = -20
settings.temperatures.supply = 90
settings.temperatures.return_ = 55



def test_feed_line_temp():
    """Test feed line temperature calculation"""
    # Optionally: calculate the supply temperatures for a given max supply t and
    # return temperature according to ambient temperature
    variable_feed_temp = hydraulic.determine_feed_line_temp(
        0,
        temp_sup_high=90,
        temp_sup_low=80,
        temp_turn_high=-14,
        temp_turn_low=6
        )
    assert variable_feed_temp == approx(83.0)


def test_regression():
    """Test thermal capacity and thermal losses regression"""
    r_thermal_cap = hydraulic.regression_thermal_capacity(settings)
    assert r_thermal_cap['a'] == approx(0.0173, rel=0.01)
    assert r_thermal_cap['b'] == approx(1153., rel=0.01)
    assert r_thermal_cap['r2'] == approx(0.9049, rel=0.01)
    # heat loss regression
    r_heat_loss = hydraulic.regression_heat_losses(
        settings, thermal_capacity=r_thermal_cap)

    # assert that each params item is within 1% of the expected value
    assert r_heat_loss['a'] == approx(4.348000e-07, rel=0.01)
    assert r_heat_loss['b'] == approx(0.02189, rel=0.01)
    assert r_heat_loss['r2'] == approx(0.658745, rel=0.01)
