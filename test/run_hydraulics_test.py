from topotherm import precalculation_hydraulic as precalc

SUPPLY = 90
RETURN = 55
AMBIENT = -20


def test_feed_line_temp():
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
    """Test regression"""
    r_thermal_cap = precalc.regression_thermal_capacity(t_supply=SUPPLY,
                                                        t_return=RETURN)
    assert abs(r_thermal_cap['params']['a'] - 0.018377) < 0.01 * 0.018377
    assert abs(r_thermal_cap['params']['b'] - 567.335) < 0.01 * 567.335
    assert abs(r_thermal_cap['params']['r2'] - 0.869778) < 0.01 * 0.869778
    # heat loss regression
    r_heat_loss = precalc.regression_heat_losses(t_supply=SUPPLY,
                                                 t_ambient=AMBIENT,
                                                thermal_capacity=r_thermal_cap)
    # assert that each params item is within 1% of the expected value
    assert abs(r_heat_loss['params']['a'] - 4.348000e-07) < 0.01 * 4.348000e-07
    assert abs(r_heat_loss['params']['b'] - 0.02189) < 0.01 * 0.02189
    assert abs(r_heat_loss['params']['r2'] - 0.658745) < 0.01 * 0.658745
