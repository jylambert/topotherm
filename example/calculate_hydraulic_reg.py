import pandas as pd

from topotherm import precalculation_hydraulic as precalc

SUPPLY = 90
RETURN = 55
AMBIENT = -20

# Optionally: calculate the supply temperatures for a given max supply t and
# return temperature according to ambient temperature
variable_feed_temp = precalc.determine_feed_line_temp(
    0,
    temp_sup_high=90,
    temp_sup_low=80,
    temp_turn_high=-14,
    temp_turn_low=6
    )
print(f'At {0}°C the supply temperature is {variable_feed_temp}°C')
print(variable_feed_temp)


# Calculate regression coefficients for the thermal capacity and heat losses
# thermal capacity regression and max power as well as part loads
r_thermal_cap = precalc.regression_thermal_capacity(t_supply=SUPPLY, t_return=RETURN)
# heat loss regression
r_heat_loss = precalc.regression_heat_losses(t_supply=SUPPLY, t_ambient=AMBIENT,
                                             thermal_capacity=r_thermal_cap)

# Save regression coefficients to csv file
df_cap = pd.DataFrame(r_thermal_cap['params'], index=[0])
print('Regression coefficients for thermal capacity:')
print(df_cap)

df_loss = pd.DataFrame(r_heat_loss['params'], index=[0])
print('Regression coefficients for heat losses:')
print(df_loss)
