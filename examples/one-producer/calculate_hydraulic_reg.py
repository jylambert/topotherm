import pandas as pd

from topotherm import precalculation_hydraulic as precalc
from topotherm.settings import Settings

settings = Settings()
settings.temperatures.ambient = -20
settings.temperatures.supply = 90
settings.temperatures.return_ = 55

# Optionally: calculate the supply temperatures for a given max supply t and
# return temperature according to ambient temperature
variable_feed_temp = precalc.determine_feed_line_temp(
    settings.temperatures.ambient,
    temp_sup_high=90,
    temp_sup_low=80,
    temp_turn_high=-14,
    temp_turn_low=6,
)
print(
    f"At {settings.temperatures.ambient} °C the supply temperature is {variable_feed_temp} °C"
)
print(variable_feed_temp)

# Calculate regression coefficients for the thermal capacity and heat losses
# thermal capacity regression and max power as well as part loads
r_thermal_cap = precalc.regression_thermal_capacity(settings)
# heat loss regression
r_heat_loss = precalc.regression_heat_losses(settings, thermal_capacity=r_thermal_cap)

# Save regression coefficients to csv file
df_cap = pd.DataFrame({key: r_thermal_cap[key] for key in ["a", "b", "r2"]}, index=[0])
print("Regression coefficients for thermal capacity:")
print(df_cap)

df_loss = pd.DataFrame({key: r_heat_loss[key] for key in ["a", "b", "r2"]}, index=[0])
print("Regression coefficients for heat losses:")
print(df_loss)
