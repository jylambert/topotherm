"""
This module contains functions for the thermo-hydraulic precalculation of the
district heating network.

Functions
---------
- ``determine_feed_line_temp`` : Calculates the supply temperature of the grid
  according to the outdoor temperature with a linear interpolation between the
  turning points.
- ``max_flow_velocity`` : Calculates the maximal flow velocity in the pipe.
- ``mass_flow`` : Calculates the maximal mass flow in the pipe.
- ``pipe_capacity`` : Calculates the maximal heat flow in the pipe.
- ``capacity_to_diameter`` : Calculates the diameter of the pipe according to
  the heat flow.
- ``thermal_resistance`` : Calculates the thermal resistance of the pipe.
- ``heat_loss_pipe`` : Calculates the heat loss of the pipe.
- ``regression_thermal_capacity`` : Calculates the regression factors for the
  linearization of the thermal capacity of the pipes.
- ``regression_heat_losses`` : Calculates the regression factors for the
  linearization of the heat losses of the pipes.
"""


import numpy as np
from scipy import stats
from scipy.optimize import root

from topotherm.settings import Settings


def determine_feed_line_temp(ambient_temperature: float,
                             temp_sup_high: float,
                             temp_sup_low: float,
                             temp_turn_high: float,
                             temp_turn_low: float) -> float:
    """
    Calculate the supply temperature of the grid according to the outdoor
    temperature using linear interpolation between the threshold points.

    Parameters
    ----------
    ambient_temperature : float
        Ambient temperature (°C).
    temp_sup_high : float
        Supply temperature high threshold (°C).
    temp_sup_low : float
        Supply temperature low threshold (°C).
    temp_turn_high : float
        Threshold for high temperature (°C).
    temp_turn_low : float
        Threshold for low temperature (°C).

    Returns
    -------
    float
        Supply temperature of the grid (°C).

    Examples
    --------
    >>> determine_feed_line_temp(0, 90, 60, 10, -10)
    75.0

    >>> determine_feed_line_temp(-15, 90, 60, 10, -10)
    90.0

    >>> determine_feed_line_temp(20, 90, 60, 10, -10)
    60.0
    """


    if ambient_temperature < temp_turn_high:
        temp_supply = temp_sup_high
    elif (ambient_temperature >= temp_turn_high) & (ambient_temperature <= temp_turn_low):
        temp_supply = temp_sup_high - ((ambient_temperature-temp_turn_high)/ \
                                       (temp_turn_low-temp_turn_high)) \
         * (temp_sup_high-temp_sup_low)
    else:
        temp_supply = temp_sup_low
    return temp_supply


def max_flow_velocity(vel_init: float,
                      diameter: float,
                      roughness: float,
                      max_spec_pressure_loss: float,
                      water_parameters: Settings) -> float:
    """
    Calculate the maximal flow velocity in a pipe based on the maximal
    specific pressure loss. All inputs must be non-negative real numbers.

    Parameters
    ----------
    vel_init : float
        Initial velocity (m/s).
    diameter : float
        Diameter of the pipe (m).
    roughness : float
        Roughness of the pipe (m).
    max_spec_pressure_loss : float
        Maximal specific pressure loss (Pa/m).
    water_parameters : Settings
        Water parameters instance.

    Returns
    -------
    float
        Maximal flow velocity (m/s).
    """

    def vel_calculation(var):
        """Function to calculate the maximal flow velocity in the pipe"""
        vel = var
        # Calculate Reynolds number
        re: float = ((water_parameters.density * vel.squeeze() * diameter)
                     / water_parameters.dynamic_viscosity)
        # Calculate friction factor f (from Haaland equation for turbulent flow)
        f = (-1.8 * np.log10((roughness / (3.7 * diameter)) ** 1.11 + 6.9 / re)) ** -2
        # Determine max. Velocity according to pressure loss and diameter
        eq = vel - np.sqrt((2 * max_spec_pressure_loss * diameter)
                           / (f * water_parameters.density))
        return eq

    sol = root(vel_calculation, vel_init, method='lm')
    if sol.success:
        vel_max = sol.x
    else:
        raise ValueError('Failed to calculate maximal flow velocity')

    return vel_max


def mass_flow(velocity: float,
              diameter: float,
              density_water: float) -> float:
    """
    Calculate the maximal mass flow in the pipe.

    Parameters
    ----------
    velocity : float
        Velocity in the pipe (m/s).
    diameter : float
        Diameter of the pipe (m).
    density_water : float
        Density of water (kg/m³).

    Returns
    -------
    float
        Maximal mass flow (kg/s).
    """
    mdot = density_water * velocity * (np.pi / 4) * diameter ** 2
    return mdot


def pipe_power(mass_flow: float,
               temperature_supply: float,
               temperature_return: float,
               cp_water: float) -> float:
    """
    Calculate the maximal heat flow in the pipe.

    Parameters
    ----------
    mass_flow : float
        Mass flow in the pipe (kg/s).
    temperature_supply : float
        Supply temperature (°C or K).
    temperature_return : float
        Return temperature (°C or K).
    cp_water : float
        Specific heat capacity of water (J/kg·K).

    Returns
    -------
    float
        Heat flow in the pipe (W).
    Examples
    --------
    >>> pipe_capacity(mass_flow=2.0,
    ...               temperature_supply=80,
    ...               temperature_return=60,
    ...               cp_water=4180)
    167200.0
    """
    p = mass_flow * cp_water * (temperature_supply - temperature_return)
    return p

#this function is it ok?
def capacity_to_diameter(thermal_power: float,
                         temperature_supply:float,
                         temperature_return: float,
                         cp_water: float) -> float:
    """
    Calculate the diameter of the pipe according to the heat flow.

    Parameters
    ----------
    thermal_power : float
        Heat power in the pipe (W).
    temperature_supply : float
        Supply temperature (°C or K).
    temperature_return : float
        Return temperature (°C or K).
    cp_water : float
        Specific heat capacity of water (J/kg·K).

    Returns
    -------
    float
        Diameter of the pipe (m).
    """
    d = thermal_power / (cp_water * (temperature_supply - temperature_return))
    return d


def thermal_resistance(diameter: float,
                       diameter_ratio: float,
                       depth: float,
                       settings: Settings) -> float:
    """
    Calculate the thermal resistance of a pipe.

    References
    ----------
    Planungshandbuch Fernwärme Version 1.3 (2021).

    Parameters
    ----------
    diameter : float
        Diameter of the pipe (m).
    diameter_ratio : float
        Outer diameter divided by inner diameter (-).
    depth : float
        Depth below ground level (m).
    settings : Settings
        Settings instance.

    Returns
    -------
    float
        Thermal resistance of the pipe (m·K/W).
    """
    outer_diameter = diameter * diameter_ratio
    thermal_resistance_ground = np.log(8 * depth / outer_diameter) / settings.ground.thermal_conductivity
    thermal_resistance_insulation = np.log(diameter_ratio) / settings.piping.thermal_conductivity
    thermal_resistance_pipe = (2 * np.pi) / (thermal_resistance_insulation + thermal_resistance_ground)
    return thermal_resistance_pipe


def heat_loss_pipe(temperature_in: float,
                   thermal_resistance_pipe: float,
                   ambient_temperature: float) -> float:
    """
    Calculate the heat loss of a pipe.

    Parameters
    ----------
    temperature_in : float
        Temperature in the pipe (°C or K).
    thermal_resistance_pipe : float
        Thermal resistance of the pipe (m·K/W).
    ambient_temperature : float
        Ambient temperature (°C or K).

    Returns
    -------
    float
        Heat loss of the pipe (W).
    """
    temp_diff_in = temperature_in - ambient_temperature
    losses = thermal_resistance_pipe * temp_diff_in

    return losses


def regression_thermal_capacity(settings: Settings) -> dict:
    """
    Calculate the regression factors for the linearization of the thermal
    capacity of the pipes. This is the first step in the thermo-hydraulic
    precalculation of the district heating network.

    Parameters
    ----------
    settings : RegressionSettings
        Regression settings instance.

    Returns
    -------
    dict
        Regression factors for the linearization (€/m).
    """
    V_INIT = 0.5  # initial velocity for hydraulic calculations

    r = {}  # results of the regression

    velocity_max = np.zeros(settings.piping.number_diameters)  # initialize array
    # iterate over all diameters
    for i, diam in enumerate(settings.piping.diameter):
        velocity_max[i] = max_flow_velocity(V_INIT,
                                            diam,
                                            settings.piping.roughness,
                                            settings.piping.max_pr_loss,
                                            settings.water).item()

    r['mass_flow_max'] = mass_flow(velocity_max,
                                   np.array(settings.piping.diameter),
                                   settings.water.density)

    r['power_flow_max'] = np.zeros([settings.piping.number_diameters])  # init

    # do the regression for each diameter
    r['power_flow_max'] = pipe_power(
        r['mass_flow_max'],
        settings.temperatures.supply,
        settings.temperatures.return_,
        settings.water.heat_capacity_cp)
    regression = stats.linregress(r['power_flow_max']/1000,
                                  np.array(settings.piping.cost))

    r['a'] = np.round(regression.slope, 6)
    r['b'] = np.round(regression.intercept, 3)
    r['r2'] = regression.rvalue**2

    # Determine maximal power flow  in kw
    r['power_flow_max_kW'] = np.round(r['power_flow_max']/1000, 3)

    # Part load according to outdoor temperature and feed line temperature
    # @TODO: refactor this
    r['power_flow_max_partload'] = 1
    return r


def regression_heat_losses(settings: Settings,
                           thermal_capacity: dict) -> dict:
    """
    Calculate the regression factors for the linearization of the heat losses
    of the pipes, based on the calculated thermal capacities
    (see ``regression_thermal_capacity``).

    Parameters
    ----------
    settings : Settings
        Settings instance.
    thermal_capacity : dict
        Thermal capacity regression factors.

    Returns
    -------
    dict
        Regression factors for the linearization of the heat losses of the pipes.
    """

    pipe_depth = np.ones(settings.piping.number_diameters) * settings.piping.depth
    ratio = np.array(settings.piping.outer_diameter) / np.array(settings.piping.middle_diameter)
    res_pipe = thermal_resistance(np.array(settings.piping.middle_diameter), ratio, pipe_depth, settings)

    heat_loss = heat_loss_pipe(settings.temperatures.supply,
                                res_pipe,
                                settings.temperatures.ambient)

    maximal_power = thermal_capacity['power_flow_max']

    regression = stats.linregress(maximal_power/1000, heat_loss/1000)

    r = {}

    r['b'] = np.round(regression.intercept, 6)
    r['a'] = np.round(regression.slope, 10)
    r['r2'] = regression.rvalue**2
    r['heat_loss'] = heat_loss

    return r
