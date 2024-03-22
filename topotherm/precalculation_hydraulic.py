""" This module contains the functions for the hydraulic precalculation of the
district heating network."""

import numpy as np
from scipy import stats
from scipy.optimize import fsolve

from topotherm.settings import Regression


def determine_feed_line_temp(ambient_temperature, temp_sup_high, temp_sup_low,
                             temp_turn_high, temp_turn_low):
    """Calculates the supply temperature of the grid according to the outdoor
    temperature with a linear interpolation between the turning points.

    Input: Ambient temperature (°C), supply temperature high (°C),
    supply temperature low (°C), turing point high (°c), turining point low
    (°C).
    Returns: supply temperature of grid (°C)"""
    if ambient_temperature < temp_turn_high:
        temp_supply = temp_sup_high
    elif (ambient_temperature >= temp_turn_high) & (ambient_temperature <= temp_turn_low):
        temp_supply = temp_sup_high - ((ambient_temperature-temp_turn_high)/ \
                                       (temp_turn_low-temp_turn_high)) \
        *(temp_sup_high-temp_sup_low)
    else:
        temp_supply = temp_sup_low
    return temp_supply


def max_flow_velocity(vel_init, diameter, roughness, max_spec_pressure_loss, water_parameters):
    """Input: initial_velocity (m/s), diameter (m), roughness(m),
    max_spec_pressure_loss (Pa/m)
    Returns: velocity_max (m/s)
    Inputs must be non-negative reals"""
    def vel_calculation(var):
        vel = var
        # Calculate Reynolds number
        Re: float = (water_parameters.density * vel.squeeze() * diameter) / water_parameters.dynamic_viscosity
        # Calculate friction factor f (from Haaland equation for turbulent flow)
        f = (-1.8 * np.log10((roughness / (3.7 * diameter)) ** 1.11 + 6.9 / Re)) ** -2
        # Determine max. Velocity according to pressure loss and diameter
        eq = vel - np.sqrt((2 * max_spec_pressure_loss * diameter) / (f * water_parameters.density))
        return eq

    vel_max = fsolve(vel_calculation, vel_init)
    return vel_max


def mass_flow(velocity, diameter, density_water):
    """Input: velocity (m/s), diameter(m)
    Returns: maximal mass flow in kg/s"""
    mass_flow = density_water * velocity * (np.pi / 4) * diameter ** 2
    return mass_flow


def pipe_capacity(mass_flow, temperature_supply, temperature_return, cp_water):
    """Input: mass_flow (kg/s), temperature_supply (°C or K) at a specific time
     step, temperature_return (°C or K)
    at a specific time step
    Returns: heat_flow_max (W)"""
    pipe_capacity = mass_flow * cp_water * (temperature_supply - temperature_return)
    return pipe_capacity


def capacity_to_diameter(pipe_capacity, temperature_supply, temperature_return, cp_water):
    mass_flow = pipe_capacity / (cp_water * (temperature_supply - temperature_return))
    return mass_flow


def thermal_resistance(diameter, diameter_ratio, depth, settings):
    """Input: diameter (m), diameter_ratio = outer diameter / diameter (-),
    depth (below ground level) (m)
    Returns: thermal_resistance_pipe (m*K/W)
    Formula: see Blommaert 2020 --> D'Eustachio 1957"""
    outer_diameter = diameter * diameter_ratio
    thermal_resistance_ground = np.log(4 * depth / outer_diameter) / (2 * np.pi * settings.ground.thermal_conductivity)
    thermal_resistance_insulation = np.log(diameter_ratio) / (2 * np.pi * settings.piping.thermal_conductivity)
    thermal_resistance_pipe = thermal_resistance_insulation + thermal_resistance_ground
    return thermal_resistance_pipe


def heat_loss_pipe(mass_flow, length, temperature_in, thermal_resistance_pipe, ambient_temperature, cp_water):
    """Input: mass_flow (kg/s), length (m), temperature_in (K or °C) at a
    specific time step, thermal_resistance_pipe (m*K/W),
    ambient_temperature (K or C°) at a specific time step

    Returns: Heat_loss_pipe (in W)"""
    temp_diff_in = temperature_in - ambient_temperature
    temp_diff_out = temp_diff_in * np.exp(
        -length / (
            mass_flow * cp_water * thermal_resistance_pipe
            )
        )
    temperature_out = temp_diff_out + ambient_temperature
    heat_loss_pipe = mass_flow * cp_water * (
        temperature_in - temperature_out) / length
    return heat_loss_pipe


def regression_thermal_capacity(settings: Regression):
    """
    Calculates the regression factors for the linearization of the thermal
    capacity of the pipes

    Input: initial velocity (m/s), inner diameter (m), roughness pipe (m),
    max specific pressure loss (Pa/m), supply temperature (°C), return
    temperature (°C)
    Returns: Regression factors for linearization (in €/m)

    Args:
        settings: Regression settings instance
    Returns:
        dict: Regression factors for linearization (in €/m)
    """
    V_INIT = 0.5  # initial velocity for hydraulic calculations

    r = {}  # results of the regression

    velocity_max = np.zeros(settings.piping.number_diameter)  # initialize array
    # iterate over all diameters
    for i, diam in enumerate(settings.piping.diameter):
        velocity_max[i] = max_flow_velocity(V_INIT, diam, settings.piping.roughness,
                                            settings.piping.max_pr_loss, settings.water).item()

    r['mass_flow_max'] = mass_flow(velocity_max, np.array(settings.piping.diameter), settings.water.density)

    r['power_flow_max'] = np.zeros([settings.piping.number_diameter])  # init

    # do the regression for each diameter
    r['power_flow_max'] = pipe_capacity(
        r['mass_flow_max'], settings.temperatures.supply, settings.temperatures.return_,
        settings.water.heat_capacity_cp)
    regression = stats.linregress(r['power_flow_max']/1000, np.array(settings.piping.cost))


    r['a'] = np.round(regression.slope, 6)
    r['b'] = np.round(regression.intercept, 3)
    r['r2'] = regression.rvalue**2

    # Determine maximal power flow  in kw
    r['power_flow_max_kW'] = np.round(r['power_flow_max']/1000, 3)

    # Part load according to outdoor temperature and feed line temperature
    # @TODO: refactor this
    r['power_flow_max_partload'] = 1
    return r


def regression_heat_losses(settings: Regression, thermal_capacity):
    """Input: mass_flow (kg/s), temperature_supply (K or °C), pipe_depth (m),
    pipe_length (m), ambient_temperature (K or °C) at a specific time step
    Returns: Heat_loss_pipe (in W/m) and regression factors for linearization
    (in kW/m)"""
    pipe_depth = np.ones(settings.piping.number_diameter)
    pipe_length = 100*np.ones(settings.piping.number_diameter)
    ratio = np.array(settings.piping.outer_diameter) / np.array(settings.piping.diameter)
    res_pipe = thermal_resistance(np.array(settings.piping.diameter), ratio, pipe_depth, settings)

    mass_flow = thermal_capacity['mass_flow_max']
    maximal_power = thermal_capacity['power_flow_max']

    heat_loss = np.zeros([settings.piping.number_diameter])

    heat_loss = heat_loss_pipe(mass_flow, pipe_length, settings.temperatures.supply, res_pipe,
                                settings.temperatures.ambient, settings.water.heat_capacity_cp)
    regression = stats.linregress(maximal_power/1000, heat_loss/1000)

    r = {}
    r['b'] = np.round(regression.intercept, 6)
    r['a'] = np.round(regression.slope, 10)
    r['r2'] = regression.rvalue**2
    r['heat_loss'] = heat_loss

    return r
