""" This module contains the functions for the hydraulic precalculation of the
district heating network."""

import numpy as np
from scipy import stats
from scipy.optimize import fsolve
from fluids.friction import Colebrook
from fluids.core import Reynolds

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
        temp_supply = temp_sup_high - ((ambient_temperature-temp_turn_high) /
                                       (temp_turn_low-temp_turn_high)) \
            * (temp_sup_high-temp_sup_low)
    else:
        temp_supply = temp_sup_low
    return temp_supply


def max_flow_velocity(vel_init, diameter, abs_roughness, max_spec_pressure_loss, water_parameters):
    """Calculate the maximum flow velocity in a pipe using the fluids library.

    Parameters:
    vel_init (float): initial guess for the flow velocity in m/s
    diameter (float): internal pipe diameter in m
    abs_roughness (float): pipe roughness in m
    max_spec_pressure_loss (float): maximum allowable pressure drop in Pa/m
    water_parameters (dict): dictionary containing 'density' and 'dynamic_viscosity' for the fluid

    Returns:
    float: maximum flow velocity in the pipe in m/s
    """
    density = water_parameters['density']
    dynamic_viscosity = water_parameters['dynamic_viscosity']
    # Reasonable initial guess for velocity
    vel = vel_init  # m/s
    rel_roughness = abs_roughness / diameter

    # Set a tolerance for convergence of the iterative method
    tolerance = 1e-6

    while True:
        # Calculate Reynold's number
        Re = Reynolds(V=vel, D=diameter, rho=density, mu=dynamic_viscosity)
        # Calculate friction factor using the Colebrook equation
        friction_factor = Colebrook(Re=Re, eD=rel_roughness)

        # Calculate the new maximum possible velocity for the given pressure loss
        vel_new = np.sqrt((2 * max_spec_pressure_loss *
                          diameter) / (friction_factor * density))

        # Check for convergence
        if abs(vel - vel_new) < tolerance:
            break
        vel = vel_new

    return vel


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
    pipe_capacity = mass_flow * cp_water * \
        (temperature_supply - temperature_return)
    return pipe_capacity


def capacity_to_diameter(pipe_capacity, temperature_supply, temperature_return, cp_water):
    mass_flow = pipe_capacity / \
        (cp_water * (temperature_supply - temperature_return))
    return mass_flow


def thermal_resistance(diameter, diameter_ratio, depth, settings):
    """Thermal Resistance of Uninsulated Buried Pipe + Soil
    Input: diameter (m), diameter_ratio = outer diameter / diameter (-),
    depth (below ground level) (m)
    Returns: thermal_resistance_pipe (m*K/W)
    Formula: see Jianguang 2018 --> 2012 ASHRAE Handbook, or similarly 2020 ASHRAE Handbook, p. 236
    estimation error is less than 1%"""
    outer_radius = diameter * diameter_ratio / 2
    if depth / outer_radius > 4:
        thermal_resistance_soil = np.log(
            2 * depth / outer_radius) / (2 * np.pi * settings.piping.thermal_conductivity)
    elif depth / outer_radius > 2:
        thermal_resistance_soil = np.log(depth / outer_radius + (
            (depth / outer_radius) ** 2 - 1) ** 0.5) / (2 * np.pi * settings.piping.thermal_conductivity)
    else:  # Formula: see Blommaert 2020 --> D'Eustachio 1957"""
        thermal_resistance_soil = np.log(
            4 * depth / (outer_radius * 2)) / (2 * np.pi * settings.ground.thermal_conductivity)
    thermal_resistance_insulation = np.log(
        diameter_ratio) / (2 * np.pi * settings.piping.thermal_conductivity)
    thermal_resistance_total = thermal_resistance_insulation + thermal_resistance_soil
    return thermal_resistance_total


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
    # V_INIT = 0.5  # initial velocity for hydraulic calculations
    # Given diameters and their corresponding maximum velocities
    hyd_diameter = np.array(
        [0.0216, 0.0285, 0.0372, 0.0431, 0.0545, 0.0703, 0.0825, 0.1071, 0.1325, 0.1603, 0.2101, 0.263, 0.3127,
         0.3444, 0.3938])  # Considered Inner diameter in m
    velocity_max_merkblatt = np.array(
        # Corresponding max velocity in m/s
        [0.6, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 1.9, 2.0, 2.5, 3.3, 3.9, 5, 5, 5])
        # """max velocity Data from Österreichisches Kuratorium für Landtechnik und Landentwicklung (ÖKL), Merkblatt-Nr. 67 – Planung von Biomasseheizwerken und Nahwärmenetzen, Wien, 3. Auflage 2016 """
    r = {}  # results of the regression

    velocity_max = np.zeros(
        settings.piping.number_diameter)  # initialize array
    # iterate over all diameters
    for i, diam in enumerate(settings.piping.diameter):
        # Calculate initial velocity for each diameter by interpolation of 
        V_INIT = np.interp(settings.piping.diameter,
                           hyd_diameter, velocity_max_merkblatt)
        velocity_max[i] = max_flow_velocity(V_INIT, diam, settings.piping.roughness,
                                            settings.piping.max_pr_loss, settings.water).item()

    r['mass_flow_max'] = mass_flow(velocity_max, np.array(
        settings.piping.diameter), settings.water.density)

    r['power_flow_max'] = np.zeros([settings.piping.number_diameter])  # init

    # do the regression for each diameter
    r['power_flow_max'] = pipe_capacity(
        r['mass_flow_max'], settings.temperatures.supply, settings.temperatures.return_,
        settings.water.heat_capacity_cp)
    regression = stats.linregress(
        r['power_flow_max']/1000, np.array(settings.piping.cost))

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
    ratio = np.array(settings.piping.outer_diameter) / \
        np.array(settings.piping.diameter)
    res_pipe = thermal_resistance(
        np.array(settings.piping.diameter), ratio, pipe_depth, settings)

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
