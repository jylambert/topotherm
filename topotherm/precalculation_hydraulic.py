# -*- coding: utf-8 -*-
import numpy as np
from scipy import stats
from scipy.optimize import fsolve

from topotherm.settings import Water, Ground, Piping, Temperatures


def determine_feed_line_temp(ambient_temperature, temp_sup_high, temp_sup_low, temp_turn_high, temp_turn_low):
    """Input: Ambient temperature (°C), supply temperature high (°C), supply temperature low (°C),
    turing point high (°c), turining point low (°C)
    Returns: supply temperature of grid (°C)"""
    if ambient_temperature < temp_turn_high:
        temp_supply = temp_sup_high
    elif (ambient_temperature >= temp_turn_high) & (ambient_temperature <= temp_turn_low):
        temp_supply = temp_sup_high - ((ambient_temperature-temp_turn_high)/(temp_turn_low-temp_turn_high))*(temp_sup_high-temp_sup_low)
    else:
        temp_supply = temp_sup_low
    return temp_supply


def init_temperatures():
    ts = {}
    # Set outdoor temperature
    ts['ambient'] = np.zeros([Piping.number_diameter]) + Temperatures.ambient
    # Determine feed line temperature according to outdoor temperature
    variable_feed_temp = determine_feed_line_temp(ts['ambient'][0], 90, 80, -14, 6)
    ts['supply'] =  variable_feed_temp * np.ones([Piping.number_diameter])
    # Set return temperature to 55 °C
    ts['return'] = np.ones([Piping.number_diameter]) * Temperatures._return
    return ts


def max_flow_velocity(vel_init, diameter, roughness, max_spec_pressure_loss):
    """Input: initial_velocity (m/s), diameter (m), roughness(m), max_spec_pressure_loss (Pa/m)
    Returns: velocity_max (m/s)
    Inputs must be non-negative reals"""
    def vel_calculation(var):
        vel = var
        reynolds: float = (Water.density * vel * diameter) / Water.dynamic_viscosity            # Calculate Reynolds number
        f = (-1.8 * np.log10((roughness / (3.7 * diameter)) ** 1.11 + 6.9 / reynolds)) ** -2    # Calculate friction factor f (from Haaland equation for turbulent flow)
        eq = vel - np.sqrt((2 * max_spec_pressure_loss * diameter) / (f * Water.density))       # Determine max. Velocity according to pressure loss and diameter
        return eq

    vel_max = fsolve(vel_calculation, vel_init)
    return vel_max


def mass_flow(velocity, diameter):
    """Input: velocity (m/s), diameter(m)
    Returns: mass_flow (kg/s)"""
    mass_flow = Water.density * velocity * (np.pi / 4) * diameter ** 2  # maximal mass flow in kg/s
    return mass_flow


def pipe_capacity(mass_flow, temperature_supply, temperature_return):
    """Input: mass_flow (kg/s), temperature_supply (°C or K) at a specific time step, temperature_return (°C or K)
    at a specific time step
    Returns: heat_flow_max (W)"""
    pipe_capacity = mass_flow * Water.heat_capacity_cp * (temperature_supply - temperature_return)
    return pipe_capacity


def capacity_to_diameter(pipe_capacity, temperature_supply, temperature_return):
    mass_flow = pipe_capacity / (Water.heat_capacity_cp * (temperature_supply - temperature_return))
    return mass_flow


def thermal_resistance(diameter, diameter_ratio, depth):
    """Input: diameter (m), diameter_ratio = outer diameter / diameter (-), depth (below ground level) (m)
    Returns: thermal_resistance_pipe (m*K/W)
    Formula: see Blommaert 2020 --> D'Eustachio 1957"""
    outer_diameter = diameter * diameter_ratio
    thermal_resistance_ground = np.log(4 * depth / outer_diameter) / (2 * np.pi * Ground.thermal_conductivity)
    thermal_resistance_insulation = np.log(diameter_ratio) / (2 * np.pi * Piping.thermal_conductivity)
    thermal_resistance_pipe = thermal_resistance_insulation + thermal_resistance_ground
    return thermal_resistance_pipe


def heat_loss_pipe(mass_flow, length, temperature_in, thermal_resistance_pipe, ambient_temperature):
    """Input: mass_flow (kg/s), length (m), temperature_in (K or °C) at a specific time step, thermal_resistance_pipe (m*K/W),
    ambient_temperature (K or C°) at a specific time step

    Returns: Heat_loss_pipe (in W)"""    
    temp_diff_in = temperature_in - ambient_temperature
    temp_diff_out = temp_diff_in * np.exp(-length / (mass_flow * Water.heat_capacity_cp * thermal_resistance_pipe))
    temperature_out = temp_diff_out + ambient_temperature
    heat_loss_pipe = mass_flow * Water.heat_capacity_cp * (temperature_in - temperature_out)/length
    return heat_loss_pipe


def regression_thermal_capacity(temperatures):
    """
    Calculates the regression factors for the linearization of the thermal capacity of the pipes

    
    Input: initial velocity (m/s), inner diameter (m), roughness pipe (m), max specific pressure loss (Pa/m),
    supply temperature (°C), return temperature (°C)
    Returns: Regression factors for linearization (in €/m)

    Args:
        temperatures (dict): dict containing the defined temperatures (°C) for supply, return and
        ambient

    Returns:
        _type_: Regression factors for linearization (in €/m)
    """
    V_INIT = 0.5  # initial velocity for hydraulic calculations

    temp_supply = temperatures['supply']
    temp_return = temperatures['return']
    # ambient_temp = temperatures['ambient']

    r = dict()  # results of the regression

    # np.zeros([Piping.number_diameter, ambient_temp.shape[1]])
    velocity_max = np.zeros(Piping.number_diameter)  # initialize array
    # itarate over all diameters
    for i, diam in enumerate(Piping.diameter):
        velocity_max[i] = max_flow_velocity(V_INIT, diam, Piping.roughness, Piping.max_pr_loss)

    r['mass_flow_max'] = mass_flow(velocity_max, Piping.diameter)

    r['power_flow_max'] = np.zeros([Piping.number_diameter])  # init

    # do the regression for each diameter
    r['power_flow_max'] = pipe_capacity(
        r['mass_flow_max'], temp_supply[0], temp_return[0])
    regression = stats.linregress(r['power_flow_max']/1000, Piping.cost)

    r['params'] = dict()
    r['params']['a'] = np.round(regression.slope, 6) 
    r['params']['b'] = np.round(regression.intercept, 3)
    r['params']['r2'] = regression.rvalue**2

    # Determine maximal power flow  in kw
    r['power_flow_max_kW'] = np.round(r['power_flow_max']/1000, 3)

    # Part load according to outdoor temperature and feed line temperature
    # @TODO: refactor this
    # r['power_flow_max_partload'] = r['power_flow_max_kW'][0, :] / r['power_flow_max_kW'][0, :].max()
    r['power_flow_max_partload'] = 1
    return r


def regression_heat_losses(temperatures, thermal_capacity):
    """Input: mass_flow (kg/s), temperature_supply (K or °C), pipe_depth (m), pipe_length (m)
    ambient_temperature (K or °C) at a specific time step
    Returns: Heat_loss_pipe (in W/m) and regression factors for linearization (in kW/m)"""
    pipe_depth = np.ones(Piping.number_diameter)
    pipe_length = 100*np.ones(Piping.number_diameter)
    res_pipe = thermal_resistance(Piping.diameter, Piping.ratio, pipe_depth)

    mass_flow = thermal_capacity['mass_flow_max']
    maximal_power = thermal_capacity['power_flow_max']

    temp_supply = temperatures['supply']
    ambient_temp = temperatures['ambient']

    heat_loss = np.zeros([Piping.number_diameter])
    # reg_factor = np.zeros([ambient_temp.shape[1], 3])
    factors = dict()

    # for col in range(ambient_temp.shape[1]):
    heat_loss = heat_loss_pipe(mass_flow, pipe_length, temp_supply, res_pipe,
                                ambient_temp)
    regression = stats.linregress(maximal_power/1000, heat_loss/1000)

    factors['b'] = np.round(regression.intercept, 6)
    factors['a'] = np.round(regression.slope, 10)
    factors['r2'] = regression.rvalue**2

    r = {}
    r['heat_loss'] = heat_loss
    r['params'] = factors
    return r
