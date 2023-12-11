# -*- coding: utf-8 -*-
import numpy as np
from scipy import stats
from scipy.optimize import fsolve
from src.settings import Water, Ground, Piping


def determine_feed_line_temp(ambient_temperature, temp_sup_high, temp_sup_low, temp_turn_high, temp_turn_low):
    """Input: Ambient temperature (°C), supply temperature high (°C), supply temperature low (°C),
    turing point high (°c), turining point low (°C)
    Returns: supply temperature of grid (°C)"""
    temp_supply = np.zeros(ambient_temperature.shape[0])
    for row in range(ambient_temperature.shape[0]):
        if ambient_temperature[row] < temp_turn_high:
            temp_supply[row] = temp_sup_high
        elif (ambient_temperature[row] >= temp_turn_high) & (ambient_temperature[row] <= temp_turn_low):
            temp_supply[row] = temp_sup_high - ((ambient_temperature[row]-temp_turn_high)/(temp_turn_low-temp_turn_high))*(temp_sup_high-temp_sup_low)
        else:
            temp_supply[row] = temp_sup_low
    return temp_supply


def calc_max_flow_velocity(vel_init, diameter, roughness, max_spec_pressure_loss):
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


def calc_mass_flow(velocity, diameter):
    """Input: velocity (m/s), diameter(m)
    Returns: mass_flow (kg/s)"""
    mass_flow = Water.density * velocity * (np.pi / 4) * diameter ** 2  # maximal mass flow in kg/s
    return mass_flow


def calc_pipe_capacity(mass_flow, temperature_supply, temperature_return):
    """Input: mass_flow (kg/s), temperature_supply (°C or K) at a specific time step, temperature_return (°C or K)
    at a specific time step
    Returns: heat_flow_max (W)"""
    pipe_capacity = mass_flow * Water.heat_capacity_cp * (temperature_supply - temperature_return)
    return pipe_capacity


def capacity_to_diameter(pipe_capacity, temperature_supply, temperature_return):
    mass_flow = pipe_capacity / (Water.heat_capacity_cp * (temperature_supply - temperature_return))
    return mass_flow


def calc_thermal_resistance(diameter, diameter_ratio, depth):
    """Input: diameter (m), diameter_ratio = outer diameter / diameter (-), depth (below ground level) (m)
    Returns: thermal_resistance_pipe (m*K/W)
    Formula: see Blommaert 2020 --> D'Eustachio 1957"""
    outer_diameter = diameter * diameter_ratio
    thermal_resistance_ground = np.log(4 * depth / outer_diameter) / (2 * np.pi * Ground.thermal_conductivity)
    thermal_resistance_insulation = np.log(diameter_ratio) / (2 * np.pi * Piping.thermal_conductivity)
    thermal_resistance_pipe = thermal_resistance_insulation + thermal_resistance_ground
    return thermal_resistance_pipe


def calc_heat_loss_pipe(mass_flow, length, temperature_in, thermal_resistance_pipe, ambient_temperature):
    """Input: mass_flow (kg/s), length (m), temperature_in (K or °C) at a specific time step, thermal_resistance_pipe (m*K/W),
    ambient_temperature (K or C°) at a specific time step
    Returns: Heat_loss_pipe (in W)"""
    temp_diff_in = temperature_in - ambient_temperature
    temp_diff_out = temp_diff_in * np.exp(-length / (mass_flow * Water.heat_capacity_cp * thermal_resistance_pipe))
    temperature_out = temp_diff_out + ambient_temperature
    heat_loss_pipe = mass_flow * Water.heat_capacity_cp * (temperature_in - temperature_out)/length
    return heat_loss_pipe


def calc_regression_thermal_capacity(velocity_init, diameter, pipe_roughness, max_pr_loss, temp_supply, temp_return, ambient_temp):
    """Input: initial velocity (m/s), inner diameter (m), roughness pipe (m), max specific pressure loss (Pa/m),
    supply temperature (°C), return temperature (°C)
    Returns: Regression factors for linearization (in €/m)"""
    np.zeros([Piping.number_diameter, ambient_temp.shape[1]])
    velocity_max = np.zeros(Piping.number_diameter)
    for i in range(Piping.number_diameter):
        velocity_max[i] = calc_max_flow_velocity(velocity_init, diameter[i], pipe_roughness, max_pr_loss)
    mass_flow_max = calc_mass_flow(velocity_max, Piping.diameter)
    power_flow_max = np.zeros([Piping.number_diameter, ambient_temp.shape[1]])
    power_flow_max_regression = np.zeros([ambient_temp.shape[1], 3])
    for col in range(ambient_temp.shape[1]):
        power_flow_max[:, col] = calc_pipe_capacity(mass_flow_max, temp_supply[:, col], temp_return[:, col])
        regression = stats.linregress(power_flow_max[:, col]/1000, Piping.cost)
        power_flow_max_regression[col, 0] = regression.intercept
        power_flow_max_regression[col, 1] = regression.slope
        power_flow_max_regression[col, 2] = regression.rvalue**2
    return mass_flow_max, power_flow_max, power_flow_max_regression


def calc_regression_heat_losses(mass_flow, temperature_supply, pipe_depth, pipe_length, ambient_temp, maximal_power):
    """Input: mass_flow (kg/s), temperature_supply (K or °C), pipe_depth (m), pipe_length (m)
    ambient_temperature (K or °C) at a specific time step
    Returns: Heat_loss_pipe (in W/m) and regression factors for linearization (in kW/m)"""
    res_pipe = calc_thermal_resistance(Piping.diameter, Piping.ratio, pipe_depth)
    heat_loss = np.zeros([Piping.number_diameter, ambient_temp.shape[1]])
    reg_factor = np.zeros([ambient_temp.shape[1], 3])
    for col in range(ambient_temp.shape[1]):
        heat_loss[:, col] = calc_heat_loss_pipe(mass_flow, pipe_length, temperature_supply[:, col], res_pipe, ambient_temp[:, col])
        regression = stats.linregress(maximal_power[:, col]/1000, heat_loss[:, col]/1000)
        reg_factor[col, 0] = regression.intercept
        reg_factor[col, 1] = regression.slope
        reg_factor[col, 2] = regression.rvalue**2
    return heat_loss, reg_factor


"""
cmap = ListedColormap(sns.color_palette("deep").as_hex())
color_list = [*cmap(np.linspace(0, 1, 10))]
global_tech_colors = [color_list[8], color_list[9], color_list[0], color_list[1], color_list[7], color_list[3]]

num_points = 100
x2 = np.linspace(0, p_max_kw[-1], num_points)
fig6, ax6 = plt.subplots()

ax6.set_xlabel('Thermal capacity in MW')
ax6.set_ylabel('Heat Losses in W/m')
ax6.set_xlim([-0.9, 70])
ax6.set_ylim([0, 55])
ax6.scatter(p_max_kw/1000, heat_loss_power_flow, marker=".", label="Real inner diameters", color =global_tech_colors[5])
ax6.plot(x2/1000, (regression_heat_loss[0][0] + regression_heat_loss[0][1]*x2)*1000, label="Fitted inner diameters", color=global_tech_colors[4])
ax6.legend(ncol=1, loc='lower right', fancybox=False, shadow=False, borderpad=1, fontsize=10, frameon=False)
fig6.set_size_inches(8.5/2, 3.5)
plt.rcParams.update({'figure.autolayout': True})
plt.savefig('linearization_paper_heat_loss.svg')
plt.close()

fig7, ax7 = plt.subplots()
ax7.set_xlabel('Thermal capacity in MW')
ax7.set_ylabel('Investment Costs in €/m ')
ax7.set_xlim([-0.9, 70])
ax7.set_ylim([0, 1850])
ax7.scatter(p_max_kw/1000, Piping.cost, marker=".", label="Real inner diameters", color= global_tech_colors[5])
ax7.plot(x2/1000, regression_capacity[0][0] + regression_capacity[0][1]*x2, label="Fitted inner diameters", color=global_tech_colors[4])
ax7.legend(ncol=1, loc='lower right', fancybox=False, shadow=False, borderpad=1, fontsize=10, frameon=False)
fig7.set_size_inches(8.5/2, 3.5)
plt.rcParams.update({'figure.autolayout': True})
plt.savefig('linearization_paper_invest.svg')
plt.close()
"""