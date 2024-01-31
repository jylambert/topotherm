"""This file contains all the settings for the optimization problem and should be modified
and adapted to each case."""

import numpy as np


class Water:
    """Water properties"""
    # For 70°C and 1 bar
    dynamic_viscosity = 4.041e-4    # in Pa*s, INSERT VALUE HERE
    density = 977.76                # in kg/m**3, INSERT VALUE HERE
    heat_capacity_cp = 4.187e3      # in J/(kg*K), INSERT VALUE HERE


class Ground:
    """Ground properties"""
    thermal_conductivity = 2.4  # in W/(m*K), for saturated soil, INSERT VALUE HERE


class Temperatures:
    """Temperatures for the regression"""
    ambient = -20  # in °C
    supply = 90  # in °C
    return_ = 55  # in °C


class Piping:
    """Piping settings and properties"""
    # Considered pipe Inner diameter in m for regression
    diameter = np.array([0.0216, 0.0285, 0.0372, 0.0431, 0.0545, 0.0703, 0.0825, 0.1071, 0.1325, 0.1603, 0.2101, 0.263, 0.3127, 0.3444, 0.3938])
    # Considered Outer diameter in m (with insulation) for regression, same order as diameter
    outer_diameter = np.array([0.09, 0.09, 0.11, 0.11, 0.125, 0.14, 0.16, 0.2, 0.225, 0.25, 0.315, 0.4, 0.45, 0.5, 0.56])
    # Ratio between the inner and outer diameter
    ratio = outer_diameter/diameter
    # Specific investment costs for pipes in €/m. Same order as diameter
    cost = np.array([390, 400, 430, 464, 498, 537, 602, 670, 754, 886, 1171, 1184, 1197, 1401, 1755])
    # Maximal specific pressure loss in Pa/m
    max_pr_loss = 250
    # Roughness for pipes in m
    roughness = 0.05e-3
    # in W/(m*K), for polyurethane
    thermal_conductivity = 0.024
    number_diameter = diameter.shape[0]


class OptSettings:
    """Settings for the optimization"""
    mip_gap = 1e-4  # MIP gap
    time_limit = 10000  # Time limit for the optimization in seconds


class Economics:
    flh = 2500
    heat_price = 120 * 10**-3  # Selling Price for heat in €/kW
    source_price = 80 * 10**-3  # Price for heat in €/kW
    c_inv_source = np.array([0])  # Investment costs for each source, same order as sources in A_p
    life_time = 40  # Number of years for deprecation
    c_irr = 0.08  # Interest rate
