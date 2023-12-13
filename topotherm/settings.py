# -*- coding: utf-8 -*-
import numpy as np


class Water:
    # For 70°C and 1 bar
    dynamic_viscosity = 4.041e-4    # in Pa*s, INSERT VALUE HERE
    density = 977.76                # in kg/m**3, INSERT VALUE HERE
    heat_capacity_cp = 4.187e3      # in J/(kg*K), INSERT VALUE HERE


class Ground:
    thermal_conductivity = 2.4      # in W/(m*K), for saturated soil, INSERT VALUE HERE


class Piping:
    diameter = np.array([0.0216, 0.0285, 0.0372, 0.0431, 0.0545, 0.0703, 0.0825, 0.1071, 0.1325, 0.1603, 0.2101, 0.263, 0.3127, 0.3444, 0.3938])   # Considered Inner diameter in m
    outer_diameter = np.array([0.09, 0.09, 0.11, 0.11, 0.125, 0.14, 0.16, 0.2, 0.225, 0.25, 0.315, 0.4, 0.45, 0.5, 0.56])     # Considered Outer diameter in m (with insulation)
    ratio = outer_diameter/diameter  # Ratio between the inner and outer diameter
    cost = np.array([390, 400, 430, 464, 498, 537, 602, 670, 754, 886, 1171, 1184, 1197, 1401, 1755])   # Specific investment costs for pipes in €/m
    max_pr_loss = 250 # Maximal specific pressure loss in Pa/m
    roughness = 0.05e-3 # Roughness for pipes in m
    thermal_conductivity = 0.024  # in W/(m*K), for polyurethane, INSERT VALUE HERE
    number_diameter = diameter.shape[0]


class OptSettings:
    mip_gap = 1e-2
    time_limit = 10000
