"""This file contains all the settings for the optimization problem and should be modified
and adapted to each case."""

from dataclasses import dataclass

@dataclass
class Water:
    """Water properties for the linearization regression."""
    # For 70°C and 1 bar. this needs to be adapted to the actual temperature
    dynamic_viscosity: float = 4.041e-4
    density: float = 977.76
    heat_capacity_cp: float = 4.187e3

@dataclass
class Ground:
    """Ground properties."""
    thermal_conductivity: float = 2.4

@dataclass
class Temperatures:
    """Temperatures for the linearization regression."""
    ambient: float = -20
    supply: float = 70
    return_: float = 55

@dataclass
class Piping:
    """Piping properties for thermal losses and investment cost linearization regression."""
    # list of all available diameters
    # list of floats of all inner diameters of the available discrete pipe sizes
    diameter: tuple = (
        0.0216, 0.0285, 0.0372, 0.0431, 0.0545, 0.0703, 0.0825, 0.1071, 0.1325, 0.1603, 0.2101,
        0.263, 0.3127, 0.3444, 0.3938
        )
    # pipe outer diameter in m
    outer_diameter: tuple = (
        0.09, 0.09, 0.11, 0.11, 0.125, 0.14, 0.16, 0.2, 0.225, 0.25, 0.315, 0.4, 0.45, 0.5, 0.56)
    # pipe costs
    cost: tuple = (
        390, 400, 430, 464, 498, 537, 602, 670, 754, 886, 1171, 1184, 1197, 1401, 1755)
    max_pr_loss: int = 250  # assumed pressure loss in Pa per meter
    roughness: float = 0.05e-3  # pipe roughtness factor
    thermal_conductivity: float = 0.024  # pipe thermal conductivity in W/mK
    number_diameter: int = len(diameter)  # number of discrete diameters


@dataclass
class OptSettings:
    """Settings for the optimization problem."""
    mip_gap: float = 1e-4  # MIP gap
    time_limit: int = 10000  # Time limit for the optimization in seconds


@dataclass
class Economics:
    flh: int = 2500  # h/kWp and y
    heat_price: float = 120 * 10**-3  # Selling Price for heat in €/kW
    source_price: float = 80 * 10**-3  # Price for heat production at supply in €/kW
    c_inv_source: tuple = (0, )  # Investment costs for each source, same dim as sources in A_p
    life_time: int = 40  # Number of years for deprecation
    c_irr: float = 0.08  # Interest rate


@dataclass
class Optimization:
    """Settings with subclasses for the optimization problem."""
    opt_settings: OptSettings = OptSettings()
    economics: Economics = Economics()
    temperatures: Temperatures = Temperatures()

class Regression:
    """Settings for the linearization of the piping."""
    ground: Ground = Ground()
    water: Water = Water()
    temperatures: Temperatures = Temperatures()
    piping: Piping = Piping()