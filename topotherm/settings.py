"""This file contains all the settings for the optimization problem and should
be modified and adapted to each case through a .yaml file (see examples)."""

import os
import warnings

from typing import List, Union

from pydantic import BaseModel, Field, model_validator
import yaml


class Water(BaseModel):
    """Water properties for the linearization of piping."""
    dynamic_viscosity: float = Field(
        4.041e-4, description="Dynamic viscosity at 70°C and 1 bar")
    density: float = Field(
        977.76, description="Density at 70°C and 1 bar")
    heat_capacity_cp: float = Field(
        4.187e3, description="Heat capacity at constant pressure")


class Ground(BaseModel):
    """Ground properties for the linearization of piping."""
    thermal_conductivity: float = Field(
        2.4, description="Thermal conductivity of the ground")


class Temperatures(BaseModel):
    """Temperatures for the linearization of piping, calculation of
    postprocessing."""
    ambient: float = Field(-20, description="Ambient temperature in °C")
    supply: float = Field(70, description="Supply temperature in °C")
    return_: float = Field(55, description="Return temperature in °C")


class Piping(BaseModel):
    """Piping properties for the linearization of piping."""
    diameter: List[float] = Field(
        default_factory=lambda: [
            0.0216, 0.0285, 0.0372, 0.0431,
            0.0545, 0.0703, 0.0825, 0.1071,
            0.1325, 0.1603, 0.2101, 0.263,
            0.3127, 0.3444, 0.3938
            ],
        description="List of all inner diameters of the available pipe sizes"
    )

    outer_diameter: List[float] = Field(
        default_factory=lambda: [
            0.09, 0.09, 0.11, 0.11,
            0.125, 0.14, 0.16, 0.2,
            0.225, 0.25, 0.315, 0.4,
            0.45, 0.5, 0.56
            ],
        description="List of all outer diameters of the available pipe sizes"
    )

    cost: List[float] = Field(
        default_factory=lambda: [
            390, 400, 430, 464,
            498, 537, 602, 670,
            754, 886, 1171, 1184,
            1197, 1401, 1755
            ],
        description="Cost of pipes"
    )

    number_diameters: int = Field(
        15, description="Number of discrete diameters")
    max_pr_loss: float = Field(
        250., description="Assumed pressure loss in Pa per meter")
    roughness: float = Field(
        0.05e-3, description="Pipe roughness factor")
    thermal_conductivity: float = Field(
        0.024, description="Pipe thermal conductivity in W/mK")

    @model_validator(mode='after')
    def check_length(self):
        """Check if the length of diameter, outer_diameter, and cost is
        consistent with the defined number diameters.
        """
        if len(self.diameter) != self.number_diameters:
            raise ValueError(
                f"""Length of diameter {len(self.diameter)} is not equal to
                number_diameters {self.number_diameters}""")
        if len(self.outer_diameter) != self.number_diameters:
            raise ValueError(
                f"""Length of outer_diameter {len(self.outer_diameter)} is
                not equal to number_diameters {self.number_diameters}""")
        if len(self.cost) != self.number_diameters:
            raise ValueError(
                f"""Length of cost {len(self.cost)} is not equal to
                number_diameters {self.number_diameters}""")
        return self


class Solver(BaseModel):
    """Solver properties for the optimization problem. Used for the
    optimization model."""
    mip_gap: float = Field(1e-4, description="MIP gap")
    time_limit: int = Field(
        10000, description="Time limit for the optimization in seconds")
    log: str = Field("solver.log", description="Log file for the solver")
    # @TODO: add more solver options, pass them to the solver flexibly
    


# @TODO: remove flh from setting and incorporate into fileio when modelling consumer specific flh
class Economics(BaseModel):
    """Economic properties for the optimization problem. Used for the
    optimization model."""
    source_price: List[List[float]] = Field(
        default_factory=lambda: [[80e-3]],
        description="Variable price for one kW of heat production at supply in €/kW")
    source_c_inv: List[float] = Field(
        default_factory=lambda: [0.],
        description="Investment costs for each source in €/kW")
    source_c_irr: List[float] = Field(
        default_factory=lambda: [0.08],
        description="Interest rate for sources")
    source_lifetime: List[float] = Field(
        default_factory=lambda: [20.],
        description="Lifetime for source investments in years")
    source_max_power: List[float] = Field(
        default_factory=lambda: [1e6],
        description="Maximum installed power for sources in kW")
    source_min_power: List[float] = Field(
        default_factory=lambda: [0],
        description="Maximum installed power for sources in kW")

    pipes_c_irr: float = Field(
        0.08, description="Interest rate for pipes")
    heat_price: float = Field(
        120e-3, description="Selling price for heat in €/kW")
    pipes_lifetime: float = Field(
        50.0, description="Lifetime for piping investments in years")


class Settings(BaseModel):
    """Class for the settings of the optimization problem which is passed
    to the regression, optimization model, and postprocessing."""
    water: Water
    ground: Ground
    temperatures: Temperatures
    piping: Piping
    solver: Solver
    economics: Economics

    def __init__(self,
                water: Water = Water(),
                ground: Ground = Ground(),
                temperatures: Temperatures = Temperatures(),
                piping: Piping = Piping(),
                solver: Solver = Solver(),
                economics: Economics = Economics()) -> None:
        super().__init__(water=water,
                        ground=ground,
                        temperatures=temperatures,
                        piping=piping, solver=solver,
                        economics=economics)


def load(file_path: Union[str, os.PathLike[str]]) -> Settings:
    """Load the settings from a yaml file."""
    with open(file_path, 'r', encoding='utf-8') as file:
        data = yaml.safe_load(file)
    # check if there are keys which are not in the model
    for key in data.keys():
        if key not in Settings.model_fields.keys():
            warnings.warn(f"Key {key} is not in the model.")
    return Settings(**data)
