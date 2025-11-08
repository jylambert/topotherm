"""This file contains all the settings for the optimization problem and should
be modified and adapted to each case through a .yaml file (see examples)."""

import logging
import os
import warnings
from typing import List, Union

from pydantic import BaseModel, Field, field_validator, computed_field
import yaml

class Water(BaseModel):
    """Water properties for the linearization of piping."""
    dynamic_viscosity: float = Field(
        default=4.041e-4, gt=0,
        description="Dynamic viscosity at 70°C and 1 bar")
    density: float = Field(
        default=977.76, gt=0,
        description="Density at 70°C and 1 bar")
    heat_capacity_cp: float = Field(
        default=4.187e3, gt=0,
        description="Heat capacity at constant pressure")


class Ground(BaseModel):
    """Ground properties for the linearization of piping."""
    thermal_conductivity: float = Field(
        default=1.2, gt=0, description="Thermal conductivity of the ground")


class Temperatures(BaseModel):
    """Temperatures for the linearization of piping, calculation of
    postprocessing."""
    ambient: float = Field(20, description="Ambient temperature in °C")
    supply: float = Field(70, description="Supply temperature in °C")
    return_: float = Field(55, description="Return temperature in °C")


class PipeDiameter(BaseModel):
    """Single pipe diameter specification, with documentation."""
    dn: int = Field(description="DN size identifier as int: 32 for DN32")
    inner: float = Field(gt=0, description="Inner diameter in m")
    outer: float = Field(gt=0, description="Outer diameter of steel pipe in m")
    jacket: float = Field(gt=0, description="Outer diameter of the jacket pipe in m")
    cost: float = Field(ge=0, description="Cost per meter")
    
    @field_validator('outer')
    @classmethod
    def check_middle_diameter(cls, v, info):
        """Ensure outer steel pipe diameter is between inner and jacket."""
        if 'inner' in info.data and 'jacket' in info.data:
            _inner = info.data['innter']
            _outer = info.data['jacket']
            if not (_inner < v < _outer):
                raise ValueError(
                    f"Outer steel pipe diameter {v} must be between "
                    f"inner {_inner} and jacket {_outer}"
                )
        return v

class Piping(BaseModel):
    """Piping properties for the linearization of piping.
    Defaults according to DIN EN 253, costs based on a literature
    and expert survey."""
    diameters: List[PipeDiameter] = Field(
        default_factory=lambda: [
            PipeDiameter(dn=25, inner=0.0291, outer=0.0337, jacket=0.09, cost=718),
            PipeDiameter(dn=32, inner=0.0372, outer=0.0424, jacket=0.11, cost=763),
            PipeDiameter(dn=40, inner=0.0431, outer=0.0483, jacket=0.11, cost=786),
            PipeDiameter(dn=50, inner=0.0545, outer=0.0603, jacket=0.125, cost=880),
            PipeDiameter(dn=65, inner=0.0703, outer=0.0761, jacket=0.14, cost=907),
            PipeDiameter(dn=80, inner=0.0825, outer=0.0889, jacket=0.16, cost=1061),
            PipeDiameter(dn=100, inner=0.1071, outer=0.1143, jacket=0.2, cost=1090),
            PipeDiameter(dn=125, inner=0.1325, outer=0.1397, jacket=0.225, cost=1256),
            PipeDiameter(dn=150, inner=0.1603, outer=0.1683, jacket=0.25, cost=1332),
            PipeDiameter(dn=200, inner=0.2101, outer=0.2191, jacket=0.315, cost=1836),
            PipeDiameter(dn=250, inner=0.263, outer=0.273, jacket=0.4, cost=2036),
            PipeDiameter(dn=300, inner=0.3127, outer=0.3239, jacket=0.45, cost=2183),
            PipeDiameter(dn=350, inner=0.3444, outer=0.3556, jacket=0.50, cost=2651),
            PipeDiameter(dn=400, inner=0.3938, outer=0.4064, jacket=0.56, cost=2902),
            PipeDiameter(dn=450, inner=0.4444, outer=0.457, jacket=0.63, cost=3345),
            PipeDiameter(dn=500, inner=0.4954, outer=0.508, jacket=0.71, cost=3580),
            PipeDiameter(dn=600, inner=0.5958, outer=0.61, jacket=0.8, cost=4507),
        ],
        description="List of available pipe sizes, sorted by inner diameter"
    )

    max_pr_loss: float = Field(
        250., gt=0, description="Assumed pressure loss in Pa per meter")

    roughness: float = Field(
        0.01e-3, gt=0, description="Pipe roughness factor")

    thermal_conductivity: float = Field(
        0.024, gt=0, description="Pipe thermal conductivity in W/mK")

    depth: float = Field(
        2, gt=0.99, description="Depth of buried steel pipes in m (measured from surface to outer)")


    def sort_diameters(self):
        """Automatically sort diameters by inner diameter."""
        self.diameters = sorted(self.diameters, key=lambda p: p.inner)
        return self

    @computed_field
    @property
    def cost(self) -> list[float]:
        """Property that returns a list of cost per meter of each pipe
        diameter, sorted by diameter size
        """
        return self.get_diameter_lists()['cost']

    @computed_field
    @property
    def inner(self) -> list[float]:
        """Property that returns a list inner steel pipe diameters in m of each pipe,
        sorted by diameter size
        """
        return self.get_diameter_lists()['inner']

    @computed_field
    @property
    def outer(self) -> list[float]:
        """Property that returns a list outer steel pipes diameter in m of each pipe,
        sorted by diameter size
        """
        return self.get_diameter_lists()['outer']

    @computed_field
    @property
    def jacket(self) -> list[float]:
        """Property that returns a list outer jacket pipe diameter in m of each pipe,
        sorted by diameter size
        """
        return self.get_diameter_lists()['jacket']

    @computed_field
    @property
    def number_diameters(self) -> int:
        """Number of discrete diameters available."""
        return len(self.diameters)

    def get_by_dn(self, dn: int) -> PipeDiameter | None:
        """Get a pipe diameter by its DN identifier as int."""
        for pipe in self.diameters:
            if pipe.dn == dn:
                return pipe
        return None

    def get_diameter_lists(self) -> dict[str, List[float]]:
        """
        Get the diameter data as sorted lists (for backward compatibility).
        
        Returns:
            Dictionary with keys: 'inner', 'outer', 'jacket', 'cost'
        """
        return {
            'inner': [p.inner for p in self.diameters],
            'outer': [p.outer for p in self.diameters],
            'jacket': [p.jacket for p in self.diameters],
            'cost': [p.cost for p in self.diameters],
        }

                
    def set_from_dict(self, data: dict[int, dict[str, float]]) -> None:
        """Set costs, inner, outer, jacket from a dictionary with DN as keys.
        If a DN is not found, it is added. If it exists, the corresponding
        property is updated and if not all are defined, the previous value is kept."""
        diameters = []
        for dn, props in data.items():
            if dn in self.dns:
                pipe = self.get_by_dn(dn)
                if pipe is not None:
                    pipe.inner = props.get('inner', pipe.inner)
                    pipe.outer = props.get('outer', pipe.outer)
                    pipe.jacket = props.get('jacket', pipe.jacket)
                    pipe.cost = props.get('cost', pipe.cost)
                    diameters.append(pipe)
                    logging.debug("Updated existing pipe DN%i with properties %s", dn, props)
            else:
                required = ['inner', 'outer', 'jacket', 'cost']
                missing = [k for k in required if k not in props]
                if missing:
                    raise ValueError(f"New DN {dn} missing required properties: {missing}")
                diameters.append(
                    PipeDiameter(
                        dn=dn,
                        inner=props['inner'],
                        outer=props['outer'],
                        jacket=props['jacket'],
                        cost=props['cost']
                    )
                )
        self.diameters = sorted(diameters, key=lambda p: p.inner)
 
    @property
    def dns(self) -> List[str | int]:
        """Get list of DN identifiers in sorted order."""
        return [p.dn for p in self.diameters]


class Solver(BaseModel):
    """Solver properties for the optimization problem. Used for the
    optimization model."""
    mip_gap: float = Field(1e-4, ge=0, description="MIP gap")
    time_limit: int = Field(
        10000, gt=0, description="Time limit for the optimization in seconds")
    log: str = Field("solver.log", description="Log file for the solver")
    # @TODO: add more solver options, pass them to the solver flexibly


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
