"""This example shows how to generated incidence matrices from geodataframes."""

from pathlib import Path

import os
import geopandas as gpd
import pyomo.environ as pyo

import topotherm as tt

current_path = Path(__file__).parent.resolve()

# read the raw geodataframe
gdf_sinks = os.path.join(current_path, "district", "sinks.gpkg")
gdf_roads = os.path.join(current_path, "district", "roads.gpkg")
gdf_supply = os.path.join(current_path, "district", "source.gpkg")


inputpaths = {"sinks": gdf_sinks,
              "roads": gdf_roads,
              "sources": gdf_supply}

outputpaths = Path(os.path.join(Path(__file__).parent.resolve(), "results"))

tt.create_matrices.from_gisfiles(
    inputpaths=inputpaths,
    outputpath=outputpaths,
    buffer=2.5,
    crs ="EPSG:25832"
)

# this results in a coherent dataset containing nodes and edges, that can be
# import to topotherm, converted to networkx graphs, etc.

mat = tt.fileio.load(Path(outputpaths))
model_sets = tt.models.sets.create(mat)

model = tt.models.single_timestep.create(
    matrices=mat,
    sets=model_sets,
)

opt = pyo.SolverFactory("gurobi")
opt.options["mipgap"] = 0.01
result = opt.solve(model, tee=True)
