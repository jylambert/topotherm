"""This example shows how to generated incidence matrices from geodataframes."""

from pathlib import Path

import geopandas as gpd
import pyomo.environ as pyo

import topotherm as tt

current_path = Path(__file__).parent.resolve()

# read the raw geodataframe
gdf_sinks = gpd.read_file(current_path / "district/sinks.gpkg")
gdf_roads = gpd.read_file(current_path / "district/roads.gpkg")
gdf_supply = gpd.read_file(current_path / "district/source.gpkg")

# these files need to be formatted in particular ways, especially the sinks one
print("Columns of sinks geodataframe:")
print(gdf_sinks.columns)

gdf_nodes, gdf_edges = tt.create_matrices.from_gdfs(
    sinks=gdf_sinks,
    roads=gdf_roads,
    sources=gdf_supply,
)

mat, gdf_nodes, gdf_edges = tt.create_matrices.create_matrices_from_gdf(gdf_nodes, gdf_edges)


a = 1

"""# this results in a coherent dataset containing nodes and edges, that can be
# import to topotherm, converted to networkx graphs, etc.
model_sets = tt.models.sets.create(incidence_matrices)
model = tt.models.single_timestep.create(
    matrices=incidence_matrices,
    sets=model_sets,
)

opt = pyo.SolverFactory("cbc")
opt.options["mipgap"] = 0.01
result = opt.solve(model, tee=True)


dummy = 1"""
