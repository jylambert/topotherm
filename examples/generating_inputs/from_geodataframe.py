"""This example shows how to generated incidence matrices from geodataframes."""

from pathlib import Path

import pyomo.environ as pyo

import topotherm as tt

# read raw geodataframes
current_path = Path(__file__).parent.resolve()
gdf_sinks = current_path / "district" / "sinks.gpkg"
gdf_roads = current_path / "district" / "roads.gpkg"
gdf_supply = current_path / "district" / "source.gpkg"

# store in dict for import function
inputpaths = {"sinks": gdf_sinks, "roads": gdf_roads, "sources": gdf_supply}

outputpaths = current_path / "results"

tt.create_matrices.from_gisfiles(
    inputpaths=inputpaths, outputpath=outputpaths, buffer=2.5, crs="EPSG:25832"
)

# this results in a coherent dataset containing nodes and edges, that can be
# imported to topotherm, converted to networkx graphs, etc.
mat = tt.fileio.load(Path(outputpaths))
model_sets = tt.models.sets.create(mat)

model = tt.models.single_timestep.create(
    matrices=mat,
    sets=model_sets,
)

opt = pyo.SolverFactory("highs")
result = opt.solve(model, tee=True)

res = tt.models.single_timestep.postprocess(model, tt.settings.Settings())
print(res.keys())
