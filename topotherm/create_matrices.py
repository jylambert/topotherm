"""Module to create incidence matrices and geo-spatial data for district heating networks
from files or dataframes."""

import logging
from pathlib import Path
from typing import Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from shapely.geometry import LineString, MultiPoint
from shapely.ops import nearest_points, snap, split

from topotherm.utils import create_dir, find_duplicate_cols

# # required information of nodes and edges input
# _COLUMNS_NODES_DF = [
#     "x",
#     "y",  # coordinates of each nodes. no overlaps
#     "nodetype",  # either "sink", "source", "internal"
#     "ts",  # at least one column containing the power demand of each sink for each time step. Example: ts_0 = 10
#     "flh",  # at least one column with the full load hours in kWh/kWp of each sink
#     # "existing"
# ]

# _COLUMN_EDGES_DF = [
#     "x_start",
#     "y_start",  # starting coordinates of the edge
#     "x_end",
#     "y_end",  # ending coordinates. crs needs to be equal to nodes
# ]

# _COLUMNS_GDF = [  # TODO implement
#     "ts_",  # columns with power demand in kW of each sink for each time step. Example: ts_0 = 10
#     "flh_",  # columns with the full load hours in kWh/kWp of each sink at time step. Example: flh_0 = 2000
# ]


def create_connection_line(point, edges):
    """
    Creates a connection line from a single point to multiple edges

    Parameters
    ----------
    point :
        Single point which should be connected
    edges :
        Geoseries of edges to which the point should be connected

    Returns
    -------
    line :
        A single linestring connecting the point to the nearest edge
    """

    # Find the nearest point on the edges
    nearest_geom = edges.distance(point).idxmin()
    nearest_point = nearest_points(point, edges.iloc[nearest_geom])[1]
    # Create a line connecting the point to the nearest point on the edge
    return LineString([point, nearest_point])


def create_nearest_point(point, edges):
    """
    Finds the nearest point on a geoseries of edges

    Parameters
    ----------
    point :
        Single point which should be projected onto an edge
    edges :
        Geoseries of edges onto which the point should be projected

    Returns
    -------
    nearest_point :
        A single nearest point object.
    """

    # Find the nearest point on the edges
    nearest_geom = edges.distance(point).idxmin()
    nearest_point = nearest_points(point, edges.iloc[nearest_geom])[1]
    return nearest_point


def create_edge_nearest_point_optimized(points1, points2, crs="EPSG:25832"):
    """
    Creates an edge for each points1 to the nearest points2

    Parameters
    ----------
    points1 :
        Series of points which should be connected to the nearest points2
    points2 :
        Series of points to which should be connected the points1

    Returns
    -------
    lines_gs :
        Geoseries of lines connecting points1 to points2
    """

    # Extract coordinates from GeoSeries
    coords1 = np.array(list(zip(points1.x, points1.y)))
    coords2 = np.array(list(zip(points2.x, points2.y)))

    # Use cKDTree for efficient nearest-neighbor search
    tree = cKDTree(coords2)
    _, indices = tree.query(coords1)

    # Create LineStrings connecting each point in points1 to the nearest point in points2
    lines = [
        LineString([points1.iloc[i], points2.iloc[indices[i]]])
        for i in range(len(points1))
    ]

    return gpd.GeoSeries(lines, crs=points1.crs if hasattr(points1, "crs") else crs)


def _validate_sinks(gdf: gpd.GeoDataFrame) -> None:
    """
    Validates if the provided GeoDataFrame contains the required columns.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        The GeoDataFrame to validate.

    Raises
    ------
    ValueError
        If any required columns are missing or invalid.
    """
    gdf_columns = set(gdf.columns)

    # Identify ts_ and flh_ columns which are obligatory
    ts_columns = [col for col in gdf_columns if col.startswith("ts_")]
    flh_columns = [col for col in gdf_columns if col.startswith("flh_")]

    # Defining only one of each is sufficient for single time step models
    if not ts_columns and "ts" not in gdf_columns:
        raise ValueError(
            "GeoDataFrame must contain at least one 'ts_' column or a 'ts' column."
        )
    if not flh_columns and "flh" not in gdf_columns:
        raise ValueError(
            "GeoDataFrame must contain at least one 'flh_' column or a 'flh' column."
        )

    if ts_columns or flh_columns:
        ts_steps = _extract_timesteps(ts_columns, "ts_")
        flh_steps = _extract_timesteps(flh_columns, "flh_")

        if ts_steps != flh_steps:
            raise ValueError(
                "The 'ts_' and 'flh_' columns must contain the same sequence of time steps."
            )

        expected_steps = list(range(len(ts_steps)))
        if ts_steps != expected_steps:
            raise ValueError(
                "Time step columns must start at 0 and increase by 1 without gaps."
            )

    if not all(gdf.geometry.type == "Point"):
        raise ValueError("All geometries in the GeoDataFrame must be of type 'Point'.")


def _extract_timesteps(columns: list[str], prefix: str) -> list[int]:
    """Extract and sort timestep numbers from column names."""
    timesteps = []
    for col in columns:
        try:
            timestep = int(col.split("_")[1])
            timesteps.append(timestep)
        except (IndexError, ValueError):
            logging.warning(
                f"Column '{col}' has invalid format. Expected: '{prefix}<timestep>'."
            )
    return sorted(timesteps)


def from_gisfiles(
    inputpaths: dict[str, str | Path],
    outputpath: str | Path,
    buffer: float = 2.5,
    crs: str | int = "EPSG:25832",
):
    """
    Creates from a set of GIS file, such as shapefiles or gepackage files.
    Needs three files as input: roads, heat sources and heat sinks.

    Parameters
    ----------
    inputpaths: dict[str, os.PathLike]
        dict of paths to the geodata files for the street layout, heat sinks and
        heat sources. Each file needs to be formatted correctly. Example:
        {"roads": "path/to/roads.shp",
        "sinks": "path/to/sinks.shp",
        "sources": "path/to/sources.gpkg"}
    outputpath: str | pathlib.Path
        Path to the result folder
    buffer: float
        Radius of the buffer layer in meter, which is used to aggregate connection lines
    crs: str | int
        Coordinate reference system to which all files are projected. Default is
        "EPSG:25832" (ETRS89 / UTM zone 32N)

    Returns
    -------
    None
        saves parquet files to the outputpath needed for the optimization of the
        district heating network and nodes and edges files.
    """
    create_dir(outputpath)

    logging.info("Loading shapefiles...")
    if not all(key in inputpaths for key in ["sinks", "roads", "sources"]):
        raise ValueError(
            f"inputpaths keys {list(inputpaths.keys())} must contain 'sinks', 'roads', and 'sources'."
        )
    sinks = gpd.read_file(inputpaths["sinks"]).to_crs(crs)
    roads = gpd.read_file(inputpaths["roads"]).to_crs(crs)
    sources = gpd.read_file(inputpaths["sources"]).to_crs(crs)

    _validate_sinks(sinks)

    logging.info("Processing geodata...")
    mat, gdf_nodes, gdf_roads = from_gdfs(sinks, roads, sources, buffer)

    logging.info("Saving results...")
    gdf_roads.to_file(outputpath + "edges.shp")
    gdf_nodes.to_file(outputpath + "nodes.shp")

    for key, val in mat.items():
        pd.DataFrame(val).to_parquet(outputpath / f"{key}.parquet")
    return


def from_dfs(
    nodes: pd.DataFrame, edges: pd.DataFrame
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Converts a pd.DataFrame containting nodes and one containing edges
    to two geodataframes that can be converted to input matrices. The matrices
    are a requirement to run the optimization models.

    If edges are not connected to the sinks and sources, they will be connected
    by shortest path distance.

    Parameters
    ----------
    nodes : pd.DataFrame
        Nodes df containing information on sinks, sources and junctions
    edges : pd.DataFrame
        Edges df containing information on how all nodes are
        connected to eachother.

    Returns
    -------
    tuple of gpd.GeoDataFrame
        converted nodes and edges as gdf.
    """

    # TODO: write script to convert from df to gdf
    return (gdf_nodes, gdf_edges)


def from_gdfs(
    sinks: gpd.GeoDataFrame,
    roads: gpd.GeoDataFrame,
    sources: gpd.GeoDataFrame,
    buffer: float = 2.5,
    crs: str | int = "EPSG:25832",
):
    """Process geodata of sinks, roads, and sources to create nodes and edges,
    as well as incidence matrices for district heating network optimization.
    Connects sources, sinks to the road network by shortest path, merges nearby
    nodes, and splits roads at connection points.

    Parameters
    ----------
    sinks : gpd.GeoDataFrame
        GeoDataFrame containing heat sink locations. Needs to have columns for
        time series data (e.g., ``ts_0``, ``ts_1``, ...) of heating demand in kW
        and full load hours (e.g., ``flh_0``, ``flh_1``, ...) in h/year.
    roads : gpd.GeoDataFrame
        GeoDataFrame containing road network.
    sources : gpd.GeoDataFrame
        GeoDataFrame containing heat source locations.
    buffer : float
        Buffer radius in meters for aggregating connection lines.
    crs : str | int
        Coordinate reference system for all geodata.

    Returns
    -------
    gdf_nodes : gpd.GeoDataFrame
        Node information.
    gdf_edges : gpd.GeoDataFrame
        Edge information.
    """
    logging.info("Processing road intersections...")
    # Optimize road intersections processing
    road_intersections = roads.boundary.explode(
        ignore_index=False, index_parts=False
    ).drop_duplicates()

    logging.info("Processing house centers...")
    # Create centroids more efficiently
    sinks_center = sinks.copy()
    sinks_center.geometry = sinks.centroid

    logging.info("Creating union...")
    # Optimize union creation
    union_src_snk = pd.concat(
        [sinks_center[["geometry"]], sources[["geometry"]]], ignore_index=True
    )

    logging.info("Creating connection lines...")
    # Create connection lines between sources and sinks to the roadnetwork
    connection_line = union_src_snk.geometry.apply(
        lambda point: create_connection_line(point, roads["geometry"])
    )
    connection_nodes = connection_line.boundary.explode(
        ignore_index=False, index_parts=False
    )

    # Optimize matching using spatial operations
    matches = connection_nodes.isin(sinks_center.geometry)
    connection_nodes = connection_nodes[~matches].drop_duplicates()

    logging.info("Merging nearby nodes...")
    # Optimize node merging
    merged = connection_nodes.buffer(buffer).union_all()
    if hasattr(merged, "geoms"):
        merged_geoms = list(merged.geoms)
    else:
        merged_geoms = [merged]

    merged = gpd.GeoDataFrame(geometry=merged_geoms, crs=crs)
    centroids = merged.centroid
    merged["x"] = centroids.x
    merged["y"] = centroids.y

    logging.info("Projecting centroids...")
    # RESTORE ORIGINAL LOGIC for centroid projection
    centroids_projected = centroids.geometry.apply(
        lambda point: create_nearest_point(point, roads["geometry"])
    )

    logging.info("Splitting roads...")
    split_points = MultiPoint(list(centroids_projected.geometry))
    roads_union = roads.geometry.union_all()

    # Snap points TO roads first (this is the key!)
    snapped_roads = snap(roads_union, split_points, tolerance=1e-6)
    roads_splitted = split(snapped_roads, split_points)

    new_roads_nodes = (
        gpd.GeoSeries(roads_splitted, crs=crs)
        .explode(ignore_index=False, index_parts=False)
        .boundary.explode(ignore_index=False, index_parts=False)
        .drop_duplicates()
    )

    logging.info("Processing internal nodes...")
    # Create a subset of all internal nodes (from roads, splitted roads and centroids)
    internal_nodes = pd.concat(
        [road_intersections, centroids_projected, new_roads_nodes]
    ).drop_duplicates()
    internal_nodes_merged = (
        internal_nodes.buffer(10e-6).fillna(internal_nodes.geometry).union_all()
    )
    internal_nodes_merged = (
        internal_nodes_merged.geoms
        if hasattr(internal_nodes_merged, "geoms")
        else internal_nodes_merged
    )
    internal_nodes_merged = gpd.GeoSeries(internal_nodes_merged, crs=crs)
    internal_nodes_merged = gpd.GeoDataFrame(geometry=internal_nodes_merged)
    internal_nodes_merged_centroids = internal_nodes_merged.centroid

    logging.info("Creating node and edge DataFrames...")
    # Create DataFrames more efficiently
    n_internal = len(internal_nodes_merged_centroids)
    n_houses = len(sinks_center)
    n_prod = len(sources)

    # Pre-allocate arrays for better performance
    node_types = ["internal"] * n_internal + ["source"] * n_prod + ["sink"] * n_houses
    geometries = (
        list(internal_nodes_merged_centroids.geometry)
        + list(sources.geometry)
        + list(sinks_center.geometry)
    )

    gdf_nodes = gpd.GeoDataFrame(
        {
            "Type": node_types,
            "Node_ID": [f"Node_{i:05d}" for i in range(len(geometries))],
        },
        geometry=geometries,
        crs=crs,
    )

    logging.info("Creating final connection lines...")
    connection_lines_final = create_edge_nearest_point_optimized(
        union_src_snk.geometry, internal_nodes
    )
    new_roads_edges = gpd.GeoSeries(roads_splitted, crs=crs).explode(
        ignore_index=False, index_parts=False
    )

    logging.info("Processing edges...")
    edges_total = pd.concat(
        [new_roads_edges, connection_lines_final], ignore_index=True
    )

    # Remove short edges
    buffer_values = np.nonzero(edges_total.length <= 10e-6)[0]
    edges_total.drop(inplace=True, index=buffer_values.tolist())
    edges_total.reset_index(inplace=True, drop=True)

    gdf_edges = gpd.GeoDataFrame(
        {
            "Length": edges_total.length,
            "Road_ID": [f"Road_{i:05d}" for i in range(len(edges_total))]
        },
        geometry=edges_total.geometry,
        crs=crs,
    )

    ts_columns = sorted(
        [col for col in sinks_center.columns if col.startswith("ts_")],
        key=lambda x: int(x.split("_")[1]),
    )
    flh_columns = sorted(
        [col for col in sinks_center.columns if col.startswith("flh_")],
        key=lambda x: int(x.split("_")[1]),
    )

    gdf_nodes["q_c"] = ""
    gdf_nodes["flh_sinks"] = ""

    gdf_nodes.loc[gdf_nodes["Type"] == "sink", "q_c"] = np.round(sinks_center[ts_columns].values.T, 2).T
    gdf_nodes.loc[gdf_nodes["Type"] == "sink", "flh_sinks"] = np.round(sinks_center[flh_columns].values.T, 2).T

    logging.info("Creating node-edge relationships...")
    # Pre-allocate u and v columns
    gdf_edges["u"] = ""
    gdf_edges["v"] = ""

    for j, _ in enumerate(gdf_edges["Road_ID"]):
        road_geom = gdf_edges.iloc[j].geometry
        boundary_points = list(road_geom.boundary.geoms)

        if len(boundary_points) >= 2:
            start_point = boundary_points[0]
            end_point = boundary_points[-1]

            u_candidates = gdf_nodes[gdf_nodes.geometry.distance(start_point) < 1e-5]
            v_candidates = gdf_nodes[gdf_nodes.geometry.distance(end_point) < 1e-5]

            if len(u_candidates) == 1 and len(v_candidates) == 1:
                u = u_candidates.iloc[0]["Node_ID"]
                v = v_candidates.iloc[0]["Node_ID"]

                gdf_edges.loc[j, ["u", "v"]] = (u, v)

    return gdf_nodes, gdf_edges


def create_matrices_from_gdf(
        gdf_nodes: gpd.GeoDataFrame,
        gdf_edges: gpd.GeoDataFrame,
):
    """Create the necessary matrices for the optimization in topotherm from GeoDataFrames.

    Parameters
    ----------
    gdf_nodes : pd.DataFrame
        Nodes df containing information on sinks, sources and junctions
    gdf_edges : pd.DataFrame
        Edges df containing information on how all nodes are
        connected to each other.
    Returns
    -------
    mat : dict
        Incidence matrices and related data.
    """

    n_nodes = len(gdf_nodes)
    n_roads = len(gdf_edges)

    mat = {}
    mat["a_i"] = np.zeros([n_nodes, n_roads], dtype="int8")
    mat["l_i"] = gdf_edges["Length"].values

    # Create lookup dictionary for faster access
    node_id_to_idx = {node_id: idx for idx, node_id in enumerate(gdf_nodes["Node_ID"])}

    for j, _ in enumerate(gdf_edges["Road_ID"]):

        u = gdf_edges.iloc[j]["u"]
        v = gdf_edges.iloc[j]["v"]

        lu = node_id_to_idx[u]
        r = node_id_to_idx[v]

        if lu != r:
            mat["a_i"][lu, j] = 1
            mat["a_i"][r, j] = -1

    logging.info("Creating producer and consumer matrices...")
    # Create matrices for the heat sources (A_p)
    gdf_nodes_src = gdf_nodes[gdf_nodes["Type"] == "source"].index
    mat["a_p"] = np.zeros([len(gdf_nodes), len(gdf_nodes_src)], dtype="int")
    mat["a_p"][gdf_nodes_src, range(len(gdf_nodes_src))] = -1

    # Create matrices for the heat sinks (A_c)
    gdf_nodes_cons = gdf_nodes[gdf_nodes["Type"] == "sink"].index
    mat["a_c"] = np.zeros([len(gdf_nodes), len(gdf_nodes_cons)], dtype="int")
    mat["a_c"][gdf_nodes_cons, range(len(gdf_nodes_cons))] = 1

    logging.info("Creating final arrays...")
    mat["positions"] = np.column_stack(
        [gdf_nodes.geometry.x.values, gdf_nodes.geometry.y.values]
    )

    # Can we somehow automate this? Easiest way -> Probably just put them into one column
    # Workaround: extract and sort all integers that start with flh_ or ts_
    gdf_nodes_sinks = gdf_nodes[gdf_nodes["Type"] == "sink"].copy()
    mat["q_c"] = gdf_nodes_sinks["q_c"].values.T
    mat["flh_sinks"] = gdf_nodes_sinks["flh_sinks"].values.T
    mat["flh_sources"] = (
        np.round(
            np.array(
                [(mat["q_c"] * mat["flh_sinks"]).sum(axis=0) / mat["q_c"].sum(axis=0)]
            ),
            2,
        ).transpose()
        * np.ones(len(gdf_nodes_src))
    ).transpose()

    duplicates = find_duplicate_cols(mat["a_i"])

    if duplicates:
        dup_arr = np.array(duplicates)
        # Compare the weights l_i at each duplicate pair
        left_vals = mat["l_i"][dup_arr[:, 0]]
        right_vals = mat["l_i"][dup_arr[:, 1]]

        # Keep the larger one -> delete the smaller
        delete_idx = np.where(left_vals > right_vals, dup_arr[:, 0], dup_arr[:, 1])
        delete_idx = np.unique(delete_idx)

        # Delete in one shot
        mat["l_i"] = np.delete(mat["l_i"], delete_idx, axis=0)
        mat["a_i"] = np.delete(mat["a_i"], delete_idx, axis=1)
        gdf_edges = gdf_edges.drop(index=delete_idx).reset_index(
            drop=True
        )

    return mat, gdf_nodes, gdf_edges
