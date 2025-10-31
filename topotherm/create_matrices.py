import logging

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.ops import nearest_points, split, snap
from shapely.geometry import LineString, MultiPoint
from scipy.spatial import cKDTree

from topotherm.utils import create_dir, find_duplicate_cols


def create_connection_line(point, edges):
    """
    Creates a connection line from a single point to multiple edges

    Parameters
    ----------
    point: Single point which should be connected
    edges: Geoseries of edges to which the point should be connected

    Returns
    -------
    line : A single linestring connecting the point to the nearest edge
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
    point: Single point which should be projected onto an edge
    edges: Geoseries of edges onto which the point should be projected

    Returns
    -------
    nearest point: A single nearest point object.
    """
    # Find the nearest point on the edges
    nearest_geom = edges.distance(point).idxmin()
    nearest_point = nearest_points(point, edges.iloc[nearest_geom])[1]
    return nearest_point


def create_edge_nearest_point_optimized(points1, points2):
    """
    Creates an edge for each points1 to the nearest points2

    Parameters
    ----------
    points1: Series of points which should be connected to the nearest points2
    points2: Series of points to which should be connected the points1

    Returns
    -------
    lines_gs: Geoseries of lines connecting points1 to points2
    """

    # Extract coordinates from GeoSeries
    coords1 = np.array(list(zip(points1.x, points1.y)))
    coords2 = np.array(list(zip(points2.x, points2.y)))

    # Use cKDTree for efficient nearest-neighbor search
    tree = cKDTree(coords2)
    distances, indices = tree.query(coords1)

    # Create LineStrings connecting each point in points1 to the nearest point in points2
    lines = [LineString([points1.iloc[i], points2.iloc[indices[i]]]) for i in range(len(points1))]

    return gpd.GeoSeries(lines, crs=points1.crs if hasattr(points1, 'crs') else "EPSG:25832")


def create_matrices(inputpath: dict, outputpath: str, buffer=2.5):
    """
    Creates from a set of shapefiles ()

    Parameters
    ----------
    inputpath: Path to the initial shape files for the street layout (roads.shp), heat sinks (heat_sinks.shp)
               and heat sources (heat_sources.shp) as a dictionary. Matrices can also be Geopackages
    outputpath: Path to the result folder
    buffer: Radius of the buffer layer in meter, which is used to aggregate connection lines

    Returns
    -------
    Parquet files, needed for the optimization of the district heating network
    Shapefiles of the nodes and edges of the district
    """
    # Create the results folder
    create_dir(outputpath)

    # Load shapefiles
    logging.info("Loading shapefiles...")
    inputpath_b = inputpath['sinks']
    inputpath_r = inputpath['roads']
    inputpath_s = inputpath['sources']

    sinks = gpd.read_file(inputpath_b).to_crs("EPSG:25832")
    roads = gpd.read_file(inputpath_r).to_crs("EPSG:25832")
    sources = gpd.read_file(inputpath_s).to_crs("EPSG:25832")

    print("Processing road intersections...")
    # Optimize road intersections processing
    road_intersections = roads.boundary.explode(ignore_index=False, index_parts=False).drop_duplicates()

    logging.info("Processing house centers...")
    # Create centroids more efficiently
    sinks_center = sinks.copy()
    sinks_center.geometry = sinks.centroid

    logging.info("Creating union...")
    # Optimize union creation
    union_src_snk = pd.concat([sinks_center[['geometry']], sources[['geometry']]], ignore_index=True)

    logging.info("Creating connection lines...")
    # Create connection lines between sources and sinks to the roadnetwork
    connection_line = union_src_snk.geometry.apply(lambda point: create_connection_line(point, roads['geometry']))
    connection_nodes = connection_line.boundary.explode(ignore_index=False, index_parts=False)

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

    merged = gpd.GeoDataFrame(geometry=merged_geoms, crs="EPSG:25832")
    centroids = merged.centroid
    merged["x"] = centroids.x
    merged["y"] = centroids.y

    logging.info("Projecting centroids...")
    # RESTORE ORIGINAL LOGIC for centroid projection
    centroids_projected = centroids.geometry.apply(lambda point: create_nearest_point(point, roads['geometry']))

    logging.info("Splitting roads...")
    split_points = MultiPoint(list(centroids_projected.geometry))
    roads_union = roads.geometry.union_all()

    # Snap points TO roads first (this is the key!)
    snapped_roads = snap(roads_union, split_points, tolerance=1e-6)
    roads_splitted = split(snapped_roads, split_points)

    new_roads_nodes = (gpd.GeoSeries(roads_splitted, crs="EPSG:25832")
                       .explode(ignore_index=False, index_parts=False)
                       .boundary
                       .explode(ignore_index=False, index_parts=False)
                       .drop_duplicates())

    logging.info("Processing internal nodes...")
    # Create a subset of all internal nodes (from roads, splitted roads and centroids)
    internal_nodes = pd.concat([road_intersections, centroids_projected, new_roads_nodes]).drop_duplicates()
    internal_nodes_merged = internal_nodes.buffer(10e-6).fillna(internal_nodes.geometry).union_all()
    internal_nodes_merged = internal_nodes_merged.geoms if hasattr(internal_nodes_merged, "geoms") else internal_nodes_merged
    internal_nodes_merged = gpd.GeoSeries(internal_nodes_merged, crs="EPSG:25832")
    internal_nodes_merged = gpd.GeoDataFrame(geometry=internal_nodes_merged)
    internal_nodes_merged_centroids = internal_nodes_merged.centroid

    logging.info("Creating node and edge DataFrames...")
    # Create DataFrames more efficiently
    n_internal = len(internal_nodes_merged_centroids)
    n_houses = len(sinks_center)
    n_prod = len(sources)

    # Pre-allocate arrays for better performance
    node_types = ['internal'] * n_internal + ['source'] * n_prod + ['sink'] * n_houses
    geometries = list(internal_nodes_merged_centroids.geometry) + list(sources.geometry) + list(sinks_center.geometry)

    gdf_nodes = gpd.GeoDataFrame({
        'Type': node_types,
        'Node_ID': [f"Node_{i:05d}" for i in range(len(geometries))]
    }, geometry=geometries, crs="EPSG:25832")

    logging.info("Creating final connection lines...")
    connection_lines_final = create_edge_nearest_point_optimized(union_src_snk.geometry, internal_nodes)
    new_roads_edges = gpd.GeoSeries(roads_splitted, crs="EPSG:25832").explode(ignore_index=False, index_parts=False)

    logging.info("Processing edges...")
    edges_total = pd.concat([new_roads_edges, connection_lines_final], ignore_index=True)

    # Remove short edges
    buffer_values = np.nonzero(edges_total.length <= 10e-6)[0]
    edges_total.drop(inplace=True, index=buffer_values.tolist())
    edges_total.reset_index(inplace=True, drop=True)

    gdf_road = gpd.GeoDataFrame({
        'Length': edges_total.length,
        'Road_ID': [f"Road_{i:05d}" for i in range(len(edges_total))]
    }, geometry=edges_total.geometry, crs="EPSG:25832")

    logging.info("Creating node-edge relationships...")
    n_nodes = len(gdf_nodes)
    n_roads = len(gdf_road)
    a_i = np.zeros([n_nodes, n_roads], dtype="int8")
    l_i = gdf_road['Length'].values

    # Create lookup dictionary for faster access
    node_id_to_idx = {node_id: idx for idx, node_id in enumerate(gdf_nodes['Node_ID'])}

    # Pre-allocate u and v columns
    gdf_road['u'] = ''
    gdf_road['v'] = ''

    for j, road_id in enumerate(gdf_road['Road_ID']):
        road_geom = gdf_road.iloc[j].geometry
        boundary_points = list(road_geom.boundary.geoms)

        if len(boundary_points) >= 2:
            start_point = boundary_points[0]
            end_point = boundary_points[-1]

            u_candidates = gdf_nodes[gdf_nodes.geometry.distance(start_point) < 1e-5]
            v_candidates = gdf_nodes[gdf_nodes.geometry.distance(end_point) < 1e-5]

            if len(u_candidates) == 1 and len(v_candidates) == 1:
                u = u_candidates.iloc[0]['Node_ID']
                v = v_candidates.iloc[0]['Node_ID']

                gdf_road.loc[j, 'u'] = u
                gdf_road.loc[j, 'v'] = v

                lu = node_id_to_idx[u]
                r = node_id_to_idx[v]

                if lu != r:
                    a_i[lu, j] = 1
                    a_i[r, j] = -1

    logging.info("Creating producer and consumer matrices...")
    # Create A_p
    gdf_nodes_prod = gdf_nodes[gdf_nodes['Type'] == 'Prod'].index
    a_p = np.zeros([len(gdf_nodes), len(gdf_nodes_prod)], dtype="int")
    a_p[gdf_nodes_prod, range(len(gdf_nodes_prod))] = -1

    # Create A_c
    gdf_nodes_cons = gdf_nodes[gdf_nodes['Type'] == 'Cons'].index
    a_c = np.zeros([len(gdf_nodes), len(gdf_nodes_cons)], dtype="int")
    a_c[gdf_nodes_cons, range(len(gdf_nodes_cons))] = 1

    logging.info("Creating final arrays...")
    # Create final arrays more efficiently
    pos = np.column_stack([gdf_nodes.geometry.x.values, gdf_nodes.geometry.y.values])

    # Extract and sort time-step columns dynamically
    ts_columns = sorted([col for col in sinks_center.columns if col.startswith("ts_")], key=lambda x: int(x.split('_')[1]))
    flh_columns = sorted([col for col in sinks_center.columns if col.startswith("flh_")], key=lambda x: int(x.split('_')[1]))

    # Can we somehow automate this? Easiest way -> Probably just put them into one column
    # Workaround: extract and sort all integers that start with flh_ or ts_
    q_c = np.round(sinks_center[ts_columns].values.T, 2).T
    flh = np.round(sinks_center[flh_columns].values.T, 2).T

    flh_prod = (np.round(np.array([(q_c*flh).sum(axis=0) / q_c.sum(axis=0)]), 2).transpose() * np.ones(2)).transpose()

    mat = {}
    mat['a_i'] = a_i
    mat['a_p'] = a_p
    mat['a_c'] = a_c
    mat['l_i'] = l_i
    mat['positions'] = pos
    mat['q_c'] = q_c
    mat['flh_sinks'] = flh
    mat['flh_sources'] = flh_prod

    duplicates = find_duplicate_cols(mat['a_i'])
    # Remove duplicate nodes
    if duplicates:
        # Convert to numpy array for easier slicing
        dup_arr = np.array(duplicates)
        # Compare the weights l_i at each duplicate pair
        left_vals = mat['l_i'][dup_arr[:, 0]]
        right_vals = mat['l_i'][dup_arr[:, 1]]

        # Keep the larger one -> delete the smaller
        delete_idx = np.where(left_vals > right_vals, dup_arr[:, 0], dup_arr[:, 1])

        # Unique + sorted indices for deletion
        delete_idx = np.unique(delete_idx)

        # Delete in one shot
        mat['l_i'] = np.delete(mat['l_i'], delete_idx, axis=0)
        mat['a_i'] = np.delete(mat['a_i'], delete_idx, axis=1)
        gdf_road = gdf_road.drop(index=delete_idx).reset_index(drop=True)

    print("Saving results...")

    # Save results
    gdf_road.to_file(outputpath + "edges.shp")
    gdf_nodes.to_file(outputpath + "nodes.shp")

    # Save matrices with compression for smaller files
    pd.DataFrame(mat['a_i']).to_parquet(outputpath + 'a_i.parquet')
    pd.DataFrame(mat['a_p']).to_parquet(outputpath + 'a_p.parquet')
    pd.DataFrame(mat['a_c']).to_parquet(outputpath + 'a_c.parquet')
    pd.DataFrame(q_c).to_parquet(outputpath + 'q_c.parquet')
    pd.DataFrame(mat['l_i']).to_parquet(outputpath + 'l_i.parquet')
    pd.DataFrame(flh).to_parquet(outputpath + 'flh_sinks.parquet')
    pd.DataFrame(flh_prod).to_parquet(outputpath + 'flh_sources.parquet')
    pd.DataFrame(mat['positions']).to_parquet(outputpath + 'positions.parquet')

    print("Processing complete!")

    return
