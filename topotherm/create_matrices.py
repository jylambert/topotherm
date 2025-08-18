import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.ops import nearest_points, split
from shapely.geometry import LineString
from geopandas.tools import sjoin
from scipy.spatial import cKDTree

from topotherm.utils import create_dir


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
    line = LineString([point, nearest_point])

    return line


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


def create_edge_nearest_point(points1, points2):
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

    # Find the nearest point in points2 for each point in points1
    distances, indices = tree.query(coords1)

    # Create LineStrings connecting each point in points1 to the nearest point in points2
    lines = [LineString([points1.iloc[i], points2.iloc[indices[i]]]) for i in range(len(points1))]

    # Convert to GeoSeries for further processing or plotting
    lines_gs = gpd.GeoSeries(lines, crs="EPSG:25832")

    return lines_gs


def create_matrices(inputpath, outputpath, buffer=2.5):
    """
    Creates from a set of shapefiles ()

    Parameters
    ----------
    inputpath: Path to the initial shape files for the street layout (street.shp), heat sinks (heat_sinks.shp)
               and heat sources (heat_sources.shp)
    outputpath: Path to the result folder
    buffer: Radius of the buffer layer in meter, which is used to aggregate connection lines

    Returns
    -------
    Parquet files, needed for the optimization of the district heating network
    Shapefiles of the nodes and edges of the district
    """
    # Create the results folder
    create_dir(outputpath)

    # Load the shape-file containing the heat sinks
    houses = gpd.read_file(inputpath + 'heat_sinks.shp')
    houses = houses.to_crs(crs="EPSG:25832")            # Ensure the same EPSG for each file

    # Load the shape-file containing the streets
    roads = gpd.read_file(inputpath + 'streets.shp')
    roads = roads.to_crs(crs="EPSG:25832")              # Ensure the same EPSG for each file

    # Load the shape-file containing the heat source
    producer = gpd.read_file(inputpath + 'heat_source.shp')
    producer = producer.to_crs(crs="EPSG:25832")        # Ensure the same EPSG for each file

    # Create the intersections at each end of the edges and remove duplicate points
    road_intersections = roads.boundary.explode(ignore_index=False, index_parts=False).drop_duplicates()

    # If needed, create a centroid for each heat sink
    houses_center = houses.to_crs(crs="EPSG:25832")
    houses_center.geometry = houses.centroid

    # Perform a union between the heat sinks and sources centroid to connect those to the streets
    union_prod_con = pd.concat([houses_center['geometry'], producer.geometry])

    # Create initial connection line to the road network
    connection_line = union_prod_con.apply(lambda point: create_connection_line(point, roads['geometry']))
    connection_nodes = connection_line.boundary.explode(ignore_index=False, index_parts=False)
    matches = connection_nodes.isin(houses_center.geometry)
    connection_nodes = connection_nodes[~matches].drop_duplicates()

    # Find near nodes and merge them (inspired by https://github.com/gboeing/osmnx/blob/main/osmnx/simplification.py)
    # Buffer nodes to passed-in distance and merge overlaps.
    # Turn merged nodes into gdf and get centroids of each cluster as x, y
    merged = connection_nodes.buffer(buffer).fillna(connection_nodes.geometry).unary_union
    merged = merged.geoms if hasattr(merged, "geoms") else merged
    merged = gpd.GeoSeries(merged, crs="EPSG:25832")
    merged = gpd.GeoDataFrame(geometry=merged)
    centroids = merged.centroid
    merged["x"] = centroids.x
    merged["y"] = centroids.y

    # Project the centroids back to the road network
    centroids_projected = centroids.geometry.apply(lambda point: create_nearest_point(point, roads['geometry']))
    roads_splitted = split(roads.geometry.unary_union, centroids_projected.buffer(1.0E-9).unary_union)
    new_roads_nodes = gpd.GeoSeries(roads_splitted, crs="EPSG:25832").explode(ignore_index=False, index_parts=False)\
        .boundary.explode(ignore_index=False, index_parts=False).drop_duplicates()

    # Create a subset of all internal nodes (from roads, splitted roads and centroids)
    internal_nodes = pd.concat([road_intersections, centroids_projected, new_roads_nodes]).drop_duplicates()
    internal_nodes_merged = internal_nodes.buffer(10e-4).fillna(internal_nodes.geometry).unary_union
    internal_nodes_merged = internal_nodes_merged.geoms \
        if hasattr(internal_nodes_merged, "geoms") else internal_nodes_merged
    internal_nodes_merged = gpd.GeoSeries(internal_nodes_merged, crs="EPSG:25832")
    internal_nodes_merged = gpd.GeoDataFrame(geometry=internal_nodes_merged)
    internal_nodes_merged_centroids = internal_nodes_merged.centroid

    # Prepare the GeoDataFrame for the nodes
    gdf_nodes = gpd.GeoDataFrame({'Type': {}}, geometry=internal_nodes_merged_centroids.geometry)
    gdf_nodes['Type'] = 'Int'

    # Prepare the GeoDataFrame for the edges
    gdf_houses_center = gpd.GeoDataFrame({'Type': {}}, geometry=houses_center.geometry)
    gdf_houses_center['Type'] = 'Cons'

    # Prepare the GeoDataFrame for the heat sources
    gdf_prod = gpd.GeoDataFrame({'Type': {}}, geometry=producer.geometry)
    gdf_prod['Type'] = 'Prod'

    # Perform a union of all nodes
    gdf_nodes = pd.concat([gdf_nodes, gdf_prod, gdf_houses_center])
    gdf_nodes['Node_ID'] = [f"Node_{i:05d}" for i in range(gdf_nodes.geometry.size)]

    # Perform a final connection line form the heat sources and sinks to the internal road nodes
    # Here the closest node is chosen for the connection line
    # Van either be an original road node or a splitted road node
    connection_lines_final = create_edge_nearest_point(union_prod_con.geometry, internal_nodes)
    new_roads_edges = gpd.GeoSeries(roads_splitted, crs="EPSG:25832").explode(ignore_index=False, index_parts=False)

    # Perform a union for all edges (connection lines and splitted roads)
    edges_total = pd.concat([new_roads_edges, connection_lines_final])
    buffer_values = np.nonzero(edges_total.length <= 10e-6)[0]
    edges_total.reset_index(inplace=True, drop=True)
    edges_total.drop(inplace=True, index=buffer_values.tolist())
    gdf_road = gpd.GeoDataFrame({'Length': edges_total.length}, geometry=edges_total.geometry)
    gdf_road['Road_ID'] = [f"Road_{i:05d}" for i in range(edges_total.size)]

    # Perform a buffer to ensure that the end point of two lines overlap
    gdf_nodes_buf = gdf_nodes.copy()
    gdf_nodes_buf.geometry = gdf_nodes.geometry.buffer(10e-6)
    gdf_nodes_buf['geo_unbuffered'] = gdf_nodes.geometry

    # Map buffered and unbuffered nodes
    intersections_fin = sjoin(gdf_nodes_buf, gdf_road, how='inner', predicate='intersects')
    edge_name = gdf_road['Road_ID'].reset_index(drop=True)
    gdf_road['u'] = {}
    gdf_road['v'] = {}

    # Create matrice A_i and vector l_i and fill them with 0
    a_i = np.zeros([len(gdf_nodes), len(gdf_road)], dtype="int")
    l_i = np.zeros(len(gdf_road))

    # Reset index to get nodes and edges starting from 0
    gdf_road = gdf_road.reset_index(drop=True)
    gdf_nodes = gdf_nodes.reset_index(drop=True)

    for j in range(edge_name.size):
        i = edge_name[j]
        inter = intersections_fin[intersections_fin['Road_ID'] == i]
        point1 = inter.iloc[0, 3]
        point2 = inter.iloc[1, 3]
        u = gdf_nodes.iloc[np.nonzero(gdf_nodes_buf['geo_unbuffered'] == point1)[0][0], 2]
        v = gdf_nodes.iloc[np.nonzero(gdf_nodes_buf['geo_unbuffered'] == point2)[0][0], 2]
        gdf_road.loc[gdf_road['Road_ID'] == i, 'u'] = u
        gdf_road.loc[gdf_road['Road_ID'] == i, 'v'] = v

        l_i[j] = gdf_road.loc[gdf_road['Road_ID'] == i, 'Length'].iloc[0]
        lu = gdf_nodes[gdf_nodes['Node_ID'] == u].index[0]
        r = gdf_nodes[gdf_nodes['Node_ID'] == v].index[0]
        if lu != r:
            a_i[lu, j] = 1
            a_i[r, j] = -1
        else:
            pass

    # Create A_p
    gdf_nodes_prod = gdf_nodes[gdf_nodes['Type'] == 'Prod'].index
    a_p = np.zeros([len(gdf_nodes), len(gdf_nodes_prod)], dtype="int")
    a_p[gdf_nodes_prod, range(len(gdf_nodes_prod))] = -1

    # Create A_c
    gdf_nodes_cons = gdf_nodes[gdf_nodes['Type'] == 'Cons'].index
    a_c = np.zeros([len(gdf_nodes), len(gdf_nodes_cons)], dtype="int")
    a_c[gdf_nodes_cons, range(len(gdf_nodes_cons))] = 1

    # Create Pos, q_c, flh and flh_prod
    pos = np.array([gdf_nodes.geometry.x, gdf_nodes.geometry.y]).transpose()
    q_c = np.array(houses_center["peak_power"].astype(float)).transpose()
    flh = np.array(houses_center["full_load"].astype(float)).transpose()
    flh_prod = np.round(np.array([(q_c*flh).sum() / q_c.sum()]), 2)

    # Save the road network with the connection lines and the nodes as a shape file
    gdf_road.to_file(outputpath + "road_fin.shp")
    gdf_nodes.to_file(outputpath + "nodes_fin.shp")

    # Create the matrices needed for the optimization
    pd.DataFrame(a_i).to_parquet(outputpath + 'A_i.parquet')
    pd.DataFrame(a_p).to_parquet(outputpath + 'A_p.parquet')
    pd.DataFrame(a_c).to_parquet(outputpath + 'A_c.parquet')
    pd.DataFrame(q_c).to_parquet(outputpath + 'Q_c.parquet')
    pd.DataFrame(l_i).to_parquet(outputpath + 'L_i.parquet')
    pd.DataFrame(flh).to_parquet(outputpath + 'flh_consumer.parquet')
    pd.DataFrame(flh_prod).to_parquet(outputpath + 'flh_source.parquet')
    pd.DataFrame(pos).to_parquet(outputpath + 'rel_positions.parquet')

    return

