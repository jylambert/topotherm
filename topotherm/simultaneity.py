"""Calculate the simultaneity factor of a given graph of the optimal district heating network.
"""
import copy

import numpy as np
import networkx as nx
import pandas as pd


def calculate(graph: nx.DiGraph) -> nx.DiGraph:
    """
    Calculate the simultaneity factor of a given graph of the optimal district
    heating network.

    Args:
        graph: A NetworkX graph object representing the network.

    Returns:
         nx.DiGraph: A graph containing the simultaneity factors of the graph
         in the attribute 'simultaneity'.
    """
    # Compute the Laplacian matrix for the graph
    laplace = compute_laplacian(graph)

    # Resolve loops in the graph and update the Laplacian matrix
    G_copy, laplace = resolve_loops(graph, laplace)

    # Resolve multiple connections between nodes and update the Laplacian matrix
    laplace = resolve_multi_connections(G_copy, laplace)

    # Identify and mark end consumers in the graph
    graph = get_n_end_consumers(graph, laplace)

    # Calculate the final simultaneity factor based on the processed graph
    graph_simultaneity = calculate_simultaneity_factor(graph)

    return graph_simultaneity


def compute_laplacian(graph: nx.DiGraph) -> np.array:
    """
    Compute the Laplacian matrix of a given graph.

    Args:
        graph: A NetworkX graph object of the solved network connections.

    Returns:
        pd.DataFrame: A pandas DataFrame representing the Laplacian matrix.
    """
    # Calculate the out-degree as a diagonal matrix
    out_degree = np.diag([graph.out_degree(node) for node in graph.nodes])

    # Convert the adjacency matrix to a numpy array
    adjacency_matrix = nx.to_numpy_array(graph, nodelist=graph.nodes)

    # Compute the Laplacian matrix (D - A)
    laplace = out_degree - adjacency_matrix

    return laplace


def resolve_loops(graph: nx.DiGraph,
                  laplace: pd.DataFrame) -> tuple[nx.DiGraph, np.array]:
    """
    Resolve loops in the graph and update the Laplacian matrix accordingly.

    Args:
        graph: A NetworkX graph object.
        laplace: The Laplacian matrix of the graph.

    Returns:
        A tuple containing the updated graph and Laplacian matrix.
    """
    def find_nodes_by_attribute(graph: nx.DiGraph, attribute: str, value) -> list:
        """
        Find nodes in a graph that have a specific attribute value.

        Args:
            graph: The NetworkX graph.
            attribute: The attribute to search for.
            value: The attribute value to match.

        Returns:
            list: A list of nodes with the specified attribute value.
        """
        return [node for node, data in graph.nodes(data=True)
                if data.get(attribute) == value]

    # Find nodes labeled as 'Heat Source'
    heat_source_nodes = find_nodes_by_attribute(graph, 'label', 'Heat Source')
    G_copy = copy.deepcopy(graph)

    # Identify and resolve simple cycles in the graph
    for loop in nx.recursive_simple_cycles(graph):
        laplace_updated = False
        for node in heat_source_nodes:
            if nx.has_path(graph, node, loop[0]) and nx.has_path(graph, node, loop[-1]):
                if nx.shortest_path_length(graph, 1, loop[-1]) > nx.shortest_path_length(graph, 1, loop[0]):
                    # Update Laplacian matrix and remove edge from the graph
                    laplace[loop[-1], loop[0]] = 0
                    laplace[loop[-1], loop[-1]] -= 1
                    G_copy.remove_edge(loop[-1], loop[0])
                else:
                    laplace[loop[0], loop[-1]] = 0
                    laplace[loop[0], loop[0]] -= 1
                    G_copy.remove_edge(loop[0], loop[-1])
                laplace_updated = True
            if laplace_updated:
                break

    return G_copy, laplace


def resolve_multi_connections(G_copy: nx.DiGraph, 
                              laplace: np.array) -> np.array:
    """
    Resolves multiple incoming connections for nodes in a directed graph by propagating 'connections' attributes 
    and updating the Laplacian matrix accordingly.

    Parameters:
    G_copy: A copy of the directed graph from the previous step.
    laplace: The Laplacian matrix to be updated.

    Returns:
    laplace: Updated Laplacian matrix with modified diagonal elements.
    """

    # Initialize 'connections {node}' attributes for nodes with in-degree > 1
    for node in G_copy.nodes:
        if G_copy.in_degree(node) > 1:
            G_copy.nodes[node][f'connections {node}'] = 1

    # Function to check for 'connections' attributes in any node
    def has_connection_attribute(graph):
        return any(
            "connections" in str(key).lower()
            for _, attrs in graph.nodes(data=True)
            for key in attrs.keys()
        )

    # Propagate 'connections' attribute
    while has_connection_attribute(G_copy):
        updates_occurred = False  # Flag
        print(G_copy.nodes.data())

        # Process nodes in normal order
        for node in G_copy.nodes:
            connection_attrs = {
                key: value for key, value in G_copy.nodes[node].items() if 'connections' in key
            }

            for key in connection_attrs.keys():
                # Count successors with the same connections attribute
                count_of_successors = sum(
                    G_copy.nodes[succ].get(key, 0) == G_copy.nodes[node].get(key, 0)
                    for succ in G_copy.successors(node)
                )

                # Update the node's connections attribute if the count is greater
                current_value = G_copy.nodes[node].get(key, 0)
                if count_of_successors > current_value:
                    G_copy.nodes[node][key] = count_of_successors
                    updates_occurred = True

            # Update predecessors if they lack the attribute
            for pred in G_copy.predecessors(node):
                for key, value in connection_attrs.items():
                    if key not in G_copy.nodes[pred]:
                        G_copy.nodes[pred][key] = 1
                        updates_occurred = True

        # Break the loop if no updates occurred
        if not updates_occurred:
            break

    for node in G_copy.nodes:
        connection_attrs = {
                key: value for key, value in G_copy.nodes[node].items() if 'connections' in key
            }
        for key, value in connection_attrs.items():
            if value > 1:
                laplace[node, node] -= (value - 1)

    return laplace


def get_n_end_consumers(graph: nx.DiGraph, laplace:np.array) -> nx.DiGraph:
    """
    Identify and mark end consumers in the graph based on the Laplacian matrix.

    Args:
        graph: A NetworkX graph object.
        laplace: The Laplacian matrix of the graph.

    Returns:
        The graph with 'n_consumers' attributes added to nodes.
    """
    laplace = pd.DataFrame(laplace, index=graph.nodes, columns=graph.nodes)

    while not laplace.empty:
        for i in laplace.index.tolist():
            laplace_single_row = laplace.loc[i]
            connected_nodes_no = laplace_single_row[i]

            laplace_connected_nodes = laplace_single_row.drop(i)
            is_isolated = not (laplace_connected_nodes != 0).any()

            if is_isolated:
                graph.nodes[i]['n_consumers'] = connected_nodes_no

                if connected_nodes_no != 0:
                    laplace_single_column = laplace[i]
                    for j in laplace_single_column.index.tolist():
                        if laplace_single_column[j] != 0:
                            laplace.loc[j, j] += connected_nodes_no - 1
                laplace.drop(index=i, columns=i, inplace=True)
    return graph


def calculate_simultaneity_factor(graph: nx.DiGraph) -> nx.DiGraph:
    """
    Calculate the simultaneity factor for edges in the graph.

    Args:
        graph: A NetworkX graph object with number of end consumers marked as 
        n_consumers attributes.

    Returns:
        nx.DiGraph: A graph that has edge names and values for the
        simultaneity factor
    """
    for u, v in graph.edges():
        end_consumers_u = graph.nodes[u].get('n_consumers', 0)
        end_consumers_v = graph.nodes[v].get('n_consumers', 0)

        max_end_consumers = max(end_consumers_u, end_consumers_v)

        # Calculate the simultaneity factor based on the maximum number of end consumers
        simultaneity_factor = (
            0.449677646267461
            + (0.551234688
               / (1 + pow((max_end_consumers / 53.84382392), 1.762743268))
               )
        )

        graph.edges[u, v]['simultaneity'] = min(simultaneity_factor, 1)

    return graph


def update_data(
        df_nodes: pd.DataFrame,
        df_edges: pd.DataFrame,
        network: nx.DiGraph
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Update the edge and node data with simultaneity attributes.

    Args:
        node_data: A DataFrame containing edge data.
        edge_data: A DataFrame containing node data.
        network_simultaneity: A NetworkX graph object with simultaneity factors,
        needs to contain power, simultaneity attributes
    
    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: Updated edge and node data DataFrames.
    """
    # p_simultaneity = power * simultaneity to the network_simultaneity object
    for u, v, d in network.edges(data=True):
        network[u][v]['power_simultaneity'] = d['simultaneity'] * d['p']

    # get simultaneity attributes from networkx object, add to edge_data
    df_edges_sim = pd.DataFrame(
        [(u, v, d['simultaneity'], d['power_simultaneity'])
         for u, v, d in network.edges(data=True)],
        columns=['start_node', 'end_node', 'simultaneity', 'power_simultaneity']
    ).set_index(['start_node', 'end_node'])

    df_edges = df_edges.join(df_edges_sim, on=['start_node', 'end_node'])

    n_consumers_dict = {
        n: d.get('n_consumers', None)
        for n, d in network.nodes(data=True)}

    df_nodes['n_consumers'] = pd.Series(n_consumers_dict)

    return df_nodes, df_edges
