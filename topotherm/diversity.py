import copy

import numpy as np
import networkx as nx
import pandas as pd

np.set_printoptions(legacy='1.25')


def get_diversity_factor(graph: nx.DiGraph) -> np.array:
    """
    Calculate the diversity factor of a given graph of the optimal district heating network.

    Args:
        graph: A NetworkX graph object representing the network.

    Returns:
        np.array: A Numpy array containing the diversity factors of the graph.
    """
    # Compute the Laplacian matrix for the graph
    laplace = compute_laplacian(graph)

    # Resolve loops in the graph and update the Laplacian matrix
    G_copy, laplace = resolve_loops(graph, laplace)

    # Resolve multiple connections between nodes and update the Laplacian matrix
    laplace = resolve_multi_connections(G_copy, laplace)

    # Identify and mark end consumers in the graph
    graph = get_end_consumers(graph, laplace)

    # Calculate the final diversity factor based on the processed graph
    graph_diversity = calculate_diversity_factor(graph)

    return graph_diversity


def compute_laplacian(graph: nx.DiGraph) -> pd.DataFrame:
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
    laplace = pd.DataFrame(out_degree - adjacency_matrix, 
                           columns=graph.nodes, 
                           index=graph.nodes)
    return laplace


def resolve_loops(graph: nx.DiGraph,
                  laplace: pd.DataFrame) -> tuple:
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
        return [
            node for node, data in graph.nodes(data=True)
            if data.get(attribute) == value
        ]

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
                    laplace.loc[loop[-1], loop[0]] = 0
                    laplace.loc[loop[-1], loop[-1]] -= 1
                    G_copy.remove_edge(loop[-1], loop[0])
                else:
                    laplace.loc[loop[0], loop[-1]] = 0
                    laplace.loc[loop[0], loop[0]] -= 1
                    G_copy.remove_edge(loop[0], loop[-1])
                laplace_updated = True
            if laplace_updated:
                break

    return G_copy, laplace


def resolve_multi_connections(G_copy, laplace):
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
                laplace.loc[node, node] -= (value - 1)

    return laplace



def get_end_consumers(graph, laplace):
    """
    Identify and mark end consumers in the graph based on the Laplacian matrix.

    Args:
        graph: A NetworkX graph object.
        laplace: The Laplacian matrix of the graph.

    Returns:
        The graph with 'end consumers' attributes added to nodes.
    """
    while not laplace.empty:
        for i in laplace.index.tolist():
            laplace_single_row = laplace.loc[i]
            connected_nodes_no = laplace_single_row[i]

            laplace_connected_nodes = laplace_single_row.drop(i)
            is_isolated = not (laplace_connected_nodes != 0).any()

            if is_isolated:
                graph.nodes[i]['end consumers'] = connected_nodes_no
                if connected_nodes_no != 0:
                    laplace_single_column = laplace[i]
                    for j in laplace_single_column.index.tolist():
                        if laplace_single_column[j] != 0:
                            laplace.loc[j, j] += connected_nodes_no - 1
                laplace.drop(index=i, columns=i, inplace=True)

    return graph


def calculate_diversity_factor(graph):
    """
    Calculate the diversity factor for edges in the graph.

    Args:
        graph: A NetworkX graph object.

    Returns:
        np.array: A Numpy array that has edge names and values for the diversity factor
    """
    for u, v in graph.edges():
        end_consumers_u = graph.nodes[u].get('end consumers', 0)
        end_consumers_v = graph.nodes[v].get('end consumers', 0)

        max_end_consumers = max(end_consumers_u, end_consumers_v)

        # Calculate the diversity factor based on the maximum number of end consumers
        diversity_factor = (
            0.449677646267461 +
            (0.551234688 / (1 + pow((max_end_consumers / 53.84382392), 1.762743268)))
        )

        if diversity_factor>1:
            diversity_factor=1
        graph.edges[u, v]['diversity factor'] = diversity_factor

    name=np.array(list(nx.get_edge_attributes(graph,'diversity factor')))
    # create a numpy array with edge names in the first column and values in the second 
    col1 = np.core.defchararray.add(name[:, 0].astype(str), ', ')
    col1 = np.core.defchararray.add(col1, name[:, 1].astype(str))
    col2 = np.array(list(nx.get_edge_attributes(graph,'diversity factor').values()))
    diversity_dataframe = pd.DataFrame({'Name': col1, 'Diversity Factor': col2})
    return diversity_dataframe

def compare(power, factors):
    
    """
    Match the diversity factors with the correct edges and multiply with powr.

    Args:
        power: An array of supplied power at each edge
        factor: An array of diversity factors for each edge.

    Returns:
        matches: An array of matched edges with updated power values
    """
    
    merged_df = pd.merge(power, factors, on='Name', how='inner')
    merged_df['revised power'] = merged_df['power'] * merged_df['Diversity Factor']
    return merged_df[['Name', 'revised power']]
