"""Sone standard plotting functions for the district heating network."""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


def district(matrices: dict,
             diameter: list[float] = [0],
             isnot_init: bool = False) -> plt.Figure:
    """
    Plot the district heating network with the given matrices.

    Parameters
    ----------
    matrices : dict
        Dictionary with the matrices of the district heating network with keys:

        - **a_i**       : Incidence matrix with rows: nodes, columns: edges.
        - **a_p**       : Adjacency matrix for the producers with rows: nodes, columns: supply ids.
        - **a_c**       : Adjacency matrix for the consumers with rows: nodes, columns: consumer ids.
        - **l_i**       : Length of edges.
        - **position**  : x, y coordinates of the nodes in the network.

    diameter : list of float
        Inner diameter of the pipes (default = [0]).

    isnot_init : bool
        True to plot after optimization, False for before (default = False).

    Returns
    -------
    plt.Figure
        Figure of the district.
    """


    G = nx.DiGraph()
    s = np.array([0, 0, 0])

    # Add the nodes to the graph
    sums = matrices['a_c'].sum(axis=1)
    prod = matrices['a_p'].T.sum(axis=0)
    ges = sums + prod

    ges = ges.tolist()
    ges = np.array(ges).squeeze()

    for q in range(matrices['a_c'].shape[0]):
        if ges[q] == 1:
            G.add_node(q, color='LightBlue', label='Consumer')
        elif ges[q] == 0:
            G.add_node(q, color='LightGrey', label='Internal Node')
        if ges[q] <= -1:
            G.add_node(q, color='Orange', label='Heat Source')

    # Add the edges to the graph
    for k in range(matrices['a_i'].shape[1]):
        for q in range(matrices['a_i'].shape[0]):
            if matrices['a_i'][q, k] == 1:
                s[0] = q
            elif matrices['a_i'][q, k] == -1:
                s[1] = q
        if isnot_init:
            G.add_edge(s[0], s[1], weight=10 * (np.round(diameter[k], 3)/0.3938), len=matrices['l_i'][k])
        else:
            G.add_edge(s[0], s[1], weight=10, len=matrices['l_i'][k])

        # Get the node colors from the color attribute
    node_colors = [G.nodes[node]['color'] for node in G.nodes()]
    node_label = [G.nodes[node]['label'] for node in G.nodes()]

    cm = 1 / 2.54
    # dhn topology plot
    fig0, ax = plt.subplots(figsize=(2*90 * cm, 2*60 * cm),
                            layout='constrained')
    if isnot_init:
        nx.draw_networkx_nodes(G, matrices['position'],
                               node_color=node_colors,
                               node_size=55,
                               ax=ax)
        for edge in G.edges(data='weight'):
            nx.draw_networkx_edges(G, matrices['position'],
                                   edgelist=[edge],
                                   width=edge[2],
                                   ax=ax,
                                   arrows=False)
    else:
        nx.draw(G, matrices['position'], node_color=node_colors, node_size=55, with_labels=False, ax=ax, width=2.5, label=node_label)

    legend_elements = [
        Line2D([0], [0],
               marker='o',
               color='w',
               label='Consumer',
               markerfacecolor='LightBlue',
               markersize=15),
        Line2D([0], [0],
               marker='o',
               color='w',
               label='Internal Node',
               markerfacecolor='LightGrey',
               markersize=15),
        Line2D([0], [0],
               marker='o',
               color='w',
               label='Heat Source',
               markerfacecolor='Orange',
               markersize=15)
    ]
    plt.legend(handles=legend_elements, loc='best', frameon=False, fontsize=22)
    plt.box(False)
    return fig0


def plot_networkx(n):
    """
    Plot a NetworkX graph.

    Parameters
    ----------
    n : networkx.Graph
        Input NetworkX graph.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure of the NetworkX graph.
    """

    fig, ax = plt.subplots(figsize=(15, 5))

    nx.draw(n, with_labels=True, font_weight='bold', ax=ax)
    plt.show()
    return fig

