# -*- coding: utf-8 -*-
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


def plot_district(a_i, a_p, a_c, title, length, diameter, positions, isnot_init):
    """Input: a_i (internal matrix), a_p (producer matrix), a_c (consumer matrix)
    title (name of the figure), length (of the pipes), diameter (inner diameter of the pipes)
    positions (positions of the nodes), isnot_init (either 0 or 1, for plot before or after the optimization)
    Returns: Figure of the district"""
    G = nx.DiGraph()
    s = np.array([0, 0, 0])

    # Add the nodes to the graph
    sums = a_c.sum(axis=1)
    prod = a_p.T.sum(axis=0)
    ges = sums + prod

    ges = ges.tolist()
    ges = np.array(ges).squeeze()

    for q in range(a_c.shape[0]):
        if ges[q] == 1:
            G.add_node(q, color='LightBlue', label='Consumer')
        elif ges[q] == 0:
            G.add_node(q, color='LightGrey', label='Internal Node')
        if ges[q] == -1:
            G.add_node(q, color='Orange', label='Heat Source')

    edge_labels = dict()
    # Add the edges to the graph
    for k in range(a_i.shape[1]):
        for q in range(a_i.shape[0]):
            if a_i[q, k] == 1:
                s[0] = q
            elif a_i[q, k] == -1:
                s[1] = q
        if isnot_init:
            G.add_edge(s[0], s[1], weight=10 * (np.round(diameter[k], 3)/0.3938), len=length[k])
        else:
            G.add_edge(s[0], s[1], weight=10, len=length[k])

        # Get the node colors from the color attribute
    node_colors = [G.nodes[node]['color'] for node in G.nodes()]
    node_label = [G.nodes[node]['label'] for node in G.nodes()]
    edge_length = [G.edges[edge]['len'] for edge in G.edges()]

    cm = 1 / 2.54
    # dhn topology plot
    #fig0, ax = plt.subplots(figsize=(3*30 * cm, 3*20 * cm), layout="constrained")
    fig0, ax = plt.subplots(figsize=(2*90 * cm, 2*60 * cm), layout="constrained")
    if isnot_init:
        nx.draw_networkx_nodes(G, positions, node_color=node_colors, node_size=55, ax=ax)
        for edge in G.edges(data='weight'):
            nx.draw_networkx_edges(G, positions, edgelist=[edge], width=edge[2], ax=ax, arrows=False)
    else:
        nx.draw(G, positions, node_color=node_colors, node_size=55, with_labels=False, ax=ax, width=2.5, label=node_label)

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Consumer', markerfacecolor='LightBlue', markersize=15),
        Line2D([0], [0], marker='o', color='w', label='Internal Node', markerfacecolor='LightGrey', markersize=15),
        Line2D([0], [0], marker='o', color='w', label='Heat Source', markerfacecolor='Orange', markersize=15)
    ]
    plt.legend(handles=legend_elements, loc='best', frameon=False, fontsize=22)
    plt.box(False)
    fig0.savefig(title + ".svg")
    plt.close()

    return