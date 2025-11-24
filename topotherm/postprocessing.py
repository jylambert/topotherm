"""
Postprocessing of the results from the optimization. This includes the
calculation of the diameter and mass flow of the pipes, and the elimination of
unused pipes and nodes.

Functions
---------
- ``calc_diam_and_velocity`` : Calculate pipe diameter and velocity depending
  on the mass flow and pipe power.
- ``sts`` : Postprocessing for the single time step model.
- ``to_networkx_graph`` : Export the postprocessed, optimal district as a
  NetworkX graph.
- ``mts`` : Postprocessing for the multiple time step model.
"""

from typing import Tuple

import networkx as nx
import numpy as np
import pandas as pd
import pyomo.environ as pyo
from scipy.optimize import root

from topotherm.settings import Settings


def diameter_and_velocity(
    v: Tuple[float, float], mass_lin: float, settings: Settings
) -> Tuple[float, float]:
    """
    Equations for calculating the diameter and velocity of a pipe based on
    mass flow and power of the pipes

    Parameters
    ----------
    v : tuple of float
        Tuple containing the velocity and diameter (``(velocity, diameter)``).
    mass_lin : float
        Mass flow of the pipe (kg/s).
    settings : Settings
        Settings object containing water and piping parameters.

    Returns
    -------
    tuple of float
        Tuple containing:

        - ``velocity`` : float
            Calculated velocity of the pipe (m/s).
        - ``diameter`` : float
            Calculated diameter of the pipe (m).
    """
    vel, d = v
    reynolds = (settings.water.density * vel * d) / settings.water.dynamic_viscosity
    # friction factor
    f = (
        -1.8
        * np.log10((settings.piping.roughness / (3.7 * d)) ** 1.11 + 6.9 / reynolds)
    ) ** -2
    # eq. for diameter
    eq1 = vel - np.sqrt(
        (2 * settings.piping.max_pr_loss * d) / (f * settings.water.density)
    )
    # eq. for velocity
    eq2 = mass_lin - settings.water.density * vel * (np.pi / 4) * d**2
    return [eq1, eq2]


def calculate_hydraulics(
    power: np.ndarray, settings: Settings
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate mass flow, diameter, and velocity for each pipe given the
    installed thermal power and the supply/return temperatures from ``settings``.

    Parameters
    ----------
    power : np.ndarray
        Installed thermal power for each pipe.
    settings : Settings
        Settings object containing temperature setpoints and water properties.

    Returns
    -------
    tuple of np.ndarray
        Tuple containing:

        - ``mass_flow`` : np.ndarray
            Mass flow for each pipe (kg/s).
        - ``diameter`` : np.ndarray
            Pipe diameter (m).
        - ``velocity`` : np.ndarray
            Flow velocity (m/s).
    """

    delta_t = settings.temperatures.supply - settings.temperatures.return_

    # Compute mass flow rate for all pipes m = P / (cp * deltaT)
    m_dot = power * 1e3 / (settings.water.heat_capacity_cp * delta_t)

    # helper function for solving per pipe
    def solve_per_pipe(m):
        sol = root(
            lambda v: diameter_and_velocity(v, m, settings), (0.5, 0.02), method="lm"
        )
        return sol.x if sol.success else (np.nan, np.nan)

    # Vectorized solving
    results = np.array([solve_per_pipe(m) for m in m_dot])
    vel, d = results.T

    return m_dot, d, vel


def sts(model: pyo.ConcreteModel, matrices: dict, settings: Settings):
    """
    Postprocessing for the single time step model. This includes the
    calculation of the diameter and velocity of the pipes and the elimination
    of unused pipes and nodes.

    Parameters
    ----------
    model : pyo.ConcreteModel
        Solved Pyomo model.
    matrices : dict
        Dictionary containing the matrices.
    settings : tt.settings.Settings
        Settings for the optimization.

    Returns
    -------
    dict
        Optimal variables and postprocessed data.
    """

    # Get the values from the model
    p_ij = np.array(pyo.value(model.P["ij", "in", :, :]))
    p_ji = np.array(pyo.value(model.P["ji", "in", :, :]))
    p_source_inst = np.array(pyo.value(model.P_source_inst[:]))
    p_source = np.array(pyo.value(model.P_source[:, :]))

    # flow direction, binary
    lambda_ij = np.around(np.array(pyo.value(model.lambda_["ij", :])), 0)
    lambda_ji = np.around(np.array(pyo.value(model.lambda_["ji", :])), 0)
    lambda_sum = lambda_ij + lambda_ji

    q_c_opt = np.zeros([matrices["a_c"].shape[1], len(model.set_t)])
    flh_c_opt = np.zeros([matrices["a_c"].shape[1], len(model.set_t)])

    for d, e in model.cons:
        if d == "ij":
            # edge in incidence matrix where pipe exits into node n (==-1)
            a_i_idx = np.where(matrices["a_i"][:, e] == -1)
            # location where a_i_idx is connected to a_c
            a_c_idx = np.where(matrices["a_c"][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError("Error in the incidence matrix!")
            # assign the heat demand to the connected consumer if lambda is 1
            q_c_opt[a_c_idx[0], :] = lambda_ij[e] * matrices["q_c"][a_c_idx[0], :]
            flh_c_opt[a_c_idx[0], :] = (
                lambda_ij[e] * matrices["flh_consumer"][a_c_idx[0], :]
            )
        elif d == "ji":
            a_i_idx = np.where(matrices["a_i"][:, e] == 1)
            a_c_idx = np.where(matrices["a_c"][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError("Error in the incidence matrix!")
            q_c_opt[a_c_idx[0], :] = lambda_ji[e] * matrices["q_c"][a_c_idx[0], :]
            flh_c_opt[a_c_idx[0], :] = (
                lambda_ji[e] * matrices["flh_consumer"][a_c_idx[0], :]
            )

    # Remove nonzero elements row-wise
    q_c_opt = q_c_opt[q_c_opt.any(axis=1)]
    flh_c_opt = np.zeros([matrices["a_c"].shape[1], len(model.set_t)])

    # Postprocessing producers
    if np.shape(matrices["a_p"])[1] == 1:
        p_source_inst_opt = p_source_inst
        p_source_opt = p_source
        flh_s_opt = matrices["flh_source"]
    else:
        # clean up the sources to exclude 0 power sources
        p_source_inst_opt = p_source_inst[p_source_inst != 0]
        p_source_opt = p_source[p_source_inst != 0]
        flh_s_opt = matrices["flh_source"][p_source_inst != 0, :]

    # Adjust Incidence Matrix for further postprocessing
    for q, _ in enumerate(lambda_ij):
        # if not active, all is 0
        if lambda_ij[q] == 0 and lambda_ji[q] == 0:
            matrices["a_i"][:, q] = 0
            matrices["l_i"][q] = 0
        # if opposite direction operational, switch a_i with -1 and switch
        # values lambda_ij for ji. This is necessary for the postprocessing.
        elif lambda_ji[q] == 1:
            matrices["a_i"][:, q] = matrices["a_i"][:, q] * (-1)

    p_lin = p_ij + p_ji  # Power of the pipes

    # drop entries with 0 in the incidence matrix to reduce size
    valid_columns = matrices["a_i"].any(axis=0)
    valid_rows = matrices["a_i"].any(axis=1)

    p_lin_opt = p_lin[valid_columns]
    pos_opt = matrices["position"][valid_rows, :]
    a_c_opt = matrices["a_c"][valid_rows, :]
    a_c_opt = a_c_opt[:, a_c_opt.any(axis=0)]
    a_p_opt = matrices["a_p"][valid_rows, :]
    a_i_opt = matrices["a_i"][valid_rows, :][:, valid_columns]
    l_i_opt = matrices["l_i"][valid_columns]

    m_lin, d_lin, v_lin = calculate_hydraulics(p_lin_opt, settings)

    res = {
        "a_i": a_i_opt,
        "a_p": a_p_opt,
        "a_c": a_c_opt,
        "q_c": q_c_opt,
        "l_i": l_i_opt,
        "d_i_0": d_lin,
        "m_i_0": m_lin,
        "position": pos_opt,
        "p": p_lin_opt,
        "flh_c_opt": flh_c_opt,
        "flh_s_opt": flh_s_opt,
        "p_s_inst_opt": p_source_inst_opt,
        "p_s_opt": p_source_opt,
        "lambda_b_orig": lambda_sum,
    }

    return res


def mts(model: pyo.ConcreteModel, matrices: dict, settings: Settings) -> dict:
    """
    Postprocessing for the multiple time step model. This includes the
    calculation of the diameter and velocity of the pipes and the elimination
    of unused pipes and nodes.

    Parameters
    ----------
    model : pyo.ConcreteModel
        Solved Pyomo model.
    matrices : dict
        Dictionary containing the matrices.
    settings : tt.settings.Settings
        Settings for the optimization.

    Returns
    -------
    dict
        Optimal variables and postprocessed data.
    """

    # Get the values from the model
    p_ij = np.reshape(
        np.array(pyo.value(model.P["ij", "in", :, :])), (-1, matrices["q_c"].shape[1])
    )
    p_ji = np.reshape(
        np.array(pyo.value(model.P["ji", "in", :, :])), (-1, matrices["q_c"].shape[1])
    )
    p_cap = np.array(pyo.value(model.P_cap[:]))
    p_source_inst = np.array(pyo.value(model.P_source_inst[:]))
    p_source = np.array(pyo.value(model.P_source[:, :]))

    # flow direction, binary
    lambda_ij = np.reshape(
        np.around(np.array(pyo.value(model.lambda_["ij", :, :])), 0),
        (-1, matrices["q_c"].shape[1]),
    )
    lambda_ji = np.reshape(
        np.around(np.array(pyo.value(model.lambda_["ji", :, :])), 0),
        (-1, matrices["q_c"].shape[1]),
    )
    # built pipes
    lambda_b = np.around(np.array(pyo.value(model.lambda_b[:])), 0)

    q_c_opt = np.zeros([matrices["a_c"].shape[1], len(model.set_t)])
    flh_c_opt = np.zeros([matrices["a_c"].shape[1], len(model.set_t)])

    # Exclude non-connected consumers in Q_c, only affects the economic case
    # Check for consumers connected in direction ij
    for d, e in model.cons:
        if d == "ij":
            # edge in incidence matrix where pipe exits into node n (==-1)
            a_i_idx = np.where(matrices["a_i"][:, e] == -1)
            # location where a_i_idx is connected to a_c
            a_c_idx = np.where(matrices["a_c"][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError("Error in the incidence matrix!")
            # assign the heat demand to the connected consumer if lambda is 1
            q_c_opt[a_c_idx[0], :] = lambda_b[e] * matrices["q_c"][a_c_idx[0], :]
            flh_c_opt[a_c_idx[0], :] = (
                lambda_b[e] * matrices["flh_consumer"][a_c_idx[0], :]
            )
        elif d == "ji":
            a_i_idx = np.where(matrices["a_i"][:, e] == 1)
            a_c_idx = np.where(matrices["a_c"][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError("Error in the incidence matrix!")
            q_c_opt[a_c_idx[0], :] = lambda_b[e] * matrices["q_c"][a_c_idx[0], :]
            flh_c_opt[a_c_idx[0], :] = (
                lambda_b[e] * matrices["flh_consumer"][a_c_idx[0], :]
            )

    # Remove nonzero elements row-wise
    q_c_opt = q_c_opt[q_c_opt.any(axis=1)]
    flh_c_opt = flh_c_opt[flh_c_opt.any(axis=1)]

    # Postprocessing producers
    if np.shape(matrices["a_p"])[1] == 1:
        p_source_inst_opt = p_source_inst
        p_source_opt = p_source
        flh_s_opt = matrices["flh_source"]
    else:
        p_source_inst_opt = p_source_inst[p_source_inst != 0]
        p_source_opt = p_source[p_source_inst != 0, :]
        flh_s_opt = matrices["flh_source"][p_source_inst != 0, :]

    # Adaption of Incidence Matrix for further postprocessing
    for q in model.set_n_i:
        # if not active, all is 0
        if lambda_b[q] == 0:
            matrices["a_i"][:, q] = 0
            matrices["l_i"][q] = 0
        # if opposite direction operational, switch a_i with -1 and switch
        # values lambda_ij for ji. This is necessary for the postprocessing.
        elif (lambda_b[q] == 1) & (lambda_ji[q, 0] == 1):
            matrices["a_i"][:, q] = matrices["a_i"][:, q] * (-1)
            lambda_ij[q, np.where(lambda_ij[q, 1:] == 0)[0]] = 1
            lambda_ji[q, np.where(lambda_ji[q, 1:] == 1)[0]] = 0

    p_lin = p_cap  # Capacity of the pipes

    # drop entries with 0 in the incidence matrix to reduce size
    valid_columns = matrices["a_i"].any(axis=0)
    valid_rows = matrices["a_i"].any(axis=1)

    p_lin_opt = p_lin[valid_columns]
    p_ij_opt = p_ij[valid_columns, :]
    p_ji_opt = p_ji[valid_columns, :]
    lambda_ij_opt = lambda_ij[valid_columns, :]
    lambda_ji_opt = lambda_ji[valid_columns, :]
    pos_opt = matrices["position"][valid_rows, :]
    a_c_opt = matrices["a_c"][valid_rows, :]
    a_c_opt = a_c_opt[:, a_c_opt.any(axis=0)]
    a_p_opt = matrices["a_p"][valid_rows, :]
    a_i_opt = matrices["a_i"][valid_rows, :][:, valid_columns]
    l_i_opt = matrices["l_i"][valid_columns]

    m_lin, d_lin, v_lin = calculate_hydraulics(p_lin_opt, settings)

    res = {
        "a_i": a_i_opt,
        "a_p": a_p_opt,
        "a_c": a_c_opt,
        "q_c": q_c_opt,
        "l_i": l_i_opt,
        "lambda_b_orig": lambda_b,
        "lambda_ij_opt": lambda_ij_opt,
        "lambda_ji_opt": lambda_ji_opt,
        "d_i_0": d_lin,
        "m_i_0": m_lin,
        "position": pos_opt,
        "p": p_lin_opt,
        "p_ij": p_ij_opt,
        "p_ji": p_ji_opt,
        "flh_c_opt": flh_c_opt,
        "flh_s_opt": flh_s_opt,
        "p_s_inst_opt": p_source_inst_opt,
        "p_s_opt": p_source_opt,
    }

    return res


def to_networkx_graph(matrices: dict) -> nx.DiGraph:
    """
    Convert incidence matrices to a directed NetworkX graph.
    Includes the nodes and edges of the district, their length, and, if
    available, installed diameter, and power.

    Parameters
    ----------
    matrices : dict
        Dictionary containing the incidence matrices as output by ``topotherm.postprocessing.sts`` or as input to the model.

    Returns
    -------
    nx.DiGraph
        NetworkX directed graph.
    """
    G = nx.DiGraph()

    # Add the nodes to the graph
    sums = matrices["a_c"].sum(axis=1)
    prod = matrices["a_p"].T.sum(axis=0)
    ges = (sums + prod).flatten()
    ges = (sums + prod).flatten()

    positions = matrices["positions"]  # Avoid redundant indexing
    colors = np.where(
        ges == 1,
        "Red",
        np.where(ges == 0, "Green", np.where(ges <= -1, "Orange", None)),
    )
    types = np.where(
        ges == 1,
        "consumer",
        np.where(ges == 0, "internal", np.where(ges <= -1, "source", None)),
    )

    for q in range(len(ges)):
        if colors[q]:  # Skip nodes that don't match any category
            G.add_node(
                q, color=colors[q], type_=types[q], x=positions[q, 0], y=positions[q, 1]
            )

    sources = np.argmax(
        matrices["a_i"] == 1, axis=0
    )  # First occurrence of 1 in each column
    targets = np.argmax(
        matrices["a_i"] == -1, axis=0
    )  # First occurrence of -1 in each column

    # Get existing to allow solved and unsolved incidence matrices
    _data = ["l_i", "d_i_0", "p"]
    available_data = [k for k in _data if k in matrices.keys()]
    available_matrices = [matrices[k].flatten() for k in available_data]
    edges = np.column_stack([sources, targets] + available_matrices)

    # hacky solution, should maybe work directly with networkx functions
    # TODO
    for row in edges:
        u, v = row[0].astype(int), row[1].astype(int)
        data = {}
        idx = 2  # start after u, v

        for key in available_data:
            val = row[idx]
            idx += 1

            if key == "l_i":
                data["l"] = val
            elif key == "d_i_0":
                data["d"] = val
            elif key == "p":
                data["p"] = val
                if val == 0:  # Skip edges with p == 0
                    break
        else:
            # Only add edge if loop didn’t break (i.e., p != 0 or p not present)
            G.add_edge(u, v, **data)

    return G


def to_dataframe(
    matrices_optimal: dict, matrices_init: dict
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Export the postprocessed, optimal district as pandas DataFrames.
    Includes the nodes and edges of the district, their length, installed
    diameter, and power.

    Parameters
    ----------
    matrices_optimal : dict
        Solved and cleaned matrices as output by
        ``topotherm.postprocessing.sts``.
    matrices_init : dict
        Original matrices as output by ``topotherm.fileio.load``.

    Returns
    -------
    nodes : pd.DataFrame
        DataFrame of nodes.
    edges : pd.DataFrame
        DataFrame of edges.
    """
    # Create a DataFrame for the nodes
    nodes = pd.DataFrame(
        index=range(matrices_optimal["a_i"].shape[0]),
        columns=["type_", "x", "y", "demand", "total_installed_power"],
        data=None,
    )

    # Calculate the sum of matrices
    sum_ac = matrices_optimal["a_c"].sum(axis=1)
    sum_ap = matrices_optimal["a_p"].T.sum(axis=0)
    total_sum = np.array(sum_ac + sum_ap).flatten()

    # Assign types based on the sum
    nodes["type_"] = [
        "consumer" if x == 1 else "internal" if x == 0 else "source" for x in total_sum
    ]
    nodes["x"] = matrices_optimal["position"][:, 0]
    nodes["y"] = matrices_optimal["position"][:, 1]

    # TODO: inherit in some way the consumer and producer id so that you don't have to search for it
    # in the original matrices, reducing the need to import them.
    sources_nodes = nodes.type_ == "source"
    positions_sources = nodes.loc[sources_nodes, ["x", "y"]].values
    matches = np.all(
        matrices_init["position"][:, None, :] == positions_sources[None, :, :], axis=2
    )
    original_source_nodes = np.where(matches)[0]
    original_source_prods = np.where(
        matrices_init["a_p"][original_source_nodes, :] == -1
    )[1]
    # assume that we want to write out the total sum of installed power for each source
    # inherent limitation of the dataframe structure, only one dimensional data possible
    nodes.loc[sources_nodes, "total_installed_power"] = matrices_optimal[
        "p_s_inst_opt"
    ][original_source_prods].sum()

    consumer_nodes = nodes[nodes.type_ == "consumer"].index
    positions_consumers = nodes.loc[consumer_nodes, ["x", "y"]].values
    matches = np.all(
        matrices_init["position"][:, None, :] == positions_consumers[None, :, :], axis=2
    )
    original_consumer_nodes = np.where(matches)[0]
    original_consumer_edges = np.where(
        matrices_init["a_c"][original_consumer_nodes, :] == 1
    )[1]
    nodes.loc[consumer_nodes, "demand"] = (
        (matrices_init["q_c"] * matrices_init["flh_consumer"])[original_consumer_edges]
        .sum(axis=1)
        .squeeze()
    )
    nodes.loc[consumer_nodes, "total_installed_power"] = (
        (matrices_init["q_c"])[original_consumer_edges].max(axis=1).squeeze()
    )
    if (
        abs(
            nodes.loc[consumer_nodes, "total_installed_power"].sum()
            - matrices_optimal["q_c"].max(axis=1).sum()
        )
        > 1e-3
    ):
        raise ValueError(
            f'Error in the incidence matrix! Demand {nodes.demand.sum()} != total demand {matrices_optimal["q_c"].sum()}'
        )

    # Create a DataFrame for the edges
    edges = pd.DataFrame()
    edges["start_node"] = np.argmax(matrices_optimal["a_i"] == 1, axis=0)
    edges["end_node"] = np.argmax(matrices_optimal["a_i"] == -1, axis=0)
    edges["x_start"] = matrices_optimal["position"][edges["start_node"], 0]
    edges["y_start"] = matrices_optimal["position"][edges["start_node"], 1]
    edges["x_end"] = matrices_optimal["position"][edges["end_node"], 0]
    edges["y_end"] = matrices_optimal["position"][edges["end_node"], 1]
    edges["length"] = matrices_optimal["l_i"]
    edges["diameter"] = matrices_optimal["d_i_0"]
    edges["power"] = matrices_optimal["p"]
    edges.loc[:, "to_consumer"] = edges["end_node"].map(nodes["type_"]).eq("consumer")
    edges.loc[:, "from_consumer"] = (
        edges["start_node"].map(nodes["type_"]).eq("consumer")
    )
    if edges.to_consumer.sum() + edges.from_consumer.sum() != len(
        matrices_optimal["q_c"]
    ):
        raise ValueError(
            f'Error in the incidence matrix! To consumer {edges.to_consumer.sum()} + from_consumer {edges.from_consumer.sum()} != total consumers {len(matrices_optimal["q_c"])}'
        )
    if edges.from_consumer.sum() > 0:
        raise ValueError(
            f"Error in the incidence matrix! From consumer {edges.from_consumer.sum()} is not 0 for single time steps"
        )
    return nodes, edges
