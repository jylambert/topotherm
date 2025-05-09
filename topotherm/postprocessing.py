"""Postprocessing of the results from the optimization. This includes the 
calculation of the diameter and mass flow of the pipes, the elimination of
unused pipes and nodes.

This module includes the following functions:
    * calc_diam_and_velocity: Equations for the calculation of the diameter and velocity of the pipes depending on the mass flow and the power of the pipes
    * sts: Postprocessing for the single time step model
    * to_networkx_graph: Export the postprocessed, optimal district as a networkx graph
    * mts: Postprocessing for the multiple time step model

"""

from typing import Tuple

import numpy as np
import pyomo.environ as pyo
from scipy.optimize import root
import networkx as nx

from topotherm.settings import Settings


def calc_diam_and_velocity(
        v: Tuple[float, float],
        mass_lin: float,
        settings: Settings) -> Tuple[float, float]:
    """Equations for the calculation of the diameter and velocity of
    the pipes depending on the mass flow and the power of the pipes.

    Args:
        v (Tuple[float, float]): Tuple containing the velocity and diameter
        mass_lin (float): mass flow of the pipe
        settings (Settings): settings for the optimization
    
    Returns:
        Tuple[float, float]: Tuple containing the velocity and diameter
    """
    vel, d = v
    reynolds = ((settings.water.density * vel * d)
                / settings.water.dynamic_viscosity)
    # friction factor
    f = (-1.8 * np.log10((settings.piping.roughness / (3.7 * d)) ** 1.11
                            + 6.9 / reynolds)
                            )**-2
    # eq. for diameter
    eq1 = vel - np.sqrt((2 * settings.piping.max_pr_loss * d)
                        / (f * settings.water.density))
    # eq. for velocity
    eq2 = mass_lin - settings.water.density * vel * (np.pi / 4) * d ** 2
    return [eq1, eq2]


def sts(model: pyo.ConcreteModel,
        matrices: dict,
        settings: Settings) -> dict:
    """Postprocessing for the single time step model. This includes the
    calculation of the diameter and velocity of the pipes, the elimination of
    unused pipes and nodes. Drops unused pipes and nodes.

    Args:
        model (pyo.ConcreteModel): solved pyomo model
        matrices (dict): dict containing the matrices
        settings (tt.settings.Settings): settings for the optimization
    
    Returns:
        dict: Optimal variables and postprocessed data
    """
    # Get the values from the model
    p_ij = np.array(pyo.value(model.P['ij', 'in', :, :]))
    p_ji = np.array(pyo.value(model.P['ji', 'in', :, :]))
    p_source_inst = np.array(pyo.value(model.P_source_inst[:]))
    p_source = np.array(pyo.value(model.P_source[:, :]))

    # flow direction, binary
    lambda_ij = np.around(np.array(pyo.value(model.lambda_['ij', :])), 0)
    lambda_ji = np.around(np.array(pyo.value(model.lambda_['ji', :])), 0)
    lambda_sum = lambda_ij + lambda_ji

    q_c_opt = np.zeros([matrices['a_c'].shape[1], len(model.set_t)])
    flh_c_opt = np.zeros([matrices['a_c'].shape[1], len(model.set_t)])

    # Exclude non-connected consumers in Q_c, only affects the economic case
    # Check for consumers connected in direction ij
    for d, e in model.cons:
        if d == 'ij':
            # edge in incidence matrix where pipe exits into node n (==-1)
            a_i_idx = np.where(matrices['a_i'][:, e] == -1)
            # location where a_i_idx is connected to a_c
            a_c_idx = np.where(matrices['a_c'][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError('Error in the incidence matrix!')
            # assign the heat demand to the connected consumer if lambda is 1
            q_c_opt[a_c_idx[0], :] = lambda_ij[e] * matrices['q_c'][a_c_idx[0], :]
            flh_c_opt[a_c_idx[0], :] = lambda_ij[e] * matrices['flh_consumer'][a_c_idx[0], :]
        elif d == 'ji':
            a_i_idx = np.where(matrices['a_i'][:, e] == 1)
            a_c_idx = np.where(matrices['a_c'][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError('Error in the incidence matrix!')
            q_c_opt[a_c_idx[0], :] = lambda_ji[e] * matrices['q_c'][a_c_idx[0], :]
            flh_c_opt[a_c_idx[0], :] = lambda_ji[e] * matrices['flh_consumer'][a_c_idx[0], :]

    # Remove nonzero elements row-wise
    q_c_opt = q_c_opt[q_c_opt.any(axis=1)]
    flh_c_opt = np.zeros([matrices['a_c'].shape[1], len(model.set_t)])

    # Postprocessing producers
    if np.shape(matrices['a_p'])[1] == 1:
        p_source_inst_opt = p_source_inst
        p_source_opt = p_source
        flh_s_opt = matrices['flh_source']
    else:
        # clean up the sources to exclude 0 power sources
        p_source_inst_opt = p_source_inst[p_source_inst != 0]
        p_source_opt = p_source[p_source_inst != 0]
        flh_s_opt = matrices['flh_source'][p_source_inst != 0, :]

    # Adjust Incidence Matrix for further postprocessing
    for q, _ in enumerate(lambda_ij):
        # if not active, all is 0
        if lambda_ij[q] == 0 and lambda_ji[q] == 0:
            matrices['a_i'][:, q] = 0
            matrices['l_i'][q] = 0
        # if opposite direction operational, switch a_i with -1 and switch
        # values lambda_ij for ji. This is necessary for the postprocessing.
        elif lambda_ji[q] == 1:
            matrices['a_i'][:, q] = matrices['a_i'][:, q] * (-1)

    p_lin = p_ij + p_ji  # Power of the pipes

    # drop entries with 0 in the incidence matrix to reduce size
    valid_columns = matrices['a_i'].any(axis=0)
    valid_rows = matrices['a_i'].any(axis=1)

    p_lin_opt = p_lin[valid_columns]
    pos_opt = matrices['position'][valid_rows, :]
    a_c_opt = matrices['a_c'][valid_rows, :]
    a_c_opt = a_c_opt[:, a_c_opt.any(axis=0)]
    a_p_opt = matrices['a_p'][valid_rows, :]
    a_i_opt = matrices['a_i'][valid_rows, :][:, valid_columns]
    l_i_opt = matrices['l_i'][valid_columns]

    a_i_shape_opt = np.shape(a_i_opt)   # (rows 0, columns 1)
    d_lin = np.zeros(a_i_shape_opt[1])  # Initialize linear diameters
    v_lin = np.zeros(a_i_shape_opt[1])  # Initialize velocities

    # Assign supply and return temperatures
    supply_temp_opt = np.ones(a_i_shape_opt[1]) * settings.temperatures.supply
    return_temp_opt = np.ones(a_i_shape_opt[1]) * settings.temperatures.return_

    # Calculate the mass flow for each pipe with m cp deltaT = P
    m_lin = (p_lin_opt * 1000
             / (settings.water.heat_capacity_cp
                * (supply_temp_opt - return_temp_opt)
                ))

    # Calculate the diameter and velocity for each pipe
    for h in range(a_i_shape_opt[1]):
        mass_lin = m_lin[h]
        sol = root(lambda v: calc_diam_and_velocity(v, mass_lin, settings),
                   (0.5, 0.02),
                   method='lm')
        if sol.success:
            v_lin[h], d_lin[h] = sol.x
        else:
            print(h, 'Warning: Failed to calculate diameter and velocity!')

    res = dict(
        a_i=a_i_opt,
        a_p=a_p_opt,
        a_c=a_c_opt,
        q_c=q_c_opt,
        l_i=l_i_opt,
        d_i_0=d_lin,
        m_i_0=m_lin,
        position=pos_opt,
        p=p_lin_opt,
        flh_c_opt=flh_c_opt,
        flh_s_opt=flh_s_opt,
        p_s_inst_opt=p_source_inst_opt,
        p_s_opt=p_source_opt,
        lambda_b_orig=lambda_sum
    )

    return res


def mts(model: pyo.ConcreteModel,
        matrices: dict,
        settings: Settings) -> dict:
    """Postprocessing for the multiple time step model. This includes the
    calculation of the diameter and velocity of the pipes, the elimination of
    unused pipes and nodes. Drops unused pipes and nodes.

    Args:
        model (pyo.ConcreteModel): solved pyomo model
        matrices (dict): dict containing the matrices
        settings (tt.settings.Settings): settings for the optimization

    Returns:
        dict: Optimal variables and postprocessed data
    """
    # Get the values from the model
    p_ij = np.reshape(np.array(pyo.value(model.P['ij', 'in', :, :])), (-1, matrices['q_c'].shape[1]))
    p_ji = np.reshape(np.array(pyo.value(model.P['ji', 'in', :, :])), (-1, matrices['q_c'].shape[1]))
    p_cap = np.array(pyo.value(model.P_cap[:]))
    p_source_inst = np.array(pyo.value(model.P_source_inst[:]))
    p_source = np.array(pyo.value(model.P_source[:, :]))

    # flow direction, binary
    lambda_ij = np.reshape(np.around(np.array(pyo.value(model.lambda_['ij', :, :])), 0), (-1, matrices['q_c'].shape[1]))
    lambda_ji = np.reshape(np.around(np.array(pyo.value(model.lambda_['ji', :, :])), 0), (-1, matrices['q_c'].shape[1]))
    lambda_b = np.around(np.array(pyo.value(model.lambda_b[:])), 0)

    q_c_opt = np.zeros([matrices['a_c'].shape[1], len(model.set_t)])
    flh_c_opt = np.zeros([matrices['a_c'].shape[1], len(model.set_t)])

    # Exclude non-connected consumers in Q_c, only affects the economic case
    # Check for consumers connected in direction ij
    for d, e in model.cons:
        if d == 'ij':
            # edge in incidence matrix where pipe exits into node n (==-1)
            a_i_idx = np.where(matrices['a_i'][:, e] == -1)
            # location where a_i_idx is connected to a_c
            a_c_idx = np.where(matrices['a_c'][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError('Error in the incidence matrix!')
            # assign the heat demand to the connected consumer if lambda is 1
            q_c_opt[a_c_idx[0], :] = lambda_b[e] * matrices['q_c'][a_c_idx[0], :]
            flh_c_opt[a_c_idx[0], :] = lambda_b[e] * matrices['flh_consumer'][a_c_idx[0], :]
        elif d == 'ji':
            a_i_idx = np.where(matrices['a_i'][:, e] == 1)
            a_c_idx = np.where(matrices['a_c'][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError('Error in the incidence matrix!')
            q_c_opt[a_c_idx[0], :] = lambda_b[e] * matrices['q_c'][a_c_idx[0], :]
            flh_c_opt[a_c_idx[0], :] = lambda_b[e] * matrices['flh_consumer'][a_c_idx[0], :]

    # Remove nonzero elements row-wise
    q_c_opt = q_c_opt[q_c_opt.any(axis=1)]
    flh_c_opt = flh_c_opt[flh_c_opt.any(axis=1)]

    # Postprocessing producers
    if np.shape(matrices['a_p'])[1] == 1:
        p_source_inst_opt = p_source_inst
        p_source_opt = p_source
        flh_s_opt = matrices['flh_source']
    else:
        p_source_inst_opt = p_source_inst[p_source_inst != 0]
        p_source_opt = p_source[p_source_inst != 0, :]
        flh_s_opt = matrices['flh_source'][p_source_inst != 0, :]

    # Adaption of Incidence Matrix for further postprocessing
    for q in model.set_n_i:
        # if not active, all is 0
        if lambda_b[q] == 0:
            matrices['a_i'][:, q] = 0
            matrices['l_i'][q] = 0
        # if opposite direction operational, switch a_i with -1 and switch
        # values lambda_ij for ji. This is necessary for the postprocessing.
        elif (lambda_b[q] == 1) & (lambda_ji[q, 0] == 1):
            matrices['a_i'][:, q] = matrices['a_i'][:, q] * (-1)
            lambda_ij[q, np.where(lambda_ij[q, 1:] == 0)[0]] = 1
            lambda_ji[q, np.where(lambda_ji[q, 1:] == 1)[0]] = 0


    p_lin = p_cap  # Capacity of the pipes

    # drop entries with 0 in the incidence matrix to reduce size
    valid_columns = matrices['a_i'].any(axis=0)
    valid_rows = matrices['a_i'].any(axis=1)

    p_lin_opt = p_lin[valid_columns]
    p_ij_opt = p_ij[valid_columns, :]
    p_ji_opt = p_ji[valid_columns, :]
    lambda_ij_opt = lambda_ij[valid_columns, :]
    lambda_ji_opt = lambda_ji[valid_columns, :]
    pos_opt = matrices['position'][valid_rows, :]
    a_c_opt = matrices['a_c'][valid_rows, :]
    a_c_opt = a_c_opt[:, a_c_opt.any(axis=0)]
    a_p_opt = matrices['a_p'][valid_rows, :]
    a_i_opt = matrices['a_i'][valid_rows, :][:, valid_columns]
    l_i_opt = matrices['l_i'][valid_columns]

    a_i_shape_opt = np.shape(a_i_opt)  # (rows 0, columns 1)
    d_lin = np.zeros(a_i_shape_opt[1])  # Initialize linear diameters
    v_lin = np.zeros(a_i_shape_opt[1])  # Initialize velocities

    # Assign supply and return temperatures
    supply_temp_opt = np.ones(a_i_shape_opt[1]) * settings.temperatures.supply
    return_temp_opt = np.ones(a_i_shape_opt[1]) * settings.temperatures.return_

    # Calculate the mass flow for each pipe with m cp deltaT = P
    m_lin = (p_lin_opt * 1000
             / (settings.water.heat_capacity_cp
                * (supply_temp_opt - return_temp_opt)
                ))

    # Calculate the diameter and velocity for each pipe
    for h in range(a_i_shape_opt[1]):
        mass_lin = m_lin[h]
        sol = root(lambda v: calc_diam_and_velocity(v, mass_lin, settings),
                   (0.5, 0.02),
                   method='lm')
        if sol.success:
            v_lin[h], d_lin[h] = sol.x
        else:
            print(h, 'Warning: Failed to calculate diameter and velocity!')

    res = dict(
        a_i=a_i_opt,
        a_p=a_p_opt,
        a_c=a_c_opt,
        q_c=q_c_opt,
        l_i=l_i_opt,
        lambda_b_orig=lambda_b,
        lambda_ij_opt=lambda_ij_opt,
        lambda_ji_opt=lambda_ji_opt,
        d_i_0=d_lin,
        m_i_0=m_lin,
        position=pos_opt,
        p=p_lin_opt,
        p_ij=p_ij_opt,
        p_ji=p_ji_opt,
        flh_c_opt=flh_c_opt,
        flh_s_opt=flh_s_opt,
        p_s_inst_opt=p_source_inst_opt,
        p_s_opt=p_source_opt
    )

    return res


def to_networkx_graph(matrices: dict) -> nx.DiGraph:
    """Export the postprocessed, optimal district as a networkx graph.
    Includes the nodes and edges of the district, their length, installed
    diameter and power.

    Args:
        matrices: a dict containing the following keys:
        - a_i (internal matrix)
        - a_p (producer matrix)
        - a_c (consumer matrix)
        - q_c (heat demand of the connected consumers)
        - l_i (of the pipes)
        - position (positions of the nodes)
        - d_i_0 (diameters of the optimal pipes)
        - m_i_0 (mass flow of the optimal pipes)
        - p (Power of the optimal pipes)
        
    Returns:
        nx.DiGraph: networkx graph
    """
    G = nx.DiGraph()

    # Add the nodes to the graph
    sums = matrices['a_c'].sum(axis=1)
    prod = matrices['a_p'].T.sum(axis=0)
    ges = sums + prod

    ges = np.array(ges).flatten()

    for q in range(matrices['a_c'].shape[0]):
        x, y = matrices['position'][q, 0], matrices['position'][q, 1]
        if ges[q] == 1:
            G.add_node(q, color='Red', type_='consumer', x=x, y=y)
        elif ges[q] == 0:
            G.add_node(q, color='Green', type_='internal', x=x, y=y)
        if ges[q] <= -1:
            G.add_node(q, color='Orange', type_='source', x=x, y=y)

    # edge_labels = dict()
    # Add the edges to the graph
    for k in range(matrices['a_i'].shape[1]):
        s = (np.where(matrices['a_i'][:, k] == 1)[0][0],
             np.where(matrices['a_i'][:, k] == -1)[0][0])   
        G.add_edge(s[0], s[1],
                   weight=matrices['l_i'][k].item(),  # important: float
                   d=matrices['d_i_0'][k],
                   p=matrices['p'][k])

    # drop all edges with p=0
    G.remove_edges_from([(u, v) for u, v, d in G.edges(data=True) if d['p'] == 0])
    return G
