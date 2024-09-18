"""Postprocessing of the results from the optimization. This includes the calculation of the
diameter and mass flow of the pipes, the elimination of unused pipes and nodes.

This module includes the following functions:
    * sts: Postprocessing for the STS model
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
        settings: Settings):
    """Postprocessing for the single time step model. This includes the
    calculation of the diameter and velocity of the pipes, the elimination of
    unused pipes and nodes.

    Args:
        model (pyo.ConcreteModel): solved pyomo model
        matrices (dict): dict containing the matrices
        settings (tt.settings.Settings): settings for the optimization
    
    Returns:
        _dict: containing the variables and postprocessed data
    """
    # Get the values from the model
    p_ij = np.array(pyo.value(model.P['ij', 'in', :, :]))
    p_ji = np.array(pyo.value(model.P['ji', 'in', :, :]))
    # flow direction, binary
    lambda_ij = np.around(np.array(pyo.value(model.lambda_['ij', :])), 0)
    lambda_ji = np.around(np.array(pyo.value(model.lambda_['ji', :])), 0)

    q_c_opt = np.zeros([matrices['a_c'].shape[1], len(model.set_t)])

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
        elif d == 'ji':
            a_i_idx = np.where(matrices['a_i'][:, e] == 1)
            a_c_idx = np.where(matrices['a_c'][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError('Error in the incidence matrix!')
            q_c_opt[a_c_idx[0], :] = lambda_ji[e] * matrices['q_c'][a_c_idx[0], :]

    # Remove nonzero elements row-wise
    q_c_opt = q_c_opt[q_c_opt.any(axis=1)]

    # Adaption of Incidence Matrix for further postprocessing
    for q, _ in enumerate(lambda_ij):
        # if not active, all is 0
        if lambda_ij[q] == 0 and lambda_ji[q] == 0:
            matrices['a_i'][:, q] = 0
            matrices['l_i'][q] = 0
        # if opposite direction operational, switch a_i with -1 and switch
        # values lambda_ij for ji. This is necessary for the postprocessing.
        elif lambda_ji[q] == 1:
            matrices['a_i'][:, q] = matrices['a_i'][:, q] * (-1)
            lambda_ij[q] = 1    # @TODO: For sts this can be removed
            lambda_ji[q] = 0    # @TODO: For sts this can be removed

    p_lin = p_ij + p_ji  # Power of the pipes

    # drop entries with 0 in the incidence matrix to reduce size
    valid_columns = matrices['a_i'].any(axis=0)
    valid_rows = matrices['a_i'].any(axis=1)

    p_lin_opt = p_lin[valid_columns]
    pos_opt = matrices['position'][valid_rows, :]
    a_c_opt = matrices['a_c'][valid_rows, :]
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
            print(h, 'failed to calculate diameter and velocity!')

    res = dict(
        a_i=a_i_opt,
        a_p=a_p_opt,
        a_c=a_c_opt,
        q_c=q_c_opt,
        l_i=l_i_opt,
        lambda_ij_opt=lambda_ij,  #@TODO: Remove lambda_ij_opt and lambda_ji_opt? -> lambda_ij_opt should only
        lambda_ji_opt=lambda_ji,  #@TODO: contain 1 and lamdba_ji_opt only 0 (for sts!).
        d_i_0=d_lin,
        m_i_0=m_lin,
        position=pos_opt,
        p=p_lin_opt
    )

    return res


def mts(model: pyo.ConcreteModel,
        matrices: dict,
        settings: Settings):
    """Postprocessing for the multiple time step model. This includes the
    calculation of the diameter and velocity of the pipes, the elimination of
    unused pipes and nodes.

    Args:
        model (pyo.ConcreteModel): solved pyomo model
        matrices (dict): dict containing the matrices
        settings (tt.settings.Settings): settings for the optimization

    Returns:
        _dict: containing the variables and postprocessed data
    """
    # Get the values from the model
    #p_ij = np.reshape(np.array(pyo.value(model.P['ij', 'in', :, :])), (-1, matrices['q_c'].shape[1]))
    p_ij = np.reshape(np.array(pyo.value(model.P['ij', 'in', :, :])), (-1, matrices['q_c'].shape[1]))
    #p_ji = np.array(pyo.value(model.P['ji', 'in', :, :]))
    p_ji = np.reshape(np.array(pyo.value(model.P['ji', 'in', :, :])), (-1, matrices['q_c'].shape[1]))
    p_cap = np.array(pyo.value(model.P_cap[:]))

    # flow direction, binary
    #lambda_ij = np.around(np.array(pyo.value(model.lambda_['ij', :, :])), 0)
    lambda_ij = np.reshape(np.around(np.array(pyo.value(model.lambda_['ij', :, :])), 0), (-1, matrices['q_c'].shape[1]))
    #lambda_ji = np.around(np.array(pyo.value(model.lambda_['ji', :, :])), 0)
    lambda_ji = np.reshape(np.around(np.array(pyo.value(model.lambda_['ji', :, :])), 0), (-1, matrices['q_c'].shape[1]))
    lambda_b = np.around(np.array(pyo.value(model.lambda_b[:])), 0)

    q_c_opt = np.zeros([matrices['a_c'].shape[1], len(model.set_t)])

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
        elif d == 'ji':
            a_i_idx = np.where(matrices['a_i'][:, e] == 1)
            a_c_idx = np.where(matrices['a_c'][a_i_idx[0], :][0] == 1)
            if len(a_i_idx) != 1 or len(a_c_idx) != 1:
                raise ValueError('Error in the incidence matrix!')
            q_c_opt[a_c_idx[0], :] = lambda_b[e] * matrices['q_c'][a_c_idx[0], :]

    # Remove nonzero elements row-wise
    q_c_opt = q_c_opt[q_c_opt.any(axis=1)]

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
            print(h, 'failed to calculate diameter and velocity!')

    res = dict(
        a_i=a_i_opt,
        a_p=a_p_opt,
        a_c=a_c_opt,
        q_c=q_c_opt,
        l_i=l_i_opt,
        lambda_ij_opt=lambda_ij_opt,
        lambda_ji_opt=lambda_ji_opt,
        d_i_0=d_lin,
        m_i_0=m_lin,
        position=pos_opt,
        p=p_lin_opt,
        p_ij=p_ij_opt,
        p_ji=p_ji_opt
    )

    return res


def to_networkx_graph(matrices):
    """Export the postprocessed, optimal district as a networkx graph. 

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
        
    Returns: Figure of the district
    """
    G = nx.DiGraph()
    s = np.array([0, 0, 0])

    # Add the nodes to the graph
    sums = matrices['a_c'].sum(axis=1)
    prod = matrices['a_p'].T.sum(axis=0)
    ges = sums + prod

    ges = np.array(ges).flatten()

    for q in range(matrices['a_c'].shape[0]):
        x, y = matrices['position'][q, 0], matrices['position'][q, 1]
        if ges[q] == 1:
            G.add_node(q, color='Red', type='consumer', x=x, y=y)
        elif ges[q] == 0:
            G.add_node(q, color='Green', type='internal', x=x, y=y)
        if ges[q] == -1:
            G.add_node(q, color='Orange', type='source', x=x, y=y)

    # edge_labels = dict()
    # Add the edges to the graph
    for k in range(matrices['a_i'].shape[1]):
        s = (np.where(matrices['a_i'][:, k] == 1)[0][0],
             np.where(matrices['a_i'][:, k] == -1)[0][0])   
        G.add_edge(s[0], s[1],
                   weight=matrices['l_i'][k],
                   d=matrices['d_i_0'][k],
                   p=matrices['p'][k])

    # drop all edges with p=0
    G.remove_edges_from([(u, v) for u, v, d in G.edges(data=True) if d['p'] == 0])
    return G
