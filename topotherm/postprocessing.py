"""Postprocessing of the results from the optimization. This includes the calculation of the
diameter and mass flow of the pipes, the elimination of unused pipes and nodes.

This module includes the following functions:
    * sts: Postprocessing for the STS model
    * mts: Postprocessing for the MTS model
"""

import numpy as np
import pyomo.environ as pyo
from scipy.optimize import root

from topotherm import settings


def postprocess(model, matrices, sets, mode, t_supply, t_return):
    """Create variables for the thermo-hydraulic coupled optimization.

    Args:
        model (pyo.ConcreteModel): pyomo model
        matrices (dict): dict containing the matrices
        sets (dict): dict containing the sets
        mode (str): sts or mts
        t_supply (float): supply temperature
        t_return (float): return temperature

    
    Returns:
        _type_: dict containing the variables
    """
    data_dict = {}

    p_21 = np.zeros([sets['a_i_shape'][1], len(model.set_t)])
    p_11 = np.zeros([sets['a_i_shape'][1], len(model.set_t)])
    lambda_dir_1 = np.zeros([sets['a_i_shape'][1], len(model.set_t)])
    lambda_dir_2 = np.zeros([sets['a_i_shape'][1], len(model.set_t)])
    lambda_built = np.zeros(sets['a_i_shape'][1])
    p_source_built = np.zeros(sets['a_p_shape'][1])
    p_cap = np.zeros(sets['a_i_shape'][1])

    for v in model.component_objects(pyo.Var, active=True):
        var_dict = {(v.name, index): pyo.value(v[index]) for index in v}
        data_dict.update(var_dict)
        if v.name == 'P_21':
            for index in v:
                p_21[index] = pyo.value(v[index])
        if v.name == 'P_11':
            for index in v:
                p_11[index] = pyo.value(v[index])
        if v.name == 'lambda_dir_1':
            for index in v:
                lambda_dir_1[index] = pyo.value(v[index])
        if v.name == 'lambda_dir_2':
            for index in v:
                lambda_dir_2[index] = pyo.value(v[index])
        if v.name == "lambda_built":
            for index in v:
                lambda_built[index] = pyo.value(v[index])
        if v.name == 'P_source_cap':
            for index in v:
                p_source_built[index] = pyo.value(v[index])
        if v.name == "P_cap":
            for index in v:
                p_cap[index] = pyo.value(v[index])

    # Round lamdda_dir and lambda_built to make sure that hey are integer
    lambda_dir_1 = np.around(lambda_dir_1, 0)
    lambda_dir_2 = np.around(lambda_dir_2, 0)
    lambda_built = np.around(lambda_built, 0)

    q_c_real = np.zeros([sets['a_c_shape'][1], len(model.set_t)])

    # Exclude non-connected consumers in Q_c, only affects the economic case
    # Check for consumers connected in direction ij
    for i in sets['connection_c_ij']:
        q_c_real[np.where(matrices['a_c'][np.where(matrices['a_i'][:, i] == -1)[0], :][0] == 1)[0], :] = \
            lambda_dir_1[i, 0] * matrices['q_c'][np.where(matrices['a_c'][np.where(matrices['a_i'][:, i] == -1)[0], :][0] == 1)[0], :]

    # Check for consumers connected in direction ji
    for i in sets['connection_c_ji']:
        q_c_real[np.where(matrices['a_c'][np.where(matrices['a_i'][:, i] == 1)[0], :][0] == 1)[0], :] = \
            lambda_dir_2[i, 0] * matrices['q_c'][np.where(matrices['a_c'][np.where(matrices['a_i'][:, i] == 1)[0], :][0] == 1)[0], :]

    q_c_real = q_c_real[q_c_real.any(axis=1)]     # Remove nonzero elements row-wise

    # Restart, Adaption of Incidence Matrix for the thermo-hydraulic coupled optimization
    if mode == "sts":
        for q, _ in enumerate(lambda_dir_1):
            if lambda_dir_1[q] == 0 and lambda_dir_2[q] == 0:
                matrices['a_i'][:, q] = 0
                matrices['l_i'][q] = 0
            elif lambda_dir_2[q] == 1:
                matrices['a_i'][:, q] = matrices['a_i'][:, q] * (-1)
    elif mode == "mts":
        for q, _ in enumerate(lambda_dir_1):
            if lambda_built[q] == 0:
                matrices['a_i'][:, q] = 0
                matrices['l_i'][q] = 0
            elif (lambda_built[q] == 1) & (lambda_dir_1[q, 0] == 0):
                matrices['a_i'][:, q] = matrices['a_i'][:, q] * (-1)
                lambda_dir_1[q, np.where(lambda_dir_1[q, 1:] == 0)[0]] = 1
                lambda_dir_2[q, np.where(lambda_dir_2[q, 1:] == 1)[0]] = 0

    lambda_dir_1 = lambda_dir_1[lambda_dir_1.any(axis=1)]
    lambda_dir_2 = lambda_dir_2[lambda_dir_2.any(axis=1)]

    if mode == "sts":
        p_lin = p_11 + p_21
    elif mode == "mts":
        p_lin = p_cap
    p_lin_opt = np.delete(p_lin, np.where(~matrices['a_i'].any(axis=0)))
    pos_opt = np.delete(matrices['position'], np.where(~matrices['a_i'].any(axis=1)), axis=0)
    a_c_opt = np.delete(matrices['a_c'], np.where(~matrices['a_i'].any(axis=1)), axis=0)
    a_p_opt = np.delete(matrices['a_p'], np.where(~matrices['a_i'].any(axis=1)), axis=0)
    a_i_opt = matrices['a_i']
    a_i_opt = np.delete(a_i_opt, np.where(~a_i_opt.any(axis=0)), axis=1)
    a_i_opt = np.delete(a_i_opt, np.where(~a_i_opt.any(axis=1)), axis=0)
    l_i_opt = matrices['l_i'][matrices['l_i'] != 0]

    # Prepare variables for the calculation of linear diameters
    a_i_shape_opt = np.shape(a_i_opt)   # (rows 0, columns 1)
    d_lin = np.zeros(a_i_shape_opt[1])  # Initialize linear diameters
    v_lin = np.zeros(a_i_shape_opt[1])  # Initialize velocities
    supply_temp_opt = np.ones(a_i_shape_opt[1]) * t_supply   # Assign constant supply temperature
    return_temp_opt = np.ones(a_i_shape_opt[1]) * t_return   # Assign constant return temperature

    def equations(v):
        vel, d = v
        reynolds = (settings.Water.density * vel * d) / settings.Water.dynamic_viscosity
        f = (-1.8 * np.log10((settings.Piping.roughness / (3.7 * d)) ** 1.11 + 6.9 / reynolds))**-2  # friction factor
        eq1 = vel - np.sqrt((2 * settings.Piping.max_pr_loss * d) / (f * settings.Water.density))  # eq. for diameter
        eq2 = mass_lin - settings.Water.density * vel * (np.pi / 4) * d ** 2    # eq. for velocity
        return [eq1, eq2]

    m_lin = (p_lin_opt*1000)/(settings.Water.heat_capacity_cp * (supply_temp_opt - return_temp_opt))

    for h in range(a_i_shape_opt[1]):
        mass_lin = m_lin[h]
        sol = root(equations, (0.5, 0.02), method='lm')
        if sol.success:
            v_lin[h], d_lin[h] = sol.x
        else:
            print(h, 'failed to calculate diameter and velocity!')

    res = dict(
        a_i=a_i_opt,
        a_p=a_p_opt,
        a_c=a_c_opt,
        q_c=matrices['q_c'],
        q_c_con=q_c_real,
        l_i=l_i_opt,
        d_i_0=d_lin,
        m_i_0=m_lin,
        lambda_dir_1=lambda_dir_1,
        lambda_dir_2=lambda_dir_2,
        position=pos_opt,
    )

    return res

