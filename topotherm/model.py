"""This module contains the optimization models for the thermo-hydraulic coupled district
heating network.

The module contains the following functions:
    * annuity: Calculate the annuity factor
    * create_sets: Create sets for the optimization based on the incidence matrices
    * sts: Create the optimization model for the thermo-hydraulic coupled with single time
    step operation
    * mts_easy: Create the optimization model for the thermo-hydraulic coupled with multiple time 
    step operation. The model is based on the STS model and implements a simplified themal
    loss equation to alleviate the computational burden.
    * mts: Create the optimization model for the thermo-hydraulic coupled with multiple time
    step operation. The model is based on the MTS easy model and implements the full thermal
    loss equation. The model is more accurate but also more computationally expensive.
"""

from collections import defaultdict

import numpy as np
import pyomo.environ as pyo

from topotherm.settings import Economics


def annuity(c_i, n):
    """Calculate the annuity factor.

    Args:
        c_i (float): Interest rate
        n (float): Number of years

    Returns:
        float: annuity factor
    """
    a = ((1 + c_i) ** n * c_i) / ((1 + c_i) ** n - 1)
    return a


def create_sets(matrices):
    """Create sets for the optimization. The sets are used to define the variables and constraints.

    Args:
        matrices (dict): Dictionary with the matrices of the district heating network with keys 
        a_i, a_p, a_c

    Returns:
        s: dictionary with the sets for the optimization
    """
    s = {}
    # Shapes of matrices
    s['a_i_shape'] = np.shape(matrices['a_i'])  # (rows 0, columns 1)
    s['a_p_shape'] = np.shape(matrices['a_p'])
    s['a_c_shape'] = np.shape(matrices['a_c'])

    # where consumers are connected to the network
    consumers = np.where(matrices['a_c'].sum(axis=1) == 1)[0]

    s['connection_c_ij'] = np.where(matrices['a_i'][consumers, :].sum(axis=0) == -1)[0]
    s['lambda_c_ij'] = np.zeros(s['a_i_shape'][1])
    s['lambda_c_ij'][s['connection_c_ij']] = 1

    s['connection_c_ji'] = np.where(matrices['a_i'][consumers, :].sum(axis=0) == 1)[0]
    s['lambda_c_ji'] = np.zeros(s['a_i_shape'][1])
    s['lambda_c_ji'][s['connection_c_ji']] = 1

    # create a dict with {row: [outgoing pipes], row: [ingoing pipes]}
    s['a_i_out'] = {}
    s['a_i_in'] = {}
    for i in range(s['a_i_shape'][0]):
        s['a_i_out'][i] = np.where(matrices['a_i'][i, :] == 1)[0]
        s['a_i_in'][i] = np.where(matrices['a_i'][i, :] == -1)[0]

    s['a_p_in'] = {}
    for i in range(s['a_p_shape'][0]):
        s['a_p_in'][i] = np.where(matrices['a_p'][i, :] == -1)[0]

    s['a_c_out'] = {}
    for i in range(s['a_c_shape'][0]):
        s['a_c_out'][i] = np.where(matrices['a_c'][i, :] == 1)[0]
    return s


def sts(matrices, sets, regression_caps, regression_losses, opt_mode):
    """Create the optimization model for the thermo-hydraulic coupled with single time
    step operation. 

    Args:
        matrices (dict): Dictionary with the matrices of the district heating network with keys 
        a_i, a_p, a_c, l_i, position, q_c
        sets (dict): Dictionary with the sets for the optimization, obtained from matrices
        regression_caps (dict): Dictionary with the regression coefficients for the thermal capacity
        regression_losses (dict): Dictionary with the regression coefficients for the heat losses
    
    Returns:
        model (pyomo.environ.ConcreteModel): pyomo model
    """
    # @TODO: Look with Amedeo if q_c can be adapted to dimensionless vector, (in theory it is possible to do
    # @TODO: a unidirectional flow formulation with multiple time step with topotherm sts)
    model = pyo.ConcreteModel()

    p_max_pipe_const = float(regression_caps['power_flow_max_kW'][-1])  # Big-M-Constraint for pipes
    p_max_source = matrices['q_c'].sum()*2  # Big-M-Constraint for source

    # Define index sets
    model.set_n_i = pyo.Set(initialize=range(sets['a_i_shape'][1]),
                            doc='Number of defined pipe connections supply/return line')    
    model.set_n_p = pyo.Set(initialize=range(sets['a_p_shape'][1]),
                            doc='Number of producers')
    model.set_n_c = pyo.Set(initialize=range(sets['a_c_shape'][1]),
                            doc='Number of consumers')
    model.set_n = pyo.Set(initialize=range(sets['a_i_shape'][0]),
                          doc='Nodes in supply/return line')
    model.set_t = pyo.Set(initialize=[0],
                          doc='Time steps')
    model.set_con_ij = pyo.Set(initialize=sets['connection_c_ij'],
                               doc='Pipes with consumer in direction ij')
    model.set_con_ji = pyo.Set(initialize=sets['connection_c_ji'],
                               doc='Pipes with consumer in direction ji')

    # Define variables
    pipe_power = dict(
        bounds=(0, p_max_pipe_const), domain=pyo.NonNegativeReals, initialize=p_max_pipe_const)
    model.P_11 = pyo.Var(model.set_n_i, model.set_t,
                         doc='Heat power at the entry of the pipe', **pipe_power)
    model.P_12 = pyo.Var(model.set_n_i, model.set_t,
                         doc='Heat power at the exit of the pipe', **pipe_power)
    model.P_21 = pyo.Var(model.set_n_i, model.set_t,
                         doc='Heat power at the entry of the pipe', **pipe_power)
    model.P_22 = pyo.Var(model.set_n_i, model.set_t,
                         doc='Heat power at the exit of the pipe', **pipe_power)
    # Building decisions of a pipe
    model.lambda_dir_1 = pyo.Var(model.set_n_i, domain=pyo.Binary, initialize=1,
                                 doc='Direction decision ij')
    model.lambda_dir_2 = pyo.Var(model.set_n_i, domain=pyo.Binary, initialize=0,
                                 doc='Direction decision ji')
    # Thermal power of the source
    source_power = dict(bounds=(0, p_max_source), domain=pyo.PositiveReals, initialize=p_max_source)
    model.P_source = pyo.Var(model.set_n_p, model.set_t,
                             doc='Thermal power of the source', **source_power)
    model.P_source_cap = pyo.Var(model.set_n_p,
                                 doc='Thermal capacity of the heat source', **source_power)

    def heat_source_cap(m, j, t):
        return m.P_source[j, t] <= m.P_source_cap[j]
    model.cons_heat_source_cap = pyo.Constraint(model.set_n_p, model.set_t, rule=heat_source_cap,
                                                doc='Investment costs for the heat source')
    # @TODO: Check if nodal power balance is the same for forced and eco (it should be the case, but testing is needed)


    def nodal_power_balance(m, j, t):
        term1 = sum(m.P_11[k, t] - m.P_22[k, t] for k in sets['a_i_out'][j]) # Sum of outgoing flows from pipes
        term2 = sum(m.P_21[k, t] - m.P_12[k, t] for k in sets['a_i_in'][j])  # Sum of incoming flows from pipes
        term3 = sum(- m.P_source[k, t] for k in sets['a_p_in'][j])           # Flows from producer
        term4 = 0
        if opt_mode == "forced":
            term4 = sum(matrices['q_c'][k, t] for k in sets['a_c_out'][j])
        elif opt_mode == "eco":
            term4 = sum((m.lambda_dir_1[sets['a_i_in'][j][0]])
                         * matrices['q_c'][k, t] for k in sets['a_c_out'][j] if len(sets['a_i_in'][j]) > 0) \
                   + sum((m.lambda_dir_2[sets['a_i_out'][j][0]])
                         * matrices['q_c'][k, t] for k in sets['a_c_out'][j] if len(sets['a_i_out'][j]) > 0)
        return term1 + term2 + term3 + term4 == 0
    model.cons_nodal_balance = pyo.Constraint(model.set_n, model.set_t, rule=nodal_power_balance,
                                              doc='Nodal Power Balance')

    def power_balance_pipe_12(m, j, t):
        term1 = m.P_11[j, t] - m.P_12[j, t]
        reg1 = regression_losses['params']['a'] * regression_caps['power_flow_max_partload'] \
            * m.P_11[j, t]
        reg2 = regression_losses['params']['b'] * m.lambda_dir_1[j]
        return term1 - (reg1 + reg2) * matrices['l_i'][j] == 0
    model.cons_power_balance_pipe_12 = pyo.Constraint(model.set_n_i, model.set_t,
                                                      rule=power_balance_pipe_12,
                                                      doc='Power balance pipe i->j')

    def power_balance_pipe_21(m, j, t):
        reg1 = regression_losses['params']['a'] * regression_caps['power_flow_max_partload'] \
            * m.P_21[j, t]
        reg2 = regression_losses['params']['b'] * m.lambda_dir_2[j]
        return m.P_21[j, t] - m.P_22[j, t] - (reg1 + reg2) * matrices['l_i'][j] == 0
    model.cons_power_balance_pipe_21 = pyo.Constraint(model.set_n_i, model.set_t,
                                                      rule=power_balance_pipe_21,
                                                      doc='Power balance pipe j->i')

    def power_max_p_11_const(m, j, t):
        return m.P_11[j, t] - p_max_pipe_const * m.lambda_dir_1[j] <= 0
    model.cons_power_max_P_11_const = pyo.Constraint(model.set_n_i, model.set_t,
                                                     rule=power_max_p_11_const,
                                                     doc='Maximum Powerflow constant i->j')

    def power_max_p_21_const(m, j, t):
        return m.P_21[j, t] - p_max_pipe_const * m.lambda_dir_2[j] <= 0
    model.cons_power_max_P_21_const = pyo.Constraint(model.set_n_i, model.set_t,
                                                     rule=power_max_p_21_const,
                                                     doc='Maximum Powerflow constant j->i')
    if opt_mode == "eco":
        def connection_to_consumer_ij(m, j):
            return m.lambda_dir_1[j] <= sets['lambda_c_ij'][j]
        msg_ = 'Constraint if houses have their own connection-pipe and set the direction (ij)'
        model.cons_connection_to_consumer_ij = pyo.Constraint(model.set_con_ij,
                                                              rule=connection_to_consumer_ij,
                                                              doc=msg_)

        def connection_to_consumer_ji(m, j):
            return m.lambda_dir_2[j] <= sets['lambda_c_ji'][j]
        msg = 'Constraint if houses have their own connection-pipe and set the direction (ji)'
        model.cons_connection_to_consumer_ji = pyo.Constraint(model.set_con_ji,
                                                              rule=connection_to_consumer_ji,
                                                              doc=msg)
    if opt_mode == "forced":
        def connection_to_consumer_ij(m, j):
            return m.lambda_dir_1[j] == sets['lambda_c_ij'][j]
        msg_ = 'Constraint if houses have their own connection-pipe and set the direction (ij)'
        model.cons_connection_to_consumer_ij = pyo.Constraint(model.set_con_ij,
                                                              rule=connection_to_consumer_ij,
                                                              doc=msg_)

        def connection_to_consumer_ji(m, j):
            return m.lambda_dir_2[j] == sets['lambda_c_ji'][j]
        msg = 'Constraint if houses have their own connection-pipe and set the direction (ji)'
        model.cons_connection_to_consumer_ji = pyo.Constraint(model.set_con_ji,
                                                              rule=connection_to_consumer_ji,
                                                              doc=msg)

    def one_pipe(m, j):
        return m.lambda_dir_1[j] + m.lambda_dir_2[j] <= 1
    model.one_pipe = pyo.Constraint(model.set_n_i, rule=one_pipe,
                                    doc='Just one Direction for each pipe')
    # @TODO: Develop total energy conservation equation for the eco mode (testing needed if beneficial)
    if opt_mode == "forced":
        def total_energy_cons(m, t):
            return sum(m.P_source[k, t] for k in m.set_n_p) \
                - sum(m.P_11[k, t] - m.P_12[k, t] for k in m.set_n_i)\
                - sum(m.P_21[k, t] - m.P_22[k, t] for k in m.set_n_i)\
                - sum(matrices['q_c'][k, t] for k in m.set_n_c) == 0
        model.cons_total_energy_cons = pyo.Constraint(model.set_t, rule=total_energy_cons,
                                                      doc='Total energy conservation')

    def objective_function(m):
        term1 = sum(
            sum(m.P_source[k, t] * Economics.source_price * Economics.flh for k in m.set_n_p)
            for t in model.set_t
                )
        term2 = sum(
            (
                ((m.P_11[k, 0] + m.P_21[k, 0]) * regression_caps['params']['a']
                    + regression_caps['params']['b'] * (m.lambda_dir_1[k]+m.lambda_dir_2[k])
                ) * annuity(Economics.c_irr, Economics.life_time) * matrices['l_i'][k]
            ) for k in m.set_n_i
        )
        term3 = sum(m.P_source_cap[k] * Economics.c_inv_source[k]
                    * annuity(Economics.c_irr, Economics.life_time) for k in m.set_n_p)

        if opt_mode == "eco":
            term4 = sum(sum(
                        sum((m.lambda_dir_1[sets['a_i_in'][j]])
                            * matrices['q_c'][k, t] for k in sets['a_c_out'][j] if len(sets['a_i_in'][j]) > 0)
                        + sum((m.lambda_dir_2[sets['a_i_out'][j]])
                            * matrices['q_c'][k, t] for k in sets['a_c_out'][j] if len(sets['a_i_out'][j]) > 0)
                        for j in model.set_n) for t in model.set_t) * Economics.flh * Economics.heat_price * (-1)
        else:
            term4 = 0

        return term1 + term2 + term3 + term4

    model.obj = pyo.Objective(rule=objective_function,
                              doc='Objective function')
    return model


# @TODO: discuss with jerry simplification strategies, since both models share a lot of equations.
# @TODO: change the flh_scaling somehow
# @TODO: implement existing pipes and sources

def mts_easy(matrices, sets, regression_caps, regression_losses, opt_mode, flh_scaling):
    """Create the optimization model for the thermo-hydraulic coupled with multiple time 
    step operation. The model is based on the STS model and implements a simplified themal
    loss equation to alleviate the computational burden.
    
    Args:
        matrices (dict): Dictionary with the matrices of the district heating network with keys 
        a_i, a_p, a_c, l_i, position, q_c
        sets (dict): Dictionary with the sets for the optimization
        regression_caps (dict): Dictionary with the regression coefficients for the thermal capacity
        regression_losses (dict): Dictionary with the regression coefficients for the heat losses
        flh_scaling (float): Scaling factor for the full load hours
    
    Returns:
        model (pyomo.environ.ConcreteModel): pyomo model
    """
    # @TODO: check how to solve the problem with the scaling of the flh
    model = pyo.ConcreteModel()

    p_max_pipe_const = float(regression_caps['power_flow_max_kW'][-1])  # Big-M-Constraint for pipes
    p_max_source = matrices['q_c'].sum()*2  # Big-M-Constraint for source
    model.flh = Economics.flh / flh_scaling

    # Define index sets
    model.set_n_i = pyo.Set(initialize=range(sets['a_i_shape'][1]),
                            doc='Number of defined pipe connections supply/return line')
    model.set_n_p = pyo.Set(initialize=range(sets['a_p_shape'][1]),
                            doc='Number of producers')
    model.set_n_c = pyo.Set(initialize=range(sets['a_c_shape'][1]),
                            doc='Number of consumers')
    model.set_n = pyo.Set(initialize=range(sets['a_i_shape'][0]),
                          doc='Nodes in supply/return line')
    model.set_t = pyo.Set(initialize=range(matrices['q_c'].shape[1]),
                          doc='Time steps')
    model.set_con_ij = pyo.Set(initialize=sets['connection_c_ij'],
                               doc='Pipes with consumer in direction ij')
    model.set_con_ji = pyo.Set(initialize=sets['connection_c_ji'],
                               doc='Pipes with consumer in direction ji')

    # Define variables
    pipe_power = {'bounds':(0, p_max_pipe_const),
                  'domain':pyo.NonNegativeReals,
                  'initialize':p_max_pipe_const}
    model.P_11 = pyo.Var(model.set_n_i, model.set_t,
                         doc='Heat power at the entry of the pipe', **pipe_power)
    model.P_12 = pyo.Var(model.set_n_i, model.set_t,
                         doc='Heat power at the exit of the pipe', **pipe_power)
    model.P_cap = pyo.Var(model.set_n_i,
                          doc='Thermal capacity of the pipe', **pipe_power)
    model.P_21 = pyo.Var(model.set_n_i, model.set_t,
                         doc='Heat power at the entry of the pipe', **pipe_power)
    model.P_22 = pyo.Var(model.set_n_i, model.set_t,
                         doc='Heat power at the exit of the pipe', **pipe_power)

    # Building decisions of a pipe
    model.lambda_built = pyo.Var(model.set_n_i, domain=pyo.Binary, initialize=1,
                                 doc='Building decision of a pipe')
    model.lambda_dir_1 = pyo.Var(model.set_n_i, model.set_t, domain=pyo.Binary, initialize=1,
                                 doc='Direction decision ij')
    model.lambda_dir_2 = pyo.Var(model.set_n_i, model.set_t, domain=pyo.Binary, initialize=0,
                                 doc='Direction decision ji')

    # Thermal power of the source
    source_power = {'bounds':(0, p_max_source),
                    'domain':pyo.PositiveReals,
                    'initialize':p_max_source}
    model.P_source = pyo.Var(model.set_n_p, model.set_t,
                             doc='Thermal power of the source', **source_power)
    model.P_source_cap = pyo.Var(model.set_n_p,
                                 doc='Thermal capacity of the heat source', **source_power)


    def heat_source_cap(m, j, t):
        return m.P_source[j, t] <= m.P_source_cap[j]
    model.cons_heat_source_cap = pyo.Constraint(model.set_n_p, model.set_t, rule=heat_source_cap,
                                                doc='Investment costs for the heat source')


    def nodal_power_balance(m, j, t):
        term1 = sum(m.P_11[k, t] - m.P_22[k, t] for k in sets['a_i_out'][j]) # Sum of outgoing flows from pipes
        term2 = sum(m.P_21[k, t] - m.P_12[k, t] for k in sets['a_i_in'][j])  # Sum of incoming flows from pipes
        term3 = sum(- m.P_source[k, t] for k in sets['a_p_in'][j])           # Flows from producer
        term4 = 0
        if opt_mode == "forced":
            term4 = sum(matrices['q_c'][k, t] for k in sets['a_c_out'][j])
        elif opt_mode == "eco":
            term4 = sum((m.lambda_dir_1[sets['a_i_in'][j][0], t])
                         * matrices['q_c'][k, t] for k in sets['a_c_out'][j] if len(sets['a_i_in'][j]) > 0) \
                   + sum((m.lambda_dir_2[sets['a_i_out'][j][0], t])
                         * matrices['q_c'][k, t] for k in sets['a_c_out'][j] if len(sets['a_i_out'][j]) > 0)
        return term1 + term2 + term3 + term4 == 0
    model.cons_nodal_balance = pyo.Constraint(model.set_n, model.set_t, rule=nodal_power_balance,
                                              doc='Nodal Power Balance')

    def power_balance_pipe_12(m, j, t):
        term1 = m.P_11[j, t] - m.P_12[j, t]
        reg1 = regression_losses['params']['a'] * regression_caps['power_flow_max_partload'] \
            * m.P_11[j, t]
        reg2 = regression_losses['params']['b'] * m.lambda_dir_1[j, t]
        return term1 - (reg1 + reg2) * matrices['l_i'][j] == 0
    model.cons_power_balance_pipe_12 = pyo.Constraint(
        model.set_n_i, model.set_t, rule=power_balance_pipe_12, doc='Power balance pipe i->j')

    def power_balance_pipe_21(m, j, t):
        reg1 = regression_losses['params']['a'] * regression_caps['power_flow_max_partload'] \
            * m.P_21[j, t]
        reg2 = regression_losses['params']['b'] * m.lambda_dir_2[j, t]
        return m.P_21[j, t] - m.P_22[j, t] - (reg1 + reg2) * matrices['l_i'][j] == 0
    model.cons_power_balance_pipe_21 = pyo.Constraint(
        model.set_n_i, model.set_t, rule=power_balance_pipe_21, doc='Power balance pipe j->i')

    def power_max_p_built_const(m, j):
        return m.P_cap[j] - p_max_pipe_const * m.lambda_built[j] <= 0
    model.cons_power_max_P_built_const = pyo.Constraint(
        model.set_n_i, rule=power_max_p_built_const, doc='Maximum Powerflow constant i->j')

    def power_max_p_11_built(m, j, t):
        return m.P_11[j, t] - regression_caps['power_flow_max_partload'] * m.P_cap[j] <= 0
    model.cons_power_max_P_11_built = pyo.Constraint(
        model.set_n_i, model.set_t, rule=power_max_p_11_built,
        doc='Maximum Powerflow according to capacity i->j')


    def power_max_p_21_built(m, j, t):
        return m.P_21[j, t] - regression_caps['power_flow_max_partload'] * m.P_cap[j] <= 0
    model.cons_power_max_P_21_built = pyo.Constraint(
        model.set_n_i, model.set_t, rule=power_max_p_21_built,
        doc='Maximum Powerflow according to capacity j->i')

    if opt_mode == "forced":
        def connection_to_consumer_ij(m, j, t):
            return m.lambda_dir_1[j, t] == sets['lambda_c_ij'][j]
        msg_ = 'Constraint if houses have their own connection-pipe and set the direction (ij)'
        model.cons_connection_to_consumer_ij = pyo.Constraint(model.set_con_ij, model.set_t,
                                                              rule=connection_to_consumer_ij,
                                                              doc=msg_)
        def connection_to_consumer_ji(m, j, t):
            return m.lambda_dir_2[j, t] == sets['lambda_c_ji'][j]
        msg = 'Constraint if houses have their own connection-pipe and set the direction (ji)'
        model.cons_connection_to_consumer_ji = pyo.Constraint(model.set_con_ji, model.set_t,
                                                              rule=connection_to_consumer_ji,
                                                              doc=msg)

    if opt_mode == "eco":
        def connection_to_consumer_ij(m, j, t):
            return m.lambda_dir_1[j, t] <= sets['lambda_c_ij'][j]
        msg_ = 'Constraint if houses have their own connection-pipe and set the direction (ij)'
        model.cons_connection_to_consumer_ij = pyo.Constraint(model.set_con_ij, model.set_t,
                                                              rule=connection_to_consumer_ij,
                                                              doc=msg_)

        def connection_to_consumer_ji(m, j, t):
            return m.lambda_dir_2[j, t] <= sets['lambda_c_ji'][j]
        msg = 'Constraint if houses have their own connection-pipe and set the direction (ji)'
        model.cons_connection_to_consumer_ji = pyo.Constraint(model.set_con_ji, model.set_t,
                                                              rule=connection_to_consumer_ji,
                                                              doc=msg)
    if opt_mode == "forced":
        def built_forced_ij(m, j):
            return m.lambda_built[j] >= 1
        model.cons_built_forced_ij = pyo.Constraint(model.set_con_ij, rule=built_forced_ij,
                                                    doc='The house connection has to be built to satisfy the demand')

        def built_forced_ji(m, j):
            return m.lambda_built[j] >= 1
        model.cons_built_forced_ji = pyo.Constraint(model.set_con_ji, rule=built_forced_ji,
                                                    doc='The house connection has to be built to satisfy the demand ji')

    def one_pipe(m, j, t):
        return m.lambda_dir_1[j, t] + m.lambda_dir_2[j, t] <= 1
    model.one_pipe = pyo.Constraint(model.set_n_i, model.set_t, rule=one_pipe,
                                    doc='Just one Direction for each pipe')

    def power_max_p_11_const(m, j, t):
        return m.P_11[j, t] - p_max_pipe_const * m.lambda_dir_1[j, t] <= 0
    model.cons_power_max_P_11_const = pyo.Constraint(
        model.set_n_i, model.set_t, rule=power_max_p_11_const,
        doc='Maximum Powerflow constant i->j')

    def power_max_p_21_const(m, j, t):
        return m.P_21[j, t] - p_max_pipe_const * m.lambda_dir_2[j, t] <= 0
    model.cons_power_max_P_21_const = pyo.Constraint(model.set_n_i, model.set_t,
                                                     rule=power_max_p_21_const,
                                                     doc='Maximum Powerflow constant j->i')
    if opt_mode == "forced":
        def total_energy_cons(m, t):
            return sum(m.P_source[k, t] for k in m.set_n_p) \
                - sum(m.P_11[k, t] - m.P_12[k, t] for k in m.set_n_i)\
                - sum(m.P_21[k, t] - m.P_22[k, t] for k in m.set_n_i)\
                - sum(matrices['q_c'][k, t] for k in m.set_n_c) == 0
        model.cons_total_energy_cons = pyo.Constraint(model.set_t, rule=total_energy_cons,
                                                      doc='Total energy conservation')

    def built_usage_mapping_help1(m, j, t):
        return m.lambda_dir_1[j, t] - m.lambda_built[j] <= 0
    model.cons_built_usage_mapping_help1 = pyo.Constraint(model.set_n_i, model.set_t,
                                                          rule=built_usage_mapping_help1)

    def built_usage_mapping_help2(m, j, t):
        return m.lambda_dir_2[j, t] - m.lambda_built[j] <= 0
    model.cons_built_usage_mapping_help2 = pyo.Constraint(model.set_n_i, model.set_t,
                                                          rule=built_usage_mapping_help2)

    def objective_function(m):

        term1 = sum(
            sum(m.P_source[k, t] * Economics.source_price * model.flh for k in m.set_n_p)
            for t in model.set_t
                )
        term2 = sum(
            (
                (m.P_cap[k] * regression_caps['params']['a']
                    + regression_caps['params']['b'] * m.lambda_built[k]
                 ) * annuity(Economics.c_irr, Economics.life_time) * matrices['l_i'][k]
            ) for k in m.set_n_i
        )
        term3 = sum(m.P_source_cap[k] * Economics.c_inv_source[k]
                    * annuity(Economics.c_irr, Economics.life_time) for k in m.set_n_p)

        if opt_mode == "eco":
            term4 = sum(sum(
                        sum((m.lambda_dir_1[sets['a_i_in'][j][0], t])
                            * matrices['q_c'][k, t] for k in sets['a_c_out'][j] if len(sets['a_i_in'][j]) > 0)
                        + sum((m.lambda_dir_2[sets['a_i_out'][j][0], t])
                            * matrices['q_c'][k, t] for k in sets['a_c_out'][j] if len(sets['a_i_out'][j]) > 0)
                        for j in model.set_n) for t in model.set_t) * model.flh * Economics.heat_price * (-1)
        else:
            term4 = 0
        return term1 + term2 + term3 + term4

    model.obj = pyo.Objective(rule=objective_function,
                              doc='Objective function')

    return model


def mts(matrices, sets, regression_caps, regression_losses, opt_mode, flh_scaling):
    """Create the optimization model for the thermo-hydraulic coupled with multiple time
    step operation. The model is based on the STS model and implements a simplified thermal
    loss equation to alleviate the computational burden.

    Args:
        matrices (dict): Dictionary with the matrices of the district heating network with keys
        a_i, a_p, a_c, l_i, position, q_c
        sets (dict): Dictionary with the sets for the optimization
        regression_caps (dict): Dictionary with the regression coefficients for the thermal capacity
        regression_losses (dict): Dictionary with the regression coefficients for the heat losses
        flh_scaling (float): Scaling factor for the full load hours
    Returns:
        model (pyomo.environ.ConcreteModel): pyomo model
    """

    model = mts_easy(matrices, sets, regression_caps, regression_losses, opt_mode, flh_scaling)
    p_max_pipe_const = float(regression_caps['power_flow_max_kW'][-1])  # Big-M-Constraint for pipes

    model.P_loss_1 = pyo.Var(model.set_n_i, model.set_t,
                             bounds=(0, p_max_pipe_const),
                             domain=pyo.NonNegativeReals,
                             initialize=0,
                             doc='Heat power at the exit of the pipe 1')
    model.P_loss_2 = pyo.Var(model.set_n_i, model.set_t,
                             bounds=(0, p_max_pipe_const),
                             domain=pyo.NonNegativeReals,
                             initialize=0,
                             doc='Heat power at the exit of the pipe 2')


    def calculation_loss_1_1(m, j, t):
        return m.P_loss_1[j, t] - p_max_pipe_const * m.lambda_dir_1[j, t] <= 0

    model.cons_calculation_loss_1_1 = pyo.Constraint(
        model.set_n_i, model.set_t, rule=calculation_loss_1_1, doc='Complex loss calc 1 to 1')

    def calculation_loss_1_2(m, j, t):
        reg1 = regression_losses['params']['a'] * regression_caps['power_flow_max_partload'] \
               * m.P_cap[j]
        reg2 = regression_losses['params']['b'] * m.lambda_dir_1[j, t]
        return (reg1 + reg2) * matrices['l_i'][j] - m.P_loss_1[j, t] \
            - p_max_pipe_const * (1 - m.lambda_dir_1[j, t]) <= 0

    model.cons_calculation_loss_1_2 = pyo.Constraint(
        model.set_n_i, model.set_t, rule=calculation_loss_1_2, doc='Complex loss calc 1 to 2')

    def calculation_loss_2_1(m, j, t):
        return m.P_loss_2[j, t] - p_max_pipe_const * m.lambda_dir_2[j, t] <= 0

    model.cons_calculation_loss_2_1 = pyo.Constraint(
        model.set_n_i, model.set_t, rule=calculation_loss_2_1)

    def calculation_loss_2_2(m, j, t):
        return (regression_losses['params']['a'] * regression_caps['power_flow_max_partload']
                * m.P_cap[j]
                + regression_losses['params']['b'] * m.lambda_dir_2[j, t]) * matrices['l_i'][j] \
            - m.P_loss_2[j, t] - p_max_pipe_const * (1 - m.lambda_dir_2[j, t]) <= 0

    model.cons_calculation_loss_2_2 = pyo.Constraint(
        model.set_n_i, model.set_t, rule=calculation_loss_2_2, doc='Complex loss calc 2 to 2')

    # delete previous power_balance constraints
    model.del_component(model.cons_power_balance_pipe_12)
    model.del_component(model.cons_power_balance_pipe_21)

    # add new power_balance constraints with complex losses
    def power_balance_pipe_12(m, j, t):
        return m.P_11[j, t] - m.P_12[j, t] - m.P_loss_1[j, t] == 0
    model.cons_power_balance_pipe_12_cmplx = pyo.Constraint(
        model.set_n_i, model.set_t, rule=power_balance_pipe_12,
        doc='Complex Power balance pipe j->i')

    def power_balance_pipe_21(m, j, t):
        return m.P_21[j, t] - m.P_22[j, t] - m.P_loss_2[j, t] == 0
    model.cons_power_balance_pipe_21_cmplx = pyo.Constraint(
        model.set_n_i, model.set_t, rule=power_balance_pipe_21,
        doc='Complex Power balance pipe j->i')

    return model

