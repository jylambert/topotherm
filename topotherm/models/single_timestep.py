"""
This module contains the optimization models for the single-timestep
district heating network design.

Functions
---------
- ``annuity`` : Calculate the annuity factor.
- ``model`` : Create the optimization model for the single time-step operation.
"""

import numpy as np
import pyomo.environ as pyo

from topotherm.settings import Economics


def annuity(c_i, n):
    """
    Calculate the annuity factor.

    Parameters
    ----------
    c_i : float
        Interest rate.
    n : float
        Number of years.

    Returns
    -------
    float
        Annuity factor.
    """
    a = ((1 + c_i) ** n * c_i) / ((1 + c_i) ** n - 1)
    return a


def model(matrices: dict,
          sets: dict,
          regression_inst: dict,
          regression_losses: dict,
          economics: Economics,
          optimization_mode: str):
    """
    Create the optimization model for the thermo-hydraulic coupled with
    single time step operation.

    Parameters
    ----------
    matrices : dict
        Dictionary with the matrices of the district heating network with keys
        ``a_i``, ``a_p``, ``a_c``, ``l_i``, ``position``, ``q_c``.
    sets : dict
        Dictionary with the sets for the optimization, obtained from ``matrices``.
    regression_inst : dict
        Dictionary with the regression coefficients for the thermal capacity.
    regression_losses : dict
        Dictionary with the regression coefficients for the heat losses.
    economics : topotherm.settings.Economics
        Object with the economic parameters.
    optimization_mode : str
        Optimization mode, either ``'economic'`` for economic or ``'forced'``
        for forced operation.

    Returns
    -------
    pyomo.environ.ConcreteModel
        pyomo model.
    """
    # Check if the optimization mode is implemented
    if optimization_mode not in ['economic', 'forced', 'sensitivity']:
        raise NotImplementedError(
            "Optimization mode %s not implemented" % optimization_mode)

    # Initialize model
    mdl = pyo.ConcreteModel()

    # Big-M-Constraint for pipes
    p_max_pipe_const = float(regression_inst['power_flow_max_kW'].max())

    # Define index sets
    mdl.set_n_i = pyo.Set(initialize=range(sets['a_i_shape'][1]),
                            doc='Number of pipe connections supply/return line')
    mdl.set_n_p = pyo.Set(initialize=range(sets['a_p_shape'][1]),
                            doc='Number of producers')
    mdl.set_n_c = pyo.Set(initialize=range(sets['a_c_shape'][1]),
                            doc='Number of consumers')
    mdl.set_n = pyo.Set(initialize=range(sets['a_i_shape'][0]),
                          doc='Nodes in supply/return line')
    mdl.set_t = pyo.Set(initialize=[0],
                          doc='Time steps')
    mdl.dirs = pyo.Set(initialize=['ij', 'ji'],
                        doc='Set of pipe directions.')
    mdl.flow = pyo.Set(initialize=['in', 'out'],
                         doc='Flow direction in the pipe')
    mdl.set_con_ij = pyo.Set(initialize=sets['connection_c_ij'],
                               doc='Pipes with consumer in direction ij')
    mdl.set_con_ji = pyo.Set(initialize=sets['connection_c_ji'],
                               doc='Pipes with consumer in direction ji')

    mdl.consumer_edges = pyo.Set(
        initialize=[(i, np.where((matrices['a_i'][np.where(matrices['a_c'][:, i] == 1)[0].item(), :] == 1) |
                                 (matrices['a_i'][np.where(matrices['a_c'][:, i] == 1)[0].item(), :] == -1)
                                 )[0].item()) for i in range(sets['a_c_shape'][1])],
        dimen=2,
        doc='Assign to each consumer the corresponding pipe'
    )

    mdl.consumer_edges_only = pyo.Set(
        initialize=[np.where((matrices['a_i'][np.where(matrices['a_c'][:, i] == 1)[0].item(), :] == 1) |
                                 (matrices['a_i'][np.where(matrices['a_c'][:, i] == 1)[0].item(), :] == -1)
                                 )[0].item() for i in range(sets['a_c_shape'][1])]
    )

    # Define the combined set for pipes with consumers in both directions
    mdl.cons = pyo.Set(
        initialize=[('ij', edge) for edge in sets['connection_c_ij']] +
                   [('ji', edge) for edge in sets['connection_c_ji']],
        dimen=2,
        doc='Pipes with consumer in both directions')

    # Define variables
    pipe_power = {'bounds': (0, p_max_pipe_const),
                  'domain': pyo.NonNegativeReals,
                  'initialize': p_max_pipe_const}
    mdl.P = pyo.Var(
        mdl.dirs, mdl.flow, mdl.set_n_i, mdl.set_t,
        doc='Heat power at the pipes',
        **pipe_power)

    # Building decisions of a pipe
    mdl.lambda_ = pyo.Var(
        mdl.dirs, mdl.set_n_i,
        initialize=1,
        domain=pyo.Binary,
        doc='Binary direction decisions')

    # Definition of thermal power for each time step and each source
    mdl.P_source = pyo.Var(
        mdl.set_n_p, mdl.set_t,
        doc='Thermal power of the source',
        domain=pyo.NonNegativeReals,
        bounds=lambda m, i, t: (0, economics.source_max_power[i]),
        initialize=lambda m, i, t: economics.source_max_power[i])

    # Definition of thermal capacity of each source
    mdl.P_source_inst = pyo.Var(
        mdl.set_n_p,
        doc='Thermal capacity of the heat source',

        domain=pyo.NonNegativeReals,
        bounds=lambda m, i: (economics.source_min_power[i], economics.source_max_power[i]),
        initialize=lambda m, i: economics.source_max_power[i])

    # Definition of constraints
    def heat_source_inst(m, j, t):
        """Never exceed the installed capacity of the heat source."""
        return m.P_source[j, t] <= m.P_source_inst[j]

    mdl.cons_heat_source_inst = pyo.Constraint(
        mdl.set_n_p, mdl.set_t,
        rule=heat_source_inst,
        doc='Upper bound for the heat source supply delivery')

    def nodal_power_balance(m, j, t):
        """REFERENCE DIRECTION: left to right
                P_ji, in            P_ji, out
        PIPE    <-------    NODE    <-------    PIPE
                ------->            -------> 
                P_ij, out           P_ij, in
        
        Energy balance system: out - in = 0
        """
        pipe_to_node = sum(m.P['ji', 'in', k, t]
                           - m.P['ij', 'out', k, t]
                           for k in sets['a_i_in'][j])
        node_to_pipe = sum(m.P['ij', 'in', k, t]
                           - m.P['ji', 'out', k, t]
                           for k in sets['a_i_out'][j])
        sources = sum(- m.P_source[k, t]
                      for k in sets['a_p_in'][j])
        sink = 0
        if len(sets['a_i_in'][j]) == 0 and len(sets['a_i_out'][j]) == 0:
            raise ValueError("Node %s is not connected. Check preprocessing." % j)
        if optimization_mode == "forced":
            sink = sum(matrices['q_c'][k, t]
                       for k in sets['a_c_out'][j])
        elif (optimization_mode == "economic") | (optimization_mode == "sensitivity"):
            sink = (
                sum(
                    (m.lambda_['ij', sets['a_i_in'][j][0]])
                    * matrices['q_c'][k, t]
                    for k in sets['a_c_out'][j] if len(sets['a_i_in'][j]) > 0)
                + sum(
                    (m.lambda_['ji', sets['a_i_out'][j][0]])
                    * matrices['q_c'][k, t]
                    for k in sets['a_c_out'][j] if len(sets['a_i_out'][j]) > 0)
                )
        expr = (node_to_pipe + pipe_to_node + sources + sink == 0)
        return expr

    mdl.cons_nodal_balance = pyo.Constraint(
        mdl.set_n, mdl.set_t,
        rule=nodal_power_balance,
        doc='Nodal Power Balance')

    def power_balance_pipe(m, d, j, t):
        """Power balance for the pipes.
        
        P_ji, out            P_ji, in
        <-------    PIPE    <-------
        ------->            ------->
        P_ij, in            P_ij, out

        """
        # flows into and out of pipe
        flows = m.P[d, 'in', j, t] - m.P[d, 'out', j, t]
        # thermal losses calculation
        variable = regression_losses['a'] * m.P[d, 'in', j, t]
        fix = regression_losses['b'] * m.lambda_[d, j]
        return flows - (variable + fix) * matrices['l_i'][j] == 0

    mdl.cons_power_balance_pipe = pyo.Constraint(
        mdl.dirs, mdl.set_n_i, mdl.set_t,
        rule=power_balance_pipe,
        doc='Power balance pipe')

    def power_bigm_P(m, d, j, t):
        lhs = m.P[d, 'in', j, t] - p_max_pipe_const * m.lambda_[d, j]
        rhs = 0
        return lhs <= rhs

    mdl.cons_power_bigm_P = pyo.Constraint(
        mdl.dirs, mdl.set_n_i, mdl.set_t,
        rule=power_bigm_P, doc='Big-M constraint for power flow')

    def connection_to_consumer_eco(m, d, j):
        return m.lambda_[d, j] <= sets[f'lambda_c_{d}'][j]

    def connection_to_consumer_fcd(m, d, j):
        return m.lambda_[d, j] == sets[f'lambda_c_{d}'][j]

    if (optimization_mode == "economic") | (optimization_mode == "sensitivity"):
        msg_ = """Constraint if houses have their own connection-pipe
            and set the direction (ij)"""
        mdl.cons_connection_to_consumer = pyo.Constraint(
            mdl.cons,
            rule=connection_to_consumer_eco,
            doc=msg_)
    elif optimization_mode == "forced":
        msg_ = """Constraint if houses have their own connection-pipe
            and set the direction (ij)"""
        mdl.cons_connection_to_consumer = pyo.Constraint(
            mdl.cons,
            rule=connection_to_consumer_fcd,
            doc=msg_)

    def one_pipe(m, j):
        return m.lambda_['ij', j] + m.lambda_['ji', j] <= 1

    mdl.one_pipe = pyo.Constraint(mdl.set_n_i,
                                    rule=one_pipe,
                                    doc='Just one Direction for each pipe')

    # @TODO: Develop total energy conservation equation for the eco mode
    # (testing needed if beneficial)
    if optimization_mode == "forced":
        def total_energy_cons(m, t):
            sources = sum(m.P_source[k, t]
                          for k in m.set_n_p)
            pipes_ij = sum(m.P['ij', 'in', k, t] - m.P['ij', 'out', k, t]
                           for k in m.set_n_i)
            pipes_ji = sum(m.P['ji', 'in', k, t] - m.P['ji', 'out', k, t]
                           for k in m.set_n_i)
            demand = sum(matrices['q_c'][k, t]
                         for k in m.set_n_c)
            return sources - pipes_ij - pipes_ji - demand == 0

        mdl.cons_total_energy_cons = pyo.Constraint(
            mdl.set_t, rule=total_energy_cons,
            doc='Total energy conservation')

    if optimization_mode == "sensitivity":
        def total_energy_qc(m):
            return sum(
                sum(m.lambda_['ij', j] * matrices['q_c'][
                    np.where(matrices['a_c'][np.where(matrices['a_i'][:, j] == -1)[0], :][0] == 1)[0], t]
                    * matrices['flh_consumer'][
                        np.where(matrices['a_c'][np.where(matrices['a_i'][:, j] == -1)[0], :][0] == 1)[0], t] for j in
                    mdl.set_con_ij)
                + sum(m.lambda_['ji', j] * matrices['q_c'][
                    np.where(matrices['a_c'][np.where(matrices['a_i'][:, j] == 1)[0], :][0] == 1)[0], t]
                      * matrices['flh_consumer'][
                          np.where(matrices['a_c'][np.where(matrices['a_i'][:, j] == 1)[0], :][0] == 1)[0], t] for j in
                      mdl.set_con_ji)
                for t in mdl.set_t) >= sets['q_c_tot']

        mdl.total_energy_qc = pyo.Constraint(rule=total_energy_qc)

        def consecutive_optimizations(m, j):
            return sets['lambda_b_previous'][j] <= m.lambda_['ij', j] + m.lambda_['ji', j]
        mdl.consecutive_opt = pyo.Constraint(mdl.consumer_edges_only,
                                             rule=consecutive_optimizations)

    mdl.revenue = pyo.Var(doc='Revenue', domain=pyo.NegativeReals)
    mdl.revenue_constr = pyo.Constraint(
        expr=mdl.revenue == (
            sum(
                sum(
                    sum(mdl.lambda_['ij', sets['a_i_in'][j].item()]
                        * matrices['flh_consumer'][k, t]
                        * matrices['q_c'][k, t]
                        for k in sets['a_c_out'][j]
                        if len(sets['a_i_in'][j]) > 0)
                    + sum(
                        (mdl.lambda_['ji', sets['a_i_out'][j].item()])
                        * matrices['flh_consumer'][k, t]
                        * matrices['q_c'][k, t]
                        for k in sets['a_c_out'][j]
                        if len(sets['a_i_out'][j]) > 0)
                    for j in mdl.set_n)
                for t in mdl.set_t) * economics.heat_price * (-1)),
        doc='Revenue constraint')
    
    mdl.opex_source = pyo.Var(doc='OPEX Source', domain=pyo.NonNegativeReals)
    mdl.opex_source_constr = pyo.Constraint(
        expr=mdl.opex_source == sum(
            sum(mdl.P_source[k, t]
                * economics.source_price[k][t]
                * matrices['flh_source'][k, t]
                for k in mdl.set_n_p)
            for t in mdl.set_t),
        doc='OPEX Source constraint')

    mdl.capex_pipes = pyo.Var(doc='CAPEX Pipe', domain=pyo.NonNegativeReals)

    # CAREFUL HARDCODED FOR 0 TIME STEPS
    def pipes_fix(k):
        return ((mdl.P['ij', 'in', k, 0] + mdl.P['ji', 'in', k, 0])
                * regression_inst['a'])

    def pipes_var(k):
        return (regression_inst['b']
                * (mdl.lambda_['ij', k] + mdl.lambda_['ji', k]))

    mdl.capex_pipe_constr = pyo.Constraint(
        expr=mdl.capex_pipes == sum(
            (pipes_fix(k) + pipes_var(k)) * matrices['l_i'][k]
            for k in mdl.set_n_i)
            * annuity(economics.pipes_c_irr, economics.pipes_lifetime),
        doc='CAPEX Pipe constraint')

    mdl.capex_source = pyo.Var(doc='CAPEX Source', domain=pyo.NonNegativeReals)
    mdl.capex_source_constr = pyo.Constraint(
        expr=mdl.capex_source == sum(
            mdl.P_source_inst[k] * economics.source_c_inv[k]
            * annuity(economics.source_c_irr[k], economics.source_lifetime[k])
            for k in mdl.set_n_p),
        doc='CAPEX Source constraint')

    mdl.obj = pyo.Objective(
        expr=mdl.capex_source + mdl.capex_pipes + mdl.opex_source + mdl.revenue,
        doc='Objective function')
    return mdl
