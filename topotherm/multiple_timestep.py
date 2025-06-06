"""This module contains the optimization models for the multiple-timestep
district heating network design.

The module contains the following functions:
    * annuity: Calculate the annuity factor
    * model: Create the optimization model for the multiple time steps operation
"""

import numpy as np
import pyomo.environ as pyo

from topotherm.settings import Economics


def annuity(c_i: float, n: float) -> float:
    """Calculate the annuity factor.

    Args:
        c_i (float): Interest rate
        n (float): Number of years

    Returns:
        float: annuity factor
    """
    a = ((1 + c_i) ** n * c_i) / ((1 + c_i) ** n - 1)
    return a


def model(matrices: dict,
          sets: dict,
          regression_inst: dict,
          regression_losses: dict,
          economics: Economics,
          optimization_mode: str) -> pyo.ConcreteModel:
    """Create the optimization model for the thermo-hydraulic coupled with
    multiple time step operation.

    Args:
        matrices (dict): Dictionary with the matrices of the district heating
            network with keys a_i, a_p, a_c, l_i, position, q_c
        sets (dict): Dictionary with the sets for the optimization, obtained
            from matrices
        regression_inst (dict): Dictionary with the regression coefficients
            for the thermal capacity
        regression_losses (dict): Dictionary with the regression coefficients
            for the heat losses
        economics (topotherm.settings.Economics): Object with the economic
            parameters
        optimization_mode (str): Optimization mode, either 'economic' for
            economic or 'forced' for forced operation

    Returns:
        pyomo.environ.ConcreteModel: multiple time step optimization model
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
    mdl.set_t = pyo.Set(initialize=range(matrices['q_c'].shape[1]),
                        doc='Time steps')
    mdl.dirs = pyo.Set(initialize=['ij', 'ji'],
                       doc='Set of pipe directions.')
    mdl.flow = pyo.Set(initialize=['in', 'out'],
                       doc='Flow direction in the pipe')

    # Define the combined set for pipes with consumers in both directions
    mdl.cons = pyo.Set(
        initialize=[('ij', edge) for edge in sets['connection_c_ij']] +
                   [('ji', edge) for edge in sets['connection_c_ji']],
        dimen=2,
        doc='Pipes with consumer in both directions')

    mdl.forced_edges = pyo.Set(
        initialize=np.concatenate([sets['connection_c_ij'], sets['connection_c_ji']]),
        doc='Pipes with a consumer in direction ij or ji'
    )

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

    # Define variables
    pipe_power = {'bounds': (0, p_max_pipe_const),
                  'domain': pyo.NonNegativeReals,
                  'initialize': p_max_pipe_const}
    mdl.P = pyo.Var(
        mdl.dirs, mdl.flow, mdl.set_n_i, mdl.set_t,
        doc='Heat power at the pipes',
        **pipe_power)

    mdl.P_cap = pyo.Var(
        mdl.set_n_i,
        doc='Thermal capacity of each pipe',
        **pipe_power
    )

    # Binaries for the chosen direction of a pipe
    mdl.lambda_ = pyo.Var(
        mdl.dirs, mdl.set_n_i, mdl.set_t,
        initialize=1,
        domain=pyo.Binary,
        doc='Binary direction decisions')

    # Binary building decision of a pipe
    mdl.lambda_b = pyo.Var(
        mdl.set_n_i,
        initialize=1,
        domain=pyo.Binary,
        doc='Binary building decision'
    )

    mdl.P_source = pyo.Var(
        mdl.set_n_p, mdl.set_t,
        doc='Thermal power of the source',
        domain=pyo.NonNegativeReals,
        bounds=lambda m, i, t: (0, economics.source_max_power[i]),
        initialize=lambda m, i, t: economics.source_max_power[i])

    mdl.P_source_inst = pyo.Var(
        mdl.set_n_p,
        doc='Thermal capacity of the heat source',
        domain=pyo.NonNegativeReals,
        bounds=lambda m, i: (economics.source_min_power[i], economics.source_max_power[i]),
        initialize=lambda m, i: economics.source_max_power[i])

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
        if optimization_mode == "forced":
            sink = sum(matrices['q_c'][k, t]
                       for k in sets['a_c_out'][j])
        elif (optimization_mode == "economic") | (optimization_mode == "sensitivity"):
            sink = sum(m.lambda_b[k] for k in sets['a_c_out_edge'][j]) \
                   * sum(matrices['q_c'][k, t] for k in sets['a_c_out'][j])
        return node_to_pipe + pipe_to_node + sources + sink == 0

    mdl.cons_nodal_balance = pyo.Constraint(
        mdl.set_n, mdl.set_t,
        rule=nodal_power_balance,
        doc='Nodal power balance for each time step')

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
        variable = (regression_losses['a']
                    * regression_inst['power_flow_max_partload']
                    * m.P[d, 'in', j, t])
        fix = regression_losses['b'] * m.lambda_[d, j, t]
        return flows - (variable + fix) * matrices['l_i'][j] == 0

    mdl.cons_power_balance_pipe = pyo.Constraint(
        mdl.dirs, mdl.set_n_i, mdl.set_t,
        rule=power_balance_pipe,
        doc='Power balance for each pipe and time step')

    def power_bigm_P(m, d, j, t):
        lhs = m.P[d, 'in', j, t] - p_max_pipe_const * m.lambda_[d, j, t]
        rhs = 0
        return lhs <= rhs
    mdl.cons_power_bigm_P = pyo.Constraint(
        mdl.dirs, mdl.set_n_i, mdl.set_t,
        rule=power_bigm_P, doc='Big-M constraint for power flow')

    def power_max_p_built_const(m, j):
        return m.P_cap[j] - p_max_pipe_const * m.lambda_b[j] <= 0
    mdl.cons_power_max_P_built_const = pyo.Constraint(
        mdl.set_n_i, rule=power_max_p_built_const, doc='Maximum Powerflow constant i->j / j->i')

    def power_max_p_built(m, d, j, t):
        return m.P[d, 'in', j, t] - m.P_cap[j] <= 0
    mdl.cons_power_max_P_built = pyo.Constraint(
        mdl.dirs, mdl.set_n_i, mdl.set_t, rule=power_max_p_built,
        doc='Maximum Powerflow according to capacity i->j / j->i')

    def built_usage_mapping(m, d, j, t):
        return m.lambda_[d, j, t] - m.lambda_b[j] <= 0
    mdl.cons_built_usage_mapping_help1 = pyo.Constraint(mdl.dirs, mdl.set_n_i, mdl.set_t,
                                                        rule=built_usage_mapping,
                                                        doc='Map lambda direction according to lambda_built')

    def one_pipe(m, j, t):
        return m.lambda_['ij', j, t] + m.lambda_['ji', j, t] <= 1

    mdl.one_pipe = pyo.Constraint(mdl.set_n_i, mdl.set_t,
                                    rule=one_pipe,
                                    doc='Just one Direction for each pipe')

    def connection_to_consumer_eco(m, d, j, t):
        return m.lambda_[d, j, t] <= sets[f'lambda_c_{d}'][j]

    def connection_to_consumer_fcd(m, d, j, t):
        return m.lambda_[d, j, t] == sets[f'lambda_c_{d}'][j]

    if (optimization_mode == "economic") | (optimization_mode == "sensitivity"):
        msg_ = """Constraint if houses have their own connection-pipe
            and set the direction (ij or ji)"""
        mdl.cons_connection_to_consumer = pyo.Constraint(
            mdl.cons, mdl.set_t,
            rule=connection_to_consumer_eco,
            doc=msg_)

    elif optimization_mode == "forced":
        msg_ = """Constraint if houses have their own connection-pipe
            and set the direction (ij or ji)"""
        mdl.cons_connection_to_consumer = pyo.Constraint(
            mdl.cons, mdl.set_t,
            rule=connection_to_consumer_fcd,
            doc=msg_)

        def connection_to_consumer_built_fcd(m, j):
            return m.lambda_b[j] >= 1
        mdl.cons_connection_forced = pyo.Constraint(
            mdl.forced_edges,
            rule=connection_to_consumer_built_fcd,
            doc='The house connection has to be built to satisfy the demand')

        def total_energy_conservation(m, t):
            return sum(m.P_source[k, t] for k in m.set_n_p) \
                - sum(m.P['ij', 'in', k, t] - m.P['ij', 'out', k, t] for k in m.set_n_i) \
                - sum(m.P['ji', 'in', k, t] - m.P['ji', 'out', k, t] for k in m.set_n_i) \
                - sum(matrices['q_c'][k, t] for k in m.set_n_c) == 0
        mdl.cons_total_energy_cons = pyo.Constraint(
            mdl.set_t,
            rule=total_energy_conservation,
            doc='Total energy conservation')

    if optimization_mode == "sensitivity":
        def total_energy_qc(m):
            return sum(sum(m.lambda_b[j] * matrices['flh_consumer'][k, t]
                       * matrices['q_c'][k, t] for k, j in m.consumer_edges)
                       for t in m.set_t) >= sets['q_c_tot']

        mdl.total_energy_qc = pyo.Constraint(rule=total_energy_qc)

        def consecutive_optimizations(m, j):
            return m.lambda_b[j] >= sets['lambda_b_previous'][j]
        mdl.consecutive_opt = pyo.Constraint(mdl.consumer_edges_only,
                                             rule=consecutive_optimizations)

    mdl.revenue = pyo.Var(doc='Revenue', domain=pyo.NegativeReals)
    mdl.revenue_constr = pyo.Constraint(
        expr=mdl.revenue == sum(
            sum(mdl.lambda_b[j] * matrices['flh_consumer'][k, t]
                * matrices['q_c'][k, t] for k, j in mdl.consumer_edges)
                for t in mdl.set_t) * economics.heat_price * (-1),
        doc='Revenue constraint')
    
    mdl.opex_source = pyo.Var(doc='OPEX Source', domain=pyo.NonNegativeReals)
    mdl.opex_source_constr = pyo.Constraint(
        expr=mdl.opex_source == sum(
            sum(mdl.P_source[k, t]
                * economics.source_price[k, t]
                * matrices['flh_source'][k, t]
                for k in mdl.set_n_p)
            for t in mdl.set_t),
        doc='OPEX Source constraint')
    
    mdl.capex_pipes = pyo.Var(doc='CAPEX Pipe', domain=pyo.NonNegativeReals)
    pipes = sum(
        ((mdl.P_cap[k] * regression_inst['a']
          + regression_inst['b'] * mdl.lambda_b[k])
          * annuity(economics.pipes_c_irr, economics.pipes_lifetime)
          * matrices['l_i'][k]) for k in mdl.set_n_i)

    mdl.capex_pipe_constr = pyo.Constraint(
        expr=mdl.capex_pipes == pipes,
        doc='CAPEX Pipe constraint')

    mdl.capex_source = pyo.Var(doc='CAPEX Source', domain=pyo.NonNegativeReals)
    mdl.capex_source_constr = pyo.Constraint(
        expr=mdl.capex_source == sum(mdl.P_source_inst[k]
                     * economics.source_c_inv[k]
                     * annuity(economics.source_c_irr[k],
                               economics.source_lifetime[k])
                     for k in mdl.set_n_p),
        doc='CAPEX Source constraint')

    mdl.obj = pyo.Objective(
        expr=mdl.capex_source + mdl.capex_pipes + mdl.opex_source + mdl.revenue,
        doc='Objective function')

    return mdl
