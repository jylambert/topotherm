"""This module contains the optimization models for the single-timestep 
district heating network design.

The module contains the following functions:
    * annuity: Calculate the annuity factor
    * model: Create the optimization model for the single time step operation
"""

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


def model(matrices: dict,
          sets: dict,
         regression_inst: dict,
         regression_losses: dict,
         economics: Economics,
         optimization_mode: str):
    """Create the optimization model for the thermo-hydraulic coupled with
    single time step operation. 

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
        mdl (pyomo.environ.ConcreteModel): pyomo model
    """

    # Check if the optimization mode is implemented
    if optimization_mode not in ['economic', 'forced']:
        raise NotImplementedError(
            "Optimization mode %s not implemented" % optimization_mode)

    # Initialize model
    mdl = pyo.ConcreteModel()

    # Big-M-Constraint for pipes
    p_max_pipe_const = float(regression_inst['power_flow_max_kW'].max())*2
    # Big-M-Constraint for source
    p_max_source = matrices['q_c'].sum()*2

    # Define index sets
    mdl.set_n_i = pyo.Set(initialize=range(sets['a_i_shape'][1]),
                            doc='N° of Pipe connections supply/return line')
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

    # Thermal power of the source
    source_power = {'bounds': (0, p_max_source),
                    'domain': pyo.PositiveReals,
                    'initialize': p_max_source}
    mdl.P_source = pyo.Var(
        mdl.set_n_p, mdl.set_t,
        doc='Thermal power of the source',
        **source_power)
    mdl.P_source_inst = pyo.Var(
        mdl.set_n_p,
        doc='Thermal capacity of the heat source',
        **source_power)

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
        elif optimization_mode == "economic":
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
        return node_to_pipe + pipe_to_node + sources + sink == 0

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
        variable = (regression_losses['a']
                * regression_inst['power_flow_max_partload']
                * m.P[d, 'in', j, t])
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

    if optimization_mode == "economic":
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
    else:
        raise NotImplementedError(
            "Optimization mode %s not implemented" % optimization_mode)

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

    def objective_function(m):
        fuel = sum(
            sum(m.P_source[k, t]
                * economics.source_price[k]
                * economics.source_flh[k]
                for k in m.set_n_p)
            for t in mdl.set_t)
        # CAREFUL HARDCODED FOR 0 TIME STEPS
        def pipes_fix(k):
            return ((m.P['ij', 'in', k, 0] + m.P['ji', 'in', k, 0])
                    * regression_inst['a'])
        def pipes_var(k):
            return (regression_inst['b']
                    * (m.lambda_['ij', k] + m.lambda_['ji', k]))
        pipes = (sum(((pipes_fix(k) + pipes_var(k))
                     * matrices['l_i'][k])
                     for k in m.set_n_i)
                     * annuity(economics.pipes_c_irr,
                               economics.pipes_lifetime))
        source = sum(m.P_source_inst[k]
                     * economics.source_c_inv[k]
                     * annuity(economics.source_c_irr[k],
                               economics.source_lifetime[k])
                     for k in m.set_n_p)

        # @TODO Implement consumer-specific flh in the economic mode
        if optimization_mode == "economic":
            revenue = (sum(
                sum(
                    sum(m.lambda_['ij', sets['a_i_in'][j].item()]
                        * matrices['q_c'][k, t]
                        for k in sets['a_c_out'][j]
                        if len(sets['a_i_in'][j]) > 0)
                    + sum(
                        (m.lambda_['ji', sets['a_i_out'][j].item()])
                        * matrices['q_c'][k, t]
                        for k in sets['a_c_out'][j]
                        if len(sets['a_i_out'][j]) > 0)
                        for j in mdl.set_n)
                for t in mdl.set_t)
                * economics.consumers_flh[0] * economics.heat_price * (-1))
        else:
            revenue = 0

        return fuel + pipes + source + revenue

    mdl.obj = pyo.Objective(rule=objective_function,
                              doc='Objective function')
    return mdl
