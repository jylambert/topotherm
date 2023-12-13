# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from collections import defaultdict
import pyomo.environ as pyo
import pickle

from topotherm.plotting import plot_district
import topotherm.utils as utils
import topotherm.import_files as imp
import topotherm.settings as settings
import topotherm.precalculation_hydraulic as precalc


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


def sts(name, time_steps, file_path, result_path, variant, mode):
    # -------------------------------- Initialization of the optimization --------------------------------
    v_init = 0.5        # Set an initial velocity for hydraulic calculations
    temp_ambient = np.zeros([settings.Piping.number_diameter, 1]) - 20  # Setting outdoor temperature to -20 °C
    supply_temp = precalc.determine_feed_line_temp(temp_ambient[0, :], 90, 80, -14, 6) * np.ones([settings.Piping.number_diameter, 1])  # Determine supply temperature to 90 °C
    return_temp = np.ones([settings.Piping.number_diameter, 1]) * 55    # Set return temperature to 55 °C

    # ------------ Perform nonlinear Calculation for thermal losses and thermal transport capacity ------------
    m_max, p_max, regression_capacity = precalc.calc_regression_thermal_capacity(v_init, settings.Piping.diameter, settings.Piping.roughness, settings.Piping.max_pr_loss, supply_temp, return_temp, temp_ambient)
    heat_loss_power_flow, regression_heat_loss = precalc.calc_regression_heat_losses(m_max, supply_temp, np.ones(settings.Piping.number_diameter), 100*np.ones(settings.Piping.number_diameter), temp_ambient, p_max)

    regression_heat_loss_a = np.round(regression_heat_loss[0][1], 10)        # Round regression parameter for heat losses
    regression_heat_loss_b = np.round(regression_heat_loss[0][0], 6)       # Round regression parameter for heat losses

    regression_capacity_a = np.round(regression_capacity[0][1], 6)          # Round regression parameter for capacity
    regression_capacity_b = np.round(regression_capacity[0][0], 3)          # Round regression parameter for capacity

    p_max_kw = np.round(p_max/1000, 3)                                         # Determine maximal power flow  in kw
    p_max_part_load = p_max_kw[0, :] / p_max_kw[0, :].max()       # Part load according to outdoor temperature and feed line temperature

    # -------------------------------- Prepare Optimization  --------------------------------
    a_i, a_p, a_c, q_c, l_i, pos = imp.load_district(name, file_path)      # Load the district
    #a_i, a_p, a_c, q_c, l_i, pos = imp.load_district_benchmark(name, file_path)

    dir_res = result_path + f"/{name}/{variant}/"       # Result Path
    utils.create_dir(dir_res)                           # Create result path
    plot_district(a_i, a_p, a_c, dir_res + name + "_plot_district_initial", l_i, 0, pos, 0) # Save initial District

    flh = 2500
    #heat_price = 110 * 10**-3       # Selling Price for heat in €/kW
    source_price = 80 * 10**-3      # Price for heat in €/kW
    c_inv_source = np.array([0])      # Investment costs
    life_time = 40                       # Number of years for deprecation
    c_irr = 0.08                    # Interest rate
    p_max_pipe_const = float(p_max_kw[-1])  # Big-M-Constraint for pipes
    p_max_source = q_c.sum()*2              # Big-M-Constraint for source

    # Shapes of matrices
    a_i_shape = np.shape(a_i)  # (rows 0, columns 1)
    a_p_shape = np.shape(a_p)
    a_c_shape = np.shape(a_c)

    connection_c_ij = np.where(a_i[np.where(a_c.sum(axis=1) == 1)[0], :].sum(axis=0) == -1)[0]
    lambda_c_ij = np.zeros(a_i_shape[1])
    lambda_c_ij[connection_c_ij] = 1

    connection_c_ji = np.where(a_i[np.where(a_c.sum(axis=1) == 1)[0], :].sum(axis=0) == 1)[0]
    lambda_c_ji = np.zeros(a_i_shape[1])
    lambda_c_ji[connection_c_ji] = 1

    set_a_i_out = defaultdict(list)
    set_a_i_in = defaultdict(list)
    row_len_a_i = range(a_i_shape[0])
    for i in row_len_a_i:
        num_pipe_out = np.where(a_i[i, :] == 1)[0]
        set_a_i_out[i].append(num_pipe_out)
        num_pipe_in = np.where(a_i[i, :] == -1)[0]
        set_a_i_in[i].append(num_pipe_in)

    set_a_p_in = defaultdict(list)
    row_len_a_p = range(a_p_shape[0])
    for i in row_len_a_p:
        num_prod_out = np.where(a_p[i, :] == -1)[0]
        set_a_p_in[i].append(num_prod_out)

    set_a_c_out = defaultdict(list)
    row_len_a_c = range(a_c_shape[0])
    for i in row_len_a_c:
        num_cons_out = np.where(a_c[i, :] == 1)[0]
        set_a_c_out[i].append(num_cons_out)

    # -------------------------------- Model Creation --------------------------------
    model = pyo.ConcreteModel()

    # Define index sets
    model.set_n_i = pyo.Set(initialize=range(a_i_shape[1]))                # Number of defined pipe connections supply/return line
    model.set_n_p = pyo.Set(initialize=range(a_p_shape[1]))                # Number of producer
    model.set_n_c = pyo.Set(initialize=range(a_c_shape[1]))                # Number of consumer
    model.set_n = pyo.Set(initialize=range(a_i_shape[0]))                  # Number of nodes in supply/return line
    model.set_t = pyo.Set(initialize=range(time_steps))                    # Number of time steps
    model.set_con_ij = pyo.Set(initialize=connection_c_ij)                 # Number of pipes with consumer in direction ij
    model.set_con_ji = pyo.Set(initialize=connection_c_ji)                 # Number of pipes with consumer in direction ji

    # Define variables
    model.P_11 = pyo.Var(model.set_n_i, model.set_t, bounds=(0, p_max_pipe_const), domain=pyo.NonNegativeReals, initialize=p_max_pipe_const)         # Heat power at the entry of the pipe
    model.P_12 = pyo.Var(model.set_n_i, model.set_t, bounds=(0, p_max_pipe_const), domain=pyo.NonNegativeReals, initialize=p_max_pipe_const)         # Heat power at the exit of the pipe
    model.P_21 = pyo.Var(model.set_n_i, model.set_t, bounds=(0, p_max_pipe_const), domain=pyo.NonNegativeReals, initialize=p_max_pipe_const)         # Heat power at the entry of the pipe
    model.P_22 = pyo.Var(model.set_n_i, model.set_t, bounds=(0, p_max_pipe_const), domain=pyo.NonNegativeReals, initialize=p_max_pipe_const)         # Heat power at the exit of the pipe                                               # Building desicion of a pipe
    model.lambda_dir_1 = pyo.Var(model.set_n_i, domain=pyo.Binary, initialize=1)        # Direction decision ij
    model.lambda_dir_2 = pyo.Var(model.set_n_i, domain=pyo.Binary, initialize=0)       # Direction decision ji
    model.P_source = pyo.Var(model.set_n_p, model.set_t, bounds=(0, p_max_source), domain=pyo.PositiveReals, initialize=p_max_source)            # Thermal power of the source
    model.P_source_cap = pyo.Var(model.set_n_p, bounds=(0, p_max_source), domain=pyo.PositiveReals, initialize=p_max_source)                     # Thermal capacity of the heat source

    # Investment costs for the heat source
    def heat_source_cap(m, j, t):
        return m.P_source[j, t] <= m.P_source_cap[j]
    model.cons_heat_source_cap = pyo.Constraint(model.set_n_p, model.set_t, rule=heat_source_cap)

    # Nodal Power Balance
    def nodal_power_balance(m, j, t):
        return sum(m.P_11[k, t] - m.P_22[k, t] for k in set_a_i_out[j][0]) \
               + sum(m.P_21[k, t] - m.P_12[k, t] for k in set_a_i_in[j][0]) \
               + sum(- m.P_source[k, t] for k in set_a_p_in[j][0]) \
               + sum(q_c[k, t] for k in set_a_c_out[j][0]) == 0
    model.cons_nodal_balance = pyo.Constraint(model.set_n, model.set_t, rule=nodal_power_balance)

    # Power balance pipe i->j
    def power_balance_pipe_12(m, j, t):
        return m.P_11[j, t] - m.P_12[j, t] - (regression_heat_loss_a * p_max_part_load[0] * m.P_11[j, t] + regression_heat_loss_b * m.lambda_dir_1[j]) * l_i[j] == 0
    model.cons_power_balance_pipe_12 = pyo.Constraint(model.set_n_i, model.set_t, rule=power_balance_pipe_12)

    # Power balance pipe j->i
    def power_balance_pipe_21(m, j, t):
        return m.P_21[j, t] - m.P_22[j, t] - (regression_heat_loss_a * p_max_part_load[0] * m.P_21[j, t] + regression_heat_loss_b * m.lambda_dir_2[j]) * l_i[j] == 0
    model.cons_power_balance_pipe_21 = pyo.Constraint(model.set_n_i, model.set_t, rule=power_balance_pipe_21)

    # Maximum Powerflow constant i->j
    def power_max_p_11_const(m, j, t):
        return m.P_11[j, t] - p_max_pipe_const * m.lambda_dir_1[j] <= 0
    model.cons_power_max_P_11_const = pyo.Constraint(model.set_n_i, model.set_t, rule=power_max_p_11_const)

    # Maximum Powerflow constant j->i
    def power_max_p_21_const(m, j, t):
        return m.P_21[j, t] - p_max_pipe_const * m.lambda_dir_2[j] <= 0
    model.cons_power_max_P_21_const = pyo.Constraint(model.set_n_i, model.set_t, rule=power_max_p_21_const)

    # Constraint if houses have their own connection-pipe and set the direction (ij)
    def connection_to_consumer_ij(m, j):
        return m.lambda_dir_1[j] == lambda_c_ij[j]
    model.cons_connection_to_consumer_ij = pyo.Constraint(model.set_con_ij, rule=connection_to_consumer_ij)

    # Constraint if houses have their own connection-pipe and set the direction (ji)
    def connection_to_consumer_ji(m, j):
        return m.lambda_dir_2[j] == lambda_c_ji[j]
    model.cons_connection_to_consumer_ji = pyo.Constraint(model.set_con_ji, rule=connection_to_consumer_ji)

    # Just one Direction for each pipe
    def one_pipe(m, j):
        return m.lambda_dir_1[j] + m.lambda_dir_2[j] <= 1
    model.one_pipe = pyo.Constraint(model.set_n_i, rule=one_pipe)

    # -------------------------------- Create some variants --------------------------------
    # Total energy conservation
    def total_energy_cons(m, t):
        return sum(m.P_source[k, t] for k in m.set_n_p) \
            - sum(m.P_11[k, t] - m.P_12[k, t] for k in m.set_n_i)\
            - sum(m.P_21[k, t] - m.P_22[k, t] for k in m.set_n_i)\
            - sum(q_c[k, t] for k in m.set_n_c) == 0
    model.cons_total_energy_cons = pyo.Constraint(model.set_t, rule=total_energy_cons)

    def objective_function(m):
        return sum(sum(m.P_source[k, t] * source_price * flh for k in m.set_n_p) for t in model.set_t) +\
            sum((((m.P_11[k, 0] + m.P_21[k, 0]) * regression_capacity_a + regression_capacity_b * (m.lambda_dir_1[k]+m.lambda_dir_2[k])) * annuity(c_irr, life_time) * l_i[k]) for k in m.set_n_i) + \
            sum(m.P_source_cap[k] * c_inv_source[k] * annuity(c_irr, life_time) for k in m.set_n_p)
    model.obj = pyo.Objective(rule=objective_function)

    # -------------------------------- Initialize Optimization --------------------------------
    # Optimization initialization
    opt = pyo.SolverFactory('gurobi')
    opt.options['mipgap'] = settings.OptSettings.mip_gap
    opt.options['timelimit'] = settings.OptSettings.time_limit
    opt.options['logfile'] = dir_res + "optimization.log"
    #opt.options['Seed'] = 56324978

    # -------------------------------- Solve the Model --------------------------------
    result = opt.solve(model, tee=True)

    # -------------------------------- Process the results --------------------------------
    # Save model results to csv
    dfres = utils.model_to_df(model)
    dfres.to_csv(dir_res + "results.csv", sep=";")

    # save solver results
    dfsol = utils.solver_to_df(result, model, solver="gurobi")
    dfsol.to_csv(dir_res + "solver.csv", sep=";")

    # Create variables for the thermo-hydraulic coupled optimization
    data_dict = {}

    p_21 = np.zeros([a_i_shape[1], time_steps])
    p_11 = np.zeros([a_i_shape[1], time_steps])
    lambda_dir_1 = np.zeros([a_i_shape[1], time_steps])
    lambda_dir_2 = np.zeros([a_i_shape[1], time_steps])
    p_source_built = np.zeros(a_p_shape[1])

    for v in model.component_objects(pyo.Var, active=True):
        var_dict = {(v.name, index): pyo.value(v[index]) for index in v}
        data_dict.update(var_dict)
        if v.name == "P_21":
            for index in v:
                p_21[index] = pyo.value(v[index])
        if v.name == "P_11":
            for index in v:
                p_11[index] = pyo.value(v[index])
        if v.name == "lambda_dir_1":
            for index in v:
                lambda_dir_1[index] = pyo.value(v[index])
        if v.name == "lambda_dir_2":
            for index in v:
                lambda_dir_2[index] = pyo.value(v[index])
        if v.name == "P_source_cap":
            for index in v:
                p_source_built[index] = pyo.value(v[index])

    lambda_dir_1 = np.around(lambda_dir_1, 0)
    lambda_dir_2 = np.around(lambda_dir_2, 0)

    # Restart, Adaption of Incidence Matrix for the thermo-hydraulic coupled optimization
    for q in range(len(lambda_dir_1)):
        if lambda_dir_1[q] == 0 and lambda_dir_2[q] == 0:
            a_i[:, q] = 0
            l_i[q] = 0
        elif lambda_dir_2[q] == 1:
            a_i[:, q] = a_i[:, q] * (-1)

    p_lin = p_11 + p_21
    p_lin_opt = np.delete(p_lin, np.where(~a_i.any(axis=0)))
    pos_opt = np.delete(pos, np.where(~a_i.any(axis=1)), axis=0)
    a_c_opt = np.delete(a_c, np.where(~a_i.any(axis=1)), axis=0)
    a_p_opt = np.delete(a_p, np.where(~a_i.any(axis=1)), axis=0)
    #a_c_opt = np.delete(a_c, np.where(~a_c.any(axis=0)), axis=1)
    a_i_opt = a_i
    a_i_opt = np.delete(a_i_opt, np.where(~a_i_opt.any(axis=0)), axis=1)
    a_i_opt = np.delete(a_i_opt, np.where(~a_i_opt.any(axis=1)), axis=0)
    l_i_opt = l_i[l_i != 0]

    a_i_shape_opt = np.shape(a_i_opt)  # (rows 0, columns 1)
    #a_p_shape_opt = np.shape(a_p_opt)
    #a_c_shape_opt = np.shape(a_c_opt)
    d_lin2 = np.zeros(a_i_shape_opt[1])
    v_lin2 = np.zeros(a_i_shape_opt[1])
    supply_temp_opt = np.ones(a_i_shape_opt[1]) * supply_temp[0, 0]
    return_temp_opt = np.ones(a_i_shape_opt[1]) * return_temp[0, 0]


    def equations(vars):
        vel, d = vars
        reynolds = (settings.Water.density * vel * d) / settings.Water.dynamic_viscosity
        f = (-1.8 * np.log10((settings.Piping.roughness / (3.7 * d)) ** 1.11 + 6.9 / reynolds)) ** -2
        eq1 = vel - np.sqrt((2 * settings.Piping.max_pr_loss * d) / (f * settings.Water.density))
        eq2 = mass_lin - settings.Water.density * vel * (np.pi / 4) * d ** 2
        return [eq1, eq2]

    m_lin = (p_lin_opt*1000)/(settings.Water.heat_capacity_cp * (supply_temp_opt - return_temp_opt))

    for h in range(a_i_shape_opt[1]):
        mass_lin = m_lin[h]
        v_lin2[h], d_lin2[h] = fsolve(equations, (0.5, 0.02))

    d_i_0 = d_lin2
    m_i_0 = m_lin

    # Saving the objects:

    with open(dir_res + name + '_matrix_preprocessed.pkl', 'wb') as h:  # Python 3: open(..., 'wb')
        pickle.dump([a_i_opt, a_p_opt, a_c_opt, q_c, l_i_opt, d_i_0, m_i_0, pos_opt], h)

    # Save figure optimized districts
    plot_district(a_i_opt, a_p_opt, a_c_opt, dir_res + name + "_plot_district_optimized", l_i_opt, d_i_0, pos_opt, 1)

    return
