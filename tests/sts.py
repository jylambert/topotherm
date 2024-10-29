import os

import pyomo.environ as pyo
import pandas as pd
from pytest import approx

import topotherm as tt

def read_regression(path, i):
    """Read the regression coefficients for the thermal capacity and heat
    losses from csv file.
    """
    # read df and force floats
    df = pd.read_csv(path, sep=',', index_col=0, header=0, dtype=float)
    r_thermal_cap = {
        "power_flow_max_kW" : df.loc[i, "power_flow_max_kW"],
        "a": df.loc[i, "capacity_a"],
        "b": df.loc[i, "capacity_b"],
        "power_flow_max_partload": df.loc[i, "power_flow_max_partload"]
    }
    r_heat_loss = {
        "a": df.loc[i, "heat_loss_a"],
        "b": df.loc[i, "heat_loss_b"]
    }
    return r_thermal_cap, r_heat_loss


def test_sts_forced(request):
    """Main function to run the optimization"""
    solver_name = request.config.getoption("--solver")
    assert solver_name in ["scip", "gurobi", "cbc"], f"Unsupported solver: {solver_name}"
    # Load the district
    current_path = os.path.dirname(os.path.abspath(__file__))
    mat = tt.fileio.load(os.path.join(current_path, 'data_sts'))

    # regression
    r_thermal_cap, r_heat_loss = read_regression(
        os.path.join(current_path, 'data_sts', 'regression.csv'), 0)

    # import settings
    settings = tt.settings.load(os.path.join(current_path, 'data_sts', 'config.yaml'))

    # modify either in code or in the config file
    settings.economics.source_c_inv = [0.]  # no investment costs for sources
    settings.temperatures.supply = 90

    model_sets = tt.sets.create(mat)
    model = tt.single_timestep.model(
        matrices=mat,
        sets=model_sets,
        regression_inst=r_thermal_cap,
        regression_losses=r_heat_loss,
        economics=settings.economics,
        optimization_mode="forced")

    # Optimization initialization
    opt = pyo.SolverFactory(solver_name)
    opt.options['mipgap'] = 0.01
    opt.options['timelimit'] = 600

    result = opt.solve(model, tee=True)
    assert result.solver.status == pyo.SolverStatus.ok
    # assert that the objective value is less than 2% away from the expected
    # value
    assert pyo.value(model.obj) == approx(519676.4358105995, rel=0.02)


def test_sts_eco(request):
    """Main function to run the optimization"""
    solver_name = request.config.getoption("--solver")
    assert solver_name in ["scip", "gurobi", "cbc"], f"Unsupported solver: {solver_name}"

    # Load the district
    current_path = os.path.dirname(os.path.abspath(__file__))
    mat = tt.fileio.load(os.path.join(current_path, 'data_sts'))

    # regression
    r_thermal_cap, r_heat_loss = read_regression(
        os.path.join(current_path, 'data_sts', 'regression.csv'), 0)

    # import settings
    settings = tt.settings.load(os.path.join(current_path, 'data_sts', 'config.yaml'))

    # modify either in code or in the config file
    settings.economics.source_c_inv = [0.]  # no investment costs for sources
    settings.temperatures.supply = 90

    model_sets = tt.sets.create(mat)
    model = tt.single_timestep.model(
        matrices=mat,
        sets=model_sets,
        regression_inst=r_thermal_cap,
        regression_losses=r_heat_loss,
        economics=settings.economics,
        optimization_mode="economic")

    # Optimization initialization
    opt = pyo.SolverFactory(solver_name)
    opt.options['mipgap'] = 0.01
    opt.options['timelimit'] = 600

    result = opt.solve(model, tee=True)
    assert result.solver.status == pyo.SolverStatus.ok
    assert pyo.value(model.obj) == approx(-4.01854e+04, rel=0.02)
    # assert (abs(pyo.value(model.obj)) - 4.01854e+04) < 0.02 * 4.01854e+04
