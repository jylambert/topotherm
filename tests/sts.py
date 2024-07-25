import os

import pyomo.environ as pyo
import pandas as pd

import topotherm as tt
from topotherm.settings import Optimization

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
    # Load the district
    current_path = os.path.dirname(os.path.abspath(__file__))
    mat = tt.fileio.load(os.path.join(current_path, 'data'))

    # regression
    r_thermal_cap, r_heat_loss = read_regression(
        os.path.join(current_path, 'data', 'regression.csv'), 0)

    model_sets = tt.model.create_sets(mat)
    model = tt.model.sts_orig(mat, model_sets, r_thermal_cap, r_heat_loss,
                         Optimization().economics, "forced")

    # Optimization initialization
    opt = pyo.SolverFactory(solver_name)
    opt.options['mipgap'] = 0.01
    opt.options['timelimit'] = 600

    result = opt.solve(model, tee=True)
    assert result.solver.status == pyo.SolverStatus.ok
    # assert that the objective value is less than 2% away from the expected
    # value
    assert (abs(pyo.value(model.obj)) - 4.6259e+06) < 0.02 * 4.6259e+06


def test_sts_eco(request):
    """Main function to run the optimization"""
    solver_name = request.config.getoption("--solver")
    # Load the district
    current_path = os.path.dirname(os.path.abspath(__file__))
    mat = tt.fileio.load(os.path.join(current_path, 'data'))

    # regression
    r_thermal_cap, r_heat_loss = read_regression(
        os.path.join(current_path, 'data', 'regression.csv'), 0)

    model_sets = tt.model.create_sets(mat)
    model = tt.model.sts(mat, model_sets, r_thermal_cap, r_heat_loss,
                         Optimization().economics, "eco")

    # Optimization initialization
    opt = pyo.SolverFactory(solver_name)
    opt.options['mipgap'] = 0.01
    opt.options['timelimit'] = 600

    result = opt.solve(model, tee=True)
    assert result.solver.status == pyo.SolverStatus.ok
    assert (abs(pyo.value(model.obj)) - 4.01854e+04) < 0.02 * 4.01854e+04
