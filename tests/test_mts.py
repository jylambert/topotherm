"""Test the multiple timestep optimization modeling workflow."""

from pathlib import Path

import numpy as np
import pyomo.environ as pyo
from pytest import approx

import topotherm as tt

SOLVERS = ["scip", "gurobi", "cbc", "highs", "cplex", "glpk"]
R_COSTS = {
    "power_flow_max_kW": np.array([6.9e04]),
    "a": 0.018377,
    "b": 567.335,
    "power_flow_max_partload": 1,
}
R_LOSS = {"a": 4.348e-07, "b": 0.02189}


def test_mts_forced(request):
    """Main function to run the optimization"""
    solver_name = request.config.getoption("--solver")
    assert solver_name in SOLVERS, f"Unsupported solver: {solver_name}"
    current_path = Path(__file__).parent
    mat = tt.fileio.load(current_path / "test_data")
    # dummy demands for mts
    mat["q_c"][:, :-1] = mat["q_c"][:, 0].reshape(-1, 1) * 0.7
    mat["flh_sinks"][:, :-1] = mat["flh_sinks"][:, 0].reshape(-1, 1) * 0.4
    mat["flh_sinks"][:, 0] = mat["flh_sinks"][:, 0] * 0.6

    settings = tt.settings.load(current_path / "test_data" / "config_mts.yaml")
    model_sets = tt.models.sets.create(mat)
    model = tt.models.multi_timestep.create(
        matrices=mat,
        sets=model_sets,
        regression_inst=R_COSTS,
        regression_losses=R_LOSS,
        economics=settings.economics,
        optimization_mode="forced",
    )

    opt = pyo.SolverFactory(solver_name)
    opt.options["mipgap"] = 0.01
    opt.options["timelimit"] = 3600
    result = opt.solve(model, tee=True)
    assert result.solver.status == pyo.SolverStatus.ok
    assert abs(pyo.value(model.obj)) == approx(683.0, rel=0.02)


def test_mts_eco(request):
    """Main function to run the optimization"""
    solver_name = request.config.getoption("--solver")
    assert solver_name in SOLVERS, f"Unsupported solver: {solver_name}"
    current_path = Path(__file__).parent
    mat = tt.fileio.load(current_path / "test_data")
    mat["q_c"][:, :-1] = mat["q_c"][:, 0].reshape(-1, 1) * 0.7
    mat["flh_sinks"][:, :-1] = mat["flh_sinks"][:, 0].reshape(-1, 1) * 0.4
    mat["flh_sinks"][:, 0] = mat["flh_sinks"][:, 0] * 0.6
    settings = tt.settings.load(current_path / "test_data" / "config_mts.yaml")

    model_sets = tt.models.sets.create(mat)
    model = tt.models.multi_timestep.create(
        matrices=mat,
        sets=model_sets,
        regression_inst=R_COSTS,
        regression_losses=R_LOSS,
        economics=settings.economics,
        optimization_mode="economic",
    )

    opt = pyo.SolverFactory(solver_name)
    opt.options["mipgap"] = 0.01
    opt.options["timelimit"] = 3600

    result = opt.solve(model, tee=True)
    assert result.solver.status == pyo.SolverStatus.ok
    assert pyo.value(model.obj) == approx(-907.7, rel=0.02)
