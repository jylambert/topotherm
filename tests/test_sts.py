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


def test_sts_forced(request):
    """Main function to run the optimization"""
    solver_name = request.config.getoption("--solver")
    assert solver_name in SOLVERS, f"Unsupported solver: {solver_name}"
    # Load the district
    datapath = Path(__file__).resolve().parent / "test_data"
    mat = tt.fileio.load(datapath)

    # import settings
    settings = tt.settings.load(datapath / "config.yaml")

    model_sets = tt.models.sets.create(mat)
    model = tt.models.single_timestep.create(
        matrices=mat,
        sets=model_sets,
        regression_inst=R_COSTS,
        regression_losses=R_LOSS,
        economics=settings.economics,
        optimization_mode="forced",
    )

    # Optimization initialization
    opt = pyo.SolverFactory(solver_name)
    opt.options["mipgap"] = 0.0001

    result = opt.solve(model, tee=True)
    assert result.solver.status == pyo.SolverStatus.ok

    tt.postprocessing.sts(model, mat, settings)
    assert pyo.value(model.obj) == approx(-990.8, rel=0.01)


def test_sts_eco(request):
    """Main function to run the optimization"""
    solver_name = request.config.getoption("--solver")
    assert solver_name in SOLVERS, f"Unsupported solver: {solver_name}"

    datapath = Path(__file__).resolve().parent / "test_data"
    mat = tt.fileio.load(datapath)
    settings = tt.settings.load(datapath / "config.yaml")

    model_sets = tt.models.sets.create(mat)
    model = tt.models.single_timestep.create(
        matrices=mat,
        sets=model_sets,
        regression_inst=R_COSTS,
        regression_losses=R_LOSS,
        economics=settings.economics,
        optimization_mode="economic",
    )

    # Optimization initialization
    opt = pyo.SolverFactory(solver_name)
    opt.options["mipgap"] = 0.0001

    result = opt.solve(model, tee=True)
    assert result.solver.status == pyo.SolverStatus.ok
    assert pyo.value(model.obj) == approx(-1436, rel=0.01)
    tt.postprocessing.sts(model, mat, settings)
