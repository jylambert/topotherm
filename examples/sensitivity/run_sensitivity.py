# -*- coding: utf-8 -*-
"""
@author: Jerry Lambert (jerry.lambert@tum.de); Amedeo Ceruti (amedeo.ceruti@tum.de)

This file is used to optimize district heating systems with the tool topotherm
of the Chair of Energy Systems. This is the example file is based on the MTS model.
The sensitvity example highlights the  'sensitivity' model, which calculates the network
configurations in steps of 10 % from 10 to 100% of the total heating demand in a district. This
way the costs can be calculated incrementally.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pyomo.environ as pyo

import topotherm as tt
from topotherm.hydraulic import regression_heat_losses, regression_thermal_capacity


DATAPATH = Path(__file__).parent / "data"
OUTPUTPATH = Path(__file__).parent / "results"
PLOTS = True  # plot districts before and after optimization
SOLVER = "gurobi"  # 'gurobi', 'cplex' or 'scip'


def main(filepath, outputpath, plots=False, solver="gurobi", mode="forced"):
    """Main function to run the optimization"""
    # Create output directory if it does not exist
    tt.utils.create_dir(outputpath)

    # Load the district
    mat = tt.fileio.load(filepath)

    if plots:
        f = tt.plotting.district(mat, isnot_init=False)  # Save initial District
        f.savefig(outputpath / "district_initial.svg", bbox_inches="tight")

    # default settings
    settings = tt.settings.Settings()
    settings.economics.heat_price = 0  # no revenue in objective function

    r_thermal_cap = regression_thermal_capacity(settings)
    r_heat_loss = regression_heat_losses(settings, r_thermal_cap)

    min_share = 0.1  # sensitivity starts at 10 % of total heating demand
    max_share = 1  # up to 100 % heating demand
    number_of_iterations = 10  # in ten steps

    share_q_c_tot = np.linspace(min_share, max_share, number_of_iterations)
    iterations = range(number_of_iterations)

    for run in iterations:
        print("============================")
        print(f"Running iteration {run} of {number_of_iterations}")
        print("============================\n")

        variant = "district_sensitivity_" + str(run)

        # Load the district
        mat = tt.fileio.load(filepath)

        model_sets = tt.models.sets.create(mat)
        # define the amount of heat that needs to be supplied, while leaving
        # open which sinks to connect
        model_sets["q_c_tot"] = (
            mat["flh_sinks"] * mat["q_c"]
        ).sum() * share_q_c_tot[run]
        model_sets["lambda_b_previous"] = np.zeros(model_sets["a_i_shape"][1])

        # The following if section is only needed, if the built-up of the district should be consecutive and
        # dependent of one another. Also comment or uncomment constraint in model file.
        if run != 0:
            variant_previous = "district_sensitivity_" + str(run - 1)
            lambda_b_previous = pd.read_parquet(
                outputpath / variant_previous / "lambda_b_orig.parquet")
            model_sets["lambda_b_previous"] = lambda_b_previous.to_numpy().reshape(-1)

        # create output dir for every sub-iteration
        outputpath_sens = outputpath / variant
        tt.utils.create_dir(outputpath_sens)

        model = tt.models.single_timestep.create(
            matrices=mat,
            sets=model_sets,
            regression_inst=r_thermal_cap,
            regression_losses=r_heat_loss,
            economics=settings.economics,
            optimization_mode=mode,
        )

        # Optimization initialization and set options
        opt = pyo.SolverFactory(solver)
        opt.options["mipgap"] = settings.solver.mip_gap
        opt.options["timelimit"] = settings.solver.time_limit
        opt.options["logfile"] = str((outputpath_sens / "optimization.log").resolve())

        result = opt.solve(model, tee=True)

        # Save model results to csv
        tt.utils.model_to_df(model).to_csv(outputpath_sens / "results.csv", sep=";")

        # save solver results
        tt.utils.solver_to_df(result, model).to_csv(outputpath_sens / "solver.csv", sep=";")

        opt_mats = tt.postprocessing.sts(model=model, matrices=mat, settings=settings)

        # iterate over opt_mats and save each matrix as parquet file
        for key, value in opt_mats.items():
            pd.DataFrame(value).to_parquet(outputpath_sens / f"{key}.parquet")

        # Save figure optimized districts
        if plots:
            f = tt.plotting.district(
                opt_mats, diameter=opt_mats["d_i_0"], isnot_init=True
            )
            f.savefig(outputpath_sens / "district_optimal.svg",
                bbox_inches="tight",
            )


if __name__ == "__main__":
    main(
        filepath=DATAPATH,
        outputpath=OUTPUTPATH,
        plots=PLOTS,
        solver=SOLVER,
        mode="sensitivity",
    )
    print(f"Finished {OUTPUTPATH}")
