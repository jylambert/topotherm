# -*- coding: utf-8 -*-
"""
@author: Jerry Lambert (jerry.lambert@tum.de); Amedeo Ceruti (amedeo.ceruti@tum.de)

This file is used to optimize district heating systems with the tool topotherm
of the Chair of Energy Systems.
"""

import os

import pandas as pd
import pyomo.environ as pyo

import topotherm as tt
import topotherm.precalculation_hydraulic as precalc
import topotherm.settings as settings


RUNID = ['Example']
DATAPATH = 'data/'
OUTPUTPATH = 'new_results/'
PLOTS = True


def regressions(temperatures):
    # thermal capacity regression and max power as well as part loads
    r_thermal_cap = precalc.regression_thermal_capacity(temperatures)
    # heat loss regression
    r_heat_loss = precalc.regression_heat_losses(temperatures, r_thermal_cap)
    return r_thermal_cap, r_heat_loss


def main(runid):
    # Create output directory if it does not exist
    resultspath = os.path.join(OUTPUTPATH, runid)
    filepath = os.path.join(DATAPATH, runid)

    tt.utils.create_dir(resultspath)

    # Load the district
    mat = tt.fileio.load(filepath)

    if PLOTS:
        f = tt.plotting.district(mat, isnot_init=False) # Save initial District
        f.savefig(os.path.join(resultspath, 'district_initial.svg'), bbox_inches='tight')

    # Run STS
    # @TODO: put this into separate notebook? and import dicts from csv
    # init temperatures
    temps = precalc.init_temperatures()

    # regression
    r_thermal_cap, r_heat_loss = regressions(temps)

    # -------------------------------- Create Model --------------------------------
    model_sets = tt.model.create_sets(mat)
    model = tt.model.sts(mat, model_sets, r_thermal_cap, r_heat_loss)

    # -------------------------------- Initialize Optimization --------------------------------
    # Optimization initialization
    opt = pyo.SolverFactory('gurobi')
    opt.options['mipgap'] = settings.OptSettings.mip_gap
    opt.options['timelimit'] = settings.OptSettings.time_limit
    opt.options['logfile'] = os.path.join(resultspath, 'optimization.log')
    #opt.options['Seed'] = 56324978

    # -------------------------------- Solve the Model --------------------------------
    result = opt.solve(model, tee=True)

        # -------------------------------- Process the results --------------------------------
    # Save model results to csv
    dfres = tt.utils.model_to_df(model)
    dfres.to_csv(os.path.join(resultspath, 'results.csv'), sep=';')

    # save solver results
    dfsol = tt.utils.solver_to_df(result, model, solver='gurobi')
    dfsol.to_csv(os.path.join(resultspath, 'solver.csv'), sep=';')

    opt_mats = tt.postprocessing.postprocess(model, mat, model_sets, temperatures=temps)

    # iterate over opt_mats and save each matrix as parquet file
    for key, value in opt_mats.items():
        pd.DataFrame(value).to_parquet(os.path.join(resultspath, key + '.parquet'))

    # Save figure optimized districts
        if PLOTS:
            f = tt.plotting.district(opt_mats, diameter=opt_mats['d_i_0'], isnot_init=True)
            f.savefig(os.path.join(resultspath, 'district_optimal.svg'), bbox_inches='tight')

if __name__ == '__main__':
    for name in RUNID:
        main(name)
        print(f'Finished {name}')
