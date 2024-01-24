# -*- coding: utf-8 -*-
"""
@author: Jerry Lambert (jerry.lambert@tum.de); Amedeo Ceruti (amedeo.ceruti@tum.de)

This file is used to optimize district heating systems with the tool topotherm
of the Chair of Energy Systems. This is the example file for the MTS-Easy model.
"""

import os

import pandas as pd
import pyomo.environ as pyo

import topotherm as tt
import topotherm.precalculation_hydraulic as precalc
import topotherm.settings as settings


RUNID = ['Example']  # name of the district to be optimized
DATAPATH = 'data/'  # path to the data
OUTPUTPATH = 'new_results/mts_eco/'  # output path for results
REGRESSION = 'regression.csv'  # regression coefficients for thermal capacity and heat losses
TIMESERIES = 'timeseries.csv'  # timeseries for heat scaling
PLOTS = True  # save plots of the district
SOLVER = 'gurobi' # 'gurobi' or 'cbc'


def read_regression(path, index):
    """Read the regression coefficients for the thermal capacity and heat losses from csv file.
    
    """
    # read df and force floats
    df = pd.read_csv(path, sep=',', index_col=0, header=0, dtype=float)
    r_thermal_cap = dict(power_flow_max_kW=[df.loc[id, "power_flow_max_kW"]],
                         params={"a": df.loc[id, "capacity_a"],
                                 "b": df.loc[id, "capacity_b"]},
                         power_flow_max_partload=df.loc[id, "power_flow_max_partload"])
    r_heat_loss = dict(params={"a": df.loc[id, "heat_loss_a"],
                                 "b": df.loc[id, "heat_loss_b"]
                                 })
    return r_thermal_cap, r_heat_loss

def regressions(temperatures):
    """Calculate the regression coefficients for the thermal capacity and heat losses"""
    # thermal capacity regression and max power as well as part loads
    r_thermal_cap = precalc.regression_thermal_capacity(temperatures)
    # heat loss regression
    r_heat_loss = precalc.regression_heat_losses(temperatures, r_thermal_cap)
    return r_thermal_cap, r_heat_loss

def main(runid):
    """Main function to run the optimization"""
    # Create output directory if it does not exist
    resultspath = os.path.join(OUTPUTPATH, runid)
    filepath = os.path.join(DATAPATH, runid)

    tt.utils.create_dir(resultspath)

    # init temperatures
    temps = precalc.init_temperatures()

    # Load the district
    mat = tt.fileio.load(filepath)
    # read in demand profile
    timeseries = pd.read_csv(os.path.join(filepath, TIMESERIES),
                             sep=';', index_col=0, header=0).iloc[7:9, :].values.squeeze() # 4:9
    # create dummy profile, q_c should already contain the timeseries of all consumer demands
    mat['q_c'] = mat['q_c'] * timeseries  # convert to timeseries
    
    if PLOTS:
        f = tt.plotting.district(mat, isnot_init=False) # Save initial District
        f.savefig(os.path.join(resultspath, 'district_initial.svg'), bbox_inches='tight')

    # regression
    #r_thermal_cap, r_heat_loss = read_regression(os.path.join(DATAPATH, runid, REGRESSION), 0)
    r_thermal_cap, r_heat_loss = regressions(temps)

    # -------------------------------- Create Model --------------------------------
    model_sets = tt.model.create_sets(mat)
    model = tt.model.mts(mat, model_sets, r_thermal_cap, r_heat_loss, "eco",
                              flh_scaling=timeseries.sum())

    # -------------------------------- Initialize Optimization --------------------------------
    # Optimization initialization
    opt = pyo.SolverFactory(SOLVER)
    opt.options['mipgap'] = settings.OptSettings.mip_gap
    opt.options['timelimit'] = settings.OptSettings.time_limit
    opt.options['logfile'] = os.path.join(resultspath, 'optimization.log')
    opt.options['Seed'] = 56324978

    # -------------------------------- Solve the Model --------------------------------
    result = opt.solve(model, tee=True)

    # -------------------------------- Process the results --------------------------------
    # Save model results to csv
    dfres = tt.utils.model_to_df(model)
    dfres.to_csv(os.path.join(resultspath, 'results.csv'), sep=';')

    # save solver results
    dfsol = tt.utils.solver_to_df(result, model, solver='gurobi')
    dfsol.to_csv(os.path.join(resultspath, 'solver.csv'), sep=';')

    opt_mats = tt.postprocessing.postprocess(model, mat, model_sets, "mts", temperatures=temps)

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
