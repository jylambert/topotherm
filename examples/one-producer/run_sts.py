# -*- coding: utf-8 -*-
"""
@author: Jerry Lambert (jerry.lambert@tum.de); Amedeo Ceruti
(amedeo.ceruti@tum.de)

This file is used to optimize district heating systems with the tool topotherm
of the Chair of Energy Systems.
"""

import os

import pandas as pd
import pyomo.environ as pyo

import topotherm as tt


DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_sts')
OUTPUTPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'results', 'sts_forced')

# regression coefficients for thermal capacity and heat losses
REGRESSION = 'regression.csv'
PLOTS = True  # plot districts before and after optimization
SOLVER = 'gurobi'  # 'gurobi', 'cplex' or 'scip'


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


def main(filepath, outputpath, plots=True, solver='gurobi', mode='economic'):
    """Main function to run the optimization"""
    # Create output directory if it does not exist
    tt.utils.create_dir(outputpath)

    # Load the district
    mat = tt.fileio.load(filepath)
    mat['q_c'] = mat['q_c'] / 1000  # convert to kW

    if plots:
        f = tt.plotting.district(mat, isnot_init=False) # Save initial District
        f.savefig(os.path.join(outputpath, 'district_initial.svg'),
                  bbox_inches='tight')

    # regression
    r_thermal_cap, r_heat_loss = read_regression(
        os.path.join(filepath, REGRESSION),
        0)

    # import settings
    settings = tt.settings.load(os.path.join(filepath, 'config.yaml'))

    # modify either in code or in the config file
    settings.economics.source_c_inv = [0.]  # no investment costs for sources

    model_sets = tt.sets.create(mat)
    model = tt.single_timestep.model(
        matrices=mat,
        sets=model_sets,
        regression_inst=r_thermal_cap,
        regression_losses=r_heat_loss,
        economics=settings.economics,
        optimization_mode=mode)

    # Optimization initialization
    opt = pyo.SolverFactory(solver)
    opt.options['mipgap'] = settings.solver.mip_gap
    opt.options['timelimit'] = settings.solver.time_limit
    opt.options['logfile'] = os.path.join(outputpath, 'optimization.log')

    result = opt.solve(model, tee=True)

    assert result.solver.termination_condition == pyo.TerminationCondition.optimal, \
        f"Optimization failed with termination condition {result.solver.termination_condition}"

    # Save model results to csv
    dfres = tt.utils.model_to_df(model)
    dfres.to_csv(os.path.join(outputpath, 'results.csv'), sep=';')

    # save solver results
    dfsol = tt.utils.solver_to_df(result, model)
    dfsol.to_csv(os.path.join(outputpath, 'solver.csv'), sep=';')

    opt_mats = tt.postprocessing.sts(
        model=model,
        matrices=mat,
        settings=settings)

    node_data, edge_data = tt.postprocessing.to_dataframe(opt_mats, mat)
    node_data.to_csv(os.path.join(outputpath, 'node_data.csv'), sep=';')
    edge_data.to_csv(os.path.join(outputpath, 'edge_data.csv'), sep=';')

    # iterate over opt_mats and save each matrix as parquet file
    for key, value in opt_mats.items():
        pd.DataFrame(value).to_parquet(os.path.join(outputpath, key + '.parquet'))

    # Save figure optimized districts
    if plots:
        f = tt.plotting.district(opt_mats, diameter=opt_mats['d_i_0'],
                                    isnot_init=True)
        f.savefig(os.path.join(outputpath, 'district_optimal.svg'),
                    bbox_inches='tight')
    
    # create networkx graph object
    network = tt.postprocessing.to_networkx_graph(opt_mats)

if __name__ == '__main__':
    main(filepath=os.path.join(DATAPATH), outputpath=os.path.join(OUTPUTPATH),
         plots=PLOTS, solver=SOLVER, mode='forced')
    print(f'Finished {OUTPUTPATH}')
