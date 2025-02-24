# -*- coding: utf-8 -*-
"""
@author: Jerry Lambert (jerry.lambert@tum.de); Amedeo Ceruti
(amedeo.ceruti@tum.de)

This file is used to optimize district heating systems with the tool topotherm
of the Chair of Energy Systems.
"""

import os

import numpy as np
import pandas as pd
import pyomo.environ as pyo

import topotherm as tt
from topotherm.settings import Settings
from topotherm.precalculation_hydraulic import regression_thermal_capacity, regression_heat_losses

DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data_sts')
OUTPUTPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'results', 'sts_forced')

# regression coefficients for thermal capacity and heat losses
REGRESSION = 'regression.csv'
PLOTS = True  # plot districts before and after optimization
SOLVER = 'gurobi'  # 'gurobi', 'cplex' or 'scip'


def main(filepath, outputpath, plots=True, solver='gurobi', mode='forced'):
    """Main function to run the optimization"""
    # Create output directory if it does not exist
    tt.utils.create_dir(outputpath)

    # Load the district
    mat = tt.fileio.load(filepath)

    if plots:
        f = tt.plotting.district(mat, isnot_init=False)  # Save initial District
        f.savefig(os.path.join(outputpath, 'district_initial.svg'), bbox_inches='tight')

    # default settings
    settings = Settings()
    settings.piping.diameter = (
        0.0229, 0.0291, 0.0372, 0.0431, 0.0545, 0.0703, 0.0825, 0.1071, 0.1325, 0.1603, 0.2101, 0.263, 0.3127, 0.3444,
        0.3938, 0.4444, 0.4954
    )

    settings.piping.outer_diameter = (
        0.09, 0.09, 0.11, 0.11, 0.125, 0.14, 0.16, 0.2, 0.225, 0.25, 0.315, 0.4, 0.45, 0.5, 0.56, 0.63, 0.71
    )

    settings.piping.cost = (
        390, 400, 425, 459, 490, 535, 599, 661, 749, 886, 1170, 1016, 1192, 1395, 1748, 2070, 2150
    )

    settings.piping.number_diameters = 17
    settings.piping.max_pr_loss = 250
    settings.piping.roughness = 0.01e-3

    settings.temperatures.supply = 70
    settings.temperatures.return_ = 40
    settings.temperatures.ambient = -20

    settings.solver.mip_gap = 2e-2
    settings.solver.time_limit = 10000

    elec_price_low = 0.080

    cop_river = 1.86
    cop_air = 1.97

    inv_river = 656.625
    inv_air = 598

    settings.economics.source_price = (elec_price_low/cop_river, elec_price_low/cop_air,
                                       elec_price_low/cop_river, elec_price_low/cop_air,
                                       elec_price_low/cop_air, elec_price_low/cop_river,
                                       elec_price_low/cop_river, elec_price_low/cop_air,)

    settings.economics.source_c_inv = (inv_river, inv_air,
                                       inv_river, inv_air,
                                       inv_air, inv_river,
                                       inv_river, inv_air)

    settings.economics.source_lifetime = (20)#, 20, 20, 20, 20, 20, 20, 20)
    settings.economics.source_c_irr = (0.08)#, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08)

    settings.economics.pipes_lifetime = 20
    settings.economics.pipes_c_irr = 0.08

    settings.economics.heat_price = 0

    settings.economics.source_max_power = (100000)#, 20000,
                                           #100000, 1500,
                                           #20000, 100000,
                                           #100000, 20000)

    # regression
    r_thermal_cap = regression_thermal_capacity(settings)
    r_heat_loss = regression_heat_losses(settings, r_thermal_cap)

    min_share = 0.1
    max_share = 1
    number_of_iterations = 10

    share_q_c_tot = np.linspace(min_share, max_share, number_of_iterations)
    iterations = range(number_of_iterations)

    for run in iterations:
        print('============================')
        print(f'Running iteration {run}')
        print('============================\n')

        variant = 'district_sensitivity_' + str(run)

        # Load the district
        mat = tt.fileio.load(filepath)

        model_sets = tt.sets.create(mat)
        model_sets['q_c_tot'] = (mat['flh_consumer'] * mat['q_c']).sum() * share_q_c_tot[run]


        outputpath_sens = os.path.join(outputpath, variant)
        tt.utils.create_dir(outputpath_sens)

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
        opt.options['logfile'] = os.path.join(outputpath_sens, 'optimization.log')

        result = opt.solve(model, tee=True)

        # Save model results to csv
        dfres = tt.utils.model_to_df(model)
        dfres.to_csv(os.path.join(outputpath_sens, 'results.csv'), sep=';')

        # save solver results
        dfsol = tt.utils.solver_to_df(result, model, solver=solver)
        dfsol.to_csv(os.path.join(outputpath_sens, 'solver.csv'), sep=';')

        opt_mats = tt.postprocessing.sts(
            model=model,
            matrices=mat,
            settings=settings)

        # iterate over opt_mats and save each matrix as parquet file
        for key, value in opt_mats.items():
            pd.DataFrame(value).to_parquet(os.path.join(outputpath_sens, key + '.parquet'))

        # Save figure optimized districts
        if plots:
            f = tt.plotting.district(opt_mats, diameter=opt_mats['d_i_0'], isnot_init=True)
            f.savefig(os.path.join(outputpath_sens, 'district_optimal.svg'), bbox_inches='tight')


if __name__ == '__main__':
    main(filepath=os.path.join(DATAPATH), outputpath=os.path.join(OUTPUTPATH),
         plots=PLOTS, solver=SOLVER, mode='economic')
    print(f'Finished {OUTPUTPATH}')
