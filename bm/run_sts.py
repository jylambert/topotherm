# -*- coding: utf-8 -*-
"""
@author: Jerry Lambert (jerry.lambert@tum.de); Amedeo Ceruti (amedeo.ceruti@tum.de)

This file is used to optimize district heating systems with the tool topotherm
of the Chair of Energy Systems.
"""

import os

import pandas as pd
import pyomo.environ as pyo
import matplotlib.pyplot as plt
import networkx as nx

import topotherm as tt
# from topotherm.settings import Settings
from topotherm.model.create_sets import main as create_sets
from topotherm.model.sts import main as sts


DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
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


def main(filepath, outputpath, plots=True, solver='gurobi', mode='forced'):
    """Main function to run the optimization"""
    # Create output directory if it does not exist
    tt.utils.create_dir(outputpath)

    # Load the district
    mat = tt.fileio.load(filepath)

    if plots:
        f = tt.plotting.district(mat, isnot_init=False) # Save initial District
        f.savefig(os.path.join(outputpath, 'district_initial.svg'),
                  bbox_inches='tight')

    # regression
    r_thermal_cap, r_heat_loss = read_regression(
        os.path.join(filepath, REGRESSION),
        0)

    # import settings
    settings = tt.settings.load('./data/config.yaml')
    print(settings)
    settings.economics.source_c_inv = [0.]  # no investment costs for sources
    settings.economics.source_flh = [2500.]  # full load hours
    settings.economics.consumers_flh = 2500.  # full load hours

    model_sets = create_sets(mat)
    model = sts(
        mat,
        model_sets,
        r_thermal_cap,
        r_heat_loss,
        settings.economics,
        mode)

    # Optimization initialization
    opt = pyo.SolverFactory(solver)
    opt.options['mipgap'] = settings.solver.mip_gap
    opt.options['timelimit'] = settings.solver.time_limit
    opt.options['logfile'] = os.path.join(outputpath, 'optimization.log')
    #opt.options['Seed'] = 56324978

    result = opt.solve(model, tee=True)

    # Save model results to csv
    dfres = tt.utils.model_to_df(model)
    dfres.to_csv(os.path.join(outputpath, 'results.csv'), sep=';')

    # save solver results
    dfsol = tt.utils.solver_to_df(result, model, solver=solver)
    dfsol.to_csv(os.path.join(outputpath, 'solver.csv'), sep=';')

    opt_mats = tt.postprocessing.sts(
        model=model,
        matrices=mat,
        settings=settings)

    # Debug prints to check shapes
    for key, value in opt_mats.items():
        print(f"{key} shape: {value.shape}")

    # iterate over opt_mats and save each matrix as parquet file
    for key, value in opt_mats.items():
        pd.DataFrame(value).to_parquet(os.path.join(outputpath, key + '.parquet'))

    # Save figure optimized districts
    if plots:
        f = tt.plotting.district(opt_mats, diameter=opt_mats['d_i_0'], isnot_init=True)
        f.savefig(os.path.join(outputpath, 'district_optimal.svg'), bbox_inches='tight')

    network = tt.postprocessing.to_networkx_graph(opt_mats)

    print(len(network.edges))
    print(mat['a_i'].shape[1])
    print(len(network.nodes))
    print(mat['a_i'].shape[0])

    fig, ax = plt.subplots(figsize=(20, 20), layout='constrained')
    node_colors = []
    node_label = []
    node_id = {}
    node_pos = []
    edges_p = []
    edges_label = {}

    for node in network.nodes(data=True):
        node_colors.append(node[1]['color'])
        node_label.append(node[1]['type'])
        node_id[node[0]] = str(node[0])
        node_pos.append([node[1]['x'], node[1]['y']])

    for edge in network.edges(data=True):
        edges_p.append(edge[2]['p'])
        edges_label[(edge[0], edge[1])] = str(edge[0]) + ' -> ' + str(edge[1])

    nx.draw_networkx_edges(network, pos=node_pos,
                            edgelist=network.edges, width=edges_p, ax=ax,
                            label=edges_label, alpha=0.3, edge_color='grey')
    nx.draw_networkx_nodes(network, pos=node_pos, node_color=node_colors,
                            ax=ax, label=node)

    nx.draw_networkx_labels(network, pos=node_pos, labels=node_id, ax=ax)#
    nx.draw_networkx_edge_labels(network, pos=node_pos, edge_labels=edges_label, ax=ax)

    fig.show()
    plt.pause(120)
    fig.savefig(os.path.join(outputpath, 'networkx.svg'), bbox_inches='tight')
    # close all figures
    plt.close('all')

if __name__ == '__main__':
    main(filepath=os.path.join(DATAPATH), outputpath=os.path.join(OUTPUTPATH),
         plots=PLOTS, solver=SOLVER, mode='forced')
    print(f'Finished {OUTPUTPATH}')
