# -*- coding: utf-8 -*-
import os
import shutil

import pandas as pd
import pyomo.environ as pyo
from pyomo.core.base.componentuid import ComponentUID


def create_dir(path: str) -> None:
    """Creates a directory if it does not exist and deletes old results.

    Args:
        path (str): path to directory
    """
    # create results directory
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        # delete old results
        for f in os.listdir(path):
            if os.path.isdir(os.path.join(path, f)):
                shutil.rmtree(os.path.join(path, f))
            elif os.path.isfile(os.path.join(path, f)):
                os.remove(os.path.join(path, f))
    return


def solver_to_df(result, model, solver):
    """Returns solver results in a dataframe. This needs to be adapted to the
    solver output (gurobi vs cplex have different naming conventions)"""

    # Useful links:
    # https://stackoverflow.com/questions/45034035/meaning-of-time-in-pyomos-results-json

    dfslvr = pd.DataFrame()
    slvr_res = result['Solver'][0]
    if solver == "cplex":
        dfslvr.loc['Termination condition', 0] = slvr_res['Termination condition']
        dfslvr.loc['Termination condition', 'unit'] = '-'
        #dfslvr.loc['Time', 0] = slvr_res['User time']
        #dfslvr.loc['Time', 'unit'] = 's'
        dfslvr.loc['Objective', 0] = result['Solution'][0]['Objective']['obj']['Value']
        dfslvr.loc['Objective', 'unit'] = '-'
    elif solver == "gurobi":
        dfslvr.loc['Termination condition', 0] = slvr_res['Termination condition']
        dfslvr.loc['Termination condition', 'unit'] = '-'
        dfslvr.loc['Wall time', 0] = slvr_res['Wall time']
        dfslvr.loc['Wall time', 'unit'] = 's'
        #dfslvr.loc['Time', 0] = slvr_res['Time']
        #dfslvr.loc['Time', 'unit'] = 's'
        dfslvr.loc['Objective', 0] = pyo.value(model.obj)
        dfslvr.loc['Objective', 'unit'] = '-'
    elif solver == "scip":
        dfslvr.loc['Termination condition', 0] = slvr_res['Termination condition']
        dfslvr.loc['Termination condition', 'unit'] = '-'
        dfslvr.loc['Termination messag', 0] = slvr_res['Termination condition']
        dfslvr.loc['Termination messag', 'unit'] = '-'
        dfslvr.loc['Time', 0] = slvr_res['Wallclock time']
        dfslvr.loc['Time', 'unit'] = 's'
        dfslvr.loc['Objective', 0] = pyo.value(model.obj)
        dfslvr.loc['Objective', 'unit'] = '-'
    else:
        raise NotImplementedError(f"Solver {solver} not implemented")
    return dfslvr


def model_to_df(model):
    """Converts a solved pyomo model to a pandas dataframe."""
    solution = {}

    # generate cuid names efficiently in bulk
    # labels = generate_cuid_names(model)
    labels = ComponentUID.generate_cuid_string_map(model)

    for var in model.component_data_objects(pyo.Var, active=True):
        solution[labels[var]] = pyo.value(var)
    for prm in model.component_data_objects(pyo.Param, active=True):
        solution[labels[prm]] = pyo.value(prm)
    for obj in model.component_data_objects(pyo.Objective, active=True):
        solution[labels[obj]] = pyo.value(obj)

    df = pd.Series(solution)
    return df
