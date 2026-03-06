# -*- coding: utf-8 -*-
import os
import shutil

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from pyomo.core.base.componentuid import ComponentUID


def create_dir(path: str) -> None:
    """
    Creates a directory if it does not exist and deletes old results.

    Parameters
    ----------
    path : str
        Path to directory.

    Returns
    -------
    None
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


def solver_to_df(result, model):
    """
    Returns solver results in a dataframe. This needs to be adapted to the
    solver output (gurobi vs cplex have different naming conventions).

    Parameters
    ----------
    result : dict
        Solver result dictionary produced by Pyomo.
    model : pyomo.core.base.PyomoModel
        Solved Pyomo model instance.

    Returns
    -------
    pandas.DataFrame or dict
        DataFrame with solver statistics (with columns 'Value' and 'Unit'),
        or raw solver output dict if parsing fails.
    """

    # Useful links:
    # https://stackoverflow.com/questions/45034035/meaning-of-time-in-pyomos-results-json

    try:
        slvr_res = result["Solver"][0]
    except (KeyError, IndexError, TypeError):
        print("Solver output not as expected. Check the solver output.")
        return result.get("Solver", {})

    try:
        # Build data dictionary for DataFrame
        data = {}

        # Termination condition
        termination = slvr_res.get("Termination condition", "Unknown")
        data["Termination condition"] = {"Value": termination, "Unit": "-"}

        # Wall time
        wall_time = slvr_res.get("Wall time")
        if wall_time is not None:
            data["Wall Time"] = {"Value": wall_time, "Unit": "s"}

        # Objective value
        try:
            obj_value = pyo.value(model.obj)
            data["Objective"] = {"Value": obj_value, "Unit": "eur/y"}
        except (AttributeError, ValueError):
            data["Objective"] = {"Value": None, "Unit": "eur/y"}

        # Create DataFrame from the data dictionary
        dfslvr = pd.DataFrame(data).T
        return dfslvr

    except Exception as e:
        print(f"Solver output not as expected. Error: {e}")
        return slvr_res


def model_to_df(model):
    """
    Converts a solved pyomo model to a pandas dataframe.

    Parameters
    ----------
    model : pyomo.core.base.PyomoModel
        Solved Pyomo model.

    Returns
    -------
    pandas.Series
        Series containing variable, parameter, and objective values.
    """

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


def find_duplicate_cols(data: np.ndarray, minoccur: int = 2) -> list:
    """
    Find duplicate columns in a numpy array.

    Parameters
    ----------
    data : np.ndarray
        Data to check for duplicates.
    minoccur : int
        Minimum number of occurrences to be considered a duplicate.

    Returns
    -------
    list
        List of indices of duplicate columns.
    """
    ind = np.lexsort(data)
    diff = np.any(data.T[ind[1:]] != data.T[ind[:-1]], axis=1)
    edges = np.where(diff)[0] + 1
    result = np.split(ind, edges)
    result = [group for group in result if len(group) >= minoccur]
    return result
