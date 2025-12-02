"""Module for reading and writing input and output data for the optimization
problem.

The input data has to be stored in parquet files and read with the function
`load`.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from topotherm.utils import find_duplicate_cols

_FILES = {
    "a_i": {"file": "a_i", "dtype": int},
    "a_p": {"file": "a_p", "dtype": int},
    "a_c": {"file": "a_c", "dtype": int},
    "l_i": {"file": "l_i", "dtype": float},
    "q_c": {"file": "q_c", "dtype": float},
    "flh_sinks": {"file": "flh_sinks", "dtype": float},
    "flh_sources": {"file": "flh_sources", "dtype": float},
    "positions": {"file": "positions", "cols": (-2, -1), "dtype": float},
}


def _load_one(path: Path, cfg: dict) -> np.ndarray:
    """Load one parquet file and convert to an array."""
    if "cols" in cfg:
        try:
            df = pd.read_parquet(path, columns=list(cfg["cols"]))
        except Exception:
            df = pd.read_parquet(path).iloc[:, list(cfg["cols"])]
    else:
        df = pd.read_parquet(path)

    arr = df.to_numpy()
    if "dtype" in cfg:
        arr = arr.astype(cfg["dtype"], copy=False)
    return arr


def load(basepath: Path, filenames: dict = None) -> dict:
    """
    Read the input data from the given path and return the matrices.

    Parameters
    ----------
    basepath : pathlib.Path
        Path to the input data.
    filenames : dict, optional
        Dict containing filenames and optional import options". Example:
        "flh_sources": {"file": "flh_sources", "dtype": float}

    Returns
    -------
    dict
        Matrices stored in the following keys:

        - ``a_i`` : Incidence matrix for the pipes (rows: nodes, columns: edges).
        - ``a_p`` : Adjacency matrix for the producers (rows: nodes, columns: supply IDs).
        - ``a_c`` : Adjacency matrix for the consumers (rows: nodes, columns: consumer IDs).
        - ``q_c`` : Heat demand of the consumers for each time step in kW (rows: consumers, columns: time steps).
        - ``l_i`` : Length of edges (rows: edges, columns: -).
        - ``positions`` : ``(x, y)`` coordinates of the nodes in the network (rows: nodes, columns: [x coordinate, y coordinate]).
        - ``flh_sources``: Full load hours of each source for each time step in h/year (rows: sources, columns: time steps)
        - ``flh_sinks``: Full load hours of each sink, for each time steps in h/year (rows: consumers, columns: time steps)
    """
    r = {}
    if filenames is not None:
        _keys = [k for k in filenames.keys() if k not in _FILES.keys()]
        if len(_keys) > 0:
            raise ValueError(
                f"Keys {_keys} passed in filenames are not valid. Accepted keys: {list(_FILES.keys())}"
            )
        for key, val in _FILES.items():
            if key not in filenames.keys():
                filenames[key] = val
    else:
        filenames = _FILES

    for key, val in filenames.items():
        p = basepath / f"{val['file']}.parquet"
        r[key] = _load_one(path=p, cfg=val)

    # @TODO Implement real warnings/errors and implement checks for full load hours
    if (r["a_i"].sum(axis=0).sum() != 0) | (
        np.abs(r["a_i"]).sum(axis=0).sum() / 2 != np.shape(r["a_i"])[1]
    ):
        print("Warning: The structure of A_i is not correct!")
    elif (-r["a_p"].sum(axis=0).sum() != np.shape(r["a_p"])[1]) | (
        np.abs(r["a_p"]).sum(axis=0).sum() != np.shape(r["a_p"])[1]
    ):
        print("Warning: The structure of A_p is not correct!")
    elif (np.abs(r["a_c"]).sum(axis=0).sum() != np.shape(r["a_c"])[1]) | (
        r["a_c"].sum(axis=0).sum() != np.shape(r["a_c"])[1]
    ):
        print("Warning: The structure of A_c is not correct!")
    elif (np.shape(r["a_i"])[0] != np.shape(r["a_p"])[0]) | (
        np.shape(r["a_i"])[0] != np.shape(r["a_c"])[0]
    ):
        print("Warning: Number of nodes doesn't correspond!")
    elif np.shape(r["q_c"])[0] != np.shape(r["a_c"])[1]:
        print("Warning: Length of Q_c doesn't match shape of A_c!")
    elif np.shape(r["l_i"])[0] != np.shape(r["a_i"])[1]:
        print("Warning: Length of Q_c doesn't match shape of A_c!")
    elif np.shape(r["positions"])[0] != np.shape(r["a_i"])[0]:
        print("Warning: Position doesn't match with the number of nodes!")
    elif len(find_duplicate_cols(r["a_i"])) != 0:
        print("Warning: There are duplicate columns in A_i, we took care of it!")
        delete_col = find_duplicate_cols(r["a_i"])
        if r["l_i"][delete_col[0][0]] > r["l_i"][delete_col[0][1]]:
            np.delete(r["l_i"], delete_col[0][0], axis=0)
            np.delete(r["a_i"], delete_col[0][0], axis=1)
        else:
            np.delete(r["l_i"], delete_col[0][1], axis=0)
            np.delete(r["a_i"], delete_col[0][1], axis=1)
    return r
