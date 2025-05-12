"""Module for reading and writing input and output data for the optimization
problem.

The input data has to be stored in parquet files and read with the function
`load`.
"""

import os

import numpy as np
import pandas as pd


def load(path: os.PathLike) -> dict:
    """Read the input data from the given path and return the matrices
        * A_i: Incidence matrix for the pipes with rows: nodes, columns: edges.
        * A_p: Adjacency matrix for the producers with rows: nodes, columns:
            supply ids.
        * A_c: Adjacency matrix for the consumers with rows: nodes, columns:
            consumer ids.
        * Q_c: Heat demand of the consumers in W.
        * L_i: Length of edges
        * rel_positions: x, y coordinates of the nodes in the network.
    
    Args:
        path (str or os.path): Path to the input data.
    
    Returns:
        dict: Matrices stored in keys `a_i`, `a_p`, `a_c`, `q_c`, `l_i`,
            and `position`.
    """
    def duplicate_columns(data: np.ndarray, minoccur: int = 2) -> list:
        """Find duplicate columns in a numpy array.
        
        Args:
            data (np.array): Data to check for duplicates.
            minoccur (int): Minimum number of occurrences to be considered a
                duplicate.
        
        Returns:
            list: List of indices of duplicate columns.
        """
        ind = np.lexsort(data)
        diff = np.any(data.T[ind[1:]] != data.T[ind[:-1]], axis=1)
        edges = np.where(diff)[0] + 1
        result = np.split(ind, edges)
        result = [group for group in result if len(group) >= minoccur]
        return result

    # Read the matrices as parquet files
    a_i = pd.read_parquet(os.path.join(path, 'A_i.parquet')).to_numpy()
    a_p = pd.read_parquet(os.path.join(path, 'A_p.parquet')).to_numpy()
    a_c = pd.read_parquet(os.path.join(path, 'A_c.parquet')).to_numpy()
    length = pd.read_parquet(
        os.path.join(path, 'L_i.parquet')).to_numpy().astype(float)
    q_c = (pd.read_parquet(
        os.path.join(path, 'Q_c.parquet')
        ).to_numpy().astype(float))
    flh_consumer = pd.read_parquet(os.path.join(path, 'flh_consumer.parquet')).to_numpy()
    flh_source = pd.read_parquet(os.path.join(path, 'flh_source.parquet')).to_numpy()
    position = pd.read_parquet(
        os.path.join(path, 'rel_positions.parquet')
        ).iloc[:,  [-2,-1]].to_numpy().astype(float)

    # @TODO Implement real warnings/errors and implement checks for full load hours
    if (a_i.sum(axis=0).sum() != 0) | (np.abs(a_i).sum(axis=0).sum()/2 != np.shape(a_i)[1]):
        print("Warning: The structure of A_i is not correct!")
    elif (-a_p.sum(axis=0).sum() != np.shape(a_p)[1]) | (np.abs(a_p).sum(axis=0).sum() != np.shape(a_p)[1]):
        print("Warning: The structure of A_p is not correct!")
    elif (np.abs(a_c).sum(axis=0).sum() != np.shape(a_c)[1]) | (a_c.sum(axis=0).sum() != np.shape(a_c)[1]):
        print("Warning: The structure of A_c is not correct!")
    elif (np.shape(a_i)[0] != np.shape(a_p)[0]) | (np.shape(a_i)[0] != np.shape(a_c)[0]):
        print("Warning: Number of nodes doesn't correspond!")
    elif np.shape(q_c)[0] != np.shape(a_c)[1]:
        print("Warning: Length of Q_c doesn't match shape of A_c!")
    elif np.shape(length)[0] != np.shape(a_i)[1]:
        print("Warning: Length of Q_c doesn't match shape of A_c!")
    elif np.shape(position)[0] != np.shape(a_i)[0]:
        print("Warning: Position doesn't match with the number of nodes!")
    elif len(duplicate_columns(a_i)) != 0:
        print("Warning: There are duplicate columns in A_i, we took care of it!")
        delete_col = duplicate_columns(a_i)
        if length[delete_col[0][0]] > length[delete_col[0][1]]:
            np.delete(length, delete_col[0][0], axis=0)
            np.delete(a_i, delete_col[0][0], axis=1)
        else:
            np.delete(length, delete_col[0][1], axis=0)
            np.delete(a_i, delete_col[0][1], axis=1)

    r = {}
    r['a_i'] = a_i
    r['a_p'] = a_p
    r['a_c'] = a_c
    r['q_c'] = q_c
    r['flh_consumer'] = flh_consumer
    r['flh_source'] = flh_source
    r['l_i'] = length
    r['position'] = position
    return r
