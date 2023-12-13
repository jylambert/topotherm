# -*- coding: utf-8 -*-
import pickle
import numpy as np
import pandas as pd


def load_district(file_name, file_path):
    """Input: file_name and file_path
    Returns: Matrices A_i A_p and A_c, Heat Demand, Length of edges, and positions"""
    def duplicate_columns(data, minoccur=2):
        ind = np.lexsort(data)
        diff = np.any(data.T[ind[1:]] != data.T[ind[:-1]], axis=1)
        edges = np.where(diff)[0] + 1
        result = np.split(ind, edges)
        result = [group for group in result if len(group) >= minoccur]
        return result

    a_i = pd.read_parquet(file_path + file_name + '/' + file_name + '_A_i.parquet', engine='pyarrow').to_numpy()
    a_p = pd.read_parquet(file_path + file_name + '/' + file_name + '_A_p.parquet', engine='pyarrow').to_numpy()
    a_c = pd.read_parquet(file_path + file_name + '/' + file_name + '_A_c.parquet', engine='pyarrow').to_numpy()
    length = pd.read_parquet(file_path + file_name + '/' + file_name + '_L_i.parquet', engine='pyarrow').to_numpy().astype(float)
    q_c = (pd.read_parquet(file_path + file_name + '/' + file_name + '_Q_c.parquet', engine='pyarrow').to_numpy().astype(float)) / 1000 #Example data is in W, optimization in kW
    position = pd.read_parquet(file_path + file_name + '/' + file_name + '_rel_positions.parquet', engine='pyarrow').loc[:, 'x_rel':'y_rel'].to_numpy().astype(float)

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
        print("Warning: There are duplicate columns in A_i, but we took care of it!")
        delete_col = duplicate_columns(a_i)
        if length[delete_col[0][0]] > length[delete_col[0][1]]:
            np.delete(length, delete_col[0][0], axis=0)
            np.delete(a_i, delete_col[0][0], axis=1)
        else:
            np.delete(length, delete_col[0][1], axis=0)
            np.delete(a_i, delete_col[0][1], axis=1)
    return a_i, a_p, a_c, q_c, length, position
