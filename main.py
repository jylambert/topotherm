# -*- coding: utf-8 -*-
"""
@author: Jerry Lambert (jerry.lambert@tum.de) and Amedeo Ceruti (amedeo.ceruti@tum.de)

This file is used to optimize district heating systems with the tool topotherm of the Chair of Energy Systems.
"""

from src.topotherm_sts_vertices import topotherm_sts_vertices
from src.topotherm_mts_easy_vertices import topotherm_mts_easy_vertices
from src.topotherm_mts_vertices import topotherm_mts_vertices


example_case = ["Example"]

file_path_opt = 'data/'


for run in example_case:
    print('============================')
    print(f'Running {run}')
    print('============================\n')
    topotherm_sts_vertices(run, 1, file_path_opt, "results/sts", "topotherm_sts_vertices", 'sts')
    topotherm_mts_easy_vertices(run, 1, file_path_opt, "results/sts", "topotherm_mts_easy_vertices", 'sts')
    topotherm_mts_vertices(run, 1, file_path_opt, "results/sts", "topotherm_mts_vertices", 'sts')

