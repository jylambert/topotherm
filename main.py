# -*- coding: utf-8 -*-
"""
@author: Jerry Lambert (jerry.lambert@tum.de) and Amedeo Ceruti (amedeo.ceruti@tum.de)

This file is used to optimize district heating systems with the tool topotherm of the Chair of
Energy Systems at Technical University of Munich.
"""

from topotherm.sts import sts_vertices
from topotherm.mts_easy_vertices import mts_easy_vertices
from topotherm.mts_vertices import mts_vertices


DISTRICTS = ["Example"]
DATAPATH = 'data/'
OUTPUTPATH = 'results/'

for run in DISTRICTS:
    print('============================')
    print(f'Running {run}')
    print('============================\n')
    sts_vertices(run, 1, DATAPATH, OUTPUTPATH, "topotherm_sts_vertices", 'sts')
    mts_easy_vertices(run, 1, DATAPATH, OUTPUTPATH, "topotherm_mts_easy_vertices", 'sts')
    mts_vertices(run, 1, DATAPATH, OUTPUTPATH, "topotherm_mts_vertices", 'sts')
