Solver
========

Gurobi
---------

The results in the paper were obtained with the commercial solver gurobi.
A free academic license is available and can be installed by following
the documentation [here](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-).

Open-source Alternatives
--------------------------

You can try the code on smaller benchmarks with several open source solvers,
such as SCIP. Other popular open-source options are COIN-OR's cbc or HiGHS. ::

 conda
 conda activate topotherm
 conda install -c conda-forge pyscipopt
