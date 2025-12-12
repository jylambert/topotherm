Installation
============

--------------------
Requirements
--------------------

topotherm requires Python 3.10 or higher. We recommend using Anaconda or Mamba to
manage the Python environment and dependencies. The required packages are listed
in the ``environment.yml`` file in the repository and the 
``requirements.txt`` file for pip installations.

----------------
Get the Code
----------------

Use git to clone topotherm to your computer or download and unzip the files.

Clone the repository::

    git clone https://github.com/jylambert/topotherm.git

----------------------
Install with Python
----------------------

Navigate into the repository and install topotherm::

    cd topotherm
    pip install .

This can also be done inside ``venv`` or equivalent.

-------------------------------
Install with Anaconda or Mamba
-------------------------------

We recommend installing dependencies using Anaconda or Mamba::

    cd topotherm
    mamba env create -f environment.yml -n topotherm
    mamba activate topotherm

------
Solver
------

topotherm is a MILP optimization model and requires a Mixed-Integer Linear Programming solver.
There are both commercial and open-source solvers available. In our experience,
commercial solvers such as gurobi are able to solve larger instances, while
open-source solvers suffice for smaller districts.

Gurobi
------

The benchmarking results were obtained with the commercial solver gurobi.
A free academic license is available at https://www.gurobi.com/.
After installing gurobi, the python interface for the corresponding version
needs to be installed, following the documentation:
https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-.

Open Source Alternatives
------------------------

You can try the code on smaller districts with several open source solvers,
such as SCIP. Other popular open-source options are COIN-OR's cbc or HiGHS.::

    mamba activate topotherm
    mamba install -c conda-forge pyscipopt
