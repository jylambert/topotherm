Install
=========

Use git to clone this repository into your computer. Then, install topotherm
with a package manager such as Anaconda, or directly with Python. ::

   git clone https://github.com/jylambert/topotherm.git

Python
---------

This can also be done with venv or equivalent. :: 
 
   cd topotherm
   python setup.py install

Anaconda or mamba
---------------------

We recommend to install the dependencies with anaconda or mamba: ::

 cd topotherm
 conda create -n topotherm python
 conda activate topotherm
 pip install -e .

Solver
---------

The results in the paper were obtained with the commercial solver gurobi.
A free academic license is available and can be installed by following
the documentation https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-.

You can try the code on smaller networks with several open source solvers,
such as SCIP. Other popular open-source options are COIN-OR's cbc or HiGHS. ::

 conda activate topotherm
 conda install -c conda-forge pyscipopt

