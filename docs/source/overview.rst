Overview
===========

Intro
------

Topotherm is a pyomo-based mixed-integer linear programming district heating
network design model which scales well into larger districts for single
and mulitple time steps.
It has been benchmarked against multiple other open-source models in a
publication available at:

> Lambert, Jerry and Ceruti, Amedeo and Spliethoff, Hartmut, Benchmark of Mixed-Integer Linear Programming Formulations for District Heating Network Design. Energy, Volume 308, 2024, 132885, ISSN 0360-5442, https://doi.org/10.1016/j.energy.2024.132885

Features
----------

* Single time-step topology and piping optimization.
* Multiple time-steps version additionally includes operation with variable heating demands.
* Greenfield optimization
* Existing network in development.
* Forced and economic expansion of the district heating network to consumers.
* Supports all solvers of pyomo, but the output helper functions in utils.py
   * might have to be rewritten (utils.py)
* Plotting functions included

Description
------------

To run the model, several incidence matrices have to be formulated. Then, the linear regression
parameters can be calculated for a given supply, ambient and return temperature of the network.
A pyomo model is then set up and solved with the solver of your choice.

Why should I use this?
-----------------------

Topotherm has the best scaling properties of multiple open-source models and
has been benchmarked and validated.

How to cite
-------------

> Lambert, Jerry and Ceruti, Amedeo and Spliethoff, Hartmut, Benchmark of Mixed-Integer Linear Programming Formulations for District Heating Network Design. Energy, Volume 308, 2024, 132885, ISSN 0360-5442, https://doi.org/10.1016/j.energy.2024.132885

Getting Started
----------------

This repository needs a PC capable to run python and its standard libraries.

Requirements
---------------

* Anaconda, mamba or venv




