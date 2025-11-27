.. _overview:

Overview
========

Topotherm is a pyomo-based mixed-integer linear programming district
heating network design model which scales well into larger districts for
single and multiple time steps.

It has been benchmarked against multiple other open-source models in a
publication available at::

    Lambert, Jerry, Amedeo Ceruti, and Hartmut Spliethoff. "Benchmark of mixed-integer linear programming formulations for district heating network design." *Energy* 308 (2024): 132885. https://doi.org/10.1016/j.energy.2024.132885

Features
----------------

* Import/export functions from GIS files, geopandas.GeoDataFrames and pandas.DataFrames.
* Thermo-hydraulic precalculation models for district heating pipes.
* Single time-step topology and piping optimization.
    * Multiple time-steps version additionally includes operation with variable
        heating demands.
    * Greenfield optimization
    * Existing network in development.
* Forced and economic expansion of the district heating network to consumers.
* Supports all solvers of pyomo
* Plotting functions included
