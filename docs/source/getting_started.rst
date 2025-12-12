Getting Started
===============

Modeling Workflow
-----------------

topotherm has a generic modeling workflow. First, generate the input incidence
matrices for the district with .parquet format. These can be optionally generated
from gis files or pandas.DataFrames. Then, modify the boundary conditions and settings
to the local district ones. Finally, the model is created, solved and postprocessed.

Examples
---------

Several examples are provided in the `./examples` folder. To run one of the provided example scripts, run::

    cd examples
    python run_sts.py

There is also an example on how to calculate the thermo-hydraulic properties of pipes
with the defined settings at::

    ./examples/one_producer/pipe_calculations.py

Lastly, there are input incidence matrices for a small district heating network
provided at::

    ./examples/generating_inputs/
