# topotherm

topotherm is a pyomo-based open-source optimization model for
district heating network design.

## Contents

- [topotherm](#topotherm)
  - [Contents](#contents)
  - [Overview](#overview)
    - [Features](#features)
    - [Why should I use this?](#why-should-i-use-this)
    - [How to cite](#how-to-cite)
  - [Install](#install)
    - [Requirements](#requirements)
    - [Get the Code](#get-the-code)
    - [Install the required Python packages](#install-the-required-python-packages)
    - [Anaconda](#anaconda)
    - [Solver](#solver)
      - [Gurobi](#gurobi)
      - [Open-source Alternatives](#open-source-alternatives)
  - [Usage](#usage)
  - [Documentation](#documentation)
  - [Contribute](#contribute)
  - [Tests](#tests)
  - [License](#license)

## Overview

topotherm is a pyomo-based mixed-integer linear programming district heating
network design model which scales well into larger districts for single
and multiple time steps.
It has been benchmarked against multiple other open-source models in a
publication available at:

>  Lambert, Jerry and Ceruti, Amedeo and Spliethoff, Hartmut, Benchmark of Mixed-Integer Linear Programming Formulations for District Heating Network Design. Energy, Volume 308, 2024, 132885, ISSN 0360-5442, https://doi.org/10.1016/j.energy.2024.132885

### Features

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

### Why should I use this?

Topotherm has the best scaling properties of multiple open-source models and
has been benchmarked and validated.

### How to cite

If you use topotherm, please cite:

>  Lambert, Jerry and Ceruti, Amedeo and Spliethoff, Hartmut, Benchmark of Mixed-Integer Linear Programming Formulations for District Heating Network Design. Energy, Volume 308, 2024, 132885, ISSN 0360-5442, https://doi.org/10.1016/j.energy.2024.132885

## Install

### Requirements

topotherm requires Python 3.10 or higher. We recommend using Anaconda or Mamba to
manage the Python environment and dependencies. The required packages are listed
in the ``environment.yml`` file in the repository and the ``requirements.txt`` file
for pip installations.

### Get the Code

Use git to clone this repository into your computer. Then, install topotherm
with a package manager such as Anaconda, or directly with Python.

```bash
git clone https://github.com/jylambert/topotherm.git
```

### Install the required Python packages

```python
cd topotherm
pip install .
```

We strongly recommend to do this with a package manager, Anaconda, mamba, venv or equivalent.

### Anaconda

To install the dependencies with anaconda or mamba:

```bash
cd topotherm
mamba env create -f environment.yml -n topotherm
mamba activate topotherm
```

### Solver

#### Gurobi

The results in the paper were obtained with the commercial solver gurobi.
A free academic license is available and can be installed by following
the documentation [here](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-).

#### Open-source Alternatives

You can also run the code with several open source solvers, such as SCIP.
Other popular open-source options are COIN-OR's cbc or HiGHS.

```bash
mamba activate topotherm
mamba install -c conda-forge pyscipopt
```

## Usage

In general, the users need to generate the required input incidence matrices for the district they want to run.
Some examples are provided in the examples folder. These can be a starting point for
custom districts.

```bash
cd examples/one_producer
python run_sts.py
```

## Documentation

topotherm has a readthedocs page available at <https://topotherm.readthedocs.io/en/latest/>.
However, the documentation can be built locally with sphinx locally by following the
steps below.

Install the optional requirements:

```bash
conda activate topotherm    # or: mamba activate topotherm / source .venv/bin/activate
pip install .[docs]
```

Build the HTML documentation from the repository root:

```bash
sphinx-build -b html docs/source docs/build/html
# then open: docs/build/html/index.html
```

## Contribute

Pull requests and any feedback regarding the code are very welcome. For major
changes, please open an issue first to discuss what you would like to change.

## Tests

Testing is done with pytest. Please pass all tests before opening a new pull request.

```python
pip install .[dev]
pytest
```

## License

[MIT](https://en.wikipedia.org/wiki/MIT_License), see LICENSE file.
