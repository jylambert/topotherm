# topotherm

topotherm is a pyomo-based open-source optimization model for
district heating network design.

## Intro

Topotherm is a pyomo-based mixed-integer linear programming district heating
network design model which scales well into larger districts for single
and mulitple time steps.
It has been benchmarked against multiple other open-source models in a
publication available at:

>  Lambert, Jerry and Ceruti, Amedeo and Spliethoff, Hartmut, Benchmark of Milp Formulations for District Heating Network Design. Available at SSRN: https://ssrn.com/abstract=4688228 or http://dx.doi.org/10.2139/ssrn.4688228 

## Feature overview

* Single time-step topology and piping optimization.
  * Multiple time-steps version additionally includes operation with variable
  heating demands.
  * Greenfield optimization
  * Existing network in development.
* Forced and economic expansion of the district heating network to consumers.
* Supports all solvers of pyomo, but the output helper functions in utils.py
might have to be rewritten (utils.py)
* Plotting functions included

## Contents

- [topotherm](#topotherm)
  - [Intro](#intro)
  - [Feature overview](#feature-overview)
  - [Contents](#contents)
  - [Description](#description)
  - [Why should I use this?](#why-should-i-use-this)
  - [How to cite](#how-to-cite)
  - [Getting Started](#getting-started)
    - [Requirements](#requirements)
  - [Install](#install)
    - [Python](#python)
    - [Anaconda or mamba](#anaconda-or-mamba)
  - [Solver](#solver)
    - [Gurobi](#gurobi)
    - [Open-source Alternatives](#open-source-alternatives)
  - [Usage](#usage)
  - [Contribute](#contribute)
  - [Tests](#tests)
  - [License](#license)

## Description

To run the model, several incidence matrices have to be formulated. Then, the linear regression
parameters can be calculated for a given supply, ambient and return temperature of the network.
A pyomo model is then set up and solved with the solver of your choice.

## Why should I use this?

Topotherm has the best scaling properties of multiple open-source models and
has been benchmarked and validated.

## How to cite

>  Lambert, Jerry and Ceruti, Amedeo and Spliethoff, Hartmut, Benchmark of Mixed-Integer Linear Programming Formulations for District Heating Network Design. Energy, Volume 308, 2024, 132885, ISSN 0360-5442, https://doi.org/10.1016/j.energy.2024.132885

## Getting Started

This repository needs a PC capable to run python and its standard libraries.

### Requirements

* Anaconda, mamba or venv

## Install

Use git to clone this repository into your computer. Then, install topotherm
with a package manager such as Anaconda, or directly with Python.

```git
git clone https://github.com/jylambert/topotherm.git
```

### Python

```Python
cd topotherm
python setup.py install
```

This can also be done with venv or equivalent.

### Anaconda or mamba

We recommend to install the dependencies with anaconda or mamba:

```conda
cd topotherm
conda create -n topotherm python
conda activate topotherm
pip install -e .
```

## Solver

### Gurobi

The results in the paper were obtained with the commercial solver gurobi.
A free academic license is available and can be installed by following
the documentation [here](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-).

### Open-source Alternatives

You can try the code on smaller benchmarks with several open source solvers,
such as SCIP. Other popular open-source options are COIN-OR's cbc or HiGHS.

```conda
conda activate topotherm
conda install -c conda-forge pyscipopt
```

## Usage

Generate the input incidence matrices for the district with .parquet format (see example).
Then, modify and run the either one of the three scripts in that folder.

```bash
cd example
python run_sts.py
```

## Contribute

Pull requests and any feedback regarding the code are very welcome. For major
changes, please open an issue first to discuss what you would like to change.

## Tests

To run the tests, use pytest.

```Python
pytest tests
```

## License

[MIT](https://en.wikipedia.org/wiki/MIT_License), see LICENSE file.
