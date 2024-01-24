# topotherm

topotherm is a pyomo-based open-source optimization model for
district heating network design.

## Intro

Topotherm is a pyomo-based mixed-integer linear programming district heating 
network design model which scales well into larger districts for single
and mulitple time steps. 
It has been benchmarked against multiple other open-source models in this
pubblication:

> Insert citation here

## Feature overview

* Single time step topology and piping optimization.
  * Greenfield optimization only for now, existing network in development.
* Forced expansion of the distric theating network to consumers.
  * Economic expansion in development.
* Multiple time steps version additionally includes operation with variable
heating demands.
* Supports all solvers of pyomo, but the output might have to be rewritten (utils.py)
* Plotting functions included

## Contents

* [Description](#description)
* [How to cite](#how-to-cite)
* [Why should I use this?](#why-should-i-use-this)
* [Getting started](#getting-started)
  * [Requirements](#requirements)
  * [Install](#install)
  * [Usage](#usage)
* [Contribute](#contribute)
* [License](#license)

## Description

To run the model, several incidence matrices have to be formulated. Then, the linear regression
parameters can be calculated for a given supply, ambient and return temperature of the network. 
A pyomo model is then set up and solved with the solver of your choice.

## Why should I use this?

Topotherm has the best scaling of multiple open-source models and has
been benchmarked and validated.

## How to cite

> Insert citation here

## Getting Started

This repository needs a PC capable to run python and its standard libraries. 

### Requirements

* Anaconda, mamba or venv

### Install

Use git to clone this repository into your computer.

```git
git clone https://github.com/jylambert/topotherm.git
```

#### Anaconda

We recommend to install the dependencies with anaconda or mamba:

```conda
conda env create -f environment.yml
conda activate topotherm
```

#### pip

Alternatively, the packages can be installed manually via pip and/or in a venv with the `requirements.txt` file.

```Python
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```

### Solver

#### Gurobi

The results in the paper were obtained with the commercial solver gurobi.
A free academic license is available and can be installed by following
the documentation [here](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-).

#### Open-source Alternatives

You can try the code on smaller benchmarks with several open source solvers,
such as GLPK. Other popular open-source options are COIN-OR's cbc, HiGHS or SCIP.

```conda
conda activate topotherm
conda install -c conda-forge glpk
```

### Usage

Generate the input incidence matrices for the district with .parquet format (see Example).
Then, modify and run the main.py script.

```bash
python run_sts.py
```

## Contribute

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://en.wikipedia.org/wiki/MIT_License), see LICENSE file.

