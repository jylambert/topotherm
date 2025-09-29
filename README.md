# topotherm

topotherm is a pyomo-based open-source optimization model for
district heating network design.

(intro)=
## Intro

Topotherm is a pyomo-based mixed-integer linear programming district heating
network design model which scales well into larger districts for single
and multiple time steps.
It has been benchmarked against multiple other open-source models in a
publication available at:

>  Lambert, Jerry and Ceruti, Amedeo and Spliethoff, Hartmut, Benchmark of Mixed-Integer Linear Programming Formulations for District Heating Network Design. Energy, Volume 308, 2024, 132885, ISSN 0360-5442, https://doi.org/10.1016/j.energy.2024.132885

(feature-overview)=
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

```{only} not html
## Contents

- topotherm
  - {ref}`Intro <intro>`
  - {ref}`Feature overview <feature-overview>`
  - {ref}`Description <description>`
  - {ref}`Why should I use this? <why-should-i-use-this>`
  - {ref}`How to cite <how-to-cite>`
  - {ref}`Getting Started <getting-started>`
    - {ref}`Requirements <requirements>`
  - {ref}`Install <install>`
    - {ref}`Python <python>`
    - {ref}`Anaconda or mamba <anaconda-or-mamba>`
  - {ref}`Solver <solver-readme>`
    - {ref}`Gurobi <gurobi>`
    - {ref}`Open-source Alternatives <open-source-alternatives>`
  - {ref}`Usage <usage>`
  - {ref}`Tests <tests>`
  - {ref}`Docs <docs>`
  

``` 

(description)=
## Description

To run the model, several incidence matrices have to be formulated. Then, the linear regression
parameters can be calculated for a given supply, ambient and return temperature of the network.
A pyomo model is then set up and solved with the solver of your choice.

(why-should-i-use-this)=
## Why should I use this?

Topotherm has the best scaling properties of multiple open-source models and
has been benchmarked and validated.

(how-to-cite)=
## How to cite

>  Lambert, Jerry and Ceruti, Amedeo and Spliethoff, Hartmut, Benchmark of Mixed-Integer Linear Programming Formulations for District Heating Network Design. Energy, Volume 308, 2024, 132885, ISSN 0360-5442, https://doi.org/10.1016/j.energy.2024.132885

(getting-started)=
## Getting Started

This repository needs a PC capable to run python and its standard libraries.

(requirements)=
### Requirements

* Anaconda, mamba or venv

(install)=
## Install

Use git to clone this repository into your computer. Then, install topotherm
with a package manager such as Anaconda, or directly with Python.

```bash
git clone https://github.com/jylambert/topotherm.git


(python)=
### Python

```python
cd topotherm
pip install .
```

This can also be done with venv or equivalent.

(anaconda-or-mamba)=
### Anaconda or mamba

We recommend to install the dependencies with anaconda or mamba:

```bash
cd topotherm
mamba env create -f environment.yml -n topotherm
mamba activate topotherm
```

(solver-readme)=
## Solver

(gurobi)=
### Gurobi

The results in the paper were obtained with the commercial solver gurobi.
A free academic license is available and can be installed by following
the documentation [here](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-).

(open-source-alternatives)=
### Open-source Alternatives

You can try the code on smaller benchmarks with several open source solvers,
such as SCIP. Other popular open-source options are COIN-OR's cbc or HiGHS.

```bash
mamba activate topotherm
mamba install -c conda-forge pyscipopt
```

(usage)=
## Usage

Generate the input incidence matrices for the district with .parquet format (see example).
Then, modify and run the either one of the three scripts in that folder.

```bash
cd examples
python run_sts.py
```

(contribute)=

## Contribute

Pull requests and any feedback regarding the code are very welcome. For major
changes, please open an issue first to discuss what you would like to change.


(tests)=
## Tests

To run the tests, use pytest.

```python
pip install .[dev]
pytest tests
```

(docs)=
## Docs

Activate the topotherm environment:

```bash
conda activate topotherm    # or: mamba activate topotherm / source .venv/bin/activate
```

Then, install the optional requirements:

```bash
pip install .[docs]
```

Build the HTML documentation from the repository root:

```bash
sphinx-build -b html docs/source docs/build/html
# then open: docs/build/html/index.html
```
(license)=
## License

[MIT](https://en.wikipedia.org/wiki/MIT_License), see LICENSE file.