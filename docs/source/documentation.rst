.. _documentation:

Documentation
==============

topotherm has a readthedocs page available at https://topotherm.readthedocs.io/en/latest/.
However, you can also build the documentation with sphinx locally by following the
steps below.

Building the Documentation
------------------------------------------

Activate the topotherm environment in your package manager::

    conda activate topotherm
    # or: mamba activate topotherm
    # or: source .venv/bin/activate

Install the documentation dependencies listed in `pyproject.toml`::

    pip install .[docs]

Then, build the HTML documentation with sphinx::

    sphinx-build -b html docs/source docs/build/html

Now you can open ``docs/build/html/index.html`` and navigate the
compiled documentation in html form.
