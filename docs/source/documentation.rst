Docs
====

Activate the topotherm environment::

    conda activate topotherm
    # or: mamba activate topotherm
    # or: source .venv/bin/activate

Install documentation dependencies::

    pip install .[docs]

Build the HTML documentation::

    sphinx-build -b html docs/source docs/build/html

Then open ``docs/build/html/index.html``.
