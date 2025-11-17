Installation
============

Use git to clone this repository into your computer. Then, install topotherm
with a package manager such as Anaconda, or directly with Python.
Clone the repository::

    git clone https://github.com/jylambert/topotherm.git

----------------
Python Installer
----------------

Navigate into the repository and install Topotherm::

    cd topotherm
    pip install .

This can also be done inside ``venv`` or equivalent.

--------------------
Anaconda or Mamba
--------------------

We recommend installing dependencies using Anaconda or Mamba::

    cd topotherm
    mamba env create -f environment.yml -n topotherm
    mamba activate topotherm
