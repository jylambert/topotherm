Install
=========

Use git to clone this repository into your computer. Then, install topotherm
with a package manager such as Anaconda, or directly with Python. ::

   git
   git clone https://github.com/jylambert/topotherm.git

Python
---------

This can also be done with venv or equivalent. :: 
 
   Python
   cd topotherm
   python setup.py install

Anaconda or mamba
---------------------

We recommend to install the dependencies with anaconda or mamba: ::

   mamba
   cd topotherm
   mamba env create -f environment.yml -n topotherm
   mamba activate topotherm
