Documentation
=========

Activate the topotherm environment: ::

 conda activate topotherm

Then, install the optional requirments if needed ::

 pip install sphinx sphinx-rtd-theme or pip install .[docs]

In the topotherm root directory run: ::

 sphinx-apidoc -o docs/source/ .

Then make the hmtl documentation files: :: 

 make html

The documentation is located under docs/build/html

Adding New Pages
-----------------

Create a new rts file under docs/source

Add the name of the created page in index.rst