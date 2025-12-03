Tests
=====

topotherm has additional dependencies for testing purposes.
To run them, first install development dependencies contained in pyproject.toml::

    pip install .[dev]

All test scripts are contained in the ``./tests`` folder. If a new feature
is added, please also add tests that cover the new behavior.
Then run all supported tests in `./tests` with pytest by running::

    pytest
