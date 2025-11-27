Tests
=====

topotherm has additional dependencies for testing purposes. First install
development dependencies contained in pyproject.toml::

    pip install .[dev]

Then run the tests with pytest by running::

    pytest

All test scripts are contained in the ``./tests`` folder. If a new feature
is added, please also add tests that cover the new behavior.
