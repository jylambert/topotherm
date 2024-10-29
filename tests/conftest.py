"""Pytest configuration file."""


def pytest_addoption(parser):
    parser.addoption(
        "--solver",
        action="store",
        default="scip",
        help="Define the solver to use for the optimization. Default is cbc.")
