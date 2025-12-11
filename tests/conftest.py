"""Pytest configuration file."""


def pytest_addoption(parser):
    parser.addoption(
        "--solver",
        action="store",
        default="highs",
        help="Define the solver to use for the optimization. Default is highs.",
    )
