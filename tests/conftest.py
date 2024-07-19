def pytest_addoption(parser):
    parser.addoption("--solver", action="store", default="cbc",
                     help="Define the solver to use for the optimization. Default is cbc.")