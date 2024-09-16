"""Set up topotherm from requirements."""

import os
from setuptools import setup, find_packages


# The directory containing this file
here = os.path.dirname(os.path.realpath(__file__))

# The text of the README file
with open(os.path.join(here, "README.md"), encoding='utf-8') as fid:
    readme = fid.read()

with open(os.path.join(here, "requirements.txt"), encoding='utf-8') as fid:
    reqs = fid.read().splitlines()

# This call to setup() does all the work
setup(
    name="topotherm",
    version="0.1.0",
    description="A package to design district heating networks",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/jylambert/topotherm",
    author="Jerry Lambert, Amedeo Ceruti",
    author_email="jerry.lambert@tum.de, amedeo.ceruti@tum.de",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: District Heating",
    ],
    packages=find_packages(exclude=("tests", "docs")),
    include_package_data=True,
    install_requires=reqs,
    keywords=["district heating", "optimization", "Linear Programming",
              "Mixed Integer Linear Programming"]
)
