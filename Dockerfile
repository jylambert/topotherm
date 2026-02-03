FROM python:3.13-slim

# define internal container file dir
WORKDIR /workspace

# copy pyproject to the WORKDIR
COPY ./pyproject.toml ./
# required for first build, will fail otherwise due to pyproj definition
COPY ./topotherm ./topotherm
# Install dependecies from the pyproj file
# just to make sure pip is up to date
RUN pip install --upgrade pip
RUN pip install --no-cache-dir .[dev,docs]

# this will be mounted in the docker compose
# COPY . .
# CMD ["/bin/bash"]
