FROM python:3.10.6-slim-bullseye AS compile-image

ENV DEBIAN_FRONTEND=noninteractive

# Install some dependencies:
#   - MulticoreTSNE: build-essential cmake
#   - git: to set the pySCENIC version
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential cmake git

# Create virtual environment.
RUN python -m venv /opt/venv

# Make sure we use the virtual environment.
ENV PATH="/opt/venv/bin:$PATH"

# Install pySCENIC dependencies with pip.
COPY requirements_docker.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Install scanpy, MulticoreTSNE and ipykernel.
#COPY requirements_docker_with_scanpy.txt /tmp/
#RUN pip install --no-cache-dir -r /tmp/requirements_docker_with_scanpy.txt

# Install pySCENIC from local copy:
COPY . /tmp/pySCENIC
RUN  cd /tmp/pySCENIC && \
     pip install . && \
     cd .. && rm -rf pySCENIC

FROM python:3.10.6-slim-bullseye AS build-image

RUN apt-get -y update && \
    apt-get -y --no-install-recommends install \
        # Need to run ps
        procps \
        libxml2 \
        less \
        # Need to run MulticoreTSNE
        libgomp1 && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

COPY --from=compile-image /opt/venv /opt/venv

# Make sure we use the virtualenv:
ENV PATH="/opt/venv/bin:$PATH"

EXPOSE 8787

