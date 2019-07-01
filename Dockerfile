FROM python:3.6.8-slim

ENV DEBIAN_FRONTEND=noninteractive
RUN BUILDPKGS="build-essential apt-utils git" && \
    apt-get update && \
    apt-get install -y $BUILDPKGS

# install dependencies:
COPY requirements.txt /tmp/
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r /tmp/requirements.txt && \
    pip install --no-cache-dir ipykernel 

# use version from argument (--build-arg version=0.9.7), or a default:
ARG version="0.9.9"
RUN pip install --no-cache-dir pyscenic==$version

RUN apt-get remove --purge -y $BUILDPKGS && \
    rm -rf /var/lib/apt/lists/*

