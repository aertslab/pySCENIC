FROM python:3.6.8-slim

RUN BUILDPKGS="build-essential apt-utils" && \
    apt-get update && \
    apt-get install -y $BUILDPKGS && \
    apt-get install -y procps && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --upgrade pyscenic==0.9.6 dask==1.0.0 pandas==0.23.4

RUN apt-get remove --purge -y $BUILDPKGS

