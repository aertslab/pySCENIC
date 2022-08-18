FROM aertslab/pyscenic:0.12.0 AS compile-image

ENV DEBIAN_FRONTEND=noninteractive

# Install some dependencies:
#   - MulticoreTSNE: build-essential cmake
RUN apt-get update && \
    apt-get install -y --no-install-recommends build-essential cmake

# Install scanpy, MulticoreTSNE and ipykernel.
COPY requirements_docker_with_scanpy.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements_docker_with_scanpy.txt

FROM aertslab/pyscenic:0.12.0  AS build-image

RUN apt-get -y update && \
    apt-get -y --no-install-recommends install \
        # Need to run MulticoreTSNE
        libgomp1 && \
    rm -rf /var/cache/apt/* && \
    rm -rf /var/lib/apt/lists/*

COPY --from=compile-image /opt/venv /opt/venv
