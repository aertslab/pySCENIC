BootStrap: docker
From: continuumio/miniconda3

%post
    . /opt/conda/etc/profile.d/conda.sh
    apt-get update && apt-get install -y build-essential
    conda install python=3.6
    conda activate
    pip install --no-cache-dir --upgrade pyscenic

