#!/bin/sh

#################################################################################################
# This BASH script explains how to set up a miniconda environment for running SCENIC and Scanpy
#################################################################################################

# Prerequisite: miniconda

# Create and set up new conda environment
conda create -n scenic_protocol python=3.6
conda activate scenic_protocol
conda install numpy pandas matplotlib seaborn
conda install -c anaconda xlrd
conda install -c anaconda cytoolz	
# Install scanpy (https://scanpy.readthedocs.io/en/latest/installation.html)
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph louvain
conda install -c conda-forge multicore-tsne
pip install scanpy
# Install pySCENIC
pip install pyscenic

# Install environment as kernel for Jupyter
pip install --user ipykernel
python -m ipykernel install --user --name=scenic_protocol

