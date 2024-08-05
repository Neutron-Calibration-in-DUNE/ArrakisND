#!/bin/bash
# This script creates a conda environment called arrakis,
# and sets up ipykernel to use the arrakis environment
# in the jupyterlab interface.

module load python

# create conda environment

conda env create -f arrakis_environment.yaml

# set up ipykernel

conda activate arrakis
python -m ipykernel install --user --name env --display-name Arrakis