#!/bin/bash
config_file=$1
pip install /local_arrakis/ArrakisND/
cd /local_scratch

# run the config file
mpirun -np 64 arrakis_nd "$config_file"