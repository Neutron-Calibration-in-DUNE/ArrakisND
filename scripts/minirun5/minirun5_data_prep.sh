#!/bin/bash
config_file=$1
cd /local_scratch

pip install /local_arrakis/ArrakisND

# run the config file
mpirun -np 64 arrakis_nd "$config_file"