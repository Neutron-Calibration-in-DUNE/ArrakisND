#!/bin/bash
config_file=$1
cd /local_scratch

# run the config file
mpirun -np 64 arrakis_nd "$config_file"