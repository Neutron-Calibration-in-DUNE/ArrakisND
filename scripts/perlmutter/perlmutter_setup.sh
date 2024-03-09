#!/bin/bash
#SBATCH -A dune                 # account to use for the job, '--account', '-A'
#SBATCH -J example              # job name, '--job-name', '-J'
#SBATCH -C cpu                  # type of job (constraint can be 'cpu' or 'gpu'), '--constraint', '-C'
#SBATCH -q shared               # Jobs requiring 1 or 2 gpus should use the shared setting, all others use 'regular'
#SBATCH -t 8:00:00              # amount of time requested for the job, '--time', 't'
#SBATCH -N 2                    # number of nodes, '--nodes', '-N'
#SBATCH -n 128                  # number of tasks '--ntasks', -n'
#SBATCH -c 1                    # number of cores per task, '--cpus-per-task', '-c'

LOCAL_SCRATCH=/pscratch/sd/${USER:0:1}/${USER}
LOCAL_BLIP=/global/cfs/cdirs/dune/users/${USER}
LOCAL_DATA=/global/cfs/cdirs/dune/users/${USER}

setfacl -m u:nobody:x /global/cfs/cdirs/dune/users/${USER}
shifter --image=docker:infophysics/nersc:latest \
        --volume="${LOCAL_SCRATCH}:/local_scratch;${LOCAL_ARRAKIS}:/local_arrakis;${LOCAL_DATA}:/local_data" \
        bash