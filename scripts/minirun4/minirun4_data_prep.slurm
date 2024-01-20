#!/bin/bash
#SBATCH -A dune                 # account to use for the job, '--account', '-A'
#SBATCH -J blip_graph_test      # job name, '--job-name', '-J'
#SBATCH -C gpu                  # type of job (constraint can be 'cpu' or 'gpu'), '--constraint', '-C'
#SBATCH -q shared               # Jobs requiring 1 or 2 gpus should use the shared setting, all others use 'regular'
#SBATCH -t 24:00:00             # amount of time requested for the job, '--time', 't'
#SBATCH -N 1                    # number of nodes, '--nodes', '-N'
#SBATCH -n 1                    # number of tasks '--ntasks', -n'
#SBATCH -c 32                   # number of cores per task, '--cpus-per-task', '-c'
#SBATCH --gpus-per-task=1       # number of gpus to be used per task
#SBATCH --gpus-per-node=1       # number of gpus per node.
#SBATCH --gpu-bind=none         # comment this out if you don't want all gpus visible to each task
#SBATCH --dependency=afterok:<minirun4_config_prep_id>      # don't run until prep job is finished
#SBATCH --array=0-<NUMBER_OF_ARRAKIS_PROCESSES>

export LOCAL_SCRATCH=/pscratch/sd/${USER:0:1}/${USER}/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID
mkdir -p $LOCAL_SCRATCH
LOCAL_BLIP=/global/cfs/cdirs/dune/users/${USER}
LOCAL_DATA=/global/cfs/cdirs/dune/users/${USER}

# find the arrakis_data associated with the array id

# use that arrakis_data as input to the minirun4_config.sh script
# along with the task_id

arrakis_data="${LOCAL_BLIP}/arrakis_data.csv"

# read in the arrakis file
if [ ! -e "$arrakis_data" ]; then
    echo "Error: Arrakis data file '$arrakis_data' not found!"
    exit 1
fi

arrakis_config=$(head -n $((SLURM_ARRAY_TASK_ID+1)) "$arrakis_data" | tail -1 | tr -d '\r')"/arrakis_config.yaml"

setfacl -m u:nobody:x /global/cfs/cdirs/dune/users/${USER}
shifter --image=docker:infophysics/blip:latest \
        --volume="${LOCAL_SCRATCH}:/local_scratch;${LOCAL_BLIP}:/local_blip;${LOCAL_DATA}:/local_data" \
        ./minirun4_data_prep.sh $arrakis_config