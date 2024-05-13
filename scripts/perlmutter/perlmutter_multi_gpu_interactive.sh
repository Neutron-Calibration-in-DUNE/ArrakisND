#!/bin/bash

#SBATCH -N 4
#SBATCH -C gpu
#SBATCH -G 16
#SBATCH -q interactive
#SBATCH -J gpu_test
#SBATCH --mail-user=ncarrara.physics@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -A dune_g
#SBATCH -t 2:0:0

# Resource allocation
salloc --constraint gpu --account dune_g --qos interactive --ntasks-per-node=4 \
       --cpus-per-task=32 --gpus-per-node=4 --time=02:00:00 \
       --image=docker:infophysics/nersc:latest 
       #--module=gpu,nccl-2.18 --reservation=dl_test -J vit-era5-mp

# Environmental settings
export FI_MR_CACHE_MONITOR=userfaultfd
export HDF5_USE_FILE_LOCKING=FALSE

# Profiling
if [ "${ENABLE_PROFILING:-0}" -eq 1 ]; then
    echo "Enabling profiling..."
    NSYS_ARGS="--trace=cuda,cublas,nvtx --cuda-graph-trace=node --kill none -c cudaProfilerApi -f true"
    NSYS_OUTPUT=${LOGDIR}/${PROFILE_OUTPUT:-"profile"}
    export PROFILE_CMD="nsys profile $NSYS_ARGS -o $NSYS_OUTPUT"
fi

export MASTER_ADDR=$(hostname)

# Reversing order of GPUs to match default CPU affinities from Slurm
export CUDA_VISIBLE_DEVICES=3,2,1,0

# Setup scratch and data directories
export LOCAL_SCRATCH=/pscratch/sd/${USER:0:1}/${USER}
export LOCAL_ARRAKIS=/global/cfs/cdirs/dune/users/${USER}
export LOCAL_DATA=/global/cfs/cdirs/dune/users/${USER}

setfacl -m u:nobody:x /global/cfs/cdirs/dune/users/${USER}

# After the salloc command completes and resources are allocated, you can manually run commands or scripts interactively
echo "Resources allocated. You can now run your commands interactively."

# To use shifter interactively after resource allocation
shifter --image=docker:infophysics/nersc:latest \
        --volume="${LOCAL_SCRATCH}:/local_scratch;${LOCAL_ARRAKIS}:/local_arrakis;${LOCAL_DATA}:/local_data" \
        bash