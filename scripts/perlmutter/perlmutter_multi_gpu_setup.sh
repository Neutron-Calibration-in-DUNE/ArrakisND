#!/bin/bash 
#SBATCH -C gpu
#SBATCH -A dune_g
#SBATCH --ntasks-per-node 4
#SBATCH --cpus-per-task 32
#SBATCH --gpus-per-node 4
#SBATCH --time=01:00:00
#SBATCH --image=docker:infophysics/nersc:latest
#SBATCH --module=gpu,nccl-2.18
#SBATCH --reservation=dl_test
#SBATCH -J vit-era5-mp
#SBATCH -o %x-%j.out

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

# if cuda graphs, use train_mp_graphs.py
set -x

LOCAL_SCRATCH=/pscratch/sd/${USER:0:1}/${USER}
LOCAL_ARRAKIS=/global/cfs/cdirs/dune/users/${USER}
LOCAL_DATA=/global/cfs/cdirs/dune/users/${USER}

setfacl -m u:nobody:x /global/cfs/cdirs/dune/users/${USER}
shifter --image=docker:infophysics/nersc:latest \
        --volume="${LOCAL_SCRATCH}:/local_scratch;${LOCAL_ARRAKIS}:/local_arrakis;${LOCAL_DATA}:/local_data" \
        bash