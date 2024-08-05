# Submit the preparation job and get its job ID

export NUMBER_OF_ARRAKIS_PROCESSES=10

minirun4_config_prep_jobid=$(sbatch minirun4_config_prep.slurm | cut -d ' ' -f 4)

# Submit the array job with a dependency on the preparation job
sbatch --dependency=afterok:$minirun4_config_prep_jobid \
       --array=0-$NUMBER_OF_ARRAKIS_PROCESSES minirun4_data_prep.slurm \
       --output=/pscratch/sd/${USER:0:1}/${USER}/arrakis_output_%A_%a.out \
       --error=/pscratch/sd/${USER:0:1}/${USER}/arrakis_error_%A_%a.err