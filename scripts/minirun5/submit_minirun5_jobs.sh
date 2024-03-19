# Submit the preparation job and get its job ID

export NUMBER_OF_ARRAKIS_PROCESSES=25

minirun5_config_prep_jobid=$(sbatch minirun5_config_prep.slurm | cut -d ' ' -f 4)

# Submit the array job with a dependency on the preparation job
sbatch --dependency=afterok:$minirun5_config_prep_jobid \
       --array=0-$NUMBER_OF_ARRAKIS_PROCESSES minirun5_data_prep.slurm \
       --output=/pscratch/sd/${USER:0:1}/${USER}/arrakis_output_%A_%a.out \
       --error=/pscratch/sd/${USER:0:1}/${USER}/arrakis_error_%A_%a.err