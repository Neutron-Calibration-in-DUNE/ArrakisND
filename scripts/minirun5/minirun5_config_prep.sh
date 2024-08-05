#!/bin/bash
cd /local_scratch

prep_config="/local_arrakis/ArrakisND/config/minirun5.yaml"

# generate arrakis configs
create_arrakis_runs "$prep_config" -arrakis_config_location=${LOCAL_SCRATCH} -number_of_processes=${NUMBER_OF_ARRAKIS_PROCESSES}