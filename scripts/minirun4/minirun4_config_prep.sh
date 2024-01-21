#!/bin/bash
cd /local_scratch

prep_config="/workspace/ArrakisND/config/minirun4.yaml"

# generate arrakis configs
create_arrakis_runs "$prep_config" -arrakis_config_location=${LOCAL_BLIP} -number_of_processes=${NUMBER_OF_ARRAKIS_PROCESSES}