#!/bin/bash

# RUN: sh scripts/run_display.sh
# OUTPUT: display app in your browser by following the generated link

pip install .

python3 - <<EOF
from arrakis_nd.utils.display.arrakis_display import ArrakisDisplay
ArrakisDisplay().run_app()
EOF