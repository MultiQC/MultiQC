#!/bin/sh
export MULTQC_IS_APP=True
source venv/bin/activate
"$PWD/venv/bin/python" scripts/multiqc "$@"