
# - [[file:workflow/scripts/fastqc_wrapper.sh][Shell script]]

#!/usr/bin/env bash

# Script variables
input="${1}"
outdir="${2}"
threads="${3}"

# Functions
fastqc  --outdir $outdir \
        --quiet \
        --threads $threads $input
