#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables
inbam=$1
outdir=$2

main(){
    macs2_wrapper \
        $inbam \
        $outdir
}

# Functions
macs2_wrapper(){
    local inbam="${1}"
    local outdir="${2}"
    #
    base=$(basename -s .bam $inbam)
    #
    macs2 callpeak \
           --extsize 150 \
          --format BAMPE \
          --gsize mm \
          --keep-dup all \
          --name $base \
          --nolambda \
          --nomodel \
          --outdir $outdir \
          -p 0.01 \
          --shift -75 \
          --treatment $inbam
}

# Run
main "$@"
