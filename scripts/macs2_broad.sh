inbam=$1
name=$2
gsize=$3
outdir=$4

macs2 callpeak -t $inbam -f BAMPE -n $name -g $gsize --broad --broad-cutoff 0.05 --keep-dup all --outdir $outdir

#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

inbam=$1
name=$2
gsize=$3
outdir=$4

macs2 callpeak --treatment $inbam \
      --bdg \
      --broad \
      --broad-cutoff 0.05 \
      --format BAMPE \
      --gsize $gsize \
      --keep-dup all \
      --name $name \
      --outdir $outdir \
      --SPMR
