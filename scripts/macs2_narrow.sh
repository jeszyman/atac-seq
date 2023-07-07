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
      --call-summits \
      --extsize 150 \
      --format BAMPE \
      --gsize $gsize \
      --keep-dup all \
      --name ${name}_multi \
      --nolambda \
      --outdir $outdir \
      -p 0.01 \
      --shift -75 \
      --SPMR \
      --nomodel

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
      --call-summits \
      --extsize 150 \
      --format BAMPE \
      --gsize $gsize \
      --keep-dup all \
      --name ${name}_multi \
      --nolambda \
      --outdir $outdir \
      -p 0.01 \
      --shift -75 \
      --SPMR \
      --nomodel
