
# - [[file:scripts/samstools_sats.sh][Base script]]

#!/usr/bin/env bash

in_bam="${1}"
out_stat="${2}"
out_flag="${3}"
threads="${4}"

samtools stats -@ $threads $in_bam > $out_stat
samtools flagstat -@ $threads $in_bam > $out_flag
