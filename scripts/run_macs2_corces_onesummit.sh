
#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables
inbam="${1}"
name="${2}"
gsize="${3}"
outdir="${4}"

main(){
    macs2_wrapper $inbam $name $gsize $outdir
}

macs2_wrapper(){
    local inbam="${1}"
    local name="${2}"
    local gsize="${3}"
    local outdir="${4}"
    #
    macs2 callpeak \
          --extsize 150 \
          --format BAMPE \
          --gsize $gsize \
          --keep-dup all \
          --name $name \
          --nolambda \
          --nomodel \
          --outdir $outdir \
          -p 0.01 \
          --shift -75 \
          --treatment $inbam
}

main "$@"
