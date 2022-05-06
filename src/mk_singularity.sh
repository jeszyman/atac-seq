#!/usr/bin/env bash
#########1#########2#########3#########4#########5#########6#########7#########8
# Notes:
#  This script is tangled from source code in the cardradbio-atac.org file
#  Therefore any changes to this compiled file in place will likely be reverted
#  by a subsequent tangle operation.
#
#  This script assumes you are running from the repository main directory in
#  order to source the appriopriate configuation file. 

# Load bash configuation
if [ -f ./config/${HOSTNAME}.sh ]; then
    source ./config/${HOSTNAME}.sh
else
    echo "No bash config found. Are you running from repo main dir?"
fi

#
singularity pull $data_dir/atac.sif docker://jeszyman/atac
cp "${data_dir}/atac.sif" "${sif_dir}/atac.sif"
