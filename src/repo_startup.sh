#!/usr/bin/env bash
repo=$1
mntpt=$2
sif_dir=$3

# Check for parameters, return usage if empty
if [ $# -ne 3 ];
then
    printf "\n usage: repo_startup.sh <REPO PATH> <RIS MOUNT PT> <SINGULARITY CONTAINER DIR>
    \n ATAC-seq repo development helper script 
    \n "
else

    # Check git file hook is read-able
    if [ -r "${repo}/.git/hooks/precommit" ]; then
        echo "Git size check is read-able"
    else
        echo
        "Git size check is not read-able"
        exit 1
    fi
    
    # Check mount point  
    if grep -qs $mntpt /proc/mounts; then
        echo "RIS storage mounted."
    else
        echo "RIS storage NOT mounted, exiting."
        exit 1
    fi

    # Check singularity container
    if [ -r $sif_dir/atac.sif ]; then
        echo "Local SIF file present"
    else
        echo "No local SIF file found"
        exit 1
    fi

    # Check singularity container up-to-date
    if [ /mnt/ris/jschwarz/cardiac-radiobiology/atac.sif -nt $sif_dir/atac.sif ]; then
        echo "Local SIF is out of date. Updating ..."
        cp /mnt/ris/jschwarz/cardiac-radiobiology/atac.sif $sif_dir/atac.sif 
    else
        echo "Local SIF file is up to date"
    fi
fi
