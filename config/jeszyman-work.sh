#!/usr/bin/env bash
#########1#########2#########3#########4#########5#########6#########7#########8

################################################################
###   Bash configuration sourcing script for jeszyman-work   ###
################################################################

# Host-local variables
data_dir=/mnt/ris/jschwarz/cardiac-radiobiology
mntpt=/mnt/ris/jschwarz
repo=/home/jeszyman/repos/cardradbio-atac
sif_dir=/home/jeszyman/sing_containers
threads=8

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
         
launch_atac() { 
    if [ -f /.dockerenv ]; then
        echo "shell already in docker, exiting";
        exit 1;
    else
        docker run --env HOME=${HOME} --hostname ${HOSTNAME} --interactive --tty --volume /home/:/home/ --volume /tmp/:/tmp/ --volume /mnt/:/mnt/ --user $(id -u ${USER}) -w "$repo" jeszyman/atac /bin/bash;
    fi
}
