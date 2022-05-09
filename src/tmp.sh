repo=$1
mntpt=$2

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
