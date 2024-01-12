
# - [[file:scripts/multiqc.sh][Shell script]]

#!/usr/bin/env bash

# Command line arguements
input="${1}"
out_name="${2}"
out_dir="${3}"

multiqc $input \
        --force \
        --outdir $out_dir \
        --filename $out_name


# - [[file:./scripts/multiqc.sh][Base script]]

multiqc_wrap()
    # Check for parameters, return usage if empty
    if [[ $# -eq 0 ]] || [[ multiqc_wrap == "h" ]] ; then
    printf "\n usage: multiqc_wrap input_dir output_dir output_prefix
           \n Wrapper for multiqc, see options in function
           \n $1 = input_dir
           \n $2 = output_dir
           \n $3 = output_dir_prefix
           \n "
    else
        multiqc $1 \
        --force \
        --dirs \
        --dirs-depth 1 \
        --outdir $2 \
        --filename atac_qc
    fi
}

# Snakemake variables
# Function
# Run command
#########1#########2#########3#########4#########5#########6#########7#########8
multiqc_wrap()
    # Check for parameters, return usage if empty
    if [[ $# -eq 0 ]] || [[ multiqc_wrap == "h" ]] ; then
    printf "\n usage: multiqc_wrap input_dir output_dir output_prefix
           \n Wrapper for multiqc, see options in function
           \n $1 = input_dir
           \n $2 = output_dir
           \n $3 = output_dir_prefix
           \n "
    else
        multiqc $1 \
        --force \
        --dirs \
        --dirs-depth 1 \
        --outdir $2 \
        --filename atac_qc
    fi
}
