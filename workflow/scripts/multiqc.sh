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
