
#!/usr/bin/env bash
#########1#########2#########3#########4#########5#########6#########7#########8
# Check for parameters, return usage if empty
if [[ "$#" -ne 3 ]];
then
    printf "\n usage: call_macs2_broad <BAM FILE> <OUTPUT BASENAME> <OUTPUT DIRECTORY>
    \n Wrapper function for calling broad beaks from ATAC-seq data with MACS2
    \n "
elif
    [[ ! -f "${1}.bai" ]]; then echo "No index for $1"
else
    macs2 callpeak \
          --bdg \
          --call-summits \
          --extsize 150 \
          --format BAMPE \
          --gsize mm \
          --keep-dup all \
          --name $2 \
          --nolambda \
          --nomodel \
          --outdir $3 \
          --SPMR \
          --treatment $1
fi
