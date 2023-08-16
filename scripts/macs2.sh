inbam=$1
name=$2
gsize=$3
outdir=$4

macs2 callpeak --treatment $inbam \
      --format BAMPE \
      --name $name \
      --gsize $gsize \
      --broad \
      --broad-cutoff 0.05 \
      --keep-dup all \
      --outdir $outdir

macs2 callpeak --treatment $inbam \
      --bdg \
      --call-summits \
      --extsize 150 \
      --format BAMPE \
      --gsize $gsize \
      --keep-dup all \
      --name ${name} \
      --nolambda \
      --outdir $outdir \
      -p 0.01 \
      --shift -75 \
      --SPMR \
      --nomodel
