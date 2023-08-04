inbam=$1
name=$2
gsize=$3
outdir=$4

macs2 callpeak -t $inbam -f BAMPE -n $name -g $gsize --broad --broad-cutoff 0.05 --keep-dup all --outdir $outdir

inbam=$1
name=$2
gsize=$3
outdir=$4

macs2 callpeak -t $inbam -f BAMPE -n $name -g $gsize --broad --broad-cutoff 0.05 --keep-dup all --outdir $outdir
