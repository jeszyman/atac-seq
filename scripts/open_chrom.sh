
input=$1
    tmp=$2
   open=$3
threads=$4

alignmentSieve --bam $input \
               --maxFragmentLength 150 \
               --numberOfProcessors $threads \
               --outFile $tmp
samtools sort -@ $threads -o $open $tmp
samtools index -@ $threads $open
