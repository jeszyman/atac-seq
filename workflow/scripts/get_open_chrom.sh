#########1#########2#########3#########4#########5#########6#########7#########8
alignmentSieve --bam $1 \
               --maxFragmentLength 150 \
               --numberOfProcessors $2 \
               --outFile $3 
samtools sort -@ $2 -o $4 $3
samtools index -@ $2 $4

#########1#########2#########3#########4#########5#########6#########7#########8
alignmentSieve --bam $1 \
               --maxFragmentLength 150 \
               --numberOfProcessors $2 \
               --outFile $3 
samtools sort -@ $2 -o $4 $3
samtools index -@ $2 $4
