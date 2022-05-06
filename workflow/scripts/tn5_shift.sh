alignmentSieve --ATACshift --bam "$1" --numberOfProcessors $2 --outFile "$3"

samtools sort -@ $2 -o $4 $3

samtools index -@ $2 $4

alignmentSieve --ATACshift --bam "$1" --numberOfProcessors $2 --outFile "$3"

samtools sort -@ $2 -o $4 $3

samtools index -@ $2 $4

alignmentSieve --ATACshift --bam "$1" --numberOfProcessors $2 --outFile "$3"

samtools sort -@ $2 -o $4 $3

samtools index -@ $2 $4

alignmentSieve --ATACshift --bam "$1" --numberOfProcessors $2 --outFile "$3"

samtools sort -@ $2 -o $4 $3

samtools index -@ $2 $4
