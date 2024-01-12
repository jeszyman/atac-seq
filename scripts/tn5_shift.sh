
# - [[file:./scripts/todo_tn5_shift.sh][Base script]]

inbam=$1
outtmp=$2
outbam=$3
threads=$4

alignmentSieve --ATACshift --bam "$1" --numberOfProcessors $4 --outFile "$2"

samtools sort -@ $4 -o $3 $2
samtools index -@ $4 $3


# - [[file:./scripts/todo_tn5_shift.sh][Base script]]

alignmentSieve --ATACshift --bam "$1" --numberOfProcessors $2 --outFile "$3"

samtools sort -@ $2 -o $4 $3

samtools index -@ $2 $4


# - [[file:./scripts/todo_tn5_shift.sh][Base script]]

alignmentSieve --ATACshift --bam "$1" --numberOfProcessors $2 --outFile "$3"

samtools sort -@ $2 -o $4 $3

samtools index -@ $2 $4
