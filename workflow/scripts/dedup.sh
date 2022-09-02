raw_bam="${1}"
threads="${2}"
dedup_bam="${3}"
samtools sort -@ $threads -n -o - $raw_bam |
    samtools fixmate -m - - |
    samtools sort -@ $threads -o - - |
    samtools markdup -@ $threads -r - $dedup_bam
samtools index $dedup_bam
