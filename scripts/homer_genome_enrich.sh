bed="${1}"
fasta="${2}"
outdir="${3}"
threads="${4}"

findMotifsGenome.pl $bed $fasta $outdir -p $threads
