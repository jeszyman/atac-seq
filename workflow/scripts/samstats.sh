#########1#########2#########3#########4#########5#########6#########7#########8
samtools stats -@ $1 $2 > $3
samtools flagstat -@ $1 $2 > $4
