string=$1
in_arr=($string)

for in_file in "${in_arr[@]}";
do
    out_file=$(echo $in_file | sed 's/atac_bams/qcpass/g')
    echo $in_file
    echo $out_file
    ln -frs $in_file $out_file
done
