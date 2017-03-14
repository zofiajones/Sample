##Bash script to get bed file from list of genes.
##An easy way to get info for gene lists in different directories

rm gene.bed
while read p;do

id=$(python get_gene.py $p | grep -oP '\<Id\>\d+\</Id\>' | grep -oP '\d+')

echo -e $p"\t"$id >> id.txt
echo -e $p"\t"$id

python get_record.py $p $id


done < ~/mine_NCBI/test_genes.txt
#done < ~/custom_2/genes.txt
#done < ~/Manuela/genes.txt
#done < ~/Panel_beds/genes.txt
#done < not_found.txt
