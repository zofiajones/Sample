##FUNCTION: To cross-referece the presence of certain variants in NGS, dbSNP and (previously processed) customer data

#use input file to get targets for tabix
grep -v '#' out_variants_3_head.vcf | awk 'BEGIN{OFS=""}{print $1,":",$2-10,"-",$2+10}' > targets.txt

#parse input file to one base per line to fine tune matching and get locatons in bed file format
cat out_variants_3_head.vcf | python get_mnp_loc.py > out_variants_4_head.vcf
grep -v '#' out_variants_4_head.vcf | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3}'  > loc.txt
#deal with asterixes
sed 's/*/\*/g' loc.txt > loc_1.txt

#grab info from vcf.gz
rm *hotspot_req.vcf
for i in $(ls *hotspot.vcf.gz);do

root=$(echo ${i%hotspot.vcf.gz})
echo $j

while read p;do

echo -e $root"\thotspot"
~/hb/tools/tabix $i $p >> $root"hotspot_req.vcf"


done < targets.txt


for i in $(ls *_hotspot_req.vcf);do

awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$4,$5,$8,$9,$10}' $i > ${i%.vcf}.bed
~/hb/tools/bedtools intersect -a request.bed -b ${i%.vcf}.bed  -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$8,$9,$10,$11,$12}' | sort | uniq | python parse_multi_ion_bed_1.py  > "request_"${i%.vcf}".txt"

done

rm tru.txt
#readin loc_1.txt
while read p;do

array=($p)

#do a dirty regex
echo $p | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | sed 's/*/\\*/g' > match.txt
matchGTB=$(grep -f match.txt request_GTB_348.txt | awk '{print $7}'  | sort | uniq | awk -F ',' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5}')

#if no variant info, get coverage from bam
if [[ -z $matchGTB  ]];then
loc=$(echo -e ${array[0]}":"${array[1]}"-"${array[2]})
DP=$(~/hb/tools/samtools view /hd1/HD701_bams/GTB_348_IonXpress_028_.bam  $loc | wc -l)
matchGTB="NA NA NA NA "$DP
fi

#get loc for further checks
matchGTB1=$(grep -f match.txt request_GTB_348.txt | awk '{print $5,$6}'  | sort | uniq | awk -F ',' 'BEGIN{OFS="@"}{print $1,$2}')
echo -e ${array[0]}"@"${array[2]}"@"$matchGTB1 | sed 's/\./\\./g' > match_base.txt


#check for rs_id (dbSNP)
~/hb/tools/tabix ~/clinvar_20160831.vcf.gz $loc | python get_mnp_loc_chr.py | awk 'BEGIN{OFS=""}{print "chr",$1,"@",$2,"@",$4,"@",$5,"@",$3}'  > rs_base.txt
rs_id=$(grep -f match_base.txt rs_base.txt | sort | uniq | awk -F '@' '{printf "%s,", $5}' )

#get list of found variant frequencies in customer data
awk 'BEGIN{OFS="@"}{print $1,$3,$4,$5}' /hd1/HD701_set/request/all_parsed.vcf > cust_match.txt
cust_info=$(grep -f match_base.txt cust_match.txt | awk -F '@' '{printf "%s", $5}' )

#again check for empty variables
if [[ -z $rs_id  ]];then
rs_id="NA"
fi

if [[ -z $cust_info  ]];then
cust_info="NA"
fi

#add to output file
echo -e $matchGTB"\t"$rs_id"\t"$cust_info >> tru.txt

done < loc_1.txt
