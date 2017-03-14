dir_ref[0]='/home/zjones/refs_ERR/'
dir_ref[1]='/home/zjones/refs_SRR/'
dir_ref[2]='/home/zjones/HG00/'
dir_ref[3]='/home/zjones/whole_batch/'
dir_ref[4]='/home/zjones/refs_ERR_77/'

for j in "${dir_ref[@]}";do

python getDist_zscore.py $j '_HG00_REF.zscore'
python ../inhouse/getDist_3_test.py $j '_HG00_REF.zscore'

python getDist_zscore.py $j '.zscore'
python getDist_zscore.py $j '_referr_REF.zscore'
python getDist_zscore.py $j '_refsrr_REF.zscore'
python getDist_zscore.py $j '_referr_77_REF.zscore'
python getDist_zscore.py $j '_53_REF.zscore'
python getDist_zscore.py $j '_whole_batch_REF.zscore'

python ../inhouse/getDist_3_test.py $j '_53_REF.zscore'
python ../inhouse/getDist_3_test.py $j '_whole_batch_REF.zscore'
python ../inhouse/getDist_3_test.py $j '_referr_REF.zscore'
python ../inhouse/getDist_3_test.py $j '_refsrr_REF.zscore'
python ../inhouse/getDist_3_test.py $j '_referr_77_REF.zscore'
python ../inhouse/getDist_3_test.py $j '.zscore'

python getDist_3c.py $j '_ruh.pickle'

done

for j in "${dir_ref[@]}";do

python ~/NIPT_pipeline/samples/Horizon_bwa_mem/getDist_3b_test.py $j '_HG00_REF.gcc'
python ~/NIPT_pipeline/samples/Horizon_bwa_mem/getDist_3b_test.py $j '_referr_77_REF.gcc'
python ~/NIPT_pipeline/samples/Horizon_bwa_mem/getDist_3b_test.py $j '_refsrr_REF.gcc'

done
