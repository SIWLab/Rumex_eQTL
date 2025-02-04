#!/bin/bash

# for each iteration 1:1000
# do 5 as a start
cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/ml

for i in {1..50}; do
# reshuffle phenotype in R
/usr/local/bin/Rscript --vanilla pheno_permute.R ${i}

# run tensorqtl
bgzip normalized_counts_mln_auto_${i}.bed && tabix -p bed normalized_counts_mln_auto_${i}.bed.gz
done

for i in {1..50}; do
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/ml/normalized_counts_mln_auto_${i}.bed.gz
prefix=ML_${i}
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1
done

for i in {1..50}; do
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/ml/normalized_counts_mln_auto_${i}.bed.gz
prefix=ML_${i}
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
done

# run MAF_null scripts in R to get ML_i_eqtl_fake.txt
#/usr/local/bin/Rscript --vanilla MAF_null.R ${i}

for i in {1..50}; do
/usr/local/bin/Rscript --vanilla MAF_null_new.R ${i}
done
cat ML_*_eqtl_fake.txt > ML_MAF_null_new.txt




for i in {1..50}; do
/usr/local/bin/Rscript --vanilla MAF_null_random_new.R ${i}
done

cat ML_*_eqtl_min_fake.txt > ML_MAF_null_min.txt
cat ML_*_eqtl_random_fake.txt > ML_MAF_null_random.txt

# new
cat ML_*_eqtl_min_fake.txt > ML_MAF_null_min_q01_filt.txt
cat ML_*_eqtl_random_fake.txt > ML_MAF_null_random_q01_filt.txt



# remove intermediate files
#rm normalized_counts_mln_auto_${i}.bed.gz*
# rm ML_${i}.cis_qtl_pairs.*.parquet
#rm ML_${i}.tensorQTL.cis_nominal.log

#mv ML_${i}.cnt cnt/
#mv ML_${i}_eqtl_fake.txt fake/
#done

###########################################################
# after the loop
cat ./fake/ML_*_eqtl_fake.txt > ML_MAF_null.txt
cat ./cnt/ML_*.cnt > ML_fake_cnt.txt



cat ./fake/MP_*_eqtl_fake.txt > MP_MAF_null.txt
cat ./cnt/MP_*.cnt > MP_fake_cnt.txt

#rm ML_*_eqtl_fake.txt