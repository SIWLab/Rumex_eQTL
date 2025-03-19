#!/bin/bash
# usage sh ml_permute.sh 1000
iterations=$1
cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/ML

run_task() {
i=$1
# reshuffle phenotype in R
/usr/local/bin/Rscript --vanilla pheno_permute.R ${i}
# run tensorqtl
bgzip normalized_counts_mln_auto_${i}.bed && tabix -p bed normalized_counts_mln_auto_${i}.bed.gz

plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/ML/normalized_counts_mln_auto_${i}.bed.gz
prefix=ML_${i}
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
--covariates ${covariates_file} \
--mode cis_nominal --window 20000

# run MAF_null scripts in R to get ML_i_eqtl_fake.txt
/usr/local/bin/Rscript --vanilla get_MAF_null.R ${i}

# remove intermediate files
rm normalized_counts_mln_auto_${i}.bed.gz*
rm ML_${i}.cis_qtl_pairs.*.parquet
rm ML_${i}.tensorQTL.cis_nominal.log
 
mv ML_${i}_random.cnt* cnt/
mv ML_${i}_eqtl_random_fake.txt fake/
echo "Running task $i"
}

export -f run_task  

# Use GNU parallel to run 15 iterations at the same time
seq 1 $iterations | parallel --tmpdir /ohta2/meng.yuan/tmp/ --joblog parallel_log.txt -j 15 run_task 