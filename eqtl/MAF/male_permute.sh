#!/bin/bash
# usage sh male_permute.sh 100
# do ML and MP permutations in one folder
iterations=$1
cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/Male

run_task() {
i=$1
########################## get fake eqtls in ML #############################
# reshuffle phenotype in R for both ML and MP
/usr/local/bin/Rscript --vanilla pheno_permute_m.R ${i}

# run tensorqtl for both
bgzip normalized_counts_mln_auto_${i}.bed && tabix -p bed normalized_counts_mln_auto_${i}.bed.gz
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/ML/normalized_counts_mln_auto_${i}.bed.gz
prefix=ML_${i}
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
--covariates ${covariates_file} \
--mode cis_nominal --window 20000

bgzip normalized_counts_mpn_auto_${i}.bed && tabix -p bed normalized_counts_mpn_auto_${i}.bed.gz
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/MP/normalized_counts_mpn_auto_${i}.bed.gz
prefix=MP_${i}
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
--covariates ${covariates_file} \
--mode cis_nominal --window 20000


# run MAF_null scripts in R to get fake eqtls
/usr/local/bin/Rscript --vanilla get_MAF_null_m.R ${i}

# remove intermediate files
rm normalized_counts_*_auto_${i}.bed.gz*
rm M*_${i}.cis_qtl_pairs.*.parquet
rm M*_${i}.tensorQTL.cis_nominal.log

mv *_fake.txt fake/ 
mv *_random.cnt cnt/

echo "Running task $i"
}

export -f run_task  

# Use GNU parallel to run 15 iterations at the same time
seq 1 $iterations | parallel --tmpdir /ohta2/meng.yuan/tmp/ --joblog parallel_log.txt -j 2 run_task 
