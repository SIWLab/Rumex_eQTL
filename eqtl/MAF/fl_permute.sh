#!/bin/bash
# usage sh fl_permute.sh 1000
iterations=$1
cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/FL

run_task() {
i=$1
# reshuffle phenotype in R
/usr/local/bin/Rscript --vanilla pheno_permute_fl.R ${i}
# run tensorqtl
bgzip normalized_counts_fln_auto_${i}.bed && tabix -p bed normalized_counts_fln_auto_${i}.bed.gz

plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/FL/normalized_counts_fln_auto_${i}.bed.gz
prefix=FL_${i}
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL.txt
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
--covariates ${covariates_file} \
--mode cis_nominal --window 20000

# run MAF_null scripts in R to get FL_i_eqtl_fake.txt
/usr/local/bin/Rscript --vanilla get_MAF_null_fl.R ${i}

# remove intermediate files
rm normalized_counts_fln_auto_${i}.bed.gz*
rm FL_${i}.cis_qtl_pairs.*.parquet
rm FL_${i}.tensorQTL.cis_nominal.log

mv FL_${i}_random.cnt cnt/
mv FL_${i}_eqtl_random_fake.txt fake/
echo "Running task $i"
}

export -f run_task  

# Use GNU parallel to run 15 iterations at the same time
seq 1 $iterations | parallel --tmpdir /ohta2/meng.yuan/tmp/ --joblog parallel_log.txt -j 15 run_task 