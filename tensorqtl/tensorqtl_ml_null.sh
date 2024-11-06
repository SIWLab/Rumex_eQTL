cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/ml

for i in {1..50}; do
file=normalized_counts_mln_auto_${i}.bed
sed -i 's/A1/1/g' ${file}
sed -i 's/A2/2/g' ${file}
sed -i 's/A3/3/g' ${file}
sed -i 's/A4/4/g' ${file}
bgzip normalized_counts_mln_auto_${i}.bed && tabix -p bed normalized_counts_mln_auto_${i}.bed.gz
done


# on slow server
for i in {5..50}; do
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/ml/normalized_counts_mln_auto_${i}.bed.gz
prefix=ML_${i}
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1
done


# on ohta
for i in {1..50}; do
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/ml/normalized_counts_mln_auto_${i}.bed.gz
prefix=ML_${i}
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
done
