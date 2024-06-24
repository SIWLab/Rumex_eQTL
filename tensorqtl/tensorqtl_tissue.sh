

plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mpn.bed.gz
prefix=MP
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_MP.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
