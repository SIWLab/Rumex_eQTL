
# phenotype file
i=l
file=normalized_counts_${i}n_auto.bed
sed -i 's/A1/1/g' ${file}
sed -i 's/A2/2/g' ${file}
sed -i 's/A3/3/g' ${file}
sed -i 's/A4/4/g' ${file}

bgzip normalized_counts_${i}n_auto.bed && tabix -p bed normalized_counts_${i}n_auto.bed.gz




# run tensorqtl without sex as covaraite
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/sex
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_ln_auto.bed.gz
prefix=nosex
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_nosex.txt
# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000


# run tensorqtl with GxSex interaction term
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/sex
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_ln_auto.bed.gz
prefix=nosex
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_nosex.txt
interactions_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/interaction_sex.txt

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix}_inter \
    --covariates ${covariates_file} \
    --interaction ${interactions_file} \
    --mode cis_nominal --window 20000

# keep only top association
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix}_inter \
    --covariates ${covariates_file} \
    --interaction ${interactions_file} \
    --best_only \
    --mode cis_nominal --window 20000

  # * 149 samples
  # * 14804 phenotypes
  # * 1 covariates
  # * 3651522 variants
  # * including 1 interaction term(s)
  #   * using 0.05 MAF threshold
  # * cis-window: ±20,000
  # * checking phenotypes: 14804/14804
  #   ** dropping 255 phenotypes without variants in cis-window

