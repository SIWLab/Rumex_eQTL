
# phenotype file
file=normalized_counts_ln_auto.bed
sed -i 's/A1/1/g' ${file}
sed -i 's/A2/2/g' ${file}
sed -i 's/A3/3/g' ${file}
sed -i 's/A4/4/g' ${file}

bgzip normalized_counts_ln_auto.bed && tabix -p bed normalized_counts_ln_auto.bed.gz

#############################################################################
# run tensorqtl with sex as covaraite
# no interaction
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/L
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_ln_auto.bed.gz
prefix=sex
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_sex.txt
# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1
#  * 143 samples
#   * 14992 phenotypes
#   * 2 covariates
#   * 3684937 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 14992/14992
#     ** dropping 237 phenotypes without variants in cis-window
#   * computing permutations
#     processing phenotype 14755/14755
#   Time elapsed: 10.31 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 14755
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.76
#   * QTL phenotypes @ FDR 0.10: 9743
#   * min p-value threshold @ FDR 0.1: 0.275477

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
# * 143 samples
#   * 14992 phenotypes
#   * 2 covariates
#   * 3684937 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 14992/14992
#     ** dropping 237 phenotypes without variants in cis-window
#   * Computing associations
#     Mapping chromosome 1
#     processing phenotype 5028/14755    time elapsed: 0.11 min
#     * writing output
#     Mapping chromosome 2
#     processing phenotype 9704/14755    time elapsed: 0.22 min
#     * writing output
#     Mapping chromosome 3
#     processing phenotype 12401/14755    time elapsed: 0.28 min
#     * writing output
#     Mapping chromosome 4
#     processing phenotype 14755/14755

#############################################################################
# run tensorqtl without sex as covaraite
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/L
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_ln_auto.bed.gz
prefix=nosex
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_L.txt
# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1
# cis-QTL mapping: empirical p-values for phenotypes
#   * 143 samples
#   * 14992 phenotypes
#   * 1 covariates
#   * 3684937 variants
#   * cis-window: ±20,000
#   * computing permutations
#   Time elapsed: 10.29 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 14755
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.75
#   * QTL phenotypes @ FDR 0.10: 9530
#   * min p-value threshold @ FDR 0.1: 0.261712

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
# * 143 samples
#   * 14992 phenotypes
#   * 1 covariates
#   * 3684937 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 14992/14992
#     ** dropping 237 phenotypes without variants in cis-window
#   * Computing associations
#     Mapping chromosome 1
#     processing phenotype 5028/14755    time elapsed: 0.86 min
#     * writing output
#     Mapping chromosome 2
#     processing phenotype 9704/14755    time elapsed: 2.17 min
#     * writing output
#     Mapping chromosome 3
#     processing phenotype 12401/14755    time elapsed: 2.44 min
#     * writing output
#     Mapping chromosome 4
#     processing phenotype 14755/14755

# run tensorqtl with GxSex interaction term
# using results from best-only
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/L
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_ln_auto.bed.gz
prefix=nosex
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_L.txt
interactions_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/interaction_sex.txt

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix}_inter \
    --covariates ${covariates_file} \
    --interaction ${interactions_file} \
    --mode cis_nominal --window 20000
  # * 143 samples
  # * 14992 phenotypes
  # * 1 covariates
  # * 3684937 variants
  # * including 1 interaction term(s)
  #   * using 0.05 MAF threshold
  # * cis-window: ±20,000
  # * checking phenotypes: 14992/14992
  #   ** dropping 237 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 1
  #   processing phenotype 5028/14755    time elapsed: 4.49 min
  #   * writing output
  #   Mapping chromosome 2
  #   processing phenotype 9704/14755    time elapsed: 8.49 min
  #   * writing output
  #   Mapping chromosome 3
  #   processing phenotype 12401/14755    time elapsed: 10.24 min
  #   * writing output
  #   Mapping chromosome 4
  #   processing phenotype 14755/14755
  #   time elapsed: 11.69 min
  #   * writing output

# keep only top association
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix}_inter \
    --covariates ${covariates_file} \
    --interaction ${interactions_file} \
    --best_only \
    --mode cis_nominal --window 20000

  # * 143 samples
  # * 14992 phenotypes
  # * 1 covariates
  # * 3684937 variants
  # * including 1 interaction term(s)
  #   * using 0.05 MAF threshold
  # * cis-window: ±20,000
  # * checking phenotypes: 14992/14992
  #   ** dropping 237 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 1
  #   processing phenotype 5028/14755    time elapsed: 3.88 min
  #   Mapping chromosome 2
  #   processing phenotype 9704/14755    time elapsed: 8.86 min
  #   Mapping chromosome 3
  #   processing phenotype 12401/14755    time elapsed: 10.87 min
  #   Mapping chromosome 4
  #   processing phenotype 14755/14755
