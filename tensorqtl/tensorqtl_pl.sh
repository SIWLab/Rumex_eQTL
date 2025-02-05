
# new phenotype files
file=normalized_counts_pl_diff.bed
sed -i 's/A1/1/g' ${file}
sed -i 's/A2/2/g' ${file}
sed -i 's/A3/3/g' ${file}
sed -i 's/A4/4/g' ${file}
bgzip normalized_counts_pl_diff.bed && tabix -p bed normalized_counts_pl_diff.bed.gz

plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_pl_diff.bed.gz
prefix=pl
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1

#   * 69 samples
#   * 16411 phenotypes
#   * 6 covariates
#   * 3288995 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 16411/16411
#     ** dropping 193 phenotypes without variants in cis-window
#   * computing permutations
#     processing phenotype 16218/16218
#   Time elapsed: 10.33 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 16218
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.45
#   * QTL phenotypes @ FDR 0.10: 2031
#   * min p-value threshold @ FDR 0.1: 0.0226032

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000

  # * 69 samples
  # * 16411 phenotypes
  # * 6 covariates
  # * 3288995 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 16411/16411
  #   ** dropping 193 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 1
  #   processing phenotype 5501/16218    time elapsed: 0.10 min
  #   * writing output
  #   Mapping chromosome 2
  #   processing phenotype 10709/16218    time elapsed: 0.21 min
  #   * writing output
  #   Mapping chromosome 3
  #   processing phenotype 13676/16218    time elapsed: 0.27 min
  #   * writing output
  #   Mapping chromosome 4
  #   processing phenotype 16218/16218
  #   time elapsed: 0.32 min

# cis independent
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --cis_output ${prefix}.cis_qtl.txt.gz \
    --mode cis_independent --window 20000 --fdr 0.1

  # * 69 samples
  # * 2031/16218 significant phenotypes
  # * 6 covariates
  # * 3288995 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 2031/2031
  # * computing independent QTLs
  #   processing phenotype 2031/2031
    
