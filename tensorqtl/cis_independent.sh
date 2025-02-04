
# male leaf
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mln_auto.bed.gz
prefix=ML
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

# cis-QTL mapping: permutations
# python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
#     --covariates ${covariates_file} \
#     --mode cis --window 20000 --fdr 0.1

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --cis_output ${prefix}.cis_qtl.txt.gz \
    --mode cis_independent --window 20000 --fdr 0.1

  # * 69 samples
  # * 2425/14777 significant phenotypes
  # * 6 covariates
  # * 3288995 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 2425/2425
  # * computing independent QTLs
  #   processing phenotype 2425/2425

# ML.cis_independent_qtl.txt.gz    

# gzip *.cis_qtl.txt
# female leaf
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_fln_auto.bed.gz
prefix=FL
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL.txt

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --cis_output ${prefix}.cis_qtl.txt.gz \
    --mode cis_independent --window 20000 --fdr 0.1
  # * 74 samples
  # * 3445/14564 significant phenotypes
  # * 5 covariates
  # * 4509904 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 3445/3445
  # * computing independent QTLs
  #   processing phenotype 3445/3445

# male pollen
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mpn_auto.bed.gz
prefix=MP
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --cis_output ${prefix}.cis_qtl.txt.gz \
    --mode cis_independent --window 20000 --fdr 0.1
  # * 69 samples
  # * 1481/13785 significant phenotypes
  # * 6 covariates
  # * 3288995 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 1481/1481
  # * computing independent QTLs
  #   processing phenotype 1481/1481
