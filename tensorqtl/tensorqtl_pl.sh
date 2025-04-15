# flip the alt allele as the nminor allele
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt2.vcf.gz
/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.MP.filt.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr  \
--set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out ML_flipped

plink --bfile ML_flipped --double-id --allow-extra-chr  \
--a2-allele ML_flipped.bim 5 2 \
--make-bed --out ML_flipped2

# double check
plink --bfile ML_flipped2 --double-id --allow-extra-chr \
--freq --out freq_check


# absolute value of difference
file=normalized_counts_pl_diff_abs.bed
sed -i 's/A1/1/g' ${file}
sed -i 's/A2/2/g' ${file}
sed -i 's/A3/3/g' ${file}
sed -i 's/A4/4/g' ${file}
bgzip normalized_counts_pl_diff_abs.bed && tabix -p bed normalized_counts_pl_diff_abs.bed.gz

plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/ML_flipped2
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_pl_diff_abs.bed.gz
prefix=pl_abs
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
#   Time elapsed: 10.36 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 16218
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.40
#   * QTL phenotypes @ FDR 0.10: 1566
#   * min p-value threshold @ FDR 0.1: 0.0161976

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
  #   processing phenotype 5501/16218    time elapsed: 0.90 min
  #   * writing output
  #   Mapping chromosome 2
  #   processing phenotype 10709/16218    time elapsed: 1.53 min
  #   * writing output
  #   Mapping chromosome 3
  #   processing phenotype 13676/16218    time elapsed: 1.89 min
  #   * writing output
  #   Mapping chromosome 4
  #   processing phenotype 16218/16218
  #   time elapsed: 2.18 min
  #   * writing output

# cis independent
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --cis_output ${prefix}.cis_qtl.txt.gz \
    --mode cis_independent --window 20000 --fdr 0.1
  # * 69 samples
  # * 1566/16218 significant phenotypes
  # * 6 covariates
  # * 3288995 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 1566/1566
  # * computing independent QTLs
  #   processing phenotype 1566/1566

# do the same allele flipping in ML and MP as well
# male leaf
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/ML_flipped2
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mln_auto.bed.gz
prefix=ML_flipped
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1
#   * 69 samples
#   * 14949 phenotypes
#   * 6 covariates
#   * 3288995 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 14949/14949
#     ** dropping 172 phenotypes without variants in cis-window
#   * computing permutations
#     processing phenotype 14777/14777
#   Time elapsed: 9.48 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 14777
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.48
#   * QTL phenotypes @ FDR 0.10: 2427
#   * min p-value threshold @ FDR 0.1: 0.0316567    

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
  # * 69 samples
  # * 14949 phenotypes
  # * 6 covariates
  # * 3288995 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 14949/14949
  #   ** dropping 172 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 1
  #   processing phenotype 5024/14777    time elapsed: 0.90 min
  #   * writing output
  #   Mapping chromosome 2
  #   processing phenotype 9733/14777    time elapsed: 1.55 min
  #   * writing output
  #   Mapping chromosome 3
  #   processing phenotype 12445/14777    time elapsed: 1.92 min
  #   * writing output
  #   Mapping chromosome 4
  #   processing phenotype 14777/14777
    
# pollen 
# genotype files are the same for male leaf and pollen
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/ML_flipped2
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mpn_auto.bed.gz
prefix=MP_flipped
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1
#   * 69 samples
#   * 13923 phenotypes
#   * 6 covariates
#   * 3288995 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 13923/13923
#     ** dropping 138 phenotypes without variants in cis-window
#   * computing permutations
#     processing phenotype 13785/13785
#   Time elapsed: 8.84 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 13785
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.35
#   * QTL phenotypes @ FDR 0.10: 1480
#   * min p-value threshold @ FDR 0.1: 0.0165735

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
  # * 69 samples
  # * 13923 phenotypes
  # * 6 covariates
  # * 3288995 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 13923/13923
  #   ** dropping 138 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 1
  #   processing phenotype 4616/13785    time elapsed: 0.84 min
  #   * writing output
  #   Mapping chromosome 2
  #   processing phenotype 9049/13785    time elapsed: 1.41 min
  #   * writing output
  #   Mapping chromosome 3
  #   processing phenotype 11586/13785    time elapsed: 1.74 min
  #   * writing output
  #   Mapping chromosome 4
  #   processing phenotype 13785/13785
