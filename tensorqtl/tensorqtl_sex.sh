# genotype file
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.filt.vcf.gz

plink --vcf $VCF --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out sex

# 65912121 variants loaded from .bim file.
# Total genotyping rate is 0.936615.
# 0 variants removed due to missing genotype data (--geno).
# --hwe: 187950 variants removed due to Hardy-Weinberg exact test.
# 62072649 variants removed due to minor allele threshold(s)
# 3651522 variants and 149 people pass filters and QC.


# covariate_sex.txt

# phenotype file
bgzip normalized_counts_ln.bed && tabix -p bed normalized_counts_ln.bed.gz


# run tensorqtl
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/sex
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_ln.bed.gz
prefix=sex
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_sex.txt
interactions_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/interaction_sex.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000

 # * 149 samples
 #  * 14804 phenotypes
 #  * 1 covariates
 #  * 3651522 variants
 #  * cis-window: ±20,000
 #  * checking phenotypes: 14804/14804
 #    ** dropping 255 phenotypes without variants in cis-window
 #  * computing permutations
 #    processing phenotype 14549/14549


plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/sex
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_ln.bed.gz
prefix=sex
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_sex.txt
interactions_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/interaction_sex.txt


# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000

  # * 149 samples
  # * 14804 phenotypes
  # * 1 covariates
  # * 3651522 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 14804/14804
  #   ** dropping 255 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome A1
  #   processing phenotype 4953/14549    time elapsed: 0.68 min
  #   * writing output
  #   Mapping chromosome A2
  #   processing phenotype 9561/14549    time elapsed: 1.13 min
  #   * writing output
  #   Mapping chromosome A3
  #   processing phenotype 12236/14549    time elapsed: 1.40 min
  #   * writing output
  #   Mapping chromosome A4
  #   processing phenotype 14549/14549


# add interaction term
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix}_interaction \
    --covariates ${covariates_file} \
    --interaction ${interactions_file} \
    --mode cis_nominal --window 20000


python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix}_interaction \
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
  # * Computing associations
  #   Mapping chromosome A1
  #   processing phenotype 4953/14549    time elapsed: 3.04 min
  #   * writing output
  #   Mapping chromosome A2
  #   processing phenotype 9561/14549    time elapsed: 5.22 min
  #   * writing output
  #   Mapping chromosome A3
  #   processing phenotype 12236/14549    time elapsed: 6.95 min
  #   * writing output
  #   Mapping chromosome A4
  #   processing phenotype 14549/14549
