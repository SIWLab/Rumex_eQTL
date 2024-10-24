############################## genotype data ##############################
# genotype data
# for autosomes, genotype files are already generated in plink_PCA.sh
# genotype data for X and PAR for females
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_X.SNP.FL.filt.vcf.gz
i=FL_X
plink --vcf ${VCF} --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out ${i}

VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.filt.vcf.gz
i=FL_PAR
plink --vcf ${VCF} --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out ${i}

############################## phenotype data ##############################
# phenotype data
grep 'maker\stranscript' merged_TX_noMatPARlarge_txanno.gtf > merged_TX_noMatPARlarge_txanno_gene_full.gtf
cut -f 1,4,5,7  merged_TX_noMatPARlarge_txanno_gene_full.gtf > merged_TX_noMatPARlarge_txanno_gene_full.bed

# already sorted 
# rename chrom names for simplicity with plink
for i in "ml" "fl" "mp"
do
file=normalized_counts_${i}n_auto.bed
sed -i 's/A1/1/g' ${file}
sed -i 's/A2/2/g' ${file}
sed -i 's/A3/3/g' ${file}
sed -i 's/A4/4/g' ${file}

bgzip normalized_counts_${i}n_auto.bed && tabix -p bed normalized_counts_${i}n_auto.bed.gz
done

# rename chrom for PAR and X so there's no issue with plink
# can rename Y for 7 later
sed -i 's/X/5/g' normalized_counts_fln_X.bed
sed -i 's/PAR/6/g' normalized_counts_fln_PAR.bed
bgzip normalized_counts_fln_X.bed && tabix -p bed normalized_counts_fln_X.bed.gz
bgzip normalized_counts_fln_PAR.bed && tabix -p bed normalized_counts_fln_PAR.bed.gz

############################## run tensorqtl ##############################
# male leaf
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/ML
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mln_auto.bed.gz
prefix=ML
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1
#   * 71 samples
#   * 14940 phenotypes
#   * 7 covariates
#   * 3524400 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 14940/14940
#     ** dropping 148 phenotypes without variants in cis-window
#   * computing permutations
#     processing phenotype 14792/14792
#   Time elapsed: 9.44 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 14792
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.49
#   * QTL phenotypes @ FDR 0.05: 1573
#   * min p-value threshold @ FDR 0.05: 0.0104994

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
# * 71 samples
#   * 14940 phenotypes
#   * 7 covariates
#   * 3524400 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 14940/14940
#     ** dropping 148 phenotypes without variants in cis-window
#   * Computing associations
#     Mapping chromosome 1
#     processing phenotype 5023/14792    time elapsed: 0.09 min
#     * writing output
#     Mapping chromosome 2
#     processing phenotype 9734/14792    time elapsed: 0.19 min
#     * writing output
#     Mapping chromosome 3
#     processing phenotype 12447/14792    time elapsed: 0.25 min
#     * writing output
#     Mapping chromosome 4
#     processing phenotype 14792/14792
#     time elapsed: 0.30 min
#     * writing output


# pollen (MP)
# genotype files are the same for male leaf and pollen
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/ML
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mpn_auto.bed.gz
prefix=MP
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
#   Time elapsed: 8.62 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 13785
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.35
#   * QTL phenotypes @ FDR 0.05: 958
#   * min p-value threshold @ FDR 0.05: 0.00537284

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
  #   processing phenotype 4616/13785    time elapsed: 0.09 min
  #   * writing output
  #   Mapping chromosome 2
  #   processing phenotype 9049/13785    time elapsed: 0.18 min
  #   * writing output
  #   Mapping chromosome 3
  #   processing phenotype 11586/13785    time elapsed: 0.24 min
  #   * writing output
  #   Mapping chromosome 4
  #   processing phenotype 13785/13785


# female leaf (autosomes)
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_fln_auto.bed.gz
prefix=FL
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1
#   * 74 samples
#   * 14634 phenotypes
#   * 5 covariates
#   * 4509904 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 14634/14634
#     ** dropping 70 phenotypes without variants in cis-window
#   * computing permutations
#     processing phenotype 14564/14564
#   Time elapsed: 9.64 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 14564
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.53
#   * QTL phenotypes @ FDR 0.05: 2226
#   * min p-value threshold @ FDR 0.05: 0.0162608

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
  # * 74 samples
  # * 14634 phenotypes
  # * 5 covariates
  # * 4509904 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 14634/14634
  #   ** dropping 70 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 1
  #   processing phenotype 4922/14564    time elapsed: 0.10 min
  #   * writing output
  #   Mapping chromosome 2
  #   processing phenotype 9557/14564    time elapsed: 0.20 min
  #   * writing output
  #   Mapping chromosome 3
  #   processing phenotype 12227/14564    time elapsed: 0.26 min
  #   * writing output
  #   Mapping chromosome 4
  #   processing phenotype 14564/14564

# X and PAR for female
# X
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL_X
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_fln_X.bed.gz
prefix=FL_X
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1


# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000

# PAR
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL_PAR
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_fln_PAR.bed.gz
prefix=FL_PAR
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1
  # * 74 samples
  # * 1102 phenotypes
  # * 0 covariates
  # * 96911 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 1102/1102
  #   ** dropping 12 phenotypes without variants in cis-window
  # * computing permutations
  #   processing phenotype 1090/1090

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
  # * 74 samples
  # * 1102 phenotypes
  # * 0 covariates
  # * 96911 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 1102/1102
  #   ** dropping 12 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 1
  #   processing phenotype 1090/1090






