############################## genotype data ##############################
# genotype data
# for autosomes, genotype files are already generated in plink_PCA.sh
# genotype data for X and PAR for females
cd /ohta2/meng.yuan/rumex/eqtl/plink
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
grep 'Ontology_term' merged_TX_noMatPARlarge.gff > merged_TX_noMatPARlarge_GO.gff
#grep 'Liftoff\sgene' merged_TX_noMatPARlarge_GO.gff > merged_TX_noMatPARlarge_gene_GO.gff # 12684

grep 'maker\sgene' merged_TX_noMatPARlarge_txanno_GO.gff > merged_TX_noMatPARlarge_txanno_gene_GO.gff # 13616
grep 'maker\sgene' merged_TX_noMatPARlarge_txanno.gff > merged_TX_noMatPARlarge_txanno_gene.gff # 37659

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
# autosomes
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mln_auto.bed.gz
prefix=ML
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
#   Time elapsed: 9.46 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 14777
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.48
#   * QTL phenotypes @ FDR 0.10: 2425
#   * min p-value threshold @ FDR 0.1: 0.0318499

# old  
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
# * 69 samples
#   * 14949 phenotypes
#   * 6 covariates
#   * 3288995 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 14949/14949
#     ** dropping 172 phenotypes without variants in cis-window
#   * Computing associations
#     Mapping chromosome 1
#     processing phenotype 5024/14777    time elapsed: 0.09 min
#     * writing output
#     Mapping chromosome 2
#     processing phenotype 9733/14777    time elapsed: 0.19 min
#     * writing output
#     Mapping chromosome 3
#     processing phenotype 12445/14777    time elapsed: 0.25 min
#     * writing output
#     Mapping chromosome 4
#     processing phenotype 14777/14777

# Y
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/Y
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mln_Y.bed.gz
prefix=Y
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000    

# pollen (MP)
# genotype files are the same for male leaf and pollen
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
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
#   Time elapsed: 8.63 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 13785
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.35
#   * QTL phenotypes @ FDR 0.10: 1481
#   * min p-value threshold @ FDR 0.1: 0.0165895

# old  
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
#     processing phenotype 9823/14564
#     processing phenotype 14564/14564
#   Time elapsed: 9.90 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 14564
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.53
#   * QTL phenotypes @ FDR 0.10: 3445
#   * min p-value threshold @ FDR 0.1: 0.0502758

# old
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
#   * 74 samples
#   * 1894 phenotypes
#   * 5 covariates
#   * 769059 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 1894/1894
#     ** dropping 24 phenotypes without variants in cis-window
#   * computing permutations
#     processing phenotype 1870/1870
#   Time elapsed: 1.30 min
# done.
#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 1870
#   * Correlation between Beta-approximated and empirical p-values: 0.9999
#   * Proportion of significant phenotypes (1-pi0): 0.77
#   * QTL phenotypes @ FDR 0.10: 1294
#   * min p-value threshold @ FDR 0.1: 0.294671

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
  # * 74 samples
  # * 1894 phenotypes
  # * 5 covariates
  # * 769059 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 1894/1894
  #   ** dropping 24 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 5
  #   processing phenotype 1870/1870
    
# PAR
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL_PAR
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_fln_PAR.bed.gz
prefix=FL_PAR
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000 --fdr 0.1

#   * 74 samples
#   * 1102 phenotypes
#   * 5 covariates
#   * 96911 variants
#   * cis-window: ±20,000
#   * checking phenotypes: 1102/1102
#     ** dropping 12 phenotypes without variants in cis-window
#   * computing permutations
#     processing phenotype 1090/1090

#   * writing output
# Computing q-values
#   * Number of phenotypes tested: 1090
#   * Correlation between Beta-approximated and empirical p-values: 1.0000
#   * Proportion of significant phenotypes (1-pi0): 0.36
#   * QTL phenotypes @ FDR 0.10: 111
#   * min p-value threshold @ FDR 0.1: 0.0163153

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000

  # * 74 samples
  # * 1102 phenotypes
  # * 5 covariates
  # * 96911 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 1102/1102
  #   ** dropping 12 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 6
  #   processing phenotype 1090/1090




