# genotype data for males
# HWE MAF 5%, missing lower 10%
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_autosomes.SNP_M.filt.vcf.gz

plink --vcf $VCF --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out autosomes_M

# 66461930 variants loaded from .bim file.
# Total genotyping rate is 0.931429.
# 3135844 variants removed due to missing genotype data (--geno).
# Warning: --hwe observation counts vary by more than 10%.  Consider using
# --geno, and/or applying different p-value thresholds to distinct subsets of
# your data.
# --hwe: 110207 variants removed due to Hardy-Weinberg exact test.
# 60041057 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 3174822 variants and 73 people pass filters and QC.

VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_autosomes.SNP_F.filt.vcf.gz

plink --vcf $VCF --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out autosomes_F

# 66461930 variants loaded from .bim file.
# 1464195 variants removed due to missing genotype data (--geno).
# Warning: --hwe observation counts vary by more than 10%.  Consider using
# --geno, and/or applying different p-value thresholds to distinct subsets of
# your data.
# --hwe: 128927 variants removed due to Hardy-Weinberg exact test.
# 61158223 variants removed due to minor allele threshold(s)
# (--maf/--max-maf/--mac/--max-mac).
# 3710585 variants and 74 people pass filters and QC.


# phenotype data
# saprate ML FL MP samples
## get tss from gff
#/ohta2/bianca.sacchi/rumex_ase_TX/genome/merged_TX_noMatPARlarge_txanno.gtf
#/ohta2/bianca.sacchi/annotation_TX/noMatPAR_liftoverAnno/merged_TX_noMatPARlarge.gff
grep 'maker\stranscript' merged_TX_noMatPARlarge_txanno.gtf > merged_TX_noMatPARlarge_txanno_gene_full.gtf
cut -f 1,4,5,7  merged_TX_noMatPARlarge_txanno_gene_full.gtf > merged_TX_noMatPARlarge_txanno_gene_full.bed

# already sorted normalized_counts_MPn.bed
bgzip normalized_counts_MPn.bed && tabix -p bed normalized_counts_MPn.bed.gz
bgzip normalized_counts_MLn.bed && tabix -p bed normalized_counts_MLn.bed.gz


# pollen
# for each chrom 
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/autosomes
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_MPn.bed.gz
prefix=autosomes_MP
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000



# cis-QTL mapping: conditionally independent QTLs
# python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
#     --covariates ${covariates_file} \
#     --cis_output ${prefix}.cis_qtl.txt.gz \
#     --mode cis_independent



# male leaf
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/autosomes
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_MLn.bed.gz
prefix=autosomes_ML
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000



