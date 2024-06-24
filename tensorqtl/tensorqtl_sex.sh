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


#
covariate_sex.txt

# phenotype file
bgzip normalized_counts_ln.bed && tabix -p bed normalized_counts_ln.bed.gz


# run tensorqtl
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/L
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_ln.bed.gz
prefix=sex
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_sex.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000

