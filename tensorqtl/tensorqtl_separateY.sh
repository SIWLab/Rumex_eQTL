# eQTL mapping for the male sex chromosomes, X, Y, PAR
# genotypes
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_Y.SNP.M.filt2.vcf.gz
eqtl_mpileup_X.SNP.M.filt2.vcf.gz
eqtl_mpileup_PAR.SNP.M.filt2.vcf.gz
eqtl_mpileup_Ys.SNP.M.filt2.vcf.gz
/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.M.filt2.vcf.gz
/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_Ys.SNP.M.filt2.vcf.gz



# need to relax these filters on the Y??? ask stephen
i=M_Y
plink --vcf ${VCF} --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out ${i}

i=M_PAR

i=M_X

# phenotypes
sed -i 's/X/5/g' normalized_counts_fln_X.bed
sed -i 's/PAR/6/g' normalized_counts_fln_PAR.bed
bgzip normalized_counts_fln_X.bed && tabix -p bed normalized_counts_fln_X.bed.gz
bgzip normalized_counts_fln_PAR.bed && tabix -p bed normalized_counts_fln_PAR.bed.gz


