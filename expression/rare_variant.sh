
bcftools query -l /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt2.vcf.gz > ML.txt
bcftools query -l /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.FL.filt.vcf.gz > FL.txt

# extract singleton sites
VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt2.vcf.gz
vcftools --gzvcf ${VCF_IN} --singletons --out ML
cut -f 1,2 ML.singletons > ML.singletons.txt

vcftools --gzvcf ${VCF_IN} --positions ML.singletons.txt --recode --stdout | bgzip -c > ML_singleton.vcf.gz
# After filtering, kept 27205538 out of a possible 64988762 Sites

VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.FL.filt.vcf.gz
vcftools --gzvcf ${VCF_IN} --singletons --out FL
cut -f 1,2 FL.singletons > FL.singletons.txt

vcftools --gzvcf ${VCF_IN} --positions FL.singletons.txt --recode --stdout | bgzip -c > FL_singleton.vcf.gz
# After filtering, kept 28553283 out of a possible 66828708 Sites

# keep only singleton not private doubleton
tabix ML_singleton.vcf.gz
bcftools view -R ML_genes_5k.bed ML_singleton.vcf.gz --threads 30 | \
bcftools query -f '%CHROM\t%POS[\t%GT]\n' > ML_singleton_5k_genotypes.tsv

tabix FL_singleton.vcf.gz
bcftools view -R FL_genes_5k.bed FL_singleton.vcf.gz --threads 30 | \
bcftools query -f '%CHROM\t%POS[\t%GT]\n' > FL_singleton_5k_genotypes.tsv

bcftools view -R MP_genes_5k.bed ML_singleton.vcf.gz --threads 30 | \
bcftools query -f '%CHROM\t%POS[\t%GT]\n' > MP_singleton_5k_genotypes.tsv


# last step
grep '0/1' ML_singleton_5k_genotypes.tsv> ML_singleton_5k_geno.tsv
grep '0/1' FL_singleton_5k_genotypes.tsv> FL_singleton_5k_geno.tsv
grep '0/1' MP_singleton_5k_genotypes.tsv> MP_singleton_5k_geno.tsv

comm -23 ML_singleton_5k_genotypes.tsv ML_singleton_5k_geno.tsv > private # 66
comm -23 FL_singleton_5k_genotypes.tsv FL_singleton_5k_geno.tsv > private2 # 252

grep '0/0' ML_singleton_5k_geno.tsv> ML_singleton_5k.tsv
grep '0/0' FL_singleton_5k_geno.tsv> FL_singleton_5k.tsv
grep '0/0' MP_singleton_5k_geno.tsv> MP_singleton_5k.tsv


comm -23 ML_singleton_5k_geno.tsv ML_singleton_5k.tsv > alt_homo # 7894
comm -23 FL_singleton_5k_geno.tsv FL_singleton_5k.tsv > alt_homo2 # 11658 remove for now

# awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1}' ML.singletons.sites > ML.singletons.bed
# bedtools intersect -a ML.singletons.bed -b ML_genes.bed > ML_sites.bed

# bedtools intersect -a ML.singletons.bed -b MP_genes.bed > MP_sites.bed


# awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1}' FL.singletons.sites > FL.singletons.bed
# bedtools intersect -a FL.singletons.bed -b FL_genes.bed > FL_sites.bed

   # 1676749 FL_sites.bed
   # 1668965 ML_sites.bed
   # 1560529 MP_sites.bed

1920203-1909879

# get gene interval regions only
# get simplifies genotype files
# ML

# VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt2.vcf.gz
# bcftools view -R ML.singletons.txt ${VCF_IN} --threads 30 | bgzip -c > ML2.vcf.gz
# tabix ML.vcf.gz
# bcftools view -R ML_genes.bed ML.vcf.gz --threads 30 | \
# bcftools query -f '%CHROM\t%POS[\t%GT]\n' > ML_rare_local_genotypes.tsv




# MAF < 5%
VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt2.vcf.gz
VCF_OUT=eqtl_mpileup_auto.SNP.ML.filt2_rare.vcf.gz
bcftools view -i 'MAF<0.05' ${VCF_IN} -Oz -o ${VCF_OUT}
bcftools view -R ML_genes.bed eqtl_mpileup_auto.SNP.ML.filt2_rare.vcf.gz --threads 30 | \
bcftools query -f '%CHROM\t%POS[\t%GT]\n' > ML_rare_local_genotypes3.tsv


VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt2.vcf.gz
bcftools view -i 'MAF<0.05' -R MP_genes.bed ${VCF_IN} --threads 30 | \
bcftools query -f '%CHROM\t%POS[\t%GT]\n' > MP_rare_local_genotypes3.tsv

VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.FL.filt.vcf.gz
bcftools view -i 'MAF<0.05' -R FL_genes.bed ${VCF_IN} --threads 30 | \
bcftools query -f '%CHROM\t%POS[\t%GT]\n' > FL_rare_local_genotypes3.tsv


# VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt2.vcf.gz
# bcftools view -R ML_sites.bed ${VCF_IN} --threads 30 | \
# bcftools query -f '%CHROM\t%POS[\t%GT]\n' > ML_rare_local_genotypes2.tsv


# MP
# bcftools view -R MP_genes.bed ML.vcf.gz --threads 30 | \
# bcftools query -f '%CHROM\t%POS[\t%GT]\n' > MP_rare_local_genotypes.tsv
# VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt2.vcf.gz
# bcftools view -R MP_sites.bed ${VCF_IN} --threads 30 | \
# bcftools query -f '%CHROM\t%POS[\t%GT]\n' > MP_rare_local_genotypes2.tsv

# FL
# VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.FL.filt.vcf.gz
# bcftools view -R FL.singletons.sites ${VCF_IN} --threads 30 | bgzip -c > FL.vcf.gz
# tabix FL.vcf.gz

# bcftools view -R FL_genes.bed FL.vcf.gz --threads 30 | \
# bcftools query -f '%CHROM\t%POS[\t%GT]\n' > FL_rare_local_genotypes.tsv
# VCF_IN=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.FL.filt.vcf.gz
# bcftools view -R FL_sites.bed ${VCF_IN} --threads 30 | \
# bcftools query -f '%CHROM\t%POS[\t%GT]\n' > FL_rare_local_genotypes2.tsv

