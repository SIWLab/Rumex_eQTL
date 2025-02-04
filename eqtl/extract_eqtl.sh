
bcftools view /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt.vcf.gz \
-r A1:274220346 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A1_274220346.vcf


VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.filt.vcf.gz
bcftools view ${VCF} -r A1:159008371 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A1_159008371.vcf



#A4:30274594 A3:13387496 A1:455264496

VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt.vcf.gz

bcftools view ${VCF} -r A4:30274594 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A4_30274594.vcf

bcftools view ${VCF} -r A3:13387496 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A3_13387496.vcf


bcftools view ${VCF} -r A1:455264496 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A1_455264496.vcf

bcftools view /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt.vcf.gz \
-r A1:274220346 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A1_274220346.vcf



VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt.vcf.gz

tabix ${VCF}
bcftools view ${VCF} -r 4:79507645 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A4_79507645.vcf

# G x Sex
bcftools view ${VCF} -r 1:70371578 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A1_70371578.vcf

bcftools view ${VCF} -r 1:159008371 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A1_159008371.vcf

bcftools view ${VCF} -r 2:14033315 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A2_14033315.vcf

bcftools view ${VCF} -r 2:250749626 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A2_250749626.vcf

bcftools view ${VCF} -r 3:10654019 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A3_10654019.vcf


# G x Lifestage
1:255437493

1:455264496
2:162241243
2:210629115

2:238366426
3:13387496
3:13442544
3:75629497
4:147510039







