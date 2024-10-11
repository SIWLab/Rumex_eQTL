
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
