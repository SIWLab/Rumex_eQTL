bcftools view /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt.vcf.gz \
-r A1:274220346 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A1_274220346.vcf


VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.filt.vcf.gz
bcftools view ${VCF} -r A1:159008371 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A1_159008371.vcf