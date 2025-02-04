
# autosomes for all inds
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.vcf.gz
VCF_OUT=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.invariant.L.vcf.gz

bcftools reheader ${VCF} --threads 20 -s sample_names.txt |\
bcftools view -S sex_147.txt | bcftools annotate --rename-chrs chr_names |\
vcftools --vcf - --max-maf 0 --max-alleles 2 --remove-indels --max-missing 0.8 \
--min-meanDP 5 --max-meanDP 20 --recode --stdout | bgzip -c > ${VCF_OUT} 

tabix ${VCF_OUT}

# remove genetic PC outliers
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.invariant.L.vcf.gz
VCF2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.invariant.L.filt.vcf.gz
bcftools view ${VCF} -s ^5aM,7bM,27eM,53bM --threads 20 | bgzip -c > ${VCF2}
tabix ${VCF2}

# combine the two VCFs using bcftools concat
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.filt2.vcf.gz
VCF2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.invariant.L.filt.vcf.gz

bcftools concat --allow-overlaps ${VCF} ${VCF2} \
-O z -o eqtl_mpileup_auto.allsites.L.vcf.gz
