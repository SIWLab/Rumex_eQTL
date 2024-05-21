
for i in "A1" "A2" "A3" "A4" "X"
do
vcf_in=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}_SNP.vcf.gz
vcf_out=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}.SNP.filt.vcf.gz

vcftools --gzvcf ${vcf_in} \
--minQ 30 \
--min-meanDP 5 \
--max-meanDP 20 \
--minGQ 30 \
--max-missing 0.8 \
--recode --stdout | bgzip -c > ${vcf_out}

tabix ${vcf_out}
done


