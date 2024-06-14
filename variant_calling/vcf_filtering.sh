
for i in "A1" "A2" "A3" "A4" "X" "Y"
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


# reheader, and separate M and F samples
# drop mismatched and questionalbe samples 35aMLD 40aFLD
# keep only sample with matching genotype-phenotype pairs 
# (ie, has both DNA and RNA samples)

# concat all autosomes together
cd /ohta2/meng.yuan/rumex/eqtl/VCF/
bcftools concat eqtl_mpileup_A1.SNP.filt.vcf.gz eqtl_mpileup_A2.SNP.filt.vcf.gz \
eqtl_mpileup_A3.SNP.filt.vcf.gz eqtl_mpileup_A4.SNP.filt.vcf.gz \
--threads 20 | bgzip -c > eqtl_mpileup_autosomes.SNP.filt.vcf.gz
tabix eqtl_mpileup_autosomes.SNP.filt.vcf.gz


VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_autosomes.SNP.filt.vcf.gz
VCF_OUT1=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_autosomes.SNP_M.filt.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_autosomes.SNP_F.filt.vcf.gz
bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S MLD.txt | bgzip -c > ${VCF_OUT1} 

bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S FLD.txt | bgzip -c > ${VCF_OUT2} 

# need to generate a VCF with both sexes and with sample ids sorted
bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S MLD_FLD.txt | bgzip -c > ${VCF_OUT} 


# use A4 as an example
CHROM="A4"
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${CHROM}.SNP.filt.vcf.gz
VCF_OUT1=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${CHROM}.SNP_M.filt.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${CHROM}.SNP_F.filt.vcf.gz

bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S MLD.txt | bgzip -c > ${VCF_OUT1} 

bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S FLD.txt | bgzip -c > ${VCF_OUT2} 

