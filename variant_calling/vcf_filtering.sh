
# concat all autosomes together
cd /ohta2/meng.yuan/rumex/eqtl/VCF/
bcftools concat eqtl_mpileup_A1_SNP.vcf.gz eqtl_mpileup_A2_SNP.vcf.gz \
eqtl_mpileup_A3_SNP.vcf.gz eqtl_mpileup_A4_SNP.vcf.gz \
--threads 20 | bgzip -c > eqtl_mpileup_auto.SNP.vcf.gz
tabix eqtl_mpileup_auto.SNP.vcf.gz

# redo genotype files separately for MP (74) ML (75) FL(74)
# prep for genotype files for leaf samples ML (75) FL(74) = 149
# 4 files
# cd /ohta2/meng.yuan/rumex/eqtl/VCF/
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.vcf.gz
VCF_OUT1=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.FL.vcf.gz
VCF_OUT3=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.MP.vcf.gz
VCF_OUT4=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.vcf.gz

bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S ML_75.txt | bgzip -c > ${VCF_OUT1} 

bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S FL_74.txt | bgzip -c > ${VCF_OUT2} 

bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S MP_73.txt | bgzip -c > ${VCF_OUT3} 

bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S sex_149.txt | bgzip -c > ${VCF_OUT4} 


tabix ${VCF_OUT1}
tabix ${VCF_OUT2}
tabix ${VCF_OUT3}
tabix ${VCF_OUT4}


for i in "ML" "FL" "MP" "L"
do 
vcf_in=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.${i}.vcf.gz
vcf_out=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.${i}.filt.vcf.gz

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


# cd /ohta2/meng.yuan/rumex/eqtl/VCF/
# X
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_X_SNP.vcf.gz
VCF_OUT=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_X.SNP.FL.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_X.SNP.FL.filt.vcf.gz

bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S FL_74.txt | bgzip -c > ${VCF_OUT} 

tabix ${VCF_OUT}

vcftools --gzvcf ${VCF_OUT} \
--minQ 30 \
--min-meanDP 5 \
--max-meanDP 20 \
--minGQ 30 \
--max-missing 0.8 \
--recode --stdout | bgzip -c > ${VCF_OUT2}

tabix ${VCF_OUT2}


# get PAR for female
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR_SNP.vcf.gz
VCF_OUT=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.filt.vcf.gz

bcftools reheader ${VCF} --threads 20 -s samples.txt \
| bcftools view -S FL_74.txt | bgzip -c > ${VCF_OUT} 

tabix ${VCF_OUT}

vcftools --gzvcf ${VCF_OUT} \
--minQ 30 \
--min-meanDP 5 \
--max-meanDP 20 \
--minGQ 30 \
--max-missing 0.8 \
--recode --stdout | bgzip -c > ${VCF_OUT2}

tabix ${VCF_OUT2}





