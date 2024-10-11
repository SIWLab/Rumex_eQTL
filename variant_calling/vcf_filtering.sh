
# concat all autosomes together
cd /ohta2/meng.yuan/rumex/eqtl/VCF/
bcftools concat eqtl_mpileup_A1.vcf.gz eqtl_mpileup_A2.vcf.gz \
eqtl_mpileup_A3.vcf.gz eqtl_mpileup_A4.vcf.gz \
--threads 20 | bgzip -c > eqtl_mpileup_auto.vcf.gz
tabix eqtl_mpileup_auto.vcf.gz

bcftools query -l eqtl_mpileup_auto.vcf.gz > sample_names.txt # 152 sample


# reheader, and separate M and F samples
# drop mismatched and questionalbe samples 35aMLD 40aFLD
# keep only sample with matching genotype-phenotype pairs 
# (ie, has both DNA and RNA samples)

# redo genotype files separately for MP (74) ML (75) FL(74)
# prep for genotype files for leaf samples ML (75) FL(74) = 149
# rename chrom names too
# X = 5, PAR = 6

# 4 files
# cd /ohta2/meng.yuan/rumex/eqtl/VCF/
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.vcf.gz
VCF_OUT1=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.FL.vcf.gz
VCF_OUT3=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.MP.vcf.gz
VCF_OUT4=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.vcf.gz

bcftools reheader ${VCF} --threads 20 -s sample_names.txt \
| bcftools view -S ML_75.txt -m2 -M2 -v snps --threads 20 |\
bcftools annotate --rename-chrs chr_names | bgzip -c > ${VCF_OUT1} 

bcftools reheader ${VCF} --threads 20 -s sample_names.txt \
| bcftools view -S FL_74.txt -m2 -M2 -v snps --threads 20 |\
bcftools annotate --rename-chrs chr_names | bgzip -c > ${VCF_OUT2} 

bcftools reheader ${VCF} --threads 20 -s sample_names.txt \
| bcftools view -S MP_73.txt -m2 -M2 -v snps |\
bcftools annotate --rename-chrs chr_names | bgzip -c > ${VCF_OUT3} 

bcftools reheader ${VCF} --threads 20 -s sample_names.txt \
| bcftools view -S sex_149.txt -m2 -M2 -v snps |\
bcftools annotate --rename-chrs chr_names | bgzip -c > ${VCF_OUT4} 


tabix ${VCF_OUT1}
tabix ${VCF_OUT2}
tabix ${VCF_OUT3}
tabix ${VCF_OUT4}

# removed outliers in males for plink later

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


VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.filt.vcf.gz
VCF2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.renamed.filt.vcf.gz
bcftools annotate --rename-chrs chr_namesY ${VCF} --threads 20 | bgzip -c > ${VCF2}


# X for female
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_X.vcf.gz
VCF_OUT=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_X.SNP.FL.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_X.SNP.FL.filt.vcf.gz

bcftools reheader ${VCF} --threads 20 -s sample_names.txt \
| bcftools view -S FL_74.txt -m2 -M2 -v snps --threads 20 |\
bcftools annotate --rename-chrs chr_namesX | bgzip -c > ${VCF_OUT} 

tabix ${VCF_OUT}

vcftools --gzvcf ${VCF_OUT} \
--minQ 30 \
--min-meanDP 5 \
--max-meanDP 20 \
--minGQ 30 \
--max-missing 0.8 \
--recode --stdout | bgzip -c > ${VCF_OUT2}

tabix ${VCF_OUT2}


# PAR for female
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_Y.vcf.gz
VCF_OUT=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.filt.vcf.gz

tabix ${VCF}
bcftools reheader ${VCF} --threads 20 -s sample_names.txt \
| bcftools view -S FL_74.txt -r Y:1-45000000 -m2 -M2 -v snps --threads 20 |\
bcftools annotate --rename-chrs chr_namesY | bgzip -c > ${VCF_OUT} 

tabix ${VCF_OUT}

vcftools --gzvcf ${VCF_OUT} \
--minQ 30 \
--min-meanDP 5 \
--max-meanDP 20 \
--minGQ 30 \
--max-missing 0.8 \
--recode --stdout | bgzip -c > ${VCF_OUT2}

tabix ${VCF_OUT2}






