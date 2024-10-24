
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

# redo genotype files separately for MP (73) ML (73) FL(74)
# prep for genotype files for leaf samples ML (73) FL(74) = 147
# rename chrom names too
# X = 5, PAR = 6


# cd /ohta2/meng.yuan/rumex/eqtl/VCF/
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.vcf.gz
VCF_OUT1=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.FL.vcf.gz
VCF_OUT3=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.vcf.gz

bcftools reheader ${VCF} --threads 20 -s sample_names.txt \
| bcftools view -S ML_73.txt -m2 -M2 -v snps --threads 20 |\
bcftools annotate --rename-chrs chr_names | bgzip -c > ${VCF_OUT1} 

bcftools reheader ${VCF} --threads 20 -s sample_names.txt \
| bcftools view -S FL_74.txt -m2 -M2 -v snps --threads 20 |\
bcftools annotate --rename-chrs chr_names | bgzip -c > ${VCF_OUT2} 

bcftools reheader ${VCF} --threads 20 -s sample_names.txt \
| bcftools view -S sex_147.txt -m2 -M2 -v snps |\
bcftools annotate --rename-chrs chr_names | bgzip -c > ${VCF_OUT3} 


tabix ${VCF_OUT1}
tabix ${VCF_OUT2}
tabix ${VCF_OUT3}


for i in "ML" "FL" "L"
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

bcftools stats /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.filt.vcf.gz \
> /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.filt.vcf.stats


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
i="Y"
# remove one sample with extremely low coverage
bcftools view /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}.vcf.gz \
-s ^NS.1594.001.IDT_i7_141---IDT_i5_141.24fMLD -r Y:1-45000000 --threads 20 -m2 -M2 -v snps \
| bgzip -c > /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR_SNP.vcf.gz


VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_Y.vcf.gz
VCF_OUT=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.vcf.gz
VCF_OUT2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.filt.vcf.gz

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


