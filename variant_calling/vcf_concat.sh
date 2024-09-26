
# concatenating VCFs for each chrom
for i in "A1" "A4" "A2" "A3" "X" "Y"
do
ls /ohta2/meng.yuan/rumex/eqtl/VCF/VCF_${i}/*.vcf.gz > VCF_${i}.txt
sort -t. -k5 -n VCF_${i}.txt -o VCF_${i}_sorted.txt

bcftools concat -f VCF_${i}_sorted.txt --threads 20 | \
 bgzip -c > /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}.vcf.gz 

tabix -p vcf /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}.vcf.gz
done

# filter for SNPs
for i in "A1" "A4" "A2" "A3" "X" 
do
# remove one sample with extremely low coverage
bcftools view /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}.vcf.gz \
-s ^NS.1594.001.IDT_i7_141---IDT_i5_141.24fMLD --threads 20 -m2 -M2 -v snps \
| bgzip -c > /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}_SNP.vcf.gz

bcftools stats /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}_SNP.vcf.gz \
> /ohta2/meng.yuan/rumex/eqtl/VCF_filtering/eqtl_mpileup_${i}_SNP.vcf.stats

# get DP field to assess its distribution
bcftools query -f '%INFO/DP\n' /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}_SNP.vcf.gz > \
/ohta2/meng.yuan/rumex/eqtl/VCF_filtering/eqtl_mpileup_${i}_SNP.vcf.DP
done

# keep only male samples for the Y
i="Y"
bcftools view /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}.vcf.gz \
-S male.txt --threads 20 -m2 -M2 -v snps \
| bgzip -c > /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}_SNP.vcf.gz

bcftools stats /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}_SNP.vcf.gz \
> /ohta2/meng.yuan/rumex/eqtl/VCF_filtering/eqtl_mpileup_${i}_SNP.vcf.stats

bcftools query -f '%POS\t%INFO/DP\n' /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}_SNP.vcf.gz > \
/ohta2/meng.yuan/rumex/eqtl/VCF_filtering/eqtl_mpileup_${i}_SNP.vcf.DP
# can separate Y and PAR for male samples later

# keep all samples for PAR (on the Y)
i="Y"
# remove one sample with extremely low coverage
bcftools view /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_${i}.vcf.gz \
-s ^NS.1594.001.IDT_i7_141---IDT_i5_141.24fMLD -r Y:1-45000000 --threads 20 -m2 -M2 -v snps \
| bgzip -c > /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR_SNP.vcf.gz

bcftools stats /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR_SNP.vcf.gz \
> /ohta2/meng.yuan/rumex/eqtl/VCF_filtering/eqtl_mpileup_PAR_SNP.vcf.stats

# get DP field to assess its distribution
bcftools query -f '%INFO/DP\n' /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR_SNP.vcf.gz > \
/ohta2/meng.yuan/rumex/eqtl/VCF_filtering/eqtl_mpileup_PAR_SNP.vcf.DP



