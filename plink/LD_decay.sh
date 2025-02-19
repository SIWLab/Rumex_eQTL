# LD decay
# use /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.L.filt2.vcf.gz instead
#SNPs 1 kb apart were ignored. 10% of SNPs were kept to calculate mean r2 between SNP pairs in 10-bp bins.
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_A4.SNP.filt.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5  \
--thin 0.1 -r2 gz --ld-window-kb 1 \
--ld-window-r2 0 \
--make-bed --out A4


VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_A1.SNP.filt.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5  \
--thin 0.1 -r2 gz --ld-window-kb 1 \
--ld-window-r2 0 \
--make-bed --out A1

VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_A2.SNP.filt.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5  \
--thin 0.1 -r2 gz --ld-window-kb 1 \
--ld-window-r2 0 \
--make-bed --out A2

VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_A3.SNP.filt.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5  \
--thin 0.1 -r2 gz --ld-window-kb 1 \
--ld-window-r2 0 \
--make-bed --out A3

python ld_decay_calc_new.py -i A1.ld.gz -o A1
python ld_decay_calc_new.py -i A2.ld.gz -o A2
python ld_decay_calc_new.py -i A3.ld.gz -o A3
python ld_decay_calc_new.py -i A4.ld.gz -o A4

