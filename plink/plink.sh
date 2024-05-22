

#for i in "A1" "A2" "A3" "A4" "X" 
# pop structure, PCA
# perform linkage pruning - i.e. identify prune sites
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_A1.SNP.filt.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out A4

# prune and create pca
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract A4.prune.in \
--make-bed --pca --out A4







# LD decay
# calc ld with plink
# A4
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5  \
--thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out A4

python ld_decay_calc.py -i A4.ld.gz -o A4


# A1
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_A1.SNP.filt.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5  \
--thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out A1

python ld_decay_calc.py -i A1.ld.gz -o A1


# A2
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_A2.SNP.filt.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5  \
--thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out A2

python ld_decay_calc.py -i A2.ld.gz -o A2


# A3
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_A3.SNP.filt.vcf.gz
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5  \
--thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out A3

python ld_decay_calc.py -i A3.ld.gz -o A3

# A4
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_A4.SNP.filt.vcf.gz

plink --vcf $VCF --recode rlist


plink --vcf $VCF --ped PEDfile_DNA_eqtl.ped --allow-extra-chr \
--set-missing-var-ids @:# \
--maf 0.01 --geno 0.1 --mind 0.5  \
--thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out A4

python ld_decay_calc.py -i A4.ld.gz -o A4




# HWE