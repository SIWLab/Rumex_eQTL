# pop structure, PCA
# LD pruning then perform PCA

cd /ohta2/meng.yuan/rumex/eqtl/plink
for i in "ML" "FL" "MP" "L"
do 	
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.${i}.filt.vcf.gz

# filtering
# HWE MAF 5%, missing lower 20%
plink --vcf $VCF --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out ${i}

# LD pruning
plink --bfile ${i} --indep-pairwise 200 50 0.1 \
--allow-extra-chr --keep-allele-order --out ${i}_ld

# pca using plink
plink --bfile ${i} --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ${i}_ld.prune.in --pca --out ${i}_ld

# get LD pruned files
plink --bfile ${i} --extract ${i}_ld.prune.in \
--allow-extra-chr --keep-allele-order --make-bed --out ${i}_ld
done



# manualy remove outlier sampled identified from PCA
for i in "ML" "MP" "L"
do 
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.${i}.filt.vcf.gz
VCF2=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.${i}.filt2.vcf.gz
bcftools view ${VCF} -s ^7bM,27eM,53b,5aM --threads 20 | bgzip -c > ${VCF2}
done


# rerun LD pruning
for i in "ML" "MP" "L"
do 	
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.${i}.filt2.vcf.gz

# filtering
# HWE MAF 5%, missing lower 20%
plink --vcf $VCF --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out ${i}

# LD pruning
plink --bfile ${i} --indep-pairwise 200 50 0.1 \
--allow-extra-chr --keep-allele-order --out ${i}_ld

# pca using plink
plink --bfile ${i} --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ${i}_ld.prune.in --pca --out ${i}_ld

# get LD pruned files
plink --bfile ${i} --extract ${i}_ld.prune.in \
--allow-extra-chr --keep-allele-order --make-bed --out ${i}_ld
done


# in EIGENSOFT, use smartpca to do PCA and TW test
for i in "ML" "FL" "MP" "L"
do 
# not marking samples as missing	
sed -i 's/\-9/1/g' ${i}_ld.fam
smartpca -p ${i}_parfile > ${i}_parfile_pca.log
done

## parfile
# genotypename: ML_ld.bed         # The input genotype file
# snpname: ML_ld.bim               # The input SNP file
# indivname: ML_ld.fam             # The input individual file
# evecoutname: ML_ld.evec          # Output file for eigenvectors (PCA results)
# evaloutname: ML_ld.eval          # Output file for eigenvalues
# twstatsoutname: ML_ld.twstats    # Output file for Tracy–Widom statistics
# numoutevec: 20                    # Number of principal components to output



# related matrix
for i in "ML" "FL" "MP" "L"
do 
plink --bfile ${i} --make-rel triangle --allow-extra-chr --out ${i}
plink --bfile ${i} -make-grm-gz no-gz --allow-extra-chr --out ${i}
done
