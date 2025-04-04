#!/bin/bash
cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/ML
mkdir cnt
mkdir fake
ln -s ../normalized_counts_mln_auto.bed ./
# make sure pheno_permute.R, get_MAF_null.R and ml_permute.sh are in this directory
vi ml_permute.sh
vi pheno_permute.R
vi get_MAF_null.R

# run permutations 1000 times
sh ml_permute.sh 1000


# mv ML_*_random.cnt cnt/
# mv ML_*_eqtl_random_fake.txt fake/

# merge results from all permutations
cat ./fake/ML_*_eqtl_random_fake.txt > ML_MAF_null.txt
cat ./cnt/ML_*_random.cnt > ML_fake_cnt.txt

# rm intermediate files
rm -rf fake
rm -rf cnt

# FL
cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/FL
mkdir cnt
mkdir fake
ln -s ../normalized_counts_fln_auto.bed ./

vi fl_permute.sh
vi pheno_permute.R
vi get_MAF_null.R

# run permutations 1000 times
sh fl_permute.sh 1000


# merge results from all permutations
cat ./fake/FL_*_eqtl_random_fake.txt > FL_MAF_null.txt
cat ./cnt/FL_*_random.cnt > FL_fake_cnt.txt

cat ./fake/FL_*_eqtl_random_fake.txt > FL_MAF_null_mar19.txt
cat ./cnt/FL_*_random.cnt1 > FL_fake_cnt1.txt
cat ./cnt/FL_*_random.cnt2 > FL_fake_cnt2.txt


# rm intermediate files
rm -rf fake
rm -rf cnt


# MP
cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/MP
mkdir cnt
mkdir fake
ln -s ../normalized_counts_mpn_auto.bed ./

vi mp_permute.sh
vi pheno_permute.R
vi get_MAF_null.R

# run permutations 1000 times
sh fl_permute.sh 1000

cat ./fake/MP_*_eqtl_random_fake.txt > MP_MAF_null.txt
cat ./cnt/MP_*_random.cnt > MP_fake_cnt.txt

# rm intermediate files
rm -rf fake
rm -rf cnt