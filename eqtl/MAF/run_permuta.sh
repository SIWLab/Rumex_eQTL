#!/bin/bash
########################### ML & MP ##########################
cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/Male
ln -s ../normalized_counts_mln_auto.bed ./
ln -s ../normalized_counts_mpn_auto.bed ./
ln -s ../ML.cis_qtl.txt ./
ln -s ../MP.cis_qtl.txt ./
mkdir cnt
mkdir fake

vi pheno_permute_m.R
vi get_MAF_null_m.R
vi male_permute.sh

# run permutations 1000 times
sh male_permute.sh 1000

# save output from all permutations
cat ./fake/ML_*_eqtl_random_fake.txt > ML_MAF_null.txt
cat ./cnt/ML_*_random.cnt > ML_fake_cnt.txt

cat ./fake/MP_*_eqtl_random_fake.txt > MP_MAF_null.txt
cat ./cnt/MP_*_random.cnt > MP_fake_cnt.txt

cat ./fake/tissue_n1_*_fake.txt > tissue_n1_fake.txt
cat ./fake/tissue_n2_*_fake.txt > tissue_n2_fake.txt
cat ./cnt/tissue_n1_n2_*_random.cnt > tissue_n1_n2_random.cnt

# rm intermediate files
rm -rf fake
rm -rf cnt

########################### FL ##########################
# FL autosomes
cd /ohta2/meng.yuan/rumex/eqtl/tensorqtl/FL
mkdir cnt
mkdir fake
ln -s ../normalized_counts_fln_auto.bed ./
ln -s ../FL.cis_qtl.txt ./

vi fl_permute.sh
vi pheno_permute_fl.R
vi get_MAF_null_fl.R

# run permutations 1000 times
sh fl_permute.sh 1000

# merge results from all permutations
cat ./fake/FL_*_eqtl_random_fake.txt > FL_MAF_null.txt
cat ./cnt/FL_*_random.cnt > FL_fake_cnt.txt

# rm intermediate files
rm -rf fake
rm -rf cnt

####################### FL X & PAR #######################
# one permutation as well



