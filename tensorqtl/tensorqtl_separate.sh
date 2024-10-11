# genotype data
# for autosomes, genotype files are already generated
# genotype data for X and PAR for females
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_X.SNP.FL.filt.vcf.gz
i=FL_X
plink --vcf ${VCF} --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out ${i}

VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_PAR.SNP.FL.filt.vcf.gz
i=FL_PAR
plink --vcf ${VCF} --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out ${i}


# phenotype data
grep 'maker\stranscript' merged_TX_noMatPARlarge_txanno.gtf > merged_TX_noMatPARlarge_txanno_gene_full.gtf
cut -f 1,4,5,7  merged_TX_noMatPARlarge_txanno_gene_full.gtf > merged_TX_noMatPARlarge_txanno_gene_full.bed

# already sorted 
bgzip normalized_counts_mln.bed && tabix -p bed normalized_counts_mln.bed.gz
bgzip normalized_counts_fln_auto.bed && tabix -p bed normalized_counts_fln_auto.bed.gz
bgzip normalized_counts_mpn.bed && tabix -p bed normalized_counts_mpn.bed.gz

# rename chrom for PAR and X so there's no issue with plink
sed -i 's/X/5/g' normalized_counts_fln_X.bed
sed -i 's/PAR/6/g' normalized_counts_fln_PAR.bed
bgzip normalized_counts_fln_X.bed && tabix -p bed normalized_counts_fln_X.bed.gz
bgzip normalized_counts_fln_PAR.bed && tabix -p bed normalized_counts_fln_PAR.bed.gz


# male leaf
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/ML
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mln.bed.gz
prefix=ML
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000


# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000


# female leaf
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_fln_auto.bed.gz
prefix=FL
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000


# X and PAR for female
# X
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL_X
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_fln_X.bed.gz
prefix=FL_X
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000


# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000

# PAR
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL_PAR
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_fln_PAR.bed.gz
prefix=FL_PAR
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000
  # * 74 samples
  # * 1102 phenotypes
  # * 0 covariates
  # * 96911 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 1102/1102
  #   ** dropping 12 phenotypes without variants in cis-window
  # * computing permutations
  #   processing phenotype 1090/1090

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000
  # * 74 samples
  # * 1102 phenotypes
  # * 0 covariates
  # * 96911 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 1102/1102
  #   ** dropping 12 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome 1
  #   processing phenotype 1090/1090


# cp covariate.txt covariate_MP.txt
# pollen
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mpn.bed.gz
prefix=MP
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_MP.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000


# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000





