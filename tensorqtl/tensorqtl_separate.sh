# genotype data for males

vcftools --vcf my_file.vcf --freq --out my_freq_out_file 

# HWE MAF 5%, missing lower 20%
#
for i in "ML" "FL" "MP"
do 
VCF=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.${i}.filt.vcf.gz

plink --vcf $VCF --double-id --allow-extra-chr  \
--keep-allele-order --set-missing-var-ids @:# \
--maf 0.05 --geno 0.2 \
--hwe 1e-6 --biallelic-only strict \
--make-bed --out ${i}

done

bcftools view /ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.SNP.ML.filt.vcf.gz \
-r A1:274220346 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' > A1_274220346.vcf


# phenotype data
grep 'maker\stranscript' merged_TX_noMatPARlarge_txanno.gtf > merged_TX_noMatPARlarge_txanno_gene_full.gtf
cut -f 1,4,5,7  merged_TX_noMatPARlarge_txanno_gene_full.gtf > merged_TX_noMatPARlarge_txanno_gene_full.bed

# already sorted 
bgzip normalized_counts_mln.bed && tabix -p bed normalized_counts_mln.bed.gz
bgzip normalized_counts_fln.bed && tabix -p bed normalized_counts_fln.bed.gz
bgzip normalized_counts_mpn.bed && tabix -p bed normalized_counts_mpn.bed.gz



# male leaf
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/ML
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mln.bed.gz
prefix=ML
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_ML_pc.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000


# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000

  # * 75 samples
  # * 14804 phenotypes
  # * 0 covariates
  # * 3664501 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 14804/14804
  #   ** dropping 146 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome A1
  #   processing phenotype 4973/14658    time elapsed: 0.40 min
  #   * writing output
  #   Mapping chromosome A2
  #   processing phenotype 9631/14658    time elapsed: 0.72 min
  #   * writing output
  #   Mapping chromosome A3
  #   processing phenotype 12331/14658    time elapsed: 0.90 min
  #   * writing output
  #   Mapping chromosome A4
  #   processing phenotype 14658/14658
  #   time elapsed: 1.06 min



# female leaf
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/FL
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_fln.bed.gz
prefix=FL
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_FL_pc.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000

  # * 74 samples
  # * 14804 phenotypes
  # * 0 covariates
  # * 4509904 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 14804/14804
  #   ** dropping 71 phenotypes without variants in cis-window
  # * computing permutations
  #   processing phenotype 14733/14733

# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000

  # * 74 samples
  # * 14804 phenotypes
  # * 0 covariates
  # * 4509904 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 14804/14804
  #   ** dropping 71 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome A1
  #   processing phenotype 4981/14733    time elapsed: 0.46 min
  #   * writing output
  #   Mapping chromosome A2
  #   processing phenotype 9668/14733    time elapsed: 0.84 min
  #   * writing output
  #   Mapping chromosome A3
  #   processing phenotype 12375/14733    time elapsed: 1.07 min
  #   * writing output
  #   Mapping chromosome A4
  #   processing phenotype 14733/14733
  #   time elapsed: 1.25 min




# cp covariate.txt covariate_MP.txt
# pollen
plink_prefix_path=/ohta2/meng.yuan/rumex/eqtl/plink/MP
expression_bed=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/normalized_counts_mpn.bed.gz
prefix=MP
covariates_file=/ohta2/meng.yuan/rumex/eqtl/tensorqtl/covariate_MP_pc.txt

# cis-QTL mapping: permutations
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis --window 20000



# cis-QTL mapping: summary statistics for all variant-phenotype pairs
python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
    --covariates ${covariates_file} \
    --mode cis_nominal --window 20000

  # * 73 samples
  # * 13916 phenotypes
  # * 0 covariates
  # * 3464540 variants
  # * cis-window: ±20,000
  # * checking phenotypes: 13916/13916
  #   ** dropping 122 phenotypes without variants in cis-window
  # * Computing associations
  #   Mapping chromosome A1
  #   processing phenotype 4621/13794    time elapsed: 0.35 min
  #   * writing output
  #   Mapping chromosome A2
  #   processing phenotype 9044/13794    time elapsed: 0.64 min
  #   * writing output
  #   Mapping chromosome A3
  #   processing phenotype 11572/13794    time elapsed: 0.80 min
  #   * writing output
  #   Mapping chromosome A4
  #   processing phenotype 13794/13794




