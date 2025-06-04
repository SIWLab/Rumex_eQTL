# calculate pi and tajD
# use 4 fold site only

# pixy
# use 4 fold site only
gff=/ohta2/bianca.sacchi/annotation_TX/noMatPAR_txanno/merged_TX_noMatPARlarge_txanno.gff
vcf=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.allsites.L.vcf.gz
ref=/ohta2/Rumex/Dovetail_XY_2023_TX_male/final_scaffolded_assembly/merged_TX_noMatPAR.fa
python3 /ohta1/meng.yuan/apps/genomics_general/codingSiteTypes.py \
-a ${gff} -f gff3 \
-v ${vcf} \
-o noMatPAR_txanno.siteTypes \
-r ${ref} --ignoreConflicts

# separate syn and nonsyn
python3 siteType.py

# get syn and nonsyn code
cut -f 1,2 syn.sites > syn.txt
sed -i 's/A//g' syn.txt

cut -f 1,2 nonsyn.sites > nonsyn.txt
sed -i 's/A//g' nonsyn.txt


vcf=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.allsites.L.vcf.gz
bcftools view -R syn.txt ${vcf} | bgzip -c > eqtl_mpileup_auto.syn.L.vcf.gz

vcf=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.allsites.L.vcf.gz
bcftools view -R nonsyn.txt ${vcf} | bgzip -c > eqtl_mpileup_auto.nonsyn.L.vcf.gz

tabix eqtl_mpileup_auto.syn.L.vcf.gz


# pixy
conda activate pixy
# vcf=/ohta2/meng.yuan/rumex/eqtl/VCF/eqtl_mpileup_auto.allsites.L.vcf.gz
# pixy --stats pi \
# --vcf ${vcf} \
# --populations pop.txt \
# --window_size 10000 \
# --n_cores 25 \
# --output_prefix allsites


vcf=/ohta2/meng.yuan/rumex/eqtl/pixy/eqtl_mpileup_auto.syn.L.vcf.gz
pixy --stats pi \
--vcf ${vcf} \
--populations pop.txt \
--window_size 10000 \
--n_cores 25 \
--output_prefix syn


cut -f 1,2,3 merged_TX_noMatPARlarge_txanno_gene_full.bed > merged_TX_noMatPARlarge_txanno_gene.bed
file=merged_TX_noMatPARlarge_txanno_gene.bed
sed -i 's/A1/1/g' ${file}
sed -i 's/A2/2/g' ${file}
sed -i 's/A3/3/g' ${file}
sed -i 's/A4/4/g' ${file}

# gene by gene stats
vcf=/ohta2/meng.yuan/rumex/eqtl/pixy/eqtl_mpileup_auto.syn.L.vcf.gz
bed=merged_TX_noMatPARlarge_txanno_gene.bed
pixy --stats pi tajima_d \
--vcf ${vcf} \
--populations pop.txt \
--n_cores 6 \
--bed_file ${bed} \
--output_prefix syn.genewise \
--output_folder syn_genewise
