# generate the bed files
samtools faidx merged_TX_noMatPAR.fa

# 2Mb windows per chromosome
python3  /ohta2/meng.yuan/apps/freebayes-1.3.6/scripts/fasta_generate_regions.py \
../merged_TX_noMatPAR.fa 2000000 --chromosomes A1 --bed chr

python3  /ohta2/meng.yuan/apps/freebayes-1.3.6/scripts/fasta_generate_regions.py \
../merged_TX_noMatPAR.fa 2000000 --chromosomes A2 --bed chr

python3  /ohta2/meng.yuan/apps/freebayes-1.3.6/scripts/fasta_generate_regions.py \
../merged_TX_noMatPAR.fa 2000000 --chromosomes A3 --bed chr

python3  /ohta2/meng.yuan/apps/freebayes-1.3.6/scripts/fasta_generate_regions.py \
../merged_TX_noMatPAR.fa 2000000 --chromosomes A4 --bed chr

python3  /ohta2/meng.yuan/apps/freebayes-1.3.6/scripts/fasta_generate_regions.py \
../merged_TX_noMatPAR.fa 2000000 --chromosomes X --bed chr

python3  /ohta2/meng.yuan/apps/freebayes-1.3.6/scripts/fasta_generate_regions.py \
../merged_TX_noMatPAR.fa 2000000 --chromosomes Y --bed chr

