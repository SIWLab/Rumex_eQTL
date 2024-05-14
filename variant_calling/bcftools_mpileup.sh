#each chrom has its bed folder

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=5:00:00
#SBATCH --job-name mpileup_A1
module load vcftools
module load samtools
module load NiaEnv/2019b gnu-parallel
ls *.bed | parallel -j 36 'bcftools mpileup -a DP,AD -B -I -f /scratch/w/wrighste/yuanmeng/genome/merged_TX_noMatPAR.fa -R {} -b /scratch/w/wrighste/yuanmeng/bam_list.txt | bcftools call -m -Oz -a GQ  -o /scratch/w/wrighste/yuanmeng/VCF/{}.vcf.gz >> /scratch/w/wrighste/yuanmeng/out/{}.out 2>&1'



#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=7:00:00
#SBATCH --job-name mpileup_A2
module load vcftools
module load samtools
module load NiaEnv/2019b gnu-parallel

ls *.bed | parallel --joblog slurm-$SLURM_JOBID.log -j 40 'bcftools mpileup -a DP,AD -B -I -f /scratch/w/wrighste/yuanmeng/genome/merged_TX_noMatPAR.fa -R {} -b /scratch/w/wrighste/yuanmeng/bam_list.txt | bcftools call -m -Oz -a GQ  -o /scratch/w/wrighste/yuanmeng/VCF_A2/{}.vcf.gz >> /scratch/w/wrighste/yuanmeng/out/{}.out 2>&1'



#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=5:00:00
#SBATCH --job-name mpileup_A3
module load vcftools
module load samtools
module load NiaEnv/2019b gnu-parallel

ls *.bed | parallel --joblog log/slurm-$SLURM_JOBID.log -j 40 'bcftools mpileup -a DP,AD -B -I -f /scratch/w/wrighste/yuanmeng/genome/merged_TX_noMatPAR.fa -R {} -b /scratch/w/wrighste/yuanmeng/bam_list.txt | bcftools call -m -Oz -a GQ  -o /scratch/w/wrighste/yuanmeng/VCF_A3/{}.vcf.gz >> /scratch/w/wrighste/yuanmeng/out/{}.out 2>&1'




#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=5:00:00
#SBATCH --job-name mpileup_A4
module load vcftools
module load samtools
module load NiaEnv/2019b gnu-parallel

ls *.bed | parallel --joblog log/slurm-$SLURM_JOBID.log -j 40 'bcftools mpileup -a DP,AD -B -I -f /scratch/w/wrighste/yuanmeng/genome/merged_TX_noMatPAR.fa -R {} -b /scratch/w/wrighste/yuanmeng/bam_list.txt | bcftools call -m -Oz -a GQ  -o /scratch/w/wrighste/yuanmeng/VCF_A4/{}.vcf.gz >> /scratch/w/wrighste/yuanmeng/out/{}.out 2>&1'



#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=5:00:00
#SBATCH --job-name mpileup_X
module load vcftools
module load samtools
module load NiaEnv/2019b gnu-parallel

ls *.bed | parallel --joblog log/slurm-$SLURM_JOBID.log -j 40 'bcftools mpileup -a DP,AD -B -I -f /scratch/w/wrighste/yuanmeng/genome/merged_TX_noMatPAR.fa -R {} -b /scratch/w/wrighste/yuanmeng/bam_list.txt | bcftools call -m -Oz -a GQ  -o /scratch/w/wrighste/yuanmeng/VCF_X/{}.vcf.gz >> /scratch/w/wrighste/yuanmeng/out/{}.out 2>&1'



#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=8:30:00
#SBATCH --job-name mpileup_Y
module load vcftools
module load samtools
module load NiaEnv/2019b gnu-parallel

ls *.bed | parallel --joblog log/slurm-$SLURM_JOBID.log -j 40 'bcftools mpileup -a DP,AD -B -I -f /scratch/w/wrighste/yuanmeng/genome/merged_TX_noMatPAR.fa -R {} -b /scratch/w/wrighste/yuanmeng/bam_list.txt | bcftools call -m -Oz -a GQ  -o /scratch/w/wrighste/yuanmeng/VCF_Y/{}.vcf.gz >> /scratch/w/wrighste/yuanmeng/out/{}.out 2>&1'


