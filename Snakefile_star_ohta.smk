from os import path

import numpy as np

import pandas as pd


#with open("RNA.samples_leaf.txt", "r") as file:
#    samples = [line.strip() for line in file if line.strip()]

with open("RNA.samples_pollen.txt", "r") as file:
     samples = [line.strip() for line in file if line.strip()]

rule all:
     input:
          idx = "STARindex/merged_TX_noMatPAR",
          aln = expand("RNAbams/{sample}Aligned.sortedByCoord.out.bam",sample = samples),
          bamidx = expand("RNAbams/{sample}Aligned.sortedByCoord.out.bam.bai", sample = samples),

rule star_index:
    input:
        fasta="merged_TX_noMatPAR.fa",
    output:
        directory("STARindex/merged_TX_noMatPAR"),
    threads: 20
    params:
        extra="",
    resources:
        tmpdir="/ohta2/bianca.sacchi/tmp"
    wrapper:
        "v2.6.1/bio/star/index"

rule star_pe_multi:
    input:
        r1 = "/ohta2/Rumex/eQTL/Transcripts/{sample}_R1.fastq.gz",
        r2 = "/ohta2/Rumex/eQTL/Transcripts/{sample}_R2.fastq.gz",
        idx = rules.star_index.output
    output:
        aln="RNAbams/{sample}Aligned.sortedByCoord.out.bam",
        log="RNAbams/{sample}Log.out",
        sj="RNAbams/{sample}SJ.out.tab",
        log_final="RNAbams/{sample}Log.final.out",
    params:
        prefix = "RNAbams/{sample}",
    threads: 4
    resources:
        tmpdir="/ohta2/bianca.sacchi/tmp"
    shell:
         """
         STAR --runThreadN {threads} \
         --readFilesCommand zcat \
         --genomeDir {input.idx} \
         --readFilesIn {input.r1} {input.r2} \
         --outSAMtype BAM SortedByCoordinate \
         --twopassMode Basic \
         --outSAMattrRGline ID:{wildcards.sample} \
         --outFileNamePrefix {params.prefix} \
         --outTmpDir {resources.tmpdir}/{wildcards.sample} \
         """

rule index_bam:
    """
    Index sorted BAM with marked duplicates
    """
    input:
        aln={rules.star_pe_multi.output.aln}
    output:
        bamidx = 'RNAbams/{sample}Aligned.sortedByCoord.out.bam.bai',
    log: 'logs/index_bam/{sample}_index_bam.rna.fullmerge.log'
    threads: 10
    resources:
          tmpdir="/ohta2/bianca.sacchi/tmp"
    shell:
        """
        samtools index -@ {threads} {input} 2> {log}
        """
