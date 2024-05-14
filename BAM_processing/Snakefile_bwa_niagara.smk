samples = ['TxF6f6', 'TxF6m23']
ref = "index/merged_TX_noMatPAR"
### workflow ###

with open("../sample_test.txt", "r") as file:
    samples = [line.strip() for line in file if line.strip()]

# This rule establishes the names of output files from other rules 
rule all:
     input:
          bwaidx = "../genome",
          bam = expand("DNAmapped/{sample}.bam", sample = samples),
          markdup = expand("DNAmapped/{sample}_merged_sorted_dupsMarked.bam", sample = samples),
          stats = expand("DNAmapped/duplication_stats/{sample}_dupStats.txt", sample = samples),
          bamidx = expand("DNAmapped/{sample}_merged_sorted_dupsMarked.bam.bai", sample = samples),

rule bwa_mem2:
    """
    Map pe reads to the reference genome using bwa-mem2. Output as BAM
    """
    input:
        r1 = "DNAreads/{sample}_R1.fastq.gz",
        r2 = "DNAreads/{sample}_R2.fastq.gz",
        bwaidx = "../genome"
    output:
        temp('../DNA_bwa/{sample}.bam')
    params:
        r"-R '@RG\tID:{sample}\tSM:{sample}'"
    log: 'logs/bwa_mem2/{sample}.pe.log'
    #threads: ????
    #resources:
    #    tmpdir=""
    shell:
        """
        module load StdEnv/2023
         ( bwa-mem2 mem -t {threads} {input.bwaidx} {input.r1} {input.r2} {params} |\
            samtools view -hb -o {output} - ) 2> {log}
        """

# This rule runs samtools fixmate, samtools sort, and samtools markduplicates

rule samtools_markdup:
    """
    Mark duplicate reads using samtools. Output sorted BAM.
    """
    input:
       'DNAmapped/{sample}.bam' 
       #rules.bwa_mem2.output
    output:
        bam = 'DNAmapped/{sample}_merged_sorted_dupsMarked.bam',
        stats = 'DNAmapped/duplication_stats/{sample}_dupStats.txt'
    log: 'logs/samtools_markdup/{sample}_samtools_markdup.log'
    threads: 8
    resources:
        tmpdir="/ohta2/bianca.sacchi/tmp"
    shell:
        """
        ( samtools fixmate --threads {threads} -m {input} - |\
            samtools sort --threads {threads} -T {wildcards.sample} -o - |\
            samtools markdup --threads {threads} -T {wildcards.sample} -f {output.stats} - {output.bam} ) 2> {log}
        """
# This rule indexes final bams
rule index_bam:
    """
    Index sorted BAM with marked duplicates
    """
    input:
        rules.samtools_markdup.output.bam
    output:
         bamidx = 'DNAmapped/{sample}_merged_sorted_dupsMarked.bam.bai',
    log: 'logs/index_bam/{sample}_index_bam.log'
    threads: 8
    resources:
          tmpdir="/ohta2/bianca.sacchi/tmp"
    shell:
        """
        samtools index -@ {threads} {input} 2> {log}
        """
