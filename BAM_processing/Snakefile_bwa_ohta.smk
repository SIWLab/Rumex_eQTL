### workflow ###

with open("DNA.samples_full.txt", "r") as file:
    samples = [line.strip() for line in file if line.strip()]

# This rule establishes the names of output files from other rules 
rule all:
     input:
          bwaidx = "index/merged_TX_noMatPAR",
          bam = expand("AnalysisReadyBams/{sample}.bam", sample = samples),
          markdup = expand("AnalysisReadyBams/{sample}_sortMarkDups.bam", sample = samples),
          bamidx = expand("AnalysisReadyBams/{sample}_sortMarkDups.bam.bai", sample = samples),

rule bwa_index:
    """
    create index
    """
    input:
        database="merged_TX_noMatPAR.fa"
    output:
        done=touch("index/merged_TX_noMatPAR")
    log:
       "logs/bwa_mem2/index.log"
    threads:
        20
    resources:
        tmpdir="/ohta2/bianca.sacchi/tmp/"
    shell:
        "bwa-mem2 index -p index/merged_TX_noMatPAR {input.database}" 


rule bwa_mem2:
    """
    Map pe reads to the reference genome using bwa-mem2. Output as BAM
    """
    input:
        r1 = "/ohta2/Rumex/eQTL/DNA/{sample}_R1.fastq.gz",
        r2 = "/ohta2/Rumex/eQTL/DNA/{sample}_R2.fastq.gz",
        bwaidx = rules.bwa_index.output
    output:
        temp('AnalysisReadyBams/{sample}.bam')
    params:
        r"-R '@RG\tID:{sample}\tSM:{sample}'"
    log: 'logs/bwa_mem2/{sample}.pe.log'
    threads: 20
    resources:
        tmpdir="/ohta2/bianca.sacchi/tmp"
    shell:
        """
         ( bwa-mem2 mem -t {threads} {input.bwaidx} {input.r1} {input.r2} {params} |\
            samtools view -hb -o {output} - ) 2> {log}
        """

# This rule runs samtools fixmate, samtools sort, and samtools markduplicates

rule samtools_markdup:
    """
    Mark duplicate reads using samtools. Output sorted BAM.
    """
    input: rules.bwa_mem2.output
    output:
        bam = 'AnalysisReadyBams/{sample}_sortMarkDups.bam',
        stats = 'duplication_stats/{sample}_dupStats.txt'
    log: 'logs/samtools_markdup/{sample}_samtools_markdup.log'
    threads: 20
    resources:
        tmpdir="/ohta2/bianca.sacchi/tmp"
    shell:
        """
        ( samtools fixmate --threads {threads} -m {input} - |\
            samtools sort --threads {threads} -T {resources.tmpdir}/{wildcards.sample} -o - |\
            samtools markdup --threads {threads} -T {resources.tmpdir}/{wildcards.sample} -f {output.stats} - {output.bam} ) 2> {log}
        """
# This rule indexes final bams
rule index_bam:
    """
    Index sorted BAM with marked duplicates
    """
    input:
        rules.samtools_markdup.output.bam
    output:
         bamidx = 'AnalysisReadyBams/{sample}_sortMarkDups.bam.bai',
    log: 'logs/index_bam/{sample}_index_bam.log'
    threads: 20
    resources:
          tmpdir="/ohta2/bianca.sacchi/tmp"
    shell:
        """
        samtools index -@ {threads} {input} 2> {log}
        """
