with open("DNA.samples_subset.txt", "r") as file:
    samples = [line.strip() for line in file if line.strip()]
reference = "merged_TX_noMatPAR.fa"

chroms = ["A1,"A2","A3","A4","X","Y"]
nchunks = 9
bamlist = "path/to/bam.list"
chunks = list(range(1,nchunks+1))

rule GenomeIndex:
    input:
        ref = reference
    output:
        idx = reference + ".fai"
    log: 
        "logs/faidx.log"
    wrapper: 
        "v0.69.0/bio/samtools/faidx"

rule GenerateFreebayesRegions:
    input:
        ref_idx = reference,
        index = reference + ".fai",
    output:
        regions = expand("resources/regions/genome.{chrom}.region.{i}.bed", chrom=chroms, i = chunks)
    log:
        "logs/GenerateFreebayesRegions.log"
    params:
        chroms = chroms,
        chunks = chunks
    conda:
        "envs/freebayes_simple.yaml"
    script:
        "~/bin/freebayes_scripts/fasta_generate_regions.py --chunks --bed resources/regions/genome --chromosome {params.chroms} {input.index} {params.chunks}"


rule VariantCallingFreebayes:
    input:
        bams = expand("AnalysisReadyBams/{sample}_sortMarkDups.bam", sample=samples),
        index = expand("AnalysisReadyBams/{sample}_sortMarkDups.bam.bai", sample=samples),
        ref = reference,
        samples = bamlist,
        regions = "resources/regions/genome.{chrom}.region.{i}.bed"
    output:
        temp("results/variants/vcfs/{chrom}/variants.{i}.vcf")
    log:
        "logs/VariantCallingFreebayes/{chrom}.{i}.log"
    conda:
        "../envs/freebayes-env.yaml"
    threads:1
    shell:	"freebayes -f {input.ref} -t {input.regions} -L {input.samples} > {output} 2> {log}"


rule ConcatVCFs:
    input:
        calls = expand("results/variants/vcfs/{{chrom}}/variants.{i}.vcf", i=chunks)
    output:
        "results/variants/vcfs/variants.{chrom}.vcf"
    log:
        "logs/ConcatVCFs/{chrom}.log"
    conda:
        "../envs/freebayes-env.yaml"
    threads:4
    shell:  
        "bcftools concat {input.calls} | vcfuniq > {output} 2> {log}"
