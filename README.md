# Mapping eQTLs in *Rumex hastatulus* 
This repository contains scripts for the *Rumex* eQTL mapping project, the manuscript is currently under review, the preprint is available here:

>Yuan M, Sacchi BM, Choudhury BI, Barrett SC, Stinchcombe JR, Wright SI. 2025. Cis-regulation of gene expression between sexes and life stages in Rumex hastatulus. [bioRxiv. 2025:2025.06.16.659834](https://www.biorxiv.org/content/10.1101/2025.06.16.659834v1).

---

Here's a brief description of the folders containing scripts for different analyses:

* `BAM_processing`: RNA and WGS raw reads alignment, BAM preprocessing.
* `variant_calling`: variant calling and SNP filtering.
* `expression`: counting RNA reads from BAMs, differential gene expression analyses, and preparing phenotype files for eQTL mapping.
* `plink`: assessing LD decay, PCA, and performing Tracy-Widom test.
* `tensorqtl`: running TensorQTL for cis-eQTL mapping.
* `eqtl`: identifying eQTLs, analyzing allele frequencies, effect sizes etc.
* `popgen`: generating neutral diversity statistics
