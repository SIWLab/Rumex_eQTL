# Mapping eQTLs in *Rumex hastatulus* 
This repository contains scripts for the *Rumex* eQTL mapping project, the manuscript is in preparation.

---

Here's a brief description of the folders containing scripts for different analyses:

* `BAM_processing`: RNA and WGS raw reads alignment, BAM preprocessing.
* `variant_calling`: variant calling and SNP filtering.
* `expression`: counting RNA reads from BAMs, differential gene expression analyses, and preparing phenotype files for eQTL mapping.
* `plink`: assessing LD decay, PCA, and performing Tracy-Widom test.
* `tensorqtl`: running TensorQTL for cis-eQTL mapping.
* `eqtl`: identifying eQTLs, analyzing allele frequencies, effect sizes etc.
* `popgen`: generating neutral diversity statistics