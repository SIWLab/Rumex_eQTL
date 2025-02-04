# keep only cis SNPs as input
data <- "/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/"
ML_A4 <- read_parquet(paste0(data,"ML.cis_qtl_pairs.4.parquet"))

data <- "/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/"
normalized_counts_mln <- read.table(paste0(data,"normalized_counts_mln_auto.bed"), header = T, check.names = F)
gff <- read.table(paste0(data,"merged_TX_noMatPARlarge_txanno_genes.gtf"))
gff <- gff[,c(13,1,4,5)]
colnames(gff) <- c("gene","chrom", "start", "end")
gff <- gff %>% filter(chrom == "A4") %>% select(gene)
colnames(gff) <- "phenotype_id"

length(unique(ML_A4$phenotype_id)) # 2332
length(unique(ML_A4$variant_id)) # 115863

ML <- read.table("/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/ML.cis_qtl.txt", header = T) 
ML <- inner_join(gff, ML, by = "phenotype_id") # 2332
ML_egene <- ML %>% filter(qval < 0.1) # 390

ML_eqtl_full <- inner_join(ML_eqtl, ML_egene[,c(1,16:19)], by = "phenotype_id") # 36592
ML_eqtl_sig <- ML_eqtl_full %>% filter(pval_nominal < pval_nominal_threshold) # 2694
length(unique(ML_eqtl_sig$variant_id)) # 2683
length(unique(ML_eqtl_sig$phenotype_id)) # 389
ML_eqtl_sig2 <- ML_eqtl_full %>% filter(pval_nominal < 0.0318507) # 7947
length(unique(ML_eqtl_sig2$variant_id)) # 7779
length(unique(ML_eqtl_sig2$phenotype_id)) # 390


t1 <- as.data.frame(unique(ML_eqtl_sig$phenotype_id))
t2 <- as.data.frame(unique(ML_egene$phenotype_id))
colnames(t1) <- "phenotype_id"
colnames(t2) <- "phenotype_id"
anti_join(t2, t1, by = "phenotype_id")

# pheno
ML_A4 <- ML_A4[1:1000,]
gene <- ML_A4 %>% select(phenotype_id)
gene <-  gene[!duplicated(gene), ]
colnames(gene) <- "geneid"
normalized_counts_mln <- normalized_counts_mln[,-c(1:3)]
colnames(normalized_counts_mln)[1]  <- "geneid"
normalized_counts_mln <- inner_join(gene, normalized_counts_mln, by = "geneid")
write.table(normalized_counts_mln, file = paste0(data,"normalized_counts_mln_test.bed"), row.names = FALSE, quote = F, sep = "\t")

normalized_counts_mln <- read.table(paste0(data, "normalized_counts_mln_test.bed"), header = T, )
gene <- normalized_counts_mln %>% select(geneid)
colnames(gff)[1] <- "geneid"
gene <- inner_join(gene, gff, by = "geneid")
colnames(gene)[2:4] <- c("chr","left", "right")
gene$chr <- factor(gene$chr)
levels(gene$chr) <- "4"
write.table(gene, file = paste0(data,"gene_pos.txt"), row.names = FALSE, quote = F, sep = "\t")

ML_A4 <- ML_A4 %>% select(variant_id)
write.table(ML_A4, file = paste0(data,"A4_variant_test.bed"), row.names = FALSE, quote = F, sep = "\t", col.names = F)

# geno
genotypes <- read.table("A4_genotypes.txt", header=T, stringsAsFactors=FALSE, check.names = F)
genotype_matrix <- genotypes[, -c(1:4)]  # Exclude chromosome, position, ID columns
genotype_matrix[genotype_matrix == "0/0"] <- 0
genotype_matrix[genotype_matrix == "0/1"] <- 1
genotype_matrix[genotype_matrix == "1/1"] <- 2
genotype_matrix[genotype_matrix == "./."] <- NA
genotype_matrix <- apply(genotype_matrix, 2, as.numeric)  # Convert all to numeric

rownames(genotype_matrix) <- paste0(genotypes[,1], "_", genotypes[,2])  # CHROM:POS
write.table(genotype_matrix, "genotype_matrix.txt", sep="\t", quote=FALSE, col.names=TRUE)
genotypes2 <- read.table("genotype_matrix.txt", header=T, stringsAsFactors=FALSE, check.names = F)
snp <- cbind(genotypes2[,1], genotypes[,1], genotypes[,2])
colnames(snp) <- c("snpid", "chr", "pos")
write.table(snp, file = "snp_pos.txt", row.names = FALSE, quote = F, sep = "\t")

############# run matrixeQTL ####################
# subset of A4
covariate_ML.txt
genotype_matrix.txt
normalized_counts_mln_test.bed

setwd("/ohta2/meng.yuan/rumex/eqtl/matrixeqtl")

# First step is to load the package:
library("MatrixEQTL")

# The toy data set files are stored with the package at the following location.
#base.dir = find.package("MatrixEQTL")

# Then we set the parameters such as selected linear model and
# names of genotype and expression data files.

useModel = modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS
# SNP_file_name = paste0(base.dir, "/data/SNP.txt")
# expression_file_name = paste0(base.dir, "/data/GE.txt")
SNP_file_name = "genotype_matrix.txt"
expression_file_name = "normalized_counts_mln_test.bed"

# A separate file may be provided with extra covariates.
# In case of no covariates set the variable covariates_file_name to character().

# covariates_file_name = paste0(base.dir, "/data/Covariates.txt")
covariates_file_name = "covariate_ML.txt"
output_file_name = tempfile()

# The p-value threshold determines which gene-SNP associations are
# saved in the output file output_file_name.
# Note that for larger datasets the threshold should be lower.
# Setting the threshold to a high value for a large dataset may
# cause excessively large output files.

pvOutputThreshold = 1e-2

# Finally, define the covariance matrix for the error term.
# This parameter is rarely used.
# If the covariance matrix is a multiple of identity, set it to numeric().

errorCovariance = numeric()

# The next section of the sample code contains three very similar parts
# loading the files with genotype, gene expression, and covariates.
# In each part one can set the file delimiter
# (i.e. tabulation "\t", comma ",", or space " "),
# the string representation for missing values,
# the number of rows with column labels, and
# the number of columns with row labels.
# Finally, one can change the number of the variables
# in a slice for the file reading procedure (do not change if not sure).

snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name )

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name )

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
cvrt$LoadFile( covariates_file_name )


genepos = SlicedData$new()
genepos$fileDelimiter = "\t"      # the TAB character
genepos$fileOmitCharacters = "NA" # denote missing values
genepos$fileSkipRows = 1          # one row of column labels
genepos$fileSkipColumns = 1       # one column of row labels
genepos$fileSliceSize = 2000      # read file in pieces of 2,000 rows
genepos$LoadFile("gene_pos.txt")

snpspos = SlicedData$new()
snpspos$fileDelimiter = "\t"      # the TAB character
snpspos$fileOmitCharacters = "NA" # denote missing values
snpspos$fileSkipRows = 1          # one row of column labels
snpspos$fileSkipColumns = 1       # one column of row labels
snpspos$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snpspos$LoadFile("snp_pos.txt")

snpspos <- read.table("snp_pos.txt", header = T)
genepos <- read.table("gene_pos.txt", header = T)

# Finally, the main Matrix eQTL function is called:

me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)



me2 = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE,
    pvOutputThreshold.cis = 0.1,
    output_file_name.cis = "cis_eqtl.txt",
    cisDist = 20000,
    snpspos = snpspos,
    genepos = genepos)


 
100.00% done, 14,729,315 eQTLs
100.00% done, 3,113,414 eQTLs

Expected number of findings > 13509625.8
Expected number of findings > 2701925.16

save(me, file="matrixeqtl-batched.rda")

me2$all$ntests # 270192516, 270192516
me2$all$ntests # 270192516

eqtl_cis = me2$cis$eqtls
min(eqtl_cis$FDR) # 0.3517423

eqtl_trans = me2$trans$eqtls
min(eqtl_trans$FDR) # 0.8151034


eqtls = dplyr::filter(me$all$eqtls, FDR < 0.1) ##select significant eQTLs

eqtl_all = me$all$eqtls
min(eqtl_all$FDR) #0.8149853



write.table(me$all$eqtls, file = "eqtl-output.txt", row.names=F, col.names=T, quote=F)



cis <- read.table("/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/cis_eqtl.txt", header = T)
ggplot(cis, aes(x=p.value)) + geom_histogram(color="black", fill="white", binwidth = 0.005) + theme_classic() + xlab("p-value")+ggtitle("matrixeqtl")

176878 cis_eqtl.txt
> nrow(eqtl_cis)
[1] 176877
> nrow(eqtl_trans)
[1] 3111437
> me2$all$ntests
[1] 270192516
