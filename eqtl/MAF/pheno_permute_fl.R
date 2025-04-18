setwd("/ohta2/meng.yuan/rumex/eqtl/tensorqtl/FL")
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])
cnt <- read.table("normalized_counts_fln_X_PAR.bed", header = T, check.names = F)
col_names <- colnames(cnt)

set.seed(1 + i) 
cnt_shuffled <- cnt[, c(1:4, sample(5:ncol(cnt)))]
colnames(cnt_shuffled) <- col_names
colnames(cnt_shuffled)[1] <- "#chr"
write.table(cnt_shuffled, paste0("normalized_counts_fln_auto_", i, ".bed"), sep = "\t", row.names = FALSE, quote = FALSE)


# might be worth merging the FL VCFs
# otherwise this:
normalized_counts_fln_X <- normalized_counts_fln %>% filter(chrom == "X") #1894
normalized_counts_fln_PAR <- normalized_counts_fln %>% filter(chrom == "PAR") #1102
normalized_counts_fln_auto <- normalized_counts_fln %>% filter(chrom != "PAR" & chrom != "X") #14634

