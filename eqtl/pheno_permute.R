setwd("/ohta2/meng.yuan/rumex/eqtl/tensorqtl/ml")
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])
cnt <- read.table("normalized_counts_mln_auto.bed", header = T, check.names = F)
col_names <- colnames(cnt)

set.seed(1 + i) 
cnt_shuffled <- cnt[, c(1:4, sample(5:ncol(cnt)))]
colnames(cnt_shuffled) <- col_names
colnames(cnt_shuffled)[1] <- "#chr"
write.table(cnt_shuffled, paste0("normalized_counts_mln_auto_", i, ".bed"), sep = "\t", row.names = FALSE, quote = FALSE)

