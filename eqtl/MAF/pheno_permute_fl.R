setwd("/ohta2/meng.yuan/rumex/eqtl/tensorqtl/FL")
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])
cnt <- read.table("normalized_counts_fln_X_PAR.bed", header = T, check.names = F)
col_names <- colnames(cnt)

set.seed(1 + i) 
cnt_shuffled <- cnt[, c(1:4, sample(5:ncol(cnt)))]
colnames(cnt_shuffled) <- col_names

cnt_shuffled$chr <- as.character(cnt_shuffled$chr)
normalized_counts_fln_X <- cnt_shuffled %>% filter(chr == "5") 
normalized_counts_fln_PAR <- cnt_shuffled %>% filter(chr == "6") 
normalized_counts_fln_auto <- cnt_shuffled %>% filter(chr != "5" & chr != "6") 

colnames(normalized_counts_fln_X)[1] <- "#chr"
colnames(normalized_counts_fln_PAR)[1] <- "#chr"
colnames(normalized_counts_fln_auto)[1] <- "#chr"
write.table(normalized_counts_fln_X, file = paste0("normalized_counts_fln_X_", i, ".bed"), row.names = FALSE, quote = F, sep = "\t")
write.table(normalized_counts_fln_PAR, file = paste0("normalized_counts_fln_PAR_", i, ".bed"), row.names = FALSE, quote = F, sep = "\t")
write.table(normalized_counts_fln_auto, file = paste0("normalized_counts_fln_auto_", i, ".bed"), row.names = FALSE, quote = F, sep = "\t")
