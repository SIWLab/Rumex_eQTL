setwd("/ohta2/meng.yuan/rumex/eqtl/tensorqtl/Male")
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

cnt1 <- read.table("normalized_counts_mln_auto.bed", header = T, check.names = F)
cnt2 <- read.table("normalized_counts_mpn_auto.bed", header = T, check.names = F)

# reshuffle column names once
col_names <- colnames(cnt1)
set.seed(1 + i)
col_names_shuffled <- c(col_names[1:4], sample(col_names[5:length(col_names)]))

# ML
# assign the reshuffled column names to the original data
colnames(cnt1) <- col_names_shuffled
# re-sort the columns so they match with the genotype file
cnt1 <- cnt1[, c(1:4, 4 + order(names(cnt1)[5:ncol(cnt1)]))]
colnames(cnt1)[1] <- "#chr"
write.table(cnt1, paste0("normalized_counts_mln_auto_", i, ".bed"), sep = "\t", row.names = FALSE, quote = FALSE)

# MP
colnames(cnt2) <- col_names
cnt2 <- cnt2[, c(1:4, 4 + order(names(cnt2)[5:ncol(cnt2)]))]
colnames(cnt2)[1] <- "#chr"
write.table(cnt2, paste0("normalized_counts_mpn_auto_", i, ".bed"), sep = "\t", row.names = FALSE, quote = FALSE)
