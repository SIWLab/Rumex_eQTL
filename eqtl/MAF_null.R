# reshuffle phenotype file
# 50 times, one one chromosome, no changes to the genotype files
setwd("path/to/your/directory")
data <- "/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/"

cnt <- read.table(paste0(data,"normalized_counts_mln_auto.bed"), header = T, check.names = F)
cnt <- cnt %>% filter(chr == "A4")
col_names <- colnames(cnt)

for (i in 1:50) {
set.seed(1 + i) 
cnt_shuffled <- cnt[, c(1:4, sample(5:ncol(cnt)))]
colnames(cnt_shuffled) <- col_names
colnames(cnt_shuffled)[1] <- "#chr"
write.table(cnt_shuffled, paste0("/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/pheno/normalized_counts_mln_auto_", i, ".bed"), sep = "\t", row.names = FALSE, quote = FALSE)
}







