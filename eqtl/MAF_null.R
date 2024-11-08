library(arrow)
library(dplyr)
library(ggplot2)
setwd("/ohta2/meng.yuan/rumex/eqtl/tensorqtl/ml")
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

ML <- read.table("ML.cis_qtl.txt", header = T) 
ML_A1 <- read_parquet(paste0("ML_", i, ".cis_qtl_pairs.1.parquet"))
ML_A2 <- read_parquet(paste0("ML_", i, ".cis_qtl_pairs.2.parquet"))
ML_A3 <- read_parquet(paste0("ML_", i, ".cis_qtl_pairs.3.parquet"))
ML_A4 <- read_parquet(paste0("ML_", i, ".cis_qtl_pairs.4.parquet"))
ML_eqtl <- rbind(ML_A1, ML_A2, ML_A3, ML_A4) 

ML_eqtl_full <- inner_join(ML_eqtl, ML[,c(1,16:19)], by = "phenotype_id") 
ML_eqtl_sig <- ML_eqtl_full %>% filter(pval_nominal < pval_nominal_threshold) 
ML_eqtl_sig <- ML_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))
ML_eqtl_sig_min <- ML_eqtl_sig %>%
    group_by(phenotype_id) %>% 
    filter(pval_nominal == min(pval_nominal)) 
ML_eqtl_sig_min_af <- ML_eqtl_sig_min %>% dplyr::select(variant_id, af, maf) %>% distinct()
ML_eqtl_sig_min_af$maf <- round(ML_eqtl_sig_min_af$maf, digits = 3)
write.table(nrow(ML_eqtl_sig_min_af), file=paste0("ML_", i, ".cnt"), row.names = FALSE, quote = FALSE,col.names = F)

p <- ggplot(ML_eqtl_sig_min_af, aes(x=maf)) + geom_histogram(binwidth=0.05, , boundary = 0.05) +
    scale_x_continuous(limits = c(0.05, 0.5))
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y/sum(p_data$y) 
p_data <- p_data %>% select(x, xmin, xmax, y, y_prop)
write.table(p_data, file = paste0("ML_", i, "_eqtl_fake.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = F) 
