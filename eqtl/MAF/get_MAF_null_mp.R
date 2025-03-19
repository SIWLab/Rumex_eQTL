library(arrow)
library(dplyr)
library(ggplot2)
setwd("/ohta2/meng.yuan/rumex/eqtl/tensorqtl/MP")
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

MP_A1 <- read_parquet(paste0("MP_", i, ".cis_qtl_pairs.1.parquet"))
MP_A2 <- read_parquet(paste0("MP_", i, ".cis_qtl_pairs.2.parquet"))
MP_A3 <- read_parquet(paste0("MP_", i, ".cis_qtl_pairs.3.parquet"))
MP_A4 <- read_parquet(paste0("MP_", i, ".cis_qtl_pairs.4.parquet"))
MP_eqtl <- rbind(MP_A1, MP_A2, MP_A3, MP_A4) 

MP_eqtl_sig <- MP_eqtl %>% filter(pval_nominal <= 0.0160989562) 
MP_eqtl_sig <- MP_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))
MP_eqtl_sig <- MP_eqtl_sig %>% group_by(phenotype_id) %>% mutate(count = n()) %>% ungroup() %>% filter(count <= 10) 
write.table(nrow(MP_eqtl_sig), file=paste0("MP_", i, "_random.cnt1"), row.names = FALSE, quote = FALSE,col.names = F)

# random
set.seed(123)
MP_eqtl_sig_random <- MP_eqtl_sig %>% group_by(phenotype_id) %>% slice_sample(n = 1) 
MP_eqtl_sig_random_af <- MP_eqtl_sig_random %>% dplyr::select(variant_id, maf) 
MP_eqtl_sig_random_af <- MP_eqtl_sig_random_af %>% group_by(variant_id) %>% mutate(count = n()) %>% ungroup() %>% filter(count == 1) 
MP_eqtl_sig_random_af$maf <- round(MP_eqtl_sig_random_af$maf, digits = 3)
write.table(nrow(MP_eqtl_sig_random_af), file=paste0("MP_", i, "_random.cnt2"), row.names = FALSE, quote = FALSE,col.names = F)


p <- ggplot(MP_eqtl_sig_random_af, aes(x=maf)) + geom_histogram(binwidth=0.05, , boundary = 0.05) +
	scale_x_continuous(limits = c(0.05, 0.5))
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y/sum(p_data$y) 
p_data <- p_data %>% select(x, xmin, xmax, y, y_prop)
write.table(p_data, file = paste0("MP_", i, "_eqtl_random_fake.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = F) 

