library(arrow)
library(dplyr)
library(ggplot2)
setwd("/ohta2/meng.yuan/rumex/eqtl/tensorqtl/FL")
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

FL_A1 <- read_parquet(paste0("FL_", i, ".cis_qtl_pairs.1.parquet"))
FL_A2 <- read_parquet(paste0("FL_", i, ".cis_qtl_pairs.2.parquet"))
FL_A3 <- read_parquet(paste0("FL_", i, ".cis_qtl_pairs.3.parquet"))
FL_A4 <- read_parquet(paste0("FL_", i, ".cis_qtl_pairs.4.parquet"))
FL_eqtl <- rbind(FL_A1, FL_A2, FL_A3, FL_A4) 

FL_eqtl_sig <- FL_eqtl %>% filter(pval_nominal <= 0.030655064) 
FL_eqtl_sig <- FL_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))

# random
set.seed(123)
FL_eqtl_sig_random <- FL_eqtl_sig %>% group_by(phenotype_id) %>% slice_sample(n = 1) 
FL_eqtl_sig_random_af <- FL_eqtl_sig_random %>% dplyr::select(variant_id, maf) 
FL_eqtl_sig_random_af <- FL_eqtl_sig_random_af %>% group_by(variant_id) %>% mutate(count = n()) %>% ungroup() %>% filter(count == 1) 
FL_eqtl_sig_random_af$maf <- round(FL_eqtl_sig_random_af$maf, digits = 3)
write.table(nrow(FL_eqtl_sig_random_af), file=paste0("FL_", i, "_random.cnt"), row.names = FALSE, quote = FALSE,col.names = F)


p <- ggplot(FL_eqtl_sig_random_af, aes(x=maf)) + geom_histogram(binwidth=0.05, , boundary = 0.05) +
scale_x_continuous(limits = c(0.05, 0.5))
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y/sum(p_data$y) 
p_data <- p_data %>% select(x, xmin, xmax, y, y_prop)
write.table(p_data, file = paste0("FL_", i, "_eqtl_random_fake.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = F) 
