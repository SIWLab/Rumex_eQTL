library(arrow)
library(dplyr)
library(ggplot2)
setwd("/ohta2/meng.yuan/rumex/eqtl/tensorqtl/FL")
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

FL <- read.table("FL.cis_qtl.txt", header = T) 
FL_egene <- FL %>% filter(qval < 0.1) 

FL_A1 <- read_parquet(paste0("FL_", i, ".cis_qtl_pairs.1.parquet"))
FL_A2 <- read_parquet(paste0("FL_", i, ".cis_qtl_pairs.2.parquet"))
FL_A3 <- read_parquet(paste0("FL_", i, ".cis_qtl_pairs.3.parquet"))
FL_A4 <- read_parquet(paste0("FL_", i, ".cis_qtl_pairs.4.parquet"))
FL_eqtl <- rbind(FL_A1, FL_A2, FL_A3, FL_A4) 

# identify eQTLs
FL_eqtl_full <- inner_join(FL_eqtl, FL_egene[,c(1,16:19)], by = "phenotype_id") 
FL_eqtl_sig <- FL_eqtl_full %>% filter(pval_nominal < pval_nominal_threshold)
cnt1 <- nrow(FL_eqtl_sig)

# remove gene with more than 10 eqtls
FL_eqtl_sig <- FL_eqtl_sig %>% group_by(phenotype_id) %>% mutate(count = n()) %>% ungroup() %>% filter(count <= 10) 
cnt2 <- nrow(FL_eqtl_sig)

# identify eqtls affecting multiple genes
FL_eqtl_multigene <- FL_eqtl_sig %>% group_by(variant_id) %>% filter(n()>1) %>% select(variant_id) %>% distinct()()

# select a random eqtl per gene
set.seed(1 + i) 
FL_eqtl_sig_random <- FL_eqtl_sig %>% group_by(phenotype_id) %>% slice_sample(n = 1) 

# remove eqtls affecting multiple genes 
FL_eqtl_sig_random <- FL_eqtl_sig_random %>% anti_join(FL_eqtl_multigene, by = c("variant_id")) %>% select(variant_id, af)

# get MAF
FL_eqtl_sig_random <- FL_eqtl_sig_random %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))
FL_eqtl_sig_random$maf <- round(FL_eqtl_sig_random$maf, digits = 3)
cnt3 <- nrow(FL_eqtl_sig_random)

p <- ggplot(FL_eqtl_sig_random, aes(x=maf)) + geom_histogram(binwidth=0.05, boundary = 0.05) + scale_x_continuous(limits = c(0.05, 0.5))
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y/sum(p_data$y) 
p_data <- p_data %>% select(x, xmin, xmax, y, y_prop)
write.table(p_data, file = paste0("FL_", i, "_eqtl_random_fake.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = F) 

cnt <- list(cnt1, cnt2, cnt3)
write.table(cnt, file=paste0("FL_", i, "_random.cnt"), row.names = FALSE, quote = FALSE, col.names = F)
