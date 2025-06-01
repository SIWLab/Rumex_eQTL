library(arrow)
library(dplyr)
library(ggplot2)
setwd("/ohta2/meng.yuan/rumex/eqtl/tensorqtl/Male")
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

get_hist_prop <- function(data, binwidth = 0.05, boundary = 0.05, x_limits = c(0.05, 0.5)) {
p <- ggplot(data, aes(x = maf)) + geom_histogram(binwidth = binwidth, boundary = boundary) + scale_x_continuous(limits = x_limits)
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y / sum(p_data$y)
p_data <- p_data %>% select(x, xmin, xmax, y, y_prop)
return(p_data)}

##################################### get eqtls in ML ###################################
ML <- read.table("ML.cis_qtl.txt", header = T) 
ML_egene <- ML %>% filter(qval < 0.1) 

ML_A1 <- read_parquet(paste0("ML_", i, ".cis_qtl_pairs.1.parquet"))
ML_A2 <- read_parquet(paste0("ML_", i, ".cis_qtl_pairs.2.parquet"))
ML_A3 <- read_parquet(paste0("ML_", i, ".cis_qtl_pairs.3.parquet"))
ML_A4 <- read_parquet(paste0("ML_", i, ".cis_qtl_pairs.4.parquet"))
ML_eqtl <- rbind(ML_A1, ML_A2, ML_A3, ML_A4) 

# identify eQTLs
ML_eqtl_full <- inner_join(ML_eqtl, ML_egene[,c(1,16:19)], by = "phenotype_id") 
ML_eqtl_sig <- ML_eqtl_full %>% filter(pval_nominal < pval_nominal_threshold)
cnt1 <- nrow(ML_eqtl_sig)

# remove gene with more than 10 eqtls
ML_eqtl_sig <- ML_eqtl_sig %>% group_by(phenotype_id) %>% mutate(count = n()) %>% ungroup() %>% filter(count <= 10) 
cnt2 <- nrow(ML_eqtl_sig)

# identify eqtls affecting multiple genes
ML_eqtl_multigene <- ML_eqtl_sig %>% group_by(variant_id) %>% filter(n()>1) %>% select(variant_id) %>% distinct()

# select a random eqtl per gene
ML_eqtl_sig <- ML_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))
ML_eqtl_sig$maf <- round(ML_eqtl_sig$maf, digits = 3)
ML_eqtl_sig_random <- ML_eqtl_sig %>% group_by(phenotype_id) %>% slice_sample(n = 1) 

# remove eqtls affecting multiple genes 
ML_eqtl_sig_random <- ML_eqtl_sig_random %>% anti_join(ML_eqtl_multigene, by = c("variant_id")) %>% select(variant_id, maf)

# get MAF
cnt3 <- nrow(ML_eqtl_sig_random)
p_data <- get_hist_prop(ML_eqtl_sig_random)
write.table(p_data, file = paste0("ML_", i, "_eqtl_random_fake.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = F) 

cnt <- list(cnt1, cnt2, cnt3)
write.table(cnt, file=paste0("ML_", i, "_random.cnt"), row.names = FALSE, quote = FALSE, col.names = F)


##################################### get eqtls in MP ###################################
MP <- read.table("MP.cis_qtl.txt", header = T) 
MP_egene <- MP %>% filter(qval < 0.1) 

MP_A1 <- read_parquet(paste0("MP_", i, ".cis_qtl_pairs.1.parquet"))
MP_A2 <- read_parquet(paste0("MP_", i, ".cis_qtl_pairs.2.parquet"))
MP_A3 <- read_parquet(paste0("MP_", i, ".cis_qtl_pairs.3.parquet"))
MP_A4 <- read_parquet(paste0("MP_", i, ".cis_qtl_pairs.4.parquet"))
MP_eqtl <- rbind(MP_A1, MP_A2, MP_A3, MP_A4) 

# identify eQTLs
MP_eqtl_full <- inner_join(MP_eqtl, MP_egene[,c(1,16:19)], by = "phenotype_id") 
MP_eqtl_sig <- MP_eqtl_full %>% filter(pval_nominal < pval_nominal_threshold)
cnt1 <- nrow(MP_eqtl_sig)

# remove gene with more than 10 eqtls
MP_eqtl_sig <- MP_eqtl_sig %>% group_by(phenotype_id) %>% mutate(count = n()) %>% ungroup() %>% filter(count <= 10) 
cnt2 <- nrow(MP_eqtl_sig)

# identify eqtls affecting multiple genes
MP_eqtl_multigene <- MP_eqtl_sig %>% group_by(variant_id) %>% filter(n()>1) %>% select(variant_id) %>% distinct()

# select a random eqtl per gene
MP_eqtl_sig <- MP_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))
MP_eqtl_sig$maf <- round(MP_eqtl_sig$maf, digits = 3)
MP_eqtl_sig_random <- MP_eqtl_sig %>% group_by(phenotype_id) %>% slice_sample(n = 1) 

# remove eqtls affecting multiple genes 
MP_eqtl_sig_random <- MP_eqtl_sig_random %>% anti_join(MP_eqtl_multigene, by = c("variant_id")) %>% select(variant_id, maf)

# get MAF
cnt3 <- nrow(MP_eqtl_sig_random)
p_data <- get_hist_prop(MP_eqtl_sig_random)
write.table(p_data, file = paste0("MP_", i, "_eqtl_random_fake.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = F) 

cnt <- list(cnt1, cnt2, cnt3)
write.table(cnt, file=paste0("MP_", i, "_random.cnt"), row.names = FALSE, quote = FALSE, col.names = F)


############################## get eqtls shared vs specific #############################
ML_eqtl_sig <- ML_eqtl_sig %>% select(phenotype_id, variant_id, maf)
MP_eqtl_sig <- MP_eqtl_sig %>% select(phenotype_id, variant_id, maf)

tissue_n2 <- inner_join(ML_eqtl_sig, MP_eqtl_sig, by = c("phenotype_id", "variant_id", "maf"))
tissue_n1 <- bind_rows(anti_join(ML_eqtl_sig, MP_eqtl_sig, by = c("phenotype_id", "variant_id", "maf")), anti_join(MP_eqtl_sig, ML_eqtl_sig, by = c("phenotype_id", "variant_id", "maf")))

# remove eqtls affecting multiple genes 
tissue_n2 <- tissue_n2 %>% anti_join(ML_eqtl_multigene, by = c("variant_id")) %>% anti_join(MP_eqtl_multigene, by = c("variant_id"))
tissue_n1 <- tissue_n1 %>% anti_join(ML_eqtl_multigene, by = c("variant_id")) %>% anti_join(MP_eqtl_multigene, by = c("variant_id"))
 
tissue_n2_random <- tissue_n2 %>% group_by(phenotype_id) %>% slice_sample(n = 1)  
tissue_n1_random <- tissue_n1 %>% group_by(phenotype_id) %>% slice_sample(n = 1) 

# get MAF
p_data <- get_hist_prop(tissue_n1_random)
write.table(p_data, file = paste0("tissue_n1_", i, "_fake.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = F) 
p_data <- get_hist_prop(tissue_n2_random)
write.table(p_data, file = paste0("tissue_n2_", i, "_fake.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = F) 

cnt <- list(nrow(tissue_n1), nrow(tissue_n2), nrow(tissue_n1_random), nrow(tissue_n2_random))
write.table(cnt, file=paste0("tissue_n1_n2_", i, "_random.cnt"), row.names = FALSE, quote = FALSE, col.names = F)
