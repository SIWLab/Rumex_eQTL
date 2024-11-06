library(arrow)
library(dplyr)
library(ggplot2)
library(patchwork)
library(qvalue)
data <- "/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/"
# real
ML <- read.table(paste0(data, "ML.cis_qtl.txt"), header = T) # 
ML_A4 <- read_parquet(paste0(data,"ML.cis_qtl_pairs.4.parquet"))
ML_eqtl <- ML_A4

qval_ML <- qvalue(ML$pval_beta)
summary(qval_ML)
p1 <- hist(qval_ML)




ML <- read.table(paste0(data, "ML_1.cis_qtl.txt"), header = T) # 
#ggplot(ML, aes(x=pval_beta)) + geom_histogram(color="black", fill="white") + theme_classic() + xlab("p-value")+ggtitle("a")

ML_A4 <- read_parquet(paste0(data,"ML_1.cis_qtl_pairs.4.parquet"))
ML_eqtl <- ML_A4
ML <- calculate_qvalues(ML)
#ggplot(ML_eqtl, aes(x=pval_nominal)) + geom_histogram(color="black", fill="white") + theme_classic() + xlab("p-value")+ggtitle("a")
ML_eqtl_full <- inner_join(ML_eqtl, ML[,c(1,16:18)], by = "phenotype_id") 
ML_eqtl_sig <- ML_eqtl_full %>% filter(qval< 0.1) 
ML_eqtl_sig <- ML_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))
ML_eqtl_sig_min <- ML_eqtl_sig %>%
    group_by(phenotype_id) %>% 
    filter(pval_nominal == min(pval_nominal)) 
ML_eqtl_sig_min_af <- ML_eqtl_sig_min %>% dplyr::select(variant_id, af, maf) %>% distinct()

p<-ggplot(ML_eqtl_sig_min_af, aes(x=maf)) + geom_histogram(binwidth=0.01) 
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y/sum(p_data$y) # 0.096



