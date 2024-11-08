library(arrow)
library(dplyr)
library(ggplot2)
data <- "/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/ml/"
# real
ML <- read.table("/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/ML.cis_qtl.txt", header = T) 
#ggplot(ML, aes(x=pval_beta)) + geom_histogram(color="black", fill="white") + theme_classic() + xlab("p-value")+ggtitle("a")
# real eqtls
ML_eqtl <- read_parquet("/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/ML.cis_qtl_pairs.4.parquet")
#ggplot(ML_eqtl, aes(x=pval_nominal)) + geom_histogram(color="black", fill="white", binwidth = 0.005) + theme_classic() + xlab("p-value")+ggtitle("a")

ML_eqtl_full <- inner_join(ML_eqtl, ML[,c(1,16:19)], by = "phenotype_id") 
ML_eqtl_sig <- ML_eqtl_full %>% filter(pval_nominal < pval_nominal_threshold) # 2812

ML_eqtl_sig <- ML_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))
ML_eqtl_sig_min <- ML_eqtl_sig %>%
    group_by(phenotype_id) %>% 
    filter(pval_nominal == min(pval_nominal))  # 476
ML_eqtl_sig_min_af <- ML_eqtl_sig_min %>% dplyr::select(variant_id, af, maf) %>% distinct()
ML_eqtl_sig_min_af$maf <- round(ML_eqtl_sig_min_af$maf, digits = 3)
    
p<-ggplot(ML_eqtl_sig_min_af, aes(x=maf)) + geom_histogram(binwidth=0.05, , boundary = 0.05) +
    scale_x_continuous(limits = c(0.05, 0.5))
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y/sum(p_data$y) 

p_data <- p_data %>% select(x, xmin, xmax, y, y_prop)
write.table(p_data, file = paste0(data, "ML_A4_eqtl_real.txt"), row.names = F, quote = F, sep = "\t")


# fake eqtls
ML_eqtl <- read_parquet(paste0(data,"ML_30.cis_qtl_pairs.4.parquet"))
# ggplot(ML_eqtl, aes(x=pval_nominal)) + geom_histogram(color="black", fill="white", binwidth = 0.005) + theme_classic() + xlab("p-value")+ggtitle("a")

# ML_eqtl_full <- inner_join(ML_eqtl, ML[,c(1,16:19)], by = "phenotype_id") 
# ML_eqtl_sig <- ML_eqtl_full %>% filter(pval_nominal < pval_nominal_threshold) 
# nrow(ML_eqtl_sig)
# 
# ML_eqtl_sig <- ML_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))
# ML_eqtl_sig_min <- ML_eqtl_sig %>%
#     group_by(phenotype_id) %>% 
#     filter(pval_nominal == min(pval_nominal)) 
# ML_eqtl_sig_min_af <- ML_eqtl_sig_min %>% dplyr::select(variant_id, af, maf) %>% distinct()
# ML_eqtl_sig_min_af$maf <- round(ML_eqtl_sig_min_af$maf, digits = 3)
# 
# p<-ggplot(ML_eqtl_sig_min_af, aes(x=maf)) + geom_histogram(binwidth=0.05) 
# p_data <- ggplot_build(p)$data[[1]]
# p_data$y_prop <- p_data$y/sum(p_data$y) 
# p_data <- p_data %>% select(x, y, y_prop)
# write.table(p_data, file = paste0(data, "ML_30_eqtl_fake.txt"), row.names = F, quote = F, sep = "\t")

# a loop with 50
fake_eqtl_cnt <- data.frame(n = integer())

for (i in 1:50) {
ML_eqtl <- read_parquet(paste0(data, "ML_", i, ".cis_qtl_pairs.4.parquet"))
ML_eqtl_full <- inner_join(ML_eqtl, ML[,c(1,16:19)], by = "phenotype_id") 
ML_eqtl_sig <- ML_eqtl_full %>% filter(pval_nominal < pval_nominal_threshold) 
nrow(ML_eqtl_sig)
ML_eqtl_sig <- ML_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))
ML_eqtl_sig_min <- ML_eqtl_sig %>%
    group_by(phenotype_id) %>% 
    filter(pval_nominal == min(pval_nominal)) 
ML_eqtl_sig_min_af <- ML_eqtl_sig_min %>% dplyr::select(variant_id, af, maf) %>% distinct()
ML_eqtl_sig_min_af$maf <- round(ML_eqtl_sig_min_af$maf, digits = 3)
fake_eqtl_cnt <- rbind(fake_eqtl_cnt, data.frame(n = nrow(ML_eqtl_sig_min_af)))

p<-ggplot(ML_eqtl_sig_min_af, aes(x=maf)) + geom_histogram(binwidth=0.05, , boundary = 0.05) +
    scale_x_continuous(limits = c(0.05, 0.5))
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y/sum(p_data$y) 
p_data <- p_data %>% select(x, xmin, xmax, y, y_prop)
write.table(p_data, file = paste0(data, "ML_", i, "_eqtl_fake.txt"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = F) 
}  

range(fake_eqtl_cnt$n)
mean(fake_eqtl_cnt$n)
median(fake_eqtl_cnt$n)

# plot the null
maf_null <- read.table(paste0(data, "ML_A4_MAF_null.txt"))
colnames(maf_null) <- c("x", "xmin","xmax","y", "y_prop")
maf_null$category <- "null"
maf_null <- maf_null %>% mutate(x_axis = paste(xmin, xmax, sep = "-"))

maf_real <- read.table(paste0(data, "ML_A4_eqtl_real.txt"), header = T)
maf_real$category <- "real"
maf_real <- maf_real %>% mutate(x_axis = paste(xmin, xmax, sep = "-"))

maf_all <- rbind(maf_real, maf_null)
maf_all$category <- factor(maf_all$category)
table(maf_all$category)


ggplot(maf_real, aes(x=x_axis, y=y_prop)) + geom_point() + theme_classic()  + labs(x="MAF", y="Proportion") +ylim(0,0.61)
ggplot(maf_null, aes(x=x_axis, y=y_prop)) + geom_point() + theme_classic()  + labs(x="MAF", y="Proportion")+ylim(0,0.61)


#ggplot(maf_all, aes(x=x_axis, y=y_prop, color=category)) + geom_point(fill = NA, shape=21) + theme_classic() + scale_color_manual(values = c("lightgrey","black")) + labs(x="MAF", y="Proportion")

ggplot(maf_null, aes(x_axis, y_prop)) + 
    geom_point(fill = NA, colour = "lightgrey")+
    geom_point(data = maf_real, color = "black", fill="black")+ 
    theme_classic()+ labs(x="MAF", y="Proportion")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 10))


