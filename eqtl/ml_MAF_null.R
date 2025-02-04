library(arrow)
library(dplyr)
library(ggplot2)
data <- "/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/"
calculate_qvalues <- function(res_df, fdr = 0.1, qvalue_lambda = NULL) {
    # Number of phenotypes tested
    num_phenotypes <- nrow(res_df)
    cat(sprintf("Number of phenotypes tested: %d\n", num_phenotypes))
    
    # Determine which p-value column to use
    if ("pval_beta" %in% colnames(res_df) && !all(is.na(res_df$pval_beta))) {
        pval_col <- "pval_beta"
        correlation <- cor(res_df$pval_perm, res_df$pval_beta, use = "complete.obs")
        cat(sprintf("Correlation between Beta-approximated and empirical p-values: %.4f\n", correlation))
    } else {
        pval_col <- "pval_perm"
        cat("WARNING: no beta-approximated p-values found, using permutation p-values instead.\n")
    }
    
    # Calculate q-values
    if (!is.null(qvalue_lambda)) {
        cat(sprintf("Calculating q-values with lambda = %.3f\n", qvalue_lambda))
        qobj <- qvalue(res_df[[pval_col]], lambda = qvalue_lambda)
    } else {
        qobj <- qvalue(res_df[[pval_col]])
    }
    
    res_df$qval <- qobj$qvalues
    pi0 <- qobj$pi0
    cat(sprintf("Proportion of significant phenotypes (1-pi0): %.2f\n", 1 - pi0))
    cat(sprintf("QTL phenotypes @ FDR %.2f: %d\n", fdr, sum(res_df$qval <= fdr)))
    
    # Determine global min(p) significance threshold and calculate nominal p-value threshold for each gene
    if (pval_col == "pval_beta") {
        lb <- sort(res_df[res_df$qval <= fdr, "pval_beta"])
        ub <- sort(res_df[res_df$qval > fdr, "pval_beta"])
        
        if (length(lb) > 0) {
            lb <- tail(lb, 1)
            if (length(ub) > 0) {
                ub <- head(ub, 1)
                pthreshold <- (lb + ub) / 2
            } else {
                pthreshold <- lb
            }
            cat(sprintf("min p-value threshold @ FDR %.2f: %.6g\n", fdr, pthreshold))
            res_df$pval_nominal_threshold <- qbeta(pthreshold, res_df$beta_shape1, res_df$beta_shape2)
        }
    }
    
    return(res_df)
}

# real
ML <- read.table("/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/ML.cis_qtl.txt", header = T) 
ggplot(ML, aes(x=pval_beta)) + geom_histogram(color="black", fill="white") + theme_classic() + xlab("p-value")+ggtitle("pval_beta")

ggplot(ML, aes(x=pval_nominal)) + geom_histogram(color="black", fill="white") + theme_classic() + xlab("p-value")+ggtitle("pval_nominal")

ggplot(ML, aes(x=pval_perm)) + geom_histogram(color="black", fill="white") + theme_classic() + xlab("p-value")+ggtitle("pval_perm")


#ML2 <- ML[,-c(18,19)]
#ML2 <- calculate_qvalues(ML2)
# all.equal(ML$qval, ML2$qval)
# all.equal(ML$pval_nominal_threshold, ML2$pval_nominal_threshold)
# ML2 <- ML2[,c(18,19)]
# ML2 <- cbind(ML2, ML[,c(18,19)])
# write.table(ML2, file=paste0(data, "ML2.txt"), , row.names = F, quote = F, sep = "\t")
ML_egene <- ML2 %>% filter(qval < 0.1) # 2425
ML_egene <- ML %>% filter(qval < 0.1) # 2425


# real eqtls
ML_eqtl <- read_parquet("/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/ML.cis_qtl_pairs.4.parquet")
#ggplot(ML_eqtl, aes(x=pval_nominal)) + geom_histogram(color="black", fill="white", binwidth = 0.005) + theme_classic() + xlab("p-value")+ggtitle("a")

ML_min <- ML_eqtl %>%
    group_by(phenotype_id) %>% 
    filter(pval_nominal == min(pval_nominal))
length(unique(ML$phenotype_id))
length(unique(ML_eqtl$phenotype_id))
length(unique(ML_min$phenotype_id))

ML_min <- ML_min %>% dplyr::select(phenotype_id, pval_nominal) %>% distinct()
ML_min <- calculate_qvalues(ML_min)

ML_eqtl_full <- inner_join(ML_eqtl, ML[,c(1,16:19)], by = "phenotype_id") 



ML$qval2 <- qvalue(ML$pval_nominal)$qvalues
thresholded_data <- ML[ML$qval2 < 0.1, ]

max(thresholded_data$pval_nominal) # 0.995526
ML_eqtl_sig2 <- ML_eqtl_full %>% filter(pval_nominal < 0.995526)



ML_eqtl_sig <- ML_eqtl_full %>% filter(pval_nominal < pval_nominal_threshold) # 2812
ML_eqtl_sig <- ML_eqtl_sig %>% mutate(maf = ifelse(af <= 0.5, af, 1 - af))


ML_eqtl_sig <-read.csv( "/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/ML_eqtl_sig.csv") # 17283
MP_eqtl_sig <-read.csv( "/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/eqtl/MP_eqtl_sig.csv") # 9835


ML_eqtl_sig_min <- ML_eqtl_sig %>%
    group_by(phenotype_id) %>% 
    filter(pval_nominal == min(pval_nominal))  # 2879
ML_eqtl_sig_min_af <- ML_eqtl_sig_min %>% dplyr::select(variant_id, af, maf) %>% distinct()
ML_eqtl_sig_min_af$maf <- round(ML_eqtl_sig_min_af$maf, digits = 3)
    
p<-ggplot(ML_eqtl_sig_min_af, aes(x=maf)) + geom_histogram(binwidth=0.05, , boundary = 0.05) +
    scale_x_continuous(limits = c(0.05, 0.5))
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y/sum(p_data$y) 

p_data <- p_data %>% select(x, xmin, xmax, y, y_prop)
write.table(p_data, file = paste0(data, "ML_eqtl_real.txt"), row.names = F, quote = F, sep = "\t")
####
MP_eqtl_sig_min <- MP_eqtl_sig %>%
    group_by(phenotype_id) %>% 
    filter(pval_nominal == min(pval_nominal))  # 2879
MP_eqtl_sig_min_af <- MP_eqtl_sig_min %>% dplyr::select(variant_id, af, maf) %>% distinct()
MP_eqtl_sig_min_af$maf <- round(MP_eqtl_sig_min_af$maf, digits = 3)

p<-ggplot(MP_eqtl_sig_min_af, aes(x=maf)) + geom_histogram(binwidth=0.05, , boundary = 0.05) +
    scale_x_continuous(limits = c(0.05, 0.5))
p_data <- ggplot_build(p)$data[[1]]
p_data$y_prop <- p_data$y/sum(p_data$y) 

p_data <- p_data %>% select(x, xmin, xmax, y, y_prop)
write.table(p_data, file = paste0(data, "MP_eqtl_real.txt"), row.names = F, quote = F, sep = "\t")

#####

# fake eqtls
# a loop with 50
data <- "/Users/yuanmeng/Library/CloudStorage/OneDrive-UniversityofToronto/Manuscripts/rumex_eqtl/ml/"
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

ML_eqtl_sig_random <- ML_eqtl_sig %>%
    group_by(phenotype_id) %>% 
    slice_sample(n = 1)

ML_eqtl_sig_min_af <- ML_eqtl_sig_min %>% dplyr::select(variant_id, af, maf) %>% distinct()
ML_eqtl_sig_min_af$maf <- round(ML_eqtl_sig_min_af$maf, digits = 3)
#fake_eqtl_cnt <- rbind(fake_eqtl_cnt, data.frame(n = nrow(ML_eqtl_sig_min_af)))

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

# maf_all <- rbind(maf_real, maf_null)
# maf_all$category <- factor(maf_all$category)
# table(maf_all$category)


# ggplot(maf_real, aes(x=x_axis, y=y_prop)) + geom_point() + theme_classic()  + labs(x="MAF", y="Proportion") 
# ggplot(maf_null, aes(x=x_axis, y=y_prop)) + geom_point() + theme_classic()  + labs(x="MAF", y="Proportion")


#ggplot(maf_all, aes(x=x_axis, y=y_prop, color=category)) + geom_point(fill = NA, shape=21) + theme_classic() + scale_color_manual(values = c("lightgrey","black")) + labs(x="MAF", y="Proportion")
maf_null <- read.table(paste0(data, "ML_MAF_null.txt"))
colnames(maf_null) <- c("x", "xmin","xmax","y", "y_prop")
maf_null$category <- "null"
maf_null <- maf_null %>% mutate(x_axis = paste(xmin, xmax, sep = "-"))

maf_real <- read.table(paste0(data, "ML_eqtl_real.txt"), header = T)
maf_real$category <- "real"
maf_real <- maf_real %>% mutate(x_axis = paste(xmin, xmax, sep = "-"))

p1 <- ggplot(maf_null, aes(x_axis, y_prop)) + 
    geom_point(fill = NA, colour = "lightgrey")+
    geom_point(data = maf_real, color = "black", fill="black")+ 
    theme_classic()+ labs(x="MAF", y="Proportion", title = "a")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 10))

maf_null <- read.table(paste0(data, "MP_MAF_null.txt"))
colnames(maf_null) <- c("x", "xmin","xmax","y", "y_prop")
maf_null$category <- "null"
maf_null <- maf_null %>% mutate(x_axis = paste(xmin, xmax, sep = "-"))

maf_real <- read.table(paste0(data, "MP_eqtl_real.txt"), header = T)
maf_real$category <- "real"
maf_real <- maf_real %>% mutate(x_axis = paste(xmin, xmax, sep = "-"))

p2 <- ggplot(maf_null, aes(x_axis, y_prop)) + 
    geom_point(fill = NA, colour = "lightgrey")+
    geom_point(data = maf_real, color = "black", fill="black")+ 
    theme_classic()+ labs(x="MAF", y="Proportion", title = "b")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 10))


p1 + p2

