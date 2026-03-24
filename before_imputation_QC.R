library(data.table)
library(RSQLite)
library(dplyr)
library(stringr)
library(ggplot2)
library(R.utils)
library(reshape2)
library(tidyr)


LAA_snp_id <- fread('/home/sysadmin/Desktop/stroke_pipeline_test/LAA/LAA444_hwe_QCFiltered_fullid.txt', header = FALSE)
head(LAA_snp_id)
SVO_snp_id <- fread('/home/sysadmin/Desktop/stroke_pipeline_test/SVO/SVO342_hwe_QCFiltered_fullid.txt', header = FALSE)
head(SVO_snp_id)

LAA_SVO_snp_id <- intersect(LAA_snp_id, SVO_snp_id)
fwrite(LAA_SVO_snp_id, "/home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA_SVO_intersect_fullid.txt",
       col.names = FALSE, quote = FALSE)


LAA_SVO_snp_id <- fread('/home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA_SVO_QCFiltered_fullid.txt', header = FALSE)
TWB1_snp_id <- fread('/home/sysadmin/Desktop/stroke_pipeline_test/TWB1/TWB1.hg19_hwe_QCFiltered_fullid.txt', header = F)

TWB1_LAA_SVO_snp_id <- intersect(LAA_SVO_snp_id, TWB1_snp_id)
fwrite(TWB1_LAA_SVO_snp_id, "/home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_intersect_fullid.txt",
       col.names = FALSE, quote = FALSE)

gwas <- fread('/home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_GWAS_nocorrected.assoc.logistic')
manhattan(gwas, chr="CHR", bp="BP", snp="SNP", p="P",
          main="Manhattan Plot",
          genomewideline=-log10(1e-8), suggestiveline=-log10(1e-6))

gwas2 <- fread('/mnt/SP-siliconpower/TWBSTROKE_test/TWB1_LAA_SVO_test/TWB1_LAA_SVO_GWAS_nocorrected.assoc.logistic')
manhattan(gwas2, chr="CHR", bp="BP", snp="SNP", p="P",
          main="Manhattan Plot2",
          genomewideline=-log10(1e-8), suggestiveline=-log10(1e-6))


het <- fread('/home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_LDprune_het.het')
head(het)


# 計算 mean 與 SD
mu <- mean(het$F, na.rm = TRUE)
sd <- sd(het$F, na.rm = TRUE)

lower <- mu - 3 * sd
upper <- mu + 3 * sd

mu
sd
lower
upper

colnames(het)[1] <- 'FID'
het_outliers <- het[F < lower | F > upper, .(FID, IID)]
nrow(het_outliers)
table(het$F < lower, het$F > upper)
fwrite(
  het_outliers,
  "/home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/het_outliers.txt",
  sep = "\t",
  col.names = FALSE
)

hist(
  het$F,
  breaks = 50,
  main = "Heterozygosity (F) distribution",
  xlab = "Inbreeding coefficient (F)"
)

fam <- fread('/home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_QCFiltered.fam')



pca <- fread("/home/sysadmin/Desktop/stroke_pipeline_test/vcf/TWB1_LAA_SVO_noPCAoutliers_PCA.eigenvec")
head(pca)
colnames(pca)[1] <- 'FID'
# 如果 .fam 有表型資訊（case/control）
fam <- fread("/home/sysadmin/Desktop/stroke_pipeline_test/vcf/TWB1_LAA_SVO_noPCAoutliers.fam")
pca$Phenotype <- factor(fam$V6, labels = c("Control", "Case"))

# 計算PC1, PC2的Z-score
pca$PC1_z <- scale(pca$PC1)
pca$PC2_z <- scale(pca$PC2)

ggplot(pca, aes(x = PC1_z, y = PC2_z, color = Phenotype)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Control" = "#1b9e77", "Case" = "#d95f02")) +
  labs(
    title = "PCA of TWB1 & Stroke",
    x = "Principal Component 1_scale",
    y = "Principal Component 2_scale",
    color = "Phenotype"
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

outlier <- abs(pca$PC1_z) > 6 | abs(pca$PC2_z) > 6
pca$Outlier <- (abs(pca$PC1_z) > 6 | abs(pca$PC2_z) > 6)
table(pca$Outlier)
pca[,c('Phenotype', 'Outlier')][pca$Outlier == 'TRUE',]


pca_no_outlier <- pca[pca$Outlier == FALSE,]
table(pca_no_outlier$Outlier)
#colnames(pca_no_outlier)[1] <- 'FID'
write.table(
  pca_no_outlier[, c("FID", "IID")],
  "/home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/no_outliers_sample_id.txt",
  quote = FALSE, row.names = FALSE, col.names = FALSE
)




#gwas_noDuplicates <- fread('/mnt/SP-siliconpower/TWBSTROKE_test/TWB1_LAA_SVO_test/TWB1_LAA_SVO_GWAS_noDuplicates.assoc.logistic')
gwas_noPCAoutliers <- fread('/home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_GWAS_noPCAoutliers.assoc.logistic')
manhattan(gwas_noPCAoutliers, chr="CHR", bp="BP", snp="SNP", p="P",
          main="Manhattan Plot:noPCAoutliers",
          genomewideline=-log10(1e-8), suggestiveline=-log10(1e-6))


gwas_noDuplicates2 <- fread('/mnt/SP-siliconpower/TWBSTROKE_test/TWB1_LAA_SVO_test/TWB1_LAA_SVO_GWAS_noDuplicates.assoc.logistic')
manhattan(gwas_noDuplicates2, chr="CHR", bp="BP", snp="SNP", p="P",
          main="Manhattan Plot2:noPCAoutliers",
          genomewideline=-log10(1e-8), suggestiveline=-log10(1e-6))



