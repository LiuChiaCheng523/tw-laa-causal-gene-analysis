library(data.table)
library(dplyr)

# 1. read gwas.bim
bim_file <- "/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/TWB1_LAA_SVO/output/chr9/TWB1_LAA_chr9_QCFiltered.bim"
bim <- fread(bim_file, fill = TRUE, header = FALSE)
head(bim)
nrow(bim)
setnames(bim, c("CHR","SNP","CM","BP","A1","A2"))
bim <- bim %>%
  filter(grepl("^[0-9]+$", CHR)) %>%    
  mutate(CHR = as.integer(CHR))             

# 2. read gwas.assoc.logistic
gwas <- fread("/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/TWB1_LAA_SVO/output/chr1/TWB1_LAA_chr1_gwas.assoc.logistic")
gwas_na_clean <- gwas %>%
  filter(!is.na(OR) & !is.na(STAT) & !is.na(P))
cat("SNP GWAS Z-score is not NA: ", nrow(gwas_na_clean), "\n")
gwas_na_clean <- gwas_na_clean %>% rename('Z'='STAT')

# 3. merge
fusion_sumstats <- merge(gwas_na_clean, bim[,c('CHR', 'BP', 'A2')], by=c('CHR', 'BP'), all=F)
sum(grepl("\\.", fusion_sumstats$SNP))
fusion_ready <- fusion_sumstats[, .(SNP, A1, A2, Z, CHR, BP)]
head(fusion_ready)
fwrite(fusion_ready,
       "/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/TWB1_LAA_SVO/output/chr9/TWB1_LAA_chr9_fusion_ready.sumstats",
       sep="\t", quote=FALSE, na="NA")


library(data.table)
library(dplyr)


chr_list <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
              'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')
for (chr_num in chr_list){
  # 1. read gwas.bim
  bim_file <- paste0('/home/sysadmin/Desktop/imputed_test/', chr_num, '/TWB1_SVO_', chr_num, '_QCFiltered.bim')
  bim <- fread(bim_file, fill = TRUE, header = FALSE)
  head(bim)
  nrow(bim)
  setnames(bim, c("CHR","SNP","CM","BP","A1","A2"))
  bim <- bim %>%
    filter(grepl("^[0-9]+$", CHR)) %>%    
    mutate(CHR = as.integer(CHR))
  
  # 2. read gwas.assoc.logistic
  gwas <- fread(paste0('/home/sysadmin/Desktop/imputed_test/', chr_num, '/TWB1_SVO_', chr_num, '_gwas.assoc.logistic'))
  gwas_na_clean <- gwas %>%
    filter(!is.na(OR) & !is.na(STAT) & !is.na(P))
  cat("SNP GWAS Z-score is not NA: ", nrow(gwas_na_clean), "\n")
  gwas_na_clean <- gwas_na_clean %>% rename('Z'='STAT')

  # 3. merge
  fusion_sumstats <- merge(gwas_na_clean, bim[,c('CHR', 'BP', 'A2')], by=c('CHR', 'BP'), all=F)
  sum(grepl("\\.", fusion_sumstats$SNP))
  fusion_ready <- fusion_sumstats[, .(SNP, A1, A2, Z, CHR, BP, P)]
  head(fusion_ready)
  fwrite(fusion_ready,
         paste0('/home/sysadmin/Desktop/dyc_lab/FUSION/GWAS/TWBSTROKE/hg19/SVO/TWB1_SVO_', chr_num, '_fusion_ready.sumstats'),
         sep="\t", quote=FALSE, na="NA")
}














LAA_SP_gwas_ss <- fread('/home/sysadmin/Desktop/dyc_lab/S_Predixcan/gwas/LAA_DR2_0.7/TWB1_LAA_DR2_07_chr1to22_spredixcan_ready.sumstats')
head(LAA_SP_gwas_ss)

LAA_FU_gwas_ss <- fread('/home/sysadmin/Desktop/dyc_lab/FUSION/GWAS/TWBSTROKE/hg19/LAA_DR2_0.7/TWB1_LAA_chr1_fusion_ready.sumstats')
head(LAA_FU_gwas_ss)

all_gwas <- {}
for (num in 1:22){
  gwas <- fread(paste0('/home/sysadmin/Desktop/dyc_lab/FUSION/GWAS/TWBSTROKE/hg19/SVO_DR2_0.7/',
                       'TWB1_SVO_chr', num, '_fusion_ready.sumstats'))
  all_gwas <- rbind(all_gwas, gwas)
}
head(all_gwas)
all_gwas <- all_gwas[,c('SNP', 'A1', 'A2', 'Z', 'CHR', 'BP')]
table(all_gwas$CHR)
fwrite(all_gwas,
       '/home/sysadmin/Desktop/dyc_lab/S_Predixcan/gwas/SVO_DR2_0.7/TWB1_SVO_DR2_07_chr1to22_spredixcan_ready.sumstats',
       sep="\t", quote=FALSE, na="NA")




