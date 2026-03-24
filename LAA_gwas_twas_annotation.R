library(data.table)
library(RSQLite)
library(dplyr)
library(stringr)
library(ggplot2)
library(R.utils)
library(reshape2)
library(tidyr)
library(qqman)
library(biomaRt)
library(ggrepel)
#BiocManager::install("rtracklayer")
library(rtracklayer)
library(ggvenn)

# load ensembl database -----
listEnsemblArchives()
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = 'GRCh37'
)
listAttributes(ensembl)

# load gencode database -----
# gencode.v26
gtf <- import("/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/GENCODE/gencode.v26.annotation.gtf")
gene_gr <- gtf[gtf$type == "gene"]

gene_dt <- data.table(
  #ensembl_gene_id = gsub("\\..*", "", gene_gr$gene_id)
  ensembl_gene_id = gene_gr$gene_id,
  gene_name       = gene_gr$gene_name,
  chr             = as.character(seqnames(gene_gr)),
  start           = start(gene_gr),
  end             = end(gene_gr),
  gene_type       = gene_gr$gene_type
)
gene_dt[, chr := sub("^chr", "", chr)]
gene_dt <- gene_dt[chr %in% c(as.character(1:22))]
head(gene_dt)
#setorder(gene_dt, chr, start)

# gencode.v19
gtf_v19 <- import("/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/GENCODE/gencode.v19.annotation.gtf")
gene_gr_v19 <- gtf_v19[gtf_v19$type == "gene"]
gene_dt_v19 <- data.table(
  #ensembl_gene_id = gsub("\\..*", "", gene_gr_v19$gene_id),  
  ensembl_gene_id = gene_gr_v19$gene_id,
  gene_name       = gene_gr_v19$gene_name,
  chr             = as.character(seqnames(gene_gr_v19)),
  start           = start(gene_gr_v19),
  end             = end(gene_gr_v19),
  gene_type       = gene_gr_v19$gene_type
)
gene_dt_v19[, chr := sub("^chr", "", chr)]
gene_dt_v19 <- gene_dt_v19[chr %in% c(as.character(1:22))]
length(unique(gene_dt_v19$ensembl_gene_id))


# load FUSION v7 wgt filename
weight_dir <- "/home/sysadmin/Desktop/dyc_lab/FUSION/WEIGHTS_v7/GTEx.ALL"
wgt_files <- list.files(
  weight_dir,
  pattern = "\\.wgt\\.RDat$",
  full.names = TRUE,
  recursive = TRUE
)

gene_map_wgt <- data.table(filename = basename(wgt_files))
gene_map_wgt[, ensembl_gene_id :=
               str_extract(filename, "ENSG[0-9]+\\.[0-9]+")
]
gene_map_wgt[, gene_id :=
               sub("\\..*$", "", ensembl_gene_id)
]
gene_map_wgt[, gene_name :=
               sub(".*ENSG[0-9]+\\.[0-9]+\\.", "", filename)
]
gene_map_wgt[, gene_name :=
               sub("\\.wgt\\.RDat$", "", gene_name)
]
gene_map_wgt <- unique(
  gene_map_wgt[, .(gene_name, ensembl_gene_id, gene_id)]
)

FU_v7_wgt_dup_symbol <- gene_map_wgt[
  , .(N = uniqueN(ensembl_gene_id)),
  by = gene_name
][N > 1]


LAA_FU_v7_annot <- merge(
  LAA_FU_v7,
  gene_map_wgt[, .(gene_name, ensembl_gene_id)],
  by="gene_name",
  all.x=TRUE
)


# LAA gwas fine-mapping (COJO) -----
annotate_one_snp <- function(chr, bp, snp_id, window, mart) {
  genes <- getBM(
    attributes = c(
      "hgnc_symbol",
      "ensembl_gene_id_version",
      "chromosome_name",
      "start_position",
      "end_position",
      "gene_biotype"
    ),
    filters = c("chromosome_name", "start", "end"),
    values = list(
      as.character(chr),
      bp - window,
      bp + window
    ),
    mart = mart
  )
  
  if (nrow(genes) == 0) return(NULL)
  
  genes$SNP    <- snp_id
  genes$SNP_bp <- bp
  #genes$distance <- pmin(
  #  abs(genes$start_position - bp),
  #  abs(genes$end_position   - bp)
  #)
  genes$distance <- ifelse(
    bp >= genes$start_position & bp <= genes$end_position,
    0,
    pmin(
      abs(genes$start_position - bp),
      abs(genes$end_position - bp)
    )
  )
  genes$location <- ifelse(
    bp >= genes$start_position & bp <= genes$end_position,
    "intragenic",
    ifelse(
      bp < genes$start_position,
      "upstream",
      "downstream"
    )
  )
  
  genes
}

chr_list <- c(1,2,3,4,6,7,9,10,11,14,15,16,17,18,19,22)
window <- 500000

LAA_all_cojo_snps <- {}
for (chr in chr_list){
  cojo_file <- paste0(
    "/mnt/data/stroke_pipeline_test/cojo_result/LAA/",
    "TWB1_LAA_chr", chr, "_DR2_0.7_cojo_p1e5.jma.cojo"
  )
  cojo_result <- fread(cojo_file)
  cojo_snps <- data.frame(
    SNP = cojo_result$SNP,
    chr = cojo_result$Chr,
    bp  = cojo_result$bp
  )
  LAA_all_cojo_snps <- rbind(LAA_all_cojo_snps,  cojo_snps)
}

LAA_all_anno_list <- list()
for (chr in chr_list) {
  message("Annotating chr", chr)
  cojo_file <- paste0(
    "/mnt/data/stroke_pipeline_test/cojo_result/LAA/",
    "TWB1_LAA_chr", chr, "_DR2_0.7_cojo_p1e5.jma.cojo"
  )
  
  if (!file.exists(cojo_file)) {
    warning("File not found: ", cojo_file)
    next
  }
  
  cojo_result <- fread(cojo_file)
  
  if (nrow(cojo_result) == 0) next
  
  cojo_snps <- data.frame(
    SNP = cojo_result$SNP,
    chr = cojo_result$Chr,
    bp  = cojo_result$bp
  )
  
  chr_anno <- lapply(
    seq_len(nrow(cojo_snps)),
    function(i) {
      annotate_one_snp(
        chr    = cojo_snps$chr[i],
        bp     = cojo_snps$bp[i],
        snp_id = cojo_snps$SNP[i],
        window = window,
        mart   = ensembl
      )
    }
  )
  
  chr_anno <- do.call(rbind, chr_anno)
  
  if (!is.null(chr_anno)) {
    LAA_all_anno_list[[paste0("chr", chr)]] <- chr_anno
  }
}

LAA_anno <- do.call(rbind, LAA_all_anno_list)
#setDT(LAA_anno)
head(LAA_anno)
unique(LAA_anno$SNP)
LAA_anno[LAA_anno$SNP == 'rs2023936',]
rs2023936
table(LAA_anno$chromosome_name)
table(LAA_anno$gene_biotype)   # protein_coding 221
length(unique(LAA_anno$ensembl_gene_id_version))
write.csv(LAA_anno, 
          '/mnt/data/stroke_pipeline_test/cojo_result/table/LAA_cojo_p1e5_annotation_genes_20260307.csv',
          row.names = F)

LAA_protein_nearest_gene <- LAA_anno %>%
  filter(
    gene_biotype == "protein_coding",
    hgnc_symbol != ""
  ) %>%
  group_by(SNP) %>%
  slice_min(distance, n = 1, with_ties = FALSE) %>%
  ungroup()


#write.csv(LAA_protein_nearest_gene, 
#          '/mnt/data/stroke_pipeline_test/cojo_result/table/LAA_cojo_p1e5_protein_coding_nearest_gene.csv',
#          row.names = F)


# LAA GWAS Manhattan plot -----
base_folder <- '/mnt/data/stroke_pipeline_test/imputed_vcf/'
all_LAA_gwas <- {}
for (num in 1:22){
  #num <- 1
  chr_num <- paste0('chr', num)
  gwas_path <- paste0(base_folder, chr_num, '/DR2_0.7/TWB1_LAA_', chr_num, '_DR2_0.7.PHENO1.glm.logistic.hybrid')
  gwas <- fread(gwas_path)
  head(gwas)
  gwas <- gwas[!is.na(gwas$P),]
  all_LAA_gwas <-  rbind(all_LAA_gwas, gwas)
}
all_LAA_gwas_plot <- all_LAA_gwas[,c('#CHROM', 'POS', 'ID', 'P')]
head(all_LAA_gwas_plot)
summary(all_LAA_gwas$A1_FREQ)
head(all_LAA_gwas)
sum(all_LAA_gwas$ALT == all_LAA_gwas$A1, na.rm = TRUE)
mean(all_LAA_gwas$ALT == all_LAA_gwas$A1, na.rm = TRUE)

gwas_df <- all_LAA_gwas_plot %>%
  rename(
    CHR = `#CHROM`,
    BP  = POS,
    SNP = ID,
    P   = P
  ) %>%
  filter(
    CHR %in% 1:22,
    !is.na(P),
    is.finite(P)
  ) %>%
  mutate(
    CHR = as.numeric(CHR),
    BP  = as.numeric(BP)
  ) %>%
  arrange(CHR, BP)

chr_max_gwas <- gwas_df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

gwas_df <- gwas_df %>%
  left_join(chr_max_gwas, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

label_df_gwas <- LAA_protein_nearest_gene %>%
  distinct(SNP, hgnc_symbol, SNP_bp) %>%
  left_join(
    gwas_df,
    by = c("SNP" = "SNP")
  ) %>%
  filter(!is.na(P))

gwas_df <- gwas_df %>%
  mutate(chr_color = factor(CHR %% 2))


p_gwas <- ggplot(gwas_df, aes(x = BP_cum, y = -log10(P))) +
  
  # ???? GWAS SNP
  geom_point(
    aes(color = chr_color),
    size = 0.4,
    alpha = 0.7
  ) +
  
  scale_color_manual(
    values = c("grey60", "black"),
    guide = "none"
  ) +
  
  geom_hline(
    yintercept = -log10(1e-5),
    color = "blue",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  
  geom_point(
    data = label_df_gwas,
    aes(x = BP_cum, y = -log10(P)),
    color = "red",
    size = 2
  ) +
  
  geom_text_repel(
    data = label_df_gwas,
    aes(label = hgnc_symbol),
    color = "red",
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  
  labs(
    title = "GWAS Manhattan Plot",
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

chr_centers_gwas <- gwas_df %>%
  group_by(CHR) %>%
  summarise(center = mean(range(BP_cum)))

p_gwas <- p_gwas +
  scale_x_continuous(
    breaks = chr_centers_gwas$center,
    labels = chr_centers_gwas$CHR
  )

p_gwas <- p_gwas +
  scale_y_continuous(
    breaks = seq(0, ceiling(max(-log10(gwas_df$P))), by = 1)
  )

ggsave("/mnt/data/stroke_pipeline_test/cojo_result/hanmanttan plot/LAA_cojo_protein_coding_manhanttan_20260308.png",
       p_gwas, width = 12, height = 6, dpi=600)



# load LAA FUSION v8 result -----
LAA_FU_v8_tissue <- fread('/home/sysadmin/Desktop/dyc_lab/FUSION/WEIGHTS/tissue_list.txt', header = F)
LAA_FU_v8_tissue_list <- LAA_FU_v8_tissue$V1

LAA_FU_v8 <-{}
for (num in 1:22){
  #num <- 1
  chr_num <- paste0('chr', num)
  for (tissue in LAA_FU_v8_tissue_list){
    base_folder <- '/home/sysadmin/Desktop/dyc_lab/FUSION/RESULTS/TWBSTROKE/hg19/LAA_DR2_0.7/'
    twas <- fread(paste0(base_folder, chr_num, '/TWAS_', chr_num, '_DR207_TWBSTROKE_LAA_GTExv8_', tissue, '.dat'))
    #head(twas)
    twas$tissue <- str_remove(twas$PANEL, "^GTExv8.EUR\\.")
    twas <- twas[,-c(1,2)]
    LAA_FU_v8 <- rbind(LAA_FU_v8, twas)
  }
}
LAA_FU_v8 <- LAA_FU_v8[order(LAA_FU_v8$TWAS.P), ]
head(LAA_FU_v8)

# load LAA FUSION v7 result -----
LAA_FU_v7_tissue <- fread('/home/sysadmin/Desktop/dyc_lab/FUSION/WEIGHTS_v7/GTEx.ALL/tissue_list_v7.txt', header = F)
LAA_FU_v7_tissue_list <- LAA_FU_v7_tissue$V1

LAA_FU_v7 <-{}
for (num in 1:22){
  #num <- 1
  chr_num <- paste0('chr', num)
  for (tissue in LAA_FU_v7_tissue_list){
    #tissue <- tissue_list_v7[1]
    base_folder <- '/home/sysadmin/Desktop/dyc_lab/FUSION/RESULTS_v7/TWBSTROKE/hg19/LAA_DR2_0.7/'
    twas <- fread(paste0(base_folder, chr_num, '/TWAS_', chr_num, '_DR207_TWBSTROKE_LAA_GTExv7_', tissue, '.dat'))
    #head(twas)
    twas$tissue <- tissue
    twas <- twas[,-c(1,2)]
    LAA_FU_v7 <- rbind(LAA_FU_v7, twas)
  }
}
LAA_FU_v7 <- LAA_FU_v7[order(LAA_FU_v7$TWAS.P), ]


# load LAA S-PrediXcan v8 result -----
LAA_SP_v8_tissue <- fread('/home/sysadmin/Desktop/dyc_lab/S_Predixcan/tissues_v8_elastic_net.txt', header = F)
LAA_SP_v8_tissue_list <- LAA_SP_v8_tissue$V1

LAA_SP_v8 <- {}
for (tissue in LAA_SP_v8_tissue_list){
  twas <- fread(paste0('/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result_v8/',
                       tissue, '_TWB1_LAA_DR2_0.7_SPrediXcan_v8_en.csv'))
  twas$tissue <- tissue
  LAA_SP_v8 <- rbind(LAA_SP_v8, twas)
}
LAA_SP_v8 <- LAA_SP_v8[order(LAA_SP_v8$pvalue),]
head(LAA_SP_v8)
nrow(LAA_SP_v8)   #277998

# load LAA S-PrediXcan v7 result -----
LAA_SP_v7_tissue <- fread('/home/sysadmin/Desktop/dyc_lab/S_Predixcan/tissues.txt', header = F)
LAA_SP_v7_tissue_list <- LAA_SP_v7_tissue$V1

LAA_SP_v7 <- {}
for (tissue in LAA_SP_v7_tissue_list){
  twas <- fread(paste0('/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result/',
                       tissue, '_TWB1_LAA_DR2_0.7_SPrediXcan.csv'))
  twas$tissue <- tissue
  LAA_SP_v7 <- rbind(LAA_SP_v7, twas)
}
LAA_SP_v7 <- LAA_SP_v7[order(LAA_SP_v7$pvalue),]
head(LAA_sp_v7)
nrow(LAA_sp_v7)   #245256


# S-PrediXcan v8 gene annotation & manhattan plot -----
#LAA_SP_v8[, gene_no_version := gsub("\\..*", "", gene)]
head(LAA_SP_v8)
nrow(LAA_SP_v8)   #277998
colnames(LAA_SP_v8)
setnames(LAA_SP_v8, "gene", "ensembl_gene_id")

LAA_SP_v8_annot <- merge(
  LAA_SP_v8,
  gene_dt,
  by.x = c("ensembl_gene_id", "gene_name"),
  by.y = c("ensembl_gene_id", "gene_name")
)
head(LAA_SP_v8_annot)
nrow(LAA_SP_v8_annot) #277998
length(unique(LAA_SP_v8_annot$ensembl_gene_id)) #21583
length(unique(LAA_SP_v8$ensembl_gene_id)) #21583
sum(is.na(LAA_SP_v8_annot$start))

LAA_SP_v8_dup_symbol <- LAA_SP_v8[
  , uniqueN(ensembl_gene_id),
  by = gene_name
][V1 > 1]

LAA_SP_v8_annot_dup_symbol <- LAA_SP_v8_annot[
  , uniqueN(ensembl_gene_id),
  by = gene_name
][V1 > 1]

LAA_SP_v8_annot <- LAA_SP_v8_annot[order(LAA_SP_v8_annot$pvalue ),]
setnames(LAA_SP_v8_annot, "chr", "CHR")
LAA_SP_v8_annot[, BP := start]

#colnames(LAA_SP_v8_annot)
LAA_SP_v8_manhattan <- LAA_SP_v8_annot %>%
  filter(CHR %in% 1:22) %>%
  mutate(P = pvalue)
head(LAA_SP_v8_manhattan)

LAA_SP_v8_manhattan <- LAA_SP_v8_manhattan %>% filter(!is.na(P))
LAA_SP_v8_manhattan <- as.data.table(LAA_SP_v8_manhattan)
LAA_SP_v8_manhattan[, FDR := p.adjust(P, method = "BH")]
LAA_SP_v8_manhattan[, FDR_tissue := p.adjust(P, method = "BH"), by = tissue]
LAA_SP_v8_manhattan[LAA_SP_v8_manhattan$FDR_tissue <= 0.1, c("gene_name", "P", "FDR", "FDR_tissue", "tissue")]

write.csv(LAA_SP_v8_manhattan,
          '/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result_table/TWBSTROKE/LAA_DR2_0.7/LAA_DR2_0.7_SPrediXcan_v8_TWAS_20260302.csv',
          row.names = F)
write.csv(LAA_SP_v8_manhattan[LAA_SP_v8_manhattan$P < 1e-4,],
          '/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result_table/TWBSTROKE/LAA_DR2_0.7/LAA_DR2_0.7_SPrediXcan_v8_TWAS_1e4_20260302.csv',
          row.names = F)
write.csv(LAA_SP_v8_manhattan[LAA_SP_v8_manhattan$FDR_tissue < 0.2,],
          '/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result_table/TWBSTROKE/LAA_DR2_0.7/LAA_DR2_0.7_SPrediXcan_v8_TWAS_FDR02_20260307.csv',
          row.names = F)

LAA_SP_v8_manhattan_1e4 <- LAA_SP_v8_manhattan[LAA_SP_v8_manhattan$P < 1e-4,]
LAA_SP_v8_manhattan_1e4[LAA_SP_v8_manhattan_1e4$gene_name == 'TWIST1',]

LAA_SP_v8_manhattan_clean <- LAA_SP_v8_manhattan %>%
  dplyr::select(CHR, BP, P, FDR_tissue, gene_name, tissue)
LAA_SP_v8_tissue_df <- data.frame(table(LAA_SP_v8_manhattan_clean$tissue))
LAA_SP_v8_tissue_df$proportion <- LAA_SP_v8_tissue_df$Freq/sum(LAA_SP_v8_tissue_df$Freq)
setnames(LAA_SP_v8_tissue_df, c('tissue', 'SP_v8_freq', 'SP_v8_proportion'))

LAA_SP_v8_manhattan_clean <- LAA_SP_v8_manhattan_clean %>%
  filter(
    !is.na(CHR),
    !is.na(BP)
  )

LAA_SP_v8_plot_df <- LAA_SP_v8_manhattan_clean %>%
  mutate(
    CHR = as.numeric(CHR),
    BP  = as.numeric(BP)
  ) %>%
  arrange(CHR, BP)

LAA_SP_v8_chr_max <- LAA_SP_v8_plot_df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

LAA_SP_v8_plot_df <- LAA_SP_v8_plot_df %>%
  left_join(LAA_SP_v8_chr_max, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

LAA_SP_v8_chr_centers <- LAA_SP_v8_plot_df %>%
  group_by(CHR) %>%
  summarise(center = mean(range(BP_cum)))

LAA_SP_v8_gene_1e4 <- LAA_SP_v8_plot_df[LAA_SP_v8_plot_df$P < 1e-4,]
LAA_SP_v8_gene_fdr <- LAA_SP_v8_plot_df[LAA_SP_v8_plot_df$FDR_tissue < 0.2,]
#LAA_SP_v8_label_genes <- unique(LAA_SP_v8_gene_1e4$gene_name)
LAA_SP_v8_label_genes <- unique(LAA_SP_v8_gene_fdr$gene_name)
LAA_SP_v8_label_df <- LAA_SP_v8_plot_df %>%
  filter(
    gene_name %in% LAA_SP_v8_label_genes,
    FDR_tissue < 0.2
  ) %>%
  group_by(gene_name) %>%
  slice_min(P, n = 1) %>%
  ungroup()

LAA_SP_v8_plot_df <- LAA_SP_v8_plot_df %>%
  mutate(chr_color = factor(CHR %% 2))

#c("grey60", "black") , values = c("#2C7BB6", "#F28E2B")
LAA_SP_v8_p <- ggplot(LAA_SP_v8_plot_df, aes(x = BP_cum, y = -log10(P))) +
  geom_point(aes(color = chr_color), size = 0.6, alpha = 0.8) +
  scale_color_manual(values = c("#2C7BB6", "#F28E2B"), guide = "none") +
  #geom_hline(yintercept = 4, color = "blue", linewidth = 0.6) +
  geom_point(
    data = LAA_SP_v8_label_df,
    aes(x = BP_cum, y = -log10(P)),
    color = "red",
    size = 2
  ) +
  geom_text_repel(
    data = LAA_SP_v8_label_df,
    aes(label = gene_name),
    color = "red",
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  labs(
    title = "S-PrediXcan GTEx_v8",
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_continuous(
    breaks = LAA_SP_v8_chr_centers$center,
    labels = LAA_SP_v8_chr_centers$CHR
  ) + scale_y_continuous(
    limits = c(0, 9),
    breaks = seq(0, 9, by = 1)
  )

ggsave("/home/sysadmin/Desktop/dyc_lab/S_Predixcan/manhantta_ plot/SPredixcan_GTEx_v8_manhanttan_fdr02_20260306.png",
       LAA_SP_v8_p, width = 12, height = 6, dpi = 600)

# S-PrediXcan v7 gene annotation & manhattan plot -----
#LAA_SP_v7[, gene_no_version := gsub("\\..*", "", gene)]
head(LAA_SP_v7)
nrow(LAA_SP_v7)   #245256
colnames(LAA_SP_v7)
setnames(LAA_SP_v7, "gene", "ensembl_gene_id")

LAA_SP_v7_annot <- merge(
  LAA_SP_v7,
  gene_dt_v19,
  by.x = c("ensembl_gene_id", "gene_name"),
  by.y = c("ensembl_gene_id", "gene_name")
)
head(LAA_SP_v7_annot)
nrow(LAA_SP_v7_annot) #245256
length(unique(LAA_SP_v7_annot$ensembl_gene_id)) #25720
length(unique(LAA_SP_v7$ensembl_gene_id)) #25720
sum(is.na(LAA_SP_v7_annot$start))

LAA_SP_v7_dup_symbol <- LAA_SP_v7[
  , uniqueN(ensembl_gene_id),
  by = gene_name
][V1 > 1]

LAA_SP_v7_annot_dup_symbol <- LAA_SP_v7_annot[
  , uniqueN(ensembl_gene_id),
  by = gene_name
][V1 > 1]

LAA_SP_v7_annot <- LAA_SP_v7_annot[order(LAA_SP_v7_annot$pvalue ),]
setnames(LAA_SP_v7_annot, "chr", "CHR")
LAA_SP_v7_annot[, BP := start]

#colnames(LAA_SP_v7_annot)
LAA_SP_v7_manhattan <- LAA_SP_v7_annot %>%
  filter(CHR %in% 1:22) %>%
  mutate(P = pvalue)
head(LAA_SP_v7_manhattan)

LAA_SP_v7_manhattan <- LAA_SP_v7_manhattan %>% filter(!is.na(P))
LAA_SP_v7_manhattan <- as.data.table(LAA_SP_v7_manhattan)
LAA_SP_v7_manhattan[, FDR := p.adjust(P, method = "BH")]
LAA_SP_v7_manhattan[, FDR_tissue := p.adjust(P, method = "BH"), by = tissue]
LAA_SP_v7_manhattan[LAA_SP_v7_manhattan$FDR_tissue <= 0.1, c("gene_name", "P", "FDR", "FDR_tissue", "tissue")]

write.csv(LAA_SP_v7_manhattan,
          '/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result_table/TWBSTROKE/LAA_DR2_0.7/LAA_DR2_0.7_SPrediXcan_v7_TWAS_20260302.csv',
          row.names = F)
write.csv(LAA_SP_v7_manhattan[LAA_SP_v7_manhattan$P < 1e-4,],
          '/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result_table/TWBSTROKE/LAA_DR2_0.7/LAA_DR2_0.7_SPrediXcan_v7_TWAS_1e4_20260302.csv',
          row.names = F)
write.csv(LAA_SP_v7_manhattan[LAA_SP_v7_manhattan$FDR_tissue < 0.2,],
          '/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result_table/TWBSTROKE/LAA_DR2_0.7/LAA_DR2_0.7_SPrediXcan_v7_TWAS_FDR02_20260307.csv',
          row.names = F)

LAA_SP_v7_manhattan_1e4 <- LAA_SP_v7_manhattan[LAA_SP_v7_manhattan$P < 1e-4,]
LAA_SP_v7_manhattan_1e4[LAA_SP_v7_manhattan_1e4$gene_name == 'ITGAV',]

LAA_SP_v7_manhattan_clean <- LAA_SP_v7_manhattan %>%
  dplyr::select(CHR, BP, P, FDR_tissue, gene_name, tissue)
LAA_SP_v7_tissue_df <- data.frame(table(LAA_SP_v7_manhattan_clean$tissue))
LAA_SP_v7_tissue_df$proportion <- LAA_SP_v7_tissue_df$Freq/sum(LAA_SP_v7_tissue_df$Freq)
setnames(LAA_SP_v7_tissue_df, c('tissue', 'SP_v7_freq', 'SP_v7_proportion'))

LAA_SP_v7_manhattan_clean <- LAA_SP_v7_manhattan_clean %>%
  filter(
    !is.na(CHR),
    !is.na(BP)
  )

LAA_SP_v7_plot_df <- LAA_SP_v7_manhattan_clean %>%
  mutate(
    CHR = as.numeric(CHR),
    BP  = as.numeric(BP)
  ) %>%
  arrange(CHR, BP)


LAA_SP_v7_chr_max <- LAA_SP_v7_plot_df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

LAA_SP_v7_plot_df <- LAA_SP_v7_plot_df %>%
  left_join(LAA_SP_v7_chr_max, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

LAA_SP_v7_chr_centers <- LAA_SP_v7_plot_df %>%
  group_by(CHR) %>%
  summarise(center = mean(range(BP_cum)))

LAA_SP_v7_gene_1e4 <- LAA_SP_v7_plot_df[LAA_SP_v7_plot_df$P < 1e-4,]
LAA_SP_v7_gene_fdr <- LAA_SP_v7_plot_df[LAA_SP_v7_plot_df$FDR_tissue < 0.2,]
#LAA_SP_v7_label_genes <- unique(LAA_SP_v7_gene_1e4$gene_name)
LAA_SP_v7_label_genes <- unique(LAA_SP_v7_gene_fdr$gene_name)
LAA_SP_v7_label_df <- LAA_SP_v7_plot_df %>%
  filter(
    gene_name %in% LAA_SP_v7_label_genes,
    FDR_tissue < 0.2
  ) %>%
  group_by(gene_name) %>%
  slice_min(P, n = 1) %>%
  ungroup()

LAA_SP_v7_plot_df <- LAA_SP_v7_plot_df %>%
  mutate(chr_color = factor(CHR %% 2))

#c("grey60", "black") , values = c("#2C7BB6", "#F28E2B")
LAA_SP_v7_p <- ggplot(LAA_SP_v7_plot_df, aes(x = BP_cum, y = -log10(P))) +
  geom_point(aes(color = chr_color), size = 0.6, alpha = 0.8) +
  scale_color_manual(values = c("#2C7BB6", "#F28E2B"), guide = "none") +
  #geom_hline(yintercept = 4, color = "blue", linewidth = 0.6) +
  geom_point(
    data = LAA_SP_v7_label_df,
    aes(x = BP_cum, y = -log10(P)),
    color = "red",
    size = 2
  ) +
  geom_text_repel(
    data = LAA_SP_v7_label_df,
    aes(label = gene_name),
    color = "red",
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  labs(
    title = "S-PrediXcan GTEx_v7",
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_continuous(
    breaks = LAA_SP_v7_chr_centers$center,
    labels = LAA_SP_v7_chr_centers$CHR
  ) + scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, by = 1)
  )

ggsave("/home/sysadmin/Desktop/dyc_lab/S_Predixcan/manhantta_ plot/SPredixcan_GTEx_v7_manhanttan_fdr02_20260306.png",
       LAA_SP_v7_p, width = 12, height = 6, dpi = 600)


# FUSION v8 gene annotation & manhattan plot -----
head(LAA_FU_v8)
nrow(LAA_FU_v8)   #295077
colnames(LAA_FU_v8)
setnames(LAA_FU_v8, "ID", "ensembl_gene_id")
LAA_FU_v8[, CHR := as.character(CHR)]

LAA_FU_v8_annot <- merge(
  LAA_FU_v8,
  gene_dt,
  by.x = c("ensembl_gene_id", "CHR"),
  by.y = c("ensembl_gene_id", "chr")
)
head(LAA_FU_v8_annot)
nrow(LAA_FU_v8_annot) #295077
length(unique(LAA_FU_v8_annot$ensembl_gene_id)) #27671
length(unique(LAA_FU_v8$ensembl_gene_id)) #27671
sum(is.na(LAA_FU_v8_annot$start))

LAA_FU_v8_annot_dup_symbol <- LAA_FU_v8_annot[
  , uniqueN(ensembl_gene_id),
  by = gene_name
][V1 > 1]

LAA_FU_v8_annot <- LAA_FU_v8_annot[order(LAA_FU_v8_annot$TWAS.P ),]
LAA_FU_v8_annot[, BP := start]

#colnames(LAA_SP_v8_annot)
LAA_FU_v8_manhattan <- LAA_FU_v8_annot %>%
  filter(CHR %in% 1:22) %>%
  mutate(P = TWAS.P)
head(LAA_FU_v8_manhattan)

LAA_FU_v8_manhattan <- LAA_FU_v8_manhattan %>% filter(!is.na(P))
LAA_FU_v8_manhattan <- as.data.table(LAA_FU_v8_manhattan)
LAA_FU_v8_manhattan[, FDR := p.adjust(P, method = "BH")]
LAA_FU_v8_manhattan[, FDR_tissue := p.adjust(P, method = "BH"), by = tissue]
LAA_FU_v8_manhattan[LAA_FU_v8_manhattan$FDR_tissue <= 0.1, c("gene_name", "P", "FDR", "FDR_tissue", "tissue")]

write.csv(LAA_FU_v8_manhattan,
          '/home/sysadmin/Desktop/dyc_lab/FUSION/result_table/GTEx_v8/TWBSTROKE/hg19_annotation/LAA_DR2_0.7/LAA_DR2_0.7_FUSION_v8_TWAS_20260302.csv',
          row.names = F)

write.csv(LAA_FU_v8_manhattan[LAA_FU_v8_manhattan$P < 1e-4,],
          '/home/sysadmin/Desktop/dyc_lab/FUSION/result_table/GTEx_v8/TWBSTROKE/hg19_annotation/LAA_DR2_0.7/LAA_DR2_0.7_FUSION_v8_TWAS_1e4_20260302.csv',
          row.names = F)

write.csv(LAA_FU_v8_manhattan[LAA_FU_v8_manhattan$FDR_tissue < 0.2,],
          '/home/sysadmin/Desktop/dyc_lab/FUSION/result_table/GTEx_v8/TWBSTROKE/hg19_annotation/LAA_DR2_0.7/LAA_DR2_0.7_FUSION_v8_TWAS_FDR02_20260307.csv',
          row.names = F)

LAA_FU_v8_manhattan_1e4 <- LAA_FU_v8_manhattan[LAA_FU_v8_manhattan$P < 1e-4,]
LAA_FU_v8_manhattan_1e4[LAA_FU_v8_manhattan_1e4$gene_name == 'DAPK2']

LAA_FU_v8_manhattan_clean <- LAA_FU_v8_manhattan %>%
  dplyr::select(CHR, BP, P, FDR_tissue, gene_name, tissue)
LAA_FU_v8_tissue_df <- data.frame(table(LAA_FU_v8_manhattan_clean$tissue))
LAA_FU_v8_tissue_df$proportion <- LAA_FU_v8_tissue_df$Freq/sum(LAA_FU_v8_tissue_df$Freq)
setnames(LAA_FU_v8_tissue_df, c('tissue', 'FU_v8_freq', 'FU_v8_proportion'))

LAA_FU_v8_manhattan_clean <- LAA_FU_v8_manhattan_clean %>%
  filter(
    !is.na(CHR),
    !is.na(BP),
  )

LAA_FU_v8_plot_df <- LAA_FU_v8_manhattan_clean %>%
  mutate(
    CHR = as.numeric(CHR),
    BP  = as.numeric(BP)
  ) %>%
  arrange(CHR, BP)

LAA_FU_v8_chr_max <- LAA_FU_v8_plot_df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

LAA_FU_v8_plot_df <- LAA_FU_v8_plot_df %>%
  left_join(LAA_FU_v8_chr_max, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

LAA_FU_v8_chr_centers <- LAA_FU_v8_plot_df %>%
  group_by(CHR) %>%
  summarise(center = mean(range(BP_cum)))

LAA_FU_v8_gene_1e4 <- LAA_FU_v8_plot_df[LAA_FU_v8_plot_df$P < 1e-4,]
LAA_FU_v8_gene_fdr <- LAA_FU_v8_plot_df[LAA_FU_v8_plot_df$FDR_tissue < 0.2,]
#LAA_FU_v8_label_genes <- unique(LAA_FU_v8_gene_1e4$gene_name)
LAA_FU_v8_label_genes <- unique(LAA_FU_v8_gene_fdr$gene_name)
LAA_FU_v8_label_df <- LAA_FU_v8_plot_df %>%
  filter(
    gene_name %in% LAA_FU_v8_label_genes,
    FDR_tissue < 0.2
  ) %>%
  group_by(gene_name) %>%
  slice_min(P, n = 1) %>%
  ungroup()

LAA_FU_v8_plot_df <- LAA_FU_v8_plot_df %>%
  mutate(chr_color = factor(CHR %% 2))

#c("grey60", "black") , values = c("#2C7BB6", "#F28E2B")
LAA_FU_v8_p <- ggplot(LAA_FU_v8_plot_df, aes(x = BP_cum, y = -log10(P))) +
  geom_point(aes(color = chr_color), size = 0.6, alpha = 0.8) +
  scale_color_manual(values = c("#2C7BB6", "#F28E2B"), guide = "none") +
  #geom_hline(yintercept = 4, color = "blue", linewidth = 0.6) +
  geom_point(
    data = LAA_FU_v8_label_df,
    aes(x = BP_cum, y = -log10(P)),
    color = "red",
    size = 2
  ) +
  geom_text_repel(
    data = LAA_FU_v8_label_df,
    aes(label = gene_name),
    color = "red",
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  labs(
    title = "FUSION GTEx_v8",
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_continuous(
    breaks = LAA_FU_v8_chr_centers$center,
    labels = LAA_FU_v8_chr_centers$CHR
  ) + scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, by = 1)
  )

ggsave("/home/sysadmin/Desktop/dyc_lab/FUSION/manhanttan plot/FUSION_GTEx_v8_manhanttan_fdr02_20260306.png",
       LAA_FU_v8_p, width = 12, height = 6, dpi = 600)


# FUSION v7 gene annotation & manhattan plot -----
colnames(LAA_FU_v7)
head(gene_dt_v19)

head(LAA_FU_v7)
nrow(LAA_FU_v7)   #81826
setnames(LAA_FU_v7, "ID", "gene_name")
LAA_FU_v7[, CHR := as.character(CHR)]

LAA_FU_v7_annot <- merge(
  LAA_FU_v7,
  gene_map_wgt[, .(gene_name, ensembl_gene_id)],
  by="gene_name",
  all.x=TRUE
)

LAA_FU_v7_annot <- merge(
  LAA_FU_v7_annot,
  gene_dt_v19,
  by.x = c("ensembl_gene_id", "gene_name", "CHR"),
  by.y = c("ensembl_gene_id", "gene_name", "chr")
)

LAA_FU_v7_annot <- LAA_FU_v7_annot[
    P0 >= start & 
    P1 <= end
]

LAA_FU_v7_dup_groups <- LAA_FU_v7_annot[
  , .N, 
  by = .(gene_name, CHR, P0, P1, tissue)
][N > 1]

LAA_FU_v7_annot_clean <- LAA_FU_v7_annot[
  !LAA_FU_v7_annot[
    , .I[.N > 1], 
    by = .(gene_name, CHR, P0, P1, tissue)
  ]$V1
]
nrow(LAA_FU_v7_annot_clean)
nrow(LAA_FU_v7_annot)
nrow(LAA_FU_v7)

head(LAA_FU_v7_annot_clean)
length(unique(LAA_FU_v7_annot_clean$gene_name)) #12001
length(unique(LAA_FU_v7$gene_name)) #12002
sum(is.na(LAA_FU_v7_annot_clean$start))

LAA_FU_v7_annot_dup_symbol <- LAA_FU_v7_annot_clean[
  , uniqueN(ensembl_gene_id),
  by = gene_name
][V1 > 1]
LAA_FU_v7_annot_dup_symbol

LAA_FU_v7_annot_clean <- LAA_FU_v7_annot_clean[order(LAA_FU_v7_annot_clean$TWAS.P ),]
LAA_FU_v7_annot_clean[, BP := start]

LAA_FU_v7_manhattan <- LAA_FU_v7_annot_clean %>%
  filter(CHR %in% 1:22) %>%
  mutate(P = TWAS.P)
head(LAA_FU_v7_manhattan)

LAA_FU_v7_manhattan <- LAA_FU_v7_manhattan %>% filter(!is.na(P))
LAA_FU_v7_manhattan <- as.data.table(LAA_FU_v7_manhattan)
LAA_FU_v7_manhattan[, FDR := p.adjust(P, method = "BH")]
LAA_FU_v7_manhattan[, FDR_tissue := p.adjust(P, method = "BH"), by = tissue]
LAA_FU_v7_manhattan[LAA_FU_v7_manhattan$FDR_tissue <= 0.15, c("gene_name", "P", "FDR", "FDR_tissue", "tissue")]

write.csv(LAA_FU_v7_manhattan,
          '/home/sysadmin/Desktop/dyc_lab/FUSION/result_table/GTEx_v7/TWBSTROKE/LAA_DR2_0.7/LAA_DR2_0.7_FUSION_v7_TWAS_20260302.csv',
          row.names = F)

write.csv(LAA_FU_v7_manhattan[LAA_FU_v7_manhattan$P < 1e-4,],
          '/home/sysadmin/Desktop/dyc_lab/FUSION/result_table/GTEx_v7/TWBSTROKE/LAA_DR2_0.7/LAA_DR2_0.7_FUSION_v7_TWAS_1e4_20260302.csv',
          row.names = F)

write.csv(LAA_FU_v7_manhattan[LAA_FU_v7_manhattan$FDR_tissue < 0.2,],
          '/home/sysadmin/Desktop/dyc_lab/FUSION/result_table/GTEx_v7/TWBSTROKE/LAA_DR2_0.7/LAA_DR2_0.7_FUSION_v7_TWAS_FDR02_20260307.csv',
          row.names = F)

LAA_FU_v7_manhattan_1e4 <- LAA_FU_v7_manhattan[LAA_FU_v7_manhattan$P < 1e-4,]
LAA_FU_v7_manhattan_1e4[LAA_FU_v7_manhattan_1e4$gene_name == 'DAPK2',]

LAA_FU_v7_manhattan_clean <- LAA_FU_v7_manhattan %>%
  dplyr::select(CHR, BP, P, FDR_tissue, gene_name, tissue) #FDR
LAA_FU_v7_tissue_df <- data.frame(table(LAA_FU_v7_manhattan_clean$tissue))
LAA_FU_v7_tissue_df$proportion <- LAA_FU_v7_tissue_df$Freq/sum(LAA_FU_v7_tissue_df$Freq)
setnames(LAA_FU_v7_tissue_df, c('tissue', 'FU_v7_freq', 'FU_v7_proportion'))

LAA_FU_v7_manhattan_clean <- LAA_FU_v7_manhattan_clean %>%
  filter(
    !is.na(CHR),
    !is.na(BP)
  )

LAA_FU_v7_plot_df <- LAA_FU_v7_manhattan_clean %>%
  mutate(
    CHR = as.numeric(CHR),
    BP  = as.numeric(BP)
  ) %>%
  arrange(CHR, BP)

LAA_FU_v7_chr_max <- LAA_FU_v7_plot_df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

LAA_FU_v7_plot_df <- LAA_FU_v7_plot_df %>%
  left_join(LAA_FU_v7_chr_max, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

LAA_FU_v7_chr_centers <- LAA_FU_v7_plot_df %>%
  group_by(CHR) %>%
  summarise(center = mean(range(BP_cum)))

LAA_FU_v7_gene_1e4 <- LAA_FU_v7_plot_df[LAA_FU_v7_plot_df$P < 1e-4,]
LAA_FU_v7_gene_fdr <- LAA_FU_v7_plot_df[LAA_FU_v7_plot_df$FDR_tissue < 0.2,]
#LAA_FU_v7_label_genes <- unique(LAA_FU_v7_gene_1e4$gene_name)
LAA_FU_v7_label_genes <- unique(LAA_FU_v7_gene_fdr$gene_name)
LAA_FU_v7_label_df <- LAA_FU_v7_plot_df %>%
  filter(
    gene_name %in% LAA_FU_v7_label_genes,
    FDR_tissue < 0.2
  ) %>%
  group_by(gene_name) %>%
  slice_min(P, n = 1) %>%
  ungroup()

LAA_FU_v7_plot_df <- LAA_FU_v7_plot_df %>%
  mutate(chr_color = factor(CHR %% 2))

#c("grey60", "black") , values = c("#2C7BB6", "#F28E2B")
LAA_FU_v7_p <- ggplot(LAA_FU_v7_plot_df, aes(x = BP_cum, y = -log10(P))) +
  geom_point(aes(color = chr_color), size = 0.6, alpha = 0.8) +
  scale_color_manual(values = c("#2C7BB6", "#F28E2B"), guide = "none") +
  #geom_hline(yintercept = 4, color = "blue", linewidth = 0.6) +
  geom_point(
    data = LAA_FU_v7_label_df,
    aes(x = BP_cum, y = -log10(P)),
    color = "red",
    size = 2
  ) +
  geom_text_repel(
    data = LAA_FU_v7_label_df,
    aes(label = gene_name),
    color = "red",
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  labs(
    title = "FUSION GTEx_v7",
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_continuous(
    breaks = LAA_FU_v7_chr_centers$center,
    labels = LAA_FU_v7_chr_centers$CHR
  ) + scale_y_continuous(
    limits = c(0, 6),
    breaks = seq(0, 6, by = 1)
  )

ggsave("/home/sysadmin/Desktop/dyc_lab/FUSION/manhanttan plot/FUSION_GTEx_v7_manhanttan_fdr02_20260306.png",
       LAA_FU_v7_p, width = 12, height = 6, dpi = 600)


# TWAS heatmap -----
library(tibble)
library(ggrepel)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

LAA_FU_v8_fdr02 <- LAA_FU_v8_manhattan %>%
  filter(FDR_tissue < 0.2)
LAA_FU_v7_fdr02 <- LAA_FU_v7_manhattan %>%
  filter(FDR_tissue < 0.2)
LAA_SP_v8_fdr02 <- LAA_SP_v8_manhattan %>%
  filter(FDR_tissue < 0.2)
LAA_SP_v7_fdr02 <- LAA_SP_v7_manhattan %>%
  filter(FDR_tissue < 0.2)


heat_df <- LAA_SP_v7_fdr02 %>%
  mutate(sig_value = -log10(FDR_tissue))

heat_mat <- heat_df %>%
  dplyr::select(gene_name, tissue, sig_value) %>%
  pivot_wider(
    names_from  = tissue,
    values_from = sig_value,
    values_fill = 0
  ) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

sig_cutoff <- -log10(0.2)   # = P < 1e-4
gene_count <- colSums(heat_mat > sig_cutoff)
tissue_count <- rowSums(heat_mat > sig_cutoff)

annotation_col <- data.frame(
  GeneCount = gene_count
)
rownames(annotation_col) <- colnames(heat_mat)

annotation_row <- data.frame(
  TissueCount = tissue_count
)
rownames(annotation_row) <- rownames(heat_mat)


colnames(heat_mat)
tissue_order <- order(colnames(heat_mat))
heat_mat <- heat_mat[, tissue_order]
gene_count <- gene_count[tissue_order]

min(heat_mat[heat_mat > 0])
median(heat_mat[heat_mat > 0])
col_fun <- colorRamp2(
  c(
    0,
    median(heat_mat[heat_mat > 0]),
    max(heat_mat)
  ),
  c('grey', "#6baed6", "#d73027")
)

ht <- Heatmap(
  heat_mat,
  name = "Significance",
  col  = col_fun,
  
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  show_row_dend   = FALSE,
  show_column_dend= FALSE,
  
  show_row_names    = TRUE,
  show_column_names = TRUE,
  
  row_names_gp    = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 60,   # ??? ?o?@??
  
  top_annotation  = HeatmapAnnotation(
    GeneCount = anno_barplot(
      gene_count,
      gp = gpar(fill = "grey40"),
      height = unit(2, "cm")
    )
  ),
  
  left_annotation = rowAnnotation(
    TissueCount = anno_barplot(
      tissue_count,
      gp = gpar(fill = "grey40"),
      width = unit(2, "cm")
    )
  )
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

png(
  filename = "/home/sysadmin/Desktop/dyc_lab/FUSION/heatmap/TWAS_heatmap_FUSION_GTEx_v7_20260306.png",
  width    = 4200,
  height   = 3000,
  res      = 300
)

png(
  filename = "/home/sysadmin/Desktop/dyc_lab/S_Predixcan/heatmap/TWAS_heatmap_SPredixcan_GTEx_v7_20260306.png",
  width    = 4200,
  height   = 3000,
  res      = 300
)

draw(
  ht,
  heatmap_legend_side    = "right",
  annotation_legend_side = "right"
)

dev.off()


# TWAS tissue summary -----
#install.packages("pheatmap")
library(pheatmap)
sum(LAA_SP_v7_tissue_df$SP_v7_freq)
sum(LAA_SP_v8_tissue_df$SP_v8_freq)

sum(LAA_FU_v7_tissue_df$FU_v7_freq)
sum(LAA_FU_v8_tissue_df$FU_v8_freq)

LAA_tissue_df <- merge(LAA_FU_v8_tissue_df, LAA_SP_v8_tissue_df, by='tissue', all=T)
LAA_tissue_df <- merge(LAA_tissue_df, LAA_FU_v7_tissue_df, by='tissue', all=T)
LAA_tissue_df <- merge(LAA_tissue_df, LAA_SP_v7_tissue_df, by='tissue', all=T)
LAA_tissue_df$all_freq <- LAA_tissue_df$FU_v8_freq + LAA_tissue_df$SP_v8_freq + 
  LAA_tissue_df$FU_v7_freq + LAA_tissue_df$SP_v7_freq
head(LAA_tissue_df)



heat_mat <- LAA_tissue_df %>%
  column_to_rownames("tissue") %>%
  dplyr::select(
    FU_v8_freq,
    SP_v8_freq,
    FU_v7_freq,
    SP_v7_freq
  ) %>%
  as.matrix()
heat_mat[is.na(heat_mat)] <- 0
colnames(heat_mat) <- c(
  "FUSION v8",
  "SPrediXcan v8",
  "FUSION v7",
  "SPrediXcan v7"
)

file.create("/mnt/data/stroke_pipeline_test/tissue_heatmap/test.txt")

png("/mnt/data/stroke_pipeline_test/twas_tissue_freq_heatmap_20260306.png",
    width = 5000,
    height = 5000,
    res = 600)

pheatmap(
  heat_mat,
  cluster_rows = FALSE,     
  cluster_cols = FALSE, 
  angle_col = 0,
  color = colorRampPalette(
    c("white", "#d0e1f2", "#74a9cf", "#2b8cbe", "#045a8d")
  )(100),                   
  main = "Number of Predictive Genes Across Tissues"
)

dev.off()

# FDR correction -----
LAA_anno
LAA_FU_v7_manhattan
LAA_FU_v8_manhattan
LAA_SP_v7_manhattan
LAA_SP_v8_manhattan

LAA_FU_v7_FDR <- LAA_FU_v7_manhattan[,c("ensembl_gene_id", "gene_name", "CHR", "P", "FDR",  "FDR_tissue", "tissue")]
LAA_FU_v8_FDR <- LAA_FU_v8_manhattan[,c("ensembl_gene_id", "gene_name", "CHR", "P", "FDR",  "FDR_tissue", "tissue")]
LAA_SP_v7_FDR <- LAA_SP_v7_manhattan[,c("ensembl_gene_id", "gene_name", "CHR", "P", "FDR",  "FDR_tissue", "tissue")]
LAA_SP_v8_FDR <- LAA_SP_v8_manhattan[,c("ensembl_gene_id", "gene_name", "CHR", "P", "FDR",  "FDR_tissue", "tissue")]

head(LAA_FU_v7_FDR)
head(LAA_FU_v8_FDR)
head(LAA_SP_v7_FDR)
head(LAA_SP_v8_FDR)
head(LAA_anno)

LAA_SP_v8_FDR[LAA_SP_v8_FDR$FDR < 0.05,]
LAA_SP_v8_FDR[LAA_SP_v8_FDR$FDR_tissue < 0.1,]
max(LAA_SP_v8_FDR$P[LAA_SP_v8_FDR$FDR_tissue < 0.1])
max(table(LAA_SP_v8_FDR$tissue))

#nrow(LAA_SP_v8_FDR[LAA_SP_v8_FDR$tissue=="Whole_Blood",])
#(3/7200)*0.1


# genes overlap by different methods (venn plot) -----
strip_version <- function(x) sub("\\..*", "", x)

get_fu_set <- function(threshold){
  fu_v7 <- unique(strip_version(
    LAA_FU_v7_FDR[FDR_tissue <= threshold, ensembl_gene_id]
  ))
  fu_v8 <- unique(strip_version(
    LAA_FU_v8_FDR[FDR_tissue <= threshold, ensembl_gene_id]
  ))
  unique(c(fu_v7, fu_v8))
}

get_sp_set <- function(threshold){
  sp_v7 <- unique(strip_version(
    LAA_SP_v7_FDR[FDR_tissue <= threshold, ensembl_gene_id]
  ))
  sp_v8 <- unique(strip_version(
    LAA_SP_v8_FDR[FDR_tissue <= threshold, ensembl_gene_id]
  ))
  unique(c(sp_v7, sp_v8))
}

fine_map_genes <- unique(strip_version(
  LAA_anno$ensembl_gene_id_version
))

get_overlap_all <- function(fu, sp, gwas){
  list(
    FU_SP     = intersect(fu, sp),
    FU_GWAS   = intersect(fu, gwas),
    SP_GWAS   = intersect(sp, gwas),
    FU_SP_GWAS = Reduce(intersect, list(fu, sp, gwas))
  )
}

thresholds <- c(0.10, 0.15, 0.20, 0.25)

all_results <- list()
all_gene_lists <- list()

for(t in thresholds){
  fu_set <- get_fu_set(t)
  sp_set <- get_sp_set(t)
  overlap <- get_overlap_all(fu_set, sp_set, fine_map_genes)
  
  all_results[[as.character(t)]] <- data.table(
    FDR_tissue = t,
    FU_n   = length(fu_set),
    SP_n   = length(sp_set),
    GWAS_n = length(fine_map_genes),
    FU_SP  = length(overlap$FU_SP),
    FU_GWAS = length(overlap$FU_GWAS),
    SP_GWAS = length(overlap$SP_GWAS),
    FU_SP_GWAS = length(overlap$FU_SP_GWAS)
  )
  
  all_gene_lists[[as.character(t)]] <- overlap
}

final_overlap_table <- rbindlist(all_results)
final_overlap_table


FU_SP_GWAS_genes <- all_gene_lists[["0.2"]]$FU_SP_GWAS
FU_SP_genes <- all_gene_lists[["0.2"]]$FU_SP
FU_GWAS_genes <- all_gene_lists[["0.2"]]$FU_GWAS
SP_GWAS_genes <- all_gene_lists[["0.2"]]$SP_GWAS


LAA_FU_gene_map <- rbind(LAA_FU_v7_FDR[,c('ensembl_gene_id', 'gene_name')], LAA_FU_v8_FDR[,c('ensembl_gene_id', 'gene_name')])
LAA_SP_gene_map <- rbind(LAA_SP_v7_FDR[,c('ensembl_gene_id', 'gene_name')], LAA_SP_v8_FDR[,c('ensembl_gene_id', 'gene_name')])
LAA_FU_gene_map$ensembl_gene_id_clean <- sub("\\..*", "", LAA_FU_gene_map$ensembl_gene_id)
LAA_SP_gene_map$ensembl_gene_id_clean <- sub("\\..*", "", LAA_SP_gene_map$ensembl_gene_id)
LAA_FU_gene_map <- unique(LAA_FU_gene_map[,c('ensembl_gene_id_clean', 'gene_name')]) #28066
LAA_SP_gene_map <- unique(LAA_SP_gene_map[,c('ensembl_gene_id_clean', 'gene_name')]) #28705

# FU & SP & COJO
LAA_FU_gene_map[ensembl_gene_id_clean %in% FU_SP_GWAS_genes]
# FU & SP 
LAA_FU_gene_map[ensembl_gene_id_clean %in% FU_SP_genes]
# FU & COJO
LAA_FU_gene_map[ensembl_gene_id_clean %in% FU_GWAS_genes]
# SP & COJO
LAA_SP_gene_map[ensembl_gene_id_clean %in% SP_GWAS_genes]


t <- 0.2
venn_fdr_input <- list(
  FUSION = get_fu_set(t),
  SPrediXcan = get_sp_set(t),
  GWAS = fine_map_genes
)

png(
  filename = "/mnt/data/stroke_pipeline_test/venn_plot/COJO_FUSION_SPrediXcan_FDR02_venn.png",
  width = 3000,
  height = 3000,
  res = 300
)

ggvenn(
  venn_fdr_input,
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"),
  stroke_size = 0.5,
  set_name_size = 4
)

dev.off()




fu_v7_1e4 <- unique(strip_version(
  LAA_FU_v7_FDR[P <= 1e-4, ensembl_gene_id]
))

fu_v8_1e4 <- unique(strip_version(
  LAA_FU_v8_FDR[P <= 1e-4, ensembl_gene_id]
))

sp_v7_1e4 <- unique(strip_version(
  LAA_SP_v7_FDR[P <= 1e-4, ensembl_gene_id]
))

sp_v8_1e4 <- unique(strip_version(
  LAA_SP_v8_FDR[P <= 1e-4, ensembl_gene_id]
))

venn_1e4_input <- list(
  FUSION = unique(c(fu_v7_1e4, fu_v8_1e4)),
  SPrediXcan = unique(c(sp_v7_1e4, sp_v8_1e4)),
  COJO = fine_map_genes
)

png(
  filename = "/mnt/data/stroke_pipeline_test/venn_plot/COJO_FUSION_SPrediXcan_1e4_venn_20260302.png",
  width = 3000,
  height = 3000,
  res = 300
)

ggvenn(
  venn_1e4_input,
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"),
  stroke_size = 0.5,
  set_name_size = 4
)

dev.off()


intersect_all <- Reduce(intersect, venn_1e4_input)
intersect_all
length(intersect_all)

intersect_FU_SP <- intersect(venn_1e4_input$FUSION,
                             venn_1e4_input$SPrediXcan)
intersect_FU_COJO <- intersect(venn_1e4_input$FUSION,
                               venn_1e4_input$COJO)
intersect_SP_COJO <- intersect(venn_1e4_input$SPrediXcan,
                               venn_1e4_input$COJO)





# LAA overlap genes pathway analysis -----
LAA_gene <- unique(rbind(LAA_FU_gene_map, LAA_SP_gene_map))
# FU & SP 
intersect_FU_SP_gene <- LAA_FU_gene_map[ensembl_gene_id_clean %in% FU_SP_genes]
# FU & COJO
intersect_FU_COJO_gene <- LAA_FU_gene_map[ensembl_gene_id_clean %in% FU_GWAS_genes]
# SP & COJO
intersect_SP_COJO_gene <- LAA_SP_gene_map[ensembl_gene_id_clean %in% SP_GWAS_genes]

LAA_FU_v7_FDR[FDR_tissue <= 0.2, ensembl_gene_id]
LAA_FU_v8_FDR[FDR_tissue <= 0.2, ensembl_gene_id]
LAA_SP_v8_FDR[FDR_tissue <= 0.2, ensembl_gene_id]
LAA_SP_v8_FDR[FDR_tissue <= 0.2, ensembl_gene_id]



#union_gene <- Reduce(union, list(LAA_FU_v7_FDR[FDR_tissue <= 0.2, gene_name],
#                                 LAA_FU_v8_FDR[FDR_tissue <= 0.2, gene_name],
#                                 LAA_SP_v7_FDR[FDR_tissue <= 0.2, gene_name],
#                                 LAA_SP_v8_FDR[FDR_tissue <= 0.2, gene_name]))

union_gene <- Reduce(union, list(intersect_FU_SP_gene,
                                 intersect_FU_COJO_gene,
                                 intersect_SP_COJO_gene))

coloc_PP4_0.7_gene <- c('ITGAV', 'RTN4', 'DAPK2', 'TWIST1')

#head(unique(LAA_gene$gene_name))
#intersection_gene <- intersect_FU_SP_gene$gene_name

gene_map <- bitr(
  coloc_PP4_0.7_gene,   #intersection_gene, union_gene, coloc_PP4_0.7_gene
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
universe_map <- bitr(
  unique(LAA_gene$gene_name),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
nrow(gene_map)
nrow(universe_map)
#setdiff(union_gene, gene_map$SYMBOL)

ego_bp <- enrichGO(
  gene          = unique(gene_map$SYMBOL),
  universe      = unique(universe_map$SYMBOL),
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.25
)

ego_union_df <- as.data.frame(ego_bp)%>%
  filter(Count > 1) %>%
  arrange(p.adjust) 
nrow(ego_union_df)

N <- 15356
K <- 28
n <- 7
k <- 2
p_value <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
p_value


plot_df <- ego_union_df %>%
  arrange(desc(RichFactor)) %>%
  slice(1:20) %>%
  mutate(
    geneID = gsub("/", ", ", geneID)  
  )

p <- ggplot(plot_df,
            aes(x = RichFactor,
                y = reorder(Description, RichFactor))) +
  
  geom_point(aes(color = p.adjust),
             size = 4) +
  
  geom_text(aes(label = geneID),
            hjust = -0.1,
            size = 3) +
  
  scale_color_gradient(
    low = "#d73027",
    high = "#4575b4",
    name = "p.adjust"
  ) +
  
  labs(
    x = "Rich Factor",
    y = "GO Biological Process",
    title = "GO Biological Process Enrichment of LAA Colocalization PP.H4 > 0.7 Candidate Genes Top20"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  ) +
  
  coord_cartesian(xlim = c(0, max(plot_df$RichFactor)*1.35))

ggsave(
  "/mnt/data/stroke_pipeline_test/pathway/GOBP/LAA/LAA_coloc_PP4_GOBP_ORA_fdr005_20260320.png",
  plot = p,
  width = 15,
  height = 8,
  dpi = 600
)

write.csv(
  ego_union_df%>%arrange(desc(RichFactor)),
  "/mnt/data/stroke_pipeline_test/pathway/GOBP/LAA/LAA_coloc_PP4_GOBP_ORA_fdr005_count2_202603020.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_bp)%>%arrange(desc(RichFactor)),
  "/mnt/data/stroke_pipeline_test/pathway/GOBP/LAA/LAA_coloc_PP4_GOBP_ORA_fdr005_202603020.csv",
  row.names = FALSE
)



ekegg <- enrichKEGG(
  gene          = unique(gene_map$ENTREZID), #ENTREZID
  universe      = unique(universe_map$ENTREZID),
  organism      = "hsa",
  pvalueCutoff  = 0.1,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.25
)

ekegg <- setReadable(
  ekegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

ekegg_union_df <- as.data.frame(ekegg) %>%
  filter(Count > 1) %>%
  arrange(desc(RichFactor)) %>%
  mutate(
    geneID = gsub("/", ", ", geneID)   # gene顯示比較好讀
  )

plot_df <- as.data.frame(ekegg) %>%
  filter(Count > 1) %>%
  arrange(desc(RichFactor)) %>%
  mutate(
    geneID = gsub("/", ", ", geneID)
  )

p_kegg <- ggplot(plot_df,
                 aes(x = RichFactor,
                     y = reorder(Description, RichFactor))) +
  
  geom_point(aes(color = p.adjust),
             size = 4) +
  
  geom_text(aes(label = geneID),
            hjust = -0.1,
            size = 3) +
  
  scale_color_gradient(
    low = "#d73027",
    high = "#4575b4",
    name = "p.adjust"
  ) +
  
  labs(
    x = "Rich Factor",
    y = "KEGG Pathway",
    title = "KEGG Pathway Enrichment of LAA Colocalization PP.H4 > 0.7 Candidate Genes"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  ) +
  
  coord_cartesian(xlim = c(0, max(plot_df$RichFactor)*1.35))

p_kegg

ggsave(
  "/mnt/data/stroke_pipeline_test/pathway/KEGG/LAA/LAA_coloc_PP4_KEGG_ORA_fdr01_20260321.png",
  plot = p_kegg,
  width = 15,
  height = 8,
  dpi = 600
)

write.csv(
  ekegg_union_df%>%arrange(desc(RichFactor)),
  "/mnt/data/stroke_pipeline_test/pathway/KEGG/LAA/LAA_coloc_PP4_KEGG_ORA_fdr01_count2_20260321.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ekegg)%>%arrange(desc(RichFactor)),
  "/mnt/data/stroke_pipeline_test/pathway/KEGG/LAA/LAA_coloc_PP4_KEGG_ORA_fdr01_20260321.csv",
  row.names = FALSE
)


