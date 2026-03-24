#!/bin/bash

# =========================================
# Sample and Variant QC Pipeline (LAA, SVO, TWB1)
# =========================================

# -------------------------
# LAA dataset
# -------------------------
cd /home/sysadmin/Desktop/stroke_pipeline_test/LAA

# Check sex distribution
awk '{print $5}' LAA444.fam | sort | uniq -c

# Remove samples with missing sex (sex = 0)
awk '$5 == 0 {print $1, $2}' LAA444.fam > LAA444_sex0.remove

plink --bfile LAA444 \
  --remove LAA444_sex0.remove \
  --make-bed \
  --out LAA444_sexFiltered

# Check sample size after filtering
wc -l LAA444_sexFiltered.fam
awk '{print $5}' LAA444_sexFiltered.fam | sort | uniq -c

# Keep only biallelic SNPs
plink2 \
  --bfile LAA444_sexFiltered \
  --max-alleles 2 \
  --make-bed \
  --out LAA444_biallelic

# Apply QC filters (missingness, MAF, HWE)
plink2 \
  --bfile LAA444_biallelic \
  --mind 0.02 \
  --geno 0.02 \
  --maf 0.01 \
  --hwe 1e-6 \
  --make-bed \
  --out LAA444_hwe_QCFiltered

# Generate SNP full IDs
awk '{print $1"_"$2"_"$4"_"$5"_"$6}' LAA444_hwe_QCFiltered.bim | sort > LAA444_hwe_QCFiltered_fullid.txt


# -------------------------
# SVO dataset
# -------------------------
cd /home/sysadmin/Desktop/stroke_pipeline_test/SVO

# Check sex distribution
awk '{print $5}' SVO342.fam | sort | uniq -c

# Remove samples with missing sex
awk '$5 == 0 {print $1, $2}' SVO342.fam > SVO342_sex0.remove

plink --bfile SVO342 \
  --remove SVO342_sex0.remove \
  --make-bed \
  --out SVO342_sexFiltered

# Check sample size
wc -l SVO342_sexFiltered.fam
awk '{print $5}' SVO342_sexFiltered.fam | sort | uniq -c

# Keep biallelic SNPs
plink2 \
  --bfile SVO342_sexFiltered \
  --max-alleles 2 \
  --make-bed \
  --out SVO342_biallelic

# Apply QC filters
plink2 \
  --bfile SVO342_biallelic \
  --mind 0.02 \
  --geno 0.02 \
  --maf 0.01 \
  --hwe 1e-6 \
  --make-bed \
  --out SVO342_hwe_QCFiltered

# Generate SNP full IDs
awk '{print $1"_"$2"_"$4"_"$5"_"$6}' SVO342_hwe_QCFiltered.bim | sort > SVO342_hwe_QCFiltered_fullid.txt


# -------------------------
# TWB1 dataset
# -------------------------
cd /home/sysadmin/Desktop/stroke_pipeline_test/TWB1

# Check sex distribution
awk '{print $5}' TWB1.hg19.fam | sort | uniq -c

# Keep biallelic SNPs
plink2 \
  --bfile TWB1.hg19 \
  --max-alleles 2 \
  --make-bed \
  --out TWB1.hg19_biallelic

# Restrict to autosomes (chr 1–22)
plink2 \
  --bfile TWB1.hg19_biallelic \
  --chr 1-22 \
  --make-bed \
  --out TWB1.hg19_autosome

# Apply QC filters
plink2 \
  --bfile TWB1.hg19_autosome \
  --mind 0.02 \
  --geno 0.02 \
  --maf 0.01 \
  --hwe 1e-6 \
  --make-bed \
  --out TWB1.hg19_hwe_QCFiltered

# Generate SNP full IDs
awk '{print $1"_"$2"_"$4"_"$5"_"$6}' TWB1.hg19_hwe_QCFiltered.bim | sort > TWB1.hg19_hwe_QCFiltered_fullid.txt


echo "QC pipeline finished successfully"


# -------------------------
# Merge LAA and SVO
# -------------------------
cd /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO

# Extract intersecting SNP rsIDs based on full ID
awk 'NR==FNR {a[$1]=1; next} ($1"_"$2"_"$4"_"$5"_"$6) in a {print $2}' \
/home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA_SVO_intersect_fullid.txt \
/home/sysadmin/Desktop/stroke_pipeline_test/LAA/LAA444_hwe_QCFiltered.bim \
> /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA_SVO_intersect_rsid.txt

# Subset LAA dataset
plink2 \
  --bfile /home/sysadmin/Desktop/stroke_pipeline_test/LAA/LAA444_hwe_QCFiltered \
  --extract /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA_SVO_intersect_rsid.txt \
  --make-bed \
  --out /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA444.intersect

# Subset SVO dataset
plink2 \
  --bfile /home/sysadmin/Desktop/stroke_pipeline_test/SVO/SVO342_hwe_QCFiltered \
  --extract /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA_SVO_intersect_rsid.txt \
  --make-bed \
  --out /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/SVO342.intersect

# Merge LAA and SVO
plink \
  --bfile /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA444.intersect \
  --bmerge /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/SVO342.intersect \
  --make-bed \
  --out /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA_SVO_merged

# Apply QC filters after merging
plink2 \
  --bfile LAA_SVO_merged \
  --mind 0.02 \
  --geno 0.02 \
  --maf 0.01 \
  --make-bed \
  --out LAA_SVO_QCFiltered

# Generate SNP full IDs
awk '{print $1"_"$2"_"$4"_"$5"_"$6}' LAA_SVO_QCFiltered.bim | sort > LAA_SVO_QCFiltered_fullid.txt


# -------------------------
# Merge TWB1 with LAA_SVO
# -------------------------
cd /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO

# Extract intersecting SNP rsIDs
awk 'NR==FNR {a[$1]=1; next} ($1"_"$2"_"$4"_"$5"_"$6) in a {print $2}' \
/home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_intersect_fullid.txt \
/home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA_SVO_QCFiltered.bim \
> /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_intersect_rsid.txt

# Subset LAA_SVO dataset
plink2 \
  --bfile /home/sysadmin/Desktop/stroke_pipeline_test/LAA_SVO/LAA_SVO_QCFiltered \
  --extract /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_intersect_rsid.txt \
  --make-bed \
  --out /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/LAA_SVO.intersect

# Subset TWB1 dataset
plink2 \
  --bfile /home/sysadmin/Desktop/stroke_pipeline_test/TWB1/TWB1.hg19_hwe_QCFiltered \
  --extract /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_intersect_rsid.txt \
  --make-bed \
  --out /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1.hg19.intersect

# Merge TWB1 with LAA_SVO
plink \
  --bfile /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/LAA_SVO.intersect \
  --bmerge /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1.hg19.intersect \
  --make-bed \
  --out /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO/TWB1_LAA_SVO_merged


echo "Merging pipeline completed successfully"


# =========================================
# GWAS and Post-QC Pipeline (TWB1_LAA_SVO)
# =========================================

cd /home/sysadmin/Desktop/stroke_pipeline_test/TWB1_LAA_SVO

# Apply basic QC filters (sample missingness, SNP missingness, MAF)
plink2 \
  --bfile TWB1_LAA_SVO_merged \
  --mind 0.02 \
  --geno 0.02 \
  --maf 0.01 \
  --make-bed \
  --out TWB1_LAA_SVO_QCFiltered


# GWAS without covariate correction
plink \
  --bfile TWB1_LAA_SVO_QCFiltered \
  --logistic hide-covar \
  --out TWB1_LAA_SVO_GWAS_nocorrected


# LD pruning for heterozygosity estimation
plink2 \
  --bfile TWB1_LAA_SVO_QCFiltered \
  --indep-pairwise 50 5 0.2 \
  --out TWB1_LAA_SVO_prune


# Estimate heterozygosity (based on pruned SNPs)
plink2 \
  --bfile TWB1_LAA_SVO_QCFiltered \
  --extract TWB1_LAA_SVO_prune.prune.in \
  --het \
  --out TWB1_LAA_SVO_het


# Remove heterozygosity outliers
plink2 \
  --bfile TWB1_LAA_SVO_QCFiltered \
  --remove het_outliers.txt \
  --make-bed \
  --out TWB1_LAA_SVO_noHetOutlier


# Remove related individuals (KING)
plink2 \
  --bfile TWB1_LAA_SVO_noHetOutlier \
  --king-cutoff 0.354

plink2 \
  --bfile TWB1_LAA_SVO_noHetOutlier \
  --keep plink2.king.cutoff.in.id \
  --make-bed \
  --out TWB1_LAA_SVO_noDuplicates


# LD pruning for PCA
plink2 \
  --bfile TWB1_LAA_SVO_noDuplicates \
  --maf 0.01 \
  --indep-pairwise 200 50 0.2 \
  --out pca_prune


# Compute principal components (top 10 PCs)
plink2 \
  --bfile TWB1_LAA_SVO_noDuplicates \
  --extract pca_prune.prune.in \
  --pca approx 10 \
  --out TWB1_LAA_SVO_PCA


# Remove PCA outliers
plink2 \
  --bfile TWB1_LAA_SVO_noDuplicates \
  --keep no_outliers_sample_id.txt \
  --make-bed \
  --out TWB1_LAA_SVO_noPCAoutliers


# GWAS after PCA outlier removal
plink \
  --bfile TWB1_LAA_SVO_noPCAoutliers \
  --logistic hide-covar \
  --out TWB1_LAA_SVO_GWAS_noPCAoutliers


echo "GWAS pipeline completed successfully"


# =========================================
# Imputation Pipeline (VCF + Beagle)
# =========================================

cd /home/sysadmin/Desktop/stroke_pipeline_test/vcf

# Convert PLINK to VCF (bgzip format)
plink \
  --bfile TWB1_LAA_SVO_noPCAoutliers \
  --recode vcf bgz \
  --out TWB_LAA_SVO_merged_final


# Normalize and sort VCF
bcftools norm -m -both -d exact TWB_LAA_SVO_merged_final.vcf.gz \
  | bcftools sort -Oz -o TWB_LAA_SVO_merged_input.vcf.gz

bcftools index TWB_LAA_SVO_merged_input.vcf.gz


# Split VCF by chromosome (1–22)
for chr in {1..22}; do
  bcftools view -r $chr TWB_LAA_SVO_merged_input.vcf.gz -Oz -o TWB_chr${chr}.vcf.gz
  bcftools index -f TWB_chr${chr}.vcf.gz
done


# Conform genotype to reference (example: chr1)
java -jar /home/sysadmin/Desktop/tools/beagle/conform-gt.24May16.cee.jar \
  ref=/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/EAS_clean/chr1_EAS_clean.vcf.gz \
  gt=TWB_chr1.vcf.gz \
  chrom=1 \
  match=POS \
  out=TWB1_chr1_conformed

# Check log
less TWB1_chr1_conformed.log


# Run Beagle imputation (chr 1–22)
for chr in {1..22}; do
  java -Xmx110g -jar /home/sysadmin/Desktop/tools/beagle/beagle.27Feb25.75f.jar \
    gt=TWB1_chr${chr}_conformed.vcf.gz \
    ref=/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/EAS_clean/chr${chr}_EAS_clean.vcf.gz \
    map=/home/sysadmin/Desktop/tools/beagle/plink.GRCh37.map/plink.chr${chr}.GRCh37.map \
    out=imputed_chr${chr} \
    nthreads=6 \
    window=30 \
    overlap=2
done


# Index imputed VCFs
for chr in {1..22}; do
  tabix -p vcf imputed_chr${chr}.vcf.gz
done


# Post-imputation QC (example: chr1)
cd /mnt/data/stroke_pipeline_test/imputed_vcf/chr1

bcftools view \
  -i 'INFO/DR2>=0.7' \
  -Oz \
  -o imputed_chr1_DR2_0.7.vcf.gz \
  imputed_chr1.vcf.gz


# Convert VCF back to PLINK format (example: chr1)
cd /mnt/data/stroke_pipeline_test/imputed_vcf/chr1/DR2_0.7

plink2 \
  --vcf imputed_chr1_DR2_0.7.vcf.gz \
  --make-bed \
  --out stroke_chr1_DR2_0.7 \
  --max-alleles 2


echo "Imputation pipeline completed successfully"


#!/bin/bash

# =========================================
# Post-imputation QC and GWAS (chr1)
# =========================================

cd /mnt/data/stroke_pipeline_test/imputed_vcf/chr1/DR2_0.7


# Step 1: Apply QC filters
plink \
  --bfile stroke_chr1_DR2_0.7 \
  --mind 0.02 --geno 0.02 --maf 0.01 --hwe 1e-6 \
  --make-bed \
  --out stroke_chr1_DR2_0.7_QC


# Step 2: Update missing rsIDs using dbSNP reference
awk 'NR==FNR{db[$1":"$3]=$4; next}
     {
       key=$1":"$4
       if($2=="." && key in db) $2=db[key]
       print
     }' \
/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/ncbi_dbsnp/chr1_snp_list.bed \
stroke_chr1_DR2_0.7_QC.bim \
> stroke_chr1_DR2_0.7_rsID.bim


# Step 3: Check missing rsID rate
awk '{total++} $2=="." {missing++} END{
  print "Total SNPs:", total
  print "Missing rsIDs:", missing
  print "Percentage:", missing/total*100 "%"
}' stroke_chr1_DR2_0.7_rsID.bim


# Step 4: Output SNPs without rsID
awk '$2=="." {print $2}' stroke_chr1_DR2_0.7_rsID.bim | sort | uniq > stroke_chr1_DR2_0.7_exclude_missingID.txt


# Step 5: Replace .bim file with updated rsIDs
cp stroke_chr1_DR2_0.7_rsID.bim stroke_chr1_DR2_0.7_QC.bim


# Step 6: Remove SNPs without rsID
plink \
  --bfile stroke_chr1_DR2_0.7_QC \
  --exclude stroke_chr1_DR2_0.7_exclude_missingID.txt \
  --make-bed \
  --out stroke_chr1_DR2_0.7_fixed


# Step 7: Restore original .fam file (after Beagle)
merged="/mnt/data/stroke_pipeline_test/TWB1_LAA_SVO_noPCAoutliers.fam"
fixed="stroke_chr1_DR2_0.7_fixed.fam"

awk '{print $1"_"$2}' $merged > merged_ids.txt
awk '{print $2}' $fixed > fixed_ids.txt

if diff -q merged_ids.txt fixed_ids.txt >/dev/null; then
    echo "IDs match. Replacing .fam file..."
    cp $merged $fixed
else
    echo "IDs do not match. No replacement performed."
fi

rm merged_ids.txt fixed_ids.txt


# Step 8: Subset LAA samples and run GWAS
plink \
  --bfile stroke_chr1_DR2_0.7_fixed \
  --keep LAA_sample_id.txt \
  --make-bed \
  --out TWB1_LAA_chr1_DR2_0.7_QCFiltered

plink2 \
  --bfile TWB1_LAA_chr1_DR2_0.7_QCFiltered \
  --glm allow-no-covars \
  --threads 8 \
  --out TWB1_LAA_chr1_DR2_0.7


# Step 9: Subset SVO samples and run GWAS
plink \
  --bfile stroke_chr1_DR2_0.7_fixed \
  --keep SVO_sample_id.txt \
  --make-bed \
  --out TWB1_SVO_chr1_DR2_0.7_QCFiltered

plink2 \
  --bfile TWB1_SVO_chr1_DR2_0.7_QCFiltered \
  --glm allow-no-covars \
  --threads 8 \
  --out TWB1_SVO_chr1_DR2_0.7


echo "Post-imputation QC and GWAS completed successfully"


# =========================================
# Fine-mapping and TWAS Pipeline
# =========================================

# -------------------------
# COJO fine-mapping (chr 1–22)
# -------------------------
for chr in {1..22}; do
  echo "Running COJO for chr${chr}"

  gcta64 \
    --bfile /mnt/data/stroke_pipeline_test/imputed_vcf/chr${chr}/DR2_0.7/TWB1_LAA_chr${chr}_DR2_0.7_QCFiltered \
    --maf 0.01 \
    --cojo-file /mnt/data/stroke_pipeline_test/imputed_vcf/chr${chr}/DR2_0.7/TWB1_LAA_chr${chr}_DR2_0.7.cojo.ma \
    --cojo-slct \
    --cojo-p 1e-5 \
    --cojo-wind 10000 \
    --diff-freq 0.2 \
    --out /mnt/data/stroke_pipeline_test/cojo_result/TWB1_LAA_chr${chr}_DR2_0.7_cojo_p1e5
done


# -------------------------
# FUSION TWAS (GTEx v8, chr 1–22)
# -------------------------
GWAS_DIR=/home/sysadmin/Desktop/dyc_lab/FUSION/GWAS/TWBSTROKE/hg19/LAA_DR2_0.7
WEIGHT_DIR=/home/sysadmin/Desktop/dyc_lab/FUSION/WEIGHTS
LDREF=/home/sysadmin/Desktop/dyc_lab/FUSION/LDREF/LDREF/1000G.EUR.
OUT_DIR=/home/sysadmin/Desktop/dyc_lab/FUSION/RESULTS/TWBSTROKE/hg19/LAA_DR2_0.7
TISSUE_LIST=${WEIGHT_DIR}/tissue_list.txt
FUSION_SCRIPT=/home/sysadmin/Desktop/dyc_lab/FUSION/fusion_twas/FUSION.assoc_test.R

for CHR in {1..22}; do

  SUMSTATS=${GWAS_DIR}/TWB1_LAA_chr${CHR}_fusion_ready.sumstats
  CHR_OUTDIR=${OUT_DIR}/chr${CHR}

  echo "Running TWAS (GTEx v8) for chr${CHR}"

  while read TISSUE; do

    POS=${WEIGHT_DIR}/GTExv8.EUR.${TISSUE}.pos

    echo "  → Tissue: ${TISSUE}"

    Rscript ${FUSION_SCRIPT} \
      --sumstats ${SUMSTATS} \
      --weights ${POS} \
      --weights_dir ${WEIGHT_DIR} \
      --ref_ld_chr ${LDREF} \
      --chr ${CHR} \
      --out ${CHR_OUTDIR}/TWAS_chr${CHR}_DR207_LAA_GTExv8_${TISSUE}.dat

  done < ${TISSUE_LIST}

done


# -------------------------
# FUSION TWAS (GTEx v7, chr 1–22)
# -------------------------
GWAS_DIR=/home/sysadmin/Desktop/dyc_lab/FUSION/GWAS/TWBSTROKE/hg19/LAA_DR2_0.7
WEIGHT_DIR=/home/sysadmin/Desktop/dyc_lab/FUSION/WEIGHTS_v7/GTEx.ALL
LDREF=/home/sysadmin/Desktop/dyc_lab/FUSION/LDREF/LDREF/1000G.EUR.
OUT_DIR=/home/sysadmin/Desktop/dyc_lab/FUSION/RESULTS_v7/TWBSTROKE/hg19/LAA_DR2_0.7
TISSUE_LIST=${WEIGHT_DIR}/tissue_list_v7.txt
FUSION_SCRIPT=/home/sysadmin/Desktop/dyc_lab/FUSION/fusion_twas/FUSION.assoc_test.R

for CHR in {1..22}; do

  SUMSTATS=${GWAS_DIR}/TWB1_LAA_chr${CHR}_fusion_ready.sumstats
  CHR_OUTDIR=${OUT_DIR}/chr${CHR}

  echo "Running TWAS (GTEx v7) for chr${CHR}"

  while read TISSUE; do

    POS=${WEIGHT_DIR}/GTEx.${TISSUE}.pos

    echo "  → Tissue: ${TISSUE}"

    Rscript ${FUSION_SCRIPT} \
      --sumstats ${SUMSTATS} \
      --weights ${POS} \
      --weights_dir ${WEIGHT_DIR} \
      --ref_ld_chr ${LDREF} \
      --chr ${CHR} \
      --out ${CHR_OUTDIR}/TWAS_chr${CHR}_DR207_LAA_GTExv7_${TISSUE}.dat

  done < ${TISSUE_LIST}

done


# -------------------------
# SPrediXcan (GTEx v7)
# -------------------------
conda activate metaxcan

TISSUE_FILE="/home/sysadmin/Desktop/dyc_lab/S_Predixcan/tissues.txt"
MODEL_DIR="/home/sysadmin/Desktop/dyc_lab/S_Predixcan/model"
GWAS_DIR="/home/sysadmin/Desktop/dyc_lab/S_Predixcan/gwas/SVO_DR2_0.7"
OUTPUT_DIR="/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result/SVO"

SPREDIXCAN="/home/sysadmin/Desktop/tools/MetaXcan/software/SPrediXcan.py"

while read tissue; do
    echo "Running SPrediXcan v7 for $tissue"

    python $SPREDIXCAN \
        --model_db_path "${MODEL_DIR}/gtex_v7_${tissue}_imputed_europeans_tw_0.5_signif.db" \
        --covariance    "${MODEL_DIR}/gtex_v7_${tissue}_imputed_eur_covariances.txt.gz" \
        --gwas_folder   "$GWAS_DIR" \
        --gwas_file_pattern "TWB1_SVO_DR2_07_chr1to22_spredixcan_ready.sumstats" \
        --snp_column SNP \
        --effect_allele_column A1 \
        --non_effect_allele_column A2 \
        --zscore_column Z \
        --output_file "${OUTPUT_DIR}/${tissue}_TWB1_SVO_DR2_0.7_SPrediXcan.csv"

done < "$TISSUE_FILE"


# -------------------------
# SPrediXcan (GTEx v8 elastic net)
# -------------------------
TISSUE_FILE="/home/sysadmin/Desktop/dyc_lab/S_Predixcan/tissues_v8_elastic_net.txt"
MODEL_DIR="/home/sysadmin/Desktop/dyc_lab/S_Predixcan/model/GTEx_v8/elastic_net_models"
GWAS_FILE="/home/sysadmin/Desktop/dyc_lab/S_Predixcan/gwas/SVO_DR2_0.7/TWB1_SVO_DR2_07_chr1to22_spredixcan_ready.sumstats"
OUTPUT_DIR="/home/sysadmin/Desktop/dyc_lab/S_Predixcan/result_v8/SVO"

SPREDIXCAN="/home/sysadmin/Desktop/tools/MetaXcan/software/SPrediXcan.py"

while read tissue; do
    echo "Running SPrediXcan v8 (elastic net) for $tissue"

    python $SPREDIXCAN \
        --model_db_path "${MODEL_DIR}/en_${tissue}.db" \
        --covariance    "${MODEL_DIR}/en_${tissue}.txt.gz" \
        --gwas_file    "$GWAS_FILE" \
        --snp_column SNP \
        --effect_allele_column A1 \
        --non_effect_allele_column A2 \
        --zscore_column Z \
        --gwas_N 22790 \
        --output_file "${OUTPUT_DIR}/${tissue}_TWB1_SVO_DR2_0.7_SPrediXcan_v8_en.csv"

done < "$TISSUE_FILE"


echo "Fine-mapping and TWAS pipeline completed successfully"
