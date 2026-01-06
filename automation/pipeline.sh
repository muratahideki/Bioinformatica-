#!/bin/bash
set -euo pipefail

# ===============================
# CONFIGURAÇÕES
# ===============================
ENV_NAME="ngs_env"

THREADS=10

RAW_DIR_IPT="data/raw/IPT"
RAW_DIR_NIPT="data/raw/NIPT"

RESULTS_IPT_SPLIT_FILE="data/raw/IPT/RESULTS"
RESULTS_NIPT_SPLIT_FILE="data/raw/NIPT/RESULTS"

TRIM_DIR_IPT="data/trimmed/IPT"
TRIM_DIR_NIPT="data/trimmed/NIPT"

QC_DIR="results/fastqc"
ALIGN_DIR="results/aligned"
COUNT_DIR="results/counts"

GENOME_INDEX="reference/hisat2/genome_index"
GTF="reference/annotation.gtf"

# ===============================
# ATIVAR CONDA
# ===============================
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# ===============================
# CRIAR DIRETÓRIOS
# ===============================
mkdir -p \
  $RAW_DIR_IPT $RAW_DIR_NIPT \
  $RESULTS_IPT_SPLIT_FILE $RESULTS_NIPT_SPLIT_FILE \
  $TRIM_DIR_IPT $TRIM_DIR_NIPT\
  $QC_DIR $ALIGN_DIR $COUNT_DIR

echo "Diretórios criados"

# ===============================
# Baixar os arquivos SRA e separar o split file 
# ===============================

pushd "$RAW_DIR_IPT"

IPT_SAMPLES= (
  SRR13447971	
  SRR13447972
  SRR13447975	
  SRR13447976
  SRR13447981	
  SRR13447982
)

for SAMPLE in "${IPT_SAMPLES[@]}"; do
  prefetch "$SAMPLE"
  fasterq-dump "$SAMPLE" \
  --split-files \
  --gzip \
  --threads 8 \
  -O "$RESULTS_IPT_SPLIT_FILE"
done 
popd

pushd "$RAW_DIR_NIPT"
NIPT_SAMPLES=(
  SRR13447973
  SRR13447974
  SRR13447977
  SRR13447978
  SRR13447979
  SRR13447980
)

for SAMPLE in "${NIPT_SAMPLES[@]}"; do
  prefetch "$SAMPLE"
  fasterq-dump "$SAMPLE" \
    --split-files \
    --gzip \
    --threads 8 \
    -O "$RESULTS_NIPT_SPLIT_FILE"
done
popd


# ===============================
# 2️⃣ TRIMMOMATIC
# ===============================
echo "✂️ Rodando Trimmomatic"

# ===============================
# TRIMMOMATIC – IPT
# ===============================


ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa"


for R1 in "$RESULTS_IPT_SPLIT_FILE"/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="$RESULTS_IPT_SPLIT_FILE/${SAMPLE}_2.fastq.gz"

    trimmomatic PE \
        -threads $THREADS \
        "$R1" "$R2" \
        "$TRIM_DIR_IPT/${SAMPLE}_1_paired.fastq.gz" "$TRIM_DIR_IPT/${SAMPLE}_1_unpaired.fastq.gz" \
        "$TRIM_DIR_IPT/${SAMPLE}_2_paired.fastq.gz" "$TRIM_DIR_IPT/${SAMPLE}_2_unpaired.fastq.gz" \
        ILLUMINACLIP:$ADAPTERS:2:30:10 \
        HEADCROP:12 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:36  
done

# ===============================
# TRIMMOMATIC – NIPT
# ===============================
echo "✂️ Trimmomatic (NIPT)"

for R1 in "$RESULTS_NIPT_SPLIT_FILE"/*_1.fastq.gz; do
  SAMPLE=$(basename "$R1" _1.fastq.gz)
  R2="$RESULTS_NIPT_SPLIT_FILE"/${SAMPLE}_2.fastq.gz

    trimmomatic PE \
      -threads $THREADS \
      "$R1" "$R2" \
      "$TRIM_DIR_NIPT/${SAMPLE}_1_paired.fastq.gz" "$TRIM_DIR_NIPT/${SAMPLE}_1_unpaired.fastq.gz" \
      "$TRIM_DIR_NIPT/${SAMPLE}_2_paired.fastq.gz" "$TRIM_DIR_NIPT/${SAMPLE}_2_unpaired.fastq.gz" \
      ILLUMINACLIP:$ADAPTERS:2:30:10 \
      HEADCROP:12 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:36  
done

echo "Removendo arquivos unpaired..."

rm -f "$TRIM_DIR_IPT"/*_unpaired.fastq.gz
rm -f "$TRIM_DIR_NIPT"/*_unpaired.fastq.gz

