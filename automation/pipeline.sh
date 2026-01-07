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
  $TRIM_DIR_IPT $TRIM_DIR_NIPT \
  $ALIGN_DIR $COUNT_DIR \
  $GENOME_INDEX

mkdir -p hisat2

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

# ===============================
# Mapeamento-  Hisat2
# ===============================

GENOME_INDEX="automation/hisat2/genome_index"
GTF="automation/annotation.gtf"

# download do genoma
curl --cookie jgi_session=/api/sessions/72c8fbec4d695cdee2d873b55b50b2e2 --output download.20260106.211503.zip -d "{\"ids\":{\"Phytozome-771\":{\"file_ids\":[\"642da257ed041a78f31f89b2\",\"642da256ed041a78f31f899c\"],\"top_hit\":\"642da256ed041a78f31f89a4\"}},\"api_version\":\"2\"}" -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/
unzip download.20260106.211503.zip

arquivos="Phytozome/PhytozomeV13/SofficinarumxspontaneumR570/v2.1"

pasta_genoma="$arquivos"/assembly
genoma="$pasta_genoma"/SofficinarumxspontaneumR570_771_v2.0.softmasked.fa.gz
gunzip $genoma
genoma="${genoma%.gz}"

#download anotation
curl --cookie jgi_session=/api/sessions/72c8fbec4d695cdee2d873b55b50b2e2 --output download.20260106.220955.zip -d "{\"ids\":{\"Phytozome-771\":{\"file_ids\":[\"642da256ed041a78f31f899a\"],\"top_hit\":\"642da256ed041a78f31f89a4\"}},\"api_version\":\"2\"}" -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/
unzip download.20260106.220955.zip
pasta_annotation="$arquivos"/annotation
annotation="$pasta_annotation"/SofficinarumxspontaneumR570_771_v2.1.gene.gff3.gz
gunzip $annotation
annotation="${annotation%.gz}"

# build 
hisat2-build \
-p $THREADS \
"$genoma" \
$GENOME_INDEX

# ===============================
# aLINHAMENTO ipt
# ===============================
echo "🧬 Alinhando IPT com HISAT2"

mkdir -p "$ALIGN_DIR/IPT"

for R1 in "$TRIM_DIR_IPT"/*_1_paired.fastq.gz; do
  SAMPLE=$(basename "$R1" _1_paired.fastq.gz)
  R2="$TRIM_DIR_IPT/${SAMPLE}_2_paired.fastq.gz"

  hisat2 \
    -p $THREADS \
    -x $GENOME_INDEX \
    -1 "$R1" \
    -2 "$R2" \
  | samtools view -@ $THREADS -bS - > "$ALIGN_DIR/IPT/${SAMPLE}.bam"

done


echo "🧬 Alinhando NIPT com HISAT2"

mkdir -p "$ALIGN_DIR/NIPT"

for R1 in "$TRIM_DIR_NIPT"/*_1_paired.fastq.gz; do 
  SAMPLE=$(basename "$R1" _1_paired.fastq.gz)
  R2="$TRIM_DIR_NIPT/${SAMPLE}_2_paired.fastq.gz"

  hisat2 \
    -p $THREADS \
    -x $GENOME_INDEX \
    -1 "$R1" \
    -2 "$R2" \
  | samtools view -@ $THREADS -bS - > "$ALIGN_DIR/NIPT/${SAMPLE}.bam"

done
