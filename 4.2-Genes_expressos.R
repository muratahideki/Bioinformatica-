# ==============================
# Script DESeq2 com featureCounts
# ==============================

# 1️⃣ Carregar pacotes
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
library(DESeq2)

# 2️⃣ Ler arquivo do featureCounts
# Substitua "counts.txt" pelo caminho do seu arquivo
counts <- read.table("counts.txt", 
                     header = TRUE, 
                     row.names = 1, 
                     comment.char = "#")

# Verificar os dados
head(counts)

# 3️⃣ Criar tabela de condições das amostras
# Substitua pelos nomes reais das condições do seu experimento
# A ordem deve corresponder às colunas do counts
coldata <- data.frame(
  row.names = colnames(counts),
  condition = c("control", "control", "treatment", "treatment")
)

# Verificar coldata
print(coldata)

# 4️⃣ Criar objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# 5️⃣ Filtrar genes com contagem muito baixa (opcional, mas recomendado)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 6️⃣ Rodar análise DESeq2
dds <- DESeq(dds)

# 7️⃣ Extrair resultados
res <- results(dds)

# Verificar resultados
head(res)

# 8️⃣ Salvar resultados em arquivo CSV
write.csv(as.data.frame(res), file = "DESeq2_results.csv")

# ==============================
# Fim do script
# ==============================
