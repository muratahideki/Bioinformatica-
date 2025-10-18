library(DESeq2)

# Ler arquivo do featureCounts
counts_raw <- read.table("counts.txt",
                         header = TRUE,
                         comment.char = "#",
                         stringsAsFactors = FALSE)

# Selecionar apenas colunas necessárias
counts <- counts_raw[, c(1, 7:ncol(counts_raw))]

# Tornar Geneid o nome das linhas
rownames(counts) <- counts$Geneid
counts <- counts[, -1]

# Criar o data frame de condições (substitua conforme seu caso)
coldata <- data.frame(
  row.names = colnames(counts),
  condition = c("control")  # altere conforme o número de amostras
)

#Cria o objeto DDS
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Filtrar genes com baixa contagem
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Rodar o DESeq2
dds <- DESeq(dds)

# Extrair resultados
res <- results(dds)

# Salvar resultados
write.csv(as.data.frame(res), file = "DESeq2_results.csv")

