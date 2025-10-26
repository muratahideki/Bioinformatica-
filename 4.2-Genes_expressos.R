#--------------------------------------
# Script DESeq2 completo com visualizações
#--------------------------------------

# 1️⃣ Instalar e carregar pacotes
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

if (!requireNamespace("pheatmap", quietly = TRUE))
    install.packages("pheatmap")

library(DESeq2)
library(pheatmap)

# Ler arquivo counts do featureCounts
counts <- read.table("counts.txt",
                     header=TRUE,
                     sep="\t",
                     comment.char="#",
                     row.names=1)

# Seleciona apenas as colunas de contagem (IPT e NIPT)
counts_matrix <- counts[, c("IPT", "NIPT")]

# Verifica
head(counts_matrix)

# Criar colData para as 2 amostras
colData <- data.frame(
  condition = factor(c("IPT","NIPT"))
)
rownames(colData) <- colnames(counts_matrix)

#  Criar objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = colData,
  design = ~1  # sem replicatas, apenas normalização
)

#  Normalizar os dados
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# Salvar contagens normalizadas
write.csv(normalized_counts, "normalized_counts.csv", row.names=TRUE)

# Calcular log2 fold-change manualmente
log2fc <- log2(normalized_counts[, "IPT"] / normalized_counts[, "NIPT"])
log2fc_df <- data.frame(Gene=rownames(normalized_counts), log2FoldChange=log2fc)
write.csv(log2fc_df, "log2FoldChange.csv", row.names=FALSE)

# Heatmap dos top 20 genes
top_genes <- head(order(abs(log2fc), decreasing=TRUE), 20)

# PNG
png("heatmap_top20_genes.png", width=1200, height=900, res=150)
pheatmap(normalized_counts[top_genes, ],
         scale="row",
         annotation_col=colData,
         main="Top 20 genes (fold-change)")
dev.off()

# PDF
pdf("heatmap_top20_genes.pdf", width=8, height=6)
pheatmap(normalized_counts[top_genes, ],
         scale="row",
         annotation_col=colData,
         main="Top 20 genes (fold-change)")
dev.off()

# MA-plot
avg_counts <- rowMeans(normalized_counts)
png("MAplot.png", width=1000, height=800, res=150)
plot(log10(avg_counts + 1), log2fc,
     pch=20, cex=0.5, col="blue",
     xlab="log10(Mean normalized counts)",
     ylab="log2 Fold Change",
     main="MA-plot")
abline(h=0, col="red")
dev.off()

pdf("MAplot.pdf", width=8, height=6)
plot(log10(avg_counts + 1), log2fc,
     pch=20, cex=0.5, col="blue",
     xlab="log10(Mean normalized counts)",
     ylab="log2 Fold Change",
     main="MA-plot")
abline(h=0, col="red")
dev.off()

#  PCA
# rlog transform
rld <- rlog(dds, blind=TRUE)

pca_data <- prcomp(t(assay(rld)))

# Cores
colors <- c("red","blue")

# PCA PNG
png("PCA.png", width=1000, height=800, res=150)
plot(pca_data$x[,1], pca_data$x[,2],
     col=colors, pch=19,
     xlab="PC1", ylab="PC2",
     main="PCA - Top 2 PCs")
text(pca_data$x[,1], pca_data$x[,2],
     labels=colnames(rld), pos=3)
dev.off()

# PCA PDF
pdf("PCA.pdf", width=8, height=6)
plot(pca_data$x[,1], pca_data$x[,2],
     col=colors, pch=19,
     xlab="PC1", ylab="PC2",
     main="PCA - Top 2 PCs")
text(pca_data$x[,1], pca_data$x[,2],
     labels=colnames(rld), pos=3)
dev.off()
