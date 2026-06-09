#! /usr/bin/env Rscript

#####################################################################
# DE annalysis using limma-voom for male and female samples

#References:
# https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html
# http://bioconductor.unipi.it/packages/3.19/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# https://rtguides.it.tufts.edu/bio/tutorials/de-seq2-in-r-ood.html


#####################################################################

# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("DEGreport")
install.packages("pheatmap")
install.packages("tidyverse")

# Call packages
library("tidyverse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("apeglm") #Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes
library("DEGreport")

# Read in feature counts matrix
cts <- as.matrix(read.csv("/home/cvillazon/mRNA_SNAKEPipes_out/featureCounts/counts_female.tsv",
                    sep="\t", row.names=1, header=TRUE))

# Read in metadata
coldata <- read.csv("/home/cvillazon/DEseq2_females/coldata.tsv", row.names=1, sep="\t")
coldata$condition <- factor(coldata$condition, levels = c("ZT0Tx", "ZT0V", "ZT12Tx", "ZT12V"))
#coldata$condition <- factor(coldata$condition)

# DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ 0 + condition)

# Pre-filtering: remove rows with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]
#Note that more strict filtering to increase power is automatically 
#applied via independent filtering on the mean of normalized counts 
#within the results function


dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "ZT0V", "ZT12Tx"))
write.table(results, file = "ZT0VvsZT12Tx_females.csv", quote = FALSE, sep = ",",
                 dec = ".", row.names = TRUE,
                 col.names = TRUE)

results <- results(dds, contrast = c("condition", "ZT0V", "ZT0Tx"))
write.table(results, file = "ZT0VvsZT0Tx_females.csv", quote = FALSE, sep = ",",
                 dec = ".", row.names = TRUE,
                 col.names = TRUE)

results <- results(dds, contrast = c("condition", "ZT0V", "ZT12V"))
write.table(results, file = "ZT0VvsZT12V_females.csv", quote = FALSE, sep = ",",
                 dec = ".", row.names = TRUE,
                 col.names = TRUE)

results <- results(dds, contrast = c("condition", "ZT0Tx", "ZT12Tx"))
write.table(results, file = "ZT0TxvsZT12Tx_females.csv", quote = FALSE, sep = ",",
                 dec = ".", row.names = TRUE,
                 col.names = TRUE)

results <- results(dds, contrast = c("condition", "ZT12V", "ZT12Tx"))
write.table(results, file = "ZT12VvsZT12Tx_females.csv", quote = FALSE, sep = ",",
                 dec = ".", row.names = TRUE,
                 col.names = TRUE)

results <- results(dds, contrast = c("condition", "ZT12V", "ZT0Tx"))
write.table(results, file = "ZT12VvsZT0Tx_females.csv", quote = FALSE, sep = ",",
                 dec = ".", row.names = TRUE,
                 col.names = TRUE)

dds <- estimateSizeFactors(dds)

rld <- rlog(dds, blind = TRUE)

write.table(assay(rld),
            file = "rlog_norm_counts_females.csv",
            quote = FALSE,
            sep = ",",
            dec = ".",
            row.names = TRUE,
            col.names = TRUE)

# PCA plot
vsd <- vst(dds, blind = TRUE)  # Variance-stabilizing transformation
svg("PCA_plot.svg")
plotPCA(vsd, intgroup = "condition")
dev.off()

# Heatmap of Sample Distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
svg("Sample_Distance_Heatmap.svg")
pheatmap(sampleDistMatrix, main = "Sample Distance Heatmap")
dev.off()


# Add significance column
results$significant <- results$padj < 0.05 & abs(results$log2FoldChange) > 1

# Volcano plot
svg("Volcano_plot.svg")
ggplot(results, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 P-value")
dev.off()