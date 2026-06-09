#! /usr/bin/env Rscript

#####################################################################
# DE annalysis using limma-voom for male and female samples

#References:
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
# https://f1000research.com/articles/5-1408/v3

#####################################################################


# install packages 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")

BiocManager::install("edgeR")

# call pakages
library(limma)
library(edgeR)
library(tidyverse)
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for better labelling of points in the plot
library(gplots) # for heatmap


# read in feature counts matrix
# original file: /home/cvillazon/mRNA_SNAKEPipes_out/featureCounts/counts.tsv
# edited file: /home/cvillazon/mRNA_SNAKEPipes_out/featureCounts/counts_edit.tsv
# females: /home/cvillazon/mRNA_SNAKEPipes_out/featureCounts/counts_female.tsv
# males: /home/cvillazon/mRNA_SNAKEPipes_out/featureCounts/counts_male.tsv
counts <- read.delim("/home/cvillazon/results/mRNA_SNAKEPipes_out/featureCounts/counts_male.tsv", row.names = 1)

d0 <- DGEList(counts) # create DGEList object
d0 <- calcNormFactors(d0) # calculate normalization factors

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] # filter lowly expressed genes

# derive experiment information from the sample names (TZ0, TZ12, Tx, V)
#group <- c("ZT0Tx", "ZT0Tx", "ZT0Tx", "ZT0V", "ZT0V", "ZT0V", "ZT0Tx", "ZT0Tx",
#            "ZT0V", "ZT0V", "ZT0V", "ZT12Tx", "ZT12Tx", "ZT12Tx", "ZT12V",
#            "ZT12V", "ZT12V", "ZT12Tx", "ZT12Tx", "ZT12Tx", "ZT12V", "ZT12V")

#males (ZT0Tx_Rep2?)
group <- factor(c("ZT0Tx", "ZT0Tx", "ZT0Tx", "ZT0V", "ZT0V", "ZT0V", "ZT12Tx", "ZT12Tx",
            "ZT12Tx", "ZT12V", "ZT12V", "ZT12V"))

# Multidimensional scaling (MDS) plot

svg('mds_group.svg')
plotMDS(d, col = as.numeric(group))
dev.off()

#####################################################
# DEA

# specify the model to be fitted
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T) # voom transformation

svg('voom.svg')
voom(d, mm, plot = T)
dev.off() # the voom-plot provides a visual check on the level of filtering performed upstream

# What is voom doing?
# 1. Counts are transformed to log2 counts per million reads (CPM), where “per million reads” 
#   is defined based on the normalization factors we calculated earlier
# 2. A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
# 3. A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression
# 4. The smoothed curve is used to obtain weights for each gene and sample that are passed 
#   into limma along with the log2 CPMs.

# Typically, the “voom-plot” shows a decreasing trend between the means and variances resulting 
# from a combination of technical variation in the sequencing experiment and biological variation 
# amongst the replicate samples from different cell populations. Experiments with high biological 
# variation usually result in flatter trends, where variance values plateau at high expression values. 
# Experiments with low biological variation tend to result in sharp decreasing trends.


# make contrasts
# Comparison between groups
contr.matrix <- makeContrasts(
   ZT0TxvsZT0V = groupZT0Tx - groupZT0V,
   ZT12TxvsZT12V = groupZT12Tx - groupZT12V, 
   ZT0TxvsZT12Tx = groupZT0Tx - groupZT12Tx,
   ZT0VvsZT12V = groupZT0V - groupZT12V,
   ZT0TxvsZT12V = groupZT0Tx - groupZT12V,
   ZT12TxvsZT12V = groupZT12Tx - groupZT12V,
   levels = colnames(mm))

# fitting linear models for comparisons of interest
fit <- lmFit(y, mm)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit)
svg('efit.svg')
plotSA(efit)
dev.off()

# Examining the number of DE genes
summary(decideTests(efit))

# For a stricter definition on significance, one may require log-fold-changes 
#(log-FCs) to be above a minimum value.
tfit <- treat(fit, lfc=0)
dt <- decideTests(tfit)
summary(dt)

# de.common contains the number of DE genes that are common between comparisons

# 1: ZT0TxvsZT12Tx-ZT0VvsZT12V
de.common <- which(dt[,3]!=0 & dt[,4]!=0)

svg('venndiagram1_ZT0TxvsZT12Tx-ZT0VvsZT12V.svg')
vennDiagram(dt[,3:4], circle.col=c("turquoise", "salmon"))
dev.off()
write.fit(tfit, dt, file="results1_ZT0TxvsZT12Tx-ZT0VvsZT12V.txt")

# 2: ZT0TxvsZT12Tx-ZT0TxvsZT12V
de.common <- which(dt[,3]!=0 & dt[,5]!=0)

svg('venndiagram2_ZT0TxvsZT12Tx-ZT0TxvsZT12V.svg')
vennDiagram(dt[,c(3,5)], circle.col=c("turquoise", "salmon"))
dev.off()
write.fit(tfit, dt, file="results2_ZT0TxvsZT12Tx-ZT0TxvsZT12V.txt")

# 4: ZT0VvsZT12V-ZT0TxvsZT12V
de.common <- which(dt[,4]!=0 & dt[,5]!=0)

svg('venndiagram4_ZT0VvsZT12V-ZT0TxvsZT12V.svg')
vennDiagram(dt[,c(4,5)], circle.col=c("turquoise", "salmon"))
dev.off()
write.fit(tfit, dt, file="results4_ZT0VvsZT12V-ZT0TxvsZT12V.txt")


de.common <- which((dt[,4]!=0 & dt[,5]!=0) & (dt[,3]!=0))

svg('venndiagram3.svg')
vennDiagram(dt[,c(3,4,5)], circle.col=c("turquoise", "salmon", "lightgreen"))
dev.off()
write.fit(tfit, dt, file="results3.txt")

# top DE genes
ZT0TxvsZT12Tx <- topTreat(tfit, coef=3, n=Inf)
ZT0VvsZT12V <- topTreat(tfit, coef=4, n=Inf)
ZT0TxvsZT12V <- topTreat(tfit, coef=5, n=Inf)




# Volcano plot
ZT0TxvsZT12Tx$diffexpressed <- "Not significant"
ZT0TxvsZT12Tx$diffexpressed[ZT0TxvsZT12Tx$adj.P.Val < 0.05 & ZT0TxvsZT12Tx$logFC >  1] <- "Upregulated"
ZT0TxvsZT12Tx$diffexpressed[ZT0TxvsZT12Tx$adj.P.Val < 0.05 & ZT0TxvsZT12Tx$logFC < -1] <- "Downregulated"



svg('volcanoplot_ZT0TxvsZT12Tx.svg')
ggplot(data = ZT0TxvsZT12Tx, aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 15), xlim = c(-8, 8)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Legend', #legend_title
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"adj.p-value")) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('ZT0Tx vs ZT12Tx')  # Plot title
dev.off()

svg("volcano.svg")
volcanoplot(tfit, coef = 1, highlight = 20, names = rownames(tfit))
dev.off()


#heatmap of top DE genes
de.genes <- rownames(ZT0TxvsZT12Tx)[ZT0TxvsZT12Tx$adj.P.Val < 0.05 &
  abs(ZT0TxvsZT12Tx$logFC) > 1]

gene.var <- apply(y$E[de.genes, ], 1, var)
topgenes <- names(sort(gene.var, decreasing = TRUE))[1:50]

mycol <- colorpanel(1000,"blue","white","red")

svg("heatmap.svg")
heatmap.2(y$E[topgenes, ], scale = "row", labRow = rownames(y$E)[i],
  labCol = group, col = mycol, trace = "none", density.info = "none",
  margins = c(8, 6), lhei = c(2, 10), dendrogram = "both")
dev.off()


heatmap.2(
  mat_scaled,
  col = mycol,
  trace = "none",
  density.info = "none",
  scale = "none",              # already scaled
  dendrogram = "both",
  ColSideColors = col_side_colors,
  labCol = group,
  labRow = rownames(mat_scaled),
  margins = c(8, 10),
  cexCol = 0.9,
  cexRow = 0.6,
  lhei = c(1.5, 10),
  lwid = c(1.5, 10),
  main = "ZT0Tx vs ZT12Tx\nTop variable DE genes"
)

# write.table(ZT0TxvsZT12Tx, file = "ZT0TxvsZT12Tx_males.csv", quote = FALSE, sep = ",",
#                 dec = ".", row.names = TRUE,
#                 col.names = TRUE)