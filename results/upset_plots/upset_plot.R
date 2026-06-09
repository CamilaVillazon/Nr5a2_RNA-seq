library(UpSetR)

df_for_upset <- read.table("/home/cvillazon/scripts/df_for_upset_filtered.csv", header = TRUE, sep = ",",
                            row.names = 1)
df_for_upset[df_for_upset == "False"] <- 0
df_for_upset[df_for_upset == "True"] <- 1
df_for_upset <- as.data.frame(lapply(df_for_upset, as.numeric))

pdf("/home/cvillazon/scripts/upset_plot_all.pdf", width = 14, height = 6)
upset(df_for_upset, nsets = ncol(df_for_upset), order.by = "freq",
        nintersects = NA, point.size = 1, line.size = 0.1,
        show.numbers = "none")
dev.off()

pdf("/home/cvillazon/scripts/upset_plot_biggest.pdf", width = 14, height = 6)
upset(df_for_upset, nsets = ncol(df_for_upset), order.by = "freq",
        nintersects = 99, point.size = 1, line.size = 0.1,
        show.numbers = "yes")
dev.off()

pdf("/home/cvillazon/scripts/upset_plot_smallest.pdf", width = 14, height = 6)
upset(df_for_upset, nsets = ncol(df_for_upset), order.by = "freq",
        point.size = 1, line.size = 0.1, decreasing = FALSE,
        nintersects = 204, show.numbers = "yes")
dev.off()


combinations <- list(c("ZT0Tx.ZT12Tx_male", "ZT0Tx.ZT12Tx_female"), 
c("ZT0Tx.ZT12Tx_male", "ZT0V.ZT0Tx_female"), 
c("ZT0Tx.ZT12Tx_male", "ZT0V.ZT12Tx_female"), 
c("ZT0Tx.ZT12Tx_male", "ZT0V.ZT12V_female"), 
c("ZT0Tx.ZT12Tx_male", "ZT12V.ZT0Tx_female"), 
c("ZT0Tx.ZT12Tx_male", "ZT12V.ZT12Tx_female"), 
c("ZT0V.ZT0Tx_male", "ZT0Tx.ZT12Tx_female"), 
c("ZT0V.ZT0Tx_male", "ZT0V.ZT0Tx_female"), 
c("ZT0V.ZT0Tx_male", "ZT0V.ZT12Tx_female"), 
c("ZT0V.ZT0Tx_male", "ZT0V.ZT12V_female"), 
c("ZT0V.ZT0Tx_male", "ZT12V.ZT0Tx_female"), 
c("ZT0V.ZT0Tx_male", "ZT12V.ZT12Tx_female"), 
c("ZT0V.ZT12Tx_male", "ZT0Tx.ZT12Tx_female"), 
c("ZT0V.ZT12Tx_male", "ZT0V.ZT0Tx_female"), 
c("ZT0V.ZT12Tx_male", "ZT0V.ZT12Tx_female"), 
c("ZT0V.ZT12Tx_male", "ZT0V.ZT12V_female"), 
c("ZT0V.ZT12Tx_male", "ZT12V.ZT0Tx_female"), 
c("ZT0V.ZT12Tx_male", "ZT12V.ZT12Tx_female"), 
c("ZT0V.ZT12V_male", "ZT0Tx.ZT12Tx_female"), 
c("ZT0V.ZT12V_male", "ZT0V.ZT0Tx_female"), 
c("ZT0V.ZT12V_male", "ZT0V.ZT12Tx_female"), 
c("ZT0V.ZT12V_male", "ZT0V.ZT12V_female"), 
c("ZT0V.ZT12V_male", "ZT12V.ZT0Tx_female"), 
c("ZT0V.ZT12V_male", "ZT12V.ZT12Tx_female"), 
c("ZT12V.ZT0Tx_male", "ZT0Tx.ZT12Tx_female"), 
c("ZT12V.ZT0Tx_male", "ZT0V.ZT0Tx_female"), 
c("ZT12V.ZT0Tx_male", "ZT0V.ZT12Tx_female"), 
c("ZT12V.ZT0Tx_male", "ZT0V.ZT12V_female"), 
c("ZT12V.ZT0Tx_male", "ZT12V.ZT0Tx_female"), 
c("ZT12V.ZT0Tx_male", "ZT12V.ZT12Tx_female"), 
c("ZT12V.ZT12Tx_male", "ZT0Tx.ZT12Tx_female"), 
c("ZT12V.ZT12Tx_male", "ZT0V.ZT0Tx_female"), 
c("ZT12V.ZT12Tx_male", "ZT0V.ZT12Tx_female"), 
c("ZT12V.ZT12Tx_male", "ZT0V.ZT12V_female"), 
c("ZT12V.ZT12Tx_male", "ZT12V.ZT0Tx_female"), 
c("ZT12V.ZT12Tx_male", "ZT12V.ZT12Tx_female"))

pdf("/home/cvillazon/scripts/upset_plot_pairs.pdf", width = 14, height = 6)
upset(df_for_upset, nsets = ncol(df_for_upset), order.by = "freq",
        point.size = 1, line.size = 0.1, show.numbers = "yes",
        intersections = combinations)
dev.off()