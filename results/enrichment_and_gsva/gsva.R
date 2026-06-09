# GSVA for pathway analysis
# source: https://www.bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html


#BiocManager::install("GSVA")
library(GSVA)

expr_for_gsva <- as.matrix(read.csv("/home/cvillazon/scripts/f_unique_tx_for_gsva.csv",
                    sep=",", row.names=1, header=TRUE))
pathway_gene_sets_Wiki <- jsonlite::fromJSON("/home/cvillazon/scripts/f_unique_tx.json")
gsvaPar_wiki <- gsvaParam(expr_for_gsva, pathway_gene_sets_Wiki)
gsva.es.wiki <- gsva(gsvaPar_wiki, verbose=FALSE, minSize=5)
write.csv(gsva.es.wiki, "/home/cvillazon/scripts/f_unique_tx_gsva_scores.csv", row.names=TRUE, quote=FALSE)


expr_for_gsva <- as.matrix(read.csv("/home/cvillazon/scripts/f_unique_v_for_gsva.csv",
                    sep=",", row.names=1, header=TRUE))
pathway_gene_sets_Wiki <- jsonlite::fromJSON("/home/cvillazon/scripts/f_unique_v.json")
gsvaPar_wiki <- gsvaParam(expr_for_gsva, pathway_gene_sets_Wiki)
gsva.es.wiki <- gsva(gsvaPar_wiki, verbose=FALSE, minSize=5)
write.csv(gsva.es.wiki, "/home/cvillazon/scripts/f_unique_v_gsva_scores.csv", row.names=TRUE, quote=FALSE)


expr_for_gsva <- as.matrix(read.csv("/home/cvillazon/scripts/m_unique_tx_for_gsva.csv",
                    sep=",", row.names=1, header=TRUE))
pathway_gene_sets_Wiki <- jsonlite::fromJSON("/home/cvillazon/scripts/m_unique_tx.json")
gsvaPar_wiki <- gsvaParam(expr_for_gsva, pathway_gene_sets_Wiki)
gsva.es.wiki <- gsva(gsvaPar_wiki, verbose=FALSE, minSize=5)
write.csv(gsva.es.wiki, "/home/cvillazon/scripts/m_unique_tx_gsva_scores.csv", row.names=TRUE, quote=FALSE)


expr_for_gsva <- as.matrix(read.csv("/home/cvillazon/scripts/m_unique_v_for_gsva.csv",
                    sep=",", row.names=1, header=TRUE))
pathway_gene_sets_Wiki <- jsonlite::fromJSON("/home/cvillazon/scripts/m_unique_v.json")
gsvaPar_wiki <- gsvaParam(expr_for_gsva, pathway_gene_sets_Wiki)
gsva.es.wiki <- gsva(gsvaPar_wiki, verbose=FALSE, minSize=5)
write.csv(gsva.es.wiki, "/home/cvillazon/scripts/m_unique_v_gsva_scores.csv", row.names=TRUE, quote=FALSE)
