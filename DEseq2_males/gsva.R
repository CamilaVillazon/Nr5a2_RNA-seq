# GSVA for pathway analysis
# source: https://www.bioconductor.org/packages/devel/bioc/vignettes/GSVA/inst/doc/GSVA.html


#BiocManager::install("GSVA")
library(GSVA)

expr_for_gsva <- as.matrix(read.csv("/home/cvillazon/DEseq2_males/expr_for_gsva.csv",
                    sep=",", row.names=1, header=TRUE))

pathway_gene_sets_GO <- jsonlite::fromJSON("/home/cvillazon/DEseq2_males/pathway_gene_sets_GO.json")
#pathway_gene_sets_GO <- as.list(pathway_gene_sets_GO)

pathway_gene_sets_Wiki <- jsonlite::fromJSON("/home/cvillazon/DEseq2_males/pathway_gene_sets_Wiki.json")


gsvaPar_go <- gsvaParam(expr_for_gsva, pathway_gene_sets_GO)
gsva.es.go <- gsva(gsvaPar_go, verbose=TRUE)
write.csv(gsva.es.go, "/home/cvillazon/DEseq2_males/gsva_es_go.csv", row.names=TRUE, quote=FALSE)


gsvaPar_wiki <- gsvaParam(expr_for_gsva, pathway_gene_sets_Wiki)
gsva.es.wiki <- gsva(gsvaPar_wiki, verbose=FALSE)
write.csv(gsva.es.wiki, "/home/cvillazon/DEseq2_males/gsva_es_wiki.csv", row.names=TRUE, quote=FALSE)