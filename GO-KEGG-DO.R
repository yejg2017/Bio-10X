library(clusterProfiler)
library(pathview)
library(topGO)
library(AnnotationHub)
library(biomaRt)
library(Rgraphviz)
library(org.Hs.eg.db)

library(argparse)
library(dplyr)

print("...configure parameters...")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--dataset",
                    type="character",
                    default=NULL,
                    help="The path of objects")

parser$add_argument("--topn",
		    type="character",
		    default="50")
args<-parser$parse_args()

print("reading markers")
markers<-readRDS(paste0("pretrain","/",args$dataset,"/","model","/","markers.rds"))
figure<-paste0("pretrain","/",args$dataset,"/","figure")
topn<-markers%>%
	group_by(cluster)%>%
	top_n(n=as.numeric(args$topn),wt=avg_logFC)
	      
print("enrich GO")
genes<-as.character(markers$gene)
symbols<-mapIds(x=org.Hs.eg.db,
		keys=genes,
		keytype="SYMBOL",
		column="ENTREZID")

head(symbols)
enrich.go.BP<-clusterProfiler::enrichGO(symbols,
					OrgDb=org.Hs.eg.db,
					keyType = "ENTREZID",
					ont="BP",
					pvalueCutoff=0.05,
					qvalueCutoff=0.2,
					readable=TRUE
					)

print("enrich GO BP Plot")
print("Barplot")
jpeg(paste0(figure,"/","BP.enrich.go.barplot.jpeg"))
barplot(enrich.go.BP,showCategory=15)
dev.off()

jpeg(paste0(figure,"/","BP.enrich.go.dotplot.jpeg"))
dotplot(enrich.go.BP,showCategory=15)
dev.off()

print("GO Graph")
jpeg(paste0(figure,"/","BP.GoGraph.jpeg"))
plotGOgraph(enrich.go.BP)
dev.off()

####################
enrich.go.CC<-clusterProfiler::enrichGO(symbols,
                                        OrgDb=org.Hs.eg.db,
                                        keyType = "ENTREZID",
                                        ont="CC",
                                        pvalueCutoff=0.05,
                                        qvalueCutoff=0.2,
                                        readable=TRUE
                                        )

print("enrich GO CC Plot")
print("Barplot")
jpeg(paste0(figure,"/","CC.enrich.go.barplot.jpeg"))
barplot(enrich.go.CC,showCategory=15)
dev.off()

jpeg(paste0(figure,"/","CC.enrich.go.dotplot.jpeg"))
dotplot(enrich.go.CC,showCategory=15)
dev.off()

print("GO Graph")
jpeg(paste0(figure,"/","CC.GoGraph.jpeg"))
plotGOgraph(enrich.go.CC)
dev.off()


##################
enrich.go.MF<-clusterProfiler::enrichGO(symbols,
                                        OrgDb=org.Hs.eg.db,
                                        keyType = "ENTREZID",
                                        ont="MF",
                                        pvalueCutoff=0.05,
                                        qvalueCutoff=0.2,
                                        readable=TRUE
                                        )

print("enrich GO MF Plot")
print("Barplot")
jpeg(paste0(figure,"/","MF.enrich.go.barplot.jpeg"))
barplot(enrich.go.CC,showCategory=15)
dev.off()

jpeg(paste0(figure,"/","MF.enrich.go.dotplot.jpeg"))
dotplot(enrich.go.CC,showCategory=15)
dev.off()

print("GO Graph")
jpeg(paste0(figure,"/","MF.GoGraph.jpeg"))
plotGOgraph(enrich.go.CC)
dev.off()

