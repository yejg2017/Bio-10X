library(Seurat)
library(stringr)
library(argparse)
library(destiny)
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
source("utils.R")


###############################
print("...configure parameters...")
parser <- ArgumentParser(description='Process some tasks')

parser$add_argument("--dataset",
                    type="character",
                    default=NULL,
                    help="The path of objects")

#parser$add_argument("--markers",
#                    type="character",
#                    nargs="+",
#                    default=NULL,
#                    help="markers to be analysized")


parser$add_argument("--topn",
                    type="character",
                    default="8",
                    help="number of markers to show in heatmap")

parser$add_argument("--markers",
                    type="character",
                    nargs="+",
                    default=NULL,
                    help="markers to be featureplot")


args <- parser$parse_args()

figure<-paste0("pretrain","/",args$dataset,"/","figure")
marker.dir<-paste0("markers","/",args$dataset,"/","figure")
model<-paste0("pretrain","/",args$dataset,"/","model","/","object.rds")
Cellmarkers<-paste0("pretrain","/",args$dataset,"/","model","/","Cellmarkers.rds")

if(!dir.exists(figure)){
        dir.create(figure,recursive=TRUE)
}

if(!dir.exists(marker.dir)){
        dir.create(marker.dir,recursive=TRUE)
}





print("---------------")
print(paste0("Loading object from : ",model))
#print(paste0("Loading markers from difference analysis : ",markers))
print(paste0("Plots are stored in : ",figure))
print("---------------")


print("Loading object")
object<-readRDS(model)


x<-subset(object,cells=colnames(object)[!str_detect(object@meta.data$CellType,"Unknown|Ambiguous")] )
Idents(x)<-x@meta.data$CellType

print("---PLOT---")
jpeg(paste0(figure,"/","TSNE.jpeg"))
DimPlot(x,reduction="tsne",pt.size=0.8,label=TRUE)
dev.off()


print("Idents table barplot")
jpeg(file=paste0(figure,"/","ReIdent-barplot.jpeg"))
bp<-barplot(table(Idents(x)),
            col = 1:length(table(Idents(x))),
            ylim = c(0,max(as.integer(table(Idents(x))))+500))
text(bp,as.integer(table(Idents(x)))+50,
     labels = as.character(as.integer(table(Idents(x)))))

dev.off()


jpeg(file=paste0(figure,"/","ReIdent-gbarplot.jpeg"))
Ident.counts.plot(x,outdir=figure)
dev.off()

print("Cluster FeaturePlot")
clusters<-Idents(x)
for(cluster in clusters){
}


print("Re Find Markers")
print("Find markers...")
if(!file.exists(Cellmarkers)){
	markers <- FindAllMarkers(x, only.pos =TRUE,
                          #features = VariableFeatures(object),
                          test.use = "wilcox",
                          min.pct = 0.2,
                          logfc.threshold = 0.25,
                          pseudocount.use = 1 )
	saveRDS(markers,file=Cellmarkers)
}
markers<-readRDS(Cellmarkers)

topn <- markers%>%
  group_by(cluster) %>%
  top_n(n =as.numeric(args$topn), wt = avg_logFC)

print(topn$gene)
jpeg(file=paste0(figure,"/","heatmap.jpeg"))
DoHeatmap(x,features=as.character(topn$gene),size=5.5,slot="data")+NoLegend()
#ggsave(file=paste0(figure,"/","heatmap.jpeg"),device="jpeg")
dev.off()



print("Cluster FeaturePlot")

for(i in 1:nrow(topn)){
	Cluster.FeaturePlot(x,
			    cluster=as.character(topn[["cluster"]][i]),
			    marker=as.character(topn[["gene"]][i]),
			    outdir=figure)
}


print("Cluster FeaturePlot:markers")
if(is.null(args$markers)){
	markers<-c("CCR7","MAL","TCF7","LTB","ILR7",
		   "AQP3","GZMK","CD8A","CMC1","MYOM2","PTGDS","CLCI3","CD8B",
		   "LINC02446","S100A9","AC020656.1","LYZ","KLRC2","GNLY","FGFBP2",
		   "GZMH","TRGC2","FOS","CXCL8","VCAN","KLRB1","IGKC","MS4A1","CD79A",
		   "FOLR3","TRGC1","CD14","TRDC","IGHM","IGHD","TCL1A","CDKN1C","LST1",
		   "AIF1","IL1B","CCL3","PPBP","PF4","NRGN","NEAT1","SOD2","GOS2","IL1R2",
		   "FCGR3B","IGHA1","IGLC2","IGKC","HLA-DRA1","HLA-DRB1","CAMP","CLC",
		   "DEFA3","TCF4","ITM2C","PTGDS","TYMS","STMN1","TUBA1B","SPINK2",
		   "AC084033.3","SOX4")
}


markers<-union(markers,args$markers)

markers<-markers[markers%in%rownames(x)]
for(cell in unique(Idents(x))){
	for(marker in markers){
		Cluster.FeaturePlot(x,
                            cluster=cell,
                            marker=marker,
                            outdir=marker.dir)
	}
}
