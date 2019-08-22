library(Seurat)
library(Signac)
library(stringr)
library(argparse)
library(destiny)
library(dplyr)
source("utils.R")

print("...configure parameters...")
parser <- ArgumentParser(description='Process some tasks')


parser$add_argument("--root",
                    type="character",
                    default="/Data/zoc/result/10X-count/PBMC/10X-VDJ-human/5RNA",
                    help="root of dataset")

parser$add_argument("--ratio",
		    type="character",
		    default="0.1",
		    help="the ratio to sample each dataset")

parser$add_argument("--dm",
		    dest="dm",
		    action="store_true")


parser$add_argument("--seed",
		    type="integer",
		    default="1234",
		    help="the seed for random sample")


parser$add_argument("--resolution",
                    type="character",
                    default="0.8",
                    help="the resolution")

args <- parser$parse_args()

print("--------------------")
print(paste0("The root of dataset is : ",args$root))
print(paste0("Whether run diffusionmap : ",args$dm))
print(paste0("The ratio of random sample for each dataset is :",args$ratio))
print(paste0("The seed for random sample is ",args$seed))
print("--------------------")

dataset<-paste0("pretrain","/","all","_ratio_",args$ratio,"_seed_",args$seed,"_resolution_",args$resolution)
model<-paste0(dataset,"/","model")
figure<-paste0(dataset,"/","figure")
if(!dir.exists(model)){
           dir.create(model,recursive=TRUE)
}

if(!dir.exists(figure)){
           dir.create(figure,recursive=TRUE)
}


print("Loading dataset")

#system("ls /Data/zoc/result/10X-count/PBMC/10X-VDJ-human/5RNA | grep -v -E '*-RNA-*' > sample.txt")
samples<-as.character(read.table("sample0.txt",header=FALSE)$V1)
classes<-as.character(read.table("sample0.txt",header=FALSE)$V4)
data<-list()
idents<-c()
class<-c()
i<-1
for(ident in samples){
   file=paste0(args$root,"/",ident,"/outs/filtered_feature_bc_matrix.h5")
   x<-Read10X_h5(file)
   n<-floor(ncol(x)*as.numeric(args$ratio))
   print(paste0("The No.",i," ",ident," sample size is : ",n))
   idents<-c(idents,rep(paste(str_split(ident,"-")[[1]][1:2],collapse="-"),n))
   class<-c(class,rep(classes[i],n))
   idx<-sample(1:ncol(x),size=n,replace=FALSE)
   x<-x[,idx]
   colnames(x)<-paste(paste(str_split(ident,"-")[[1]][1:2],collapse="-"),colnames(x),sep="_")
   data[[i]]<-x
   i<-i+1
}

counts<-do.call(cbind,data)
print(dim(counts))
print(table(idents))

print("remove MT- genes")
genes.use<-rownames(counts)[!str_detect(rownames(counts),"^MT-+")]
counts<-counts[genes.use,]

print("Create Seurat Object")
idents<-data.frame("ident"=idents,"class"=class)
rownames(idents)<-colnames(counts)
print(dim(idents))
object<-CreateSeuratObject(counts= counts,
                       assay = "RNA",
                       project ="scRNA",
		       meta.data=idents,
                       min.cells=250,
                       min.features=100)


object@meta.data$disease<-unlist(lapply(object@meta.data$ident,function(i){return(str_split(i,"-")[[1]][1])}))
#object@meta.data$class<-ifelse(str_detect(object@meta.data$disease,"UBQ|UBR"),"BD","VHK")
print(head(object@meta.data))

print("Data preprocessing")
object[["percent.mt"]] <- PercentageFeatureSet(object,pattern = "^MT+")
object[["percent.cd"]] <- PercentageFeatureSet(object,pattern = "^CD+")

print("Subset...")
jpeg(paste0(figure,"/VlnPlot.jpeg"))
p1<-VlnPlot(object, features = c("nFeature_RNA"),group.by="disease",pt.size=0.3)
p2<-VlnPlot(object, features = c("percent.cd"),group.by="disease",pt.size=0.3)
p3<-VlnPlot(object, features = c("percent.mt"),group.by="disease",pt.size=0.3)

CombinePlots(plots = list(p1, p2,p3))
dev.off()


object<- subset(object, subset = nFeature_RNA > 500 & nFeature_RNA < 3800& percent.mt <25&percent.cd<2)


object <- FindVariableFeatures(object, selection.method = "vst",
                            nfeatures = 2000,verbose = TRUE)

object<-NormalizeData(object,normalization.method = "LogNormalize",verbose = TRUE)
object<-ScaleData(object,model.use = "linear",
               vars.to.regress = c("ident",
                                   "nFeature_RNA",
                                   "percent.mt"),verbose = TRUE)



print("Run LSI...")
object <- RunLSI(object, n = 50, scale.max = NULL)
print("Run UMAP...")
object <- RunUMAP(object, reduction = "lsi", dims = 1:30)

print("Run PCA...")
object<-RunPCA(object,assay = "RNA",npcs = 50)

print("Run TSNE...")
object<-RunTSNE(object,reduction="lsi",dims=1:30)


if(args$dm){
	print("Run Diffusionmap...")
        x<-Seurat2Monocle(object)
        x<-DiffusionMap(data = x,k = floor(sqrt(ncol(x))),n_eigs=3)

        diffusionmap.mat<-x@eigenvectors

        jpeg(paste0(figure,"/DiffusionMap.jpeg"))
        plot(x)
        dev.off()

        colnames(diffusionmap.mat)<-paste("DM_",1:ncol(diffusionmap.mat),sep = "")
        rownames(diffusionmap.mat)<-colnames(object)

        object[["dm"]]<-CreateDimReducObject(embeddings = diffusionmap.mat,
                                  key = "DM_",
                                  assay = DefaultAssay(object)
                                  )
	saveRDS(x,paste0(model,"/","DM.rds"))
}


object<-FindNeighbors(object,reduction = "lsi",dims = 1:30)
object<-FindClusters(object,resolution = as.numeric(args$resolution))


print("Find markers...")
markers <- FindAllMarkers(object, only.pos = FALSE,
                          #features = VariableFeatures(object),
                          test.use = "wilcox",
                          min.pct = 0.2,
                          logfc.threshold = 0.25,
                          pseudocount.use = 1 )


print("Saving Model :")
saveRDS(object,paste0(model,"/","object.rds"))
saveRDS(markers,paste0(model,"/","markers.rds"))
print("Successfully Done")



