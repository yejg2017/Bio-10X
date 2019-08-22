library(Seurat)
library(Signac)
library(stringr)
library(destiny)
library(argparse)
library(dplyr)
source("utils.R")

print("...configure parameters...")
parser <- ArgumentParser(description='Process some tasks')


parser$add_argument("--root",
                    type="character",
                    default="/Data/zoc/result/10X-count/10X-VDJ-human",
                    help="VDJ root")


parser$add_argument("--VHK",
                    type="character",
                    default="UVQ-7",
                    help="dataset 1 to be analysized")

parser$add_argument("--BD",
                    type="character",
                    default="UBQ-4",
                    help="dataset 2 to be analysized")




parser$add_argument("--resolution",
                    type="character",
                    default="0.8",
                    help="the resolution")

args <- parser$parse_args()

##################################
add_clonotype <- function(folder){
    x <- read.csv(paste0(folder,"filtered_contig_annotations.csv"))
    # Remove the -1 at the end of each barcode.
    # Subsets so only the first line of each barcode is kept,
    # as each entry for given barcode will have same clonotype.
    x$barcode <- gsub("-1", "", x$barcode)
    x <- x[!duplicated(x$barcode), ]

    # Only keep the barcode and clonotype columns.
    # We'll get additional clonotype info from the clonotype table.
    x <-x[,c("barcode", "raw_clonotype_id")]
    names(x)[names(x) == "raw_clonotype_id"] <- "clonotype_id"

    # Clonotype-centric info.
    clono <- read.csv(paste0(folder,"clonotypes.csv"))

    # Slap the AA sequences onto our original table by clonotype_id.
    x <- merge(x, clono[, c("clonotype_id", "cdr3s_aa")])

    # Reorder so barcodes are first column and set them as rownames.
    x <- x[, c(2,1,3)]
    rownames(x) <- paste(x[,1],"-1",sep="")
    x[,1] <- NULL

    # Add to the Seurat object's metadata.
    #clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
    return(x)
}
###################################

#dataset<-paste0(
#		paste(str_split(args$BD,"-")[[1]][1:2],collapse="-"),
#		"_",
#		paste(str_split(args$VHK,"-")[[1]][1:2],collapse="-"))

dataset<-paste0("pretrain","/",paste0(args$BD,"_",args$VHK))
if(!dir.exists(dataset)){
           dir.create(dataset,recursive=TRUE)
}

model<-paste0(dataset,"/","model")
figure<-paste0(dataset,"/","figure")

dir.create(model,recursive=TRUE)
dir.create(figure,recursive=TRUE)

###########################################

BD<-paste0(args$root,"/","5RNA","/",
	    args$BD,"-5RNA","/",
	    "outs/filtered_feature_bc_matrix.h5")

bcr.bd<-paste0(args$root,"/","BCR","/",
            args$BD,"-BCR","/",
            "outs","/")

tcr.bd<-paste0(args$root,"/","TCR","/",
            args$BD,"-TCR","/",
            "outs","/")
############################################
VHK<-paste0(args$root,"/","5RNA","/",
            args$VHK,"-5RNA","/",
            "outs/filtered_feature_bc_matrix.h5")

bcr.vhk<-paste0(args$root,"/","BCR","/",
            args$VHK,"-BCR","/",
            "outs","/")

tcr.vhk<-paste0(args$root,"/","TCR","/",
            args$VHK,"-TCR","/",
            "outs","/")
#############################################

print("---------------------")
print(paste0("Loading object from root ",args$root))
print(paste0("Loading data from data 1 :",BD))
print(paste0("Loading BCR data 1 from : ",bcr.bd))
print(paste0("Loading TCR data 1  from : ",tcr.bd))

print("---------------------")
print(paste0("Loading data from data 2 :",VHK))
print(paste0("Loading BCR data 2 from : ",bcr.vhk))
print(paste0("Loading TCR data 2  from : ",tcr.vhk))


print(paste0("Save Model in  : ",model))
print(paste0("Save Plot in : ",figure))
print("---------------------")



print("Loading object")
print("Loading object BD")
BD<-Read10X_h5(BD)
metadata.bd<-add_clonotype(tcr.bd)
BD <- CreateSeuratObject(
  counts = BD,
  assay = 'RNA',
  project = str_split(args$BD,"-")[[1]][1],
  min.cells = 200,
  min.features=1000,
  meta.data = metadata.bd
)

DefaultAssay(BD) <- "RNA"
BD<- FindVariableFeatures(BD,nfeatures=3000)
BD<-NormalizeData(BD,verbose=TRUE)
BD<-ScaleData(BD,vars.to.regress=c("nFeature_RNA"))

BD<-RenameCells(BD,add.cell.id=str_split(args$BD,"-")[[1]][1])

###############################################

print("Loading object VHK")
VHK<-Read10X_h5(VHK)
metadata.vhk<-add_clonotype(tcr.vhk)
VHK <- CreateSeuratObject(
  counts = VHK,
  assay = 'RNA',
  project =str_split(args$VHK,"-")[[1]][1],
  min.cells = 200,
  min.features=1000,
  meta.data = metadata.vhk
)

DefaultAssay(VHK) <- "RNA"
VHK<- FindVariableFeatures(VHK,nfeatures=3000)
VHK<-NormalizeData(VHK,verbose=TRUE)
VHK<-ScaleData(VHK,vars.to.regress=c("nFeature_RNA"))

VHK<-RenameCells(VHK,add.cell.id=str_split(args$VHK,"-")[[1]][1])


print("-------------------")
print("Merge VHK and BD on CCA")

f1<-VariableFeatures(BD)
f2<-VariableFeatures(VHK)
features<-intersect(f1,f2)


if(length(features)==0){
        stop("features are NULL,invalid")
}else{
        print(paste0("Use : ",length(features)," genes for CCA"))
}


object<-RunCCA(BD,VHK,
               features=features,
               num.cc=50)

print("Data preprocessing")
object[["percent.mt"]] <- PercentageFeatureSet(object,pattern = "^MT+")
object[["percent.cd"]] <- PercentageFeatureSet(object,pattern = "^CD+")


object <- FindVariableFeatures(object, selection.method = "vst",
                            nfeatures = 2000,verbose = TRUE)

object<-NormalizeData(object,normalization.method = "LogNormalize",verbose = TRUE)
object<-ScaleData(object,model.use = "linear",
               vars.to.regress = c("orig.ident",
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

print("Run Diffusionmap...")
x<-Seurat2Monocle(object)
x<-DiffusionMap(data = x,k = floor(sqrt(ncol(x))),n_eigs=3)

diffusionmap.mat<-x@eigenvectors

jpeg(paste0(dataset,"/","figure","/DiffusionMap.jpeg"))
plot(x)
dev.off()

colnames(diffusionmap.mat)<-paste("DM_",1:ncol(diffusionmap.mat),sep = "")
rownames(diffusionmap.mat)<-colnames(object)

object[["dm"]]<-CreateDimReducObject(embeddings = diffusionmap.mat,
                                  key = "DM_",
                                  assay = DefaultAssay(object)
                                  )


object<-FindNeighbors(object,reduction = "lsi",dims = 1:30)
object<-FindClusters(object,resolution = as.numeric(args$resolution))


print("Find markers...")
markers <- FindAllMarkers(object, only.pos = FALSE,
                          features = VariableFeatures(object),
                          test.use = "wilcox",
                          min.pct = 0.2,
                          logfc.threshold = 0.25,
                          pseudocount.use = 1 )


print("Saving Model :")
saveRDS(x,paste0(dataset,"/","model","/","DM.rds"))
saveRDS(object,paste0(dataset,"/","model","/","object.rds"))
saveRDS(markers,paste0(dataset,"/","model","/","markers.rds"))
print("Successfully Done")

