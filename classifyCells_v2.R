library(Seurat)
library(monocle)
library(stringr)
library(argparse)
source("utils.R")


print("...configure parameters...")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--dataset",
                    type="character",
                    default=NULL,
                    help="The path of objects")

args<-parser$parse_args()

updateSeurat<-function(object,meta.data=NULL){
	object<-AddMetaData(object,
			    metadata=meta.data)
	return(object)
}


model<-paste0("pretrain","/",args$dataset,"/","model","/","object.rds")
print(paste0("Loading object from : ",args$dataset))
object<-readRDS(model)

print("Convert  Seurat into Monocle")
monocle.object<-Seurat2Monocle(object)


print("Add Cell Type Message")

cth <- newCellTypeHierarchy()
cth <-addCellType(cth,cell_type_name="B Cells",
		  classify_func=function(x){ x["PTPRC",]>0 & x["CD19",] > 0 },
		  parent_cell_type_name="root")

cth <-addCellType(cth,cell_type_name="T Cells",
		  classify_func=function(x){ x["PTPRC",] > 0 & x["CD3D",] > 0 },
		  parent_cell_type_name="root")


cth <-addCellType(cth,cell_type_name="CD8+ T Cells",
                  classify_func=function(x){ x["CD8A",] > 0 },
                  parent_cell_type_name="T Cells")

cth <-addCellType(cth,cell_type_name="CD4+ T Cells",
                  classify_func=function(x){ x["CD4",] > 0 },
                  parent_cell_type_name="T Cells")



cth <-addCellType(cth,cell_type_name="Treg Cells",
                  classify_func=function(x){ x["FOXP3",] > 0},
                  parent_cell_type_name="CD4+ T Cells")

cth <-addCellType(cth,cell_type_name="Myeloid Cells",
                  classify_func=function(x){ x["PTPRC",] > 0 & x["CD19",]==0 & x["CD3D",]==0 & x["NCAM1",]==0 & x["ITGAM",] > 0 },
                  parent_cell_type_name="root")



monocle.object<-classifyCells(monocle.object,cth)
print("Update Seurat object")
new.meta.data<-pData(monocle.object)[,"CellType",drop=FALSE]
object<-updateSeurat(object,new.meta.data)
saveRDS(object,model)
print("Successfully Done")



