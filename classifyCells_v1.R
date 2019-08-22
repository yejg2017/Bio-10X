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


print(paste0("Loading object from : ",args$dataset))
object<-readRDS(args$dataset)

print("Convert  Seurat into Monocle")
monocle.object<-Seurat2Monocle(object)


print("Add Cell Type Message")

cth <- newCellTypeHierarchy()
cth <-addCellType(cth,cell_type_name="T Cells",
		  classify_func=function(x){ x["CCR6",] > 0 },
		  parent_cell_type_name="root")

cth <-addCellType(cth,cell_type_name="B Cells",
		  classify_func=function(x){ x["CD28",] > 0 & x["CD38",] > 0 },
		  parent_cell_type_name="root")


cth <-addCellType(cth,cell_type_name="CD8+ T Cells",
                  classify_func=function(x){ x["CD8A",] > 0 },
                  parent_cell_type_name="T Cells")

cth <-addCellType(cth,cell_type_name="CD4+ T Cells",
                  classify_func=function(x){ x["CD4",] > 0 & x["CD38",] > 0 },
                  parent_cell_type_name="T Cells")



cth <-addCellType(cth,cell_type_name="TCR alpha-gama T Cells",
                  classify_func=function(x){ x["CXCR5",] > 0 & x["HLA-DRA",] > 0 },
                  parent_cell_type_name="T Cells")

cth <-addCellType(cth,cell_type_name="CD8+ Memory T Cells",
                  classify_func=function(x){ x["HLA-DRA",] > 0 },
                  parent_cell_type_name="CD8+ T Cells")

cth <-addCellType(cth,cell_type_name="CD4+ Effector Memory T Cells",
                  classify_func=function(x){ x["HLA-DRA",] > 0 & x["CCR6",] > 0 },
                  parent_cell_type_name="CD4+ T Cells")


cth <-addCellType(cth,cell_type_name="CD4+ Navie/Effector T Cells",
                  classify_func=function(x){ x["ICOS",] > 0 },
                  parent_cell_type_name="CD4+ T Cells")

cth <-addCellType(cth,cell_type_name="CD4+ Treg-type Cells",
                  classify_func=function(x){ x["CCR6",] > 0 },
                  parent_cell_type_name="CD4+ T Cells")


cth <-addCellType(cth,cell_type_name="CD8+ Tc17-type Cells",
                  classify_func=function(x){ x["CXCR3",] > 0 },
                  parent_cell_type_name="CD8+ T Cells")

cth <-addCellType(cth,cell_type_name="NKT Cells",
                  classify_func=function(x){ x["CXCR3",] > 0 & x["CCR6",] > 0 },
                  parent_cell_type_name="root")

cth <-addCellType(cth,cell_type_name="Plasmablasts",
                  classify_func=function(x){ x["CCR6",] > 0 & x["CD28",] > 0 },
                  parent_cell_type_name="root")


cth <-addCellType(cth,cell_type_name="Basophils",
                  classify_func=function(x){ x["CCR7",] > 0 & x["LILRB1",] > 0 },
                  parent_cell_type_name="root")


cth <-addCellType(cth,cell_type_name="CXCR3+ Cells",
                  classify_func=function(x){ x["CXCR3",] > 0 },
                  parent_cell_type_name="root")


cth <-addCellType(cth,cell_type_name="CXCR3+ Dendritic Cells",
                  classify_func=function(x){ x["CXCR5",] > 0 & x["CD38",] > 0 },
                  parent_cell_type_name="CXCR3+ Cells")



monocle.object<-classifyCells(monocle.object,cth)
print("Update Seurat object")
new.meta.data<-pData(monocle.object)[,"CellType",drop=FALSE]
object<-updateSeurat(object,new.meta.data)
saveRDS(object,file=args$dataset)
print("Successfully Done")



