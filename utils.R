library(Seurat)
library(ggplot2)
library(ggpubr)
library(stringr)


# PBMC markers
pbmc.human.markers<-c("CD19",# B Cells",
                      "CD3","CD56", # T Cells",
                      "CD25","CD69",# Actived T Cells",
                      "CD1C","CD83","CD141","CD209", #Dendritic Cells* ",
                      "CD123","CD303","CD304", #Plasmacytoid Dendritic Cells* ",
                      "CD42B", #Platelets (resting)",
                      "CD62P", # Platelets (active)",
                      "CD3","CD16","CD56", #Natural Killer Cells",
                      "CD34","CD90",#Hematopoietic Stem Cell",
                      "CD11B","CD68","CD163",#Macrophage* ",
                      "CD14","CD16","CD64", #Monocyte* ",
                      "CD138",#Plasma Cells",
                      "CD235A",# Red Blood Cells",
                      "CD15","CD16","CD49D", # Neutrophils",
                      "CD117","CD123","CD203C",#Basophils",
                      "CD11B","CD193","SIGLEC-8" #Eosinophils"
)

pbmc.mouse.markers<-c("CD19",# B Cells",
                      "CD3","CD49B", # T Cells",
                      "CD25","CD69",# Actived T Cells",
                      "CD11C", #Dendritic Cells* ",
                      "CD11C","CD317", #Plasmacytoid Dendritic Cells* ",
                      "CD41", #Platelets (resting)",
                      "CD62P", # Platelets (active)",
                      "CD3","CD49B", #Natural Killer Cells",
                      "CD48","CD117","CD150","SCA-1",#Hematopoietic Stem Cell",
                      "F4","F80","CD68",#Macrophage* ",
                      "CD11B","CD115","GR-1", #Monocyte* ",
                      "CD138",#Plasma Cells",
                      "TER-119",# Red Blood Cells",
                      "CD11B","CD115","LY6G","GR-1", # Neutrophils",
                      "CD200R3",#Basophils",
                      "CD11B","CD193","SIGLEC-F","F4","F80" #Eosinophils"
)

must.markers<-c("CD3","CD4","CD8A","CD19","CD45RA","CD45RO")


cell.surface.markers<-c("HLA","DR","CD11C","CD14","CD16",
                "CD123","CD11B","CD33","CXCR3",
                "CD3E","CD3G","CD3D","CD3","CD19",
                "LGD","CD27","CD38","CD24","CD36",
                "CD4","CD8","CD45RA","CCR7","CD25","CD127",
                "CXCR5","CD279","CD161","CD57","CD56",
                "CD94")


###############################################
NK.Cells<-c("NCAM1","NCR1","NKG2")
Cytotoxic.T<-c("GNLY","PFN1","GZMA","GZMB","GMZM","GZMH")
Exhausted.T<-c("FOXP3", "CTLA4", "TIGIT","TNFRSF4", "LAG3", "PDCD1")
T.Cells<-c("CD8A","CD4")
T.Cells.CD4<-"CD4"
T.Cells.CD8<-"CD8A"
Naive.T.cells<-c("IL7R")
B.Cells<-c("CD19")
Mast.Cells<-c("ENPP3","KIT")
Plasmacytoid.DC<-c("IL3RA","LILRA4")
Monocytic.Lineage<-c("HLA-DRB1","FCGR3A", "CD68", 
                     "ANPEP", "ITGAX", "CD14","ITGAM", "CD33")
Dendritic.Cell<-c("FCER1A")
Memory.B.Cells<-c("CD38","MS4A1")

CD11B<-"ITGAM"
CD11C<-"ITGAX"
CD20<-"MS4A1"

############################################# caculate gene variance
require(ggplot2)
Clusters.CV<-function(object=NULL,use.ident=TRUE,cv=TRUE,
                      genes=NULL,slot="data",eps=1e-4,return.cv=FALSE){
  
  if(slot=="count"){
    data<-object[["RNA"]]@counts
  }
  else if(slot=="data"){
    data<-object[["RNA"]]@data
  }
  else{
    data<-object[["RNA"]]@scale.data
  }
  
  if(is.null(genes)){
    genes<-VariableFeatures(object)
  }
  data<-data[genes,]
  print(dim(data))
  
  if(use.ident){
    idents<-unique(Idents(object))
    cv<-lapply(idents,function(i){
      cells<-WhichCells(object,idents=i)
      mat<-data[,cells]
      u<-apply(mat,1,mean)
      std<-apply(mat,1,sd)
      if(cv){
        x<-data.frame(std/(u+eps))
        colnames(x)<-paste0("Cluster_",i)
        return(x)
      }else{
        x<-data.frame(std)
        colnames(x)<-paste0("Cluster_",i)
        return(x)
      }
      
    })}
  else{
    idents<-unique(object@meta.data$orig.ident)
    cv<-lapply(idents,function(i){
      cells<-colnames(object)[object@meta.data$orig.ident==i]
      mat<-data[,cells]
      u<-apply(mat,1,mean)
      std<-apply(mat,1,sd)
      if(cv){
        x<-data.frame(std/(u+eps))
        colnames(x)<-paste0("Cluster_",i)
        return(x)
      }else{
        x<-data.frame(std)
        colnames(x)<-paste0("Cluster_",i)
        return(x)
      }
    })
  }
  cv<-as.data.frame(cv)
  cv<-reshape2::melt(cv)
  
  print(ggplot(data=cv,aes(x=value,fill=variable))+
    geom_density()+theme_bw()+xlab("gene")+ylab("pdf")+ggtitle("Gene Expression Variance")+
    theme(axis.text = element_text(face = "bold",size = 10),
                                    axis.title = element_text(face = "bold",size=14),
                                    legend.title = element_blank(),
                                    legend.text = element_text(face = "bold",size=14)))
  
  if(return.cv){
   return(cv) 
  }else{
    return(NULL)
  }
  
  
}

Cell.Chart.Ratio<-function(object,idents=NULL){
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  n.cells<-ncol(object)
  if(is.null(idents)){
    idents<-unique(Idents(object))
  }
  ratio.idents<-lapply(idents,function(i){
    len.ident<-list(round(length(WhichCells(object,idents =i))/n.cells,3))
    names(len.ident)<-paste0("Cluster_",i)
    return(len.ident)
  })
  
  count.data<-melt(as.data.frame(ratio.idents))
  
  # Add label positio
  count.data <- count.data %>%
    arrange(desc(variable)) %>%
    mutate(lab.ypos = cumsum(value) - 0.5*value)
  
  
  ggplot(count.data, aes(x = "", y = round(value,3), fill = variable)) +
    geom_bar(width = 100, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    geom_text(aes(y = lab.ypos, label =value), color = "black",check_overlap = TRUE)+
    theme(legend.title = element_blank())
  
}


StackCSV<-function(csv.dir,sequence="TCR",type="contig",object=NULL){
  files<-list.files(csv.dir,recursive = TRUE,include.dirs = FALSE,full.names = TRUE)
  
  if(type=="contig"){
    target=paste0(sequence,"_","all_contig_annotations.csv")
  }
  else if(type=="consensus"){
    target=paste0(sequence,"_","consensus_annotations.csv")
  }
  else{
    target=paste0(sequence,"_","clonotypes.csv")
  }
  csv<-data.frame()
  print("...Stack...")
  for(file in files){
    if(target==basename(file)){
      df<-read.csv(file,header = TRUE,encoding = "utf-8")
      csv<-rbind(csv,df)
    }
  }
  print("...Done...")
  
  if(type=="contig"||type=="consensus"){
    
    csv$cell<-unlist(
      lapply(csv$barcode,function(v){
        return(str_split(v,"-")[[1]][1])
      })
    )
  }
  if(!is.null(object)){
    csv<-csv[csv$cell%in%colnames(object),]
  }
  
  return(csv)
}


Add.Feature<-function(object,tcr,feature="chain"){
  if(!feature%in%colnames(tcr)){
    stop("Feature is not valid!")
  }
  add<-c()
  n<-nrow(tcr)
  
  feature<-tcr[,feature]
  if(is.factor(feature)){
    feature<-as.character(feature)
  }else{
    feature<-as.numeric(feature)
  }
  tcr.barcode<-tcr[,"cell"]
  barcode<-colnames(object)
  for(i in 1:length(barcode)){
    f<-ifelse(barcode[i]%in%tcr.barcode,feature[i],
              ifelse(is.character(feature[i]),"Unknow",0))
    add<-c(add,f)
  }
  #object@meta.data[,feature]<-add
  return(add)
}



require(iCellR)
vdj<-function(data.dir,type="BCR",feature="UBQ-3"){
  data.dir<-paste(data.dir,feature,paste0(type,"_","all_contig_annotations.csv"),sep = "/")
  vdj<-prep.vdj(vdj.data =data.dir, cond.name =feature)
  return(vdj)
}

sample.idx<-function(df,ratio=0.5){
  n.features<-dim(df)[2]
  n.sample<-floor(ratio*n.features)
  idx<-sample(1:n.features,size=n.sample,replace=F)
  return(idx)
}


Barcode<-function(vdj.data,sep="-"){
  barcode=vdj.data$barcode
  barcode=unlist(
    lapply(barcode,function(b){
      bs=str_split(b,"-")[[1]]
      bs=paste(bs[1],bs[2],sep=sep)
      return(bs)}))
  
  #rownames(vdj.data)<-barcode
  return(barcode)
}




Seurat.Clono.Plot<-function(object,vdj.data,plot.data.type="tsne",
                            clono=1,interactive=FALSE,cell.size = 1, 
                            cell.colors = c("red", "gray"), box.cell.col = "black",
                            back.col = "white", cell.transparency = 0.5,out.name="plot"){
  if(plot.data.type=="tsne"){
    MyTitle<-"TSNE Plot"
    DATA<-as.data.frame(object@reductions$tsne@cell.embeddings)
  }
  if(plot.data.type=="pca"){
    MyTitle<-"PCA Plot"
    DATA<-as.data.frame(object@reductions$pca@cell.embeddings)
  }
  if(plot.data.type=="umap"){
    MyTitle<-"UMAP Plot"
    DATA<-as.data.frame(object@reductions$umap@cell.embeddings)
  }
  
  colono <- unique(vdj.data[,c("raw_clonotype_id","barcode")])
  rownames(colono)<-Barcode(colono,sep = ".")
  
  orig.ident<-as.character(object@meta.data$orig.ident)
  orig.ident<-unlist(lapply(orig.ident,function(x){
    x<-str_split(x,"-")[[1]]
    x<-paste(x[1],x[2],sep = ".")
  }))
  barcode<-rownames(DATA)
  
  barcode<-unlist(lapply(1:length(barcode),function(i){
    return(paste(orig.ident[i],barcode[i],sep = "_"))
  }))
  

  rownames(DATA)<-barcode
  
  colono$raw_clonotype_id <- gsub("clonotype", " ", colono$raw_clonotype_id)
  colono <- colono[1]
  colnames(colono) <- c("Clonotypes")

  colonoData <- merge(DATA, colono, by = "row.names", all.x = T,
                      all.y = F)
  
  #clono<-paste0("clonotype",clono)
  colonoData$Clonotypes <- gsub(" ", "", colonoData$Clonotypes)
  colonoData$Clonotypes[is.na(colonoData$Clonotypes)] <- "NA"
  
  #top.clono<-names(table(colonoData$Clonotypes)[order(table(colonoData$Clonotypes),decreasing = TRUE)][1:5])

  colonoData$Clonotypes[colonoData$Clonotypes != clono] <- "NA"
  colonoData$Clonotypes[colonoData$Clonotypes == clono] <- clono
  
  #print(table(colonoData$Clonotypes))
  
  DATA <- colonoData
  row.names(DATA) <- DATA$Row.names
  DATA <- DATA[, -1]
  DATA <- (DATA[order(DATA$Clonotypes, decreasing = T), ])
  DATA$Clonotypes<-as.factor(DATA$Clonotypes)
  
  if (interactive == F) {
    myPLOT <- ggplot(DATA, aes(DATA[, 1], y = DATA[,
                                                   2], col = Clonotypes, text = row.names(DATA))) +
      geom_point(size = cell.size, alpha = cell.transparency) +
      xlab("Dim1") + ylab("Dim2") + ggtitle(paste(MyTitle,
                                                  "(clonotype", clono, ")")) + theme(legend.position = "none") +
      scale_color_manual(values = cell.colors) + theme(panel.background = element_rect(fill = back.col,
                                                                                       colour = "black"), panel.grid.major = element_blank(),
                                                       panel.grid.minor = element_blank(), legend.key = element_rect(fill = back.col))
  }
  else {
    myPLOT <- ggplot(DATA, aes(DATA[, 1], y = DATA[,
                                                   2], text = row.names(DATA), color = Clonotypes)) +
      geom_point(size = cell.size, alpha = cell.transparency) +
      scale_color_manual(values = cell.colors) + xlab("Dim1") +
      ylab("Dim2") + ggtitle(MyTitle) + theme_bw()
  }

  if (interactive == T) {
    OUT.PUT <- paste(out.name, ".html", sep = "")
    htmlwidgets::saveWidget(ggplotly(myPLOT), OUT.PUT)
  }
  else {
    return(myPLOT)
  }
}


require(factoextra)
opt.clust.num<-function (x = NULL, max.clust = 20, gap.stat.nboot = 100, verbose = TRUE, 
          clust.type = "tsne", opt.method = "silhouette") 
{
  if ("scSeqR" != class(x)[1] && "iCellR" != class(x)[1]) {
    stop("x should be an object of class scSeqR")
  }
  if (clust.type == "tsne") {
    df <- (x@tsne.data)[1:2]
  }
  if (clust.type == "pca") {
    df <- (x@pca.data)[1:2]
  }
  if (clust.type == "distance") {
    df <- x@dist.data
  }
  if (opt.method == "all") {
    Elbow = fviz_nbclust(df, kmeans, method = "wss", k.max = max.clust) + 
      geom_vline(xintercept = 4, linetype = 2) + labs(subtitle = "Elbow method")
    Silhouette = fviz_nbclust(df, kmeans, method = "silhouette", 
                              k.max = max.clust) + labs(subtitle = "Silhouette method")
    set.seed(123)
    Gap.statistic = fviz_nbclust(df, kmeans, nstart = 25, 
                                 method = "gap_stat", nboot = gap.stat.nboot, k.max = max.clust) + 
      labs(subtitle = "Gap statistic method")
    return(grid.arrange(Elbow, Silhouette, Gap.statistic, 
                        nrow = 3))
  }
  if (opt.method == "elbow.wss") {
    Elbow = fviz_nbclust(df, kmeans, method = "wss", k.max = max.clust) + 
      geom_vline(xintercept = 4, linetype = 2) + labs(subtitle = "Elbow method")
    return(Elbow)
  }
  if (opt.method == "silhouette") {
    geom_vline(xintercept = 4, linetype = 2)
    Silhouette = fviz_nbclust(df, kmeans, method = "silhouette", 
                              k.max = max.clust) + labs(subtitle = "Silhouette method")
    return(Silhouette)
  }
  if (opt.method == "gap.stat") {
    geom_vline(xintercept = 4, linetype = 2)
    Gap.statistic = fviz_nbclust(df, kmeans, nstart = 25, 
                                 method = "gap_stat", nboot = gap.stat.nboot, k.max = max.clust) + 
      labs(subtitle = "Gap statistic method")
    return(Gap.statistic)
  }
}


cluster.cond.info<-function (x = NULL, plot.type = "pie",split=TRUE,
                             normalize.ncell = TRUE) 
{
  if ("Seurat" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }

  MYConds<-unique(as.character(x@meta.data[,"orig.ident"]))
  if(split){
    MYConds<-unlist(lapply(MYConds,function(s){return(str_split(s,"-")[[1]][1])}))
  }
  if(length(MYConds)<2){
    stop("You need more then one condition/sample to run this function")
  }

  DATA <-x@meta.data[,"seurat_clusters",drop=FALSE]
  rownames(DATA)<-rownames(x@meta.data)
  colnames(DATA)<-"clusters"
  Conds <- as.character(x@meta.data[,"orig.ident"])
  if(split){
    Conds<-unlist(lapply(Conds,function(s){return(str_split(s,"-")[[1]][1])}))
  }
  
  ForNorm1 <- as.data.frame(table(Conds))
  ForNorm <- min(ForNorm1$Freq)
  clusts <- (as.data.frame(DATA$clusters))
  cond.clust <- cbind(Conds, clusts)
  colnames(cond.clust) <- c("conditions", "clusters")
  Conds <- as.character(ForNorm1$Conds)
  My.Conds.data <- cond.clust
  for (i in Conds) {
    NameCol <- paste("My_Cond", i, sep = "_")
    myDATA <- head(subset(My.Conds.data, My.Conds.data$conditions == 
                            i), ForNorm)
    eval(call("<-", as.name(NameCol), myDATA))
  }
  filenames <- ls(pattern = "My_Cond_")
  datalist <- mget(filenames)
  NormDATA <- do.call(rbind.data.frame, datalist)
  if (normalize.ncell == T) {
    cond.clust <- NormDATA
  }
  DATA <- as.data.frame(table(cond.clust))
  Freq <- DATA$Freq
  if (normalize.ncell == T) {
    colnames(DATA) <- c("conditions", "clusters", "NormalizedFreq")
  }
  myBP <- ggplot(DATA, aes(y = Freq, x = clusters, fill = conditions)) + 
    geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90))
  myPIE <- ggplot(DATA, aes(y = Freq, x = "", fill = conditions)) + 
    geom_bar(stat = "identity", position = "fill") + theme_bw() + 
    facet_wrap(~clusters) + theme(axis.title.y = element_blank(), 
                                  axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    coord_polar(theta = "y")
  write.table((DATA), file = "clust_cond_freq_info.txt", sep = "\t", 
              row.names = F)
  print("clust_cond_freq_info.txt file has beed generated.")
  if (plot.type == "bar") {
    return(myBP)
  }
  if (plot.type == "pie") {
    return(myPIE)
  }
}


cond.cluster.info<-function(x=NULL,plot.type="bar",split=TRUE,normalize.ncell=TRUE){
  if ("Seurat" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  clusters<-Idents(x)
  MYConds<-x@meta.data$orig.ident
  
  if(split){
    MYConds<-unlist(lapply(MYConds,function(s){return(str_split(s,"-")[[1]][1])}))
  }
  if(length(unique(MYConds))<2){
    stop("You need more then one condition/sample to run this function")
  }
  
  DATA <-x@meta.data[,"seurat_clusters",drop=FALSE]
  rownames(DATA)<-rownames(x@meta.data)
  colnames(DATA)<-"clusters"
  Conds<-MYConds
  
  ForNorm1 <- as.data.frame(table(clusters))
  ForNorm <- min(ForNorm1$Freq)
  clusts <- (as.data.frame(DATA$clusters))
  cond.clust <- cbind(Conds, clusts)
  colnames(cond.clust) <- c("conditions", "clusters")
  
  clusters <- as.character(ForNorm1$clusters)
  My.Conds.data <- cond.clust
  
  for (i in clusters) {
    NameCol <- paste("My_Cluster", i, sep = "_")
    myDATA <- head(subset(My.Conds.data, My.Conds.data$clusters == 
                            i), ForNorm)
    eval(call("<-", as.name(NameCol), myDATA))
  }
  filenames <- ls(pattern = "My_Cluster_")
  datalist <- mget(filenames)
  NormDATA <- do.call(rbind.data.frame, datalist)
  if (normalize.ncell == T) {
    cond.clust <- NormDATA
  }
  DATA <- as.data.frame(table(cond.clust))
  Freq <- DATA$Freq
  if (normalize.ncell == T) {
    colnames(DATA) <- c("conditions", "clusters", "NormalizedFreq")
  }
  myBP <- ggplot(DATA, aes(y = Freq, x = conditions, fill = clusters)) + 
    geom_bar(stat = "identity",position = position_dodge(width = 0.8)) + theme_bw() + theme(axis.text.x = element_text(angle = 90))
  myPIE <- ggplot(DATA, aes(y = Freq, x = "", fill = clusters)) + 
    geom_bar(stat = "identity", position = "fill") + theme_bw() + 
    facet_wrap(~conditions) + theme(axis.title.y = element_blank(), 
                                  axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    coord_polar(theta = "y")
  #write.table((DATA), file = "clust_cond_freq_info.txt", sep = "\t", 
  #            row.names = F)
  #print("clust_cond_freq_info.txt file has beed generated.")
  if (plot.type == "bar") {
    return(myBP)
  }
  if (plot.type == "pie") {
    return(myPIE)
  }
}


gene.stats<-function (x = NULL, which.data = "raw.data", each.cond = F) 
{
  if ("Seurat" != class(x)[1]) {
    stop("x should be an object of class iCellR")
  }
  if (which.data == "raw") {
    DATA <- x@assays$RNA@counts
  }
  if (which.data == "data") {
    DATA <- x@assays$RNA@data
  }
  if (each.cond == T) {
    do<-x@meta.data$orig.ident
    if(split){
     do<-unlist(lapply(conds,function(s){return(str_split(s,"-")[[1]][1])}))
     }
    #do <- data.frame(do.call("rbind", strsplit(as.character(colnames(DATA)), 
    #                                         "_", fixed = TRUE)))[1]
    Myconds <- as.character(as.matrix(unique(do)))
    if (length(Myconds) > 1) {
      for (i in Myconds) {
        #ha <- subset(do, do == i)
        #dim2 <- max(as.numeric(rownames(ha)))
        #dim1 <- min(as.numeric(rownames(ha)))
        #myDATA <- DATA[dim1:dim2]
        myDATA<-DATA[,do==i]
        #mymat = as.matrix(myDATA)
        SDs <- apply(myDATA, 1, function(mymat) {
          sd(mymat)
        })
        
        Means<- apply(myDATA, 1, function(mymat) {
          mean(mymat)
        })
        
        nozero.nums<-apply(myDATA,1,function(mymat){
          num<-sum(as.numeric(mymat)!=0)
          return(num)
        })
        Table <- list(row.names(myDATA),nozero.nums,dim(myDATA)[2], 
                      nozero.nums/dim(myDATA)[2]* 100, as.numeric(Means), 
                      as.numeric(SDs))
        names(Table) <- c("genes", "numberOfCells", "totalNumberOfCells", 
                          "percentOfCells", "meanExp", "SDs")
        Table <- as.data.frame(Table)
        Table <- cbind(Table, condition = i)
        NameCol = paste("MyGeneStat", i, sep = "_")
        eval(call("<-", as.name(NameCol), Table))
      }
    }
  }
  #mymat = as.matrix(DATA)
  SDs <- apply(DATA, 1, function(mymat) {
    sd(mymat)
  })
  
  
  Means<- apply(DATA, 1, function(mymat) {
    mean(mymat)
  })
  
  nozero.nums<-apply(DATA,1,function(mymat){
    num<-sum(as.numeric(mymat)!=0)
    return(num)
  })
  
  Table <- list(row.names(DATA),nozero.nums,dim(DATA)[2], 
                  nozero.nums/dim(DATA)[2]* 100, as.numeric(Means), 
                  as.numeric(SDs))
  
  names(Table) <- c("genes", "numberOfCells", "totalNumberOfCells", 
                    "percentOfCells", "meanExp", "SDs")
  Table <- as.data.frame(Table)
  MyGeneStat_all <- cbind(Table, condition = "all")
  #filenames <- ls(pattern = "MyGeneStat_")
  #datalist <- mget(filenames)
  #Table <- do.call(rbind.data.frame, datalist)
  #row.names(Table) <- NULL
  #attributes(x)$gene.data <- Table
  return(Table)
}

Cov.Heatmap<-function(x=NULL,genes=NULL,
                      idents=0,slot="data",type="lower",
                      label=FALSE,p.mat=FALSE){
  if("Seurat"!=class(x)){
    stop("x must be a Seurat object!")
  }
  cells<-WhichCells(x,idents =idents)
  if (slot == "raw") {
    DATA <- x@assays$RNA@counts
  }
  if (slot== "data") {
    DATA <- x@assays$RNA@data
  }
  if(is.null(genes)){
    genes<-rownames(x)[str_detect(rownames(x),"^CD+")]
  }
  DATA<-t(as.matrix(DATA[genes,cells]))
  cormat <- round(cor(DATA),2)

  #require(reshape2)
  #require(ggplot2)

  require(ggcorrplot)
  if(p.mat){
    p.mat <- cor_pmat(DATA) ## Barring the no significant coefficient
  }else{
    p.mat<-NULL
  }
  p<-ggcorrplot(cormat,method = "square",
                type = type,
                outline.col = "white",
                hc.order = TRUE,
                show.diag = TRUE,
                lab = label,
                p.mat = p.mat,insig = "pch")+theme_bw()+
    theme(legend.text    = element_text(size = 12,face = "bold"),
          legend.title = element_text(size=12,face = "bold"),
          axis.text.y = element_text(size=12,face = "bold"),
          axis.text.x =element_text(size = 12,face = "bold",angle = 90),
          axis.title = element_blank(),
          panel.grid = element_blank(),panel.border = element_blank())
  print(p)
}

Seurat2Monocle<-function(seurat){
  #Load Seurat object
  require(monocle)
  if(class(seurat)[[1]]=="Seurat"){
    seurat_object<-seurat
  }else{
    seurat_object <- readRDS(seurat)
  }
  
  #Extract data, phenotype data, and feature data from the SeuratObject
  #data <- as(as.matrix(seurat_object@assays$RNA@data), 'sparseMatrix')
  data<-as.sparse(seurat_object@assays$RNA@data)
  pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)
  
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  
  #Construct monocle cds
  monocle_cds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.1,
                                expressionFamily = negbinomial.size())
  #print("save monocle object")
  #saveRDS(monocle_cds,file=monocle.path)
  print("done")
  return(monocle_cds)
}


violin.plot<-function(x=NULL,features=NULL,
                      group.by=NULL,
                      split=TRUE,
                      idents=NULL,
                      log=FALSE){
  if(!is.null(group.by)){
    if(group.by!="orig.ident"){
      stop("group.by is not accurate")
    }
    conds<-as.character(x@meta.data[,group.by])
    if(split){
      conds<-unlist(lapply(conds,function(c){
        return(str_split(c,"-")[[1]][1])
      }))
      x@meta.data$new.ident<-conds
      p<-VlnPlot(x,features = features,
                 group.by ="new.ident",
                 idents = idents,
                 log = log)
    }else{
      p<-VlnPlot(x,features = features,
                 group.by =group.by,
                 idents = idents,
                 log = log)
    }
  }else{
    p<-VlnPlot(x,features = features,idents = idents,log = log)
  }

  p<-p+stat_summary(fun.y=mean, geom="point", shape=23, size=5)
  p<-p+theme(legend.position = "none",
             axis.title = element_text(size = 18,face = "bold"),
             plot.title = element_text(size=20,face = "bold"))+
    xlab(NULL)
  return(p)
}

Cluster.VariableFeaturePlot<-function(object,ident=NULL,points=NULL,pt.size=1,
                                      cols=c("black", "red")){

  if(is.nul(ident)){
    obj<-object
  }else{
    obj<-subset(object,idents=ident)
  }

  p<-VariableFeaturePlot(obj,cols = cols,pt.size = pt.size)
  p<-LabelPoints(plot = p,points = points)
  p<-p+ggtitle(paste0("Cluster : ",ident))
  return(p)
}



Ident.counts.plot<-function(x,outdir,feature="disease"){
  features<-x@meta.data[[feature]]
  mat<-as.data.frame(as.matrix(table(Idents(x),features)))
  colnames(mat)<-c("cell","feature","Freq")

  p<-ggplot(data=mat,aes(feature))+
    geom_bar(aes(weight=Freq,fill=feature))+
    facet_wrap(~cell)+theme(legend.title = element_blank(),
                            legend.text = element_text(size = 12,face = "bold"),
                            axis.title.x = element_blank(),
                            axis.title.y = element_text(size=15,face="bold"),
                            axis.text.x = element_text(size = 12,face = "bold",angle = 90),
                            panel.grid.major = element_blank(),
                            strip.text = element_text(size=15,face = "bold")
                              )
  jpeg(paste0(outdir,"/","Ident.count.by","_",feature,".jpeg"))
  print(p)
  dev.off()
}


Cluster.FeaturePlot<-function(x,cluster,marker="CD4",outdir="."){
  x<-subset(x,idents = cluster)
  plot_list<-list()
  feature<-unique(x@meta.data[["disease"]])
  for(i in 1:length(feature)){
      cells<-Cells(x)[x@meta.data[["disease"]]==feature[i]]
      p<-FeaturePlot(x,features = marker,
                     cells = cells,
                     reduction = "tsne",
                     pt.size = 1.0)+ggtitle(label="")
        #ggtitle(paste0(feature[[i]],"(",cluster," : ",marker,")"))
      plot_list[[i]]<-p
  }
  p1<-CombinePlots(plot_list,ncol = 2)
  #jpeg(paste0(outdir,"/",cluster,"_",marker,"_","scatter.jpeg"))
  #print(p1)
  #dev.off()

  expr.value<-GetAssayData(x,slot = "data")[marker,]
  expr.value<-data.frame("expr"=expr.value,
                         "group"=x@meta.data[["disease"]])

  comparisions<-list(c("CAT","UBQ"),
                     c("CAT","UBR"),
                     c("CAT","UVQ"),
                     c("CAT","UVR"))
  p2<-ggplot(data=expr.value,aes(x=group,y=expr,fill=group))+
    geom_boxplot()+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),
                       method = "t.test",
                       comparisons = comparisions)+
    theme_bw()+NoGrid()+theme(axis.title.x = element_blank(),
                              axis.text = element_text(size=10,face = "bold"),
                              axis.text.y = element_text(size=12,face="bold"),
                              axis.title.y = element_text(size=12,face="bold"),
                              legend.title = element_blank(),
                              legend.text = element_text(size=12,face = "bold"))+
    ggtitle(paste(cluster," (",marker,") "))

  jpeg(file=paste0(outdir,"/",cluster,"_",marker,".jpeg"),height=720,width=720)
  print(CombinePlots(list(p1,p2),ncol = 2))
  dev.off()
}
