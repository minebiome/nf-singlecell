#!/usr/bin/env Rscript

library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)


if(F){
    seuratObjPath <- "results/seurat/pmbc_1683691139445_seurat_filter_seurat.rds"
    markerPath <- "results/marker/marker.rds"
    name <- str_replace(basename(seuratObjPath),".rds","")
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Usage: bin/featurePlot.R  results/seurat/pmbc_1683691139445_seurat_filter_seurat.rds results/marker/marker.rds", call.=FALSE)
}
seuratObjPath <- args[1]
markerPath <- args[2]
name <- str_replace(basename(seuratObjPath),".rds","")



seuratObj <- readRDS(seuratObjPath)
markers <- readRDS(markerPath)




lapply(names(markers), function(marker){
    cellType <- make.names(marker)
    len <- length(marker)
    fp <- FeaturePlot(seuratObj,features=markers[[marker]],ncol = 3,cols = c("lightgrey", "red"))
    if(len<=3){
        ggsave(filename=paste0(name,"/",cellType,"_featurePlot.png"),plot=fp,width = 300,height =100,units = "mm")
        ggsave(filename=paste0(name,"/",cellType,"_featurePlot.pdf"),plot=fp,width = 300,height =100,units = "mm")
    }
    
    if(len <=6 && len > 3){
        ggsave(filename=paste0(name,"/",cellType,"_featurePlot.png"),plot=fp,width = 300,height =200,units = "mm")
        ggsave(filename=paste0(name,"/",cellType,"_featurePlot.pdf"),plot=fp,width = 300,height =200,units = "mm")
    }
    
    if(len<=9 &&len>6){
       ggsave(filename=paste0(name,"/",cellType,"_featurePlot.png"),plot=fp,width = 300,height =300,units = "mm")
       ggsave(filename=paste0(name,"/",cellType,"_featurePlot.pdf"),plot=fp,width = 300,height =300,units = "mm")
    }
    
    if(len<=12&&len>9){
        ggsave(filename=paste0(name,"/",cellType,"_featurePlot.png"),plot=fp,width = 300,height =400,units = "mm")
        ggsave(filename=paste0(name,"/",cellType,"_featurePlot.pdf"),plot=fp,width = 300,height =400,units = "mm")
    }
    
})


