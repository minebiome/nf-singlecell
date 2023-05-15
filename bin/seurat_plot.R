#!/usr/bin/env Rscript
if(F){
    seuratObjPath <- "results/seurat/pmbc_1683691139445_seurat_filter_seurat.rds"
    cellassignObjPath <- "results/assign/pmbc_1683691139445_seurat_filter_cellassign.rds"
    name <- str_replace(basename(seuratObjPath),".rds","")
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Usage: bin/seurat_plot.R  pmbc_1683691139445.rds", call.=FALSE)
}
seuratObjPath <- args[1]
cellassignObjPath <- args[2]
name <- str_replace(basename(seuratObjPath),".rds","")

library(Seurat)
library(tidyverse)
library(patchwork)

reticulate::use_python("/ssd1/wy/.conda/envs/python3.9/bin")
library(cellassign)

library(SingleCellExperiment)
library(scran)


seuratObj <- readRDS(seuratObjPath)
cellassignObj <- readRDS(cellassignObjPath)

dp <- DimPlot(seuratObj, reduction = "umap",label=T)
ggsave(filename=paste0(name,"_umap.png"),plot=dp,width=150,units="mm",height=150)


seuratObj@meta.data$cellassign <- celltypes(cellassignObj)

dp <- DimPlot(seuratObj, reduction = "umap",group.by="cellassign",label=T)
ggsave(filename=paste0(name,"_cellassign_umap.png"),plot=dp,width=150,units="mm",height=150)


delta <- cellassignObj$mle_params$delta
feature <- rownames(delta)
dp <- DotPlot(seuratObj, features =feature)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(bg = "white",filename=paste0(name,"_cellassign_dotplot.png"),plot=dp,width=200,units="mm",height=150)

fp <- FeaturePlot(obj, features =feature, reduction="umap") + plot_annotation(title = "")
file = paste0("output3/E14.5","/umap_feature.png")
ggsave(filename=file,plot=fp,width=230,units="mm",height=150)