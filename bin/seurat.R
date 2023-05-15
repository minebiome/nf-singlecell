#!/usr/bin/env Rscript
if(F){
    rds <- "pmbc_1683691139445.rds"
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Usage: bin/seurat.R pmbc_1683691139445.rds", call.=FALSE)
}
rds <- args[1]
message(rds)


library(Seurat)
library(tidyverse)
library(patchwork)

name <- str_replace(basename(rds),".rds","")

obj <- readRDS(rds)
# obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 50)
# vp <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ggsave(filename=paste0(name,"_VlnPlot.png"),plot=vp,width=230,units="mm",height=150)

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)
obj <- RunPCA(obj)

ggsave(filename=paste0(name,"_eblow.png"), plot=ElbowPlot(obj,ndims=50))

obj <- RunUMAP(obj, reduction = "pca", dims = 1:20,min.dist=0.3)
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.4)

saveRDS(obj,file=paste0(name,"_seurat.rds"))

# ggsave(filename=paste0(name,"_umap.png"),plot=dp,width=150,units="mm",height=150)
# ggsave(filename=paste0(name,"_umap.pdf"),plot=dp,width=150,units="mm",height=150)

