#!/usr/bin/env Rscript

reticulate::use_python("/ssd1/wy/.conda/envs/python3.9/bin")
library(cellassign)

# library(tensorflow)
# tensorflow::tf_config()

library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(scran)


if(F){

    rds <- "/data/wangyang/yachi_singlecell_GSE189381/work/08/d7fa4efcd60b9844c84e5df722ebcb/E16.5_seurat.rds"
    name <- str_replace(basename(rds),".rds","")
    marker <- "/data/wangyang/yachi_singlecell_GSE189381/work/08/d7fa4efcd60b9844c84e5df722ebcb/yachi.rds"
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Usage: bin/cellassign.R  results/seurat/pmbc_1683691139445_seurat.rds pbmc.rds pbmc.rds", call.=FALSE)
}
rds <- args[1]
marker <- args[2]
message(rds)
name <- str_replace(basename(rds),".rds","")






# data(example_sce)
# data(example_marker_mat)

# fit <- em_result <- cellassign(example_sce[rownames(example_marker_mat),],
#     marker_gene_info = example_marker_mat,
#     s = colSums(SummarizedExperiment::assay(example_sce, "counts")),
#     learning_rate = 1e-2,
#     shrinkage = TRUE,
#     verbose = FALSE)





marker_gene_list <- readRDS(marker)
marker_gene_mat <- marker_list_to_mat(marker_gene_list)
head(marker_gene_mat)


obj <- readRDS(paste0(rds))

mat <- obj@assays$RNA@counts
intersect_gene <-  intersect(rownames(mat), rownames(marker_gene_mat)) 
dim(marker_gene_mat)
setdiff(rownames(marker_gene_mat), rownames(mat))
marker_gene_mat <- marker_gene_mat[intersect_gene,]

# mat[rownames(marker_gene_mat),] |> as.matrix()

sceset <- SingleCellExperiment(assays=list(counts=mat))
qclust  <- quickCluster(sceset)
sceset <-  computeSumFactors(sceset,clusters=qclust)
s <- sizeFactors(sceset)


fit  <- cellassign(sceset[rownames(marker_gene_mat),],
                marker_gene_info = marker_gene_mat,
                s =s,
                learning_rate = 1e-2,
                shrinkage = TRUE,
                verbose = FALSE,threads=30)
saveRDS(fit,file=paste0(name,"_cellassign.rds"))