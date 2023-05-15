#!/usr/bin/env Rscript
if(F){
    txt <- "yachi_marker.txt"
}

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    stop("Usage:bin/marker.R yachi_marker.txt", call.=FALSE)
}
txt <- args[1]
message(txt)


library(tidyverse)

lines <- readLines(txt)

cell <- sapply(lines, function(x){
    items <- str_split(x,",")[[1]]
    items[1]
})
marker <- lapply(lines, function(x){
    items <- str_split(x,",")[[1]]
    items[2:length(items)]
})

names(marker) <- cell
saveRDS(marker, file="marker.rds")
