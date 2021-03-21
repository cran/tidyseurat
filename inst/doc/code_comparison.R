## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(tidyseurat)

tidy_object = pbmc_small %>% tidy()
seurat_object = pbmc_small


## -----------------------------------------------------------------------------
tidy_object %>%
  add_count(file) %>%
  filter(PC_1 > 0 & n > 40)

## -----------------------------------------------------------------------------
pca_emb = Embeddings(object = seurat_object, reduction = "pca")
cell_pca = rownames(pca_emb[pca_emb[,1]>0,])
n =
  seurat_object@meta.data %>%
  add_count(file) %>%
  pull(n)
seurat_object = AddMetaData( object = seurat_object, metadata = n, col.name = 'n')
seurat_object %>%
  subset(cells = cell_pca) %>%
  subset( n < 40)


