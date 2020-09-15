## ----include=FALSE------------------------------------------------------------
# Set path to plotly screenshot. We don't run the plotly code chunk as most servers do not have javascript libraries needed for interactive plotting
screenshot <- "../man/figures/plotly.png"

# The chunk below uses Rmd in man/fragments to avoid duplication, as the content is shared with the vignette and README. As suggested here: https://www.garrickadenbuie.com/blog/dry-vignette-and-readme/

## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(knitr)
knitr::opts_chunk$set(warning = FALSE, message = FALSE)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("tidyseurat")

## ---- eval=FALSE--------------------------------------------------------------
#  devtools::install_github("stemangiola/tidyseurat")

## -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(Seurat)
library(tidyseurat)

## -----------------------------------------------------------------------------
pbmc_small_tidy <- tidyseurat::pbmc_small %>% tidy()

## -----------------------------------------------------------------------------
pbmc_small_tidy

## -----------------------------------------------------------------------------
pbmc_small_tidy@assays

## -----------------------------------------------------------------------------
pbmc_small_tidy$file[1:5]

## -----------------------------------------------------------------------------
# Create sample column
pbmc_small_polished <-
  pbmc_small_tidy %>%
  extract(file, "sample", "../data/([a-z0-9]+)/outs.+", remove = FALSE)
# Reorder to have sample column up front
pbmc_small_polished %>%
  select(sample, everything())

## -----------------------------------------------------------------------------
# Use colourblind-friendly colours
if (requireNamespace("dittoSeq", quietly = TRUE)) {
      friendly_cols <- dittoSeq::dittoColors()
   } else {
      friendly_cols <- c("red", "blue", "green", "purple")
   }

# Set theme
my_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )

## ----plot1--------------------------------------------------------------------
pbmc_small_polished %>%
  tidyseurat::ggplot(aes(nFeature_RNA, fill = groups)) +
  geom_histogram() +
  my_theme

## ----plot2--------------------------------------------------------------------
pbmc_small_polished %>%
  tidyseurat::ggplot(aes(groups, nCount_RNA, fill = groups)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  my_theme

## -----------------------------------------------------------------------------
pbmc_small_polished %>%
  join_transcripts(transcripts = c("HLA-DRA", "LYZ")) %>%
  ggplot(aes(groups, abundance_RNA + 1, fill = groups)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(size = nCount_RNA), alpha = 0.5, width = 0.2) +
  scale_y_log10() +
  my_theme

## ----preprocess---------------------------------------------------------------
pbmc_small_pca <-
  pbmc_small_polished %>%
  SCTransform(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

pbmc_small_pca

## ----pc_plot------------------------------------------------------------------
pbmc_small_pca %>%
  as_tibble() %>%
  select(contains("PC"), everything()) %>%
  GGally::ggpairs(columns = 1:5, ggplot2::aes(colour = groups)) +
  my_theme

## ----cluster------------------------------------------------------------------
pbmc_small_cluster <-
  pbmc_small_pca %>%
  FindNeighbors(verbose = FALSE) %>%
  FindClusters(method = "igraph", verbose = FALSE)

pbmc_small_cluster

## ----cluster count------------------------------------------------------------
pbmc_small_cluster %>%
  tidyseurat::count(groups, seurat_clusters)

## -----------------------------------------------------------------------------
# Identify top 10 markers per cluster
markers <-
  pbmc_small_cluster %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_logFC)

# Plot heatmap
pbmc_small_cluster %>%
  DoHeatmap(
    features = markers$gene,
    group.colors = friendly_cols
  )

## ----umap---------------------------------------------------------------------
pbmc_small_UMAP <-
  pbmc_small_cluster %>%
  RunUMAP(reduction = "pca", dims = 1:15, n.components = 3L, )

## ----umap plot, eval=FALSE----------------------------------------------------
#  pbmc_small_UMAP %>%
#    plot_ly(
#      x = ~`UMAP_1`,
#      y = ~`UMAP_2`,
#      z = ~`UMAP_3`,
#      color = ~seurat_clusters,
#      colors = friendly_cols[1:4]
#    )

## ----eval=FALSE---------------------------------------------------------------
#  # Get cell type reference data
#  blueprint <- celldex::BlueprintEncodeData()
#  
#  # Infer cell identities
#  cell_type_df <-
#    pbmc_small_UMAP@assays[["SCT"]]@counts %>%
#    log1p() %>%
#    Matrix::Matrix(sparse = TRUE) %>%
#    SingleR::SingleR(
#      ref = blueprint,
#      labels = blueprint$label.main,
#      method = "single"
#    ) %>%
#    as.data.frame() %>%
#    as_tibble(rownames = "cell") %>%
#    select(cell, first.labels)

## -----------------------------------------------------------------------------
# Join UMAP and cell type info
pbmc_small_cell_type <-
  pbmc_small_UMAP %>%
  left_join(cell_type_df, by = "cell")

# Reorder columns
pbmc_small_cell_type %>%
  tidyseurat::select(cell, first.labels, everything())

## -----------------------------------------------------------------------------
pbmc_small_cell_type %>%
  count(seurat_clusters, first.labels)

## -----------------------------------------------------------------------------
pbmc_small_cell_type %>%

  # Reshape and add classifier column
  pivot_longer(
    cols = c(seurat_clusters, first.labels),
    names_to = "classifier", values_to = "label"
  ) %>%

  # UMAP plots for cell type and cluster
  ggplot(aes(UMAP_1, UMAP_2, color = label)) +
  geom_point() +
  facet_wrap(~classifier) +
  my_theme

## -----------------------------------------------------------------------------
pbmc_small_cell_type %>%

  # Add some mitochondrial abundance values
  mutate(mitochondrial = rnorm(n())) %>%

  # Plot correlation
  join_transcripts(transcripts = c("CST3", "LYZ"), shape = "wide") %>%
  ggplot(aes(CST3 + 1, LYZ + 1, color = groups, size = mitochondrial)) +
  geom_point() +
  facet_wrap(~first.labels, scales = "free") +
  scale_x_log10() +
  scale_y_log10() +
  my_theme

## -----------------------------------------------------------------------------
pbmc_small_nested <-
  pbmc_small_cell_type %>%
  filter(first.labels != "Erythrocytes") %>%
  mutate(cell_class = if_else(`first.labels` %in% c("Macrophages", "Monocytes"), "myeloid", "lymphoid")) %>%
  nest(data = -cell_class)

pbmc_small_nested

## -----------------------------------------------------------------------------
pbmc_small_nested_reanalysed <-
  pbmc_small_nested %>%
  mutate(data = map(
    data, ~ .x %>%
      FindVariableFeatures(verbose = FALSE) %>%
      RunPCA(npcs = 10, verbose = FALSE) %>%
      FindNeighbors(verbose = FALSE) %>%
      FindClusters(method = "igraph", verbose = FALSE) %>%
      RunUMAP(reduction = "pca", dims = 1:10, n.components = 3L, verbose = FALSE)
  ))

pbmc_small_nested_reanalysed

## -----------------------------------------------------------------------------
pbmc_small_nested_reanalysed %>%

  # Convert to tibble otherwise Seurat drops reduced dimensions when unifying data sets.
  mutate(data = map(data, ~ .x %>% as_tibble())) %>%
  unnest(data) %>%

  # Define unique clusters
  unite("cluster", c(cell_class, seurat_clusters), remove = FALSE) %>%

  # Plotting
  ggplot(aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point() +
  facet_wrap(~cell_class) +
  my_theme

## ---- eval = FALSE------------------------------------------------------------
#  library(SingleCellSignalR)
#  
#  pbmc_small_nested_interactions <-
#    pbmc_small_nested_reanalysed %>%
#  
#    # Unnest based on cell category
#    unnest(data) %>%
#  
#    # Create unambiguous clusters
#    mutate(integrated_clusters = first.labels %>% as.factor() %>% as.integer()) %>%
#  
#    # Nest based on sample
#    tidyseurat::nest(data = -sample) %>%
#    tidyseurat::mutate(interactions = map(data, ~ {
#  
#      # Produce variables. Yuck!
#      cluster <- .x@meta.data$integrated_clusters
#      data <- data.frame(.x[["SCT"]]@data)
#  
#      # Ligand/Receptor analysis using SingleCellSignalR
#      data %>%
#        cell_signaling(genes = rownames(data), cluster = cluster) %>%
#        inter_network(data = data, signal = ., genes = rownames(data), cluster = cluster) %$%
#        `individual-networks` %>%
#        map_dfr(~ bind_rows(as_tibble(.x)))
#    }))
#  
#  pbmc_small_nested_interactions %>%
#    select(-data) %>%
#    unnest(interactions)

## ---- echo = FALSE------------------------------------------------------------
tidyseurat::pbmc_small_nested_interactions

## -----------------------------------------------------------------------------
sessionInfo()

