# Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

set.seed(42)

# Set project paths
project_root <- Sys.getenv("SCRNA_PROJECT_ROOT", unset = getwd())
data_root <- file.path(project_root, "output", "data")
result_root <- file.path(project_root, "output", "result")

################ Figure 1b: Before QC ################
# Load pre-QC metadata
meta_before <- readRDS(
  file.path(data_root, "plot_data", "qc_meta_before.rds")
)

qc_before <- rbind(
  data.frame(
    sample = meta_before$sample,
    metric = "nFeature_RNA",
    value = meta_before$nFeature_RNA
  ),
  data.frame(
    sample = meta_before$sample,
    metric = "nCount_RNA",
    value = meta_before$nCount_RNA
  ),
  data.frame(
    sample = meta_before$sample,
    metric = "max_umi_ratio_clean",
    value = meta_before$max_umi_ratio_clean
  )
)

qc_before$metric <- factor(
  qc_before$metric,
  levels = c(
    "nFeature_RNA",
    "nCount_RNA",
    "max_umi_ratio_clean"
  )
)

qc_before$sample <- factor(
  qc_before$sample,
  levels = c("ner_1", "ner_2", "ten_1", "ten_2")
)

# Plot QC metric distributions
p1b <- ggplot(
  qc_before,
  aes(x = sample, y = value, fill = sample)
) +
  geom_violin(
    scale = "width",
    linewidth = 0.3
  ) +
  facet_wrap(
    ~metric,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "ner_1" = "#5b99c5",
      "ner_2" = "#faa256",
      "ten_1" = "#badaba",
      "ten_2" = "#f0afaf"
    )
  ) +
  scale_y_continuous(
    labels = scales::label_comma(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Before QC",
    x = "Sample",
    y = NULL
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    strip.background = element_blank(),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    axis.ticks.length = grid::unit(0.08, "cm"),
    legend.position = "none"
  )

p1b

rm(meta_before, qc_before, p1b)

gc()

################ Figure 1c: After QC ################
# Load QC-filtered Seurat object
seob_qc <- readRDS(
  file.path(data_root, "seurat", "seob_qc.rds")
)

meta_after <- seob_qc@meta.data

qc_after <- rbind(
  data.frame(
    sample = meta_after$sample,
    metric = "nFeature_RNA",
    value = meta_after$nFeature_RNA
  ),
  data.frame(
    sample = meta_after$sample,
    metric = "nCount_RNA",
    value = meta_after$nCount_RNA
  ),
  data.frame(
    sample = meta_after$sample,
    metric = "max_umi_ratio_clean",
    value = meta_after$max_umi_ratio_clean
  )
)

qc_after$metric <- factor(
  qc_after$metric,
  levels = c(
    "nFeature_RNA",
    "nCount_RNA",
    "max_umi_ratio_clean"
  )
)

qc_after$sample <- factor(
  qc_after$sample,
  levels = c("ner_1", "ner_2", "ten_1", "ten_2")
)

# Plot QC metric distributions after QC
p1c <- ggplot(
  qc_after,
  aes(x = sample, y = value, fill = sample)
) +
  geom_violin(
    scale = "width",
    linewidth = 0.3
  ) +
  facet_wrap(
    ~metric,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "ner_1" = "#5b99c5",
      "ner_2" = "#faa256",
      "ten_1" = "#badaba",
      "ten_2" = "#f0afaf"
    )
  ) +
  scale_y_continuous(
    labels = scales::label_comma(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "After QC",
    x = "Sample",
    y = NULL
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    strip.background = element_blank(),
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    axis.ticks.length = grid::unit(0.08, "cm"),
    legend.position = "none"
  )

p1c

rm(seob_qc, meta_after, qc_after, p1c)

gc()

################ Figure 1d: Cell counts before and after QC ################
# Load QC summary
qc_summary <- read.csv(
  file.path(data_root, "seurat", "qc_summary.csv")
)

qc_counts <- qc_summary %>%
  dplyr::select(sample, before, after) %>%
  tidyr::pivot_longer(
    cols = c(before, after),
    names_to = "stage",
    values_to = "cells"
  )

qc_counts$sample <- factor(
  qc_counts$sample,
  levels = c("ner_1", "ner_2", "ten_1", "ten_2")
)

qc_counts$stage <- factor(
  qc_counts$stage,
  levels = c("before", "after"),
  labels = c("Before QC", "After QC")
)

# Plot cell counts before and after QC
p1d <- ggplot(
  qc_counts,
  aes(x = stage, y = cells, fill = sample)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "grey30",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = scales::comma(cells)),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 3.5,
    family = "Arial"
  ) +
  scale_fill_manual(
    values = c(
      "ner_1" = "#5b99c5",
      "ner_2" = "#faa256",
      "ten_1" = "#badaba",
      "ten_2" = "#f0afaf"
    )
  ) +
  scale_y_continuous(
    labels = scales::label_comma(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Cell counts before and after QC",
    x = NULL,
    y = "Cell count",
    fill = "Sample"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 10,
      face = "bold"
    ),
    axis.ticks.length = grid::unit(0.08, "cm"),
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    legend.text = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "right"
  )

p1d

rm(qc_summary, qc_counts, p1d)

gc()

################ Figure 1e: Singlet and doublet counts ################
# Load doublet detection results
doublet_results <- readRDS(
  file.path(data_root, "seurat", "doublet_results.rds")
)

doublet_counts <- as.data.frame(
  table(
    sample = factor(
      doublet_results$sample,
      levels = c("ner_1", "ner_2", "ten_1", "ten_2")
    ),
    classification = factor(
      doublet_results$class,
      levels = c("singlet", "doublet")
    )
  )
)

names(doublet_counts)[3] <- "cells"

doublet_counts$classification <- factor(
  doublet_counts$classification,
  levels = c("singlet", "doublet"),
  labels = c("Singlet", "Doublet")
)

cell_totals <- aggregate(
  cells ~ sample,
  data = doublet_counts,
  FUN = sum
)

doublet_labels <- subset(
  doublet_counts,
  classification == "Doublet"
)

doublet_labels$total <- cell_totals$cells[
  match(doublet_labels$sample, cell_totals$sample)
]

doublet_labels$label <- paste0(
  scales::comma(doublet_labels$cells),
  " (",
  sprintf(
    "%.1f",
    doublet_labels$cells / doublet_labels$total * 100
  ),
  "%)"
)

print(
  xtabs(cells ~ sample + classification, data = doublet_counts)
)

# Plot singlet and doublet counts
p1e <- ggplot(
  doublet_counts,
  aes(x = sample, y = cells, fill = classification)
) +
  geom_col(
    position = position_stack(reverse = TRUE),
    width = 0.7,
    color = "grey30",
    linewidth = 0.3
  ) +
  geom_text(
    data = subset(doublet_counts, classification == "Singlet"),
    aes(label = scales::comma(cells)),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    size = 3.5,
    family = "Arial"
  ) +
  geom_text(
    data = doublet_labels,
    aes(x = sample, y = total, label = label),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 3.5,
    family = "Arial",
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Singlet" = "#3C77AF",
      "Doublet" = "#D1352B"
    )
  ) +
  scale_y_continuous(
    labels = scales::label_comma(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Singlet and doublet counts",
    x = "Sample",
    y = "Cell count",
    fill = "Classification"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      face = "bold"
    ),
    axis.ticks.length = grid::unit(0.08, "cm"),
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    legend.text = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "right"
  )

p1e

rm(
  doublet_results, doublet_counts,
  cell_totals, doublet_labels, p1e
)

gc()

################ Figure 2a: Before Harmony ################
# Load the clustered Seurat object
seob <- readRDS(
  file.path(data_root, "seurat", "seob_clustered.rds")
)

# Set the pre-Harmony UMAP cache path
umap_pca_file <- file.path(data_root, "plot_data", "umap_pca.rds")

if (file.exists(umap_pca_file)) {
  seob[["umap.pca"]] <- readRDS(umap_pca_file)
} else {
  seob <- RunUMAP(
    seob,
    reduction = "pca",
    dims = 1:50,
    reduction.name = "umap.pca"
  )
  
  saveRDS(seob[["umap.pca"]], umap_pca_file)
}

umap_df <- Embeddings(seob, "umap.pca")

x_range <- range(umap_df[, 1])
y_range <- range(umap_df[, 2])

x0 <- x_range[1] - diff(x_range) * 0.03
y0 <- y_range[1] - diff(y_range) * 0.03

x_len <- diff(x_range) * 0.06
y_len <- diff(y_range) * 0.06

# Plot UMAP before Harmony integration
p2a <- DimPlot(
  seob,
  reduction = "umap.pca",
  group.by = "sample",
  pt.size = 0.02,
  shuffle = TRUE,
  cols = c(
    "ner_1" = "#5b99c5",
    "ner_2" = "#faa256",
    "ten_1" = "#badaba",
    "ten_2" = "#f0afaf"
  )
) +
  scale_x_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  coord_equal(clip = "off") +
  labs(
    title = "Before Harmony",
    x = NULL,
    y = NULL,
    color = "Sample"
  ) +
  
  annotate(
    "segment",
    x = x0, y = y0,
    xend = x0 + x_len, yend = y0,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = x0, y = y0,
    xend = x0, yend = y0 + y_len,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = x0 + x_len * 0.55,
    y = y0 - y_len * 0.50,
    label = "UMAP 1",
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +
  annotate(
    "text",
    x = x0 - x_len * 0.50,
    y = y0 + y_len * 0.55,
    label = "UMAP 2",
    angle = 90,
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +
  
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    legend.text = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "right",
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

p2a

rm(
  seob, umap_pca_file, umap_df,
  x_range, y_range,
  x0, y0, x_len, y_len, p2a
)

gc()

################ Figure 2b: After Harmony ################
# Load the clustered Seurat object
seob <- readRDS(
  file.path(data_root, "seurat", "seob_clustered.rds")
)

umap_df <- Embeddings(seob, "umap.harmony")

x_range <- range(umap_df[, 1])
y_range <- range(umap_df[, 2])

x0 <- x_range[1] - diff(x_range) * 0.03
y0 <- y_range[1] - diff(y_range) * 0.03

x_len <- diff(x_range) * 0.06
y_len <- diff(y_range) * 0.06

# Plot UMAP after Harmony integration
p2b <- DimPlot(
  seob,
  reduction = "umap.harmony",
  group.by = "sample",
  pt.size = 0.02,
  shuffle = TRUE,
  cols = c(
    "ner_1" = "#5b99c5",
    "ner_2" = "#faa256",
    "ten_1" = "#badaba",
    "ten_2" = "#f0afaf"
  )
) +
  scale_x_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  coord_equal(clip = "off") +
  labs(
    title = "After Harmony",
    x = NULL,
    y = NULL,
    color = "Sample"
  ) +

  annotate(
    "segment",
    x = x0, y = y0,
    xend = x0 + x_len, yend = y0,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = x0, y = y0,
    xend = x0, yend = y0 + y_len,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +

  annotate(
    "text",
    x = x0 + x_len * 0.55,
    y = y0 - y_len * 0.50,
    label = "UMAP 1",
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +
  annotate(
    "text",
    x = x0 - x_len * 0.50,
    y = y0 + y_len * 0.55,
    label = "UMAP 2",
    angle = 90,
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +
  
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    legend.text = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "right",
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

p2b

rm(
  seob, umap_df,
  x_range, y_range,
  x0, y0, x_len, y_len, p2b
)

gc()

################ Figure 2c: Cell clusters ################
# Load the clustered Seurat object
seob <- readRDS(
  file.path(data_root, "seurat", "seob_clustered.rds")
)

umap_df <- data.frame(
  UMAP_1 = Embeddings(seob, "umap.harmony")[, 1],
  UMAP_2 = Embeddings(seob, "umap.harmony")[, 2],
  cluster = factor(seob$harmony_clusters)
)

label_df <- aggregate(
  cbind(UMAP_1, UMAP_2) ~ cluster,
  data = umap_df,
  FUN = median
)

base_colors <- c(
  "#D2EBC8", "#3C77AF",
  "#7DBFA7", "#AECDE1",
  "#EE934E", "#D1352B",
  "#9B5B33", "#F5CFE4",
  "#B383B9", "#8FA4AE",
  "#FCED82", "#F5D2A8",
  "#BBDD78"
)

cluster_colors <- grDevices::colorRampPalette(base_colors)(
  nlevels(umap_df$cluster)
)

x_range <- range(umap_df$UMAP_1)
y_range <- range(umap_df$UMAP_2)

x0 <- x_range[1] - diff(x_range) * 0.03
y0 <- y_range[1] - diff(y_range) * 0.03

x_len <- diff(x_range) * 0.06
y_len <- diff(y_range) * 0.06

# Plot UMAP of 33 clusters
p2c <- ggplot(
  umap_df,
  aes(x = UMAP_1, y = UMAP_2, color = cluster)
) +
  geom_point(
    size = 0.12,
    alpha = 0.9
  ) +
  geom_text(
    data = label_df,
    aes(label = cluster),
    color = "black",
    size = 3.2,
    fontface = "bold",
    family = "Arial"
  ) +
  scale_color_manual(
    values = cluster_colors
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  coord_equal(clip = "off") +
  labs(
    title = "Cell clusters",
    x = NULL,
    y = NULL
  ) +
  
  annotate(
    "segment",
    x = x0, y = y0,
    xend = x0 + x_len, yend = y0,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = x0, y = y0,
    xend = x0, yend = y0 + y_len,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +

  annotate(
    "text",
    x = x0 + x_len * 0.55,
    y = y0 - y_len * 0.50,
    label = "UMAP 1",
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +
  annotate(
    "text",
    x = x0 - x_len * 0.50,
    y = y0 + y_len * 0.55,
    label = "UMAP 2",
    angle = 90,
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    legend.position = "none",
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

p2c

rm(
  seob, umap_df, label_df,
  base_colors, cluster_colors,
  x_range, y_range,
  x0, y0, x_len, y_len, p2c
)

gc()

################ Figure 3a: Cell types ################
# Load the annotated Seurat object
seob <- readRDS(
  file.path(data_root, "seurat", "seob_celltype.rds")
)

# Define cell type colors
cell_type_colors <- c(
  "Neurons" = "#3C77AF",
  "Mesenchymal cells" = "#7DBFA7",
  "Epithelial cells" = "#EE934E",
  "Muscle cells" = "#D1352B",
  "Immune cells" = "#9B5B33",
  "Ciliated cells" = "#B383B9",
  "Proliferating cells" = "#BBDD78",
  "Uncertain" = "#8FA4AE"
)

umap_df <- Embeddings(seob, "umap.harmony")

x_range <- range(umap_df[, 1])
y_range <- range(umap_df[, 2])

x0 <- x_range[1] - diff(x_range) * 0.03
y0 <- y_range[1] - diff(y_range) * 0.03

x_len <- diff(x_range) * 0.06
y_len <- diff(y_range) * 0.06

# Plot UMAP by cell type
p3a <- DimPlot(
  seob,
  reduction = "umap.harmony",
  group.by = "cell_type",
  pt.size = 0.12,
  shuffle = TRUE,
  cols = cell_type_colors
) +
  scale_x_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  coord_equal(clip = "off") +
  labs(
    title = "Cell types",
    x = NULL,
    y = NULL,
    color = "Cell type"
  ) +

  annotate(
    "segment",
    x = x0, y = y0,
    xend = x0 + x_len, yend = y0,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +
  annotate(
    "segment",
    x = x0, y = y0,
    xend = x0, yend = y0 + y_len,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +

  annotate(
    "text",
    x = x0 + x_len * 0.55,
    y = y0 - y_len * 0.50,
    label = "UMAP 1",
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +
  annotate(
    "text",
    x = x0 - x_len * 0.50,
    y = y0 + y_len * 0.55,
    label = "UMAP 2",
    angle = 90,
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    legend.text = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "right",
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

p3a

rm(
  seob, cell_type_colors, umap_df,
  x_range, y_range,
  x0, y0, x_len, y_len, p3a
)

gc()

################ Figure 3b: Cell type marker heatmap ################
# Load the annotated Seurat object
seob <- readRDS(
  file.path(data_root, "seurat", "seob_celltype.rds")
)

DefaultAssay(seob) <- "RNA"

# Load the top 30 cluster markers
top30_markers <- read.csv(
  file.path(result_root, "top30_markers.csv")
)

cluster_annotation <- unique(
  seob@meta.data[, c("harmony_clusters", "cell_type")]
)

cluster_to_type <- setNames(
  as.character(cluster_annotation$cell_type),
  as.character(cluster_annotation$harmony_clusters)
)

cell_types <- levels(seob$cell_type)

cell_type_colors <- c(
  "Neurons" = "#3C77AF",
  "Mesenchymal cells" = "#7DBFA7",
  "Epithelial cells" = "#EE934E",
  "Muscle cells" = "#D1352B",
  "Immune cells" = "#9B5B33",
  "Ciliated cells" = "#B383B9",
  "Proliferating cells" = "#BBDD78",
  "Uncertain" = "#8FA4AE"
)

# Select 20 markers per cell type
n_genes_type <- 20

type_marker_pool <- top30_markers %>%
  mutate(
    cell_type = unname(
      cluster_to_type[as.character(cluster)]
    )
  ) %>%
  filter(cell_type %in% cell_types) %>%
  group_by(cell_type, gene) %>%
  summarise(
    max_avg_log2FC = max(avg_log2FC),
    .groups = "drop"
  ) %>%
  group_by(cell_type) %>%
  slice_max(
    order_by = max_avg_log2FC,
    n = n_genes_type,
    with_ties = FALSE
  ) %>%
  ungroup()

genes_info <- type_marker_pool %>%
  arrange(
    factor(cell_type, levels = cell_types),
    desc(max_avg_log2FC)
  ) %>%
  mutate(
    row_id = paste(cell_type, gene, sep = "__")
  )

genes <- genes_info$gene

# Save selected marker genes
write.csv(
  genes_info,
  file.path(result_root, "heatmap_cell_type_top20_genes.csv"),
  row.names = FALSE
)

# Calculate average expression by cell type
avg_expr <- AverageExpression(
  seob,
  assays = "RNA",
  group.by = "cell_type",
  layer = "data"
)$RNA

avg_expr <- avg_expr[
  genes,
  cell_types,
  drop = FALSE
]

rownames(avg_expr) <- genes_info$row_id

# Calculate row Z-scores
mat <- t(scale(t(as.matrix(avg_expr))))

mat[is.na(mat)] <- 0

heatmap_colors <- colorRampPalette(
  c("#D9E6F2", "#FFFFFF", "#C186A8")
)(100)

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  "Cell type" = factor(
    cell_types,
    levels = cell_types
  ),
  col = list(
    "Cell type" = cell_type_colors
  ),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    "Cell type" = list(
      title = "Cell type",
      title_gp = grid::gpar(
        fontsize = 11,
        fontface = "bold",
        fontfamily = "Arial"
      ),
      labels_gp = grid::gpar(
        fontsize = 10,
        fontfamily = "Arial"
      )
    )
  )
)

# Construct the heatmap
ht <- ComplexHeatmap::Heatmap(
  mat,
  name = "Row Z-score",
  col = heatmap_colors,
  top_annotation = top_annotation,

  column_split = factor(
    cell_types,
    levels = cell_types
  ),

  column_gap = grid::unit(0.4, "mm"),
  border = TRUE,
  border_gp = grid::gpar(
    col = "#CCCCCC",
    lwd = 0.35
  ),

  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,

  show_row_names = FALSE,
  show_column_names = FALSE,

  row_title = NULL,
  column_title = NULL,

  rect_gp = grid::gpar(col = NA),

  use_raster = FALSE,

  heatmap_legend_param = list(
    title = "Row Z-score",
    at = c(-2, -1, 0, 1, 2),
    title_gp = grid::gpar(
      fontsize = 11,
      fontface = "bold",
      fontfamily = "Arial"
    ),
    labels_gp = grid::gpar(
      fontsize = 10,
      fontfamily = "Arial"
    )
  )
)

p3b <- grid::grid.grabExpr(
  {
    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(
          nrow = 2,
          ncol = 1,
          heights = grid::unit.c(
            grid::unit(10, "mm"),
            grid::unit(1, "null")
          )
        )
      )
    )

    grid::grid.text(
      "Cell type marker expression",
      x = 0.5,
      y = 0.5,
      vp = grid::viewport(
        layout.pos.row = 1
      ),
      gp = grid::gpar(
        fontsize = 14,
        fontface = "bold",
        fontfamily = "Arial"
      )
    )

    grid::pushViewport(
      grid::viewport(
        layout.pos.row = 2
      )
    )

    ComplexHeatmap::draw(
      ht,
      newpage = FALSE,
      merge_legends = TRUE,
      heatmap_legend_side = "right",
      annotation_legend_side = "right"
    )

    grid::popViewport(2)
  },
  width = 230 / 25.4,
  height = 165 / 25.4
)

grid::grid.newpage()
grid::grid.draw(p3b)

rm(
  seob, top30_markers,
  cluster_annotation, cluster_to_type, cell_types,
  cell_type_colors, n_genes_type, type_marker_pool,
  genes_info, genes, avg_expr, mat,
  heatmap_colors, top_annotation, ht, p3b
)

gc()

################ Figure 3c: Cell type annotation markers ################
# Load the annotated Seurat object
seob <- readRDS(
  file.path(data_root, "seurat", "seob_celltype.rds")
)

DefaultAssay(seob) <- "RNA"

cell_types <- levels(seob$cell_type)

# Define 21 manually reviewed marker genes
marker_genes <- data.frame(
  cell_type = rep(
    cell_types[1:7],
    each = 3
  ),
  Gene_Symbol = c(
    "RASGRF1", "SYT2", "ELAVL2",
    "COL5A2", "TIMP3", "COL9A1",
    "MUC2", "KLF5", "MUC5AC",
    "MYH3", "KLHL41", "MYLK",
    "NLRC4", "MRC1", "CTSL",
    "ODF3L2", "PIFO", "TEKT4",
    "PRC1", "BUB1", "TPX2"
  ),
  gene = c(
    "evm.TU.Chr16.324",
    "evm.TU.Chr18.241",
    "evm.TU.Chr7.455",
    
    "evm.TU.Chr10.483",
    "evm.TU.Chr5.197",
    "evm.TU.Chr3.701",
    
    "evm.TU.Chr20.42",
    "evm.TU.Chr13.666",
    "evm.TU.Chr20.41",
    
    "evm.TU.Chr8.472",
    "evm.TU.Chr10.310",
    "evm.TU.Chr23.447",
    
    "evm.TU.Chr18.561",
    "evm.TU.Chr7.710",
    "evm.TU.Chr10.347",
    
    "evm.TU.Chr16.757",
    "evm.TU.Chr14.555",
    "evm.TU.Chr8.740",
    
    "evm.TU.Chr17.555",
    "evm.TU.Chr12.408",
    "evm.TU.Chr9.575"
  )
)

# Check marker gene IDs
missing_genes <- setdiff(
  marker_genes$gene,
  rownames(seob)
)

if (length(missing_genes) > 0) {
  stop(
    "The following marker gene IDs were not found in the Seurat object: ",
    paste(missing_genes, collapse = ", ")
  )
}

gene_labels <- setNames(
  marker_genes$Gene_Symbol,
  marker_genes$gene
)

feature_order <- marker_genes$gene

# Plot cell type annotation markers
p3c <- DotPlot(
  seob,
  features = feature_order,
  group.by = "cell_type",
  cols = c("#E8F1F6", "#3C77AF"),
  dot.scale = 4,
  scale.by = "size",
  col.min = -2,
  col.max = 2
) +
  scale_x_discrete(
    labels = gene_labels
  ) +
  scale_y_discrete(
    limits = rev(cell_types)
  ) +
  guides(
    size = guide_legend(
      order = 1,
      ncol = 1,
      title.position = "left",
      theme = theme(
        legend.title = element_text(
          size = 11,
          face = "bold",
          angle = 90,
          hjust = 0.5,
          vjust = 0.5
        )
      )
    ),
    color = guide_colorbar(
      order = 2,
      direction = "vertical",
      title.position = "left",
      barheight = grid::unit(30, "mm"),
      barwidth = grid::unit(4, "mm"),
      theme = theme(
        legend.title = element_text(
          size = 11,
          face = "bold",
          angle = 90,
          hjust = 0.5,
          vjust = 0.5
        )
      )
    )
  ) +
  labs(
    title = "Cell type annotation markers",
    x = NULL,
    y = NULL,
    color = "Scaled average expression",
    size = "Percent expressed"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",

    axis.text.x = element_text(
      size = 10,
      face = "bold.italic",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),

    axis.text.y = element_text(
      size = 10,
      face = "bold"
    ),

    axis.line = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(0.08, "cm"),

    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.spacing.y = grid::unit(4, "mm"),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

p3c

rm(
  seob, cell_types,
  marker_genes, missing_genes,
  gene_labels, feature_order, p3c
)

gc()

################ Figure 4a: Cell types by tissue ################
# Load the annotated Seurat object
seob <- readRDS(
  file.path(data_root, "seurat", "seob_celltype.rds")
)

# Define cell type colors
cell_type_colors <- c(
  "Neurons" = "#3C77AF",
  "Mesenchymal cells" = "#7DBFA7",
  "Epithelial cells" = "#EE934E",
  "Muscle cells" = "#D1352B",
  "Immune cells" = "#9B5B33",
  "Ciliated cells" = "#B383B9",
  "Proliferating cells" = "#BBDD78",
  "Uncertain" = "#8FA4AE"
)

umap <- Embeddings(seob, "umap.harmony")
meta <- seob@meta.data[rownames(umap), ]

plot_df <- data.frame(
  UMAP_1 = umap[, 1],
  UMAP_2 = umap[, 2],
  tissue = factor(
    ifelse(
      meta$sample %in% c("ner_1", "ner_2"),
      "Nerve ring",
      "Tentacles"
    ),
    levels = c("Nerve ring", "Tentacles")
  ),
  cell_type = meta$cell_type
)

x_range <- range(plot_df$UMAP_1)
y_range <- range(plot_df$UMAP_2)

x0 <- x_range[1] - diff(x_range) * 0.03
y0 <- y_range[1] - diff(y_range) * 0.03

x_len <- diff(x_range) * 0.06
y_len <- diff(y_range) * 0.06

axis_df <- data.frame(
  tissue = factor(
    "Nerve ring",
    levels = c("Nerve ring", "Tentacles")
  ),
  x0 = x0,
  y0 = y0,
  x_len = x_len,
  y_len = y_len
)

# Plot cell types by tissue
p4a <- ggplot(
  plot_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = cell_type
  )
) +
  geom_point(
    size = 0.12,
    alpha = 0.9
  ) +
  facet_wrap(
    ~tissue,
    nrow = 1
  ) +
  scale_color_manual(
    values = cell_type_colors,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.07, 0.05))
  ) +
  coord_equal(clip = "off") +

  geom_segment(
    data = axis_df,
    aes(
      x = x0,
      y = y0,
      xend = x0 + x_len,
      yend = y0
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +

  geom_segment(
    data = axis_df,
    aes(
      x = x0,
      y = y0,
      xend = x0,
      yend = y0 + y_len
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(0.10, "cm")),
    linewidth = 0.5
  ) +

  geom_text(
    data = axis_df,
    aes(
      x = x0 + x_len * 0.55,
      y = y0 - y_len * 0.50
    ),
    label = "UMAP 1",
    inherit.aes = FALSE,
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +

  geom_text(
    data = axis_df,
    aes(
      x = x0 - x_len * 0.50,
      y = y0 + y_len * 0.55
    ),
    label = "UMAP 2",
    inherit.aes = FALSE,
    angle = 90,
    size = 2.6,
    fontface = "bold",
    family = "Arial"
  ) +
  
  labs(
    title = "Cell types by tissue",
    x = NULL,
    y = NULL,
    color = "Cell type"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 3,
        alpha = 1
      )
    )
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    plot.title.position = "plot",
    
    strip.text = element_text(
      size = 12,
      face = "bold",
      margin = margin(t = 4, b = 8)
    ),
    strip.background = element_blank(),
    
    panel.spacing.x = grid::unit(8, "mm"),
    
    legend.title = element_text(
      size = 11,
      face = "bold"
    ),
    legend.text = element_text(
      size = 10,
      face = "bold"
    ),
    legend.position = "right",
    legend.key.size = grid::unit(5, "mm"),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.margin = margin(
      t = 6,
      r = 6,
      b = 6,
      l = 6
    )
  )

p4a

rm(
  seob, cell_type_colors,
  umap, meta, plot_df, x_range, y_range,
  x0, y0, x_len, y_len, axis_df, p4a
)

gc()

################ Figure 4b: Cell type composition across samples ################
# Load the annotated Seurat object
seob <- readRDS(
  file.path(data_root, "seurat", "seob_celltype.rds")
)

sample_colors <- c(
  "ner_1" = "#5B99C5",
  "ner_2" = "#FAA256",
  "ten_1" = "#BADABA",
  "ten_2" = "#F0AFAF"
)

meta <- seob@meta.data

meta$sample <- factor(
  meta$sample,
  levels = names(sample_colors)
)

# Calculate cell type counts by sample
composition <- as.data.frame(
  table(
    sample = meta$sample,
    cell_type = meta$cell_type
  )
)

names(composition)[3] <- "cells"

sample_totals <- aggregate(
  cells ~ sample,
  data = composition,
  FUN = sum
)

composition$total <- sample_totals$cells[
  match(composition$sample, sample_totals$sample)
]

	# Calculate cell type proportions
composition$percentage <- (
  composition$cells / composition$total
) * 100

composition$label <- ifelse(
  composition$cells == 0,
  "0%",
  ifelse(
    composition$percentage < 1,
    sprintf("%.2f%%", composition$percentage),
    sprintf("%.1f%%", composition$percentage)
  )
)

composition$label_x <- ifelse(
  composition$cells == 0,
  0.005,
  composition$percentage
)

composition$label_hjust <- ifelse(
  composition$cells == 0,
  0,
  -0.15
)

sample_y <- c(
  "ner_1" = 2.45,
  "ner_2" = 2.05,
  "ten_1" = 1.65,
  "ten_2" = 1.25
)

bar_half_height <- 0.145

composition$sample_y <- unname(
  sample_y[as.character(composition$sample)]
)

print(composition)

# Plot cell type composition across samples
p4b <- ggplot(
  composition,
  aes(fill = sample)
) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = percentage,
      ymin = sample_y - bar_half_height,
      ymax = sample_y + bar_half_height
    )
  ) +
  geom_text(
    aes(
      x = label_x,
      y = sample_y,
      label = label,
      hjust = label_hjust
    ),
    size = 3,
    fontface = "bold",
    family = "Arial"
  ) +
  facet_wrap(
    ~cell_type,
    nrow = 2,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = sample_colors
  ) +
  scale_y_continuous(
    breaks = unname(sample_y),
    labels = names(sample_y),
    limits = c(0.98, 2.72),
    expand = expansion(mult = 0)
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.28))
  ) +
  labs(
    title = "Cell type proportions across samples",
    x = "Proportion of cells (%)",
    y = "Sample"
  ) +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",
    
    strip.text = element_text(
      size = 12,
      face = "bold",
      margin = margin(b = 3)
    ),
    strip.background = element_blank(),
    
    axis.title = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text = element_text(
      size = 10,
      face = "bold"
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.25
    ),
    axis.ticks.length = grid::unit(0.08, "cm"),
    
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.80
    ),
    panel.spacing.x = grid::unit(5, "mm"),
    panel.spacing.y = grid::unit(4, "mm"),
    
    legend.position = "none",
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

p4b

rm(
  seob, sample_colors,
  meta, composition, sample_totals, sample_y,
  bar_half_height, p4b
)

gc()

################ Figure 5: KEGG pathway enrichment ################
# Load significant KEGG enrichment results
kegg_sig <- read.csv(
  file.path(result_root, "kegg", "kegg_significant.csv")
)

cell_types <- c(
  "Neurons",
  "Mesenchymal cells",
  "Epithelial cells",
  "Muscle cells",
  "Immune cells",
  "Ciliated cells",
  "Proliferating cells"
)

# Select the top 10 significant pathways per cell type
fig5_data <- kegg_sig %>%
  mutate(
    cell_type = factor(cell_type, levels = cell_types)
  ) %>%
  group_by(cell_type) %>%
  slice_min(
    order_by = p.adjust,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(cell_type, p.adjust) %>%
  mutate(
    Description = factor(
      Description,
      levels = rev(unique(Description))
    ),
    log10_fdr = -log10(p.adjust)
  )

table(fig5_data$cell_type)

nrow(fig5_data)
n_distinct(fig5_data$Description)

# Plot KEGG enrichment across cell types
p5 <- ggplot(
  fig5_data,
  aes(
    x = cell_type,
    y = Description,
    size = Count,
    color = log10_fdr
  )
) +
  geom_point(
    alpha = 0.95
  ) +

  scale_size_continuous(
    range = c(2, 6.5),
    breaks = c(10, 30, 50, 70)
  ) +

  scale_color_gradientn(
    colours = c(
      "#785199",
      "#B383B9",
      "#EE934E",
      "#FCED82"
    ),
    transform = "log10",
    breaks = c(2, 5, 10, 20, 50)
  ) +

  scale_x_discrete(
    expand = expansion(add = 0.6)
  ) +
  scale_y_discrete(
    expand = expansion(add = 0.6)
  ) +

  guides(
    color = guide_colorbar(
      order = 1,
      title.position = "left",
      barheight = grid::unit(35, "mm"),
      barwidth = grid::unit(5, "mm")
    ),
    size = guide_legend(
      order = 2,
      title.position = "left"
    )
  ) +

  labs(
    title = "KEGG pathway enrichment across cell types",
    x = "Cell type",
    y = "KEGG pathway",
    size = "Gene count",
    color = "-log10(adjusted P value)"
  ) +

  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.title.position = "plot",

    axis.title = element_text(
      size = 12,
      face = "bold"
    ),

    axis.text.x = element_text(
      size = 10,
      face = "bold",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),

    axis.text.y = element_text(
      size = 10,
      face = "bold"
    ),

    axis.ticks.y = element_blank(),
    axis.ticks.length = grid::unit(0.08, "cm"),

    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),

    legend.title = element_text(
      size = 11,
      face = "bold",
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.text = element_text(
      size = 10
    ),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = grid::unit(6, "mm"),
    
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )

p5

rm(kegg_sig, cell_types, fig5_data, p5)

gc()
