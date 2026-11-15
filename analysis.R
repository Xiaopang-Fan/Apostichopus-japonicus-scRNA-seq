# Load packages
library(Seurat)
library(harmony)
library(ggplot2)
library(tidyverse)
library(scDblFinder)
library(clusterProfiler)

set.seed(42)

# Set project paths
project_root <- Sys.getenv("SCRNA_PROJECT_ROOT", unset = getwd())
input_root <- file.path(project_root, "input")
data_root <- file.path(project_root, "output", "data")
result_root <- file.path(project_root, "output", "result")

# Create output directories
for (dir_path in c(
  file.path(data_root, "plot_data"),
  file.path(data_root, "seurat"),
  result_root,
  file.path(result_root, "kegg")
)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

################ Read input data ################
# Define sample metadata and 10x data directories
sample_info <- data.frame(
  sample = c(
    "ner_1", "ner_2", "ten_1", "ten_2"
  ),
  tissue = c(
    "ner", "ner", "ten", "ten"
  ),
  data_dir = c(
    file.path(input_root, "ner_1"),
    file.path(input_root, "ner_2"),
    file.path(input_root, "ten_1"),
    file.path(input_root, "ten_2")
  )
)

# Initialize a list for individual Seurat objects
seob_list <- list()

# Read 10x data and create Seurat objects
for (i in seq_len(nrow(sample_info))) {
  
  scrna_data <- Read10X(
    data.dir = sample_info$data_dir[i]
  )
  
  seob_sample <- CreateSeuratObject(
    counts = scrna_data,
    min.cells = 3,
    min.features = 200
  )
  
  seob_sample$sample <- sample_info$sample[i]
  seob_sample$tissue <- sample_info$tissue[i]
  
  seob_list[[sample_info$sample[i]]] <- seob_sample
}

# Check gene and cell counts for each sample
sapply(seob_list, dim)

# Merge four samples with sample-specific cell prefixes
seob <- merge(
  x = seob_list[[1]],
  y = seob_list[-1],
  add.cell.ids = names(seob_list)
)

# Check merged object dimensions
dim(seob)

# Load gene annotations
gene_info <- read_delim( file.path(input_root, "out.emapper.annotations"), 
                        delim = "\t", 
                        escape_double = FALSE, 
                        col_names = FALSE, 
                        comment = "#", 
                        trim_ws = TRUE) %>%
  dplyr::select(
    Gene_Id = X1,
    Gene_Symbol = X9,
    GO = X10,
    Ko = X12,
    Pathway = X13,
    COG = X7,
    Description = X8,
    Domain = X21)

gene_info <- gene_info %>%
  mutate(Gene_Id = Gene_Id %>%
           gsub("model", "TU", .) %>%
           str_replace("(\\.\\d+)\\.\\d+$", "\\1") %>%
           str_replace_all("_", "-")) %>%
  distinct(Gene_Id, .keep_all = TRUE)

# Remove intermediate objects
rm(
  seob_list, scrna_data, seob_sample, i
)

gc()

################ Quality control ################
# Merge RNA count layers
seob[["RNA"]] <- JoinLayers(seob[["RNA"]])

counts_matrix <- LayerData(
  seob,
  assay = "RNA",
  layer = "counts"
)

# Identify highly expressed genes
gene_stats <- data.frame(
  gene_id = rownames(counts_matrix),
  mean_exp = Matrix::rowMeans(counts_matrix),
  detection_rate = Matrix::rowMeans(counts_matrix > 0)
)


gene_stats <- gene_stats[
  order(-gene_stats$mean_exp),
]

head(gene_stats, 60)

# Exclude FTH1 and ND5 from the QC metric calculation
genes_to_exclude <- c(
  "evm.TU.Chr14.277",     # FTH1
  "evm.TU.Scaffold6.1"    # ND5
)

counts_clean <- counts_matrix[
  !rownames(counts_matrix) %in% genes_to_exclude,
]

# Calculate the maximum UMI proportion for each cell
max_umi <- numeric(ncol(counts_clean))

for (i in seq_len(ncol(counts_clean))) {
  
  start <- counts_clean@p[i] + 1L
  end <- counts_clean@p[i + 1L]
  
  if (start <= end) {
    max_umi[i] <- max(counts_clean@x[start:end])
  }
}

seob$max_umi_ratio_clean <- max_umi / Matrix::colSums(counts_clean)

probs <- c(
  seq(0.01, 0.10, 0.01),
  0.25, 0.50, 0.75,
  seq(0.90, 1.00, 0.01)
)

# Examine QC metric distributions by sample
for (s in unique(seob$sample)) {
  
  cat("\n=====", s, "=====\n")
  
  qc_quantile <- rbind(
    
    nFeature_RNA = sprintf(
      "%.0f",
      quantile(
        seob$nFeature_RNA[seob$sample == s],
        probs
      )
    ),
    
    nCount_RNA = sprintf(
      "%.0f",
      quantile(
        seob$nCount_RNA[seob$sample == s],
        probs
      )
    ),
    
    max_umi_ratio_clean = sprintf(
      "%.5f",
      quantile(
        seob$max_umi_ratio_clean[seob$sample == s],
        probs
      )
    )
  )
  
  colnames(qc_quantile) <- paste0(probs * 100, "%")
  
  print(qc_quantile, quote = FALSE)
}

qc_keep <- logical(ncol(seob))
names(qc_keep) <- colnames(seob)

qc_summary <- data.frame()

# Apply sample-specific QC thresholds
for (s in unique(seob$sample)) {
  
  idx <- seob$sample == s
  
  nFeature_max <- as.numeric(
    quantile(
      seob$nFeature_RNA[idx],0.99)
  )
  
  nCount_max <- as.numeric(
    quantile(
      seob$nCount_RNA[idx],0.99)
  )
  
  max_ratio_min <- as.numeric(
    quantile(
      seob$max_umi_ratio_clean[idx],0.05)
  )
  
  sample_keep <- (
    seob$nFeature_RNA[idx] > 200 &
      seob$nFeature_RNA[idx] < nFeature_max &
      seob$nCount_RNA[idx] < nCount_max &
      seob$max_umi_ratio_clean[idx] > max_ratio_min
  )
  
  qc_keep[idx] <- sample_keep
  
  qc_summary <- rbind(
    qc_summary,
    data.frame(
      sample = s,
      before = sum(idx),
      nFeature_max = nFeature_max,
      nCount_max = nCount_max,
      max_ratio_min = max_ratio_min,
      nFeature_low = sum(seob$nFeature_RNA[idx] <= 200),
      nFeature_high = sum(seob$nFeature_RNA[idx] >= nFeature_max),
      nCount_high = sum(seob$nCount_RNA[idx] >= nCount_max),
      max_ratio_low = sum(seob$max_umi_ratio_clean[idx] <= max_ratio_min),
      after = sum(sample_keep)
    )
  )
}

qc_summary$removed <- qc_summary$before - qc_summary$after

qc_summary$removed_percent <- round(
  qc_summary$removed / qc_summary$before * 100, 2)

qc_summary$nFeature_max <- round(
  qc_summary$nFeature_max, 2)

qc_summary$nCount_max <- round(
  qc_summary$nCount_max, 2)

qc_summary$max_ratio_min <- round(
  qc_summary$max_ratio_min, 5)

qc_summary

# Save pre-QC metadata for visualization
meta_before <- seob@meta.data

saveRDS(
  meta_before[, c(
    "sample", "nFeature_RNA", "nCount_RNA",
    "max_umi_ratio_clean"
  )],
  file.path(data_root, "plot_data", "qc_meta_before.rds")
)

# Filter cells based on QC criteria
seob <- subset(
  seob,
  cells = colnames(seob)[qc_keep]
)

dim(seob)

table(seob$sample)

# Save QC results
saveRDS(
  seob,
  file = file.path(data_root, "seurat", "seob_qc.rds")
)

write.csv(
  qc_summary,
  file = file.path(data_root, "seurat", "qc_summary.csv"),
  row.names = FALSE
)

# Remove intermediate objects
rm(
  sample_info, counts_matrix, gene_stats, genes_to_exclude,
  counts_clean, max_umi, start, end, probs, qc_quantile,
  qc_keep, qc_summary, meta_before, idx, nFeature_max,
  nCount_max, max_ratio_min, sample_keep, s, i
)

gc()

################ Doublet detection and filtering ################
# Extract RNA counts 
counts <- LayerData(
  seob,
  assay = "RNA",
  layer = "counts"
)

sample_ids <- seob@meta.data[colnames(counts), "sample"]

# Identify doublets for each sample using scDblFinder
sce <- scDblFinder(
  counts,
  samples = sample_ids,
  dbr.per1k = 0.008,
  BPPARAM = BiocParallel::SerialParam(RNGseed = 42)
)

doublet_results <- data.frame(
  barcode = colnames(sce),
  sample = seob@meta.data[colnames(sce), "sample"],
  score = sce$scDblFinder.score,
  class = as.character(sce$scDblFinder.class),
  row.names = colnames(sce)
)

doublet_results <- doublet_results[colnames(seob), ]

seob$scDblFinder.score <- doublet_results$score
seob$scDblFinder.class <- doublet_results$class

# Summarize doublet detection results
doublet_summary <- as.data.frame.matrix(
  table(
    sample = seob$sample,
    class = seob$scDblFinder.class
  )
)

doublet_summary$before <- rowSums(doublet_summary)
doublet_summary$after <- doublet_summary$singlet
doublet_summary$doublet_percent <- round(
  doublet_summary$doublet / doublet_summary$before * 100, 2
)

doublet_summary[, c(
  "before", "doublet", "doublet_percent", "after"
)]

table(seob$scDblFinder.class)

# Save doublet detection results
saveRDS(
  doublet_results,
  file = file.path(data_root, "seurat", "doublet_results.rds")
)

# Retain predicted singlets
seob <- subset(
  seob,
  cells = colnames(seob)[seob$scDblFinder.class == "singlet"]
)

dim(seob)

table(seob$sample)

# Save the filtered Seurat object
saveRDS(
  seob,
  file = file.path(data_root, "seurat", "seob_singlet.rds")
)

rm(
  counts, sample_ids, sce,
  doublet_results, doublet_summary
)

gc()

################ Data preprocessing, integration, and clustering ################
# Normalize RNA expression
seob <- NormalizeData(seob)

Layers(seob[["RNA"]])

seob <- FindVariableFeatures(
  seob,
  selection.method = "vst",
  nfeatures = 2000
)

seob <- ScaleData(seob)

seob <- RunPCA(seob)

# Integrate samples using Harmony
seob <- RunHarmony(
  seob,
  group.by.vars = "sample"
)

seob <- FindNeighbors(
  seob,
  reduction = "harmony",
  dims = 1:50
)

# Evaluate clustering resolutions
resolutions <- c(
  0.3, 0.4, 0.5, 0.6, 0.7,
  0.8, 0.9, 1.0, 1.1, 1.2, 
  1.3, 1.4, 1.5
)

for (res in resolutions) {
  seob <- FindClusters(
    seob,
    graph.name = "RNA_snn",
    resolution = res,
    random.seed = 42
  )
  
  col_name <- paste0("harmony_res_", res)
  seob@meta.data[[col_name]] <- seob$seurat_clusters
}

harmony_res_cols <- paste0("harmony_res_", resolutions)

sapply(
  harmony_res_cols,
  function(x) length(unique(seob@meta.data[[x]]))
)

# Select the final clustering resolution
seob$harmony_clusters <- seob$harmony_res_0.8

Idents(seob) <- "harmony_clusters"

table(seob$harmony_clusters)

# Run UMAP using Harmony embeddings
seob <- RunUMAP(
  seob,
  reduction = "harmony",
  dims = 1:50,
  reduction.name = "umap.harmony"
)

# Save the clustered Seurat object
saveRDS(
  seob,
  file = file.path(data_root, "seurat", "seob_clustered.rds")
)


############### Marker identification ###############
DefaultAssay(seob) <- "RNA"

# Identify positive markers for each cluster
cluster_markers <- FindAllMarkers(
  seob, only.pos = TRUE,
  min.pct = 0.25, logfc.threshold = 0.25
)

# Add gene annotations
cluster_markers <- cluster_markers %>%
  left_join(gene_info, by = c("gene" = "Gene_Id"))

# Select the top 30 significant markers per cluster
top30_markers <- cluster_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC, n = 30,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    cluster, gene,
    Gene_Symbol, Description,
    avg_log2FC, pct.1, pct.2, p_val_adj
  )

# Select the top 30 annotated markers per cluster
top30_named_markers <- cluster_markers %>%
  filter(
    p_val_adj < 0.05,
    !is.na(Gene_Symbol),
    Gene_Symbol != "",
    Gene_Symbol != "-"
  ) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 30) %>%
  ungroup() %>%
  select(
    cluster, gene,
    Gene_Symbol, Description,
    avg_log2FC, pct.1, pct.2, p_val_adj
  )

# Save marker results
write.csv(
  cluster_markers, row.names = FALSE, 
  file = file.path(result_root, "cluster_markers.csv")
)

write.csv(
  top30_markers, row.names = FALSE, 
  file = file.path(result_root, "top30_markers.csv")
)

write.csv(
  top30_named_markers, row.names = FALSE, 
  file = file.path(result_root, "top30_named_markers.csv")
)

################ Cell type annotation ################
# Define cluster-to-cell-type assignments
cell_type_clusters <- list(
  "Neurons" = c(0, 7, 9, 18, 25, 26, 29, 30, 31),
  "Mesenchymal cells" = c(1, 3, 4, 5, 6, 13, 16, 17),
  "Epithelial cells" = c(2, 8, 10, 11, 12, 21, 27),
  "Muscle cells" = c(14, 22),
  "Immune cells" = c(15, 23, 24, 28),
  "Ciliated cells" = 19,
  "Proliferating cells" = 20,
  "Uncertain" = 32
)

cluster_to_type <- rep(
  names(cell_type_clusters),
  lengths(cell_type_clusters)
)

names(cluster_to_type) <- as.character(
  unlist(cell_type_clusters, use.names = FALSE)
)

# Assign cell type annotations
seob$cell_type <- factor(
  unname(cluster_to_type[as.character(seob$harmony_clusters)]),
  levels = names(cell_type_clusters)
)

table(seob$cell_type, useNA = "ifany")

# Save the annotated Seurat object
saveRDS(
  seob,
  file = file.path(data_root, "seurat", "seob_celltype.rds")
)

############### Cell type marker analysis ###############
DefaultAssay(seob) <- "RNA"

Idents(seob) <- "cell_type"

# Identify positive markers for each cell type
cell_type_markers <- FindAllMarkers(
  seob,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Add gene annotations
cell_type_markers <- cell_type_markers %>%
  left_join(gene_info, by = c("gene" = "Gene_Id"))

# Select significant markers for seven cell types
cell_type_markers_sig <- cell_type_markers %>%
  filter(
    cluster != "Uncertain",
    p_val_adj < 0.05
)

table(cell_type_markers_sig$cluster)

# Save cell type marker results
write.csv(
  cell_type_markers,
  file = file.path(result_root, "cell_type_markers.csv"),
  row.names = FALSE
)

# Save significant markers for KEGG enrichment
write.csv(
  cell_type_markers_sig,
  file = file.path(result_root, "cell_type_markers_sig.csv"),
  row.names = FALSE
)

nrow(cell_type_markers_sig)

################ KEGG pathway enrichment analysis ################
# Build KEGG pathway-to-gene mapping
pathway2gene <- gene_info %>%
  dplyr::select(
    path_id = Pathway,
    gene_id = Gene_Id
  ) %>%
  separate_rows(path_id, sep = ",") %>%
  filter(str_detect(path_id, "^ko[0-9]{5}$")) %>%
  mutate(path_id = str_remove(path_id, "^ko")) %>%
  distinct()

# Retrieve KEGG pathway names
get_path2name <- function() {
  keggpathid2name.df <- clusterProfiler:::kegg_list("pathway")
  keggpathid2name.df[, 1] <- gsub("map", "", keggpathid2name.df[, 1])
  colnames(keggpathid2name.df) <- c("path_id", "path_name")
  return(keggpathid2name.df)
}

pathway2name <- get_path2name()

# Define the background using detected genes
background_genes <- rownames(seob)[
  Matrix::rowSums(
    SeuratObject::LayerData(
      seob,
      assay = "RNA",
      layer = "counts"
    ) > 0
  ) > 0
]

cell_type_order <- setdiff(
  levels(seob$cell_type),
  "Uncertain"
)

kegg_genes <- split(
  cell_type_markers_sig$gene,
  factor(
    cell_type_markers_sig$cluster,
    levels = cell_type_order
  ),
  drop = TRUE
)

sapply(
  kegg_genes,
  function(x) length(intersect(x, pathway2gene$gene_id))
)

length(intersect(background_genes, pathway2gene$gene_id))

# Perform KEGG enrichment for seven cell types
kegg_list <- lapply(kegg_genes, function(genes) {
  clusterProfiler::enricher(
    gene = genes,
    universe = background_genes,
    TERM2GENE = pathway2gene,
    TERM2NAME = pathway2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
})

kegg_results <- dplyr::bind_rows(
  lapply(names(kegg_list), function(type) {
    result <- as.data.frame(kegg_list[[type]])
    if (nrow(result) > 0) result$cell_type <- type
    result
  })
)

kegg_summary <- data.frame(
  cell_type = names(kegg_list),
  Total = sapply(
    kegg_list,
    function(x) nrow(as.data.frame(x))
  ),
  P_0.05 = sapply(
    kegg_list,
    function(x) sum(as.data.frame(x)$pvalue < 0.05)
  ),
  FDR_0.05 = sapply(
    kegg_list,
    function(x) sum(as.data.frame(x)$p.adjust < 0.05)
  ),
  row.names = NULL
)

kegg_summary

# Filter significant KEGG pathways
kegg_sig <- kegg_results %>%
  filter(p.adjust < 0.05)

nrow(kegg_sig)

n_distinct(kegg_sig$Description)

# Save KEGG enrichment results
write.csv(
  kegg_results,
  file = file.path(result_root, "kegg", "kegg_results.csv"),
  row.names = FALSE
)

write.csv(
  kegg_summary,
  file = file.path(result_root, "kegg", "kegg_summary.csv"),
  row.names = FALSE
)

write.csv(
  kegg_sig,
  file = file.path(result_root, "kegg", "kegg_significant.csv"),
  row.names = FALSE
)

rm(
  pathway2gene, pathway2name, get_path2name,
  background_genes, cell_type_order,
  kegg_genes, kegg_list
)

gc()
