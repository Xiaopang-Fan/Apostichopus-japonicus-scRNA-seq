# Single-cell RNA-seq analysis of the nerve ring and tentacles in *Apostichopus japonicus*

This repository contains the R scripts used to analyze and visualize single-cell RNA-seq data from the nerve ring and tentacles of the sea cucumber *Apostichopus japonicus*.

The dataset comprises four samples (`ner_1`, `ner_2`, `ten_1`, and `ten_2`). The analysis includes quality control, doublet removal, sample integration, cell clustering, cell type annotation, marker gene identification, and KEGG pathway enrichment.

## Scripts

- **analysis.R** — Performs single-cell data preprocessing, quality control, doublet detection using scDblFinder, Seurat and Harmony integration, clustering, marker identification, cell type annotation, and KEGG enrichment analysis.
- **plot.R** — Generates visualizations of quality control metrics, doublet detection, UMAP embeddings, marker gene expression, cell type composition, and KEGG enrichment results.

## Requirements

The scripts require R and the following packages:

- Seurat
- harmony
- ggplot2
- tidyverse
- scDblFinder
- clusterProfiler
- ComplexHeatmap
- circlize
- BiocParallel
- Matrix
- SeuratObject

## Input data

The analysis requires four 10x Genomics gene expression matrices and the corresponding gene annotation file.

The expected input directory structure is:

```text
input/
├── ner_1/
├── ner_2/
├── ten_1/
├── ten_2/
└── out.emapper.annotations
```

Each sample directory should contain a 10x-formatted gene expression matrix. The annotation file contains gene annotations generated using eggNOG-mapper.

These input files are not included in this GitHub repository.

## Usage

Run the scripts from the project root directory:

```r
source("analysis.R")
source("plot.R")
```

The scripts use the current working directory as the project root by default. Alternatively, the project root can be specified using the `SCRNA_PROJECT_ROOT` environment variable.

`analysis.R` creates the required output directories and saves intermediate Seurat objects and analysis results under `output/`.

`plot.R` uses these generated files to reproduce the visualizations. The analysis script should therefore be run before the plotting script when starting from raw count matrices.

## Data availability

Raw sequencing data have been deposited in the Genome Sequence Archive (GSA) and NCBI Sequence Read Archive (SRA). Accession numbers will be added upon public release.

The processed Seurat object and top 30 marker genes for each cluster have been deposited in Figshare (DOI: [10.6084/m9.figshare.35013703](https://doi.org/10.6084/m9.figshare.35013703)).

## License

The code in this repository is distributed under the MIT License. See the `LICENSE` file for details.
