# Pseudobulk Expression Profile Input

**Script:** `bin/preprocess/input_seurat.R`  

## Overview

This script takes a Seurat single-cell RNA-seq object and generates **pseudobulk expression profiles** per donor, both per cluster and aggregated across all clusters for a given cell type. It is intended to be run **prior to the eQTL pipeline** to produce the expression input files required by `Snakefile`.

For each cluster (and for the full aggregated cell type), it computes:

- **Mean gene expression** per donor (pseudobulk profile)
- **Proportion of expressing cells** per gene per donor
- **Cell counts** per donor

Donors with fewer than 10 cells in a given cluster are automatically excluded.

---

## Requirements

### R packages

| Package | Purpose |
|---|---|
| `Seurat` | Single-cell object handling and average expression |
| `Matrix` | Sparse matrix operations for expression frequency |
| `data.table` | Fast file writing |
| `optparse` | Command-line argument parsing |

Install with:

```r
install.packages(c("Matrix", "data.table", "optparse"))
# Seurat: https://satijalab.org/seurat/articles/install.html
```

---

## Usage

```bash
Rscript input_seurat.R \
  -s <seurat_object.RDS> \
  -d <donor_annotation.csv> \
  -c <cell_type_label> \
  -r <resolution_column> \
  -o <output_path>
```

### Arguments

| Flag | Long form | Description |
|---|---|---|
| `-s` | `--seuratobj` | Path to the `.RDS` Seurat object |
| `-d` | `--donorannot` | Path to donor annotation CSV (see format below) |
| `-c` | `--cell` | Cell type label (used for naming output directories) |
| `-r` | `--resolutioncolumn` | Metadata column in the Seurat object containing cluster assignments (e.g. `RNA_snn_res.0.2`) |
| `-o` | `--outputpath` | Root output directory where results will be saved |

---

## Input file formats

### Seurat object (`.RDS`)
A standard Seurat object with:

- Raw or normalized counts in the `RNA` assay
- A metadata column matching `--resolutioncolumn` containing cluster labels
- A `barcode` column-compatible row index in `@meta.data`

### Donor annotation CSV (`--donorannot`)
A CSV file with at minimum the following columns:

| Column | Description |
|---|---|
| `barcode` | Cell barcode matching the Seurat object |
| `full.donor.id.tag` | Numeric donor ID (will be zero-padded to 4 digits and prefixed with `DLCP`, e.g. `DLCP0042`) |
| `population.tag` | Cell population label |

> **Note:** Barcodes labeled `Doublet` in `full.donor.id.tag` are automatically removed.

---

## Output structure

For each cluster defined by `--resolutioncolumn`, and for the full aggregated cell type (`AGGR`), the script creates a subdirectory:

```
<outputpath>/
  <cell>/
    <cluster>/
      gene_expression.txt    # Genes × donors mean expression matrix (pseudobulk)
      freq_cell.txt          # Genes × donors proportion of expressing cells
      cellsxdonor.txt        # Number of cells per donor in this cluster
      donors.txt             # List of donor IDs (one per line)
    AGGR/
      gene_expression.txt    # Same files but aggregated across all clusters
      freq_cell.txt
      cellsxdonor.txt
      donors.txt
```


### Output file details

**`gene_expression.txt`** — tab-separated, genes as rows, donors as columns. Contains the mean expression value per gene per donor (pseudobulk). This is the primary input for the eQTL pipeline (`expFile` in the dataset table).

**`freq_cell.txt`** — tab-separated, genes as rows, donors as columns. Each value is the proportion of cells from that donor where the gene has expression > 0. Used as the `prop` (frequency table) input in the pipeline.

**`cellsxdonor.txt`** — two-column table (`donor`, `n.cells`). Donors with fewer than 10 cells are excluded before this is written.

**`donors.txt`** — plain text list of donor IDs in the same order as the columns of `gene_expression.txt`. Used as the `donorFile` input in the pipeline.

---

## Important notes

- **Cluster merging:** If a cluster label contains an underscore (e.g. `2_3`), the script interprets it as a merge of clusters 2 and 3 and pools those cells together before computing statistics.
- **Donor ID formatting:** Donor IDs are zero-padded to 4 digits and prefixed with `DLCP` automatically (e.g. `42` → `DLCP0042`). Make sure your downstream pipeline dataset table uses this same format.
- **Minimum cell filter:** Any donor contributing fewer than 10 cells to a given cluster is dropped from that cluster's outputs. This threshold is hardcoded.
- **`AGGR` folder:** Always generated regardless of cluster structure — represents the full cell type without cluster subdivision.