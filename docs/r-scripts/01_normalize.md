# 01 · Normalize

**Script:** `scripts/01_normalize.R`  
**Rule:** `normalize` in `workflow/rules/normalize.smk`

## What it does

Reads the raw count matrix produced by featureCounts, applies TMM normalization via `edgeR`, and writes log-CPM values to disk.

## Inputs & outputs

| | Path |
|--|------|
| **Input** | `results/counts_raw.csv` |
| **Output** | `results/counts_norm.csv` |
| **Log** | `logs/normalize.log` |

## Script walkthrough

```r
# Access Snakemake inputs/outputs via the snakemake object
counts_path <- snakemake@input[["counts"]]
out_path    <- snakemake@output[["norm"]]

library(edgeR)

# Load raw counts (genes × samples)
counts <- read.csv(counts_path, row.names = "gene_id")

# Build DGEList and compute normalization factors
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method = snakemake@config[["normalization"]])

# Export log-CPM matrix
norm <- cpm(dge, log = TRUE)
write.csv(norm, out_path)
```

!!! tip
    The `snakemake@config[["normalization"]]` call reads the `normalization` key directly from `config/config.yaml`, so you never hardcode it in the script.

## Changing the normalization method

Edit `config/config.yaml`:

```yaml
normalization: RLE   # or TMM, upperquartile
```
