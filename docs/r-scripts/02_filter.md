# 02 · Filter

**Script:** `scripts/02_filter.R`  
**Rule:** `filter` in `workflow/rules/filter.smk`

Removes low-expression genes using the thresholds defined in `config.yaml` (`min_counts` and `min_samples`).

```r
counts_path <- snakemake@input[["norm"]]
out_path    <- snakemake@output[["filtered"]]

min_counts  <- snakemake@config[["min_counts"]]
min_samples <- snakemake@config[["min_samples"]]

library(dplyr)

counts  <- read.csv(counts_path, row.names = "gene_id")
keep    <- rowSums(counts >= min_counts) >= min_samples
filtered <- counts[keep, ]

message(sprintf("Kept %d / %d genes after filtering", sum(keep), nrow(counts)))
write.csv(filtered, out_path)
```
