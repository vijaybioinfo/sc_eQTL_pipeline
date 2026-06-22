# 03 · QC report

**Script:** `scripts/03_qc_report.R`  
**Rule:** `qc_report` in `workflow/rules/report.smk`

Generates an HTML QC report using `ggplot2` and `rmarkdown`, summarizing sample distributions before and after normalization and filtering.

## Output

`results/qc_report.html` — a self-contained HTML file you can share with collaborators directly.

## Plots included

- Raw vs normalized count distributions (boxplot per sample)
- PCA of normalized counts colored by sample group
- Gene filtering summary (how many genes were removed and why)
