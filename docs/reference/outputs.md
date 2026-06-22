# Output files

All outputs are written to `results/` (configurable via `results_dir` in `config.yaml`).

| File | Produced by | Description |
|------|-------------|-------------|
| `results/counts_raw.csv` | featureCounts (external) | Raw gene × sample count matrix |
| `results/counts_norm.csv` | `01_normalize.R` | Log-CPM normalized counts |
| `results/counts_filtered.csv` | `02_filter.R` | Filtered count matrix |
| `results/qc_report.html` | `03_qc_report.R` | Interactive QC summary report |
| `logs/*.log` | Each rule | Per-rule log files for debugging |
