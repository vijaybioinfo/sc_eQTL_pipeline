# Preprocessing Pipeline

A reproducible **Snakemake + R** pipeline for sc-eQTL analysis, from single cell counts to finemapping results.

## What this pipeline does

This pipeline takes raw single cell RNA-seq data counts and genotype data as input, and performs eQTL analysis on pseudobulk-profiles following three different technical (and consistent) tool combinations:

1. **PEER + MatrixeQTL + EigenMT** (default) 
2. **PEER + MatrixeQTL + permutation tests**
3. **PEER + FastQTL + permutation tests**

It summarizes the results in a single table, filters significant hits, and performs finemapping with SuSiE. The pipeline is modular and can be easily adapted to different datasets and analysis needs.

## Requirements

- Snakemake `> 6.0, < 8.0` (Compatible with Snakemake 7.x)
- R ≥ 4.3
- [PEER](https://github.com/pmbio/peer)
- [Plink](https://www.cog-genomics.org/plink/) `> 1.9`
- [EigenMT](https://github.com/joed3/eigenMT)
- bcftools
- Conda or Mamba (recommended)
- R packages: `arrow`, `data.table`, `edgeR`, `ggplot2`, `Matrix`, `MatrixEQTL`, `matrixStats`, `qvalue`, `Seurat`, `susieR`, `tidyverse`

