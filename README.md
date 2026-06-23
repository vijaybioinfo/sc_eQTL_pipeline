# sc-eQTL Pipeline

> A reproducible Snakemake + R pipeline for single-cell expression quantitative trait loci (sc-eQTL) analysis, developed by the [DICE](https://www.dice-database.org/) team at La Jolla Institute for Immunology (LJI).

[![Documentation](https://img.shields.io/badge/docs-mkdocs--material-blue)](https://vijaybioinfo.github.io/sc_eQTL_pipeline)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen)](https://snakemake.readthedocs.io)

---

## Overview

This pipeline identifies **expression quantitative trait loci (eQTLs)** — genetic variants that influence gene expression — across multiple tissue, cell type, and subset combinations derived from single-cell RNA-seq data.

It implements **three parallel QTL mapping strategies**:

- **MatrixEQTL + permutation-based FDR correction**
- **MatrixEQTL + EigenMT multiple testing correction**
- **FastQTL permutation pass** (per-gene beta-approximated p-values)

And includes **Bayesian fine-mapping** via SuSiE to identify credible causal variant sets.

---

## Pipeline stages

1. **Preprocessing** — Expression data is normalized and formatted. Genotype data is subset to relevant donors, LD-pruned, and used to compute genotype PCs. PEER factors are estimated as latent covariates. All covariates are merged into a single matrix.

2. **MatrixEQTL mapping** — Cis-eQTL associations are tested per chromosome. A null distribution is built via permutations and FDR is calculated using either EigenMT or permutation-based methods.

3. **FastQTL mapping** — Nominal pass (all pairs) and permutation pass (per-gene p-values) run in parallel.

4. **Fine-mapping** — Significant eQTLs are fine-mapped using a credible sets approach, outputting credible sets and their constituent variants with posterior inclusion probabilities (PIPs).

---

## Requirements

### Software

| Tool | Version | Purpose |
|---|---|---|
| [Snakemake](https://snakemake.readthedocs.io) | ≥ 7.0 | Workflow manager |
| [R](https://www.r-project.org/) | ≥ 4.1 | QTL mapping and statistics |
| [Python](https://www.python.org/) | ≥ 3.8 | EigenMT and FastQTL wrapper |
| [PLINK](https://www.cog-genomics.org/plink/) | 1.9 | LD pruning and PCA |
| [bcftools](https://samtools.github.io/bcftools/) | ≥ 1.12 | Genotype subsetting |
| [FastQTL](http://fastqtl.sourceforge.net/) | — | FastQTL mapping |
| [tabix / bgzip](http://www.htslib.org/) | — | BED file indexing |

### Conda environments

The pipeline uses two conda environments defined in the repo:

- `DLCP_v3` — R-based rules (MatrixEQTL, FDR, fine-mapping)
- `pyEigenMT` — Python-based EigenMT correction

---

## Installation

```bash
# Clone the repository
git clone https://github.com/vijaybioinfo/sc_eQTL_pipeline.git
cd sc_eQTL_pipeline

# Create conda environments
conda env create -f envs/DLCP_v3.yaml
conda env create -f envs/pyEigenMT.yaml
```

---

## Quick start

### 1. Generate pseudobulk input files

Before running the pipeline, generate per-donor pseudobulk expression profiles from your Seurat object:

```bash
Rscript bin/input_seurat.R \
  -s path/to/seurat_object.RDS \
  -d path/to/donor_annotation.csv \
  -c <cell_type> \
  -r <resolution_column> \
  -o path/to/output/
```

This script was designed to our own Seurat objects, but can be adapted to other single-cell RNA-seq datasets.

### 2. Configure the pipeline

Edit `snake_conf.yaml` to set your input files, parameters, and working directory. See the [Configuration docs](https://vijaybioinfo.github.io/sc_eQTL_pipeline/pipeline/configuration/) for a full description of all options.

---

## Input files

| File | Description |
|---|---|
| Seurat object (`.RDS`) | Single-cell RNA-seq object per cell type |
| Donor annotation (`.csv`) | Maps barcodes to donor IDs |
| Genotype VCF (`.vcf.gz`) | Phased or unphased genotypes |
| SNP files table | Per-chromosome SNP and SNP location files |
| Gene annotation | Gene coordinates for cis-window definition |
| Covariates file | Additional covariates (e.g. sex, age, batch) |
| Dataset table (`.csv`) | Sample sheet linking tissue/cell/subset to input files |

---
> [!NOTE]
> ## Example Data
>
> The example dataset provided in this repository is **simulated** and does not
> represent real patient or donor information. It was generated synthetically for
> demonstration purposes only, to allow users to test and explore the pipeline
> without requiring access to the original data.
>
> No real genomic or personal information from study participants 
> is included or can be inferred from these files.
>
> To reproduce the results reported in the paper, access to the original dataset
> is required.

---

## Output files

For each `tissue / cell type / subset` combination:

```
results/
  matrix_eqtl/eQTL/<tissue>/<cell>/<subset>/
    credible_sets.txt               # Fine-mapped credible sets per eQTL
    credible_sets_variants.txt      # Per-variant PIPs within credible sets
    eQTL_sigpairs.csv               # Significant eQTL pairs
    all_pairs_adjust_pvalue.rds     # All pairs with adjusted p-values

  fastQTL/eQTL/<tissue>/<cell>/<subset>/
    Output_all.allpairs.txt.gz      # All nominal pairs
    Output_all.genes.txt.gz         # Per-gene permutation results
```

---

## Documentation

Full documentation is available at:

**[https://vijaybioinfo.github.io/sc_eQTL_pipeline](https://vijaybioinfo.github.io/sc_eQTL_pipeline)**

---

## Citation

If you use this pipeline in your work, please cite:

> Benjamin J. Schmiedel et al., Single-cell eQTL analysis of activated T cell subsets reveals activation and cell type–dependent effects of disease-risk variants. Sci. Immunol.7,eabm2508(2022).DOI:10.1126/sciimmunol.abm2508
> Benjamin J. Schmiedel et al., Tissue-resident immune cells drive genetic risk in autoimmune and lung diseases. Nat. immunol. 2026 (in press).

---

## Contact

For questions or issues, please open a [GitHub Issue](https://github.com/vijaybioinfo/sc_eQTL_pipeline/issues) or contact us at **[cgonzalez@lji.org](mailto:cgonzalez@lji.org)**.