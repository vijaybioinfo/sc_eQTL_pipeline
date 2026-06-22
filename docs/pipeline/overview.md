# Pipeline overview

The pipeline is defined in a single `Snakefile` and runs four sequential rules. Dependencies between rules are resolved automatically by Snakemake based on input/output file paths.
It identifies expression quantitative trait loci (eQTLs) across multiple tissue/cell type/subset combinations. It is design to run three parallel QTL mapping approaches (MatrixEQTL + permutation correction, MatrixEQTL + EigenMT correction and FastQTL permutation correction) and includes statistical fine-mapping to identify credible causal variants.

The pipeline has four major stages:

1. **Preprocessing:** Expression data is normalized and formatted for both MatrixEQTL and FastQTL. Genotype data is subset to relevant donors, LD-pruned, and used to compute genotype PCs. PEER factors are estimated as latent technical covariates. All covariates are merged into a single matrix.

2. **MatrixEQTL mapping:** Cis-eQTL associations are tested per chromosome. A null distribution is built via permutations, and FDR is calculated using either EigenMT or permutation-based methods to account for multiple testing. Significant pairs are compiled into a summary table with allele frequencies and gene annotations.

3. **Fine-mapping:** Significant eQTLs are fine-mapped using a credible sets approach (SuSiER), outputting credible sets and their constituent variants. Pipeline only generate finemapping results from MatrixEQTL + EigenMT approach, further improvement is needed for the other two approaches.

4. **FastQTL mapping:** Run in parallel: nominal pass (all pairs) and permutation pass (per-gene beta-approximated p-values), producing gzipped outputs for downstream use.

## DAG of rules

<img src="../eqtl_pipeline_overview.svg"/>
