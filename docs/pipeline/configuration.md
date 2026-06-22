# Running the pipeline

Our pipeline was run on a high-performance computing cluster using the SLURM workload manager. 

We used the following command to run the pipeline:

```bash
cd /path/to/working/directory
mkdir -p logs
ln -s /path/to/bin bin # Create a symbolic link to the bin directory containing the scripts
ln -s /path/to/Snakefile Snakefile # Create a symbolic link to the Snakefile
ln -s /path/to/config.yaml config.yaml # Create a symbolic link to the configuration file
mamba activate DLCP_v3 # Make sure to activate the appropriate conda environment

snakemake -s Snakefile --profile /path/to/profiles/slurm --configfile /path/to/config.yaml --stats logs/snakemake.stats >& logs/snakemake.log # Add -np options for a dry run to check for any potential issues before actual execution

```

For specific details on the SLURM profile used, please refer to the `profiles/slurm` directory in the repository. The `--stats` option generates a statistics file that provides insights into the execution of the pipeline, while the `-np` option allows for a dry run to check for any potential issues before actual execution. For further details on the available options and configurations, please refer to the Snakemake documentation.
# Configuration

All pipeline parameters live in `snake_config.yaml`. Edit this file before running.

## Full config reference

```yaml
# snake_config.yaml

inpFiles:
  vcf_genotype: /path/to/genotype.vcf.gz
  datasets: /path/to/datasets.csv
  covariates: /path/to/covariates.txt
  snpfiles: /path/to/snpfiles.txt
  geneannot: /path/to/geneannot.txt
  geneannot_complete: /path/to/geneannot_complete.txt
  rstables: /path/to/rstables/

fastqtl_parameters:
  chunks: 100
  ma_sample: 10
  cis_window: 1000000
  permute_down: 1000
  permute_up: 10000

parameters:
  MAF: 0.05
  PCS: 1
  Factors: 2
  permutations: 1
  covariates: ['Gender', 'Age']
  run_var_genes: False ## True or False
  n_var_genes: NA ### if 'run_var_genes' == True, 'n_var_genes' int > 0 & < number genes in the matrix. If 'run_var_genes' == False, 'n_var_genes' == NA

config:
  workdir: /path/to/workdir/
  peertool: /path/to/peertool
```

## Parameter descriptions

| Parameter | Default | Description |
|-----------|---------|-------------|
| `inpFiles/vcf_genotype` | /path/to/genotype.vcf.gz | Path to the input VCF file containing genotype data. The VCF file should be indexed and compressed with bgzip. |
| `inpFiles/datasets` | /path/to/datasets.csv | Path to a CSV file containing the list of datasets to analyze. The CSV file should have a column named 'dataset' with the names of the datasets. |
| `inpFiles/covariates` | /path/to/covariates.txt | Path to a tab-delimited text file containing covariate information for the samples. The file should have a header row with column names corresponding to the covariates (e.g., 'Gender', 'Age'). |
| `inpFiles/snpfiles` | /path/to/snps_chrs_files_maf0.05.txt | Path to a text file containing the list of SNP files for each chromosome. |
| `inpFiles/geneannot` | /path/to/geneAnnot.txt | Path to a text file containing gene annotation information with matrixeQTL format. Because our analysis is set around TSS, the start and end positions should be specified relative to the transcription start site. |
| `inpFiles/geneannot_complete` | /path/to/geneAnnot_complete.txt | Path to a text file containing complete gene annotation information with matrixeQTL format. This file should contain the full gene coordinates (start and end positions) relative to the reference genome. |
| `inpFiles/rstables` | /path/to/rstables/ | Path to a directory containing rsID tables for each chromosome. It is use for final annotation of significant eQTLs with rsIDs. |
| `parameters/MAF` | 0.05 | Minor allele frequency threshold for filtering variants. |
| `parameters/PCS` | 1 | Number of genotype principal components to include as covariates in the eQTL analysis. |
| `parameters/Factors` | 2 | Number of expression PEER factors to include as covariates in the eQTL analysis. |
| `parameters/permutations` | 1 | Number of permutations to perform for significance testing. Set to 1 to skip permutation testing. |
| `parameters/covariates` | ['Gender', 'Age'] | List of covariates to include in the eQTL analysis. This must match the column names in the covariates file. |
| `parameters/run_var_genes` | False | Whether to run the variable gene selection step. Set to True to select the top variable genes for eQTL analysis. |
| `parameters/n_var_genes` | NA | If `run_var_genes` is True, this specifies the number of top variable genes to select for eQTL analysis. Must be an integer greater than 0 and less than the total number of genes in the expression matrix. If `run_var_genes` is False, this should be set to NA. |
| `config/workdir` | /path/to/workdir/ | Path to the working directory where all intermediate and final results will be stored. |
| `config/peertool` | /path/to/peertool | Path to the PEER tool executable. This tool is used for inferring hidden confounding factors in gene expression data. |


> **__NOTE:__** Minimal example files for the input parameters can be found in the [example_files](https://github.com/vijaybioinfo/example_files/) directory of this repository.

## Example input files

### datasets CSV file

Datasets CSV file should have the following format:

| cell | tissue | seurat.obj.file | resolution | subset | expFile | donorFile | clusterID | cluster_def |
|------|--------|-----------------|------------|--------|---------|-----------|-----------|-------------|
| CD8 | LUNG | /path/to/seurat_object.RDS | 0.2 | 0 | /path/to/pseudobulk/profile/gene_expression.txt | /path/to/donors.txt | Naive CD8 T cells (TN) | Naive CD8 T cells (TN) |
| CD8 | LUNG | /path/to/seurat_object.RDS | 0.2 | 1 | /path/to/pseudobulk/profile/gene_expression.txt | /path/to/donors.txt | Central memory CD8 T cells (TCM) | Central memory CD8 T cells (TCM) |

- **expFile:** Gene expression file should be a tab-delimited text file with donors as columns and genes as rows. The first column should contain gene identifiers (e.g., Ensembl IDs or gene symbols). [expression files](https://github.com/vijaybioinfo/example_files/expression_files/).
- **donorFile:** Donor file should be list of donors in the same order as the columns in the gene expression file.

### Covariates file

Covariates file should be a tab-delimited text file with samples as rows and covariates as columns. The first row should contain the column names for the covariates (e.g., 'Gender', 'Age').

| id | SAMPLE1 | SAMPLE2 | SAMPLE3 | SAMPLE4 | SAMPLE5 | SAMPLE6 | SAMPLE7 | SAMPLE8 | SAMPLE9 | SAMPLE10 | SAMPLE11 | SAMPLE12 | SAMPLE13 | SAMPLE14 | SAMPLE15 | SAMPLE16 | SAMPLE17 | SAMPLE18 | SAMPLE19 | SAMPLE20 | SAMPLE21 | SAMPLE22 | SAMPLE23 | SAMPLE24 | SAMPLE25 | SAMPLE26 | SAMPLE27 | SAMPLE28 | SAMPLE29 | SAMPLE30 |
|----|---------|---------|---------|---------|---------|---------|---------|---------|---------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| Age  | 69 | 27 | 70 | 63 | 50 | 47 | 28 | 6838 | 36 | 25 | 64 | 45 | 30 | 50 | 24 | 5038 | 69 | 42 | 48 | 21 | 39 | 60 | 58 | 5261 | 37 | 45 | 65 | 33 | 49 | 59 |
| Gender | 0 | 1 | 1 | 1 | 1 | 0 | 1 | 101 | 1 | 0 | 0 | 1 | 0 | 0 | 0 | 100 | 0 | 0 | 1 | 1 | 0 | 0 | 0 | 101 | 1 | 1 | 0 | 0 | 0 | 1 |

### SNP files

SNP files contain three  columns representing chromosome, SNP identifier, and SNP location. SNP and SNP_LOC columns should be consistent with the VCF file used for genotypes. Both files should have matrixeQTL format. (SNP files)[https://github.com/vijaybioinfo/example_files/snp_files/]

The SNP files should be tab-delimited text files with the following format:

| CHR | SNP | SNP_LOC |
|-----|-----|---------|
| 20 | example_files/snp_files/snp/snp_20.txt | example_files/snp_files/snp_loc/snppos_20.txt |
| 21 | example_files/snp_files/snp/snp_21.txt | example_files/snp_files/snp_loc/snppos_21.txt |

### Gene annotation files

Gene annotation files should have similar matrixeQTL format with three columns representing gene identifier, chromosome, and gene location. The gene location should be specified relative to the transcription start site (TSS).

| geneid | chr | s1 | s2 |
|--------|-----|----|----|
| ENSG00000178591 20 | 87250 | 87251 |
| ENSG00000125788 20 | 142590 | 142591 |

### Gene annotation complete files

Gene annotation complete file should have extra columns with the full gene coordinates (start and end positions) relative to the reference genome. This file is used for final annotation of significant eQTLs with gene names.

### rsID tables

RsID tables should have columns representing chromosome, position, reference allele, alternate allele, ID used in the VCF file, and the rsID. These tables are used for final annotation of significant eQTLs with rsIDs. (rstables)[https://github.com/vijaybioinfo/example_files/rs_tables/]