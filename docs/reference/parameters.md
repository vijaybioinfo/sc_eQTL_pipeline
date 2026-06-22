# Parameters

Full reference for all keys in `config/config.yaml`.

## `samples`
- **Type:** list of strings
- **Required:** yes
- **Example:** `["sample_01", "sample_02"]`

Sample IDs. Each ID must match a file named `{sample}.fastq.gz` in `data_dir`.

## `data_dir`
- **Type:** string
- **Default:** `data/raw`

Path to the directory containing raw FASTQ files, relative to the repo root.

## `normalization`
- **Type:** string — one of `TMM`, `RLE`, `upperquartile`
- **Default:** `TMM`

Normalization method passed to `edgeR::calcNormFactors()`.

## `min_counts`
- **Type:** integer
- **Default:** `10`

Minimum count value a gene must have in at least `min_samples` samples to be retained.

## `min_samples`
- **Type:** integer
- **Default:** `3`

Minimum number of samples in which a gene must pass the `min_counts` threshold.

## `results_dir`
- **Type:** string
- **Default:** `results`

Directory where all output files are written.
