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
