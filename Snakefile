import pandas as pd
import os
import re

configfile: "snake_conf.yaml"

workdir: config['config']['workdir']
#report: "report/workflow.rst"

class Table:
    def __init__(self, file):
        self.data = pd.read_csv(file, dtype=str)
    def extract(self, wildcards):
        #sel = datasets[ datasets["cell"] == wildcards.cell ]
        sel = self.data[ self.data["cell"] == wildcards.cell ]
        sel = sel[ sel["subset"].astype(str) == wildcards.subset ]
        sel = sel[ sel["tissue"] == wildcards.tissue ]
        return sel.to_dict(orient="list")

    def donor(self, wildcards):
        dic = self.extract(wildcards)
        return {"donorFile": dic["donorFile"]}

    def expression(self, wildcards):
        dic = self.extract(wildcards)
        return {"expFile": dic["expFile"]}

    def frequency(self, wildcards):
        dic = self.extract(wildcards)
        return {"freqtable": dic["prop"]}

def snpfiles(wildcards):
    dic = SNP_FILES[ SNP_FILES['CHR'] == wildcards.chrom ].to_dict(orient = 'list')
    return {"snps": dic["SNP"], "snpsloc": dic["SNP_LOC"]}


#### Load dataset file
dataset_file = config["inpFiles"]["datasets"]
datasets = Table(dataset_file)

TISSUES = datasets.data.to_dict(orient='list')['tissue']
CELLS = datasets.data.to_dict(orient='list')['cell']
SUBSETS = datasets.data.to_dict(orient='list')['subset']
PERMUTATIONS = int(config["parameters"]["permutations"])
#######
SNP_FILES = pd.read_table(config["inpFiles"]["snpfiles"])
CHROMS =  SNP_FILES['CHR'].to_list()
###
PEER = config["config"]["peertool"]
###
covariates = config['parameters']['covariates']
covariates = '-'.join(covariates)

run_var_genes = config['parameters']['run_var_genes']
n_var_genes = config['parameters']['n_var_genes']
print(run_var_genes)
print(n_var_genes)
mst_var_gns = '_VarGenes' + str(run_var_genes) + '_Ngenes_' + str(n_var_genes)
print(mst_var_gns)

if run_var_genes:
    if type(n_var_genes) != int or n_var_genes < 1:
        sys.exit("ERROR: if 'run_var_genes' == True, 'n_var_genes' int > 0 & < number genes in the matrix")


GPREFIX = "MAF_" + str(config['parameters']['MAF']) + "_covariates" + covariates + "_genPCs" + str(config['parameters']['PCS']) + "_expPEER" + str(config['parameters']['Factors']) + mst_var_gns
######
rule all:
    input:
        expand("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/credible_sets.txt", zip, tissue = TISSUES, cell = CELLS, subset = SUBSETS),
        expand("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/credible_sets_variants.txt", zip, tissue = TISSUES, cell = CELLS, subset = SUBSETS),
        expand("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/eQTL_sigpairs.csv", zip, tissue = TISSUES, cell = CELLS, subset = SUBSETS),
        expand("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/all_pairs_adjust_pvalue.rds", zip, tissue = TISSUES, cell = CELLS, subset = SUBSETS),
        expand("results/fastQTL/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_all.genes.txt.gz", zip, tissue = TISSUES, cell = CELLS, subset = SUBSETS),
        expand("results/fastQTL/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_all.allpairs.txt.gz", zip, tissue = TISSUES, cell = CELLS, subset = SUBSETS)

rule expression_preprocess:
    input:
        unpack(datasets.donor),
        unpack(datasets.expression),
        annot_file = config["inpFiles"]["geneannot"]
    output:
        exp_file = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.tab",
        fast_exp = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.fastqtl.txt.gz",
        fast_exp_idx = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.fastqtl.txt.gz.tbi"
    params:
        script = "bin/general_run/process_expression.R",
        fast_exp = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.fastqtl.txt"
    conda:
        "DLCP_v3"
    shell:
        """
        Rscript {params.script} -e {input.expFile} -d {input.donorFile} -q -o {output.exp_file} --fastqtl --annotFile {input.annot_file} 
        bgzip {params.fast_exp}
        tabix -f -p bed {output.fast_exp}
        """

rule peer_calc:
    input:
        expression = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.tab"
    output:
        "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/PEERanalysis/factors.txt"
    params:
        script = PEER,
        peerfactors_n = config["parameters"]["Factors"],
        odir = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/PEERanalysis/",
        factors = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/PEERanalysis/X.csv",
        run_var_genes = run_var_genes,
        n_genes = n_var_genes
    run:
        print(int(params.peerfactors_n))
        if  int(params.peerfactors_n) < 1:
            shell("{params.script} -f {input.expression}  --has_header --has_rownames -o {params.odir} -n 1 --transpose;")
            shell("cat <(seq 1 {params.peerfactors_n}|sed 's/^/Factor/'|sed -z 's/\\n/,/'|sed 's/^/,/' ) <(paste -d ',' <(head {input.expression} -n1|sed 's/\\t/\\n/g') {params.factors}) > {output}")
        else:
            shell("{params.script} -f {input.expression}  --has_header --has_rownames -o {params.odir} -n {params.peerfactors_n} --transpose;")
            shell("cat <(seq 1 {params.peerfactors_n}|sed 's/^/Factor/'|sed -z 's/\\n/,/g'|sed 's/^/,/'|sed 's/$/\\n/'  ) <(paste -d ',' <(head {input.expression} -n1|sed 's/\\t/\\n/g') {params.factors}) > {output}")

rule subset_geno:
    input:
        unpack(datasets.donor),
        genotype = config["inpFiles"]["vcf_genotype"]
    output:
        temp("results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.vcf.gz")
    shell:
        "bcftools view -S {input.donorFile} -Oz -o {output}  {input.genotype}"

rule LD_prune:
    input:
        vcf = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.vcf.gz"
    output:
        prunein = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.prune.in",
        pruneout = temp("results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.prune.out")
    params:
        "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype"
    shell:
        "plink --vcf {input.vcf} --indep-pairwise 200 100 0.1 --out {params}"

rule pca:
    input:
        genotype = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.vcf.gz",
        variants = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.prune.in"
    output:
        eigenvec = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.eigenvec",
        eigenval = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.eigenval",
        glog = temp("results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.log"),
        gnsx = temp("results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.nosex")
    params:
        "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype"
    shell:
        "plink --vcf {input.genotype} --extract {input.variants} --pca header tabs --out {params}"

rule merge_cov:
    input:
        covariates = config["inpFiles"]["covariates"],
        peer_file = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/PEERanalysis/factors.txt",
        eigenvec = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.eigenvec"
    output:
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates.txt"
    params:
        script = "bin/general_run/merge_covariates.R",
        pcs = config["parameters"]["PCS"],
        factors = config["parameters"]["Factors"],
        covar = covariates
    run:
        if not params.covar:
            shell("Rscript {params.script} --pcafile {input.eigenvec} --factorfile {input.peer_file} --covariatesFile {input.covariates} --factors {params.factors} --pcs {params.pcs} --output {output.covariates} ")
        else:
            shell("Rscript {params.script} --pcafile {input.eigenvec} --factorfile {input.peer_file} --covariatesFile {input.covariates} --covariates {params.covar} --factors {params.factors} --pcs {params.pcs} --output {output.covariates} ")

rule matrix_eQTL:
    input:
        # snps = config["inpFiles"]["snps"],
        # snpsloc = config["inpFiles"]["snpsloc"],
        unpack(snpfiles),
        expression = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.tab",
        geneloc = config["inpFiles"]["geneannot"],
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates.txt"
    output:
        cis = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_cis.txt",
        trans = temp("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_tra.txt"),
        dgfree = temp("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_df.txt"),
        cis_all = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_all_cis.txt",
        cis_all_rds = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_all_cis.rds",
        qqplot = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_qqplot.png"
    params:
        script = "bin/general_run/Matrix_eQTL.R",
        null = "FALSE",
        MAF = config["parameters"]["MAF"],
        outputdir = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/",
        prefix = "Output"
    conda:
        "DLCP_v3"
    shell:
        "Rscript {params.script}  --nullDist {params.null} --snpfile {input.snps} \
         --snplocation {input.snpsloc} --expressionfile {input.expression} \
         --covariates {input.covariates} --genelocation {input.geneloc} --MAF {params.MAF} \
         --output {params.prefix} --outdir {params.outputdir} "

rule eigenMT:
    input:
        unpack(datasets.donor),
        unpack(snpfiles),
        qtl = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_all_cis.txt",
        geneloc = config["inpFiles"]["geneannot"]
    output:
        eigen = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_eigen_cis_sig.tsv"
    params:
        script = "bin/general_run/eigenMT.py",
        chromosome = "{chrom}"
    conda:
        "pyEigenMT"
    shell:
        "python {params.script} --QTL {input.qtl} --GEN {input.snps} --GENPOS {input.snpsloc} --PHEPOS {input.geneloc} --OUT {output.eigen} --sample_list {input.donorFile} --CHROM {params.chromosome} "

rule null_distribution:
    input:
        # snps = config["inpFiles"]["snps"],
        # snpsloc = config["inpFiles"]["snpsloc"],
        unpack(snpfiles),
        expression = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.tab",
        geneloc = config["inpFiles"]["geneannot"],
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates.txt"
    output:
        cis = temp("temp/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_{permutation}_cis.txt"),
        qqplot = temp("temp/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/Output_{permutation}_qqplot.png")
    params:
        script = "bin/general_run/Matrix_eQTL.R",
        null = "TRUE",
        MAF = config["parameters"]["MAF"],
        outputdir = "temp/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/{chrom}/",
        prefix = "Output_{permutation}"
    conda:
        "DLCP_v3"
    shell:
        "Rscript {params.script}  --nullDist {params.null} --snpfile {input.snps} "
        " --snplocation {input.snpsloc} --expressionfile {input.expression} "
        " --covariates {input.covariates} --genelocation {input.geneloc} --MAF {params.MAF} "
        " --output {params.prefix} --outdir {params.outputdir} "


rule merge_null:
    input:
        null = expand("temp/matrix_eqtl/null_distr/{{tissue}}/{{cell}}/{{subset}}/" + GPREFIX + "/{chrom}/Output_{permutation}_cis.txt", chrom = CHROMS, permutation = list(range(1,PERMUTATIONS + 1)))
    output:
        null = temp("results/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/null_distr.tsv")
    params:
        Rfile = "bin/general_run/Merge_Results_null.R"
    conda:
        "DLCP_v3"
    shell:
        "Rscript {params.Rfile} {input.null} {output.null}"


rule calculate_FDR:
    input:
        cis = expand("results/matrix_eqtl/eQTL/{{tissue}}/{{cell}}/{{subset}}/" + GPREFIX + "/{chrom}/Output_cis.txt", chrom = CHROMS),
        tra = expand("results/matrix_eqtl/eQTL/{{tissue}}/{{cell}}/{{subset}}/" + GPREFIX + "/{chrom}/Output_tra.txt", chrom = CHROMS),
        null = "results/matrix_eqtl/null_distr/{tissue}/{cell}/{subset}/" + GPREFIX + "/null_distr.tsv",
        eigen = expand("results/matrix_eqtl/eQTL/{{tissue}}/{{cell}}/{{subset}}/" + GPREFIX + "/{chrom}/Output_eigen_cis_sig.tsv", chrom = CHROMS)
    output:
        cis = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_All_cis_sig.tsv",
        tra = temp("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_All_tra_sig.tsv"),
        qqcis = temp("results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_All_cis_sig_QQplot.png")
    params:
        Rfile = "bin/general_run/Calculate_FDR_merge.R",
        prefix = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output",
        cis_results = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/",
        cis_name = "Output_cis.txt",
        tra_results = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/",
        tra_name = "Output_tra.txt",
        eigen_results = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/",
        eigen_name = "Output_eigen_cis_sig.tsv"
    conda:
        "DLCP_v3"
    shell:
        "Rscript {params.Rfile} -c {params.cis_results} -t {params.tra_results} -i {params.cis_name} -r {params.tra_name} -n {input.null} -o {params.prefix} -m {params.eigen_results} -e {params.eigen_name} "

rule summary_table:
    input:
        unpack(datasets.donor),
        unpack(datasets.expression),
        unpack(datasets.frequency),
        geneannot = config["inpFiles"]["geneannot_complete"],
        snptable = config["inpFiles"]["snpfiles"],
        sigFile = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_All_cis_sig.tsv"
    output:
        summaryfile = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/summary_table.csv"
    params:
        script = "bin/general_run/summary_table.R",
        rstables = config["inpFiles"]["rstables"]
    shell:
        """
        Rscript {params.script} --sigFile {input.sigFile} \
            --geneannot {input.geneannot} \
            --snpTable {input.snptable} \
            --rsdatabase {params.rstables} \
            --expFile {input.expFile} \
            --donorFile {input.donorFile} \
            --freqtable {input.freqtable} \
            --outputfile {output.summaryfile}
        """

rule sigpairs:
    input:
        summaryfile = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/summary_table.csv",
    output:
        sigpairs = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/eQTL_sigpairs.csv",
        allpairs = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/all_pairs_adjust_pvalue.rds"
    params:
        script = "bin/general_run/sig_pairs_adj_pvalues.R",
        assocpath = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/"
    shell:
        """
        Rscript {params.script} --summfile {input.summaryfile} \
            --assocpath {params.assocpath} \
            --sigoutput {output.sigpairs} \
            --allrdsoutput {output.allpairs}
        """

rule finemapping:
    input:
        unpack(datasets.donor),
        exp_file = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.tab",
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates.txt",
        summaryfile = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/summary_table.csv",
        genotypefile = config["inpFiles"]["snpfiles"]
    output:
        credible_sets = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/credible_sets.txt",
        finemap_variants = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/credible_sets_variants.txt"
    params:
        script = "bin/general_run/credible_sets_stats_v3.R",
        associations = "results/matrix_eqtl/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX
    conda:
        "DLCP_v3"
    resources:
        walltime = 400,
        mem_gb = 80
    shell:
        """
        Rscript {params.script} --expFile {input.exp_file} \
            --donorFile {input.donorFile} \
            --output {output.credible_sets} \
            --covariates {input.covariates} \
            --sfile {input.summaryfile} \
            --associations {params.associations} \
            --genotypefile {input.genotypefile}
        """

# ## Fastqtl rules:
rule merge_cov_fasqtl:
    input:
        covariates = config["inpFiles"]["covariates"],
        peer_file = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/PEERanalysis/factors.txt",
        eigenvec = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/genotype.eigenvec"
    output:
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates_fastqtl.txt.gz"
    params:
        script = "bin/general_run/merge_covariates.R",
        pcs = config["parameters"]["PCS"],
        factors = config["parameters"]["Factors"],
        covar = covariates,
        covarfile = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates_fastqtl.txt"
    run:
        if not params.covar:
            shell("Rscript {params.script} --pcafile {input.eigenvec} --factorfile {input.peer_file} --covariatesFile {input.covariates} --factors {params.factors} --pcs {params.pcs} --output {params.covarfile} --fastqtl")
        else:
            shell("Rscript {params.script} --pcafile {input.eigenvec} --factorfile {input.peer_file} --covariatesFile {input.covariates} --covariates {params.covar} --factors {params.factors} --pcs {params.pcs} --output {params.covarfile} --fastqtl ")


rule fastqtl_nom:
    input:
        genotype = config["inpFiles"]["vcf_genotype"],
        exp_file = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.fastqtl.txt.gz",
        exp_file_idx = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.fastqtl.txt.gz.tbi",
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates_fastqtl.txt.gz"
    output:
        qtl = "results/fastQTL/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_all.allpairs.txt.gz"
    params:
        maf = config["parameters"]["MAF"],
        chunks = config["fastqtl_parameters"]["chunks"],
        ma_sample = config["fastqtl_parameters"]["ma_sample"],
        cis_window = config["fastqtl_parameters"]["cis_window"],
        out_prefix = "results/fastQTL/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_all"
    threads: 10
    shell:
        "/home/cgonzalez/tools/fastqtl/python/run_FastQTL_threaded.py {input.genotype} {input.exp_file} {params.out_prefix} --covariates  {input.covariates} --window {params.cis_window} --ma_sample_threshold {params.ma_sample} --maf_threshold {params.maf} --chunks {params.chunks} --threads {threads} "

rule fastqtl_perm:
    input:
        genotype = config["inpFiles"]["vcf_genotype"],
        exp_file = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.fastqtl.txt.gz",
        exp_file_idx = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/GE.fastqtl.txt.gz.tbi",
        covariates = "results/matrix_eqtl/covariates/{tissue}/{cell}/{subset}/" + GPREFIX + "/covariates_fastqtl.txt.gz",
	    cis = "results/fastQTL/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_all.allpairs.txt.gz"
    output:
        qtl = "results/fastQTL/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_all.genes.txt.gz"
    params:
        maf = config["parameters"]["MAF"],
        chunks = config["fastqtl_parameters"]["chunks"],
        ma_sample = config["fastqtl_parameters"]["ma_sample"],
        cis_window = config["fastqtl_parameters"]["cis_window"],
        permute_down = config["fastqtl_parameters"]["permute_down"],
        permute_up = config["fastqtl_parameters"]["permute_up"],
        out_prefix = "results/fastQTL/eQTL/{tissue}/{cell}/{subset}/" + GPREFIX + "/Output_all"
    threads: 20
    resources:
        walltime=540
    shell:
        "/home/cgonzalez/tools/fastqtl/python/run_FastQTL_threaded.py {input.genotype} {input.exp_file} {params.out_prefix} --covariates {input.covariates} --permute {params.permute_down} {params.permute_up} --window {params.cis_window} --ma_sample_threshold {params.ma_sample} --maf_threshold {params.maf} --chunks {params.chunks} --threads {threads} "


