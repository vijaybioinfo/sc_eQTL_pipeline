###############################################################################
## Script to get summary table for each dataset: eGenes, summary stats, mean
## expression per genotype group, proportion of expressing cells, etc.
## Author: Cristian Gonzalez-Colin
###############################################################################
# Load libraries
library(data.table)
library(Seurat)
library(arrow)
library(tidyverse)

########
# Get command line arguments
suppressPackageStartupMessages(require(optparse))

option_list = list(
  optparse::make_option(c("--sigFile"), type="character", default=NULL, help="Path to the file with significant associations from matrix eQTL", metavar="character"),
  optparse::make_option(c("--geneannot"), type="character", default=NULL, help="Path to the file with gene annotation", metavar="character"),
  optparse::make_option(c("--snpTable"), type="character", default=NULL, help="Path to the file with the list of genotype files for each chromosome", metavar="character"),
  optparse::make_option(c("--expFile"), type="character", default=NULL, help="Path to the file with gene expression data", metavar="character"),
  optparse::make_option(c("--donorFile"), type="character", default=NULL, help="Path to the file with the list of donors in the dataset", metavar="character"),
  optparse::make_option(c("--freqtable"), type="character", default=NULL, help="Path to the file with the frequency of expressing cells per gene and donor", metavar="character"),
  optparse::make_option(c("--rsdatabase"), type="character", default=NULL, help="Path to the database with rsnumbers for each variant", metavar="character"),
  optparse::make_option(c("--outputfile"), type="character", default=NULL, help="Path to the output file where the summary table will be saved", metavar="character")
  )
opt = parse_args(OptionParser(option_list=option_list))

## function to get the mean expression per genotype group and the mean
## frequency of expressing cells per genotype group for each SNP-gene pair
get_metrics <- function(expdf, geno, donors, gs_pair, freqCounts){
  ## get gene and snp names from the pair
  genes <- gsub(",.+$", "", gs_pair)
  snps <- gsub("^[^.]+,", "", gs_pair)
  expdf <- expdf[genes,]
  freqCounts <- freqCounts[gene %in% genes]
  dfmean <- data.table(rowMeans(expdf))
  colnames(dfmean) <- "Mean_Exp"
  dfmean$gene <- rownames(expdf)
  ######
  expdf$gene <- rownames(expdf)
  expdf <- data.table(expdf)
  #####
  dfpair <- data.table(gene = genes, snp = snps, gene_snp = gs_pair)
  #####
  expdf <- melt(expdf, id.vars = "gene")
  setnames(expdf, "value", "expression")
  ##expression dataframe with SNP_gene pairs
  dfpair_exp <- merge(dfpair, expdf)
  dfpair_exp[, snpdonor := paste(snp, variable, sep = ",")]
  ###freq dataframe with SNP_gene pairs
  freqCounts <- melt(freqCounts, id.vars = "gene")
  setnames(freqCounts, "value", "frequency")
  # merge with the SNP-gene pairs
  dfpair_frq <- merge(dfpair, freqCounts)
  dfpair_frq[, snpdonor := paste(snp, variable, sep = ",")]

  ## get mean expression per genotype group for each SNP-gene pair
  dfmeanAlleles <- geno %>%  filter(ID %in% snps) %>%
         select(ID, donors) %>%
         pivot_longer(cols = -ID, names_to = "donors", values_to = "genotype") %>%
         mutate(snpdonor = paste(ID, donors, sep = ",")) %>%
         right_join(dfpair_exp, by = c("snpdonor" = "snpdonor")) %>%
         group_by(gene_snp, genotype) %>%
         summarise(meanGroup = mean(expression)) %>%
         ungroup() %>%
         pivot_wider(names_from = "genotype", values_from = "meanGroup", names_prefix = "Mean_Allele") %>%
         mutate(gene = gsub(",.+$", "", gene_snp), snp = gsub("^[^.]+,", "", gene_snp))

  ## get mean frequency of expressing cells per genotype group for each SNP-gene pair
  dfFreqAlleles <- geno %>%  filter(ID %in% snps) %>%
         select(ID, donors) %>%
         pivot_longer(cols = -ID, names_to = "donors", values_to = "genotype") %>%
         mutate(snpdonor = paste(ID, donors, sep = ",")) %>%
         right_join(dfpair_frq, by = c("snpdonor" = "snpdonor")) %>%
         group_by(gene_snp, genotype) %>% 
         summarise(meanGroup = mean(frequency)) %>%
         ungroup() %>%
         pivot_wider(names_from = "genotype", values_from = "meanGroup", names_prefix = "Freq_Allele") %>%
         mutate(gene = gsub(",.+$", "", gene_snp), snp = gsub("^[^.]+,", "", gene_snp)) %>%
         select(gene_snp, Freq_Allele0, Freq_Allele1, Freq_Allele2)

  ##### merge everything in one table
  dfall <- merge( dfmean, dfmeanAlleles, by = "gene")
  dfall <- merge(dfall, dfFreqAlleles, by = "gene_snp")
  dfall$gene <- NULL
  dfall$snp <- NULL
  return(dfall)
}
########
## DEBUG
if (FALSE){
  opt <- list(sigFile = "results/matrix_eqtl/eQTL/rarecells/DC/0/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/Output_All_cis_sig.tsv",
              geneannot = "/mnt/bioadhoc/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/data/geneAnnotation/DLCP_geneAnnot_hg38.complete.txt",
              snpTable = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_GALAXY/eQTL/data/genotype/matrix_eQTL/snps_chrs_files.txt", 
              expFile = "/mnt/cephfs/vd-vijay/Cristian/DICE_GALAXY/eQTL/data/sc_gene_expression/DC/0/gene_expression.txt",
              donorFile = "/mnt/cephfs/vd-vijay/Cristian/DICE_GALAXY/eQTL/data/sc_gene_expression/DC/0/donors.txt",
              freqtable = "/mnt/cephfs/vd-vijay/Cristian/DICE_GALAXY/eQTL/data/sc_gene_expression/DC/0/freq_cell.txt",
              rsdatabase = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/genomes/dbsnp/GRCh38/bed/extendIDs/",
              outputfile = "results/matrix_eqtl/eQTL/rarecells/DC/0/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/summary_table.csv"
  )
}
## MAIN
####
# Load data
gene_info <- data.table::fread(opt$geneannot)

# load genotype data
snpTable <- data.table::fread(opt$snpTable)
genotype <- data.table::rbindlist(lapply(snpTable$SNP, function(fileid){
  iFile <- data.table::fread(fileid, header = TRUE)
}))

if(!"ID" %in% colnames(genotype)){
  if("V1" %in% colnames(genotype)){
    setnames(genotype, "V1", "ID")
  }else if("id" %in% colnames(genotype)){
    setnames(genotype, "id", "ID")
  } else {
    stop("No ID column found in genotype data")
  }
}
# load dbsnp database with rsnumbers
ds <- open_dataset(opt$rsdatabase, format = "tsv")
###########
cat("\tGet proportion data\n")
freqtable <- data.table::fread(opt$freqtable, header = TRUE)
# get proportion of donors expressin x% of cells for each gene
propdata <- data.table::data.table(gene = freqtable$gene,
                                    prop.0.01 = rowSums(freqtable[,-1] > 0.01) / (ncol(freqtable)-1),
                                    prop.0.05 = rowSums(freqtable[,-1] > 0.05) / (ncol(freqtable)-1),
                                    prop.0.10 = rowSums(freqtable[,-1] > 0.10) / (ncol(freqtable)-1)
)

###
cat("\tLoading eQTL data\n")
eqtlData <- data.table::fread(opt$sigFile)
eqtlData$chr <- NULL

eqtlData[, neg.log.10.FDR := -log10(BF.FDR)]
######
allsnps <- unique( eqtlData$snps )
QTLpairs <- unique( paste( eqtlData$gene, eqtlData$snps, sep = "," ) )
###
cat("\tGet mean expression per group")
expFile <- read.table(opt$expFile, check.names = FALSE, row.names = 1, header = TRUE)
donFile <- readLines(opt$donorFile)

dfmeangroup <- get_metrics(expFile, genotype, donFile, QTLpairs, freqtable)
###
cat("\tGetting rsnumbers\n")
snpdf <- ds %>% select(ID, variant2, REF, ALT)  %>%
    filter(variant2 %in% allsnps) %>% collect()

snpdf <- snpdf %>% group_by(variant2) %>%
        mutate(ID = paste(ID, collapse = "|")) %>% unique()
################
cat("\tMerging data\n")
dfall <- merge(eqtlData, propdata, by = "gene", all.y = T)
#######
dfall[, gene_snp := paste(gene, snps, sep = ",")]
dfall <- merge(dfmeangroup, dfall, by = "gene_snp", all.y = T)
dfall <- merge( snpdf, dfall, by.y = "snps", by.x = "variant2", all.y = T )
dfall <- merge(gene_info, dfall, by.x = "Geneid", by.y = "gene")
#######remove genes in MT and ribosomal
dfall <- dfall[Chr %in% paste0( "chr", c(1:22,'X') ) ]
dfall <- dfall[ grep("^MT-|^RP", dfall$gene_name, invert = T)]
dfall$Mean_AlleleNA <- NULL
######
if("pvalue.x" %in% names(dfall)){
  dfall$pvalue.y <- NULL
  setnames(dfall, "pvalue.x", "pvalue")
}
dfall[, is.eGene_noMAC := prop.0.01 > 0.05 & BF.FDR < 0.05 & pvalue < 1e-4]

## add MAC threshold
# remove ID column to calculate allele count
# subet genotype to the donors in the dataset
allelecount <-  rowSums(genotype[, -1][, ..donFile], na.rm = TRUE)
minor_allele_count <- pmin(allelecount, 2*(length(donFile)) - allelecount)

dfall$MAC <- minor_allele_count[match(dfall$variant2, genotype$ID)]
dfall[, is.eGene := prop.0.01 > 0.05 & BF.FDR < 0.05 & pvalue < 1e-4 & MAC >= 10]
######
cat("\tSaving data\n")
write.csv(dfall, opt$outputfile)

