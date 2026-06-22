library(data.table)

suppressPackageStartupMessages(require(optparse))

option_list <- list(
  optparse::make_option(c("--summfile"), type="character", default = NULL,
                        help="Path to the file with significant associations from matrix eQTL",
                        metavar = "character"),
  optparse::make_option(c("--assocpath"), type="character", default = NULL,
                        help="Path to the folder with the association files for each chromosome",
                        metavar = "character"),
  optparse::make_option(c("--sigoutput"), type="character", default = NULL,
                        help="Path to the output file where the significant eQTLs will be saved",
                        metavar = "character"),
  optparse::make_option(c("--allrdsoutput"), type="character", default = NULL,
                        help="Path to the output file where all associations will be saved",
                        metavar = "character")
)
opt = parse_args(OptionParser(option_list=option_list))

## DEBUG
if(FALSE){
  opt = list(
    summfile = "results/matrix_eqtl/eQTL/rarecells/rarecells/DC/3/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/summary_table.csv",
    assocpath = "results/matrix_eqtl/eQTL/rarecells/rarecells/DC/3/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/",
    sigoutput = "results/matrix_eqtl/eQTL/rarecells/rarecells/DC/3/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/eQTL_sigpairs.csv",
    allrdsoutput = "results/matrix_eqtl/eQTL/rarecells/rarecells/DC/3/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/all_pairs_adjust_pvalue.rds"
  )
}

summfile <- data.table::fread(opt$summfile, header = TRUE)
# load association files to include the N.TESTs column in all associations
# (not only the significant ones)
chrvec <- c(1:22, "X")

alldata <- lapply(chrvec, function(chr){
  # load all cis associations for the chromosome
  afile <- readRDS(file.path(opt$assocpath, chr, "/Output_all_cis.rds"))
  afile <- data.table::as.data.table(afile)
  # load eigenMT file for the chromosome to get the number of tests per gene
  efile <- data.table::fread(file.path(opt$assocpath, chr,
                                       "/Output_eigen_cis_sig.tsv"),
                             header = TRUE)
  ## add number of tests to the association file
  mfile <- merge(afile, efile[, .(gene, TESTS)], by.x = "gene",
                 by.y = "gene", all.x = TRUE)

  mfile[, pvalue_genelevel_corrected := pvalue * TESTS]
  mfile[, pvalue_genelevel_corrected := ifelse(pvalue_genelevel_corrected > 1,
                                               1, pvalue_genelevel_corrected)]

  return(mfile)
})
alldata <- data.table::rbindlist(alldata)
# filter for significant associations at gene level
eqtlsig <- alldata[pvalue_genelevel_corrected < 0.05]
eqtlsig <- eqtlsig[gene %in% summfile[is.eGene == TRUE]$Geneid]
eqtlsig[, FDR := NULL]
eqtlsig[, se := beta/statistic]
eqtlsig[, c("chr", "snp.pos", "REF", "ALT") := tstrsplit(snps, ":")]
# save results
write.csv(eqtlsig, file = opt$sigoutput, row.names = FALSE)
saveRDS(alldata, file = opt$allrdsoutput)