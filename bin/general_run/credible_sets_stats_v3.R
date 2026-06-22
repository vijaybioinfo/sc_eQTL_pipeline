library(data.table)
library(susieR)
library(Matrix)

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("--expFile"), action="store", default=NULL, type='character',
              help="Initial gene expression file."),
  make_option(c("--donorFile"), action="store", default=NULL, type='character',
              help="File with donors for the analysis in the given gene expression."),
  make_option(c("--output"), action="store", default=NULL, type='character',
              help="Output file to save normalize gene expression."),
  make_option(c("--covariates"), action="store", default=NULL, type='character',
              help="Covariates file used for eQTL analysis."),
  make_option(c("--sfile"), action="store", default=NULL, type='character',
              help="File with summary eQTL statistics."),
  make_option(c("--associations"), action="store", default=NULL, type='character',
              help="General path to eQTL associations."),
  make_option(c("--genotypefile"), action="store", default=NULL, type='character',
              help="File with path to genotype files.")
)
opt = parse_args(OptionParser(option_list=option_list))
####
if( FALSE ){
  opt = list(
    expFile = "results/matrix_eqtl/covariates/rarecells/DC/2/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/GE.tab",
    donorFile = "/mnt/cephfs/vd-vijay/Cristian/DICE_GALAXY/eQTL/data/sc_gene_expression/DC/2/donors.txt",
    output = "results/matrix_eqtl/eQTL/rarecells/DC/2/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/credible_sets.txt",
    covariates = "results/matrix_eqtl/covariates/rarecells/DC/2/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/covariates.txt",
    sfile = "results/matrix_eqtl/eQTL/rarecells/DC/2/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA/summary_table.csv",
    associations = "results/matrix_eqtl/eQTL/rarecells/DC/2/MAF_0.05_covariatesGender-Age_genPCs4_expPEER6_VarGenesFalse_Ngenes_NA",
    genotypefile = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_GALAXY/eQTL/data/genotype/matrix_eQTL/snps_chrs_files.txt"
  )
}
####
egene_genotype <- function(egene, genotype, aFile, donors){
  egene_asso <- aFile[aFile$gene == egene,]
  rownames(egene_asso) <- 1:nrow(egene_asso)
  egene_asso <- as.data.table(egene_asso)
  setkey(egene_asso, snps)
  genotype <- genotype[ID %in% egene_asso$snps]
  setkey(genotype, ID)
  genotype <- genotype[,..donors]
  return(list(geno = data.table::transpose(genotype) , asso = egene_asso))
}
## covariates regress out: https://stephenslab.github.io/susieR/articles/finemapping.html#a-note-on-covariate-adjustment
remove.covariate.effects <- function (X, Z, y) {
  # include the intercept term
  if (any(Z[,1]!=1)) Z = cbind(1, Z)
  A   <- forceSymmetric(crossprod(Z))
  SZy <- as.vector(solve(A,c(y %*% Z), tol = 1e-30))
  SZX <- as.matrix(solve(A,t(Z) %*% X, tol = 1e-30))
  y <- y - c(Z %*% SZy)
  X <- X - Z %*% SZX
  return(list(X = X,y = y,SZy = SZy,SZX = SZX))
}

####
# expFile <- read.csv(opt$expFile, row.names = 1)
expFile <- read.table(opt$expFile, row.names = 1, header = TRUE, sep = "\t",
                      check.names = FALSE)
donorFile <- readLines(opt$donorFile)
donorFile <- intersect(donorFile, names(expFile))

sfile <- data.table::fread(opt$sfile)
sfile <- sfile[is.eGene == TRUE]
chromosomes <- sort(gsub("chr", "", unique(sfile$Chr)))
covariates <- read.table(opt$covariates, sep = "\t", row.names = 1,
                         check.names = FALSE)
covariates <- covariates[, donorFile]
covariates <- t(covariates)


expFile <- expFile[, donorFile]
expFile <- t(expFile)
n <- length(donorFile)

## load genotype
snpfiles <- data.table::fread(opt$genotypefile, header = TRUE)

all_cs <- lapply(chromosomes, function(chr){
# all_cs <- lapply(c(6), function(chr){
  cat("Analyzing chromosome:", chr, "\n")
  stmp <- sfile[Chr == paste0("chr", chr)]
  egenes <- stmp$Geneid
  ####
  genoFile <- data.table::fread(snpfiles[CHR == chr]$SNP, header = TRUE)
  assoFile <- readRDS(paste(opt$associations, chr, "/Output_all_cis.rds", sep = "/"))
  #### Correct genotype 'ID' column if needed
  if(!"ID" %in% colnames(genoFile)){
    if("V1" %in% colnames(genoFile)){
      setnames(genoFile, "V1", "ID")
    }else if("id" %in% colnames(genoFile)){
      setnames(genoFile, "id", "ID")
    } else {
      stop("No ID column found in genotype data")
    }
  }
  ####
  genes_sets <- lapply(egenes, function(gene){
    cat("\t", gene, "\n")
    gene_data <- egene_genotype(egene = gene, genotype = genoFile, aFile = assoFile, donors = donorFile)
    # extra check
    # SNPS tested do not have same genotype for all donors
    idxs <- sapply(gene_data$geno, function(col) length(unique(col)) > 1)
    gene_data$geno <- gene_data$geno[, idxs, with = FALSE]
    gene_data$asso <- gene_data$asso[idxs,]
    ## regress out covariates in genotype and expression data
    out <-  remove.covariate.effects(X = as.matrix(gene_data$geno), Z = covariates, y = expFile[, gene])
    ####
    gene_data$asso$beta_se = gene_data$asso$beta / gene_data$asso$statistic
    # if statistic is 0 it will triger a x/0 -> NA (not defined), change beta_se to small values
    gene_data$asso$statistic[which(gene_data$asso$statistic==0)] <- 0.00000000001
    # change beta_se values of 0 to close to 0
    gene_data$asso$beta_se[which(gene_data$asso$beta_se==0 | is.na(gene_data$asso$beta_se))] <- 0.00000000001
    # change beta_se values of 0 to close to 0
    gene_data$asso$beta[which(gene_data$asso$beta==0)] <- 0.00000000001
    z_scores <- gene_data$asso$beta / gene_data$asso$beta_se
    ####
    # fitted_rss <- tryCatch({ susie_rss( bhat = gene_data$asso$beta, shat = gene_data$asso$beta_se,
    #                                     R = cov2cor(crossprod(out$X)), n = n, L = 10,
    #                                     coverage = 0.95, min_abs_corr = 0.5,
    #                                     max_iter = 200, estimate_residual_variance = TRUE)
    #                         }, error = function(e) {
    #                           susie_rss( bhat = gene_data$asso$beta, shat = gene_data$asso$beta_se,
    #                                     R = cov2cor(crossprod(out$X)), n = n,
    #                                     L = 10, coverage = 0.95, min_abs_corr = 0.5,
    #                                     max_iter = 200, estimate_residual_variance = FALSE)
    #                                     })
    # Estimate CS for small clusters were really challenging.
    # For this reason we integrate multiple strategies below.
    # One common issue we could not avoid was when 'estimated prior variance is unreasonably large'.
    # Even after inspect the whole data. If this error persist, we just ignore this specific error.
    fitted_rss <- tryCatch({
                            # Attempt 1: out$X with residual variance
                            susie_rss(
                              bhat = gene_data$asso$beta,
                              shat = gene_data$asso$beta_se,
                              R = cov2cor(crossprod(out$X)),
                              n = n,
                              L = 10,
                              coverage = 0.95,
                              min_abs_corr = 0.5,
                              max_iter = 200,
                              estimate_residual_variance = TRUE
                            )
                          }, error = function(e1) {
                            # Attempt 2: out$X without residual variance
                            tryCatch({
                              susie_rss(
                                bhat = gene_data$asso$beta,
                                shat = gene_data$asso$beta_se,
                                R = cov2cor(crossprod(out$X)),
                                n = n,
                                L = 10,
                                coverage = 0.95,
                                min_abs_corr = 0.5,
                                max_iter = 200,
                                estimate_residual_variance = FALSE
                              )
                            }, error = function(e2) {
                              # Attempt 3: fallback to geno-derived LD
                              tryCatch({
                                susie_rss(
                                  bhat = gene_data$asso$beta,
                                  shat = gene_data$asso$beta_se,
                                  R = cov2cor(crossprod(as.matrix(gene_data$geno))),
                                  n = n,
                                  L = 10,
                                  coverage = 0.95,
                                  min_abs_corr = 0.5,
                                  max_iter = 200,
                                  estimate_residual_variance = FALSE
                                )
                              }, error = function(e3) {
                                msg2 <- conditionMessage(e3)
                                if (!grepl("estimated prior variance is unreasonably large", msg2)) {
                                  stop("[ERROR] SuSiE failed even after fallback")
                                }
                              })
                            })
                          })

    ## Check if previous function retrieve null due to large variance
    if(is.null(fitted_rss)){
      return()
    }

    ##
    cs <- susie_get_cs(fitted_rss, coverage=0.95, min_abs_corr=0.5, X = out$X, n =120)
    ## No credible sets detected
    if(is.null(cs$cs)){
      return()
    }
    best_pip_per_cs <- sapply(cs$cs, function(idx){
       max(fitted_rss$pip[idx])
    })
    ####
    dfout <- data.table::data.table(egene = gene, cs = paste0("CS_", 1:length(cs$cs)), n.variants = sapply(cs$cs, length), best.pip = best_pip_per_cs)
    ###
    # loop through all sets and retrieve variants on them
    variantsdf <- data.table::rbindlist(lapply(cs$cs, function(cset){
       data.table::data.table(egene = gene, rsID = gene_data$asso$snps[cset], pip = fitted_rss$pip[cset])
    }))
    len_cs <- unlist(lapply(cs$cs, length))

    variantsdf[, cs := paste0("CS_", rep(1:length(len_cs), len_cs))]
    return(list(dfout, variantsdf))
  })
  genes_csets <- data.table::rbindlist(lapply(genes_sets, '[[', 1))
  variants_cset <- data.table::rbindlist(lapply(genes_sets, '[[', 2))
  return(list(genes_csets, variants_cset))
})
dfout <- data.table::rbindlist(lapply(all_cs, '[[', 1))
cs_variants <- data.table::rbindlist(lapply(all_cs, '[[', 2))
#####

cat("Saving output files...")
write.csv(dfout, opt$output, row.names = FALSE)
write.csv(cs_variants, gsub(".txt", "_variants.txt", opt$output), row.names = FALSE)
