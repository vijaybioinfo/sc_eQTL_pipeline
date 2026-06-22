library(Matrix)
library(Seurat)
library(data.table)
#############

getsummarystats <- function(obj_srt, identifier = "donor.id.tag"){
  #####Number of cells per donor
  dfdonors <- data.table::data.table(table(obj_srt@meta.data[, identifier ]))
  names(dfdonors) <- c("donor", "n.cells")
  #remove donors with less than 10 cells
  dfdonors <- dfdonors[dfdonors$n.cells >= 10,]
  ###########filter donors from above
  if( nrow(dfdonors) != length( unique(obj_srt@meta.data[, identifier ]) ) | any(unique(inpObj@meta.data$full.donor.id.tag) == "Doublet")){
    cat("Number of donors", nrow(dfdonors), "\n")
    obj_srt <- obj_srt[, obj_srt@meta.data[[identifier]] %in% dfdonors[donor != "Doublet"]$donor]
  }
  #####summary donors
  ## Mean
  Seurat::Idents(obj_srt) <- identifier
  df.cell <- Seurat::AverageExpression(obj_srt)
  dfmean <- df.cell$RNA

  ## Freq/prop cells
  donor_filter <- lapply(unique( unlist(obj_srt@meta.data[identifier]) ), function(donor.id){
    #donor.id <-unique( unlist(obj_srt@meta.data[identifier]) )[1]
    cat("Looking at donorID: ", donor.id, "\n")
    cells.to.keep <- Cells(obj_srt)[obj_srt@meta.data[, identifier ] == donor.id]
    t.obj <- subset(obj_srt, cells=cells.to.keep)
    ##proportion of expressing cells
    df.out <- data.frame(Matrix::rowSums(GetAssayData(object = t.obj) > 0) / ncol(t.obj))
    names(df.out) <- donor.id
    df.out$gene <- rownames(df.out)
    df.out <- data.table::data.table(df.out)
    return(df.out)
  })
  dffreq <- Reduce(function(x,y) merge(x,y, by = "gene", all = T), donor_filter)

  return(list(mean = dfmean, freq = dffreq, cellsxdonor = dfdonors))
}

############
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-s", "--seuratobj"), action="store", default=NA, type='character',
              help="Seurat object to process."),
  make_option(c("-d", "--donorannot"), action="store", default=NA, type='character',
              help="Donor annot object to get final id."),
  make_option(c("-c", "--cell"), action="store", default=NA, type='character',
              help="Cell type to analyse."),
  make_option(c("-r", "--resolutioncolumn"), action="store", default=NA, type='character',
              help="Resolution column to filter seurat object."),
  make_option(c("-o", "--outputpath"), action="store", default=NA, type='character',
              help="Output path to save filter results.")
)
opt = parse_args(OptionParser(option_list=option_list))
### DEBUG
if(FALSE){
  opt = list( seuratobj = "/mnt/hpcscratch/vfajardo/R24/seurat_analysis/R24_Cancer_NK/R24_Cancer_Batches-1-to-20_NK_ALL_CD16neg-2/R24_Cancer_Batches-1-to-20_NK_ALL_CD16neg-2_05-09-2024_qc-nk-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/seurat_objects/SeuratObjectForPrjR24_Cancer_Batches-1-to-20_NK_ALL_CD16neg-2_WithArgs_NoPCs_30.RDS",
              donorannot = "/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/paper_items/labels_for_transfer/CD8_v1.csv",
              cell = "NKneg",
              resolutioncolumn = "RNA_snn_res.0.2",
              outputpath = "/mnt/BioAdHoc/Groups/vd-vijay/Cristian/DICE_LungCancer/tissue_tumor_combine/data/sc_gene_expression/"
  )
}
###
inpObj <- base::readRDS(opt$seuratobj)
metadata <- inpObj@meta.data
metadata$barcode <- rownames(metadata)
donorannot <- read.csv(opt$donorannot)
# Add pre leading zeros on column full.donor.id.tag if needed (4 digits)
donorannot$full.donor.id.tag <- paste0("DLCP", sprintf("%04d", donorannot$full.donor.id.tag))
donorannot$full.donor.id.tag[donorannot$full.donor.id.tag == "DLCP  NA"] <- NA

newmetadata <- merge(metadata, donorannot[, c("population.tag", "full.donor.id.tag", "barcode")], by = "barcode", all.x = TRUE)
rownames(newmetadata) <- newmetadata$barcode
newmetadata$barcode <- NULL
newmetadata <- newmetadata[rownames(inpObj@meta.data),]
inpObj@meta.data <- newmetadata
##### cluster data object
clusters <- levels(inpObj@meta.data[[opt$resolutioncolumn]])
# clusters <- c(0,2:8,10)

lapply(clusters, function(cluster.use){
  # cluster.use <- 0
  cat(opt$cell, " cluster: ", cluster.use, "\n")
  ##
  if( length(grep("_", cluster.use)) > 0 ){
    clusterObj <- inpObj[, inpObj@meta.data[[opt$resolutioncolumn]] %in% strsplit(cluster.use, "_")[[1]]]
  }else{
    clusterObj <- inpObj[, inpObj@meta.data[[opt$resolutioncolumn]] == cluster.use]
  }
  #####
  clusterstats <- getsummarystats(clusterObj, identifier = "full.donor.id.tag")
  clusterdonors <- as.vector(gsub("g", "", colnames(clusterstats$mean)))
  meandf <- as.data.frame(clusterstats$mean)
  colnames(meandf) <- clusterdonors
  dirpath <- paste(opt$outputpath, opt$cell, cluster.use, sep = "/")
  if(!dir.exists(dirpath)){
    dir.create(dirpath, recursive = TRUE)
  }
  data.table::fwrite(meandf, paste(opt$outputpath, opt$cell, cluster.use, "gene_expression.txt", sep = "/"), quote = FALSE, sep = "\t",  row.names = TRUE)
  data.table::fwrite(clusterstats$freq, paste(opt$outputpath, opt$cell, cluster.use, "freq_cell.txt", sep = "/"), quote = FALSE, sep = "\t",  row.names = FALSE)
  data.table::fwrite(clusterstats$cellsxdonor, paste(opt$outputpath, opt$cell, cluster.use, "cellsxdonor.txt", sep = "/"), quote = FALSE, sep = "\t",  row.names = FALSE)
  writeLines(clusterdonors, paste(opt$outputpath, opt$cell, cluster.use, "donors.txt", sep = "/"))
})


cellstats <- getsummarystats(inpObj, identifier = "full.donor.id.tag")
cellsdonors <- as.vector(gsub("g", "", colnames(cellstats$mean)))
meandf <- as.data.frame(cellstats$mean)
colnames(meandf) <- cellsdonors

dirpath <- paste(opt$outputpath, opt$cell, "AGGR", sep = "/")
if(!dir.exists(dirpath)){
  dir.create(dirpath, recursive = TRUE)
}
data.table::fwrite(meandf, paste(opt$outputpath, opt$cell, "AGGR", "gene_expression.txt", sep = "/"), quote = FALSE, sep = "\t",  row.names = TRUE)
data.table::fwrite(cellstats$freq, paste(opt$outputpath, opt$cell, "AGGR", "freq_cell.txt", sep = "/"), quote = FALSE, sep = "\t",  row.names = FALSE)
data.table::fwrite(cellstats$cellsxdonor, paste(opt$outputpath, opt$cell, "AGGR", "cellsxdonor.txt", sep = "/"), quote = FALSE, sep = "\t",  row.names = FALSE)
writeLines(cellsdonors, paste(opt$outputpath, opt$cell, "AGGR", "donors.txt", sep = "/"))
