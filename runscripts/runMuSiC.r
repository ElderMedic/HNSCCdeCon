library(MuSiC,quietly = T)
library(Biobase,quietly = T)
library(ggplot2,quietly = T)
library(dplyr,quietly = T)
library(ggcorrplot,quietly = T)
library(reshape,quietly = T)
library(corrplot,quietly = T)
library(tidyr,quietly = T)
library(caret,quietly = T)

args = commandArgs(trailingOnly=TRUE)
print("running BayesPrism with following args:")

if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)==6) {
  # default subtype location
  args[7] = args[3]
}
print(args)

path_scRNA <- args[1]
path_bulk <- args[2]
path_label <- args[3]
path_marker <- args[4]
FS_setup <- args[5]
path_out <- args[6]
path_label_subtype <- args[7]

df_Puram <- read.table(path_scRNA,sep='\t',header=TRUE,row.names = 1,check.names = FALSE)
df_Puram = t(df_Puram)

df_bulk <- read.table(path_bulk,sep='\t',header=TRUE,row.names = 1,check.names = FALSE)
df_bulk <- t(df_bulk)

# for LOT setup, copy a test pseudobulk sample!
if (grepl('LOT',FS_setup)) {
    df_bulk_new <- cbind(df_bulk,df_bulk[,colnames(df_bulk)[1]])
    colnames(df_bulk_new) <- c(colnames(df_bulk),0)
    df_bulk <- df_bulk_new
}


label <- read.table(paste0(path_label,"cellcategory_simple.csv"),sep=',',header=TRUE,row.names = 1)
label_subtype <- read.table(paste0(path_label_subtype,"cellcategory_subtype.csv"),sep=',',header=TRUE,row.names = 1)

label$subtype <- label_subtype$cell_category
label$sampleID <- rownames(label)
label <- subset(label, row.names(label) %in% colnames(df_Puram))
label_subtype <- subset(label_subtype, row.names(label_subtype) %in% colnames(df_Puram))

if (path_marker!='noFS') {
    marker_genes <- scan(path_marker,character(),sep=',')
    }

# generate expressionset object and run MuSiC main function, process result cell type fraction
run_MuSiC <- function(df_Puram_filtered,df_bulk,label){
    # prepare scRNA-seq ref data, as expressionset object
    metadata <- data.frame(labelDescription= c("sampleID","cell_category", "subtype"), row.names=c("sampleID","cell_category", "subtype"))
    Puram.eset = ExpressionSet(assayData = data.matrix(df_Puram_filtered), phenoData =  new("AnnotatedDataFrame", data = label, varMetadata = metadata) )
    bulk_raw.eset = ExpressionSet(assayData = data.matrix(df_bulk))
    # Estimate cell type proportions
    Est.prop.tcga_raw = music_prop(bulk.eset = bulk_raw.eset, sc.eset = Puram.eset, clusters = 'cell_category',samples = 'sampleID')
    MuSiC_res <- data.matrix(Est.prop.tcga_raw$Est.prop.weighted)
    NNLS_res <- data.frame(Est.prop.tcga_raw$Est.prop.allgene)
    colnames(MuSiC_res)[5] <- 'other'
    MuSiC_res <- MuSiC_res[, order(colnames(MuSiC_res))]
    colnames(MuSiC_res)[9] <- 'T-cell'
    colnames(MuSiC_res)[1] <- 'B-cell'
    
    colnames(NNLS_res)[5] <- 'other'
    NNLS_res <- NNLS_res[, order(colnames(NNLS_res))]
    colnames(NNLS_res)[9] <- 'T-cell'
    colnames(NNLS_res)[1] <- 'B-cell'
    return(list(MuSiC_res,NNLS_res))
}

if (path_marker!='noFS'){
    df_Puram_filtered <- subset(df_Puram, row.names(df_Puram) %in% marker_genes)
} else {
    df_Puram_filtered <- df_Puram
}

list_res <- run_MuSiC(df_Puram_filtered,df_bulk,label)
MuSiC_res <- list_res[1]
NNLS_res <- list_res[2]
write.csv(MuSiC_res, paste0(path_out,"celltypefrac_MuSiC_",FS_setup,".csv"))
write.csv(NNLS_res, paste0(path_out,"celltypefrac_NNLS_",FS_setup,".csv"))


