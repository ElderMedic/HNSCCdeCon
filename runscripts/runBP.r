#!/usr/bin/Rscript
library(TED)

args = commandArgs(trailingOnly=TRUE)
print("running BayesPrism with following args:")


if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)==7) {
  # default subtype location
  args[8] = args[3]
}
print(args)

path_scRNA <- args[1]
path_bulk <- args[2]
path_label <- args[3]
path_marker <- args[4]
FS_setup <- args[5]
path_out <- args[6]
path_signature_out <- args[7]
path_label_subtype <- args[8]

df_Puram <- read.table(path_scRNA,sep='\t',header=TRUE,row.names = 1,check.names = FALSE)

df_bulk <- read.table(path_bulk,sep='\t',header=TRUE,row.names = 1,check.names = FALSE)
df_bulk = data.matrix(df_bulk)

label <- read.table(paste0(path_label,"cellcategory_simple.csv"),sep=',',header=TRUE,row.names = 1)
label_subtype <- read.table(paste0(path_label_subtype,"cellcategory_subtype.csv"),sep=',',header=TRUE,row.names = 1)

label <- subset(label, row.names(label) %in% row.names(df_Puram))
label_subtype <- subset(label_subtype, row.names(label_subtype) %in% row.names(df_Puram))

if (path_marker!='noFS') {
    marker_genes <- scan(path_marker,character(),sep=',')
    }

Get_filtered_df <- function(marker_genes,df_Puram,FS_setup,path_signature_out){
    df_Puram_filtered = df_Puram[,marker_genes]
    #writing signature for CIBERSORTx
    write.csv(t(df_Puram_filtered),paste0(path_signature_out,FS_setup,"_signature.tsv"))
    return(df_Puram_filtered)
}

run_BP <- function(df_Puram_filtered,df_bulk,label,label_subtype){
    ref.data.filtered <- cleanup.genes(ref.dat=data.matrix(df_Puram_filtered), gene.type=c("RB","chrM","chrX","chrY"), species="hs", input.type="scRNA",exp.cells=5)
    #run BayesPrism main func
    res <- run.Ted(ref.dat=ref.data.filtered, X=df_bulk,cell.type.labels=label$cell_category,cell.subtype.labels=label_subtype$cell_category,
        tum.key='tumor',input.type="scRNA",n.cores=10,pdf.name="PuramPseudobulkRaw")
    return(res)
}

get_BP_res <- function(res){
    BP_res <- data.frame(res$res$final.gibbs.thet)
    colnames(BP_res)[8] <- 'T-cell'
    colnames(BP_res)[3] <- 'B-cell'
    colnames(BP_res)[5] <- 'other'
    BP_res <- BP_res[, order(names(BP_res))]
    return(BP_res)
}

if (path_marker!='noFS'){
    df_Puram_filtered <- Get_filtered_df(marker_genes,df_Puram,FS_setup,path_signature_out)
} else {
    df_Puram_filtered <- df_Puram
}
res <- run_BP(df_Puram_filtered,df_bulk,label,label_subtype)
BP_res <- get_BP_res(res)

# ?run.Ted
saveRDS(res,paste0(path_out,FS_setup,"_BP_res.rds"))
write.csv(BP_res, paste0(path_out,"celltypefrac_BP_",FS_setup,".csv"))
