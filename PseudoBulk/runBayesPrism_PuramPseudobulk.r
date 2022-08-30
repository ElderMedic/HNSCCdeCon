library(TED)

df_Puram <- read.table("/home/cke/Puram/Puram_pseudobulk_scRNAtrain.tsv",sep='\t',header=TRUE,row.names = 1,check.names = FALSE)

df_pseudobulk_raw <- read.table("/home/cke/Puram/Puram_pseudobulk_fromraw_test.tsv",sep='\t',header=TRUE,row.names = 1,check.names = FALSE)
df_pseudobulk_raw = data.matrix(df_pseudobulk_raw)

# read in marker genes
marker_genes_100DEGs <- scan("/home/cke/Puram/top100DEGs_pseudobulk.txt",character(),sep=',')
marker_genes_100 <- scan("/home/cke/Puram/top100markers_de_cor_symbol.txt",character(),sep=',')
marker_genes_50 <- scan("/home/cke/Puram/top50markers_de_cor_symbol.txt",character(),sep=',')
marker_genes_20 <- scan("/home/cke/Puram/top20markers_de_cor_symbol.txt",character(),sep=',')

label <- read.table("/home/cke/Puram/PuramHNSCC_cell_category_mappedCelltypes_simple.csv",sep=',',header=TRUE,row.names = 1)
label_subtype <- read.table("/home/cke/Puram/PuramHNSCC_cell_category_mappedCelltypes.csv",sep=',',header=TRUE,row.names = 1)

label <- subset(label, row.names(label) %in% row.names(df_Puram))
label_subtype <- subset(label_subtype, row.names(label_subtype) %in% row.names(df_Puram))

Get_filtered_df <- function(marker_genes,df_Puram,num){
    df_Puram_filtered = df_Puram[,marker_genes]
    write.csv(t(df_Puram_filtered),paste0("/home/cke/Puram/Puram_pseudobulk_scRNAref_",num,"_signature.tsv"))
    return(df_Puram_filtered)
}

run_BP <- function(df_Puram_filtered,df_pseudobulk_raw,label){
    ref.data.filtered <- cleanup.genes(ref.dat=data.matrix(df_Puram_filtered), gene.type=c("RB","chrM","chrX","chrY"), species="hs", input.type="scRNA",exp.cells=5)
    #run BayesPrism main func
    res <- run.Ted(ref.dat=ref.data.filtered, X=df_pseudobulk_raw,cell.type.labels=label$cell_category,cell.subtype.labels=label_subtype$cell_category,
        tum.key='tumor',input.type="scRNA",n.cores=10,pdf.name="PuramPseudobulkRaw",seed=71822)
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

df_Puram_filtered_100DEGs <- Get_filtered_df(marker_genes_100DEGs,df_Puram,"100DEGs")
res_100DEGs <- run_BP(df_Puram_filtered_100DEGs,df_pseudobulk_raw,label)
BP_res_100DEGs <- get_BP_res(res_100DEGs)
saveRDS(res_100DEGs,"/home/cke/BayesPrism/BP_Puram_withsub_PseudobulkRaw_res_100DEGs.rds")
write.csv(BP_res_100DEGs, "/home/cke/BayesPrism/celltypefrac_BP_Puramfiltered_pseudobulk_100DEGs.csv")

df_Puram_filtered_100 <- Get_filtered_df(marker_genes_100,df_Puram,"100")
res_100 <- run_BP(df_Puram_filtered_100,df_pseudobulk_raw,label)
BP_res_100 <- get_BP_res(res_100)
saveRDS(res_100,"/home/cke/BayesPrism/BP_Puram_withsub_PseudobulkRaw_res_100.rds")
write.csv(BP_res_100, "/home/cke/BayesPrism/celltypefrac_BP_Puramfiltered_pseudobulk_100.csv")

df_Puram_filtered_50 <- Get_filtered_df(marker_genes_50,df_Puram,"50")
res_50 <- run_BP(df_Puram_filtered_50,df_pseudobulk_raw,label)
BP_res_50 <- get_BP_res(res_50)
saveRDS(res_50,"/home/cke/BayesPrism/BP_Puram_withsub_PseudobulkRaw_res_50.rds")
write.csv(BP_res_50, "/home/cke/BayesPrism/celltypefrac_BP_Puramfiltered_pseudobulk_50.csv")

df_Puram_filtered_20 <- Get_filtered_df(marker_genes_20,df_Puram,"20")
res_20 <- run_BP(df_Puram_filtered_20,df_pseudobulk_raw,label)
BP_res_20 <- get_BP_res(res_20)
saveRDS(res_20,"/home/cke/BayesPrism/BP_Puram_withsub_PseudobulkRaw_res_20.rds")
write.csv(BP_res_20, "/home/cke/BayesPrism/celltypefrac_BP_Puramfiltered_pseudobulk_20.csv")