library(TED)

df_Puram<-read.table("/home/cke/Puram/Puram_pseudobulk_scRNA.tsv",sep='\t',header=TRUE,row.names = 1)

df_pseudobulk_raw <- read.table("/home/cke/Puram/Puram_pseudobulk_fromraw.tsv",sep='\t',header=TRUE,row.names = 1)

# df_pseudobulk_harmony <- read.table("/home/cke/Puram/Puram_pseudobulk_fromharmony.tsv",sep='\t',header=TRUE,row.names = 1)

# read in marker genes
# marker_genes <- scan("/home/cke/Puram/top100DEGs.txt",character(),sep=',')
# marker_genes <- scan("/home/cke/Puram/top100markers_de_cor.txt",character(),sep=',')

label <- read.table("/home/cke/Puram/PuramHNSCC_cell_category_mappedCelltypes_simple.csv",sep=',',header=TRUE,row.names = 1)
label_subtype <- read.table("/home/cke/Puram/PuramHNSCC_cell_category_mappedCelltypes.csv",sep=',',header=TRUE,row.names = 1)

label <- subset(label, row.names(label) %in% row.names(df_Puram))
label_subtype <- subset(label_subtype, row.names(label_subtype) %in% row.names(df_Puram))

# fixing the problem summing characters in df
# class(df_Puram)<-'numeric'
# class(df_pseudobulk_raw)<-'numeric'

# df_Puram_t_filtered <- df_Puram_t[,marker_genes]

# class(df_Puram_t)
ref.data.filtered <- cleanup.genes(ref.dat=df_Puram, gene.type=c("RB","chrM","chrX","chrY"), species="hs", input.type="scRNA",exp.cells=5)

#run BayesPrism main func
res <- run.Ted(ref.dat=ref.data.filtered, X=df_pseudobulk_raw,cell.type.labels=label$cell_category,cell.subtype.labels=label_subtype$cell_category,
      tum.key='tumor',input.type="scRNA",n.cores=10,pdf.name="PuramPseudobulkRaw",seed=62422)

# ?run.Ted
saveRDS(res,"/home/cke/BayesPrism/BP_Puram_withsub_PseudobulkRaw_res.rds")
#final est cell type fraction in each bulk samples
saveRDS(res$res$final.gibbs.thet,"/home/cke/BayesPrism/BP_Puram_withsub_PseudobulkRaw_res_celltypeprop.rds")

