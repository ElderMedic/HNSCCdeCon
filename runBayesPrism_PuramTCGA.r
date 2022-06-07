library(TED)

df_Puram<-read.table("/home/cke/Puram/HNSCC2PuramGSE103322_HNSCC_exp.tsv",sep='\t',header=TRUE)
df_TCGA <- read.table("/home/cke/TCGA-HNSC.htseq_counts_exp2.tsv",sep='\t',header=TRUE)

df_TCGA_t<-t(df_TCGA)
colnames(df_TCGA_t)<-df_TCGA_t[1,]
df_TCGA_t <- df_TCGA_t[-1,]

df_Puram_t <- t(df_Puram)
colnames(df_Puram_t)<-df_Puram_t[1,]
df_Puram_t <- df_Puram_t[-1,]

# read in marker genes
marker_genes <- scan("/home/cke/Puram/top100DEGs.txt",character(),sep=',')

label <- read.table("/home/cke/Puram/PuramHNSCC_cell_category_mappedCelltypes_simple.csv",sep=',',header=TRUE)

# head(celltypelab)

# fixing the problem summing characters in df
class(df_Puram_t)<-'numeric'
class(df_TCGA_t)<-'numeric'

df_Puram_t_filtered <- df_Puram_t[,marker_genes]

# class(df_Puram_t)
ref.data.filtered <- cleanup.genes(ref.dat=df_Puram_t_filtered, gene.type=c("RB","chrM","chrX","chrY"), species="hs", input.type="scRNA",exp.cells=5)

#run BayesPrism main func
res <- run.Ted(ref.dat=ref.data.filtered, X=df_TCGA_t,cell.type.labels=label$cell_category,
               tum.key='tumor',input.type="scRNA",n.cores=20,pdf.name="PuramTCGA",seed=60722)

# ?run.Ted
saveRDS(res,"/home/cke/BayesPrism/BP_Puramfiltered_TCGA_res.rds")
#final est cell type fraction in each bulk samples
saveRDS(res$res$final.gibbs.thet,"/home/cke/BayesPrism/BP_Puramfiltered_TCGA_res_celltypeprop.rds")

