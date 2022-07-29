#!/usr/bin/env python
# coding: utf-8

# In[113]:


import sys, os
import argparse
import time
from Deconvolution.BLADE import Framework
import numpy as np
from numpy import transpose as t
import itertools
import pickle
from scipy.optimize import nnls
from scipy.stats import gaussian_kde
from matplotlib.colors import LogNorm
from sklearn.svm import SVR
from sklearn.svm import NuSVR

from sklearn.metrics import mean_squared_error as mse
import pandas as pd
from tqdm import trange,tqdm
# modules for visualization
import qgrid
from matplotlib import pyplot as plt
import seaborn as sns
import anndata as ad
import scanpy as sc
import scanorama
from sklearn.model_selection import LeaveOneOut,StratifiedKFold,KFold
import cycler

import warnings
warnings.filterwarnings('ignore')


# In[114]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, facecolor='white')


# In[115]:


def getloclist(wd,keyword=["top","marker","DEG"]):
    loc_list = []
    for root, dirs, files in os.walk(wd):
        for file in files:
            for key in keyword:
                if key in file:
                    loc = os.path.join(root, file)
                    loc_list.append(loc) #get file location
                    break
    return loc_list


# In[116]:


def record_time():
    t = time.localtime()
    current_time = time.strftime("%Y-%m-%d %H:%M:%S", t)
    print(current_time)


# Run BayesPrism first to get signature.csv file in place


def prepare_CIBERSORTx(path_label,path_signature):
    labels = pd.read_csv(path_label+"cellcategory_simple.csv",index_col=0)
    dict_label = labels.to_dict()['cell_category']
    scRNAref =  pd.read_csv(path_signature,index_col=0)
    scRNAref.columns = scRNAref.columns.map(dict_label)
    scRNAref.to_csv(path_signature.split(".")[0]+"_toCIBERSORTx.tsv",sep='\t')
    print("upload file:",path_signature.split(".")[0]+"_toCIBERSORTx.tsv","to CIBERSORTx as input")
    # see preapare_BLADE for test set, signature file from runBP 
    # run this func in notebook! after runBP done.


# ## File preparation
# ### Pseudobulk: Cross-validation of dataset


def prepare_BLADE(final_adata_mean,final_adata,mode,out):
    # Leave-one-out CV
#prepare var and mean signature matrix for BLADE
    merge_celltype = pd.merge(final_adata_mean.to_df(),final_adata.obs,left_index=True,right_index=True)
    if mode == 'real':
        counts_Puram_mean = merge_celltype.groupby(['Cell_type']).aggregate(np.mean).fillna(0)[final_adata_mean.to_df().columns]
        counts_Puram_std = merge_celltype.groupby(['Cell_type']).aggregate(np.std).fillna(0)[final_adata_mean.to_df().columns]
        counts_Puram_mean.T.to_csv(out+mode+"_mean.tsv",sep='\t')
        #save signature matrices for TCGA deconv
        counts_Puram_std.T.to_csv(out+mode+"_std.tsv",sep='\t')
        list_LOT = []

    elif mode == 'pseudobulk':
        merge_celltype_pseudobulk = merge_celltype.groupby(['batch']).sum()
        sample = merge_celltype_pseudobulk.index.tolist()
        list_LOT = []
        for train_index, test_index in LeaveOneOut().split(sample):
            print("LOT TRAIN:", train_index, "TEST:", test_index)
            train_sample = merge_celltype_pseudobulk.iloc[train_index,].index.tolist()
            train = merge_celltype[merge_celltype['batch'].isin(train_sample)]
            test_sample = merge_celltype_pseudobulk.iloc[test_index,].index.tolist()
            list_LOT.append("".join(str(i) for i in test_sample))
            print("Leaving out: ","".join(str(i) for i in test_sample))
            test = merge_celltype_pseudobulk[merge_celltype_pseudobulk.index.isin(test_sample)]
            counts_Puram_mean = train.groupby(['Cell_type']).aggregate(np.mean).fillna(0)[final_adata_mean.to_df().columns]
            counts_Puram_std = train.groupby(['Cell_type']).aggregate(np.std).fillna(0)[final_adata_mean.to_df().columns]
            counts_Puram_mean.T.to_csv(out+mode+"_LOT"+"".join(str(i) for i in test_sample)+"_mean.tsv",sep='\t')
            counts_Puram_std.T.to_csv(out+mode+"_LOT"+"".join(str(i) for i in test_sample)+"_std.tsv",sep='\t')
            # save the leftout testset for all methods
            test.iloc[:,:final_adata.to_df().shape[1]].to_csv(out+mode+"_LOT"+"".join(str(i) for i in test_sample)+"_test.tsv",sep='\t')
            # for CIBERSORTx
            test.iloc[:,:final_adata.to_df().shape[1]].T.to_csv(out+mode+"_LOT"+"".join(str(i) for i in test_sample)+"_test_transpose.tsv",sep='\t')
    return list_LOT

# In[121]:


def prepare_BLADE_kfoldCV(final_adata_mean,final_adata,mode,out,n_splits=5):
    #stratified n-fold CV
#prepare var and mean signature matrix for BLADE
    merge_celltype = pd.merge(final_adata_mean.to_df(),final_adata.obs,left_index=True,right_index=True)
    if mode == 'real':
        counts_Puram_mean = merge_celltype.groupby(['Cell_type']).aggregate(np.mean).fillna(0)[final_adata_mean.to_df().columns]
        counts_Puram_std = merge_celltype.groupby(['Cell_type']).aggregate(np.std).fillna(0)[final_adata_mean.to_df().columns]
        counts_Puram_mean.T.to_csv(out+mode+"_mean.tsv",sep='\t')
        #save signature matrices for TCGA deconv
        counts_Puram_std.T.to_csv(out+mode+"_std.tsv",sep='\t')
        list_CV = []

    elif mode == 'pseudobulk':
        singlecells = merge_celltype.index.tolist()
        batch = merge_celltype['batch'].tolist()
        list_CV = []
        count = 1
        for train_index, test_index in StratifiedKFold(n_splits=n_splits).split(singlecells,batch):
            print("CV TRAIN:",count , len(train_index), "TEST:", len(test_index))
            train_sc = merge_celltype.iloc[train_index,].index.tolist()
            train = merge_celltype.loc[train_sc,]
            test_sc = merge_celltype.iloc[test_index,].index.tolist()
            list_CV.append(str(count))
            test = merge_celltype.loc[test_sc,].groupby(['batch']).sum()
            counts_Puram_mean = train.groupby(['Cell_type']).aggregate(np.mean).fillna(0)[final_adata_mean.to_df().columns]
            counts_Puram_std = train.groupby(['Cell_type']).aggregate(np.std).fillna(0)[final_adata_mean.to_df().columns]
            counts_Puram_mean.T.to_csv(out+mode+"_CV"+str(count)+"_mean.tsv",sep='\t')
            counts_Puram_std.T.to_csv(out+mode+"_CV"+str(count)+"_std.tsv",sep='\t')
            # save the leftout testset for all methods
            test.iloc[:,:final_adata.to_df().shape[1]].to_csv(out+mode+"_CV"+str(count)+"_test.tsv",sep='\t')
            # for CIBERSORTx
            test.iloc[:,:final_adata.to_df().shape[1]].T.to_csv(out+mode+"_CV"+str(count)+"_test_transpose.tsv",sep='\t')
            count+=1
    return list_CV


def prepare_Rdeconv(final_adata,mode,out):
    merge_sample = pd.merge(final_adata.to_df(),final_adata.obs,left_index=True,right_index=True)
    if mode == 'real':
        scRNA_input = final_adata.to_df().loc[:,final_adata.to_df().columns]
        scRNA_input.to_csv(out+mode+"_scRNAref.tsv",sep='\t')
        list_LOT = []
        
    elif mode == 'pseudobulk':
        merge_sample_pseudobulk = merge_sample.groupby(['batch']).sum()
        sample = merge_sample_pseudobulk.index.tolist()
        list_LOT = []
        for train_index, test_index in LeaveOneOut().split(sample):
            print("LOT TRAIN:", train_index, "TEST:", test_index)
            train_sample = merge_sample_pseudobulk.iloc[train_index,].index.tolist()
            train = merge_sample[merge_sample['batch'].isin(train_sample)]
            test_sample = merge_sample_pseudobulk.iloc[test_index,].index.tolist()
            list_LOT.append("".join(str(i) for i in test_sample))
            print("Leaving out: ","".join(str(i) for i in test_sample))
            test = merge_sample_pseudobulk[merge_sample_pseudobulk.index.isin(test_sample)]
            train.iloc[:,:final_adata.to_df().shape[1]].to_csv(out+mode+"_LOT"+"".join(str(i) for i in test_sample)+"_scRNAtrain.tsv",sep='\t')
            test.iloc[:,:final_adata.to_df().shape[1]].to_csv(out+mode+"_LOT"+"".join(str(i) for i in test_sample)+"_test.tsv",sep='\t')
            # for CIBERSORTx
            test.iloc[:,:final_adata.to_df().shape[1]].T.to_csv(out+mode+"_LOT"+"".join(str(i) for i in test_sample)+"_test_transpose.tsv",sep='\t')
    return list_LOT



def prepare_Rdeconv_kfoldCV(final_adata,mode,out,n_splits=5):
    merge_sample = pd.merge(final_adata.to_df(),final_adata.obs,left_index=True,right_index=True)
    if mode == 'real':
        scRNA_input = final_adata.to_df().loc[:,final_adata.to_df().columns]
        scRNA_input.to_csv(out+mode+"_scRNAref.tsv",sep='\t')
        list_CV = []
        
    elif mode == 'pseudobulk':
        singlecells = merge_celltype.index.tolist()
        batch = merge_celltype['batch'].tolist()
        list_CV = []
        count = 1
        for train_index, test_index in StratifiedKFold(n_splits=n_splits).split(singlecells,batch):
            print("CV TRAIN:",count , len(train_index), "TEST:", len(test_index))
            train_sc = merge_celltype.iloc[train_index,].index.tolist()
            train = merge_celltype.loc[train_sc,]
            test_sc = merge_celltype.iloc[test_index,].index.tolist()
            list_CV.append(str(count))
            test = merge_celltype.loc[test_sc,].groupby(['batch']).sum()
            train.iloc[:,:final_adata.to_df().shape[1]].to_csv(out+mode+"_CV"+str(count)+"_scRNAtrain.tsv",sep='\t')
            test.iloc[:,:final_adata.to_df().shape[1]].to_csv(out+mode+"_CV"+str(count)+"_test.tsv",sep='\t')
            # for CIBERSORTx
            test.iloc[:,:final_adata.to_df().shape[1]].T.to_csv(out+mode+"_CV"+str(count)+"_test_transpose.tsv",sep='\t')
            count+=1
    return list_CV



def run_cmd(mode,out,out_res,dict_FS,folder_marker,methods,
            test_sample=False,path_label="home/cke/Puram/scRNAlabels/",
            path_bulk=False,name="unnamed_job"):
    # "out" is folder where input files prepared in last steps
    if mode == 'real':
        for FS_setup, marker_file in tqdm(dict_FS.items()):
            time.sleep(30)
            print(mode,"mode: now with feature selection setup: \r",FS_setup)
            record_time()
            if 'MuSiC' in methods:
                cmd_MuSiC = "".join(["Rscript /home/cke/runscripts/runMuSiC.r ", 
                                 out,mode,"_scRNAref.tsv ", # scRNA signature matrix
                                path_bulk," ", # real bulk RNAseq matrix
                                 path_label," ", 
                              marker_file," ",
                               name+"_"+mode+"_"+FS_setup+" ", # job name FS_setup
                              out_res+"MuSiC/"," &", # output frac folder
                                ])
                print("now running in cmd: \r",cmd_MuSiC)
                os.system(cmd_MuSiC)
            if 'BP' in methods:
                cmd_BP = "".join(["Rscript /home/cke/runscripts/runBP.r ", 
                                 out,mode,"_scRNAref.tsv ", # signature matrix
                                path_bulk," ", # testset pseudobulk matrix
                                 path_label," ",
                              marker_file," ",
                               name+"_"+mode+"_"+FS_setup+" ", # job name FS_setup
                              out_res+"BayesPrism/"," ", # output frac folder
                              out_res," &" # output CIBERSORTx prelim signature
                                ])
                print("now running in cmd: \r",cmd_BP)
                os.system(cmd_BP)
        # for BLADE, feature selection is done within wrapper
        # R is lame, plz use python to develop new tool :)
        if 'BLADE' in methods:
            cmd_BLADE = "".join(["python /home/cke/runscripts/runBLADE.py ", 
                             out,mode,"_std.tsv ", # std signature matrix
                             out,mode,"_mean.tsv ", # mean signature matrix
                            path_bulk," ", # testset pseudobulk matrix
                             out_res+"BLADE/"," ", # output folder
                             "--folder_marker ",folder_marker," ", # marker folder
                             "--name ",name+"_"+mode," &" # job name, background run
                            ])
            print("now running in cmd: \r",cmd_BLADE)
            os.system(cmd_BLADE)
    elif mode == "pseudobulk":
        for FS_setup, marker_file in tqdm(dict_FS.items()):
            time.sleep(30)
            record_time()
            print(mode,"mode: now with feature selection setup: \r",FS_setup)
            if 'MuSiC' in methods:
                cmd_MuSiC = "".join(["Rscript /home/cke/runscripts/runMuSiC.r ", 
                                 out,mode,"_LOT",str(test_sample),"_scRNAtrain.tsv ", # signature matrix
                                out,mode,"_LOT",str(test_sample),"_test.tsv ", # testset pseudobulk matrix
                                 path_label," ",
                              marker_file," ",
                               name+"_"+mode+"_"+FS_setup+"_"+str(test_sample)+"LOT ", # job name, FS_setup
                              out_res+"MuSiC/"," &", # output frac folder
                                ])
                print("now running in cmd: \r",cmd_MuSiC)
                os.system(cmd_MuSiC)
            if 'BP' in methods:   
                cmd_BP = "".join(["Rscript /home/cke/runscripts/runBP.r ", 
                                 out,mode,"_LOT",str(test_sample),"_scRNAtrain.tsv ", # signature matrix
                                out,mode,"_LOT",str(test_sample),"_test.tsv ", # testset pseudobulk matrix
                                 path_label," ",
                              marker_file," ",
                               name+"_"+mode+"_"+FS_setup+"_"+str(test_sample)+"LOT ", # job name, FS_setup
                              out_res+"BayesPrism/"," ", # output frac folder
                              out_res," &" # output CIBERSORTx prelim signature
                                ])
                print("now running in cmd: \r",cmd_BP)
                os.system(cmd_BP)
        # for BLADE, feature selection is done within wrapper
        # R is lame, plz use python to develop new tool :)
        if 'BLADE' in methods:
            cmd_BLADE = "".join(["python /home/cke/runscripts/runBLADE.py ", 
                             out,mode,"_LOT",str(test_sample),"_std.tsv ", # std signature matrix
                             out,mode,"_LOT"+str(test_sample),"_mean.tsv ", # mean signature matrix
                            out,mode,"_LOT",str(test_sample),"_test.tsv ", # testset pseudobulk matrix
                             out_res+"BLADE/"," ", # output folder
                             "--folder_marker ",folder_marker," ", # marker folder
                             "--name ",name+"_"+mode+"_"+"_LOT"+str(test_sample)," &" # job name, background run
                            ])
            print("now running in cmd: \r",cmd_BLADE)
            os.system(cmd_BLADE)

# In[132]:


def run_cmd_noFS(mode,out,out_res,methods,test_sample=False,
                 path_label="home/cke/Puram/scRNAlabels/",
                 path_bulk=False,name="unnamed_noFS_job"):
    # "out" is folder where input files prepared in last steps
    marker_file = "noFS"
    if mode == 'real':
        print(mode,"mode: no feature selection\r")
        record_time()
        if 'MuSiC' in methods:
            cmd_MuSiC = "".join(["Rscript /home/cke/runscripts/runMuSiC.r ", 
                             out,mode,"_scRNAref.tsv ", # scRNA signature matrix
                            path_bulk," ", # real bulk RNAseq matrix
                             path_label," ", 
                          marker_file," ",
                           name+"_"+mode+"_"+"noFS"+" ", # job name, FS_setup
                          out_res+"MuSiC/"," &", # output frac folder
                            ])
            print("now running in cmd: \r",cmd_MuSiC)
            os.system(cmd_MuSiC)
        if 'BP' in methods: 
            cmd_BP = "".join(["Rscript /home/cke/runscripts/runBP.r ", 
                             out,mode,"_scRNAref.tsv ", # signature matrix
                            path_bulk," ", # testset pseudobulk matrix
                             path_label," ",
                          marker_file," ",
                           name+"_"+mode+"_"+"noFS"+" ", # job name FS_setup
                          out_res+"BayesPrism/"," ", # output frac folder
                          out_res," &" # output CIBERSORTx prelim signature
                            ])
            print("now running in cmd: \r",cmd_BP)
            os.system(cmd_BP)
        # for BLADE, feature selection is done within wrapper
        # R is lame, plz use python to develop new tool :)
        if 'BLADE' in methods:
            cmd_BLADE = "".join(["python /home/cke/runscripts/runBLADE.py ", 
                             out,mode,"_std.tsv ", # std signature matrix
                             out,mode,"_mean.tsv ", # mean signature matrix
                            path_bulk," ", # testset pseudobulk matrix
                             out_res+"BLADE/"," ", # output folder
                             "--name ",name+"_"+mode,"_noFS &" # job name, background run
                            ])
            print("now running in cmd: \r",cmd_BLADE)
            os.system(cmd_BLADE)
    elif mode == "pseudobulk":
        print(mode,"mode: no feature selection\r")
        record_time()
        if 'MuSiC' in methods:
            cmd_MuSiC = "".join(["Rscript /home/cke/runscripts/runMuSiC.r ", 
                             out,mode,"_LOT",str(test_sample),"_scRNAtrain.tsv ", # signature matrix
                            out,mode,"_LOT",str(test_sample),"_test.tsv ", # testset pseudobulk matrix
                             path_label," ",
                          marker_file," ",
                           name+"_"+mode+"_"+"noFS"+"_"+str(test_sample)+"_LOT ", # job name, FS_setup
                          out_res+"MuSiC/"," &", # output frac folder
                            ])
            print("now running in cmd: \r",cmd_MuSiC)
            os.system(cmd_MuSiC)
        if 'BP' in methods: 
            cmd_BP = "".join(["Rscript /home/cke/runscripts/runBP.r ", 
                             out,mode,"_LOT",str(test_sample),"_scRNAtrain.tsv ", # signature matrix
                            out,mode,"_LOT",str(test_sample),"_test.tsv ", # testset pseudobulk matrix
                             path_label," ",
                          marker_file," ",
                           name+"_"+mode+"_"+"noFS"+"_"+str(test_sample)+"_LOT ", # job name, FS_setup
                          out_res+"BayesPrism/"," ", # output frac folder
                          out_res," &" # output CIBERSORTx prelim signature
                            ])
            print("now running in cmd: \r",cmd_BP)
            os.system(cmd_BP)
        # for BLADE, feature selection is done within wrapper
        # R is lame, plz use python to develop new tool :)
        if 'BLADE' in methods:
            cmd_BLADE = "".join(["python /home/cke/runscripts/runBLADE.py ", 
                             out,mode,"_LOT",str(test_sample),"_std.tsv ", # std signature matrix
                             out,mode,"_LOT"+str(test_sample),"_mean.tsv ", # mean signature matrix
                            out,mode,"_LOT",str(test_sample),"_test.tsv ", # testset pseudobulk matrix
                             out_res+"BLADE/"," ", # output folder
                             "--name ",name+"_"+mode+"_"+"LOT"+str(test_sample)+"_noFS"," &" # job name, background run
                            ])
            print("now running in cmd: \r",cmd_BLADE)
            os.system(cmd_BLADE)

# In[133]:


def main(path_adata,mode,out,out_res,path_label,path_bulk=False,folder_marker=False,
         name="unnamed_job",methods=['MuSiC','BP','BLADE'],keyword=["top","marker","DEG"]):
    # "out" is folder where input files prepared in last steps
#     print(folder_marker)
    record_time()
    final_adata = sc.read_h5ad(path_adata)
    final_adata.uns['log1p']["base"] = None
    final_adata_mean = final_adata.copy()
    sc.pp.log1p(final_adata_mean)
    if folder_marker:
        print("running Feature selection!\r")
        list_markers = getloclist(folder_marker,keyword)
        dict_FS = {}
        for marker_file in list_markers:
#             marker_genes =  pd.read_csv(marker_file,header=None).iloc[0,:]
            dict_FS[os.path.split(marker_file)[1].split("_")[0]] = marker_file
            
        if mode == 'real':
            list_CV = prepare_BLADE(final_adata_mean,final_adata,mode,out)
            list_CV = prepare_Rdeconv(final_adata,mode,out)
            run_cmd(mode,out,out_res,dict_FS,folder_marker,methods,
                    path_label=path_label,path_bulk=path_bulk,name=name)
            
        elif mode == 'pseudobulk':
            list_CV = prepare_BLADE(final_adata_mean,final_adata,mode,out)
            list_CV = prepare_Rdeconv(final_adata,mode,out)
            for test_sample in tqdm(list_CV):
                time.sleep(30)
                run_cmd(mode,out,out_res,dict_FS,folder_marker,methods,
                        test_sample=test_sample,path_label=path_label,name=name)
#                 if True: 
#                     prepare_CIBERSORTx(path_label,path_signature)
    else:
        print("no Feature selection!\r")
        if mode == 'real':
            list_CV = prepare_BLADE(final_adata_mean,final_adata,mode,out)
            list_CV = prepare_Rdeconv(final_adata,mode,out)
            run_cmd_noFS(mode,out,out_res,methods,
                    path_label=path_label,path_bulk=path_bulk,name=name)
        elif mode == 'pseudobulk':
            list_CV = prepare_BLADE(final_adata_mean,final_adata,mode,out)
            list_CV = prepare_Rdeconv(final_adata,mode,out)
            for test_sample in tqdm(list_CV):
                time.sleep(60)
                run_cmd_noFS(mode,out,out_res,methods,
                        test_sample=test_sample,path_label=path_label,name=name)
# In[129]:


def parse_args():
    """
        Parses inputs from the commandline.
        :return: inputs as a Namespace object
    """
    parser = argparse.ArgumentParser(description='Generates pipeline')
    # Arguments
    parser.add_argument('path_adata', help='directory of preprocessed raw scRNA anndata object')
    parser.add_argument('mode', help='scheme for data processing', choices=['pseudobulk','real'])
    parser.add_argument('out', help='output CV input directory')
    parser.add_argument('out_res', help='output of decon methods directory')
    parser.add_argument('path_label', help='labels of single-cell type identity directory')
    parser.add_argument('--path_bulk', help='bulk rnaseq data directory',default=False)
    parser.add_argument('--folder_marker', help='the folder where markers is stored',default=False)
    parser.add_argument('--name', help='give this job a name to help remember',default='unnamed_job')
    parser.add_argument('--methods', help='which methods you want to use',default=['MuSiC','BP','BLADE'])
    parser.add_argument('--keyword', help='keyword in marker file name to identify them',default=["top","marker","DEG"])
    return parser.parse_args()


# In[ ]:


if __name__ == "__main__":
    args = parse_args()
    path_adata = args.path_adata
    mode = args.mode
    out = args.out
    out_res = args.out_res
    path_label = args.path_label
    path_bulk = args.path_bulk
    name = args.name
    folder_marker = args.folder_marker
    methods = args.methods
    keyword = args.keyword
    if name == "unnamed_job":
        print("You did not name this job!\r")
    else:
        print(name," - pipeline initiated! Welcome, contact author for support: kechanglin1998@hotmail.com\r")

    if mode == "real":
        if path_bulk == False:
            raise ValueError("no bulk RNAseq data input! LOAD UP YOUR AMMO!\r")
        main(path_adata,mode,out,out_res,path_label,path_bulk=path_bulk,folder_marker=folder_marker,
         name=name,methods=methods,keyword=keyword)
    elif mode == 'pseudobulk':
        main(path_adata,mode,out,out_res,path_label,path_bulk=path_bulk,folder_marker=folder_marker,
         name=name,methods=methods,keyword=keyword)
    print("job finished! details of arguments used in this run: ",args)



# In[ ]:
