#!/usr/bin/env python
# coding: utf-8

# In[4]:


import sys, os
from Deconvolution.BLADE import Framework
import numpy as np
from numpy import transpose as t
import itertools
import pickle
from scipy.optimize import nnls
from sklearn.svm import SVR
from sklearn.svm import NuSVR
from sklearn.metrics import mean_squared_error as mse
import pandas as pd

# modules for visualization
import qgrid
from matplotlib import pyplot as plt
import seaborn as sns
import logging

# logging  
#LOG = "/home/cke/runBLADE.log"                                                     
#logging.basicConfig(filename=LOG, filemode="w", level=logging.DEBUG,format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %H:%M:%S')  

# console handler  
#console = logging.StreamHandler()  
#console.setLevel(logging.ERROR)  
#logging.getLogger("").addHandler(console)

# ### Run BLADE with TCGA bulk and Puram scRNA-seq reference
# 

# #### Application of deconvolution methods
# 
# From here, we will apply the following three methods for further performance comparison:
# 1. BLADE (estimation of cellular fraction + group-mode/high-resolution-mode purification)
# 2. NNLS (estimation of fraction)
# 3. SVR followed by NNLS (estimation of fraction + group-mode purification) - similar to CIBERSORTx
# 
# 
# ##### 1. Application of BLADE
# 
# These are the key parameters used in BLADE (note that there is default settings of these parameters, if not specified):
# - Hyperparameters (`hyperpars`): `Alpha`, `Alpha0`, `Kappa0` and `SigmaY`, each of which can be defined as a list of options. BLADE takes an empirical Bayes approach to find the optimal parameter set given the all possible combinations. 
# - `Nrep`: Number of repeat for evaluating each parameter configuration.
# - `Nrepfinal`: Number of repeated optimizations for the final parameter set.
# - `Njob`: Number of parallel jobs.

# In[17]:


hyperpars = {
    'Alpha': [1, 10],
    'Alpha0': [0.1, 0.5, 1, 5, 10],
    'Kappa0': [1, 0.5, 0.1],
    'SY': np.sqrt([0.1, 0.5, 1, 1.5, 2])
}

Nrep=3
Nrepfinal=10
Njob=20


# In[10]:

# read in marker genes
marker_genes = pd.read_csv("/home/cke/Puram/top100markers_de_cor.txt",header=None).iloc[0,:]

# df_Puram_std = pd.read_csv("/home/cke/Puram/HNSCC2PuramGSE103322_HNSCC_exp_std.tsv",sep='\t',index_col=0)
# df_Puram_mean = pd.read_csv("/home/cke/Puram/HNSCC2PuramGSE103322_HNSCC_exp_mean.tsv",sep='\t',index_col=0)

# merged all tumor cell types
df_Puram_std = pd.read_csv("/home/cke/Puram/HNSCC2PuramGSE103322_HNSCC_exp_std_simple.tsv",sep='\t',index_col=0)
df_Puram_mean = pd.read_csv("/home/cke/Puram/HNSCC2PuramGSE103322_HNSCC_exp_mean_simple.tsv",sep='\t',index_col=0)

df_Puram_std_filtered = df_Puram_std.loc[marker_genes,:]
df_Puram_mean_filtered = df_Puram_mean.loc[marker_genes,:]

df_TCGA = pd.read_csv("/home/cke/TCGA-HNSC.htseq_counts_exp2.tsv",sep='\t',index_col=0)
df_Puram_mean_log2 = np.log2(df_Puram_mean_filtered+1)

merge_genes_mean = pd.merge(df_Puram_mean_log2,df_TCGA,left_index=True,right_index=True,how='inner')
merge_genes_std = pd.merge(df_Puram_std_filtered,df_TCGA,left_index=True,right_index=True,how='inner')

print("Get mean and std exp!")

# 21706 genes in common
# df_TCGA_shared = merge_genes_mean.iloc[:,24:]
# df_shared_mean = merge_genes_mean.iloc[:,:24]
# df_shared_std = merge_genes_std.iloc[:,:24]

#simple tumor cell type setup
df_TCGA_shared = merge_genes_mean.iloc[:,10:]
df_shared_mean = merge_genes_mean.iloc[:,:10]
df_shared_std = merge_genes_std.iloc[:,:10]

print("Get common genes! ",df_shared_mean.shape[0])
print("cell types: ",df_shared_mean.shape[1])
print("bulk samples: ",df_TCGA_shared.shape[1])


# Given the configuration above, BLADE is applied to each of the simulation dataset created previously.  
# 
# BLADE produce several outcomes:
# - `final_obj`: final BLADE object with optimized variational parameters
# - `best_obj`: BLADE object trained with the best parameter set found by the Empirical Bayes framework. Empirical Bayes framework is applied after selecting a subset of samples (5 samples; indicated by `Ind_sample` below), and thus the outcome contains only 5 samples. If `Nsample` <= 5, `final_obj` is identical to `best_obj`.
# - `best_set`: Best parameter set defined by Empirical Bayes framework.
# - `outs`: Outcome of BLADE for every possible combination of hyperparameters, used in the Empirical Bayes framework. 
# 

# - There are nan in mean and std matrix! NAs are filled with 0?
# 
# full tumor type setup:
# - ngenes = 21706 common genes
# - ncells = 24, including all 16 tumor types
# - nsample = 546
# 
# simple tumor type setup:
#- ngenes = 21706 common genes
#- ncells = 10, all tumor types are merged, including one NA type?
#- nsample = 546
# - marker genes = 900 (including 9 genes not shared)

print("start BLADE!")
Y = df_TCGA_shared.to_numpy()
mean = df_shared_mean.to_numpy() 
sd = df_shared_std.to_numpy() 

outfile = './BLADE/data/Puramfiltered_TCGA_corDEmarkers_BLADEout.pickle'

final_obj, best_obj, best_set, outs = Framework(
    mean, sd, Y,
    Alphas=hyperpars['Alpha'], Alpha0s=hyperpars['Alpha0'], 
    Kappa0s=hyperpars['Kappa0'], SYs=hyperpars['SY'],
    Nrep=Nrep, Njob=Njob, Nrepfinal=Nrepfinal)

pickle.dump(
    {
        'final_obj': final_obj,
        'best_obj': best_obj,
        'best_set': best_set,
        'outs' : outs
    }, open(outfile, 'wb')
    )

# In[ ]:




