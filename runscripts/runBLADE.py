#!/usr/bin/env python
# coding: utf-8

# In[45]:


import sys, os
import argparse

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
from tqdm import trange,tqdm
# modules for visualization
import qgrid
from matplotlib import pyplot as plt
import seaborn as sns

from sklearn.model_selection import LeaveOneOut


# ## BLADE runsrcipt 
# input:
# output:
# options:

# In[ ]:





# In[46]:


def run_BLADE(marker_genes, df_Puram_std, df_Puram_mean, df_bulk,celltype=10):
    if marker_genes is not None:
        marker_genes = marker_genes.drop_duplicates()
        df_Puram_std_filtered = df_Puram_std.loc[marker_genes,:]
        df_Puram_mean_filtered = df_Puram_mean.loc[marker_genes,:]

        merge_genes_mean = pd.merge(df_Puram_mean_filtered,df_bulk,left_index=True,right_index=True,how='inner')
        merge_genes_std = pd.merge(df_Puram_std_filtered,df_bulk,left_index=True,right_index=True,how='inner')
    else:
        merge_genes_mean = pd.merge(df_Puram_mean,df_bulk,left_index=True,right_index=True,how='inner')
        merge_genes_std = pd.merge(df_Puram_std,df_bulk,left_index=True,right_index=True,how='inner')

    print("Get mean and std exp!")

    #simple tumor cell type setup, there are 10 annotated cell types
    df_bulk_shared = merge_genes_mean.iloc[:,celltype:]
    df_shared_mean = merge_genes_mean.iloc[:,:celltype]
    df_shared_std = merge_genes_std.iloc[:,:celltype]

    print("Get common genes! ",df_shared_mean.shape[0])
    print("cell types: ",df_shared_mean.shape[1])
    print("bulk samples: ",df_bulk_shared.shape[1])
    return df_bulk_shared, df_shared_mean, df_shared_std


# In[40]:


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


# In[ ]:


def get_result(final_obj,df_bulk_shared,df_shared_mean,path_out,name,FS_setup="noFS"):
    obj = final_obj
    outcomes = {
        'BLADE': {
            'Fraction': t(obj.ExpF(obj.Beta)), 
            'Signature': np.mean(obj.Nu, 0), #group mode purification
            'HighRes': obj.Nu                #highresolution mode purification
        }}
    filtered_celltypefrac_BLADE = pd.DataFrame(outcomes['BLADE']['Fraction'])
    filtered_celltypefrac_BLADE.columns = df_bulk_shared.columns
    filtered_celltypefrac_BLADE.index = df_shared_mean.columns
    outfile = path_out+name+"_celltypefrac_BLADEout_"+FS_setup+'.csv'
    filtered_celltypefrac_BLADE.T.to_csv(outfile)
    return filtered_celltypefrac_BLADE


# In[43]:


def main(path_std,path_mean,path_bulk,path_out,folder_marker=False,name="unnamed_job",keyword=["top","marker","DEG"]):
    df_Puram_std = pd.read_csv(path_std,sep='\t',index_col=0)
    df_Puram_mean = pd.read_csv(path_mean,sep='\t',index_col=0)
    df_bulk = pd.read_csv(path_bulk,sep='\t',index_col=0).T
    
    hyperpars = {
        'Alpha': [1, 10],
        'Alpha0': [0.1, 1, 5],
        'Kappa0': [1, 0.5, 0.1],
        'SY': [1,0.3,0.5],
    }

    Nrep=3
    Nrepfinal=10
    Njob=20
    
    if folder_marker:
        list_markers = getloclist(folder_marker,keyword)
        dict_FS = {}
        for marker_file in list_markers:
            marker_genes =  pd.read_csv(marker_file,header=None).iloc[0,:]
            dict_FS[os.path.split(marker_file)[1].split("_")[0]] = marker_genes

        for FS_setup, marker_genes in tqdm(dict_FS.items()):
            print("now with feature selection setup: ",FS_setup)
            df_bulk_shared, df_shared_mean, df_shared_std = run_BLADE(marker_genes, df_Puram_std, df_Puram_mean, df_bulk)
            print("start BLADE!")
            Y = df_bulk_shared.to_numpy()
            mean = df_shared_mean.to_numpy() 
            sd = df_shared_std.to_numpy() 
            outfile = path_out+name+"_BLADEout_"+FS_setup+'.pickle'
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
            print("export to: ",outfile)
            filtered_celltypefrac_BLADE = get_result(final_obj,df_bulk_shared,df_shared_mean,path_out,name,FS_setup=FS_setup)    
            
    else:
        print("no feature selection on BLADE, be advised: You might be waiting for the end of the world!")
        df_bulk_shared, df_shared_mean, df_shared_std = run_BLADE(None, df_Puram_std, df_Puram_mean, df_bulk)
        print("start BLADE!")
        Y = df_bulk_shared.to_numpy()
        mean = df_shared_mean.to_numpy() 
        sd = df_shared_std.to_numpy() 
        outfile = path_out+name+"_BLADEout_"+"noFS"+'.pickle'
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
        print("export to: ",outfile)
        filtered_celltypefrac_BLADE = get_result(final_obj,df_bulk_shared,df_shared_mean,path_out,name)
        


# In[44]:


def parse_args():
    """
        Parses inputs from the commandline.
        :return: inputs as a Namespace object
    """
    parser = argparse.ArgumentParser(description='Generates pipeline')
    # Arguments
    parser.add_argument('path_std', help='variability signature directory')
    parser.add_argument('path_mean', help='mean signature directory')
    parser.add_argument('path_bulk', help='bulk rnaseq data directory')
    parser.add_argument('path_out', help='output cell type fractions directory')
    parser.add_argument('--folder_marker', help='the folder where markers is stored')
    parser.add_argument('--name', help='give this job a name to help remember')
    return parser.parse_args()


# In[ ]:


if __name__ == '__main__':
    args = parse_args()
    path_std = args.path_std
    path_mean = args.path_mean
    path_bulk = args.path_bulk
    path_out = args.path_out
    folder_marker = args.folder_marker
    name = args.name
    main(path_std,path_mean,path_bulk,path_out,folder_marker,name,keyword=["top","marker","DEG"])

