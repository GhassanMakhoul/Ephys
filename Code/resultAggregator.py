#I/O  Setup
import os
import sys
import getopt
import re
import h5py
import glob


#data and math packages
import numpy as np
import pandas as pd
from scipy.io import loadmat

def agg_sesh_df(h5file, str):
    sesh_df =[]
    with h5py.File(h5file, 'r') as f:
        for key in f.keys():
            df = entry_to_df(key, f[key])
            sesh_df.append(df)
    sesh_df = pd.concat(sesh_df)
    return sesh_df

def entry_to_df(key, resp_h5):
    """assembles df with 
    1. resp_region
    2. alphas
    3. TR
    """
    alphas = resp_h5['alphas']
    TR = resp_h5.attrs['Tr']
    df = pd.DataFrame(data=alphas, column='alphas') #TODO check shape
    df['TR'] = TR
    df['resp_reg'] = key.split("_")[-1] #messy but keyshould be 'response_RH14' for example
    return df



def get_stim_folders(subj: str, res_folder: str):
    """Returns fodlers of each stim trial

    Args:
        subj (str): subject ID, e.g. 'Spat30',
        res_folder (str): path to result from crp pipeline
    """
    inp_pattern = os.path.join(res_folder,subj,'*-*_*mA' )
    folders = glob.glob(inp_pattern)

    #TODO sanity checks
    return folders


def get_sesh_params(folder:str):
    """Returns stim contacts and ma

    Args:
        folder (str): folder should be of form "path/to/Con1_ConP2_#ma"
    """
    stim_ma = folder.split("/")[0]
    return stim_ma.split("_")

def agg_responses(subj: str, h5file: str, stim_folders: list, pathout: str):
    """aggregates all responses across all stim trials in to one
    mega dataframe

    Args
        subj (str): subject id
        h5file (str): name of hdf5 file should, eg 'stim_resp.hdf5'
        stim_folders (list): list of stim response dirs, full path
        pathout (str): where to save large df
    """
    dfs = []
    for folder in stim_folders:
        h5f = os.path.join(folder, h5file)
        df = agg_sesh_df(h5f)
        df['stim_reg'], df['ma'] = get_sesh_params(folder)
        dfs.append(df)
    agg_df = pd.concat(dfs)
    agg_df['subj'] = subj
    agg_df.to_csv(os.path.join(pathout, f"{subj}_stim.csv"))

def main(argv):
    subj = ''
    res_folder = '/mnt/ernie_main/Ghassan/ephys/data'
  
    opts, _ = getopt.getopt(argv,"s:i:p:",["subj=",'inpf=','pathout'])
    for opt, arg in opts:
        if opt in ("-i", 'inpf'):
            ma = arg
        elif opt in ("-s", "--subj"):
            subj = arg
        elif opt in ('-p', '--pathout'):
            pathout = arg
    stim_folders = get_stim_folders(subj, res_folder)
    agg_responses(subj, 'stim_resp.hdf5', stim_folders, pathout)
