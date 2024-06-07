import logging
import os
import re
from scipy.io import loadmat
import mat73
import glob
import pdb
from multiprocessing import Pool


from collections import Counter
import pandas as pd
import numpy as np
import mne

import seaborn as sns
import matplotlib.pyplot as plt
import tqdm
from utils import *

#TODO: work through and log relevant output for mass run
#TODO: make more modular for ANY connectivity metric, not just PDC
#TODO: refine assemble_net_conn to iterate through regions
#TODO: make assemble_net_conn more modular so it can be extended to exclude certain within soz connections

TIMELINE_F = '/mnt/ernie_main/000_Data/SEEG/SEEG_Periictal/data/Extracted_Per_Event_Interictal/all_time_data_01042023_212306.csv'
SEEG_FOLDER = '/mnt/ernie_main/000_Data/SEEG/SEEG_Entire_EMU_Downloads/data/'
BANDS = ['delta', 'theta', 'alpha', 'beta','gamma_l', 'gamma_H']
PERIOD = ['inter','pre','ictal','post']

#configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename='../logs/conn_dyn.log', level=logging.WARNING)
def load_mat(f):
    try:
        return loadmat(f)
    except NotImplementedError:
        try:
            return mat73.loadmat(f)
        except:
            print(f"Problem Loading {f}")
            return None
    

def get_pat_conn_labels(pat_df, regions, subj_id):
    """Returns an array of channels as NZ, SOZ, PZ, etc

    Args:
        pat_df (_type_): dataframe with subj and labels
    """
    assert len(set(pat_df.subj)) == 1, "Multiple subjects not allowed!"
    assert pat_df.subj.values[0] == subj_id, "incorrect subject dataframe"
    contact_designation = pat_df.label.values
    bipoles = pat_df.bipole.values
    bip_dict = dict(zip(bipoles, contact_designation))
    pat_conn_labels = np.array([bip_dict[b] for b in regions])
    return pat_conn_labels

def get_reg_inds(pat_conn_labels):
    soz_inds = np.where(pat_conn_labels == 'SOZ')[0]
    pz_inds = np.where(pat_conn_labels == 'PZ')[0]
    nz_bool = pat_conn_labels == 'NIZ'
    iz_bool = pat_conn_labels =='IZ'
    nz_inds = np.where(np.logical_or(nz_bool, iz_bool))[0]
    return {'soz':soz_inds, 'pz' : pz_inds, 'nz': nz_inds}

def get_conn_dict(conn_obj,key='pdc', periods=PERIOD, filt_dist=0,**kwargs):
    #TODO fix final keys
    inter_conn = read_conn_struct(conn_obj, key,'PDC_interictal')
    pre_conn = read_conn_struct(conn_obj, key, 'PDC_pre')
    ictal_conn = read_conn_struct(conn_obj, key, 'PDC_ictal')
    post_conn = read_conn_struct(conn_obj, key, 'PDC_post')
    if filt_dist > 0:
        dist_mat = conn_obj['pdc']['seizure']['dist_mat']
        inter_conn = filter_dist(inter_conn, dist_mat, filt_dist)
        pre_conn = filter_dist(pre_conn, dist_mat, filt_dist)
        ictal_conn = filter_dist(ictal_conn, dist_mat, filt_dist)
        post_conn = filter_dist(post_conn, dist_mat, filt_dist)
    conn_dict = dict(zip(periods, [inter_conn, pre_conn, ictal_conn, post_conn]))
    return conn_dict

def read_conn_struct(conn_obj, metric_key, data_key):
    try:
        conn_data = conn_obj[metric_key]['seizure'][data_key]
    except KeyError:
        conn_data = conn_obj['seizure'][data_key]
    return conn_data

def filter_dist(conn_mat:np.ndarray, dist_mat:np.ndarray, filt_dist:float)-> np.ndarray:
    """removes connections from a functional connectivity matrix if the euclidean
    distance is under a FILT_DIST threshold

    Args:
        conn_mat (np.ndarray): functional connectivity matrix, can be shape NxN or K_BAND x N x N
        filt_dist (np.ndarray): N x N euclidean distance matrix_

    Returns:
        np.ndarray: filtered connectivity matrix    
    """
    d_x, d_y= np.where(dist_mat < filt_dist)
    if len(conn_mat.shape) == 3:
        assert conn_mat.shape[1] > conn_mat.shape[0], "Weird shape order, check connectivity matrix!"
        conn_mat[:, d_x, d_y] = np.nan
        return conn_mat
    conn_mat[d_x, d_y] = np.nan
    return conn_mat

def assemble_net_conn(subj_id, pdc_dict,soz_inds,pz_inds,nz_inds,bands=BANDS):
    """Assemble net directed connectivity across periods of interest for all 
    frequency bands

    Args:
        subj_id (_type_): string of subject id
        pdc_dict (_type_): dictionary mapping period (str) to pdc (np.array)
        soz_inds (_type_): indices corresponding to soz nodes
        pz_inds (_type_): pz node indices
        nz_inds (_type_): nz node indices

    Returns:
        _type_: Dataframe of next connectivity for each region.
    """
    net_df = pd.DataFrame(columns=['subj', 'period', 'region', 'net_pdc', 'freq_band'])
    ind = 0 
    for period, pdc in pdc_dict.items():
        for b in range(len(bands)):
            z_pdc_in = z_score_conn(pdc[b, :,:],direction='col')
            z_pdc_out = z_score_conn(pdc[b,:,:], direction='row')

            band = bands[b]
            soz_in = z_pdc_out[:,soz_inds]
            soz_out = z_pdc_in[soz_inds,:]
            net_soz = np.nanmean(soz_in) - np.nanmean(soz_out)
            net_df.loc[ind] = [subj_id,period,'soz',net_soz,band]
            ind += 1

            pz_in = z_pdc_out[:,pz_inds]
            pz_out = z_pdc_in[pz_inds,:]
            net_pz = np.nanmean(pz_in) - np.nanmean(pz_out)
            net_df.loc[ind] = [subj_id,period,'pz',net_pz,band]
            ind +=1


            nz_in = z_pdc_out[:,nz_inds]
            nz_out = z_pdc_in[nz_inds,:]
            net_nz = np.nanmean(nz_in) - np.nanmean(nz_out)
            net_df.loc[ind] = [subj_id,period,'nz',net_nz,band]
            ind +=1
    return net_df    

def assemble_file(subj_id: str, conn_f: str, label_df: pd.DataFrame, **kwargs)->pd.DataFrame:
    """For a given subject's connectivity struct, return the dataframe containing net
    connectivity for each frequency band, over relevant periods
    NOTE: this will only work with DIRECTED connecitivity

    Args:
        subj_id (str): subject id ex Epat02
        conn_f (str): .mat struct with connectivity matrices
        label_df (df): df of channel labels -> {SOZ_inds,NZ_inds, PZ_inds}

    Returns:
        pd.DataFrame: net connecitivity
    """

    
    conn_obj = load_mat(conn_f)
    
    conn_dict, label_inds = prep_conn(subj_id, label_df, conn_obj, **kwargs)
    soz_inds, pz_inds, nz_inds = label_inds['soz'], label_inds['pz'], label_inds['nz']
    return assemble_net_conn(subj_id, conn_dict, soz_inds, pz_inds,nz_inds)

def prep_conn(subj_id, label_df, conn_obj, **kwargs):
    if conn_obj == None:
        raise ValueError("Issue with load_mat")
    conn_dict = get_conn_dict(conn_obj, **kwargs)
    pat_conn_labels = get_regions(label_df, conn_obj,subj_id)
    label_inds = get_reg_inds(pat_conn_labels)
    return conn_dict, label_inds

def get_regions(label_df, conn_obj, subj_id):
    regions = [reg[0] for reg in read_conn_struct(conn_obj, 'pdc', 'bip_labels_used')]
    
    regions = format_bipoles(regions)
    pat_conn_labels = get_pat_conn_labels(label_df, regions,subj_id)
    return pat_conn_labels

def agg_patient(subj_id, conn_folder, label_df, **kwargs) -> pd.DataFrame:
    """Aggregate all connecivity files within a folder and return as one df
    NOTE: assumes tht all .mat file in this folder are connectivity structs to analyze
    """
    files = glob.glob(os.path.join(conn_folder, "*.mat"))

    
    conn_dfs = []
    for f in files:
        try:
            net_df = assemble_file(subj_id, f, label_df, **kwargs)
            net_df['conn_file'] = os.path.basename(f)
            conn_dfs.append(net_df)
        except KeyError as e:
            logger.debug(f"Issue hanndling keys in bipole labels for {subj_id} in this file: {f}\n\t\tMessage: {e}")
        except ValueError as e:
            logger.debug(f"Issue with load_mat for {subj_id} in file: {f}\n\t\tMessage{e}")
            
    if len(conn_dfs) == 0:
        return pd.DataFrame()
    return pd.concat(conn_dfs)

def agg_subjects(folders: list[str], subject_ids: list[str], bipole_label_df:pd.DataFrame, **kwargs ) -> pd.DataFrame:
    """Given a list of folders and subjects aggregate connectivity matrices across subjects

    Args:
        folders (list[str]): list of subject folders
        subject_ids (list[str]): list of subject IDS NOTE: subj_ids must be aligned with subj folders!
        bipole_label_df (pd.DataFrame): df with bipole labels should have following columns: 'subj',
        'bipole', 'label'. Used to designate SOZ NZ, ets

    Returns:
        pd.DataFrame: _description_
    """
    #first check that there are .mat files in the folders
    folders = check_folders(folders)
    conn_dfs = []

    for i, conn_dir in enumerate(folders):
        subj_id = subject_ids[i]
        label_df  = bipole_label_df[bipole_label_df.subj == subj_id]
        df = agg_patient(subj_id, conn_dir, label_df, **kwargs)
        if df.shape[0] > 0:
            conn_dfs.append(df)
    return pd.concat(conn_dfs)
        
    
def check_folders(folders:list[str])-> list[str]:
    """Standard input checker to make sure that each folder has a .mat file in it

    Args:
        folders (list[str]): list of folders containing .mat structs of connecivity

    Returns:
        list[str]: filtered list with empty folders filtered out
    """
    
    return [f for f in folders if len(glob.glob(os.path.join(f, "*.mat"))) != 0 ]


def map_num_soz(l:int, remap_iz =False)->str:
    """Returns label for each number following this convention
            # 0 - NIZ 
            # 1 - SOZ 
            # 2 - PZ
            # 3 - IZ
    Args:
        l (int): soz designation code found in structs

    Returns:
        str: string designation of 'SOZ', 'NIZ', etc
    """
    match l:
        case 0:
            return "NIZ"
        case 1:
            return "SOZ"
        case 2: 
            return "PZ"
        case 3:
            if remap_iz:
                return "NIZ"
            return "IZ"
def struct_to_pat_df(struct, sub_ids:list[str], filt_dist=0)->pd.DataFrame:
    """Takes the ISH cohort struct and converts to dataframe similar to agg_subjects
    NOTE: this method assumes struct organization for accessing data is the same as 
    found in the PCA project original 81 structs folder
    Args:
        struct (_type_): matlab struct
        sub_ids (list[str]): list of subject ids

    Returns:
        pd.DataFrame: pdc connectivity net in out perf freq, per region
    """
    
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename='conn_dyn.log')
    subj_dfs = []
    for i, pat_obj in enumerate(struct['pats'][0]):
        subj_id = sub_ids[i]
        # map numbers to SOZ IZ etc 
        node_labels = np.array([map_num_soz(l, remap_iz=True) for l in pat_obj['SOZ']])
        
        label_inds = get_reg_inds(node_labels)
        soz_inds, pz_inds, nz_inds = label_inds['soz'], label_inds['pz'], label_inds['nz']
        if pz_inds.shape[0] == 0:
            print(f"Subj {subj_id} has no PZ designation")
            logger.warning(f"Subj {subj_id} has no PZ zone!")
        conn_mat = long_to_mat(pat_obj['long_Z'] )
        if filt_dist >0:
            dist_mat = pat_obj['dist_mat']
            conn_mat = filter_dist(conn_mat, dist_mat,filt_dist)
        pdc_dict = {'eyes_closed_inter':conn_mat}
        df = assemble_net_conn(subj_id, pdc_dict,soz_inds,pz_inds,nz_inds,bands=BANDS[1:])
        subj_dfs.append(df)
    return pd.concat(subj_dfs)

def long_to_mat(long_mat:np.ndarray)-> np.ndarray:
    """convert an array (or list ) of 2D matrices to a proper 3 D matrix

    Args:
        long_mat (_type_): 

    Returns:
        np.ndarray: 3D matrix with N_BANDS x K_NODES x K_Nodes
    """
    k_nodes, _ = long_mat[0][0].shape
    n_bands = long_mat.shape[0]
    conn_mat = np.zeros((n_bands, k_nodes,k_nodes))
    for b in range(n_bands):
        conn_mat[b,:,:] = long_mat[b][0]
    return conn_mat

def conn_to_flow_df(conn_mat:np.ndarray, reg_inds:dict)->pd.DataFrame:
    """Converts a connectivity matrix to a df containing generalized 
    directed connectivity in a source -> target , value format

    Args:
        conn_mat (np.ndarray): NxN connectivity matrix
        reg_inds (dict): description of regions of interest to summarize long (SOZ,PZ, etc)

    Returns:
        pd.DataFrame: dataframe with source to target, and values as a ratio of total connectivity
    """
    flow_targets = [(src, trgt) for src in reg_inds.keys() for trgt in reg_inds.keys()]
    flow_df = pd.DataFrame(columns=['source','target','value'])
    total_conn = np.nansum(conn_mat)
    for i, (src, trgt) in enumerate(flow_targets):
        src_inds = reg_inds[src]
        trgt_inds = reg_inds[trgt]
        src_rows = conn_mat[src_inds,:]
        src_trgt_conns = np.nansum(src_rows[:,trgt_inds])/total_conn * 100 #calc proportion of total connections
        flow_df.loc[i] = [src,trgt,src_trgt_conns]
    return flow_df

def map_subject_to_flow(subj_id, pat_files, label_df, filt_dist=0, **kwargs ):
    """_summary_

    Args:
        subj_id (_type_): _description_
        conn_obj (_type_): _description_
        filt_dist (int, optional): _description_. Defaults to 0.
    """
    flow_dfs = []
    pat_structs = load_structs(pat_files, **kwargs)
    for pat_obj in pat_structs:
        conn_dict, label_inds = prep_conn(subj_id, label_df, pat_obj,filt_dist=filt_dist, **kwargs)
        if label_inds['pz'].shape[0] == 0:
            print(f"Subj {subj_id} has no PZ designation")
            logger.warning(f"Subj {subj_id} has no PZ zone!")
        for b, band in enumerate(BANDS):
            for period in PERIOD:
                conn_mat = conn_dict[period][b,:,:]
                df = conn_to_flow_df(conn_mat, label_inds)
                df['subj'] = subj_id
                df['band'] = band
                df['period'] = period
                df['seizure'] = pat_obj['sz_type']
                flow_dfs.append(df)
    return pd.concat(flow_dfs)

def load_structs(file_list, cores=12)->list[dict]:
    """loading .mat structs is the limiting step in most of these pipelines
    so let's paralellize this. With a 10G port, this should easily speed things up.
    Will slow down if querying ernie from diff network (maybe a VPN)

    Args:
        file_list (_type_): list of .mats to load
        cores (int, optional) Defaults to 12.
    """
    p = Pool(cores)
    return p.map(load_mat,file_list)

def map_cohort_to_flow(subj_ids, folders,label_df,**kwargs):
    flow_df = []
    for i,subj in enumerate(subj_ids):
        sub_files = glob.glob(os.path.join(folders[i], '*mat'))
        sub_label_df = label_df[label_df.subj==subj]
        df  = map_subject_to_flow(subj, sub_files, sub_label_df, **kwargs)
        flow_df.append(df)
    return pd.concat(flow_df)