import logging
import os
import sys
from scipy.io import loadmat
import mat73
import glob
import pdb
import yaml
from multiprocessing import Pool
from functools import partial
import getopt


from collections import Counter, defaultdict
import pandas as pd
import numpy as np
import copy
from pandarallel import pandarallel


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

PERIOD = ['inter','pre','ictal','post']

pandarallel.initialize()


#configure logging
deflogger= logging.getLogger(__name__)
logging.basicConfig(filename='../logs/conn_dyn.log', level=logging.WARNING)
from loguru import logger

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

def get_node_labels(conn_obj):
    """returns the designation of each node in the connectivity matrix, 
    PZ, SOZ, NIZ, etc on a per seizure basis

    Args:
        conn_obj (_type_): 
    """
    soz_designation = read_conn_struct(conn_obj, 'pdc','soz_per_seizure')
    return [map_label(ch) for ch in soz_designation]

def get_chan_names(conn_obj):
    """returns the name of each bipolar channel

    Args:
        conn_obj (_type_): 
    """
    chan_names = read_conn_struct(conn_obj, 'pdc','bip_labels_used')
    return np.array([chan[0] for chan in chan_names])

def get_reg_inds(pat_conn_labels):
    if type(pat_conn_labels) != np.ndarray:
        pat_conn_labels = np.array(pat_conn_labels)
    soz_inds = np.where(pat_conn_labels == 'SOZ')[0]
    pz_inds = np.where(pat_conn_labels == 'PZ')[0]
    nz_bool = pat_conn_labels == 'NIZ'
    iz_bool = pat_conn_labels =='IZ'
    nz_inds = np.where(np.logical_or(nz_bool, iz_bool))[0]
    return {'soz':soz_inds, 'pz' : pz_inds, 'nz': nz_inds}


def get_conn_dict_full(conn_obj,key='pdc', filt_dist=0,**kwargs):
    """Reads the patient seizure struct and returns a connectivity dictionary with keys
    as the window and the values as the connectivity matrix

    Args:
        conn_obj (_type_): matlab struct read in as dictionary 
        key (str, optional): key for reading struct (somewhat depends on how it was saved). Defaults to 'pdc'.
        filt_dist (int, optional): _description_. Defaults to 0.

    Returns:
        _type_: dictionary of win_number -> PDC_matrix, note that window state is not tracked here and should be
        mapped elsewhere 
    """
    #TODO fix final keys
    pdc_all_windowed = read_conn_struct(conn_obj, key, 'PDC_all_windowed') #returns a 6 x N_win x N_ch x N_ch

    if filt_dist > 0:
        dist_mat = read_conn_struct(conn_obj, 'pdc', 'dist_mat')
        pdc_all_windowed = filter_dist(pdc_all_windowed, dist_mat, 20)
    n_win = pdc_all_windowed.shape[1]
    return dict(zip([i for i in range(n_win)],[pdc_all_windowed[:,i,:,:] for i in range(n_win) ] ))

def get_conn_dict_summary(conn_obj,key='pdc', periods=PERIOD, filt_dist=0,**kwargs):
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
    if len(conn_mat.shape) == 4:
        assert conn_mat.shape[-1] == conn_mat.shape[-2], "weird shape order check conn matrix!"
        conn_mat[:,:, d_x, d_y] = np.nan
        return conn_mat
    if len(conn_mat.shape) == 3:
        assert conn_mat.shape[1] > conn_mat.shape[0], "Weird shape order, check connectivity matrix!"
        conn_mat[:, d_x, d_y] = np.nan
        return conn_mat
    conn_mat[d_x, d_y] = np.nan
    return conn_mat

def assemble_net_conn_selector(pdc_dict, soz_inds, pz_inds, nz_inds, chan_names, bands=BANDS, ref_stat=defaultdict(lambda: None), period_meta=defaultdict(lambda :''), verbose=False, **kwargs):
    """Select function used to assemble net connectivity based on verbose value.

    Args:
        subj_id (_type_): string of subject id
        pdc_dict (_type_): dictionary mapping period (int) to pdc (np.array)
        soz_inds (_type_): indices corresponding to soz nodes
        pz_inds (_type_): pz node indices
        nz_inds (_type_): nz node indices

    Returns:
        _type_: Dataframe of next connectivity for each region.
    """

    if verbose:
        # outputs a csv row for every bipole
        net_df = assemble_net_conn_verbose(pdc_dict, soz_inds, pz_inds, nz_inds, chan_names, ref_stat=ref_stat, period_meta=period_meta)
    else:
        # averages SOZ, PZ, and NIZ electrodes within a sz event
        net_df = assemble_net_conn(pdc_dict, soz_inds, pz_inds, nz_inds, ref_stat=ref_stat, period_meta=period_meta)

    return net_df

@logger.catch
def assemble_net_conn_verbose(pdc_dict, soz_inds, pz_inds, nz_inds, chan_names, bands=BANDS, ref_stat=defaultdict(lambda: None), period_meta=defaultdict(lambda :'')):
    """Assemble net directed connectivity across periods of interest for all 
    frequency bands without averaging by region

    Args:
        subj_id (_type_): string of subject id
        pdc_dict (_type_): dictionary mapping period (str) to pdc (np.array)
        soz_inds (_type_): indices corresponding to soz nodes
        pz_inds (_type_): pz node indices
        nz_inds (_type_): nz node indices

    Returns:
        _type_: Dataframe of next connectivity for each region.
    """
    cols = ['period', 'region', 'bip', 'in_pdc', 'out_pdc', 'freq_band', 'window_designations']
    
    net_dict = {}

    for period, pdc in pdc_dict.items():
        period_designations = period_meta[period]
        for b in range(len(bands)):
            #TODO add z_score to interictal period and LOCK the reference!
            mu_in, std_in = ref_stat['mu_in'], ref_stat['std_in']
            z_pdc_in = z_score_conn(pdc[b, :,:],mu=mu_in[b,:],std=std_in[b,:],direction='col')
            mu_out, std_out = ref_stat['mu_out'], ref_stat['std_out']
            z_pdc_out = z_score_conn(pdc[b,:,:], mu=mu_out[:,b], std=std_out[:,b], direction='row')

            band = bands[b]
            
            if len(soz_inds) != 0:
                soz_ins = np.nanmean(z_pdc_out[:,soz_inds], axis=0)
                soz_outs = np.nanmean(z_pdc_in[soz_inds,:], axis=1)
                for soz_in, soz_out, chan in zip(soz_ins, soz_outs, chan_names[soz_inds]):
                    for key,val in zip(cols,[period,'soz',chan,soz_in,soz_out,band,period_designations]):
                        net_dict = update_dict(net_dict,key,val)

            if len(pz_inds) != 0:
                pz_ins = np.nanmean(z_pdc_out[:,pz_inds], axis=0)
                pz_outs = np.nanmean(z_pdc_in[pz_inds,:], axis=1)
                for pz_in, pz_out, chan in zip(pz_ins, pz_outs, chan_names[pz_inds]):
                    for key,val in zip(cols,[period,'pz',chan,pz_in,pz_out,band,period_designations]):
                        net_dict = update_dict(net_dict,key,val)

            nz_ins = np.nanmean(z_pdc_out[:,nz_inds], axis=0)
            nz_outs = np.nanmean(z_pdc_in[nz_inds,:], axis=1)

            for nz_in, nz_out, chan in zip(nz_ins, nz_outs, chan_names[nz_inds]):
                for key,val in zip(cols,[period,'nz',chan,nz_in,nz_out,band,period_designations]):
                    net_dict = update_dict(net_dict,key,val)

    net_df = pd.DataFrame.from_dict(net_dict)
    net_df.insert(loc=5,column='net_pdc',value=net_df.in_pdc - net_df.out_pdc)
    
    return net_df

@logger.catch
def assemble_net_conn(pdc_dict, soz_inds, pz_inds, nz_inds, bands=BANDS, ref_stat=defaultdict(lambda: None), period_meta=defaultdict(lambda :'')):
    """Assemble net directed connectivity across periods of interest for all 
    frequency bands with an average per region (SOZ, PZ, NIZ)

    Args:
        subj_id (_type_): string of subject id
        pdc_dict (_type_): dictionary mapping period (int) to pdc (np.array)
        soz_inds (_type_): indices corresponding to soz nodes
        pz_inds (_type_): pz node indices
        nz_inds (_type_): nz node indices

    Returns:
        _type_: Dataframe of next connectivity for each region.
    """
    cols = ['period', 'region', 'net_pdc', 'in_pdc', 'out_pdc', 'freq_band', 'window_designations']
    
    net_dict = {}

    for period, pdc in pdc_dict.items():
        period_designations = period_meta[period]
        for b in range(len(bands)):
            #TODO add z_score to interictal period and LOCK the reference!
            mu_in, std_in = ref_stat['mu_in'], ref_stat['std_in']
            z_pdc_in = z_score_conn(pdc[b, :,:],mu=mu_in[b,:],std=std_in[b,:],direction='col')
            mu_out, std_out = ref_stat['mu_out'], ref_stat['std_out']
            z_pdc_out = z_score_conn(pdc[b,:,:], mu=mu_out[:,b], std=std_out[:,b], direction='row')

            band = bands[b]
            
            soz_in = np.nanmean(z_pdc_out[:,soz_inds])
            soz_out = np.nanmean(z_pdc_in[soz_inds,:])
            net_soz = soz_in - soz_out

            for key,val in zip(cols,[period,'soz',net_soz,soz_in, soz_out, band,period_designations]):
                net_dict = update_dict(net_dict,key,val)

            if len(pz_inds) != 0:
                pz_in = np.nanmean(z_pdc_out[:,pz_inds])
                pz_out = np.nanmean(z_pdc_in[pz_inds,:])
                net_pz = pz_in - pz_out
                for key,val in zip(cols,[period,'pz',net_pz,pz_in, pz_out, band, period_designations]):
                    net_dict = update_dict(net_dict,key,val)

            nz_in = np.nanmean(z_pdc_out[:,nz_inds])
            nz_out = np.nanmean(z_pdc_in[nz_inds,:])
            net_nz = nz_in - nz_out
            for key,val in zip(cols,[period,'nz',net_nz,nz_in, nz_out, band, period_designations]):
                net_dict = update_dict(net_dict,key,val)

    net_df = pd.DataFrame.from_dict(net_dict)
    
    return net_df

def load_assemble_obj(path, **kwargs):
    """Loads the .mat object specified in path
    and runs the assemble_obj routine
    
    Returns a pd.dataframe or None (if unable to load object) 
    """  
    struct_obj = load_mat(path)
    if struct_obj == None:
        return pd.DataFrame()
    return assemble_obj(struct_obj, **kwargs)

def assemble_obj(conn_obj, wintype='full' ,**kwargs)->pd.DataFrame:
    """For a given subject's connectivity struct, return the dataframe containing net
    connectivity for each frequency band, over relevant periods
    NOTE: this will only work with DIRECTED connectivity

    Args:
        conn_f (str): .mat struct with connectivity matrices
        label_df (df): df of channel labels -> {SOZ_inds,NZ_inds, PZ_inds}

    Returns:
        pd.DataFrame: net connectivity
    """

    conn_dict, label_inds, chan_names = prep_conn(conn_obj, wintype, **kwargs)
    soz_inds, pz_inds, nz_inds = label_inds['soz'], label_inds['pz'], label_inds['nz']

    #get mu and std to Z-score against interictal period
    
    if wintype == 'full':
        window_dict = prep_window_dict(conn_obj)
        win_size = kwargs['win_size'] if "win_size" in kwargs.keys() else 5
        ref_stats= score_period(list(conn_dict.values()), window_dict, agg_win="interictal", win_size=win_size, **kwargs)
        conn_df = assemble_net_conn_selector(conn_dict, soz_inds, pz_inds, nz_inds, chan_names, ref_stat=ref_stats, period_meta=window_dict, **kwargs)
    else:
        conn_df = assemble_net_conn_selector(conn_dict, soz_inds, pz_inds, nz_inds, chan_names, **kwargs)

    conn_df['eventID'] = read_conn_struct(conn_obj, 'pdc', 'eventID')
    conn_df['patID'] = read_conn_struct(conn_obj, 'pdc', 'patID')
    conn_df['sz_type'] = read_conn_struct(conn_obj, 'pdc', 'sz_type')
    return conn_df

def score_period(conn_mats, window_dict, agg_win="interictal", buffer=60, win_size=5, directed=True, stats='full',bands=BANDS,**kwargs):
    """
    score a set of connectivity matrices, filtered by period criterion. 
    Establishes the set of reference values to compare other periods against.
    Often use interictal as baseline to get mu and std to see how other periods compare
    conn_mats should be a list of connectivity matrices within a band.
        Each matrix has shape BAND x N_ch x N_ch
    Returns a dictionary of reference stats for each band of interest
    #in = 'col'
    """
    conn_mats = filter_periods(conn_mats, window_dict, agg_win, buffer, win_size)
    conn_mats = np.array(conn_mats)
    if not directed:
        
        means = [np.nanmean(conn_mats[:,b,:,:]) for b in bands]
        stds = [np.nanstd(conn_mats[:,b,:,:]) for b in bands]
        return {"mu": means,"std":stds}
    else:
        stat_dict = {}
        num_win, _, d,_ = conn_mats.shape
        stat_dict['mu_in'] = np.zeros((len(bands),d)) 
        stat_dict['mu_out'] = np.zeros((d,len(bands)))
        stat_dict['std_in'] = np.ones((len(bands),d))
        stat_dict['std_out'] = np.ones((d,len(bands)))
        for b in range(len(bands)):
            if stats == 'mu' or stats == 'full':
                mu_in = np.array([np.nanmean(conn_mats[i,b,:,:],axis=0) for i in range(num_win)])
                mu_in = np.mean(mu_in, axis=0)
                stat_dict['mu_in'][b,:] = mu_in
                
                mu_out = np.array([np.nanmean(conn_mats[i,b,:,:],axis=1) for i in range(num_win)])
                mu_out = np.mean(mu_out, axis=0)
                stat_dict['mu_out'][:,b] = mu_out
            if stats == 'std' or stats == 'full':
                sigma_in = np.array([np.nanstd(conn_mats[i,b,:,:],axis=0) for i in range(num_win)])
                sigma_in = np.mean(sigma_in, axis=0)
                stat_dict['std_in'][b,:] = sigma_in

                sigma_out = np.array([np.nanstd(conn_mats[i,b,:,:],axis=1) for i in range(num_win)])
                sigma_out = np.mean(sigma_out, axis=0)
                stat_dict['std_out'][:,b] = sigma_out 

        #Additional code in to allow for preserving original pdc values. For example, not removing the mean but just 
        # compressing the range based on the standard deviations, use stats = 'std'       
        return stat_dict

def filter_periods(conn_mats, window_dict, agg_win, buffer, win_size):
    """
    filters through the window designations and returns only connectivity matrices
    that belong to the specified AGG_WIN and are outside of the buffer period
    between window transitions
    NOTE: assumes that dictionary WINDOW_DICT is ordered as in, the periods go from
    interictal -> ictal -> post-ictal
    """
    keys = list(window_dict.keys())
    designations = list(window_dict.values())
    if agg_win == 'interictal':
        transition_win = find_transition(designations, ['0.0_0.0_1.0', '0.0_1.0_1.0'])
        num_buffer_win = buffer // win_size
        if transition_win ==  -1:
            logger.warning(f"Error Finding transition from interictal to ictal state, set of all window designations: {set(designations)}")
            return conn_mats
        filtered_wins = keys[0: transition_win-num_buffer_win]
        return [conn_mats[k] for k in filtered_wins]

    elif agg_win == 'postictal':
        transition_win = find_transition(designations, ['1.0_1.0_2.0', '1.0_2.0_2.0'])
        num_buffer_win = buffer // win_size
        if transition_win ==  -1:
            logger.warning(f"Error Finding transition from ictal to post_ictal state, set of all window designations: {set(designations)}")
            return conn_mats
        filtered_wins = keys[transition_win+num_buffer_win:]
        return [conn_mats[k] for k in filtered_wins]    
    return conn_mats    

def prep_window_dict(conn_obj):
    """Given a patient struct, returns a dictionary of window_designations mapped
    keyed by window number. Should allow for finding transistions in seizure onset and 
    centering the seizure period.

    Args:
        conn_obj (_type_): loaded matlab struct from peri-ictal dataset. 

    Returns:
        _type_: dictionary mapping period_number -> window_designation.
    """
    w_end = read_conn_struct(conn_obj, 'pdc', 'window_end_state')
    w_start = read_conn_struct(conn_obj, 'pdc', 'window_start_state')
    w_mid = read_conn_struct(conn_obj, 'pdc' ,'window_middle_state')
    num_win = len(w_mid)
    keys = [i for i in range(num_win)]
    designations = [ f"{w_start[i]}_{w_mid[i]}_{w_end[i]}" for i in range(num_win)] 
    return dict(zip(keys, designations))

def prep_conn( conn_obj, wintype='full', **kwargs):
    if conn_obj == None:
        raise ValueError("Issue with load_mat")
    if wintype =='full':
        conn_dict = get_conn_dict_full(conn_obj, **kwargs)
    else:
        conn_dict = get_conn_dict_summary(conn_obj, **kwargs)
    pat_conn_labels = get_node_labels(conn_obj)
    label_inds = get_reg_inds(pat_conn_labels)
    chan_names = get_chan_names(conn_obj)
    return conn_dict, label_inds, chan_names


def get_regions( conn_obj):
    regions = [reg[0] for reg in read_conn_struct(conn_obj, 'pdc', 'bip_labels_used')]
    
    # regions = format_bipoles(regions)
    # pat_conn_labels = get_pat_conn_labels(label_df, regions,subj_id)
    pat_conn_labels = get_reg_inds(regions)
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
    
    # logger = logging.getLogger(__name__)
    # logging.basicConfig(filename='conn_dyn.log')
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

def conn_to_flow_dict_selector(conn_mat, reg_inds, window_dict, period, chan_names, band, flow_dict, verbose=False, **kwargs):
    """Selects function for creating flow dictionary

    Args:
        conn_mat (np.ndarray): NxN connectivity matrix
        reg_inds (dict): description of regions of interest to summarize long (SOZ,PZ, etc)
        window_dict (dict): window designations associated with each period
        period (int): time period
        flow_dict (dict): master dictionary to update

    Returns:
        dict: keys include source, target, mean connectivity, period, and window_designation
    """

    if verbose:
        # outputs a csv row for every bipole
        flow_dict = conn_to_flow_dict_verbose(conn_mat, reg_inds, window_dict, period, chan_names, band, flow_dict)
    else:
        # averages SOZ, PZ, and NIZ electrodes within a sz event
        flow_dict = conn_to_flow_dict(conn_mat, reg_inds, window_dict, period, band, flow_dict)

    return flow_dict

def conn_to_flow_dict(conn_mat:np.ndarray, reg_inds:dict, window_dict, period, band, flow_dict)->dict:
    """Converts a connectivity matrix to a dict containing generalized 
    directed connectivity in a source -> target , value format per region

    Args:
        conn_mat (np.ndarray): NxN connectivity matrix
        reg_inds (dict): description of regions of interest to summarize long (SOZ,PZ, etc)
        window_dict (dict): window designations associated with each period
        period (int): time period
        flow_dict (dict): master dictionary to update

    Returns:
        dict: keys include source, target, mean connectivity, period, and window_designation
    """
    # only includes categories with at least one electrode
    flow_targets = [(src, trgt) for src in reg_inds.keys() for trgt in reg_inds.keys() if (reg_inds[src].size & reg_inds[trgt].size)]
    keys = ['source','target','value','period','window_designations', 'freq_band']

    for (src, trgt) in flow_targets:
        src_inds = reg_inds[src]
        trgt_inds = reg_inds[trgt]

        src_rows = conn_mat[src_inds,:]
        src_trgt_conns = np.nanmean(src_rows[:,trgt_inds])

        for key, val in zip(keys,(src,trgt,src_trgt_conns,period,window_dict[period],band)):
            flow_dict = update_dict(flow_dict, key, val)

    return flow_dict

def conn_to_flow_dict_verbose(conn_mat:np.ndarray, reg_inds:dict, window_dict, period, chan_names, band, flow_dict)->dict:
    """Converts a connectivity matrix to a dict containing generalized 
    directed connectivity in a source -> target , value format per bipolar channel

    Args:
        conn_mat (np.ndarray): NxN connectivity matrix
        reg_inds (dict): description of regions of interest to summarize long (SOZ,PZ, etc)
        window_dict (dict): window designations associated with each period
        period (int): time period
        flow_dict (dict): master dictionary to update

    Returns:
        dict: keys include source, target, mean connectivity, period, and window_designation
    """
    # only includes categories with at least one electrode
    flow_targets = [(src, trgt) for src in reg_inds.keys() for trgt in reg_inds.keys() if (reg_inds[src].size & reg_inds[trgt].size)]
    keys = ['source','target','src_bip','value','period','window_designations','freq_band']

    for (src, trgt) in flow_targets:
        src_inds = reg_inds[src]
        trgt_inds = reg_inds[trgt]

        for src_ind in src_inds:
            src_rows = conn_mat[src_ind,:]
            src_trgt_conns = np.nanmean(src_rows[trgt_inds])
            src_bip = chan_names[src_ind]
            for key, val in zip(keys,(src,trgt,src_bip,src_trgt_conns,period,window_dict[period],band)):
                flow_dict = update_dict(flow_dict, key, val)

    return flow_dict

def gen_flow_dfs(pathout, inp_paths, num_cores =1, **kwargs)->None:
    """generates a tri-node summary or a given subjects connectivity matrices
    across all periods. 

    Args:
        pathout (str): folder to save csv's
        inp_paths (list[str]): list of .mat files for a given subject

        NOTE: for now this assumes that all files belong to the same subject!
    """
    subj_dfs = []
    if num_cores == 1:
        for i,inp_path in enumerate(inp_paths):
            df  = map_subject_to_flow(inp_path, **kwargs)
            subj_dfs.append(df)
    else:
        logger.info("Running on parallel track")
        p = Pool(min(num_cores, len(inp_paths)))
        subj_dfs = p.imap(partial(map_subject_to_flow, **kwargs), inp_paths)
        p.close()
    
    subj_dfs = pd.concat(subj_dfs)
    subj = subj_dfs.patID.values[0]
    if kwargs['verbose']:
        fname = f"peri_ictal_flow_verbose_{subj}.csv"
    else:
        fname = f"peri_ictal_flow_{subj}.csv"
    
    subj_dfs.to_csv(os.path.join(pathout, fname), index=False)
    logger.success(f'Saving {subj} 3-node view of {len(inp_paths)} seizures to {pathout} as {fname}')

def map_subject_to_flow(pat_file, filt_dist=0, **kwargs):
    """_summary_

    Args:
        subj_id (_type_): _description_
        conn_obj (_type_): _description_
        filt_dist (int, optional): _description_. Defaults to 0.
    """

    flow_dict = {}
    pat_obj = load_mat(pat_file)

    subj_id = read_conn_struct(pat_obj, 'pdc', 'patID')
    # NOTE: need to add functionality for retaining all channels
    pdc_dict, label_inds , chan_names = prep_conn(pat_obj, wintype='full',filt_dist=filt_dist, **kwargs)
    window_dict = prep_window_dict(pat_obj)
    if label_inds['pz'].shape[0] == 0:
        print(f"Subj {subj_id} has no PZ designation for the following struct {pat_file}")
        logger.warning(f"Subj {subj_id} has no PZ channels!")
    if label_inds['soz'].shape[0] == 0:
        print(f"Subj {subj_id} has no SOZ designation for the following struct {pat_file}")
        logger.warning(f"Subj {subj_id} has no SOZ channels!")
    
    for period, pdc in pdc_dict.items():
        for i, band in enumerate(BANDS):
            conn_mat = pdc[i,:,:]
            flow_dict = conn_to_flow_dict_selector(conn_mat, label_inds, window_dict, period, chan_names, band, flow_dict, **kwargs)
    
    flow_df = pd.DataFrame.from_dict(flow_dict)
    flow_df['eventID'] = read_conn_struct(pat_obj, 'pdc', 'eventID')
    flow_df['patID'] = subj_id
    flow_df['sz_type'] = read_conn_struct(pat_obj, 'pdc', 'sz_type')
    return flow_df
    

def load_structs(file_list, cores=12, **kwargs)->list[dict]:
    
    """loading .mat structs is the limiting step in most of these pipelines
    so let's paralellize this. With a 10G port, this should easily speed things up.
    Will slow down if querying ernie from diff network (maybe a VPN)

    Args:
        file_list (_type_): list of .mats to load
        cores (int, optional) Defaults to 12.

    NOTE: if pipeline is going to use parlalleization after load step, it may may sense to load
    lazily. 
    WARNING: any paralellization used after this point will copy working memory into each workers process memmory
    so if you just loaded 10gb of data and each worker onlly needs 1gb, then this will unnecessarily copy 9gb
    per worker
    """
    p = Pool(cores)
    return p.map(load_mat,file_list)


def assemble_peri_obj_para(struct_paths:list[str], cores =12, **kwargs):
    """Assemble connectivity in parallelized fashion"""
    p = Pool(cores)

    dfs =  p.imap(partial(load_assemble_obj, **kwargs), struct_paths)
    dfs = [df for df in dfs if not df.empty]
    return pd.concat(dfs)

def center_onset(peri_df: pd.DataFrame, win_size=5, stride=1, center_designations=["0.0_0.0_1.0", "0.0_1.0_1.0"], **kwargs)->pd.DataFrame:
    """Given a dataframe with each row matching net connectivity per period
    and a string representation of the beginning_middle_end of each window designation, 
    this method returns a modified dataframe with a "win_sz_centered" column that counts up to 
    the seizure transition. The new window counts will have transitions at zero. This
    will allow you to align all peri-ictal dynamics according to the transition point of seizures
    This method will also add a "period_label" column with the following labels
    "Interictal, pre_ictal, transition_into_sz, ictal, transition out, post_ictal"

    Pre_ictal will be the minute before a seizure and post_ictal will be the minute after, assuming
    WIN_SIZE is in seconds

    Args:
        peri_df (pd.DataFrame): dataframe from assemble_obj subroutine. should have net/in/out connectivity
        across all periods 
        win_size (int) : length of window size in seconds
        center_designations (list, optional): _description_. Defaults to ["0_0_1", "0_1_1"].

    Returns:
        pd.DataFrame: df of peri-ictal connectivity with period_count centered at 
        the transition point into the ictal state
    """
    centered_dfs = []
    subj = peri_df.patID.values[0]
    peri_df.eventID = peri_df.eventID.astype(str) # fixes bug where some patients have 3 and "3" as separate events for example
    for event in set(peri_df.eventID): #NOTE: could this be a groupby apply?
        event_df = peri_df[peri_df.eventID == event]
        try:
            centered_event_df = center_event_df(win_size, stride, center_designations, event_df, **kwargs)
            centered_dfs.append(centered_event_df)
        except IndexError as e:
             logger.warning(f"Issue centering {subj} on event: {event}.\nMore details: {e}")
        except ValueError as e:
             logger.warning(f"Issue centering {subj} on event: {event}.\nMore details: {e}")
    
    return pd.concat(centered_dfs)

def center_event_df(win_size, stride, center_designations, event_df, start_buffer=10, end_buffer=10, mid_sz_length=10)-> pd.DataFrame:
    """center an event_level data frame: see center_onset above

    Args:
        win_size (_type_): _description_
        stride (_type_): _description_
        center_designations (_type_): _description_
        event_df (_type_): _description_

    Returns:
        pd.DataFrame: _description_
    """
    ##This bit of code will help to reduce the redundancy of the event_df
    ## We will take advantage of the fact that dictionaries
    ## are guaranteed to be ordered in order to remove redundant rows
    # 
    p_d_dict = dict(zip(event_df.period, event_df.window_designations))
    periods, window_designations = np.array(list(p_d_dict.keys())), np.array(list(p_d_dict.values()))
    centered_dict = center_windows(window_designations, periods,center_designations=center_designations)
    centered_windows = [centered_dict[period] for period in event_df.period.values]
    event_df.insert(loc =0, column='win_sz_centered',value = centered_windows)
    if event_df.win_sz_centered.isna().all():
        event_df.insert(loc=0,column='sz_end',value= np.nan)
        event_df.insert(loc=0,column='win_sz_st_end', value = np.nan)
        event_df.insert(loc=0,column='win_label', value = np.nan)
    else:
        sz_end =  get_sz_end(event_df)
        event_df.insert(loc=0,column='sz_end', value =sz_end)
        sz_st_end = sample_seizures(event_df, start_buffer=start_buffer, end_buffer=end_buffer, mid_sz_length=mid_sz_length, win_size=win_size)
        event_df.insert(loc=0,column='win_sz_st_end', value = sz_st_end)
        labelled_window = event_df.parallel_apply(label_window, args=[win_size, stride], axis=1)
        event_df.insert(loc=0,column='win_label', value = labelled_window)
    centered_event_df = event_df
    return centered_event_df

def sample_seizures(peri_df, start_buffer, end_buffer, mid_sz_length, win_size=5, stride=1):
    """
    Returns index tracking all windows with zero centered at seizure onset and all seizures
    ending at the same time. Seizure window is determined by length of start_buffer + end_buffer + mid_sz_length
    Aligns all seizures in time while preserving transition windows into seizure 
    (start_buffer) and windows out of seizure (end_buffer. Seizures that do not 
    have the minumum number of windows for start and end will have np.nan values returned 
    in their final window counts.

    NOTE: buffers and mid_sz_length are in seconds and win_size is number of seconds per window
    """
    assert len(set(peri_df.eventID)) == 1, "Can only pass one seizure per call to get_sz_end!"
    assert 'win_sz_centered' in peri_df.columns, "Need to center seizures on start before sampling seizure windows!"
    assert "window_designations" in peri_df.columns, "Unable to track seizure state without window designations!"
    assert "sz_end" in peri_df.columns, "need to find end of szr first!"

    start_buffer =  start_buffer//stride
    end_buffer = end_buffer//stride
    mid_sz_length = mid_sz_length//stride
    # Get window number that seizure ends in
    sz_end = peri_df.sz_end.values[0]

    # Buffers define the minimally acceptable seizure length
    if sz_end < (start_buffer + mid_sz_length + end_buffer):
        return [np.nan for _ in range(peri_df.shape[0])]

    # pull out all windows leading up to (or after) seizure and the transition windows that we will preserve
    all_wins = np.unique(peri_df.win_sz_centered)
    pre_wins = np.unique(peri_df[peri_df.win_sz_centered < start_buffer].win_sz_centered)
    post_wins = np.unique(peri_df[peri_df.win_sz_centered > sz_end-end_buffer].win_sz_centered)
    mid_wins = np.setdiff1d(all_wins, pre_wins) 
    mid_wins = np.setdiff1d(mid_wins, post_wins)

    # resample middle windows such that all windows are compressed into mid_sz_length
    # example: 10-30 with mid_sz_length 5 is mapped to 10,10,10,10,11,11,11,11,etc.
    mid_sample = np.array([i+start_buffer for i, arr in enumerate(np.array_split(mid_wins,mid_sz_length)) for _ in range(arr.size)])

    # shift post_wins back to begin at start_buffer + mid_sz_length
    post_resample = (post_wins - sz_end + end_buffer - 1 + start_buffer + mid_sz_length)
    resampled_wins = np.hstack((pre_wins, mid_sample, post_resample))
    remap_periods = dict(zip(all_wins,resampled_wins))
    resampled_periods = peri_df.win_sz_centered.parallel_apply(lambda x: remap_periods[x])

    return resampled_periods
    
def get_sz_end(peri_df, use_col='win_sz_centered'):
    """Given a calculated peri-ictal connectivity df, returns the period to reference 
    as seizure termination

    Args:
        peri_df (pd.DataFrame): compiles ,peri_connectivity dataframe. Assumes only one event (seizure ) per df 
        use_col (str) :  column to use as time column, may be centered based on seizure (win_sz_centered) or the 
                        same period defined from original aggregation
    Returns
        sz_end time as int
    """
    assert len(set(peri_df.eventID)) == 1, "Can only pass one seizure per call to get_sz_end!"
    assert use_col in peri_df.columns, f"column used to index time is not in df! {use_col}"
    assert "window_designations" in peri_df.columns, "Unable to track seizure state without window designations!"

    sz_end = np.where("2.0_2.0_2.0" == peri_df.window_designations)[0]
    sz_end = sz_end[0]
    return peri_df[use_col].values[sz_end]


def label_window(df_row, win_size, stride=1):
    """Labels a row of the data frame using the window_designation column
    and the win_sz_centerd column to assign states. 

    Current logic:
            Interictal - > 1 min prior to a seizure
            Pre-ictal - 1 min prior to the seizure transition (uses win_sz_centered)
            Early Ictal - all labels corresponding to 0_0_1 up through 0_1_1 that occur after 
                        the designated transition (win_sz_centerd > 0)
            Ictal - all periods with the 1_1_1 designation
            Late ictal - all periods with 1_1_2 up to 1_2_2 
            Early post-ictal - all periods with designation 2_2_2 occuring in first minute after end of seizure
            Post ictal - all periods designated 2_2_2 after the minute mark of the sezure end

    NOTE: in the future can further subdivide the ictal states
    Args:
        df_row (_type_): _description_
        win_size (int) : length of window size in seconds, important for calculating pre- early post-ictal times
        sz_end (int) : period where seizure ends 
            
    """
    win_designation =    df_row['window_designations']
    centered_period = df_row['win_sz_centered']
    sz_end = df_row['sz_end']
    match win_designation:
        case "0.0_0.0_0.0":
            if centered_period * stride < -60:
                return "interictal"
            return "pre_ictal"
        case "0.0_0.0_1.0" | "0.0_1.0_1.0":
            if centered_period >= 0:
                return "early_ictal"
            return "pre_ictal"
        case "1.0_1.0_1.0": 
            return "ictal"
        case "1.0_1.0_2.0" | "1.0_2.0_2.0":
            return "late_ictal"
        case "2.0_2.0_2.0":
            #NOTE: magic number at 60, in future make early post ictal period designation more modular
            if (centered_period - sz_end)* stride <= 60:
                return "early_post_ictal" 
            return "post_ictal"
    return 

def center_windows(window_designations, periods, center_designations=["0.0_0.0_1.0", "0.0_1.0_1.0"]):
    """Returns an array of window counts with the 0th window
    occuring at the transition. If non of the windows in CENTER_DESIGNATION are found, then method 
    will default to defining the transition window as the first window containing 1's after a window of all "0_0_0"
    NOTE: this method assumes only one "seizure" per series
    Args:
        window_designations (_type_): series of strings labelling the beginning_middle_end of a window.
                                        0 -> before a seizure
                                        1 - > during seizure
                                        2 - > after seizure
        periods : should track the period 
        center_designations (list): defaults to 2 window types to center at. If the first window type is not available
        will default to second.ljl
    Returns:
        dictionary: original period number mapped to count of windows with zero marking the transition into the seizure state. THis is important for aligning across groups
    """
    transition_ind= find_transition(window_designations, center_designations)
    if transition_ind == -1:
        centered_wins = [np.nan for _ in periods]
        logger.warning(f"Problem with finding transition")
    else:
        trans_period = periods[transition_ind]
        centered_wins = periods - trans_period
    return dict(zip(periods, centered_wins))

def find_transition(window_designations, center_designations):
    """finds the point in an ORDERED series of windows where window designation changes 
    from one state to another. For example most commonly used to find the end of the interictal period
    to the beginning of the ictal period

    Args:
        window_designations (_type_): ordered list of string designations that correspond to the state of the recording. 
            often a 3 element string with the first designating window start, mid-window, and third is end 
            of window label                        
        center_designations (_type_): transition window

    Returns:
        int : index where window types change
    """
    center_transitions= copy.copy(center_designations)
    if type(window_designations) != np.ndarray:
        window_designations = np.array([w for w in window_designations])
    try:
        transition = np.where(window_designations == center_transitions[0])[0]
    except IndexError:
        pdb.set_trace()
    if len(transition) == 0 and len(center_transitions) > 1:
        while len(transition) ==0 and len(center_transitions) >=1:
            transition = np.where(window_designations == center_transitions[0])[0]
            center_transitions.pop(0)
    try:
        assert len(transition) != 0, f"Check center_designation {center_designations} and df! No transitions found"
    except AssertionError:
        if len(set(window_designations)) < 3:
            return -1
    return transition[0]

def peri_net_pipeline(pathout, paths, num_cores=16, **kwargs):
    count = 0
    peri_dfs = []
    #NOTE THIS assumes that the subject list is the same for all runs
    for path_chunk in chunker(paths, num_cores):
        conn_df = assemble_peri_obj_para(path_chunk, num_cores, **kwargs)
        peri_dfs.append(conn_df)
        count += len(path_chunk)
        logger.success(f"Finished  {count} seizures",)
    peri_dfs = pd.concat(peri_dfs)
    subj = peri_dfs.patID.values[0]
    assert len(set(peri_dfs.patID)) ==1, "mixing subjects, For now that's bad!"

    if kwargs['verbose']:
        fname = f"peri_ictal_network_verbose_{subj}.csv"
    else:
        fname = f"peri_ictal_network_{subj}.csv"
    peri_dfs.to_csv(os.path.join(pathout, fname),index=False)
    logger.success(f"Saving {subj} net periconnectivity in {pathout} as {fname} ")



def main(argv):
    opts, _ = getopt.getopt(argv,"d:p:c:l:",["datadir=",'pathout=','config=','logdir='])
    for opt, arg in opts:
        if opt in ("-d", 'datadir'):
            datadir = arg
        elif opt in ('-p', '--pathout'):
            pathout = arg
        elif opt in ("-c", '--config'):
            config_f = arg
        elif opt in ("-l", "--logdir"):
            logdir = arg
    

    logger.add(os.path.join(logdir, "connectivity_dynamics_run.log"), enqueue=True,level=40)
    with open(config_f, 'r') as f:
            config =  yaml.safe_load(f)

    kwargs = config['peri_para']
    print(kwargs)
    pipeline = "net" if 'pipeline' not in config.keys() else config['pipeline']

    paths = glob.glob(os.path.join(datadir, "*PDC.mat"))
    match pipeline:
        case 'net':
            peri_net_pipeline(pathout, paths, **kwargs)
        case 'trinode':
            gen_flow_dfs(pathout, paths, **kwargs)
        case "power_spectral":
            return NotImplementedError("Need to Implement Straight UP PSD")
        case "e_i_bal":
            return NotImplementedError("have not yet implemented The E/I Autocorr")
        
    


if __name__ == "__main__":
    with logger.catch():
        main(sys.argv[1:])
