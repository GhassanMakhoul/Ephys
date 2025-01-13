
#Systems I/O
import os
from pathos.pools import ProcessPool
import glob
from scipy.io import loadmat
import mat73
import re
import h5py
import pdb
from loguru import logger
#data/stats
import numpy as np
import pandas as pd
from scipy.stats import f_oneway, ttest_ind
from collections import Counter
np.random.seed(10555)
#specialty
"""Originally coded during the dynamic ISH project. Should continue updating and reusing over all projects in the Ephys repository
Author: Ghassan S. Makhoul 
Last modified: 1/10/2025
"""

#GLOBAL Variables
BANDS = ['delta', 'theta','alpha', 'beta', 'gamma_low','gamma_high']
DIR = '/mnt/ernie_main/000_Data/SEEG/SEEG_EyesClosed_RestingState/results/Graham_81pats/PDC_RestingState/'
DTYPE = h5py.special_dtype(vlen=str)
COLOR_MAP= {'pz':'#D95319', 'soz': "#A2142F", "nz":"#0072BD", "niz":"#0072BD", 
            'pz_soz': '#D95319', 'soz_soz': "#A2142F", "nz_soz":"#0072BD"}
COLOR_MAP['niz_soz'] = COLOR_MAP['soz']
COLOR_MAP['nz_soz'] = COLOR_MAP['soz']
COLOR_MAP['niz_pz'] = COLOR_MAP['pz']
COLOR_MAP['nz_pz'] = COLOR_MAP['pz']
COLOR_MAP['niz_niz'] = COLOR_MAP['nz']
COLOR_MAP['nz_nz'] = COLOR_MAP['nz']
FLOWMAP  = {'nz_soz_False':"#00bdaa", 'pz_soz_False':"#bdb400", 'pz_nz_True':"#00bd58", 
            'nz_pz_False':"#9700bd", 'pz_nz_False':"#00bd2c", 'nz_nz_True':"#9400bd", 'soz_pz_False':"#bd6e00",
            'pz_pz_False':"#d98919", 'pz_soz_True': "#bdb000", 'soz_soz_False':"#bd0078", 'soz_nz_True':"#bd6e00",
            'nz_soz_True': "#0072BD", 'nz_nz_False' : "#2000bd", 'soz_soz_True':"#A2142F", 'nz_pz_True':"#00bd09", 
            'pz_pz_True':"#D95319", 'soz_pz_True':"#bd9400", 'soz_nz_False':"#00bd45"}

ctx_keys = {'nz_nz_True_outCTX', 'nz_soz_False_outCTX', 'nz_soz_False_ctx', 'nz_pz_False_ctx', 
            'nz_pz_True_ctx', 'nz_pz_False_outCTX', 'nz_soz_True_ctx', 'nz_nz_False_outCTX',
            'nz_nz_True_ctx', 'nz_nz_False_ctx', 'nz_soz_True_outCTX', 'nz_pz_True_outCTX',
            'soz_nz_False_outCTX', 'soz_pz_False_ctx', 'soz_pz_True_ctx', 'soz_nz_False_ctx', 
            'soz_nz_True_ctx', 'soz_soz_False_outCTX', 'soz_soz_False_ctx', 'soz_nz_True_outCTX',
            'soz_soz_True_ctx', 'soz_pz_True_outCTX', 'soz_soz_True_outCTX', 'soz_pz_False_outCTX'}


#nz_nz_True:"#9400bd"
#nz_nz_False : "#2000bd"


def map_p(p_value) -> 'str':
    """Returns a p_value text mapping to convey level of significant 
                - "#" for  p < .0005
                - "+" for p < .005
                - "*" for p < .05
    Args:
        p_value (_type_): _description_

    Returns:
        str: _description_
    """
    if p_value < .0005:
        return "#"
    if p_value < .005:
        return "+"
    if p_value <= .05:
        return "*"
    return "NS"
    

for key in ctx_keys:
    og_key = "_".join(key.split("_")[0:3])
    hex_color = FLOWMAP[og_key]
    if "outCTX" in key:
        hex_color = hex_color[0:-2] + "00"
    FLOWMAP[key] = hex_color
keys = list(COLOR_MAP.keys())
for key in keys:
    COLOR_MAP[key.upper()] = COLOR_MAP[key]

COLOR_MAP['SOZ_False'] = "#c28893"
COLOR_MAP['SOZ_True'] = COLOR_MAP['SOZ']

COLOR_MAP['PZ_False'] = "#dba38a"
COLOR_MAP['PZ_True'] = COLOR_MAP['PZ']

COLOR_MAP['NIZ_False'] = "#84a2b5"
COLOR_MAP['NIZ_True'] = COLOR_MAP['NIZ']

BANDS = ['delta', 'theta', 'alpha', 'beta','gamma_l', 'gamma_H']
#NOTE: the distinctino between gamma low and gamma high is kind of arbitrary
BAND_RANGES = dict(zip(BANDS, [(0.01, 4), (4, 8), (8, 12), (12, 30),(30,60),(60,120) ]))

def load_mat(f):
    """Loads structs and attmpts to use scipy's loadmat functionality or 
    mat73. Some structs from older matlab work with scipy while others need mat73
    """
    try:
        return loadmat(f)
    except NotImplementedError:
        try:
            return mat73.loadmat(f)
        except:
            print(f"Problem Loading {f}")
            return None
        

def setup_logs(logger, pathout: str, runID: str):
    """Set up logging, adds dir if not exist, and adds sink

    Args:
        pathout (str): DATA_DIR/subj_stim_ma/
    """
    logdir = os.path.join(pathout, 'logs/')
    if not os.isdir(logdir):
        os.mkdir(logdir)
    logf = os.path.join(logdir, f"run_{runID}.log")
    logger.add(logf)
    return logger

def split_bipole(bip_df: pd.DataFrame):
    """splits the bipole column of a bipole df
    duplicates rows

    """
    assert 'bipole' in bip_df.columns, "Need bipole column!"

    contact1 = bip_df.bipole.apply(lambda x: x.split("-")[0].strip(" "))
    contact2 = bip_df.bipole.apply(lambda x: x.split("-")[1].strip(" "))
# python resultAggregator.py -s 'Epat38' -p '/mnt/ernie_main/Ghassan/ephys/data/Epat38'
    df2 = bip_df.copy(deep=True)

    bip_df['contact'] = contact1
    df2['contact'] = contact2


    return pd.concat([bip_df,df2])


def resample_bipole_df(verbose_df :pd.DataFrame, subgroup_col: str, bal_col : str, resamp_schema='balance',**kwargs):
    """Given a dataframe with imbalanced samples, returns a resampled AND summarized dataframe, accounting
    for imbalances in subgroups. For example, if the dataframe is a contact level net connectivity dataframe, there
    may be 200 contacts total and 10 of them may be SOZ's vs 190 NIZ's. A resampled df for boostrapping may then 
    balance SOZ's and NIZ's such that there are 10 NIZ's and 10 SOZ's. Then the aggregator would generate summary
    statistics (mean, variance, median, etc.) over this balanced dataframe. A proper boostrapping pipeline will
    call this method many many times in order to proper sample the dataset.

    Args:
        verbose_df (pd.DataFrame): dataframe with imbalanced subgroups, rows may contain repeat measures which are accounted for
        in the repeat measures keyword arg. Currently pipeline works for one event at a time. 
        subgroup_col (str): column to balance dataset against
        bal_col (str) : column that should be balanced (e.g. bipole level balancing)
        unique_measures (list[str], optional) :  List of columns of different measures
            for single data point , most often frequency bands. Used to make a UID per row.
              Defaults to '', indicating no repeat measures .
        resamp_schema (str, optional): resampling schema defaults to perfectly balancing the dataset
                    also allowed to send in proportions dictionary. Defaults to 'balance'.
    """
    # def concat_entries(row, cooncat_cols):
    #     str_cols = [str(row[col]) for col in cooncat_cols]
    #     return "_".join(str_cols)
    # uids = verbose_df.paralell_apply(lambda x: concat_entries(x, unique_measures), axis=1)
    assert 'bip' in verbose_df.columns, "Need bip "
    assert 'freq_band' in verbose_df.columns, "Need freq band for now"
    assert 'period' in verbose_df.columns, "Need periods of time for now"
    #NOTE: these specific requirements for certain columns help to narrow down the dataframe to just 
    # one set of the bipoles. This is important because we want to know the sampling schema and pull out 
    # a balanced set of bipoles per event. TODO: in the future add modularity such that you can pass in 
    # most df's and rebalance no matter the format. 
    resamp_dfs = []
    
    
    for event in set(verbose_df.eventID):

        event_df = verbose_df[verbose_df.eventID == event]
        sub_df = event_df[event_df.freq_band == 'alpha']
        sub_df = sub_df[sub_df.period==0]
        n_samp = np.min(sub_df.groupby('region').count())
        sub_df = sub_df.groupby(subgroup_col).sample(n=n_samp, replace=True)
        event_df = event_df[event_df[bal_col].isin(sub_df[bal_col])]
        resamp_dfs.append(event_df)
    return pd.concat(resamp_dfs)

def agg_verbose_df(verbose_df: pd.DataFrame, measure_cols:list, categorical_cols:list[str], **kwargs)->pd.DataFrame:
    """Given verbose dataframe with many repeated measures (bipole level and freq band level), returns an aggregated dataframe
    along the groups of interest (specified in categorical_cols)

    Args:
        verbose_df (pd.DataFrame): long form dataframe containing entries from connectivity or eibal pipelines. may be at various
        levels of verbosity. For example may contain a row for every Bipole X frequency_band
        measure_cols (list): Numeric columns to apply mean to (connectivity measure, ei_bal, etc)
        categorical_cols (list[str]): metadata to preserve in agged df, could be period designation, freq band etc
    """
    v_df = verbose_df[measure_cols+categorical_cols]
    v_df = v_df.groupby(categorical_cols).mean().reset_index()
    return v_df


def map_label(label):
    label = int(label)
    match label:
        case 0:
            return "NIZ"
        case 1:
            return 'SOZ'
        case 2:
            return 'PZ'
        case 3:
            return 'IZ'
        
#TODO refactor subsample

# def subsample_df(df:pd.DataFrame, factor:int, random_state=42, groups=['win_label'])->pd.DataFrame:
#     """Returns a subsampled dataframe that selects samples along groupings 
#     specified by GROUPS arg.

#     Args:
#         df (pd.DataFrame): dataframe to resample, often a peri-ictal verbose df
#         factor (int): amount to refactor sample by 2 -> return half the sample, 5-> 1/5 samples returned
#         groups (list, optional): stratifications to segment df by,. Defaults to ['win_label'].

#     Returns:
#         pd.DataFrame: returns subsampled df
#     """

#     group_df = df.groupby(by=groups)
#     count_df= group_df.count()
#     min_grp_size = min(count_df.values)
#     resamp_df = group_df.sample(min_grp_size//factor, random_state=random_state)
#     return resamp_df

def format_soz(soz_labels):
    return [map_label(l) for l in soz_labels]

def format_bipole(bipole):
    bipole = bipole.strip()
    bipole = bipole.replace("-","")
    p = "[0-9][a-zA-Z]"
    bip_match = re.search(p, bipole)

    assert bip_match != None, f"Bipole {bipole}, is not captures by regex pattern {p}"
    s, _ = bip_match.span()
    return bipole[0:s+1] + "-" + bipole[s+1:]

def format_bipoles(char_list):
    return [format_bipole(c) for c in char_list]

def paried_region_significance(region_vals:dict, permutations=10000)-> np.ndarray:
    """Given a dictionary mapping region name to connectivity values
    returns the pairwise N_reg x N_reg significance matrix

    Args:
        region_vals (dict): dictionary with keys as region names and values as region
        connectivity values

        permutations : number of times to run permutation over dataset for significance test
                        recommend at least 10,0000 or NONE

    Returns:
        list[str], np.ndarray: region index names, N_reg x N _reg significance matrix
        making all pairwise comparisons
    """
    reg_keys = list(region_vals.keys())
    n = len(reg_keys)
    sig_mat = np.ones((n,n))
    for i in range(n):
        dist_a = region_vals[reg_keys[i]]
        for j in range(n):
            if i ==j:
                continue
            dist_b = region_vals[reg_keys[j]]
            sig_mat[i,j] = ttest_ind(dist_a, dist_b, permutations=permutations)[1]
    return reg_keys, sig_mat


def load_pdc(pdc_struct):
    """Loads pdc from the .mat struct that in SEEG, eyes closed resting

    
    NOTE: OG struct organization was as follows 
    Index to selects col in struct
                            |
    Index to select patients  |
                        |  |  
                        |  |  
                        |  |  
    epat_26_alpha = tst[0][0][2][0][0][2]
    Columns in Struct
    0 - subID
    1 - labels
    2 - "long" -> 6 x N x N PDC matrix (N= number of bipoles)
    3 - pat_ID_clean - NOTE: use this one instead of labels! 
    4 - SOZ, labels for each bipole -> 0 NIZ, 1 SOZ, 2 PZ, 3 NIZ
    5 - AVG_SOZ  - ??
    6 - long_Z - ?
   
    Args:
        pdc_struct (list of arrays?): lists of arrays from matlab struct

    Returns:
        list: PDC adjacency matrices as list, should be 6, 1 per band
    """
    pdc_all = [pdc for pdc in pdc_struct[0][0]]
    return dict(zip(BANDS, pdc_all))

def save_pdc_to_h5(h5f, subj_id, tissue_labels, contact_labels, pdc_mats):
    subj_id = f"/{subj_id}/"
    with h5py.File(h5f, 'a') as f:
        grp = f.require_group(subj_id)
        dset = grp.require_dataset('tissue_label',len(tissue_labels), DTYPE)
        dset[:] = tissue_labels
        dset.attrs['description'] = "Labels for contacts as SOZ, NIZ, etc, derived from SOZ label in OG matlab struct \
            \nNOTE: str will be saved as ascii. \
            \nFor utf-8 (the normal python format), use .decode() per entry)\n"

        dset = grp.require_dataset('bipole_labels', len(contact_labels), DTYPE)
        dset[:] = contact_labels
        dset.attrs['description'] = "Contact labels indicating where they were implanted \
            \nNOTE: str will be saved as ascii. \
            \nFor utf-8 (the normal python format), use .decode() per entry)\n"

        pdc_group =f.require_group(os.path.join(subj_id,'pdc'))

        for k,v in pdc_mats.items():
            dset = pdc_group.require_dataset(k, v.shape, float)
            dset[:] = v

def z_score_conn(conn_mat:np.ndarray, direction='none', mu=np.array([]), std=np.array([]))->np.ndarray:
    """Z scores connectivity matrices, can handle directed, or non directed
    Zscoring directed connectivity matrices requires specifying the DIRECTION parameter.
    Z scoring along the column ('col') will z-score inward connections. Z scoring the rows
    will score against outward connections. This may be used when you want to normalize against highly connected regions.
    
    NOTE:For example, imagine that the amygdala recieves a lot of connections from every location (a bright column along the amygdala entry)
    if you are interested in the average outward connectivity of a cortical area (the row) these outward statistics will likely be dominated
    by the amygdala outward connectivity. This may obscure the specific connectivity of the cortical area. The idea is to normalize the 
    contributions per area. 

    

    NOTE: Reminder for axes operations
            >>> a = np.array([[1,2],[3,4]])
            >>> print(a)
            [[1 2]
            [3 4]] 
            >>> # axis=0 -> collapse along all rows and returns the average of each column
            >>> np.mean(a, axis=0)
            [2. 3.] 
            >>> #axis=1 -> collapse along columns and return the average of a row
            >>> np.mean(a, axis=1)
            [1.5, 
             3.5]a

    NOTE: reminder for estimation of standard deviation. Matlab uses a degrees of freedom factor to correct for estimating
    variance from a sample of larger distribution. Numpy has ddof=0 by default. We set ddof=1 to correct for this difference
            
    Args:
        conn_mat (np.ndarray): symmetric connectivity matrix
        direction (str, optional): determine direction of z_score. options: 'col', 'row', 'none' Defaults to 'None'.
                                    Use col for generating an average out ward connectivity to score
        mu (np.ndarray, optional): pre_specified mean value to z-score against, often used if comparing to a known baseline, like the interictal period Defaults to None.
        std (np.ndarray, optional): option to pre-specify standard deviation to z-score current connectivity matrix against. Often used when known baseline to compare to such as the intericta
                            period Defaults to None.
    Returns:
        np.ndarray: z_scored connectivity matrix
    """
    n_rows, m_cols = conn_mat.shape
    match direction:
        case 'none':
            mu_mat =mu if len(mu) > 0 else np.nanmean(conn_mat)
            std_mat = std if len(std) > 0 else np.nanstd(conn_mat) 
            return (conn_mat -  mu_mat)/ std_mat
        case 'outward': #columnwise
            inward_mean = mu if len(mu) > 0  else np.nanmean(conn_mat, axis=0).reshape(1,m_cols)
            inward_std = std if len(std) > 0  else np.nanstd(conn_mat, axis=0, ddof=1).reshape(1,m_cols)
            diff = np.subtract(conn_mat, inward_mean)
            return np.divide(diff, inward_std)
        
        case 'inward': #row-wise
            outward_mean = mu if len(mu) > 0 else  np.nanmean(conn_mat, axis =1).reshape(n_rows,1)
            outward_std =  std if len(std) > 0 else np.nanstd(conn_mat, axis =1,ddof=1).reshape(n_rows,1)
            diff = np.subtract(conn_mat, outward_mean)
            return np.divide(diff, outward_std)
        
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def update_dict(net_dict, key, val):
    
    net_dict.setdefault(key, []).append(val)

    return net_dict