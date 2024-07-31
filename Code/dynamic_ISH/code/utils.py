
#Systems I/O
import os
from pathos.pools import ProcessPool
import glob
from scipy.io import loadmat
import mat73
import re
import h5py
import pdb

#data/stats
import numpy as np
import pandas as pd
from scipy.stats import f_oneway, ttest_ind

#specialty

#GLOBAL Variables
BANDS = ['delta', 'theta','alpha', 'beta', 'gamma_low','gamma_high']
DIR = '/mnt/ernie_main/000_Data/SEEG/SEEG_EyesClosed_RestingState/results/Graham_81pats/PDC_RestingState/'
DTYPE = h5py.special_dtype(vlen=str)

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

def paried_region_significance(region_vals:dict)-> np.ndarray:
    """Given a dictionary mapping region name to connectivity values
    returns the pairwise N_reg x N_reg significance matrix

    Args:
        region_vals (dict): dictionary with keys as region names and values as region
        connectivity values

    Returns:
        list[str], np.ndarray: region index names, N_reg x N _reg significance matrix
        making all pairwise comparisons
    """
    return None


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
             3.5]
            
    Args:
        conn_mat (np.ndarray): symmetric connectivity matrix
        direction (str, optional): determine direction of z_score. options: 'col', 'row', 'none' Defaults to 'None'.
        mu (np.ndarray, optional): pre_specified mean value to z-score against, often used if comparing to a known baseline, like the interictal period Defaults to None.
        std (np.ndarray, optional): option to pre-specify standard deviation to z-score current connectivity matrix against. Often used when known baseline to compare to such as the intericta
                            period Defaults to None.
    Returns:
        np.ndarray: z_scored connectivity matrix
    """
    d, _ = conn_mat.shape
    match direction:
        case 'none':
            mu_mat =mu if len(mu) > 0 else np.nanmean(conn_mat)
            std_mat = std if len(std) > 0 else np.nanstd(conn_mat) 
            (conn_mat -  mu_mat)/ std_mat
        case 'col':
            col_mean = mu if len(mu) > 0  else np.nanmean(conn_mat, axis=0).reshape(1,d)
            col_std = std if len(std) > 0  else np.nanstd(conn_mat, axis=0).reshape(1,d)
            diff = np.subtract(conn_mat, col_mean)
            return np.divide(diff, col_std)
        
        case 'row':
            row_mean = mu if len(mu) > 0 else  np.nanmean(conn_mat, axis =1).reshape(d,1)
            row_std =  std if len(std) > 0 else np.nanstd(conn_mat, axis =1).reshape(d,1)
            diff = np.subtract(conn_mat, row_mean)
            return np.divide(diff, row_std)
        
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def update_dict(net_dict, key, val):
    
    net_dict.setdefault(key, []).append(val)

    return net_dict