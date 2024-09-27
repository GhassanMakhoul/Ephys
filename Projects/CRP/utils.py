
#Systems I/O
import os
from pathos.pools import ProcessPool
import glob
import mat73
from scipy.io import loadmat as lm
import re
import h5py

#data/stats
import numpy as np
import pandas as pd

#specialty

#GLOBAL Variables
BANDS = ['delta', 'theta','alpha', 'beta', 'gamma_low','gamma_high']
DIR = '/mnt/ernie_main/000_Data/SEEG/SEEG_EyesClosed_RestingState/results/Graham_81pats/PDC_RestingState/'
DTYPE = h5py.special_dtype(vlen=str)

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
    p = "[0-9][a-zA-Z]"
    bip_match = re.search(p, bipole)

    assert bip_match != None, f"Bipole {bipole}, is not captures by regex pattern {p}"
    s, _ = bip_match.span()
    return bipole[0:s+1] + "-" + bipole[s+1:]

def format_bipoles(char_list):
    return [format_bipole(c) for c in char_list]

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
