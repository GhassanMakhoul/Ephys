import sys
import os
import re
import logging
import warnings
import getopt
from tqdm import tqdm


import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from utils import *
from connectivity_dynamics import *

SPECTRAL_KEYS = [ 'pwelch_ictal',
                'pwelch_interictal',
                'pwelch_post',
                'pwelch_pre']

def gen_global_peri_psd(subj_obj, spectral_keys)-> pd.DataFrame:
    """Generate a dataframe with all psd's per recording

    Args:
        subj_obj (_type_): matlab struct loaded from periictal data
        spectral_keys (list, optional): keys for struct to pull out specific periods. Defaults to [].

    Returns:
        pd.DataFrame: dataframe with columns for freq, power, and period
    """

    spectral_dfs = []

    freqs = subj_obj['pwelch_freqs']
    for key in spectral_keys:
        decomps = subj_obj[key]
        avg_decomp = np.mean(decomps, axis=0)
        df = pd.DataFrame()
        df['freq'] = freqs
        df['power'] = np.log(avg_decomp)
        df['period'] = key
        spectral_dfs.append(df)
    return pd.concat(spectral_dfs)

def gen_contact_peri_psd(subj_obj, spectral_keys, contact_labels)->pd.DataFrame:
    """Returns a datafrme with the contact level power spectral decompositions for each 
    period in Spectral keys

    Args:
        subj_obj (_type_): matlab strcut loaded from peri-ictal data
        spectral_keys (_type_): keys for struct to pull out specific periods
        contact_labels (_type_): designation of channels as 'SOZ', 'PZ', or 'NIZ', 

    Returns:
        pd.DataFrame: dataframe with a contact column. NOTE: this DF will be built up over many 
        contact level dfs. Meaning, this will have many many rows for all the measurements
    """
    spectral_dfs = []
    freqs = subj_obj['pwelch_freqs']
    #NOTE: using the subj_obj dictionary like this is not the most general way to code this pipeline
    #TODO: return to this and remove hardcoded key referencing

    for key in spectral_keys:
        decomps = subj_obj[key]
        num_ch, num_freq = decomps.shape
        #These list comps are crazy -> repeats the label for number of frequency decomps
        ch_labels = [l for label in contact_labels for l in [label]*num_freq ] 
        #repeats frequency so that there are x values per contact psd
        freq_rep = [f for fs in [freqs for i in range(num_ch)] for f in fs]
        flat_decomp = decomps.reshape(-1,1)
        df = pd.DataFrame()
        df['freq'] = freq_rep
        # pdb.set_trace()
        df['power'] = np.log(flat_decomp)
        df['labels'] = ch_labels
        df['period'] = key
        spectral_dfs.append(df)
    return pd.concat(spectral_dfs)



def get_reg_ei(subj_obj, spectral_keys='all')->pd.DataFrame:
    """Loops through peiods and pulls out the psd per period then calculates
    the E/I index of the period using the get_ei_slope method. Builds a dataframe of each 
    periods values on a per contact level.

    Args:
        subj_obj (_type_): matlab struct 
        spectral_keys (_type_): _keys for struct to pull out specific periods' psd
        contact_label (_type_): labels of the channels to be used for boolean indexing 

    Returns:
        pd.DataFrame: dataframe with columns for period, E/I slope, region
    """
    if 'pdc' in subj_obj.keys():
        subj_obj = subj_obj['pdc']['seizure']
    if spectral_keys == 'all':
        spectral_keys = SPECTRAL_KEYS

    freqs = subj_obj['pwelch_freqs']
    ei_dfs = []
    contact_label = format_soz(subj_obj['soz_per_seizure'])

    for key in spectral_keys:
        decomps = subj_obj[key]
        for reg in set(contact_label):
            reg_inds = np.where(np.array(contact_label) == reg)
            # import pdb
            # pdb.set_trace()
            reg_psds = decomps[reg_inds]
            ei_vals = [get_ei_slope(freqs, reg_psds[i,:]) for i in range(reg_psds.shape[0])]
            df = pd.DataFrame()
            df['e_i'] = ei_vals
            df['region'] = reg
            df['period'] = key
            ei_dfs.append(df)
    return pd.concat(ei_dfs)
        


def get_reg_ei_para(sub_objects:list, sub_list, sz_list, cores =12, spectral_keys='all')->pd.DataFrame:
    """Similar to get_reg_ei just parallellized

    Args:
        sub_objects (list): list of struct objects to apply get_reg_ei to
        spectral_keys (str, optional): _description_. Defaults to 'all'.

    Returns:
        pd.DataFrame: dataframe with columns for period, E/I slope, region
    """
    p = Pool(cores)
    dfs =  p.map(get_reg_ei, sub_objects)
    for i, df in enumerate(dfs):
        df['subj'] = sub_list[i]
        df['sz_type'] = sz_list[i]
    return pd.concat(dfs)

def gen_patient_psd(subj_id, pat_structs,spectral_keys, level='gen'):
    """should generate all peri psd's per patient's seizure

    """
    
    psd_dfs = []
    for subj_obj in pat_structs:
        subj_obj = subj_obj['pdc']['seizure']
        # assert subj_id == subj_obj['patID'],f"Subject {subj_id} \
        #     doesnt match with struct f{subj_obj['patID']}!"
        
        if level == 'gen':
            df = gen_global_peri_psd(subj_obj, spectral_keys)
        elif level =='channel':
            soz_labels = format_soz(subj_obj['soz_per_seizure'])
            df = gen_contact_peri_psd(subj_obj, spectral_keys, soz_labels)
        else:
            return #TODO: raise error?
        # df['seizure'] = subj_obj['sz_type']
        psd_dfs.append(df)
    return pd.concat(psd_dfs)


    #NOTE: using soz_per_seizure
    map_label



def get_ei_slope(freqs:np.ndarray, power:np.ndarray, ei_range=(30,50))->float:
    """given a spectral decomposistion (POWER), determines the slope of a frequency range of interest
    indexed by FREQ and determined by EI_RANGE. Fits a linear model using using sklearn's linear regression
    NOTE: consider returning confidence interval, p_value of fit, etc
    TODO: revisit algorithm in the future

    Args:
        freqs (np.ndarray): list of frequency bins estimated by psd, make sure array is 2D!
        power (np.ndarray): list of power values should track freqs, make sure array is 2D!
        ei_range (tuple, optional): excitation inhibition range. Defaults to (30,50).

    Returns:
        float: value of slope of PSD at frequency range of interest
    """

    assert len(freqs) == len(power), "Freqs array needs to track power array! check array sizes."
    f_inds = np.logical_and(freqs > ei_range[0], freqs < ei_range[1])
 
    power_range = power[f_inds]
    freq_range = freqs[f_inds]
    lr = LinearRegression()
    # import pdb
    # pdb.set_trace()
    lr.fit(freq_range.reshape(-1,1), power_range.reshape(-1,1))
    slope = lr.coef_
    return slope[0][0]

def gen_ei_dfs(data_dir, out_dir):
    """generates the ei_df in a paralellized fashion that saves output. Saves intermediate csvs too
    and returns final df with pd.concat

    Args:
        data dir: data with  structs /patID/struct.mat
        out_put : where to save final csv and intermediate data
    """
    from tqdm import tqdm

    sub_paths= glob.glob(os.path.join(data_dir, "*pat*", "*.mat"))
    num_cores = 20
    ei_dfs = []
    count = 1
    for f_paths in tqdm(chunker(sub_paths, num_cores)):
        structs = load_structs(f_paths,num_cores)
        incl_inds = [i for i in range(len(structs)) if structs[i] != None]
        structs = [structs[i] for i in incl_inds]

        sub_list = [sub_path.split("/")[-2] for sub_path in f_paths]
        sub_list = [sub_list[i] for i in incl_inds] #NOTE: god this is messy. TODO: fix struct chars for ID and sztype
        sz_list = [sub_path.split("/")[-1].split("_")[-2] for sub_path in f_paths]
        sz_list = [sz_list[i] for i in incl_inds]

        res_dfs = get_reg_ei_para(structs, sub_list, sz_list, cores=num_cores)
        res_dfs.to_csv(os.path.join(out_dir,f'ei_tmp_{count}.csv'),index=False)
        print(f"Saved first {count*num_cores} seizures to folder: {out_dir}")
        ei_dfs.append(res_dfs)
        
        count += 1

    return pd.concat(ei_dfs)


def main(argv):
    opts, _ = getopt.getopt(argv,"d:p:c:",["datadir=",'pathout=','config='])
    for opt, arg in opts:
        if opt in ("-d", 'datadir'):
            datadir = arg
        elif opt in ('-p', '--pathout'):
            pathout = arg
        elif opt in ("-c", '--config'):
            config_f = arg
    #TODO use yamls and configs
    # with open(config_f, 'r') as f:
    #     config =  yaml.safe_load(f)
    ei_df = gen_ei_dfs(datadir, pathout)
    ei_df.to_csv(os.path.join(pathout, "ei_bal.csv"),index=False)

if __name__ == '__main__':
    main(sys.argv[1:])    