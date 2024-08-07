import sys
import os
import re
import logging
from loguru import logger
import warnings
import getopt
from tqdm import tqdm

import pdb

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
    window_dict = prep_window_dict(subj_obj)
    pwelch_all_windows = read_conn_struct(subj_obj,'pdc','pwelch_all_windowed')
    for key, window_designation in window_dict.items():
        decomps = pwelch_all_windows[key,:,:]
        avg_decomp = np.mean(decomps, axis=0)
        df = pd.DataFrame()
        df['freq'] = freqs
        df['power'] = np.log(avg_decomp)
        df['window_designations'] = window_designation
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
    freqs = read_conn_struct(subj_obj, 'pdc','pwelch_freqs')
  
    window_dict = prep_window_dict(subj_obj)
    pwelch_all_windows = read_conn_struct(subj_obj,'pdc','pwelch_all_windowed')
    #NOTE: using the subj_obj dictionary like this is not the most general way to code this pipeline
    #TODO: return to this and remove hardcoded key referencing

    for key, window_designation in window_dict.items():
        pwelch = pwelch_all_windows[key, :, :]
        num_ch, num_freq = pwelch.shape
        #These list comps are crazy -> repeats the label for number of frequency decomps
        ch_labels = [l for label in contact_labels for l in [label]*num_freq ] 
        #repeats frequency so that there are x values per contact psd
        freq_rep = [f for fs in [freqs for i in range(num_ch)] for f in fs]
        flat_decomp = pwelch.reshape(-1,1)
        df = pd.DataFrame()
        df['freq'] = freq_rep
        # pdb.set_trace()
        df['power'] = np.log(flat_decomp)
        df['labels'] = ch_labels
        df['period'] = key
        df['window_deisgnations'] = window_designation
        spectral_dfs.append(df)
    return pd.concat(spectral_dfs)


def score_period_power(power_decomps: np.ndarray, freq_inds, window_dict, agg_win = "interictal", buffer=60, win_size=5, stride=1,stats='full', **kwargs)->dict:
    """Generate reference statistics to compute z-score PSD's against

    Args:
        power_decomps (np.ndarray): raw pwelch should be T X N_ch X M_freqs
        freq_inds (np.ndarray): array that is M_freqs long and should track the frequency bins of p_welch
        window_dict (_type_): dictionary mapping indices of windows to window deisgnations that 
                            tracks the beginning middle and end of a recording 0_0_0 = interictal_interictal_interictal
        agg_win (str, optional): windows to create references against. Defaults to "interictal".
        buffer (int, optional): exclusion time before window transitions  . Defaults to 60.
        win_size (int, optional): length in s of window. Defaults to 5.
        stride (int, optional): stride of window in s. Defaults to 1.
        stats (str, optional): will probable remove. Defaults to 'full'.

    Returns:
        dict: _description_
    """

    power_decomps = filter_periods(power_decomps, window_dict, agg_win, buffer, win_size)
    power_decomps = np.log(np.array(power_decomps))
    mu_decomps = np.mean(power_decomps, axis=0)
    mu_auc = [get_power(mu_decomps, freq, freq_inds) for freq in BANDS ]
    mu_dict = dict(zip(BANDS, mu_auc))
    #std dev is going to be a bit trickier
    std_dict = {}
    widths = np.diff(freq_inds)
    widths = np.append(widths, widths[-1]) #make up for losing one entry
    for freq, freq_ranges in BAND_RANGES.item():
        lo, hi = freq_ranges
        inds_h = np.where(freq_inds <= hi)
        inds_l = np.where(freq_inds > lo)
        inds  = np.intersect1d(inds_h, inds_l)
        
        # AUC's (N_periods x CH x N_freq_inds) * Widths (N_freq_inds,1)
        # AUC diff (N_periods x CH x 1) - CH x 1 #NOTE: should I make a repmat?
        # Check shapes
        auc_diff =  np.dot(power_decomps[:,:,inds], widths[inds]) - mu_dict[freq]
        auc_var = np.dot(auc_diff, auc_diff)
        auc_std = np.sqrt(np.mean(auc_var, axis=0))
        std_dict[freq] = auc_std
    ref_stats = dict(zip(BANDS, [(mu_dict[b], std_dict[b]) for b in BANDS]))
    return ref_stats

def assemble_psd_verbose(subj_obj, sz_band='beta',window='full'):
    
    freqs = read_conn_struct(subj_obj, 'pdc','pwelch_freqs')
    if window == 'full':
        window_dict = prep_window_dict(subj_obj)
        pwelch_all_windows = read_conn_struct(subj_obj,'pdc','pwelch_all_windowed')

    psd_dfs = []
    bip_names = get_chan_names(subj_obj)
    soz_per_szr = read_conn_struct(subj_obj, 'pdc', 'soz_per_seizure')
    contact_label = format_soz(soz_per_szr)
    freqs = read_conn_struct(subj_obj, 'pdc', 'pwelch_freqs')
    subj = read_conn_struct(subj_obj,'pdc','patID')
    sz_type = read_conn_struct(subj_obj,'pdc', 'sz_type')
    eventID = read_conn_struct(subj_obj, 'pdc','eventID')
    ref_stats = score_period_power(pwelch_all_windows, freqs, window_dict)
    for key, window_designation in window_dict.items():
        pwelch = pwelch_all_windows[key, :,:]
        #NOTE: TAKE LOG of psd before calculating power AUC
        pwelch = np.log(pwelch)
        band_pow = get_power(pwelch, sz_band, freqs)
        #TODO: add iteration over multiple bands
        df = pd.DataFrame(data=band_pow,columns=[f"power_{sz_band}"])
        z_bands = z_score(df[f"power_{sz_band}"], ref_stats[sz_band])
        df.insert(loc=1, column=f"z_{sz_band}", value=z_bands)
        df.insert(loc=1, column='region', value=contact_label)
        df.insert(loc=1, column='bip', value=bip_names)
        df.insert(loc=1, column='period', value=key)
        df.insert(loc=1, column='window_designations', value=window_designation)
        #Z-score

        #NOTE: is there a better way to z-score? Maybe I repmat and just use over aggregated psd_df

        psd_dfs.append(df)
    psd_dfs = pd.concat(psd_dfs)
    psd_dfs.insert(loc=1, column='patID', value = subj)
    psd_dfs.insert(loc=1, column='sz_type', value=sz_type)
    psd_dfs.insert(loc=1, column='eventID',value= eventID)
    return psd_dfs


def z_score(samps, mu, std):
    return (samps - mu)/std
        
def get_power(psd, freq, freq_inds):
    """Return the area under the curve of the PSD in the frequency range of interes
    NOTE: THIS ASSUMES that you are passing the logarithm of the raw pwelch.
    Args:
        psd (np.ndarray): 2D array of PSD's shape should be N_ch x n_freq
        freq (str): frequency range of interest, 
        freq_inds (np.npdarray) : array of frequency bins, each entry should correspond to a
                    frequency bin measured by PSD
    Returns
        np.ndarray of AUC for frequency of interest for N_ch's
    """
    #psd = np.log(psd) #assumes that the log has already been calculated
    lo, hi = BAND_RANGES[freq]
    inds_h = np.where(freq_inds <= hi)
    inds_l = np.where(freq_inds > lo)
    inds  = np.intersect1d(inds_h, inds_l)
    widths = np.diff(freq_inds) #NOTE: that this may be off by 1!

    return np.dot(psd[:,inds], widths[inds])


def get_reg_ei(subj_obj, window='full')->pd.DataFrame:
    """Loops through peiods and pulls out the psd per period then calculates
    the E/I index of the period using the get_ei_slope method. Builds a dataframe of each 
    periods values on a per contact level.

    Args:
        subj_obj (_type_): matlab struct 
        window_type (_type_): _keys for struct to pull out specific periods' psd
        contact_label (_type_): labels of the channels to be used for boolean indexing 

    Returns:
        pd.DataFrame: dataframe with columns for period, E/I slope, region
    """
    
    freqs = read_conn_struct(subj_obj, 'pdc','pwelch_freqs')
    if window == 'full':
        window_dict = prep_window_dict(subj_obj)
        pwelch_all_windows = read_conn_struct(subj_obj,'pdc','pwelch_all_windowed')

    ei_dfs = []
    soz_per_szr = read_conn_struct(subj_obj, 'pdc', 'soz_per_seizure')
    contact_label = format_soz(soz_per_szr)
    freqs = read_conn_struct(subj_obj, 'pdc', 'pwelch_freqs')
    subj = read_conn_struct(subj_obj,'pdc','patID')
    sz_type = read_conn_struct(subj_obj,'pdc', 'sz_type')
    eventID = read_conn_struct(subj_obj, 'pdc','eventID')
    for key, window_designation in window_dict.items():
        pwelch = pwelch_all_windows[key, :,:]
        for reg in set(contact_label):
            reg_inds = np.where(np.array(contact_label) == reg)
            # import pdb
            # pdb.set_trace()
            reg_psds = pwelch[reg_inds]
            ei_vals = [get_ei_slope(freqs, reg_psds[i,:]) for i in range(reg_psds.shape[0])]
            df = pd.DataFrame()
            df['e_i'] = ei_vals
            df['region'] = reg
            df['period'] = key
            df['window_designations'] = window_designation
            ei_dfs.append(df)
    ei_dfs = pd.concat(ei_dfs)
    ei_dfs ['patID'] = subj
    ei_dfs['sz_type'] = sz_type
    ei_dfs['eventID'] = eventID
    return ei_dfs
        




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
    

def get_reg_ei_para(sub_objects:list, cores =12, spectral_keys='all')->pd.DataFrame:
    """Similar to get_reg_ei just parallellized

    Args:
        sub_objects (list): list of struct objects to apply get_reg_ei to
        spectral_keys (str, optional): _description_. Defaults to 'all'.

    Returns:
        pd.DataFrame: dataframe with columns for period, E/I slope, region
    """
    p = Pool(cores)
    dfs =  p.map(get_reg_ei, sub_objects)
    return pd.concat(dfs)


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
    log_power = np.log(power)
    power_range = log_power[f_inds]
    freq_range = freqs[f_inds]
    lr = LinearRegression()
    # import pdb
    # pdb.set_trace()
    lr.fit(freq_range.reshape(-1,1), power_range.reshape(-1,1))
    slope = lr.coef_
    return slope[0][0]

def gen_ei_dfs(data_dir, out_dir,num_cores=20):
    """generates the ei_df in a paralellized fashion that saves output. Saves intermediate csvs too
    and returns final df with pd.concat
    NOTE: THIS assumes that the directory you're passing in as DATA_DIR is the top level directory
    for ONE AND ONLY ONE subject/patient at a time. The code will run the same with a differnt directory
    organization tbh, but the book keeping may be off somewhere downstream/things may be slower too.

    Args:
        data dir: data with  structs /patID/struct.mat
        out_put : where to save final csv and intermediate data
    """
    from tqdm import tqdm

    sub_paths= glob.glob(os.path.join(data_dir, "*PDC.mat"))
    
    ei_dfs = []
    count = 1
    assert len(sub_paths) > 0, f"No files to load in {sub_paths}, check {data_dir}"
    for f_paths in tqdm(chunker(sub_paths, num_cores)):
        structs = load_structs(f_paths,num_cores)
        incl_inds = [i for i in range(len(structs)) if structs[i] != None]
        structs = [structs[i] for i in incl_inds]
        res_dfs = get_reg_ei_para(structs, cores=num_cores)
        subj = res_dfs.patID.values[0]
        res_dfs.to_csv(os.path.join(out_dir,f'ei_{subj}_tmp_{count}.csv'),index=False)
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
    with open(config_f, 'r') as f:
        config =  yaml.safe_load(f)
        config = config['ei_bal']
    logger.info(f"Running eibal on {datadir} with the following config:\n{config}")
    ei_df = gen_ei_dfs(datadir, pathout)
    subj = ei_df.patID.values[0]
    ei_df.to_csv(os.path.join(pathout, f"ei_bal_{subj}.csv"),index=False)

if __name__ == '__main__':
    main(sys.argv[1:])    