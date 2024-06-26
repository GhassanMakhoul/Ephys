import os
import re
import logging
import warnings

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


def gen_global_peri_psd(subj_obj, spectral_keys)-> pd.DataFrame:
    """Generate a dataframe with all psd's per recording

    Args:
        subj_obj (_type_): matlab struct loaded from periictal data
        spectral_keys (list, optional): keys for struct to pull out specific periods. Defaults to [].

    Returns:
        pd.DataFrame: dataframe with columns for freq, power, and period
    """

    spectral_dfs = []

    freqs = subj_obj['pdc']['seizure']['pwelch_freqs']
    for key in spectral_keys:
        decomps = subj_obj['pdc']['seizure'][key]
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
    freqs = subj_obj['pdc']['seizure']['pwelch_freqs']
    #NOTE: using the subj_obj dictionary like this is not the most general way to code this pipeline
    #TODO: return to this and remove hardcoded key referencing

    for key in spectral_keys:
        decomps = subj_obj['pdc']['seizure'][key]
        num_ch, _ = decomps.shape
        #These list comps are crazy -> repeats the label for number of frequency decomps
        ch_labels = [l for label in contact_labels for l in [label]*num_ch ] 
        #repeats frequency so that there are x values per contact psd
        freq_rep = [f for fs in [freqs for i in range(num_ch)] for f in fs]
        flat_decomp = np.reshape(-1,1)
        df = pd.DataFrame()
        df['freq'] = freq_rep
        df['power'] = flat_decomp
        df['labels'] = ch_labels
        df['period'] = key
        spectral_dfs.append(df)
    return pd.concat(spectral_dfs)

    
def gen_patient_psd():
    """should generate all peri psd's per patient's seizure
    
    """


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
    power_range = power[f_inds[:,0],:]
    freq_range = freqs[f_inds[:,0],:]
    lr = LinearRegression()
    lr.fit(freq_range, power_range)
    slope = lr.coef_
    return slope





    