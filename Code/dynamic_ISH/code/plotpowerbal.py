import logging

logging.getLogger('mat73').setLevel(logging.CRITICAL)
import os
import re
from scipy.io import loadmat
import mat73
logging.getLogger('mat73').setLevel(logging.CRITICAL)

import warnings

from collections import Counter
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from sklearn import linear_model

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from utils import *
from connectivity_dynamics import *
from eibal import *


def get_involved_inds(event_df, band="beta", threshold=2):
    """Given an event_df, which contains the dynamics of all bipoles for one seizure event,
    thresholds all areas which achieve an average band power above threshold 
    and returns a column detailing which bipoles are involved as BANDS_increased"""
    ictal_df = event_df[event_df.win_label.isin(['early_ictal','ictal','late_ictal'])]
    ictal_bip_df = ictal_df[['bip', f'z_{band}',]].groupby('bip').mean().reset_index()
    # get regions with > threshold power
    thresh_bool = ictal_bip_df[f'z_{band}'].values > threshold
    power_dict = dict(zip(ictal_bip_df.bip, thresh_bool))
    all_bip_thresh = event_df.bip.apply(lambda x: power_dict[x])
    event_df.insert(loc=len(event_df.columns), column=f"{band}_involved", value=all_bip_thresh)
    return event_df


def merge_flow__power(flow_df: pd.DataFrame, power_df: pd.DataFrame, band='beta'):
    """merge bipole level peri_ictal flow_df to beta power channel
    assumes that eventID is same types"""
    power_flow_dfs = []
    for event in set(power_df.eventID.values):
        pow_event_df = center_onset(power_df[power_df.eventID == event])
        pow_event_df = get_involved_inds(pow_event_df)
        flow_event_df = center_onset(flow_df[flow_df.eventID == event])
   
        df = pow_event_df[['bip','z_beta','win_sz_centered', f'{band}_involved']].merge(\
            flow_event_df,how='right',
            right_on=['win_sz_centered','src_bip'],
            left_on=['win_sz_centered','bip'])
        df['region_involved'] = df.apply( \
            lambda x: f"{x['source']}_{x['target']}_{x['beta_involved']}",\
            axis=1)
        power_flow_dfs.append(df)
    return pd.concat(power_flow_dfs)

def load_dfs(flow_fname, power_fname):
    """Loads the channel level connectivity and power dataframes"""
    power_df = pd.read_csv(power_fname)
    power_df = center_onset(power_df)
    power_df['eventID'] = power_df.eventID.apply(int)

    flow_df = pd.read_csv(flow_fname)
    flow_df = center_onset(flow_df)
    flow_df['eventID'] = flow_df.eventID.apply(int)
    assert len(set(flow_df.eventID)) == len(set(power_df.eventID)), "need same number of events!"
    return flow_df, power_df

def plot_flow_power(flow_power_df, fname):
     with sns.plotting_context("paper"):
            grid = sns.FacetGrid(flow_power_df, row='source',row_order=['nz','soz'],
                                col='win_label', 
                                col_order=['interictal', 'pre_ictal','early_ictal','ictal','late_ictal','early_post_ictal','post_ictal'],
                                ) 
            ax = grid.map_dataframe(sns.scatterplot, y='value',x='z_beta', hue='ictal_involvement')
            grid.add_legend()
            grid.figure.suptitle("Pre-Ictal PDC",y=1.01)
            plt.savefig(f"../viz/{fname}.pdf", transparent=True)


def main(argv):
    DATA_DIR = '/mnt/ernie_main/Ghassan/ephys/data/ei_bal'
    sub_paths= glob.glob(os.path.join(DATA_DIR, "*power*"))

    path = sub_paths[1]
    print(f"Loading sub_path {path}")
    tst_df = pd.read_csv(path)
    tst_df = center_onset(tst_df)



if __name__ == '__main__':
     main(sys.argv[1:])