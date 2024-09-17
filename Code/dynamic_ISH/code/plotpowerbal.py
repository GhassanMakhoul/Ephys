import logging

logging.getLogger('mat73').setLevel(logging.CRITICAL)
import os
import sys
import getopt
import re
from scipy.io import loadmat
import yaml
logging.getLogger('mat73').setLevel(logging.CRITICAL)
from multiprocessing import Pool
from functools import partial

import warnings

from collections import Counter
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from sklearn import linear_model

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from utils import *
from connectivity_dynamics import center_onset


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


def merge_flow_power(flow_df: pd.DataFrame, power_df: pd.DataFrame, subsample=1, band='beta',**kwargs):
    """merge bipole level peri_ictal flow_df to beta power channel
    assumes that eventID is same types and that there is an equal number
    of events between the flow_df and the power_df"""
    power_flow_dfs = []
    for event in set(power_df.eventID.values):
        pow_event_df = power_df[power_df.eventID == event]
        pow_event_df = get_involved_inds(pow_event_df)
        flow_event_df = flow_df[flow_df.eventID == event]
        assert "win_sz_centered" in flow_event_df.columns, "Need to center flow first!"
        assert "win_sz_centered" in pow_event_df.columns, "Need to center flow first!"
   
        df = pow_event_df[['bip','z_beta','win_sz_centered', f'{band}_involved']].merge(\
            flow_event_df,how='right',
            right_on=['win_sz_centered','src_bip'],
            left_on=['win_sz_centered','bip'])
        df['region_involved'] = df.apply( \
            lambda x: f"{x['source']}_{x['target']}_{x['beta_involved']}",\
            axis=1)
        if subsample > 1:
            df = subsample_df(df, subsample, groups=['source','win_label','region_involved'])
        power_flow_dfs.append(df)
    return pd.concat(power_flow_dfs)

def subsample_df(df:pd.DataFrame, factor:int, random_state=42, groups=['win_label'])->pd.DataFrame:
    """Returns a subsampled dataframe that selects samples along groupings 
    specified by GROUPS arg.

    Args:
        df (pd.DataFrame): dataframe to resample, often a peri-ictal verbose df
        factor (int): amount to refactor sample by 2 -> return half the sample, 5-> 1/5 samples returned
        groups (list, optional): stratifications to segment df by,. Defaults to ['win_label'].

    Returns:
        pd.DataFrame: returns subsampled df
    """
    n_rows = len(df)
    logger.info(f"Original data frame has {n_rows} row")
    group_df = df.groupby(by=groups)
    count_df= group_df.count()
    min_grp_size = min(count_df.min().values)
    n_resamp = max(n_rows//factor, min_grp_size)
    resamp_df = group_df.sample(n_resamp, random_state=random_state,replace=True)
    resamp_df = resamp_df.drop_duplicates()
    logger.info(f"New dataframe has {len(resamp_df)} rows with minimum {max(min_grp_size//factor,min_grp_size)} samples per group.\
            \n\t\t\t\t\tSubsample factor : {n_rows/len(resamp_df)} \n\t\t\t\t\tOriginal supsample factor: {factor}")
    return resamp_df


def load_dfs(flow_fname, power_fname):
    """Loads the channel level connectivity and power dataframes"""
    power_df = pd.read_csv(power_fname)
    power_df['eventID'] = power_df.eventID.apply(str)

    flow_df = pd.read_csv(flow_fname)
    flow_df['eventID'] = flow_df.eventID.apply(str)
    assert len(set(flow_df.eventID)) == len(set(power_df.eventID)), "need same number of events!"
    return flow_df, power_df

def kde_flow_power(flow_power_df, pltname, x='value',y='z_beta', **kwargs):
     grp_flowpow_df = flow_power_df.groupby(['win_label','source'])
     with sns.plotting_context("paper"):
            fig, axs = plt.subplots(2,7, sharex=True, sharey=True):
            for source, i in enumerate(['nz', 'soz']):
                for period, t, in enumerate(['interictal', 'pre_ictal','early_ictal','ictal','late_ictal','early_post_ictal','post_ictal']):
                    ax = axs[i,t]
                    plot_df = grp_flowpow_df.get_group((period, source))
                    sns.kdeplot(data=plot_df, x=x, y=y, hue='region_involved', legend=True, ax=ax, palette=FLOWMAP)

            plt.savefig(f"../viz/{pltname}.pdf", transparent=True, bbox_inches="tight")

def scatter_flow_power(flow_power_df, pltname, x='value', y='z_beta', **kwargs):
     with sns.plotting_context("paper"):
            grid = sns.FacetGrid(flow_power_df, row='source',row_order=['nz','soz'],
                                col='win_label', 
                                col_order=['interictal', 'pre_ictal','early_ictal','ictal','late_ictal','early_post_ictal','post_ictal'],
                                legend_out=True
                                ) 
            ax = grid.map_dataframe(sns.scatterplot, y=y,x=x, hue='region_involved')
            grid.add_legend()
            grid.figure.suptitle(f"Peri-ictal {x} vs {y} Power",y=1.0)
            plt.savefig(f"../viz/{pltname}.pdf", transparent=True, bbox_inches="tight")


def load_merge_subject(subject, power_dir="/mnt/ernie_main/Ghassan/ephys/data/ei_bal/", flow_dir="/mnt/ernie_main/Price/ephys/data/periconnectivity", **kwargs):
    """Given a subject, coordinates dataframe loading
    and merging, then generates plots for power vs connectivity

    Args:
        subjects (str): str of subject id  to load 
        power_dir (str): directory leading to channel wise, peri-ictal power csvs
        flow_dir (str): directory leading to channel-wise peri-ictal connecivity csvs
    """
    try:
        pow_fname = os.path.join(power_dir, f"power_bal_{subject}_centered.csv")
        assert os.path.exists(pow_fname)
        flow_fname = os.path.join(flow_dir, f"peri_ictal_flow_verbose_{subject}_centered.csv")
        assert os.path.exists(flow_fname)
    except AssertionError:
        logger.warning(f"Subject {subject} not found in flow and power!")
        return pd.DataFrame()
    #maybe assert that there is a file in each glob before indexing?
    flow_df, power_df = load_dfs(flow_fname, pow_fname)
    
    return merge_flow_power(flow_df, power_df, **kwargs)

def center_seizure_onset(subjects, power_dir="/mnt/ernie_main/Ghassan/ephys/data/ei_bal/", flow_dir="/mnt/ernie_main/Price/ephys/data/periconnectivity", **kwargs):
    pre_centered = 0
    for subject in subjects:
        centered_pow_dfs = []
        centered_flow_dfs = []
        # for event in events ids
        pow_fname = os.path.join(power_dir, f"power_bal_{subject}.csv")
        flow_fname = os.path.join(flow_dir, f"peri_ictal_flow_verbose_{subject}.csv")

        pow_c_fname = pow_fname.replace(".csv", "_centered.csv")
        flow_c_fname = flow_fname.replace('.csv', '_centered.csv')
        if os.path.exists(pow_c_fname) and os.path.exists(flow_c_fname):
            pre_centered +=1
            continue
        #maybe assert that there is a file in each glob before indexing?
        flow_df, power_df = load_dfs(flow_fname, pow_fname)
        for event in set(power_df.eventID.values):
            pow_event_df = power_df[power_df.eventID == event]
            flow_event_df =flow_df[flow_df.eventID == event]
            centered_ev_flow = center_onset(flow_event_df)
            centered_flow_dfs.append(centered_ev_flow)
        
            centered_ev_power = center_onset(pow_event_df)
            centered_pow_dfs.append(centered_ev_power)
        
        centered_flow_dfs = pd.concat(centered_flow_dfs)
        centered_pow_dfs = pd.concat(centered_pow_dfs)
        
        centered_flow_dfs.to_csv(flow_c_fname)
        centered_pow_dfs.to_csv(pow_c_fname)
    logger.info(f"Centered: {len(subjects)- pre_centered} subjects' flow and power csvs to seizure onset")
    
def verify_subjlist(subjlist, flow_dir, pow_dir):
    """returns a list of subjects that has both verbose flow and verbose power csvs"""
    verify_list = []
    for sub in subjlist:
        flow_f = os.path.join(flow_dir, f"peri_ictal_flow_verbose_{sub}.csv")
        pow_f = os.path.join(pow_dir,  f"power_bal_{sub}.csv")
        if os.path.exists(flow_f) and os.path.exists(pow_f):
            verify_list.append(sub)
    diff = set(subjlist).difference(verify_list)
    if len(diff) > 0:
        logger.warning(f"The following patients were not aligned across power and flow csvs: {diff}")
    return verify_list


def load_merge_all(subjects, num_cores=20, merged_fname="merged_df", **kwargs):
    #may need to have spec test for certain centering 
    # if len(glob.glob(os.path.join(kwargs['power_dir'], "*_centered.csv"))) < len(subjects):
    #     logger.info("Centering Seizure events first!")
    #     center_seizure_onset(subjects, **kwargs)
    logger.info(f"Loading subjects {len(subjects)}")
    if len(subjects) <2 :
        #for debugging purposes
        merged_dfs =[ load_merge_subject(subj, **kwargs) for subj in subjects]
    elif len(subjects) >1:
        num_cores = min(len(subjects), num_cores)
        p = Pool(num_cores)
        merged_dfs = p.map(partial(load_merge_subject, **kwargs), subjects)
        p.close()
        p.join()
        merged_dfs = [m_df for m_df in merged_dfs if not m_df.empty]
    merged_dfs = pd.concat(merged_dfs)
    logger.success(f"Merged and loaded {len(subjects)}")
    merged_dfs.to_csv(merged_fname)
    logger.info(f"Saved merged df to {merged_fname}")

def plot_subjects(merged_dfs, pltname='flow_pow.pdf', band='beta', plot_type='kde', **kwargs):   
    #TODO maybe make x more modular
    x = f"z_{band}"
    logger.info(f"About to plot peri-ictal {plot_type}")
    match plot_type:
        case "kde":
            kde_flow_power(merged_dfs, pltname,x=x, **kwargs)
        case "scatter":
            scatter_flow_power(merged_dfs, pltname, x=x, **kwargs)

    logger.success(f"plotting successful! Saved plots to {pltname}")


def main(argv):
    logger.info(argv)
    kwargs = {}
    opts, _ = getopt.getopt(argv,"s:f:p:t:c:l:",["subjlist=","flowdir=",'powdir=',"plotname=",'config=','logdir='])
    for opt, arg in opts:
        if opt in ("-s", "subjlist="):
            subjlist = arg
        elif opt in ('-f', '--flowdir='):
            flowdir = arg
        elif opt in ("-p", '--powdir='):
            powdir = arg
        elif opt in ("-t","--plotname="):
            pltname = arg
        elif opt in ("-c", "--config="):
            config_f = arg
            with open(config_f, 'r') as f:
                    config =  yaml.safe_load(f)
            kwargs = config['plt_flow_power']
        elif opt in ("-l", "--logdir="):
            logdir = arg
    logger.add(os.path.join(logdir, "pow_flow_plot_run.log"), enqueue=True,level=40)
    kwargs['pltname'] = pltname
    kwargs['power_dir'] = powdir
    kwargs['flow_dir'] = flowdir
    pipeline = kwargs['pipeline'] if 'pipeline' in kwargs.keys() else 'plot_center'
    logger.info(f"Running pipeline with following kwargs:\n\t\t\t\t{kwargs}")
    subjlist = pd.read_csv(subjlist,index_col=False, header=None)
    subjlist = [s[0] for s in subjlist.values] ## parse out 
    logger.info(f"Runinng pipeline on {len(subjlist)} subjects")
    subjlist_verified = verify_subjlist(subjlist, flowdir, powdir)
    
    subjlist_verified = subjlist[0:5]
    match pipeline:
        case "full_plot":
            center_seizure_onset(subjlist_verified, **kwargs)
            merged_dfs = load_merge_all(subjlist_verified, **kwargs)
            plot_subjects(merged_dfs, **kwargs)
        case "merge_plot":
            merged_dfs = load_merge_all(subjlist_verified, **kwargs)
            plot_subjects(merged_dfs, **kwargs)
        case "center":
            logger.info("Centering Subjects")
            center_seizure_onset(subjlist_verified, **kwargs)
        case "merge":
            load_merge_all(subjlist_verified, **kwargs)
        case "plot":
            assert kwargs['pre_merged'], "Need to merge flow and power csv files first!"
            merged_dfs = pd.read_csv(kwargs['merged_fname'])
            plot_subjects(merged_dfs, **kwargs)

if __name__ == '__main__':
    with logger.catch():
        main(sys.argv[1:])