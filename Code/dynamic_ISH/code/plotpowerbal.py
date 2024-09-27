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

from collections import Counter, defaultdict
import pandas as pd
import numpy as np
from scipy.stats import f_oneway
from sklearn import linear_model
from statsmodels.formula.api import ols
import statsmodels.api as sms

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from utils import *
from connectivity_dynamics import center_onset

def stim_to_df(stim_struct):
    """Returns a flattened dataframe with
    stim site and response site as well as banded
    change in power from pre/post stim z-scoring
    for all gray matter response sites

    Args:
        stim_struct (_type_): stim struct with banded response values
        should be passed through at the level of pat_spes

    Returns:
        _type_: _data frame of stim responses
    """
    patID = stim_struct['patID_clean']
    stim_labels = stim_struct['stim_labels']
    resp_labels = stim_struct['response_labels']
    z_freq = stim_struct['data_RBP_Z_AllDists']
    resp_dict = defaultdict(lambda : [])
    for i, stim_site in enumerate(stim_labels):
        for j, band in enumerate(BANDS[2:]):
            resp_dict["stim_bip"].extend(stim_site*len(resp_labels))
            resp_dict['resp_bip'].extend([r[0] for r in resp_labels])
            resp_dict[f'change_z'].extend(z_freq[i,:,j])
            resp_dict['band'].extend([band]*len(resp_labels))
    stim_df = pd.DataFrame.from_dict(resp_dict)
    stim_df['patID'] = patID
    return stim_df

def merge_flow_spes(flow_df: pd.DataFrame, power_df:pd.DataFrame ,spes_df: pd.DataFrame, subsample=1, band='beta',**kwargs):
    """merge bipole level peri_ictal flow_df to banded spes channel-wise df
    For seizure involvement aggregates deisgnations across all flow_dfs"""
    spes_flow_dfs = []
    sz_bips_df = get_overall_involved(power_df)

    df = pow_event_df[['bip',f'z_{band}','win_sz_centered', f'{band}_involved',"z_delta","anat_region"]].merge(\
            flow_event_df,how='right',
            right_on=['win_sz_centered','src_bip'],
            left_on=['win_sz_centered','bip'])
    df['region_involved'] = df.apply( \
            lambda x: f"{x['source']}_{x['target']}_{x[f'{band}_involved']}",\
            axis=1)
    if subsample > 1:
        df = subsample_df(df, subsample, groups=['source','win_label','region_involved'])
        spes_flow_dfs.append(df)
    return

def get_involved_inds(event_df, band="beta", threshold=2):
    """Given an event_df, which contains the dynamics of all bipoles for one seizure event,
    thresholds all areas which achieve an average band power above threshold 
    and returns a column detailing which bipoles are involved as BANDS_increased"""
    ictal_df = event_df[event_df.win_label.isin(['early_ictal','ictal','late_ictal'])]
    ictal_bip_df = ictal_df[['bip', f'z_{band}',]].groupby('bip').max().reset_index()
    # get channels that on average achieve a z_score above  threshold power
    thresh_bool = ictal_bip_df[f'z_{band}'].values > threshold
    power_dict = dict(zip(ictal_bip_df.bip, thresh_bool))
    all_bip_thresh = event_df.bip.apply(lambda x: power_dict[x])
    event_df.insert(loc=len(event_df.columns), column=f"{band}_involved", value=all_bip_thresh)
    return event_df

def get_involved_dict(event_df, band="beta", threshold=2):
    """Given an event_df, which contains the dynamics of all bipoles for one seizure event,
    thresholds all areas which achieve an average band power above threshold 
    and returns a DICTIONARY detailing which bipoles are involved as BANDS_increased"""
    ictal_df = event_df[event_df.win_label.isin(['early_ictal','ictal','late_ictal'])]
    ictal_bip_df = ictal_df[['bip', f'z_{band}',]].groupby('bip').max().reset_index()
    # get channels that on average achieve a z_score above  threshold power
    thresh_bool = ictal_bip_df[f'z_{band}'].values > threshold
    power_dict = dict(zip(ictal_bip_df.bip, thresh_bool))
    return power_dict


def merge_flow_power(flow_df: pd.DataFrame, power_df: pd.DataFrame, subsample=1, band='beta',**kwargs):
    """merge bipole level peri_ictal flow_df to band power channel
    assumes that eventID is same types and that there is an equal number
    of events between the flow_df and the power_df"""
    power_flow_dfs = []
    for event in set(power_df.eventID.values):
        pow_event_df = power_df[power_df.eventID == event]
        pow_event_df = get_involved_inds(pow_event_df, band=band)
        flow_event_df = flow_df[flow_df.eventID == event]
        assert "win_sz_centered" in flow_event_df.columns, "Need to center flow first!"
        assert "win_sz_centered" in pow_event_df.columns, "Need to center flow first!"
   
        df = pow_event_df[['bip',f'z_{band}','win_sz_centered', f'{band}_involved',"z_delta","anat_region"]].merge(\
            flow_event_df,how='right',
            right_on=['win_sz_centered','src_bip'],
            left_on=['win_sz_centered','bip'])
        df['region_involved'] = df.apply( \
            lambda x: f"{x['source']}_{x['target']}_{x[f'{band}_involved']}",\
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
    logger.info(f"Original data frame has {n_rows} rows")
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



#TODO separate out connection types into certain plots
# 1. first plot Nz_soz_true, nz_soz_false, nz_nz_true
# 2. Z-score PDC against interictal period

def joint_kde_flow_power_plotter(flow_power_df, pltname, relationships='all', sz_types=[] ,anatomical =False, **kwargs):
    #NOTE holy shit recursion is useful sometimes
    if len(sz_types) > 1:
        sz_type_df = flow_power_df[flow_power_df.sz_type == sz_types[0]]
        sz_pltname = f"{sz_types[0]}_{pltname}"
        joint_kde_flow_power_plotter(sz_type_df, sz_pltname, relationships=relationships, sz_types=[],anatomical=anatomical, **kwargs)
        joint_kde_flow_power_plotter(flow_power_df, pltname, relationships=relationships, sz_types=sz_types[1:],anatomical=anatomical, **kwargs) # do subtyping here
    elif len(sz_types) == 1:
        flow_power_df = flow_power_df[flow_power_df.sz_type == sz_types[0]]
        sz_pltname = f"{sz_types[0]}_{pltname}"
    if relationships == 'all':
         #NOTE make more modular to accomodate any grouping ?
       joint_kde_flow_power(flow_power_df, pltname,anatomical=anatomical, **kwargs)
    else:
        for rel in relationships:
            #assyne that list of lists is passed for this
            plot_df = flow_power_df[flow_power_df.region_involved.isin(rel)]
            rel_pltname = f"{pltname}_{rel}"
            # assumes following formatting: "src_trgt_szBoolean"
            # TODO reference the sources from leftover df? check df in live sesh
            sources = list(plot_df.source.unique())
            joint_kde_flow_power(plot_df, rel_pltname, sources=sources,anatomical=anatomical, **kwargs)
    # if want to segment by seizure type and relatinoships
    # call joint_kde_flow_power on one seizure recursively, then call back to joint_kde_power with remaining list?

def kde_flow_power(flow_power_df, pltname, y='value',x='z_beta', title="", sources= ['nz', 'soz'],relationships="all", **kwargs):
     
     grp_flowpow_df = flow_power_df.groupby(['win_label','source'])
     with sns.plotting_context("paper"):
            fig, axs = plt.subplots(2,7, sharex=True, sharey=True)
            fig.set_figheight(15)
            fig.set_figwidth(30)
            for i, source in enumerate(sources):
                for t, period, in enumerate(['interictal', 'pre_ictal','early_ictal','ictal','late_ictal','early_post_ictal','post_ictal']):
                    ax = axs[i][t]
                    plot_df = grp_flowpow_df.get_group((period, source))
                    
                    if source == 'nz' and period == 'post_ictal':
                        ax = sns.kdeplot(data=plot_df, x=x, y=y, hue='region_involved', legend=True, ax=ax, palette=FLOWMAP)
                        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
                    elif source =='soz' and period =='post_ictal':
                        ax = sns.kdeplot(data=plot_df, x=x, y=y, hue='region_involved', legend=True, ax=ax, palette=FLOWMAP)
                        sns.move_legend(ax, "lower left", bbox_to_anchor=(2,2) )
                    else:
                        ax = sns.kdeplot(data=plot_df, x=x, y=y, hue='region_involved', legend=False, ax=ax, palette=FLOWMAP)
                    ax.title.set_text(f"{period} -Flow from {source} vs {x}")
            plt.suptitle(title)
            plt.savefig(f"../viz/{pltname}.pdf", transparent=True, bbox_inches="tight")

def joint_kde_flow_power(flow_power_df, pltname, y='value',x='z_beta', sources=['nz','soz'], title="", anatomical=False, **kwargs):
    if anatomical:
        logger.info("Adding Anatomical Designation")
        #kind of hacky should TODO: add this to merge flow part of pipeline
        flow_power_df.region_involved = flow_power_df.apply(lambda row: f"{row['region_involved']}_{'ctx' if 'tx' in row['anat_region']else 'outCTX'}",axis=1)
    grp_flowpow_df = flow_power_df.groupby(['win_label','source'])
    with sns.plotting_context("paper"):
            #TODO consider removing enumeration
            for i, source in enumerate(sources):
                for t, period, in enumerate(['interictal', 'pre_ictal','early_ictal','ictal','late_ictal','early_post_ictal','post_ictal']):
                    plt.figure(figsize=(15,30))
                    plot_df = grp_flowpow_df.get_group((period, source))
                    count_df = plot_df.groupby(by="region_involved").count().reset_index()

                    ax = sns.jointplot(data=plot_df, x=x, y=y, kind='kde', hue='region_involved', legend=True, palette=FLOWMAP, **kwargs)
                        #sns.move_legend(og_ax, "upper left", bbox_to_anchor=(1, 1))
                    # Add a text box
                    text_box = dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.5)
                    logger.info(f"Count_df for plots\n\n\t\t{count_df}")
                    
                    plt.xlim(-25, 40)
                    plt.ylim(0,0.5)
                    plt.suptitle(f"{period} -Flow from {source} vs {x}")
                    out_f = f"../viz/{pltname}_source_{source}_{period}.pdf"
                    plt.savefig(out_f, transparent=True, bbox_inches="tight")
                    logger.success(f"Saved {out_f}")
                    plt.close()


def scatter_flow_power(flow_power_df, pltname, y='value', x='z_beta', **kwargs):
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

#NOTE: I hate that I have to do this, should have z_scored verbose calculations
# at the initial compilation but for now, I'll post-hoc z-score
#TODO in future re-analysis add z-score verbose option

def  load_score_merge_subject(subject, power_dir="/mnt/ernie_main/Ghassan/ephys/data/ei_bal/", flow_dir="/mnt/ernie_main/Price/ephys/data/periconnectivity", **kwargs):
    """same as above except includes a step for z_scoring connectivity matrix at biipole level against
        the interictal period    """
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
    ## INSERT Z score here
    z_flow_df = z_score_all_events(flow_df)
    logger.success(f"Z-scored all {len(flow_df.eventID.unique())} events for {subject}")
    return merge_flow_power(z_flow_df, power_df, **kwargs)

def z_score_all_events(flow_df):
    """Iterate through all events and z=scure each
    bipoles connectivity values against the interictal period

    Args:
        flow_df (pd.DataFrame): df of verbose flow, channel level connectivity (in a source target format)
    """
    #for each event select the interictal period
    #       then group by bip get mean of 'value' (the PDC)
    #       get mean,std of bip connectivity values
    #       apply z_score to the remaining periods at the bipole level
    z_events = []
    for event in flow_df.eventID.unique():
        event_df = flow_df[flow_df.eventID == event]
        interictal_flow_df = event_df[event_df.win_label =='interictal']
        bip_df = interictal_flow_df[['src_bip', 'freq_band','value']].groupby(by=['src_bip', 'freq_band'])
        mu_df = bip_df.mean().reset_index()
        std_df = bip_df.std().reset_index()
        mu_dict, std_dict = defaultdict(lambda : np.nan), defaultdict(lambda : np.nan)
        # creates a dictionary mapping a (src bipole, connectivity_band) - > mu (or std) thus the double zip
        mu_dict.update(dict(zip(zip(mu_df.src_bip, mu_df.freq_band), mu_df.value))) # MESSY!
        std_dict.update(dict(zip(zip(std_df.src_bip, std_df.freq_band), std_df.value))) 
        # MESSY Again! but this is necessary to prevent edge cases from failing
        # some patients don't have entries for all bips at every frequency band

        apply_z = lambda row : (row['value'] - mu_dict[(row['src_bip'], row['freq_band'])])/std_dict[(row['src_bip'], row['freq_band'])]
        z_pdc = event_df.apply(apply_z, axis=1)
        event_df.insert(loc=len(event_df.columns), column='z_pdc', value=z_pdc)
        z_events.append(event_df)

    return pd.concat(z_events)


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
        #maybe assert that there is a file : UserWaruse joinning: Ignoring `ax`; jointplot is a figure-level function.
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


def load_merge_all(subjects, num_cores=20, merged_fname="merged_df", out_dir="../data", **kwargs):
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
    merged_fname = os.path.join(out_dir, merged_fname)
    merged_dfs.to_csv(merged_fname)
    logger.info(f"Saved merged df to {merged_fname}")

def run_stats(merged_df):
    """
    run 2 way anova"""
    logger.info("RUNNING STATS")
    cols = ['patID', 'bip',  'win_sz_centered','beta_involved', "win_label", 'freq_band','region_involved',  "value", 'z_beta']
    group_df =  merged_df[cols].groupby(cols[0:-2]).mean().reset_index()
    group_df['target'] = group_df.region_involved.apply(lambda x: x.split("_")[1])
    #TODO make bands more modular in future
    group_alpha = group_df[group_df.freq_band =='alpha']
    #TODO make modular for regionwise
    regions =  ['nz_nz_True', 'nz_nz_False', 'nz_soz_True', 'nz_soz_False']
    stats_df = group_alpha[group_alpha.region_involved.isin(regions)]
    for window in stats_df.win_label.unique():
        logger.info(f"Running {window}")
        win_stats_df = stats_df[stats_df.win_label == window]
        anova2way_model = ols("value~beta_involved+target", data=win_stats_df).fit()
        stats_res = sms.stats.anova_lm(anova2way_model, type=2)
        logger.success(f"Computed stats for {window} window. Results:\n{stats_res}")
    return
def load_score_merge_all(subjects, num_cores=20, merged_fname="merged_df", **kwargs):
    #may need to have spec test for certain centering 
    # if len(glob.glob(os.path.join(kwargs['power_dir'], "*_centered.csv"))) < len(subjects):
    #     logger.info("Centering Seizure events first!")
    #     center_seizure_onset(subjects, **kwargs)
    logger.info(f"Loading subjects {len(subjects)}")
    if len(subjects) <2 :
        #for debugging purposes
        merged_dfs =[ load_score_merge_subject(subj, **kwargs) for subj in subjects]
    elif len(subjects) >1:
        num_cores = min(len(subjects), num_cores)
        p = Pool(num_cores)
        merged_dfs = p.map(partial(load_score_merge_subject, **kwargs), subjects)
        p.close()
        p.join()
        merged_dfs = [m_df for m_df in merged_dfs if not m_df.empty]
    merged_dfs = pd.concat(merged_dfs)
    logger.success(f"Merged and loaded {len(subjects)}")
    merged_dfs.to_csv(merged_fname)
    logger.info(f"Saved merged df to {merged_fname}")

def plot_subjects(merged_df, pltname='flow_pow.pdf', band='beta', plot_type='kde', conn_band="all",plot_kwargs={}, norm_szrs=False, **kwargs):   
    kwargs.update(plot_kwargs)
    ## by default exclude unknown awareness
    merged_df = merged_df[merged_df.sz_type.isin(['FIAS','FBTC','FAS'])]
    #TODO maybe make x more modular
    if 'x' not in kwargs.keys():
        x = f"z_{band}"
    if conn_band != "all":
        merged_df = merged_df[merged_df.freq_band ==conn_band]
        pltname = f"{pltname}_{conn_band.upper()}"
    logger.info(f"About to plot peri-ictal {plot_type}")
    if norm_szrs:
        #NOTE should I add a groupby or just manually
        #get single seizure representation
        #TODO this is a brittle point, could make more modular in the future
        logger.info("Averaging within patient")
        cols = ['patID', 'bip',  'win_sz_centered','beta_involved', 'source', 'win_label', 'freq_band','region_involved',  "value", 'z_beta']
        merged_df = merged_df[cols].groupby(cols[0:-2]).mean().reset_index()
    match plot_type:
        case "kde":
            kde_flow_power(merged_df, pltname,x=x, **kwargs) 
        case "scatter":
            scatter_flow_power(merged_df, pltname, x=x, **kwargs)
        case "joint_kde":
            joint_kde_flow_power_plotter(merged_df, pltname, **kwargs)

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
    
    subjlist_verified = subjlist_verified
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

        case "zscore_merge":
            load_score_merge_all(subjlist_verified, **kwargs)
        case "plot":
            assert kwargs['pre_merged'], "Need to merge flow and power csv files first!"
            merged_dfs = pd.read_csv(kwargs['merged_fname'])
            plot_subjects(merged_dfs, **kwargs)
        case "stats":
            assert kwargs['pre_merged'], "Need to merge flow and power csv files first!"
            merged_dfs = pd.read_csv(kwargs['merged_fname'])
            run_stats(merged_dfs)

if __name__ == '__main__':
    with logger.catch():
        main(sys.argv[1:])