#system I/O
import os
import h5py
from pathos.pools import ProcessPool

#data packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

CONTEXT = 'paper'

#specialty 
from crp import reparam_trial
import resultAggregator as agg
import time

#Debugging
from loguru import logger
import pdb
import warnings
warnings.filterwarnings('ignore')



        #plot
        # print("Plotting data") 
        # plot_cross_project(S, pathout,subj, ma, stim, contact)
        
        # trial_reparam_df = reparam_trial(V_tr, canonical_response, tr_win)
        # out_f = os.path.join(pathout, f'figs/{subj}_reparam_trials_{contact}_stim_{stim}_{ma}.pdf')
        # plot_reparam_trials(trial_reparam_df, k, out_f)
        # out_f = os.path.join(pathout, f'figs/{subj}_reparam_agg_{contact}_stim_{stim}_{ma}.pdf')
        # plot_reparam_agg(trial_reparam_df,out_f)

        # ## on norm data
        # norm = np.linalg.norm(V_tr, axis=1)
        # V_norm =V_tr/ norm[:,None]
        # #TODO look at this tomorrow
        # 
        # norm_reparam_df = reparam_trial(V_norm, canonical_response, tr_win)
        # out_f = os.path.join(pathout, f'figs/{subj}_reparam_NORM_agg_{contact}_stim_{stim}_{ma}.pdf')
        # plot_reparam_agg(norm_reparam_df,out_f, proc='NORMED')
PLOT_COLS = ["fname","plot_type","key","out_fname","notes"]


def plot_channels(spes_df, channel_list, out_f=''):

    nrows = len(channel_list)
    if nrows >10:
        nrows = 10
    channel_list = gen_plot_channels(channel_list, nrows)

    with sns.plotting_context(CONTEXT):
        fig, axes = plt.subplots(nrows=nrows, ncols=1,sharex=True)
        plt.subplots_adjust(left=0.3,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=1)

        for i,ch in enumerate(channel_list):
            ax = axes[i]
            sns.lineplot(x=spes_df.index, y=spes_df[ch], ax=ax)
        fig.suptitle("Average SP Resp")
        plt.yticks(rotation=40)
        if out_f != "":
            plt.savefig(out_f, transparent=True)
    plt.close()


def gen_plot_channels(ch_list:list, num_samps:int):
    """Checks channels from spes_df for any errant columns (trial for example)
    If more than 10 channels requestes, automatically randomly samples to 10

    Args:
        ch_list (list): channels to plot, should be specified in spes datafram
        num_samps (int): number of desired samples, max 10

    Returns:
        np.array: array of strings, resampling channels
    """
    ch_list = list(ch_list)
    if "trial" in ch_list:
        ch_list.remove("trial")
    n_ch = len(ch_list)
    if n_ch <= num_samps:
        return ch_list
    inds = np.random.randint(0, n_ch, num_samps)
    ch_array = np.array(ch_list)
    return ch_array[inds]

def plot_cross_project(S, out_f, ma, stim, contact, tr_win):

    with sns.plotting_context(CONTEXT):
        ax = sns.lineplot(data=S, x="win_size", y="cross_proj")
        ax.axvspan(tr_win[0], tr_win[1], color='orange', alpha=0.5)
        plt.title(f"Cross Proj for {contact} resp to {stim}_{ma} mA")
        plt.legend( bbox_to_anchor=[0.15, 0.5], loc='right')
        plt.savefig(out_f, transparent=True)
        plt.close()



def plot_reparam_trials(trial_reparam_df, k,out_f):
    my_colors = ['k','r','g']

    with sns.plotting_context(CONTEXT):
        fig, axes = plt.subplots(nrows=k, ncols=1,sharex=True)
        for i in range(k):
            axis=axes[i]
            cols = [f'raw_{i}', f'proj_{i}', f'epsilon_{i}','time']
            ax = trial_reparam_df[cols].plot(kind='line',x='time',color=my_colors, ax=axis,legend=i==k)
            plt.ylabel('voltage')
            plt.xlabel("time")
            ax.get_legend()
        plt.legend( bbox_to_anchor=[1.15, 0.5], loc='center', labels=['Raw','Proj','Epsilon'])
        fig.suptitle("CRP Reparamaterization")
        plt.savefig(out_f,transparent=True)
        plt.close()

def plot_reparam_agg(trial_reparam_df, out_f,resp,stim ):
    my_colors = ['r']*10+['g']*10+['k']*10 +['y'] #TODO fix magic number

    with sns.plotting_context(CONTEXT):
        plt. figure(figsize=(10, 12)) 
        ax = trial_reparam_df.plot(kind='line',x='time',  color=my_colors)
        plt.legend( bbox_to_anchor=[1.15, 0.5], loc='center')
        plt.title(f"Reparam {resp} response, stim: {stim}")
        plt.savefig(out_f,transparent=True)
        plt.close()


def plot_row(row):

    fname = row['fname']
    plot_type = row['plot_type']
    key = row['key']
    out_f = row['out_fname']
    start = time.time()
    match plot_type:
        case "raw":
            df = pd.read_csv(fname,index_col=0)
            plot_channels(df, df.columns, out_f)
        case "reparam-agg":
            reparam_df = get_reparam(fname, key)
            resp = key.split("_")[-1]
            stim = fname.split("/")[-2]
            plot_reparam_agg(reparam_df,out_f,resp, stim)
        case "cross":
            S, args = get_crossproject(fname, key)
            plot_cross_project(S, out_f, *args)
    t = time.time() - start
    t = "%.3f" % t
    logger.info(f"Plotting {plot_type} Took :{t}s")


def get_reparam(fname, key):
    with h5py.File(fname, 'r') as h5:
        trial = h5[key]
        V_tr = trial['V_tr'][:]
        crp = trial['crp'][:]
        tr_win = trial.attrs['tr_win']
    return reparam_trial(V_tr, crp, tr_win)

def get_crossproject(fname, key):
    with h5py.File(fname, 'r') as h5:
        trial = h5[key] # key should get you to the stimtrial/resp-reg level
        stim = trial.attrs['stim']
        ma = trial.attrs['ma']
        contact = trial.attrs['contact']
        tr_win = trial.attrs['tr_win']
    S_key = os.path.join(key, "cross_proj")
    S = pd.read_hdf(fname,S_key)
    return S, [ma, stim, contact, tr_win]


#@logger.catch
def verify_df(df: pd.DataFrame):
    """data validation routine that ensures the plot_df is formatted correctly
    example df: 
                |    fname       |  plot_type   |        key            |     out_fname      |  note     |
                |'trial.csv'     | 'raw'        | ''                    |'avg-STIM-resp.pdf' |  non-sig  |
                | /path/to/h.hdf5| 'reparam-agg'| 'stim-sesh/resp-reg/' | 'agg-reparam.pdf   |  artifact |



    Args:
        df (pd.DataFrame): datagrame of files to plot and necessary information to call plotting
        functions. Must conform to example in description
    """
    assert (df.columns == PLOT_COLS).all(),'Improper column formatting'
    #check file names
    fnames = df['fname']
    fnames_check = fnames.str.contains('hdf5') | fnames.str.contains('.csv')
    assert fnames_check.all(), "Improper input filenames: must be csv or hdf5"


    plots = df['plot_type']
    plots_check = plots.str.contains('raw') | plots.str.contains('reparam-agg') | plots.str.contains('reparam') | plots.str.contains('cross')
    assert plots_check.all(), "Improper plot type specified!"


    out_fnames = df['out_fname']
    out_fnames_check = out_fnames.str.contains('pdf')
    assert out_fnames_check.all(), "Plots must be saved to .pdf"

    h5_entries = df[df.fname.str.contains(".hdf5")]
    if h5_entries.shape[0] == 0: 
        return #nothing to check!
    h5_check = h5_entries.keys != ""
    assert h5_check, "H5 files need keys to plot!"
    

def visualize_pipeline(plot_file:str)->None:
    """Given a file formatted properly, visualize pipeline will coordinate plotting derivatives and results
    from the CRP pipeline

    Args:
        plot_file (str): csv file formatted to load other csv intermediates or h5 files
    """
    assert ".csv" in plot_file, "Must be a csv file"

    plot_df = pd.read_csv(plot_file)
    verify_df(plot_df)
    # for _, row in plot_df.iterrows():
    #     plot_row(row)
    # plot_df.apply(plot_row, axis=1)
    rows = [row for _,row in plot_df.iterrows()]
    pool = ProcessPool(nodes=7)
    pool.map(plot_row, rows)


def gen_plot_df(subj: str, h5file: str, plot_types=['reparam-agg'], **kwargs) -> pd.DataFrame:
    """Helper function that generates a .csv with necessary information for crplot to generate plots
    Useful for significance plotting, artifact detection and result chcking

    Args:
        subj (str) : subject id 
        h5file (str): hdf5 file where keys can be found for generating plots
        kwargs mostly for plotting options -sig testing and artifact testing
        plot_types (list[str]): list of strings for plot types to specify
    Returns:
        pd.DataFrame : df with each row in correct plotting format for selected keys
    """
    artifact_test = kwargs['artifact_test'] if 'artifact_test' in kwargs.keys() else False
    sig_test = kwargs['sig_test'] if 'sig_test' in kwargs.keys() else False
    plot_path = kwargs['plot_path'] if "plot_path" in kwargs.keys() else ""
    plot_dfs = []
    with h5py.File(h5file, 'r') as f:
        notes = ""
        for key in f.keys():
            if sig_test:
                sig = agg.get_sig(key, h5file, f[key])
                notes = f"sig:{sig}-"
            if artifact_test:
                artifact_test = False #agg.artifact_test
                notes = f"{notes}artifact:{artifact_test}" #feel like I should make a JSON for this
            out_fnames = [gen_fname(p, subj, h5file, key,plot_path) for p in plot_types]
            rows = [[h5file, plot_types[i],key,out_f,notes] for i, out_f in enumerate(out_fnames)]
            df = pd.DataFrame(columns=PLOT_COLS, data=rows)
            plot_dfs.append(df)
    plot_dfs = pd.concat(plot_dfs)
    return plot_dfs

def gen_fname(plt_type:str,subj:str, file:str, key:str,plot_path='')->str:
    """Generate an output file name given information provided

    Args:
        plt_type (str): Type of plot desired, raw, cross project etc
        subj (str): Subject ID
        file (str): filename of hdf5 should contain stim info
        key (str): key for response region should have response contact info

    Returns:
        str: save file name with .pdf
    """
    stim = file.split("/")[-2]
    resp = key.split('_')[-1].replace(" ", "")
    match plt_type:
        case 'cross':
            plot = 'cross-proj'
        case 'reparam-agg':
            plot = 'reparam-agg'
        case 'reparam-trial':
            plot = 'reparam-trial'
        case 'raw':
            plot = 'raw-averagedERP'
            file = f"{subj}-{stim}-{plot}.pdf"
    file = f"{subj}-{stim}-{resp}-{plot}.pdf"
    return os.path.join(plot_path, file)