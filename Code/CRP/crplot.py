#system I/O
import os
import h5py
from loguru import logger

#data packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#specialty 
from crp import reparam_trial


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
        # norm_reparam_df = reparam_trial(V_norm, canonical_response, tr_win)
        # out_f = os.path.join(pathout, f'figs/{subj}_reparam_NORM_agg_{contact}_stim_{stim}_{ma}.pdf')
        # plot_reparam_agg(norm_reparam_df,out_f, proc='NORMED')


def plot_channels(spes_df , channel_list,out_f=''):

    nrows = len(channel_list)
    if nrows >10:
        nrows = 10
    channel_list = gen_plot_channels(channel_list, nrows)

    with sns.plotting_context("poster"):
        fig, axes = plt.subplots(nrows=nrows, ncols=1,sharex=True)
        for i,ch in enumerate(channel_list):
            ax = axes[i]
            sns.lineplot(x=spes_df.index, y=spes_df[ch], ax=ax)
        fig.suptitle("Average SP Resp")

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
    n = len(ch_list)
    if "trial" in ch_list:
        ch_list.remove("trial")
    if num_samps > n:
        return ch_list
    inds = np.random.randint(0, n, num_samps)
    ch_array = np.array(ch_list)
    return ch_array[inds]

def plot_cross_project(S, out_f, ma, stim, contact, tr_win):

    with sns.plotting_context("poster"):
        ax = sns.lineplot(data=S, x="win_size", y="cross_proj")
        ax.axvspan(tr_win[0], tr_win[1], color='orange', alpha=0.5)
        plt.title(f"Cross Proj for {contact} resp to {stim}_{ma} mA")
        plt.legend( bbox_to_anchor=[0.15, 0.5], loc='right')
        plt.savefig(out_f, transparent=True)
        plt.close()



def plot_reparam_trials(trial_reparam_df, k,out_f):
    my_colors = ['k','r','g']

    with sns.plotting_context("poster"):
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

def plot_reparam_agg(trial_reparam_df, out_f, proc='RAW'):
    my_colors = ['r']*10+['g']*10+['k']*10 +['y'] #TODO fix magic number

    with sns.plotting_context("poster"):
        trial_reparam_df.plot(kind='line',x='time',  color=my_colors)
        plt.legend( bbox_to_anchor=[1.15, 0.5], loc='center')
        plt.title("Reparamaterization on {proc} voltage")
        plt.savefig(out_f,transparent=True)
        plt.close()


def plot_row(row):
    fname = row['fname']
    plot_type = row['plot_type']
    keys = row['keys']
    out_f = row['out_fname']

    match plot_type:
        case "raw":
            df = pd.read_csv(fname,index_col=0)
            plot_channels(df, df.columns, out_f)
        case "reparam-agg":
            reparam_df = get_reparam(fname, keys)
            plot_reparam_agg(reparam_df,out_f)
        case "cross":
            S, args = get_crossproject(fname, keys)
            plot_cross_project(S, out_f, *args)

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


@logger.catch
def verify_df(df: pd.DataFrame):
    """data validation routine that ensures the plot_df is formatted correctly
    example df: 
                |    fname       |  plot_type   |        key           |     out_fname      |
                |'trial.csv'     | 'raw'        | ''                    |'avg-STIM-resp.pdf' |
                | /path/to/h.hdf5| 'reparam-agg'| 'stim-sesh/resp-reg/' | 'agg-reparam.pdf   |



    Args:
        df (pd.DataFrame): datagrame of files to plot and necessary information to call plotting
        functions. Must conform to example in description
    """
    assert (df.cols == ['fname', 'plot_type','key','out_fname']).all(),'Improper column formatting'
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
    h5_check =h5_entries.keys != ""
    assert h5_check.all(), "H5 files need keys to plot!"
    

def visualize_pipeline(plot_file:str)->None:
    """Given a file formatted properly, visualize pipeline will coordinate plotting derivatives and results
    from the CRP pipeline

    Args:
        plot_file (str): csv file formatted to load other csv intermediates or h5 files
    """
    assert ".csv" in plot_file, "Must be a csv file"

    plot_df = pd.read_csv(plot_file)
    verify_df(plot_df)
    plot_df.apply(plot_row)