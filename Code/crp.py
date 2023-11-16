import mne
import numpy as np
import pandas as pd
import numpy as np
from scipy.io import loadmat
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import re
import pdb

DATA_FOLDER = '/mnt/ernie_main/000_Data/SPES/data/preprocessed/'

def construct_spes_df(spes_trains, contact_labels, fs,reref=None):
    assert len(contact_labels) == min(spes_trains.shape)
    if np.argmax(spes_trains.shape) != 0:
        spes_trains = spes_trains.T
    nsamps,ch = spes_trains.shape
    t = np.arange(nsamps)/fs
    df = pd.DataFrame(columns=contact_labels, index =t, data=spes_trains)
    if reref != None:
        reref(df)
    return  df

def CAR(df):
    regions = get_regions(df.columns)
    #TODO look for bad channels or exclude
    n_samps = max(df.shape)
    for reg in regions:
        #get data
        reg_cols = [c for c in df.columns if reg == re.sub('[0-9]+','',c)]
        #reassign
        mean = df[reg_cols].mean(1).values
        df[reg_cols] = df[reg_cols].values - np.reshape(mean, (n_samps,1))

    return df

def plot_channels(spes_df , channel_list):
    nrows = len(channel_list)
    fig, axes = plt.subplots(nrows=nrows, ncols=1,sharex=True)
    for i,ch in enumerate(channel_list):
        ax = axes[i]
        sns.lineplot(x=spes_df.index, y=spes_df[ch], ax=ax)
    
def t_ix(t,fs=512):
    return round(t*fs)
    
def assemble_trial(subj, stim_pair, ma):
    files = []
    spes_dfs = []
    for pulse in range(1,11):
        file = f'{subj}/{subj}_CCEP_single_pulses/{subj}_{stim_pair}_{ma}_pulse_{pulse}.mat'
        files.append(file)
        spes_trial = loadmat(os.path.join(DATA_FOLDER, file))
        fs = spes_trial['fs'][0][0]
        full_train = spes_trial["pulse"]
        labels  = [l[0] for l in spes_trial['labels'][0]]
        df = construct_spes_df(full_train, labels,fs, CAR) #TODO find more pipeline way
        df['trial'] =pulse
        spes_dfs.append(df)
    spes_dfs = pd.concat(spes_dfs)
    return spes_dfs

def get_regions(contacts):
    return set([re.sub("[0-9]+",'', c) for c in contacts])



def main(argv):
    subj = ''
    pathout =''
    ma = ''
    stim = ''
    opts, _ = getopt.getopt(argv,"s:st:m:p",["subj=",'stim=',"ma=",'pathout'])
    for opt, arg in opts:
        if opt in ("-m", 'ma'):
            ma = arg
        elif opt in ("-s", "--subj"):
            subj = arg
        elif opt in ("-st", "--stim"):
            stim = arg
        elif opt in ('-p', '--pathout'):
            pathout = arg
    spes_df = assemble_trial(subj, stim, ma)


if __name__ == "__main__":
    ## load files
    ## bipolar montage them
    ## assemble trial per relationship
    main(sys.argv[1:])