import mne
import numpy as np
import pandas as pd


def construct_spes_df(spes_trains, contact_labels, fs):
    assert len(contact_labels) == min(spes_trains.shape)
    if np.argmax(spes_trains.shape) != 0:
        spes_trains = spes_trains.T
    nsamps,ch = spes_trains.shape
    t = np.arange(nsamps)/fs
    df = pd.DataFrame(columns=contact_labels, index =t, data=spes_trains)
    return  df

def plot_channels(spes_df , channel_list):
    nrows = len(channel_list)
    fig, axes = plt.subplots(nrows=nrows, ncols=1,sharex=True)
    for i,ch in enumerate(channel_list):
        ax = axes[i]
        sns.lineplot(x=spes_df.index, y=spes_df[ch], ax=ax)
    
def t_ix(t,fs=512):
    return round(t*fs)
    
def assemble_trial():
     files = []
    spes_dfs = []
    for pulse in range(1,11):
        file = f'Epat26/Epat26_CCEP_single_pulses/Epat26_LAC3-LAC4_5mA_pulse_{pulse}.mat'
        files.append(file)
        spes_trial = loadmat(os.path.join(DATA_FOLDER, file))
        fs = spes_trial['fs'][0][0]
        full_train = spes_trial["pulse"]
        labels  = [l[0] for l in spes_trial['labels'][0]]
        df = construct_spes_df(full_train, labels,fs)
        df['trial'] = pulse
        spes_dfs.append(df)
    stim_contact = get_stim(files[0])
    spes_dfs = pd.concat(spes_dfs)

if '__name__' == __main__():
    ## load files
    ## bipolar montage them
    ## assemble trial per relationship