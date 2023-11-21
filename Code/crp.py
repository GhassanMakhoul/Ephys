#I/O  Setup
import os
import sys
import getopt
import re
import h5py

#data and math packages
import numpy as np
import pandas as pd
from scipy.io import loadmat
#plotting
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

# clean up pipeline IO
import warnings
warnings.filterwarnings('ignore')
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

def plot_channels(spes_df , channel_list,save=False,fout=''):
    nrows = len(channel_list)
    fig, axes = plt.subplots(nrows=nrows, ncols=1,sharex=True)
    for i,ch in enumerate(channel_list):
        ax = axes[i]
        sns.lineplot(x=spes_df.index, y=spes_df[ch], ax=ax)
    fig.suptitle("Average SP Resp")

    if save:
        plt.savefig(fout, transparent=True)
    plt.close()
    
def t_ix(n_samp,fs=512):
    return round(fs*n_samp)

#TODO remove IO/ file prep from main pipelin and add to runner
def assemble_trial(subj, stim_pair, ma):
    files = []
    spes_dfs = []
    for pulse in range(1,11):
        file = f'{subj}/{subj}_CCEP_single_pulses/{subj}_{stim_pair}_{ma}_pulse_{pulse}.mat'
        if not os.path.isfile(file):
            continue #check to see if file exists then skip if not
        print(f"Loading {file}")
        spes_trial = loadmat(os.path.join(DATA_FOLDER, file))
        fs = spes_trial['fs'][0][0]
        full_train = spes_trial["pulse"]
        labels  = [l[0] for l in spes_trial['labels'][0]]
        df = construct_spes_df(full_train, labels,fs, CAR) #TODO find more pipeline way
        df['trial'] =pulse
        spes_dfs.append(df)
    spes_dfs = pd.concat(spes_dfs)
    return spes_dfs, fs

def get_regions(contacts):
    return set([re.sub("[0-9]+",'', c) for c in contacts])

def cross_project_trial(V_full, fs, step=2):
    s =0 #paper assumes clean stim after 15ms
    n_samps, k = V_full.shape
    assert n_samps >k, 'incorrect data shape or not enough samples!'
    cross_dfs = []
    import pdb

    #TODO consider adding area relationships
    for e in np.arange(10,n_samps,step):
        #e = t_ix(t)
        win_len = e/fs
        V_raw = V_full[s:e,:]
        # norm length of time
        norm = np.linalg.norm(V_raw, axis=0)
        V_norm =V_raw/ norm[None,:]

        full_crossproj = V_norm.T@V_raw
        full_crossproj = full_crossproj - np.diag(full_crossproj)*np.identity(k)
        full_crossproj /= np.sqrt(fs)

        df = pd.DataFrame(data=full_crossproj.flatten(), columns=['cross_proj'])
        df['win_size'] = win_len
        cross_dfs.append(df)
    cross_proj_df = pd.concat(cross_dfs)
    return cross_proj_df

def flatten_df(og_df, cols, flat_name, const_col):
    n,_ = og_df.shape
    flat_data = og_df[cols].values.flatten()
    flat_df = pd.DataFrame(data=flat_data,columns=['voltage'])
    flat_cols = []
    for c in cols:
        flat_cols = flat_cols + [c]*n
    flat_df[flat_name] = flat_cols

    flat_df[const_col] = np.tile(og_df[const_col],len(cols))
    return flat_df


def plot_reparam(reparam_df,k):
    assert "proj_0" in reparam_df.columns and 'epsilon_0' in reparam_df.columns, "Incorrect columns!"
    fig, axes = plt.subplots(nrows=k, ncols=1, sharex=True)
    for i in range(k):
        cols = [f'proj_{i}', f'epsilon_{i}', 'crp']
        ax = axes[i]
        df = flatten_df(reparam_df,cols,'type','time')
        sns.lineplot(data =df,y='voltage', x='time', hue='type', ax=ax)

def reparam(V_data,crp):
    n = crp.shape[0]
    d = min(V_data.shape)
    alphas = crp.T@V_data
    proj = crp.reshape(n, 1)@alphas.reshape(1,d)
    ep = V_data - proj
    return alphas, proj, ep

def lin_PCA(V):
    """returns the result of lin kernal PCA"""
    [S2,eigen_vectors] = np.linalg.eig(V.T@V,)
    S = np.sqrt(S2)
    sort_inds = np.argsort(S)[::-1]
    S = S[sort_inds]
    eigen_vectors = eigen_vectors[:,sort_inds]
    V_kernel = V@eigen_vectors
    return V_kernel / S[None,:]

def plot_cross_project(S, pathout, subj, ma, stim, contact):

    with sns.plotting_context("paper"):
        sns.lineplot(data=S, x="win_size", y="cross_proj")
        fout = os.path.join(pathout, f'figs/{subj}_cross_proj_{contact}_{stim}_{ma}.pdf')
        plt.title(f"Cross Proj for {contact} resp to {stim}")
        plt.legend( bbox_to_anchor=[1.15, 0.5], loc='right')
        plt.savefig(fout, transparent=True)
    plt.close()
def get_tr(cross_proj_df):
    m = cross_proj_df.groupby(by='win_size').mean()
    max_val = np.max(m['cross_proj'].values)

    tr_win = m[m.cross_proj == max_val].index[0] #win_size will be index here
    #TODO change back
    tr_ind = t_ix(tr_win) #t_ix(tr_win+.15)
    return tr_ind, tr_win

def reparam_trial(V_tr, canonical_response, tr_win ):
    n_t, k = V_tr.shape 
    spes_trial = pd.DataFrame(data=V_tr, columns=[f'raw_{i}'for i in range(k)])
    t = np.linspace(0,tr_win ,n_t)
    spes_trial['time'] = t

    alphas, proj, ep = reparam(V_tr,canonical_response)

    df1 = pd.DataFrame(data=proj,columns=[f'proj_{i}'for i in range(k)])
    df1['time'] = t

    df2 = pd.DataFrame(data=ep, columns=[f'epsilon_{i}' for i in range(k)])
    df2['time'] = t

    trial_reparam_df = df1.merge(df2,on='time')
    trial_reparam_df = trial_reparam_df.merge(spes_trial, on='time')
    trial_reparam_df['crp'] = canonical_response
    return trial_reparam_df


def plot_reparam_trials(trial_reparam_df, k,fout):
    my_colors = ['k','r','g']
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
    plt.savefig(fout,transparent=True)
    plt.close()

def plot_reparam_agg(trial_reparam_df, fout, proc='RAW'):
    my_colors = ['r']*10+['g']*10+['k']*10 +['y']
    trial_reparam_df.plot(kind='line',x='time',  color=my_colors)
    plt.legend( bbox_to_anchor=[1.15, 0.5], loc='center')
    plt.title("Reparamaterization on {proc} voltage")
    plt.savefig(fout,transparent=True)
    plt.close()

def run_crp_pipeline(subj, pathout, ma, stim):
    spes_df, fs = assemble_trial(subj, stim, ma)
    spes_df.to_csv(os.path.join(pathout,f'spes_{stim}_{ma}.csv'))
    
    print("Saving Average Responses")
    #for reg in tqdm(get_regions(spes_df.columns[0:-1])):
    #    cols = [c for c in spes_df.columns if reg == re.sub('[0-9]+','',c)]
    #    fout = os.path.join(pathout,f"figs/{subj}_avg_sp_resp_{stim}_stim_{reg}.pdf")
    #    plot_channels(spes_df, cols,save=True,fout=fout)

    print("Running Cross Projection!")
    for contact in tqdm(spes_df.columns[0:-1]):
        #calc
        V_trial = spes_df.pivot( columns='trial', values=contact)
        S = cross_project_trial(V_trial.values,fs)
        tr_ind, tr_win = get_tr(S)
        V_tr = V_trial.values[0:tr_ind,:]
        eigenV = lin_PCA(V_tr)
        n_t, k = eigenV.shape
        canonical_response = eigenV[:,0] #get top PC for crp

        #plot
        # print("Plotting data") 
        # plot_cross_project(S, pathout,subj, ma, stim, contact)
        
        # trial_reparam_df = reparam_trial(V_tr, canonical_response, tr_win)
        # fout = os.path.join(pathout, f'figs/{subj}_reparam_trials_{contact}_stim_{stim}_{ma}.pdf')
        # plot_reparam_trials(trial_reparam_df, k, fout)
        # fout = os.path.join(pathout, f'figs/{subj}_reparam_agg_{contact}_stim_{stim}_{ma}.pdf')
        # plot_reparam_agg(trial_reparam_df,fout)

        # ## on norm data
        # norm = np.linalg.norm(V_tr, axis=1)
        # V_norm =V_tr/ norm[:,None]
        # #TODO look at this tomorrow
        # norm_reparam_df = reparam_trial(V_norm, canonical_response, tr_win)
        # fout = os.path.join(pathout, f'figs/{subj}_reparam_NORM_agg_{contact}_stim_{stim}_{ma}.pdf')
        # plot_reparam_agg(norm_reparam_df,fout, proc='NORMED')

        #save:
        print("Packaging out Results and Derivative Data")
        ## package out derivatives and data 
        
        writeout(eigenV,S, V_tr, fs, tr_win, pathout, subj, stim, contact, ma)

def writeout(basis, S, V_tr, fs, tr_win, pathout: str, subj: str, stim :str, contact :str, ma :str):
    """writeout deriviative data to an hdf5 file. Follows the following hierearchy
    filename is subject, the n create groups

    hdf5['/stim_sesh/resp_contact/...]
        group attrs : {fs: sampling rate, tr_win: (ix1, ix2), Tr: resp (s)}

        datasets:
        .../V_tr 
        .../S {fs: sampling rate, tr_win: (ix1, ix2), Tr: resp(s)}
        .../crp
        .../alphas
        .../epsilons
        .../alphas_prime


    Args:
        basis ( nd.array): canonicl responses  from linear kernel pca
        S (numpy.array): seminormalized cross projections of resp trials
        V_tr (numpy.array): windowed trials matrix TR (numsamps) x k (num trials)
        fs (float) : sampling rate
        t_win (numpy.array): indices marking response duration
        athout (str): folder to saveh hdf5 a
        subj (str): 
        stim (str): 
        contact (str):
        ma (str): _
    """
    crp = basis[:,0]
    alphas, proj, ep = reparam(V_tr, crp)
    h5f = os.path.join(pathout,'stim_resp.hdf5')
    key = f'/response_{contact}/'
    try:
        S.to_hdf(h5f, key=os.path.join(key,'cross_proj'), mode='a')
    except:
        print(f"Going to fail on this key: {key}")
        S.to_hdf(h5f, key=os.path.join(key,'cross_proj'), mode='a')
    #     raise(ValueError, f"failed to save cross proj of {key}")

    with h5py.File(h5f, 'a') as f:
        grp  =f[key]      
        
        grp.attrs['fs'] = fs
        grp.attrs['tr_win'] = (0,tr_win) #assumes starting at 0 with clean stim artifact
        grp.attrs['Tr'] = t_ix(tr_win,fs)


        dset = grp.require_dataset('V_tr',V_tr.shape, float)
        dset[:] = V_tr
        dset = grp.require_dataset('crp',crp.shape,float)
        dset[:] = crp
        dset =  grp.require_dataset('alphas', alphas.shape, float)
        dset[:] = alphas
        dset = grp.require_dataset('epsilons', ep.shape, float)
        dset[:] = ep
        dset = grp.require_dataset("projections", proj.shape, float)
        dset[:] = proj

def main(argv):
    subj = ''
    pathout =''
    ma = ''
    stim = ''
    opts, _ = getopt.getopt(argv,"s:t:m:p:",["subj=",'stim=',"ma=",'pathout'])
    for opt, arg in opts:
        if opt in ("-m", 'ma'):
            ma = arg
        elif opt in ("-s", "--subj"):
            subj = arg
        elif opt in ("-t", "--stim"):
            stim = arg
        elif opt in ('-p', '--pathout'):
            pathout = arg

    run_crp_pipeline(subj, pathout, ma, stim)


if __name__ == "__main__":
    ## load files
    ## bipolar montage them
    ## assemble trial per relationship
    main(sys.argv[1:])
