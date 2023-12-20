#I/O  Setup
import os
import sys
import getopt
import re
import h5py
import glob
from tqdm import tqdm
from loguru import logger
import yaml


#data and math packages
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import ttest_1samp

#specialty

def agg_sesh_df(h5file, **kwargs):
    """aggregate across stim sesh, each key in for loop is diff stim session

    Args:
        h5file (str): hdf5 file for stim sesh

    Returns:
        pandas.DataFrame: concatenation of all stim trials
    """
    sig_test = kwargs['sig_test'] if "sig_test" in kwargs.keys() else False
    sesh_df =[]
    with h5py.File(h5file, 'r') as f:
        for key in f.keys():
            if sig_test and not get_sig(key, h5file, f[key]):
               continue
            df = entry_to_df(key, f[key])
            sesh_df.append(df)
        rej =len(f.keys()) - len(sesh_df)
        if rej >0 :
            logger.info(f"\t\t\tRejected {rej} for resp: {key}")

    sesh_df = pd.concat(sesh_df)
    return sesh_df


def get_explained_var(pulse_trial) -> np.number:
    """Compute explained variance of CRP parameterization
        1 - (e.T@e)/V_tr.T @ V_tr

    Args:
        pulse_trial (h5 group): response_reg group should contain, V_tr, and epsilon

    Returns:
        np.number: dimensionless representation of explained 
        variance of the first canonical response
    """
    V_tr = pulse_trial['V_tr'][:]
    e = pulse_trial["epsilons"][:]

    ev = 1 - np.divide(e.T @ e, V_tr.T @ V_tr)
    return ev.diagonal()



def get_sig(key:str, file_path:str, resp_h5:h5py.File) -> bool:
    """Computes significance using cross projection and returns TRUE/FALSE 
    depending on if cross proj is significant across Tr

    Args:
        key : str - stim resp region for loading S
        resp_h5 (h5py.File): h5 file obj spec to response

    Returns:
        bool: whether response is significant within TR
    """
    TR = resp_h5.attrs['Tr']
    # import pdb
    # pdb.set_trace()
    key = os.path.join(key, 'cross_proj')
    S = pd.read_hdf(file_path,key)
    S_tr = S[S.win_size ==TR]
    stats_inds = get_stats_inds(S_tr)
    vals = S_tr.iloc[stats_inds].cross_proj
    stats = ttest_1samp(vals, 0)
    return stats[1] < .05

def get_stats_inds(S):
    N,_ = S.shape
    inds = np.arange(0,N,2)
    return inds #TODO double check this, may be diff with flattened implementation

def entry_to_df(key, resp_h5):
    """assembles df with 
    1. resp_region
    2. alphas
    3. TR
    4. get explained variance
    """
    alphas = resp_h5['alphas']
    TR = resp_h5.attrs['Tr']
    fs = resp_h5.attrs['fs']
    df = pd.DataFrame(data=alphas, columns=['alphas']) #TODO check shape
    df['TR'] = TR
    df['resp_reg'] = key.split("_")[-1] #messy but keyshould be 'response_RH14' for example
    df['alpha_prime'] = alphas/(TR*fs) #TODO check on fs
    df['explained_variance'] = get_explained_var(resp_h5)
    return df



def get_stim_folders(subj: str, res_folder: str):
    """Returns fodlers of each stim trial using glob.glob's fuzzy pattern matching

    Args:
        subj (str): subject ID, e.g. 'Spat30',
        res_folder (str): path to result from crp pipeline
    """
    inp_pattern = os.path.join(res_folder,subj,'*-*_*mA' )
    folders = glob.glob(inp_pattern)

    #TODO sanity checks
    return folders


def get_sesh_params(folder:str):
    """Returns stim contacts and ma

    Args:
        folder (str): folder should be of form "path/to/Con1_ConP2_#ma"
    """
    stim_ma = folder.split("/")[-1]
    return stim_ma.split("_")

def agg_responses(subj: str, h5file: str, stim_folders: list, pathout: str, **kwargs):
    """aggregates all responses across all stim trials in to one
    mega dataframe

    Args
        subj (str): subject id
        h5file (str): name of hdf5 file should, eg 'stim_resp.hdf5'
        stim_folders (list): list of stim response dirs, full path
        pathout (str): where to save large df
    """
    dfs = []
    logger.info(f"Aggregating subj{subj} over {len(stim_folders)} folders")

    for folder in tqdm(stim_folders):
        h5f = os.path.join(folder, h5file)
        df = agg_sesh_df(h5f, **kwargs)
        stim_reg, ma  = get_sesh_params(folder)
        df['stim_reg'], df['ma'] = stim_reg, ma
        # logger.info(f"\tAgg: {stim_reg} at {ma}")
        dfs.append(df)
    agg_df = pd.concat(dfs)
    agg_df['subj'] = subj
    agg_df.to_csv(os.path.join(pathout, f"{subj}_stim.csv"),index=False)

def verify_pathout(pathout:str)->None:
    """Ensures that the path exists.

    Args:
        pathout (str): path to save files to
    """
    if not os.path.exists(pathout):
        os.makedirs(pathout)
    
def gen_plot_file(subj: str, h5file:str, stim_folders: list, pathout:str, **kwargs):
    """Creates a csv with plotting commands that crplot.py can use. kwargs determines whether to 
    employ artifact detection and significance detection. If artifact/significance detection enabled, a rejection file will be generated.  

    Args:
        subj (str): Subject's SPES files to aggregate
        h5file (str): name of CRP file
        stim_folders (list): folder with stim sesh
        kwargs (dict, optional): Plotting commands options, controls number of channels to plot(rand select
        Whether to detect artifatcts and insignificant files).
    """
    from crplot import gen_plot_df #avoids circular import
    verify_pathout(pathout)
    plot_dfs = []
    for folder in stim_folders:
        h5f = os.path.join(folder, h5file)
        df = gen_plot_df(subj, h5f, **kwargs)
        plot_dfs.append(df)
        #TODO gen plot max filter
    plot_dfs= pd.concat(plot_dfs)
    plot_dfs.to_csv(os.path.join(pathout, f"{subj}_plots.csv"),index=False) 
    

@logger.catch
def main(argv):
    subj = ''
    config_f = 'config_agg.yml'
    opts, _ = getopt.getopt(argv,"s:i:p:c:",["subj=",'inpf=','pathout=','config='])
    for opt, arg in opts:
        if opt in ("-i", 'inpf'):
            inpf = arg
        elif opt in ("-s", "--subj"):
            subj = arg
        elif opt in ('-p', '--pathout'):
            pathout = arg
        elif opt in ("-c", '--config'):
            config_f = arg
    with open(config_f, 'r') as f:
        config =  yaml.safe_load(f)


    res_folder = config['general']['res_folder']
    inpf = config['general']['h5filename']
    logdir = config['general']['logdir']
    
    logger.add(os.path.join(logdir, f'agg_{subj}.log'))
    stim_folders = get_stim_folders(subj, res_folder)
    agg_kwargs = config['agg'] if 'agg' in config.keys() else {}
    agg_responses(subj, inpf, stim_folders, pathout, **agg_kwargs)
    logger.success(f"SUCCESSFULLY aggregated responses for {subj}!")

    if config['plot']['gen_plots']:
        logger.info(f"Generating plot file for {subj}")
        plot_kwargs = config['plot']   #TODO should I feed yaml directly as kwargs? seems safer to process input first
        plot_kwargs['plot_path'] = plot_kwargs['plot_path'].replace("SUBJ", subj)

        gen_plot_file(subj, inpf, stim_folders, pathout,**plot_kwargs)
        logger.success(f"Generated Plot file for {subj}")


if __name__ == '__main__':
    main(sys.argv[1:])
