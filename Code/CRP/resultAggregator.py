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
from scipy.signal import find_peaks

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
            df = entry_to_df(key, f[key],**kwargs)
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


def get_times_to_peak(curve: np.array, fs: int, n_peaks =1, dist_time = .05) -> list[float]:
    """Using a peakfinding approach, this returns the index of the peak
    We are implementing this to characterize delays in response timing.
    The matsumoto2017 review on SPES posits that delays in peak response 
    characterize the SOZ

    Args:
        curve (np.array): any SPES curve, could be the CRP, could be an average, or a single 
        pulse, but needs to be an np.array
        
        fs (int): sampling rate - samples/second
        n_peaks (int) : number of peaks to examine. Defaults to 1.
        dist_time (float) : peaks need to be this many ms apart

        #TODO - implement a windowing period may help peakfinder with known N1, N2
        #TODO - clarify if we want to only find positive deflections. Given the nature of 
        non cortical SPES we may have some flipped responses.


    Returns:
        list[float]: list of times correlating to number of peaks
    """
    start_ind= int(.01*fs) #exclude first 10ms 
    peak_inds = np.argmax(curve[start_ind:])
    return [peak_inds/fs + .01]
    dists =   dist_time * fs #50 ms (default) = .010s * fs samp/s = n_samps
    # peak_inds, _ = find_peaks(curve,distance=dists, threshold=min_height)
    
    # if n_peaks < len(peak_inds):
    #     diff = n_peaks - len(peak_inds)
    #     np.append(peak_inds, [np.nan]*diff)
    # peak_inds = peak_inds[0:n_peaks]
    # peak_times =  [ind/fs  for ind in peak_inds if ind/fs]

    #return peak_times


def get_timing(resp_h5,n_peaks=1):
    fs = resp_h5.attrs['fs']
    trial_dim = np.argmin(resp_h5['projections'].shape)
    curve = np.mean(resp_h5['projections'],axis=trial_dim)
    return get_times_to_peak(curve, fs, n_peaks=n_peaks)

def add_timing(df, resp_h5, **kwargs):
    """Return dataframe with columns for timing of peaks"""

    n_peaks = kwargs['n_peaks'] if 'n_peaks' in kwargs.keys() else 1
    peak_times = get_timing(resp_h5,n_peaks=n_peaks)
    for i,t in enumerate(peak_times):
        df[f't_peak_{i}'] = t
    return df

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
    stats = ttest_1samp(vals, 0,alternative="greater")
    return stats[1] < .05

def get_stats_inds(S):
    N,_ = S.shape
    n = np.sqrt(N)
    inds = np.arange(0,N-n,2)
    if n %2 ==1:
        inds = tri_checkerboard(n)
    return inds #TODO double check this, may be diff with flattened implementation

def tri_checkerboard(n):
    """Returns every other index of the upper and lower triangular indices. This satisfies 
    2 properties of the stats sample for cross projection:
        1. Half of the samples reflect normalized trial projected onto un-norm
        2. No double counting of normed trials
        To satisfy this, the upper and lower triangular indices of the matrix need to
        to be masked with a checkerboard pattern. The checkerboard pattern is offset
        differntly for the upper and lower triangular matrices 


    Args:
        n (_type_): _description_
    """
    n = int(n)
    raw_inds = np.arange(0, n**2).reshape(n,n)

    u_inds = tri_mask(n)

    checker_upper = checkerboard((n,n), 1)
    checker_inds = checker_upper[u_inds]
    try:
        upper_samps = raw_inds[u_inds][checker_inds]
    except IndexError:
        import pdb
        pdb.set_trace()
    

    l_inds = tri_mask(n, side='lower')
    checker_lower = checkerboard((n,n), 0)
    checker_inds = checker_lower[l_inds]
    lower_samps = raw_inds[l_inds][checker_inds]
    return np.append(upper_samps, lower_samps)


def checkerboard(shape, offset):
    """Creates a checkerboard boolean indices

    Args:
        shape (int): shape of one dimenstion
        offset (int): 0 or 1 if 0 then 0,0 will be true

    Returns:
        np.ndarray : boolean matrix

    example:
    in : checkerboard((9,9), 0)
    
    out: array([[ True, False,  True, False,  True, False,  True, False,  True],
       [False,  True, False,  True, False,  True, False,  True, False],
       [ True, False,  True, False,  True, False,  True, False,  True],
       [False,  True, False,  True, False,  True, False,  True, False],
       [ True, False,  True, False,  True, False,  True, False,  True],
       [False,  True, False,  True, False,  True, False,  True, False],
       [ True, False,  True, False,  True, False,  True, False,  True],
       [False,  True, False,  True, False,  True, False,  True, False],
       [ True, False,  True, False,  True, False,  True, False,  True]])

    """
    try:
        return np.indices(shape).sum(axis=0) %2 == offset
    except TypeError:
        logger.warning("INDS SHOULD BE INTS")
        shape = ( int(shape[0]), int(shape[1]) )
        return np.indices(shape).sum(axis=0) %2 == offset

def tri_mask(n, side='upper'):
    """n should be the shape of the matrix

    Args:
        n (int): shape of
        side (str, optional): _description_. Defaults to 'upper'.

    Returns:
        np.ndarray: boolean array of upper triangular
    """
    r = np.arange(n)
    if side == 'upper':
        mask = r[:, None] > r
    else:
        mask = r[:, None] < r
    return mask


def entry_to_df(key, resp_h5, **kwargs):
    """assembles df with 
    1. resp_region
    2. alphas
    3. TR
    4. get explained variance
    5. add time delay
    """
    alphas = resp_h5['alphas']
    TR = resp_h5.attrs['Tr']
    fs = resp_h5.attrs['fs']
    df = pd.DataFrame(data=alphas, columns=['alphas']) #TODO check shape
    df['TR'] = TR
    df['resp_reg'] = key.split("_")[-1] #messy but keyshould be 'response_RH14' for example
    df['alpha_prime'] = alphas/(TR*fs) #TODO check on fs
    df['explained_variance'] = get_explained_var(resp_h5)

    df = add_timing(df, resp_h5,**kwargs)
    return df



def get_stim_folders(subj: str, res_folder: str):
    """Returns folders of each stim trial using glob.glob's fuzzy pattern matching

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


def agg_crp(subj: str, h5file: str, stim_folders: list, pathout: str, **kwargs): #parallel agg_reponses
    # open H5File (represents a stim setting) (all H5Files have same name but in separate folders)
    # for each stim response region -> pull out CRP
        # zero pad from 0 - 1000ms *f3
    # add whole CRP to dataframe
    # add stim_df to subj_df
    # Save out subj_df
    subj_df = []
    for folder in tqdm(stim_folders):
        h5f = os.path.join(folder, h5file)
        stim_reg, ma = get_sesh_params(folder)
        #create df that is fully zero padded
        stim_df = agg_crp_df(h5f, **kwargs ) # write function that adds crp data, leaves zeros if nothing to add
        subj_df = pd.concat(stim_df)

    return subj_df #a df for one subject with all crps for all significant stim-resp pairs as rows


def agg_crp_df(stim_reg: str, ma: str, hd5: str, **kwargs) -> pd.DataFrame: #parallel agg_sesh_df
    """Returns a df with the crp saved in each row,

    Args:
        stim_sesh (str): string of where stim occurred and how much
        ma (str) : mA of stim sesh
        hd5 (str): filename with the full path for opening the hdf5 file associated with this stim

    Returns:
        pd.DataFrame: N_stim_k_resp X M_samples one CRP per row
    """
    sig_test = kwargs['sig_test'] if "sig_test" in kwargs.keys() else False
    stim_df = [] #zero pad here?
    with h5py.File(f, 'r') as f:
        for key in f.keys():
            if sig_test and not get_sig(key, hd5, f[key]):
               continue
            response = f[key]
            crp = response[crp][:] #do timestamp/second conversions, 1 column = 1 stamp, add metadata for sampling rate
            stim_df.append(crp)
        rej = len(f.keys()) - len(stim_df)
        if rej >0 :
            logger.info(f"\t\t\tRejected {rej} for resp: {key}")

    stim_df = pd.concat(stim_df) #indexed by stim_resp pairs
    return stim_df #a df with all crps for one stim and all its resp pairs (i.e., one HD5 file), consult AGG SESH 

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
