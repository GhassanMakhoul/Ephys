#!/usr/bin/env python3

import logging
import os
logging.getLogger('mat73').setLevel(logging.CRITICAL)

import pandas as pd

from utils import *
from connectivity_dynamics import *


if __name__ == '__main__':

    ## get FPAC regions and store in lut_data
    fpac_lut = pd.read_csv('FPAC_lut.txt')
    lut_data = fpac_lut.region.to_numpy()

    DATA_DIR = '/mnt/ernie_main/Price/ephys/data/periconnectivity'

    # load multiple patients
    agg_df_lst = []
    stats_dict = {}
    stats_keys = ['patID','eventID','FPAC_bips','other_bips']
    files = glob.glob(os.path.join(DATA_DIR,"*flow*pat*.csv"))
    files = [f for f in files if ('verbose' in f)]

    for f in files[0:2]:
        flow_df = pd.read_csv(os.path.join(DATA_DIR,f))
        flow_df = flow_df[flow_df.freq_band == 'alpha']
        flow_df = flow_df[flow_df.source == 'nz']

        # get patient's bip to region dictionary
        patID = flow_df.patID.iloc[0]
        mat_paths = glob.glob(f'/mnt/ernie_main/Ghassan/ephys/data/connectivity/{patID}/{patID}*PDC.mat')
        obj = load_mat(mat_paths[0])
        regions = get_atlas_regions(obj)
        chans = get_chan_names(obj)
        bip_to_region_dict = dict(zip(chans, regions))

        flow_df['region'] = flow_df.parallel_apply(lambda x: bip_to_region_dict[x['src_bip']], axis=1)
        flow_df['FPAC'] = flow_df.parallel_apply(lambda x: (x.region in lut_data), axis=1)
        flow_df = center_onset(flow_df, **{'mid_sz_length':20})
        flow_df.dropna()

        # collect stats
        for event in flow_df.eventID.unique():
            stat_df = flow_df.copy()
            stat_df = stat_df[(flow_df.period == 0) & (flow_df.eventID == event) & (flow_df.target == 'nz')]
            n_fpac = stat_df[stat_df.FPAC].shape[0]
            n_other = stat_df[~stat_df.FPAC].shape[0]
            for key, val in zip(stats_keys,(patID,event,n_fpac,n_other)):
                stats_dict.setdefault(key,[]).append(val)
        del stat_df

        agg_df = agg_verbose_df(
            flow_df,
            measure_cols = ['value'],
            categorical_cols = ['target','win_sz_st_end','eventID','patID','sz_type','FPAC']
        )

        agg_df = agg_verbose_df(
            agg_df,
            measure_cols = ['value'],
            categorical_cols = ['patID','sz_type','target','win_sz_st_end','FPAC']
        )

        del flow_df
        agg_df_lst.append(agg_df)

    grp_flow_df = pd.concat(agg_df_lst)
    grp_flow_df['src_trgt'] = grp_flow_df.parallel_apply(lambda x : "nz_"+x['target'], axis=1)
    grp_flow_df.to_csv(f'{DATA_DIR}/verbose_grp_flow.csv', index=False)

    stats_df = pd.DataFrame.from_dict(stats_dict)
    stats_df.to_csv(f'{DATA_DIR}/stats_verbose_grp_flow.csv', index=False)