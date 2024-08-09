#!/usr/bin/env python3

import argparse
from loguru import logger
import logging
import os
import sys
logging.getLogger('mat73').setLevel(logging.CRITICAL)

import pandas as pd

from utils import *
from connectivity_dynamics import *


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='create summary of flow from/to a specified source/target'
                                                 '. user can only specify source or target. defaults to '
                                                 'target = soz')
    parser.add_argument("-s", "--source", choices=['soz','pz','nz'], default='', help="source of flow")
    parser.add_argument("-t", "--target", choices=['soz','pz','nz'], default='', help="target of flow")
    args = parser.parse_args()

    if (len(args.source) > 0) & (len(args.target) > 0):
        # user cannot specify both source and target
        logger.error('Must specify either source or target, not both')
        sys.exit()
    elif (len(args.source) == 0) & (len(args.target) == 0):
        # if user doesn't specify source or target, defaults to target=soz
        args.target = 'soz'

    if args.source:
        retain_col = 'target'
        logger.info(f'Source set to {args.source}')
    else:
        retain_col = 'source'
        logger.info(f'Target set to {args.target}')

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

    for f in files:
        flow_df = pd.read_csv(os.path.join(DATA_DIR,f))
        flow_df = flow_df[flow_df.freq_band == 'alpha']
        if args.source:
            flow_df = flow_df[flow_df.source == args.source]
        else:
            flow_df = flow_df[flow_df.target == args.target]

        # get patient's bip to region dictionary
        patID = flow_df.patID.iloc[0]
        logger.info(f'Running patient {patID}')
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
            logger.info(f'+Working on event {event}+')
            stat_df = flow_df.copy()
            stat_df = stat_df[(stat_df.period == 0) & (stat_df.eventID == event)]
            if args.source:
                stat_df = stat_df[stat_df.target == 'nz']
            else:
                stat_df = stat_df[stat_df.source == 'nz']

            n_fpac = stat_df[stat_df.FPAC].shape[0]
            n_other = stat_df[~stat_df.FPAC].shape[0]
            for key, val in zip(stats_keys,(patID,event,n_fpac,n_other)):
                stats_dict.setdefault(key,[]).append(val)
        del stat_df

        agg_df = agg_verbose_df(
            flow_df,
            measure_cols = ['value'],
            categorical_cols = [retain_col,'win_sz_st_end','eventID','patID','sz_type','FPAC']
        )

        agg_df = agg_verbose_df(
            agg_df,
            measure_cols = ['value'],
            categorical_cols = ['patID','sz_type',retain_col,'win_sz_st_end','FPAC']
        )

        del flow_df
        agg_df_lst.append(agg_df)

    grp_flow_df = pd.concat(agg_df_lst)

    if args.source:
        grp_flow_df['src_trgt'] = grp_flow_df.parallel_apply(lambda x : f"{args.source}_"+x['target'], axis=1)
        file_suffix = f'source_{args.source}'
    else:
        grp_flow_df['src_trgt'] = grp_flow_df.parallel_apply(lambda x : x['source']+f"_{args.target}", axis=1)
        file_suffix = f'target_{args.target}'

    grp_flow_df.to_csv(f'{DATA_DIR}/verbose_grp_flow_{file_suffix}.csv', index=False)

    logger.success(f'Created summary file at {DATA_DIR}/verbose_grp_flow_{file_suffix}.csv')

    stats_df = pd.DataFrame.from_dict(stats_dict)
    stats_df.to_csv(f'{DATA_DIR}/stats_verbose_grp_flow_{file_suffix}.csv', index=False)
    logger.success(f'Created stats file at {DATA_DIR}/stats_verbose_grp_flow_{file_suffix}.csv')