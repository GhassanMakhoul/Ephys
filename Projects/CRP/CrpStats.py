import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from statistics import NormalDist
import glob
from collections import defaultdict
from scipy.stats import ttest_ind, f_oneway

from utils import *


def merge_label(subj_df: pd.DataFrame, label_df: pd.DataFrame, leftcol: str, rightcol: str ) -> pd.DataFrame:
    """Merges contact label into the subj_df for both stim and response regions

    Args:
        subj_df (pd.DataFrame): dataframe with SPES data, could be CRP, raw SPES, etc
        as long as each row has an entry for a bipole (or monopolar) contact to merge on. 
        ASSUMES 'subj' is a column for sanity checks.
       
        label_df (pd.DataFrame): labels designating 'SOZ', 'NIZ', 'PZ', should be originating 
        from the 'all_pats_bipole.csv', may contain more than one subject, but ASSUMES 'subj' is a column.

        leftcol (str): column of subject bipoles, could be stim, or resp, will be renamed
        'leftcol_label'

       merge_labeleturns:
        pd.DataFrame: subj df with columns for labels
    """
    assert 'subj' in subj_df.columns, "Need a 'subj' column! in subj_df"
    assert 'subj' in label_df.columns, "Needa 'subj' column in label_df"
    assert 'label' in label_df.columns, "Need a label to merge into subj_df, check label_df!"
    subj = subj_df.subj.values[0]
    assert subj in set(label_df.subj.values) ,f"Subject: {subj} missing from label df!"
    assert len(set(subj_df.subj)) == 1, 'Can only merge one subject at a time!'
    og_rows = subj_df.shape[0]
    label_df = label_df[label_df.subj == subj]
    reg_map = defaultdict(lambda: "UNLABELED")
    reg_map.update({reg:label for reg,label in label_df[[rightcol,'label']].values})
    subj_df[ f'{leftcol}_label'] = subj_df[leftcol].apply(lambda x: reg_map[x])  # subj_df = subj_df.merge(label_df[[rightcol, 'label']], left_on=leftcol, right_on=rightcol)
    assert og_rows >= subj_df.shape[0], f"For {subj}\n\t\texpected at most {og_rows}, after merge: {subj_df.shape[0]} rows detected"
    return subj_df


import pdb
def agg_subject_results(result_files: list[str], label_df: pd.DataFrame) -> pd.DataFrame:
    dfs = []
    expected_rows = 0
    for f in result_files:
        res_df = pd.read_csv(f)
        subj = res_df.subj.values[0]
        res_df.resp_reg = res_df.resp_reg.apply(lambda x: x.replace(" ", ""))
        expected_rows += res_df.shape[0]
        tmp = merge_label(res_df,label_df, 'resp_reg', 'bipole')
        merge_df = merge_label(tmp, label_df, 'stim_reg', 'bipole')        
        dfs.append(merge_df)
    dfs = pd.concat(dfs)
    dfs = dfs[dfs.stim_reg_label != "UNLABELED"]
    dfs = dfs[dfs.resp_reg_label != "UNLABELED"]
    assert expected_rows >= dfs.shape[0], f"Expected at most {expected_rows}, got {dfs.shape[0]}\n\t\t diff: {expected_rows - dfs.shape[0]}"
    print(f"{len(set(dfs.subj))} subjects, total of {dfs.shape[0]} trials, dropped: {expected_rows - dfs.shape[0]}")
    return dfs

def score_adj_matrix(adj_mat, by='col'):
    """z-scores a directed adjacency matrix against in-degree or out-degree

    Args:
        mat (np.array): NxN matrix of directed connections
        by (str, optional): normalize by row or column. If specify col then this will center 
        all in degrees (column-wise normalize). Defaults to 'col'.
    NOTE: when choosing a column wise operations, the dim will = 0, this is because
    the columnwise norm collapses rows. numpy.nanmeans asks which dimension to summarize,
    and thus specifying 0 means we wish to collapse along the row dimension
    Returns:
        _type_: _description_
    """
    assert by== 'col' or by =='row', "select proper dimension to summarize, only supports 2D"
    mu, std = get_stats(adj_mat, by)
    centered = np.subtract(adj_mat, mu)
    z_scored = np.divide(centered, std)
    return z_scored

def get_stats(adj_mat, by='col'):
    assert by== 'col' or by =='row', "select proper dimension to summarize, only supports 2D"

    dim = 0 if by == 'col' else 1
    n = adj_mat.shape[0] 

    mu = np.nanmean(adj_mat, axis=dim)
    std = np.nanstd(adj_mat, axis=dim)

    if by == "col": #reshape so that broadcast on subtract does element-wise along row/col
        mu = mu.reshape(1,n)
        std = std.reshape(1,n)
    else:
        mu = mu.reshape(n,1)
        std = std.reshape(n,1)
    return mu,std