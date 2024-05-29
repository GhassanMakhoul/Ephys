import os
import re
from scipy.io import loadmat
import mat73

from collections import Counter
import pandas as pd
import numpy as np
import mne

import seaborn as sns
import matplotlib.pyplot as plt

from utils import *

TIMELINE_F = '/mnt/ernie_main/000_Data/SEEG/SEEG_Periictal/data/Extracted_Per_Event_Interictal/all_time_data_01042023_212306.csv'
SEEG_FOLDER = '/mnt/ernie_main/000_Data/SEEG/SEEG_Entire_EMU_Downloads/data/'
BANDS = ['delta', 'theta', 'alpha', 'beta','gamma_l', 'gamma_H']
PERIOD = ['inter','pre','ictal','post']

def load_mat(f):
    try:
        return loadmat(f)
    except NotImplementedError:
        return mat73.loadmat(f)

