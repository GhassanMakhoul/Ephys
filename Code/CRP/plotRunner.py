#I/O  Setup
import os
import sys
import getopt
import re
import h5py
import glob
from loguru import logger
import yaml


#data and math packages
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.stats import ttest_1samp

#specialty
import crplot


def main(argv):
    subj = ''
    config_f = 'config_plot.yml'
    filt = {}
    opts, _ = getopt.getopt(argv,"i:c:",["subj=",'config='])
    for opt, arg in opts:
        if opt in ("-i", 'inpf'):
            inpf = arg
        elif opt in ("-c", '--config'):
            config_f = arg
    with open(config_f, 'r') as f:
        config =  yaml.safe_load(f)
        filt = config['filter']
        plot_opts = config['plot_options']
    crplot.visualize_pipeline(inpf, filt, **plot_opts)

if __name__ == '__main__':
    main(sys.argv[1:])
