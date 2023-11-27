import os
import sys
import h5py

import unittest

sys.path.append("../")
from resultAggregator import *
#from icecream import ic
#ic.configureOutput(prefix='ic UnitTest| -> ')

import pandas as pd

class TestResultAgg(unittest.TestCase):

    def setUp(self):
        filepath = '/mnt/ernie_main/Ghassan/ephys/data/Epat26/LA1-LA2_1mA/stim_resp.hdf5'
        self.res_folder = '/mnt/ernie_main/Ghassan/ephys/data/'
        self.subj = 'Epat26'
        self.filepath = filepath
        #TODO modularize test file
        # open h5File
        self.h5 = h5py.File(filepath,'r')
        print(f"Loaded h5 {self.h5}\n")
        self.keys = [k for k in self.h5.keys()]


    def tearDown(self):
        self.h5.close()
        print("Closed h5")


    def test_entryToDF(self):
        key1 = self.keys[0]
        resp1 = self.h5[key1]
        print(f'Testing {key1}')
        df = entry_to_df(key1, resp1)
        self.assertIn('TR', df.columns)
        self.assertIn('alphas', df.columns)
        self.assertIn('resp_reg', df.columns)
    
    def test_getStimFolders(self):
        folders = get_stim_folders(self.subj, self.res_folder)
        self.assertNotIn("derivatives", folders)
        self.assertNotIn("figs", folders)
        self.assertEqual(len(folders), 243) #Based on running find . -name '*mA' -type d | wc -l
        #TODO make it test each folder by running this command from python

    def test_getSeshParams(self):
        
        folders = get_stim_folders(self.subj, self.res_folder)
        reg, ma = get_sesh_params(folders[0])
        self.assertRegexpMatches(reg, '[A-Z]+[0-9]+\-[A-Z]+[0-9]+')
        self.assertRegexpMatches(ma, '[0-9]+mA')
    

    def test_aggSesh(self):
        df = agg_sesh_df(self.h5)
        self.assertEqual(len(self.keys), len(set(df['resp_reg'].values)))
    
    def test_aggResponses(self):
        """run the aggregation once as a sanity check, should surface any errors
        """
        folders = get_stim_folders(self.subj, self.res_folder)
        agg_responses(self.subj, self.filepath, \
                      folders, '/mnt/ernie_main/Ghassan/ephys/test/')
if __name__ == '__main__':
    unittest.main()