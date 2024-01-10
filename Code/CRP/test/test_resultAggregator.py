import os
import sys
import h5py
import gc
import yaml
import unittest
import time

sys.path.append("../")
import resultAggregator as ragg
import numpy as np
import pandas as pd
#from icecream import ic
#ic.configureOutput(prefix='ic UnitTest| -> ')
#Debugging
import pdb
from loguru import logger



class TestResultAgg(unittest.TestCase):

    def setUp(self):
        self.h5fname = 'stim_resp_bipole.hdf5'
        filepath = '/mnt/ernie_main/Ghassan/ephys/data/Epat26/LA1-LA2_1mA/stim_resp_bipole.hdf5'
        self.res_folder = '/mnt/ernie_main/Ghassan/ephys/data/'
        self.subj = 'Epat26'
        self.filepath = filepath
        self.pathout = '/mnt/ernie_main/Ghassan/ephys/data/test/'
        #TODO modularize test file
        # open h5File
        self.h5 = h5py.File(filepath,'r')
        print(f"Loaded h5 {self.h5}\n")
        self.keys = [k for k in self.h5.keys()]
        self.startTime = time.time()

    
    def tearDown(self):
        t = time.time() - self.startTime
        logger.info("Test %s: took %.3f s to run" % (self.id(), t))
        self.h5.close()
        print("Closed h5")
        print("Garbage Collecting")
        gc.collect()

    def test_getSig(self):
        key = 'response_LAM9 - LAM10'
        sig = ragg.get_sig(key, self.filepath, self.h5[key])
        self.assertIsInstance(sig, np.bool_) #TODO make stronger test

    def test_entryToDF(self):
        key1 = self.keys[0]
        resp1 = self.h5[key1]
        print(f'Testing {key1}')
        df = ragg.entry_to_df(key1, resp1)
        self.assertIn('TR', df.columns)
        self.assertIn('alphas', df.columns)
        self.assertIn('resp_reg', df.columns)
        self.assertIn('explained_variance', df.columns)
    
    def test_getStimFolders(self):
        folders = ragg.get_stim_folders(self.subj, self.res_folder)
        self.assertNotIn("derivatives", folders)
        self.assertNotIn("figs", folders)
        self.assertEqual(len(folders), 243) #Based on running find . -name '*mA' -type d | wc -l
        #TODO make it test each folder by running this command from python

    def test_getSeshParams(self):
        
        folders = ragg.get_stim_folders(self.subj, self.res_folder)
        reg, ma = ragg.get_sesh_params(folders[0])
        self.assertRegex(reg, '[A-Z]+[0-9]+\-[A-Z]+[0-9]+', "Stim region not found!")
        self.assertRegex(ma, '[0-9]+mA', "No mA in sesh params!")


    def test_explainedVariance(self):
        k = self.keys[0]
        pulse_trial = self.h5[k]
        V_tr = pulse_trial['V_tr'][:]
        num_pulses = min(V_tr.shape)
        ev = ragg.get_explained_var(pulse_trial)
        self.assertEqual(len(ev.shape), 1)
        self.assertEqual(ev.shape[0],num_pulses)
        self.assertTrue(np.all(ev < 1)) # explained var needs to be less than 1
        self.assertFalse(np.all(np.isnan(ev))) # no NaN values allowed!
    
    def test_genplot(self):
        with open("../config_agg.yml", 'r') as f:
            config =  yaml.safe_load(f)
        config['plot']['plot_path'] =  config['plot']['plot_path'].replace("SUBJ", self.subj)
        stim_folders = ragg.get_stim_folders(self.subj, self.res_folder)
        plot_kwargs = config['plot']
        ragg.gen_plot_file(self.subj,self.h5fname,stim_folders,self.pathout,**plot_kwargs )
        plot_file = os.path.join(self.pathout,self.subj+ "_plots.csv")
        self.assertTrue(os.path.exists(plot_file))
        # run gen plot file
        # check for notes ? entries ?
    
    def test_aggSesh(self):
        print("Test agg sesh")
        df = ragg.agg_sesh_df(self.filepath)
        self.assertTrue(len(self.keys) >= len(set(df['resp_reg'].values)))
        
    
    
    # def test_aggResponses(self):
    #     """run the aggregation once as a sanity check, should surface any errors
    #     """
    #     print("Test Agg Resp")
    #     folders = ragg.get_stim_folders(self.subj, self.res_folder)
    #     ragg.agg_responses(self.subj, self.filepath, \
    #                   folders, '/mnt/ernie_main/Ghassan/ephys/test/')
    #     df = pd.read_csv("/mnt/ernie_main/Ghassan/ephys/test/Epat26_stim.csv")
    #     self.assertEqual(len(set(df.subj)), 1)
        
    # def test_main(self):
    #    print("testing main method")
    #    config = 'config_tst_agg.yml'
    #    pathout = '/mnt/ernie_main/Ghassan/ephys/data/test/'
    #    ragg.main(['-s', 'Epat26', '-p', pathout, '-c', config])
    #    res_file = os.path.join(pathout,self.subj + "_stim.csv")
    #    self.assertTrue(os.path.exists(res_file))


if __name__ == '__main__':
    logger.add('logs/tst_res_agg.log', level=20)
    logger.info("Running Result aggregator")
    unittest.main()