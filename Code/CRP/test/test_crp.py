import os
import sys
import h5py
import gc

import unittest

sys.path.append("../")
import crp
#from icecream import ic
#ic.configureOutput(prefix='ic UnitTest| -> ')

import pandas as pd
from loguru import logger

class TestResultAgg(unittest.TestCase):

    def setUp(self):

        self.old_subj = {'subj':'Epat27', 'stim_pair':'LA1-LA2', 'mA':'1mA'}
        self.new_subj = {'subj':'Spat52', 'stim_pair':'LI5-LI6', 'mA':'5mA'}

    
    def test_assembleTrial(self):
        #test mat73 loading
        spes_old, fs_old = crp.assemble_trial(self.old_subj['subj'], self.old_subj['stim_pair'], self.old_subj['mA'])
        self.assertIsInstance(spes_old, pd.DataFrame)
        
    

        #test scipy loadmat
        spes_new, fs_new = crp.assemble_trial(self.new_subj['subj'], self.new_subj['stim_pair'], self.new_subj['mA'])
        self.assertIsInstance(spes_new, pd.DataFrame)
        try:
            self.assertEqual(type(fs_new), type(fs_old))
        except AssertionError:
            import pdb
            pdb.set_trace()
            print(fs_old)
            print(fs_new)
        t_old, d_old = spes_old.shape
        t_new, d_new = spes_new.shape
        self.assertAlmostEqual(t_new, t_old, delta=10)
        self.assertAlmostEqual(d_new, d_old, delta=10)

        
        

    def test_crossproject(self):

        spes_old, fs_old = crp.assemble_trial(self.old_subj['subj'], self.old_subj['stim_pair'], self.old_subj['mA'])
        contact =spes_old.columns[0]
        V_trial = spes_old.pivot( columns='trial', values=contact)
        S_old = crp.cross_project_trial(V_trial.values,fs_old)

        #test scipy loadmat
        spes_new, fs_new = crp.assemble_trial(self.new_subj['subj'], self.new_subj['stim_pair'], self.new_subj['mA'])
        contact =spes_new.columns[0]
        V_trial = spes_new.pivot( columns='trial', values=contact)
        S_new = crp.cross_project_trial(V_trial.values,fs_new)
       
        self.assertEqual(type(S_new), type(S_old))


if __name__ == '__main__':
    logger.add("../logs/tst_crp.log", level=20)
    logger.info("Running tests of pipeline")
    unittest.main()