import os
import sys
import h5py
import gc
from functools import partial

import unittest
import time

sys.path.append("../")
import crplot as cr
import numpy as np
#from icecream import ic
#ic.configureOutput(prefix='ic UnitTest| -> ')

#Debugging
from loguru import logger
import pdb
import pandas as pd
import warnings
warnings.simplefilter(action='ignore')



class TestPlotPipeline(unittest.TestCase):

    def setUp(self):
        self.plot_file = 'ex_plot.csv'
        plot_df = pd.read_csv('ex_plot.csv')
        self.plot_df = plot_df

        self.bad_f = "ex_plot_error.csv"
        self.bad_df = pd.read_csv('ex_plot_error.csv')

        crossplot_df = plot_df[plot_df.plot_type =='cross']
        self.cross_f = crossplot_df.fname.values[0]
        self.cross_k = crossplot_df.key.values[0]
        self.cross_out_f = crossplot_df.out_fname.values[0]

        df = plot_df[plot_df.plot_type == 'raw']
        self.spes_df = pd.read_csv(df.fname.values[0])
        self.spes_out_f = df.out_fname.values[0]

        reparam_df = plot_df[plot_df.plot_type =='reparam-agg']
        self.reparam_f = reparam_df.fname.values[0]
        self.reparam_key = reparam_df.key.values[0]
        self.reparam_out_f = reparam_df.out_fname.values[0]
        self.startTime = time.time()
    
    def tearDown(self):
        t = time.time() - self.startTime
        logger.info('Test %s: took %.3f s to run' % (self.id(), t))




    def test_Verify(self):
        self.assertRaises(AssertionError,cr.verify_df, self.bad_df) 
         #TODO consider adding more failure conditions systematically to the df
         #If this does not raise an exception, test fails
        cr.verify_df(self.plot_df)
        #If this runs, and no exceptions raised -> code works!
    
    def test_getCrossProject(self):
        
        S, [ma, stim, contact, tr_win] = cr.get_crossproject(self.cross_f,self.cross_k)
        warnings.simplefilter(action='ignore')
        self.assertTrue((S.columns == ['cross_proj', 'win_size']).all())
        self.assertLess(tr_win[1], S.win_size.values[-1])
        self.assertAlmostEqual(tr_win[1], 0.109375) #manually opened

    def test_plotCrossProject(self):
        S, args = cr.get_crossproject(self.cross_f, self.cross_k)
        cr.plot_cross_project(S,self.cross_out_f, *args)
        self.assertTrue(os.path.exists(self.cross_out_f))

    def test_plotGetReparam(self):
        f = self.reparam_f
        k = self.reparam_key
        reparam_df = cr.get_reparam(f,k)
        self.assertIsInstance(reparam_df,pd.DataFrame)
        

    def test_plotReparamAgg(self):
        f = self.reparam_f
        k = self.reparam_key 
        resp = k.split("_")[-1]
        stim = f.split("/")[-2]
        out_f = self.reparam_out_f
        reparam_df = cr.get_reparam(f,k)
        cr.plot_reparam_agg(reparam_df,out_f,resp,stim)
        self.assertTrue(os.path.exists(out_f))
    
    def test_resampleChannels(self):
        chs = cr.gen_plot_channels(self.spes_df.columns, 10)
        self.assertNotIn('trial',chs)
        self.assertEqual(len(chs), 10)
    
    def test_genPlotDf(self):
        subj = 'Epat27'
        h5file = '/mnt/ernie_main/Ghassan/ephys/data/Epat27/LA3-LA4_2mA/stim_resp_bipole.hdf5'
        plot_df = cr.gen_plot_df(subj, h5file)
        with h5py.File(h5file, 'r') as f:
            keys = f.keys()
            num_keys = len(keys)
        plot_df.to_csv("Epat_27_tst_plot.csv",index=False)
        self.assertEqual(plot_df.shape[0],num_keys)
    
    def test_filter(self):
        conditions = {}
        lg_df = pd.read_csv('/mnt/ernie_main/Ghassan/ephys/data/test/Epat26_plots.csv')
        logger.info(f"OG SHAPE { lg_df.shape }")
        conditions['notes'] = {'sig': {'sig:True':10, 'sig:False':5}}
        plot_df = cr.filter_plot_df(lg_df, conditions)

        self.assertEqual(15, plot_df.shape[0])

    def test_plotRaw(self):
        warnings.simplefilter(action='ignore')

        chs = cr.gen_plot_channels(self.spes_df.columns, 5000)
        cr.plot_channels(self.spes_df,chs, self.spes_out_f)
        self.assertTrue(os.path.exists(self.spes_out_f))

    def test_visualizePipe(self):
        cr.visualize_pipeline(self.plot_file)
        files = self.plot_df.out_fname
        for f in files:
           self.assertTrue(os.path.exists(f))


        
        
if __name__ == '__main__':
    logger.add("logs/tst_plot.log", level=20)
    logger.info("Running tests of pipeline")
    unittest.main()

