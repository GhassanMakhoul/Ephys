import os
import sys
import h5py
import gc

import unittest

sys.path.append("../")
import crplot as cr
import numpy as np
#from icecream import ic
#ic.configureOutput(prefix='ic UnitTest| -> ')
import pdb

import pandas as pd

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



    def test_Verify(self):
        with self.assertRaises(AssertionError) as ae:
            cr.verify_df(self.bad_df)#TODO consider adding more failure conditions systematically to the df
            #If this does not raise an exception, test fails
        cr.verify_df(self.plot_df)
        #If this runs, and no exceptions raised -> code works!
    
    def test_getCrossProject(self):
        
        S, [ma, stim, contact, tr_win] = cr.get_crossproject(self.cross_f,self.cross_k)

        self.assertAlmostEqual(S.columns, ['cross_proj', 'win_size'])
        self.assertLess(tr_win[1], S.win_size.values[-1])
        self.assertAlmostEqual(tr_win[1], .527)

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
        out_f = self.reparam_out_f
        reparam_df = cr.get_reparam(f,k)
        cr.plot_reparam_agg(reparam_df,out_f,proc="BIPOLE")
        self.assertTrue(os.path.exists(out_f))
    
    def test_resampleChannels(self):
        chs = cr.gen_plot_channels(self.spes_df.columns, 10)
        self.assertNotIn('trial',chs)
        self.assertEqual(len(chs), 10)

        
        chs = cr.gen_plot_channels(self.spes_df.columns, 5000)
        self.assertNotIn('trial',chs)
        self.assertEqual(len(chs), len(self.spes_df.columns)-1)

    def test_plotRaw(self):
        cr.plot_channels(self.spes_df,chs, self.spes_out_f)
        self.assertTrue(os.path.exists(self.spes_out_f))

    def test_visualizePipe(self):
        cr.visualize_pipeline(self.plot_file)
        files = self.plot_df.out_fname
        for f in files:
            self.assertTrue(os.path.exists(f))


        
        
if __name__ == '__main__':
    logger.add("../logs/tst_crp.log", level=20)
    logger.info("Running tests of pipeline")
    unittest.main()

