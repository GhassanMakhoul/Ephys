import os
import sys
import h5py
import gc

import unittest

sys.path.append("../")
import resultAggregator as ragg
import numpy as np
#from icecream import ic
#ic.configureOutput(prefix='ic UnitTest| -> ')
import pdb

import pandas as pd

class TestPlotPipeline(unittest.TestCase):

    def setUp(self):
        