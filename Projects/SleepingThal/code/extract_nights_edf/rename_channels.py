# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 06:55:53 2023

@author: Derek
"""

from pyedflib import highlevel

clabel_file = highlevel.read_edf_header("Z:/000_Data/SEEG/SEEG_Sleep_Staging/data/nights_of_sleep/Spat36/Spat36_Night0-0_08252020.edf")
normal_file = highlevel.read_edf_header("Z:/000_Data/SEEG/SEEG_Sleep_Staging/data/nights_of_sleep/Spat36/Spat36_Night2-0_08272020.edf")

clabel_ch = clabel_file["channels"]
correct_ch = normal_file["channels"]

#%%

res = {clabel_ch[i]: correct_ch[i] for i in range(len(clabel_ch))}

#%%

highlevel.rename_channels("Z:/000_Data/SEEG/SEEG_Sleep_Staging/data/nights_of_sleep/Spat36/Spat36_Night0-0_08252020.edf", res, verbose=True)