# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 11:15:55 2023

@author: Derek
"""

import mne
import datetime
import pandas as pd

ideal_starttime_hours = 22 # 2200
ideal_endtime_hours = 6 # 0600


edf_file = mne.io.read_raw_edf("Y:/SEEG_sleep/Spat36_08292020_07005005.EDF")

info = edf_file.info

edf_file_starttime = info['meas_date'].replace(tzinfo=None)
edf_file_endtime = edf_file_starttime + datetime.timedelta(seconds=edf_file._last_time)

#%%
# all channels from the edf file
edf_ch_names = edf_file.ch_names

import pandas as pd
channels_used = pd.read_excel("Z:/000_Data/SEEG/SEEG_EyesClosed_RestingState/data/Spat36/Spat36_Channels_Used.xlsx", header=None)

bip_names = channels_used[0].to_list()

first_ch, second_ch = list(zip(*(k.split(" - ") for k in bip_names)))

first_ch = list(first_ch)
second_ch = list(second_ch)
all_used_channels = first_ch + second_ch

# remove duplicates
all_used_channels = [*set(all_used_channels)]

ch_to_drop = list(set(edf_ch_names).difference(all_used_channels))

#%% CROP the file

sleep_starttime = datetime.datetime(year=2020, month=8, day=29, hour=22, minute=0, second=0)
sleep_endtime = datetime.datetime(year=2020, month=8, day=30, hour=6, minute=0, second=0)

starttime_in_seconds = (sleep_starttime - edf_file_starttime).seconds
endtime_in_seconds = (sleep_endtime - edf_file_starttime).seconds

edf_file = edf_file.crop(tmin=starttime_in_seconds, tmax=endtime_in_seconds)

#%% remove unneeded channels
edf_file = edf_file.drop_channels(ch_to_drop)

#%%
edf_file.export("Y:/SEEG_sleep/test_v2.edf")