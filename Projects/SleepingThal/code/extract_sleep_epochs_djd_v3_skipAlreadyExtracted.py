# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 13:57:39 2022

@author: dossd
"""

#%% Imports
import pandas as pd
import datetime
import numpy as np
import glob
import math
import os
import mne
from scipy.io import savemat

from pyedflib import highlevel

from tqdm import tqdm

import multiprocessing

from joblib import Parallel, delayed

import traceback

data_from_root = "Z:/000_Data/SEEG/SEEG_Entire_EMU_Downloads/data/"

# save root
save_root = "Z:/000_Data/SEEG/SEEG_Periictal/data/Extracted_Per_Event_Interictal/"

csv_file_path = "Z:/000_Data/SEEG/SEEG_Periictal/data/Extracted_Per_Event_Interictal/all_time_data_01092023_112957.csv"

ictal_sheet_fpath = "Z:/000_Data/SEEG/SEEG_Periictal/notes/Master_Ictal_Event_Timestamps_SeizureType_SOZ_PZ_batch3_preprocessing_pats.xlsx"

patient_batch_sheet = "Z:/000_Data/SEEG/SEEG_Periictal/notes/all_preprocessing_pats.xlsx"

all_time_df = pd.read_csv(csv_file_path,sep='\t')

sleep_epoch_dir = "Z:/000_Data/SEEG/SEEG_Sleep_Staging/data/extracted_sleep_epochs_v2/"
nights_of_sleep_dir = "Z:/000_Data/SEEG/SEEG_Sleep_Staging/data/nights_of_sleep/"

# sleep categories that I definitely want. N1 is a bit buggy, so it's not completely necessary
confident_sleep_categories = ["R", "W", "N2", "N3"]

# per seizure interictal duration
interictal_duration = 5 # minutes
interictal_timedelta = datetime.timedelta(minutes=interictal_duration)

# duration of windows that will be calculated
window_duration = 5 # seconds
window_timedelta = datetime.timedelta(seconds=window_duration)

# buffer times between seizures for interictal determination (minutes)
ictal_buffer_time = 60
ictal_buffer_time_hours = int(ictal_buffer_time/60)

# optimum time before sizure onset
optimum_time_before = 4.00 + (interictal_duration/60) #hours


pat_list_excel = pd.read_excel(patient_batch_sheet)

patID_list = pat_list_excel['patIDs'].tolist()
pat_names = patID_list

pat_nights_of_sleep_dirs = glob.glob(nights_of_sleep_dir+"*pat*")
pat_names = list()

already_extracted_pat_dirs = glob.glob(sleep_epoch_dir + "*pat*")
pats_to_rem = list()
for i in range(len(already_extracted_pat_dirs)):
    # get all the monopole files in the directory
    tmp_pat_epoch_dir = already_extracted_pat_dirs[i].replace("\\", "/")
    tmp_file_types = confident_sleep_categories.copy()
    for tmp_monopole_file in glob.iglob(tmp_pat_epoch_dir + "/*_monopole.mat"):
        tmp_monopole_file = tmp_monopole_file.replace("\\","/").split("/")[-1]
        tmp_sleep_type = tmp_monopole_file.split("_")[1]
        # could throw an error if N1 is found and it isn't in the list
        if tmp_sleep_type != "N1":
            tmp_file_types.remove(tmp_sleep_type)
    
    if len(tmp_file_types) == 0:
        pats_to_rem.append(already_extracted_pat_dirs[i].replace("\\", "/").split("/")[-1])



for i in range(len(pat_nights_of_sleep_dirs)):
    tmp_pat_name = pat_nights_of_sleep_dirs[i].replace("\\", "/").split("/")[-1]
    if tmp_pat_name not in pats_to_rem:
        pat_names.append(tmp_pat_name)
        
        

clean_sleep_df = pd.read_csv("Z:/000_Data/SEEG/SEEG_Sleep_Staging/notes/cleaned_sleep_times_06122023.csv",
                             delimiter="\t")

sleep_categories = ["R", "W", "N1", "N2", "N3"]
sleep_epoch_duration_s = 300 # 5 minutes

#%%

def append_error_to_master_error(fname):
    with open(save_root+"error.txt",'a') as f:
        f.write(fname+"\n")
        
def append_custom_error_to_master_error(fname, custom_err_message):
    with open(save_root+"error.txt",'a') as f:
        f.write(fname + custom_err_message + "\n")

def change_str_to_datetime(df, cols_to_loop_through):
    for j in range(len(cols_to_loop_through)):
    
        tmp_list = list()
        for i in range(df.shape[0]):
    
            if df.iloc[i][cols_to_loop_through[j]] == 'None' or df.iloc[i][cols_to_loop_through[j]] == 'Unknown'or pd.isna(df.iloc[i][cols_to_loop_through[j]]):
                tmp_list.append(df.iloc[i][cols_to_loop_through[j]])
                continue
            
            split_datestr = df.iloc[i][cols_to_loop_through[j]].split(".")
            if len(split_datestr) == 2:
                clean_datestr = split_datestr[0] + ":" + split_datestr[1].ljust(6,"0")
                tmp_start_datetime =  datetime.datetime.strptime(clean_datestr, "%Y-%m-%d %H:%M:%S:%f")
            else:
                clean_datestr = split_datestr[0]
                tmp_start_datetime =  datetime.datetime.strptime(clean_datestr, "%Y-%m-%d %H:%M:%S")
            tmp_list.append(tmp_start_datetime)
        df[cols_to_loop_through[j]] = tmp_list
    return df
        
def convert_csv_datetime_to_python_datetime(datetime_string):
    split_datestr = datetime_string.split(".")
    if len(split_datestr) == 2:
        clean_datestr = split_datestr[0] + ":" + split_datestr[1].ljust(6,"0")
        tmp_start_datetime =  datetime.datetime.strptime(clean_datestr, "%Y-%m-%d %H:%M:%S:%f")
    else:
        clean_datestr = split_datestr[0]
        tmp_start_datetime =  datetime.datetime.strptime(clean_datestr, "%Y-%m-%d %H:%M:%S")
    return tmp_start_datetime

def fix_clabel_file(file_path, file_list):
        
    edf_file = mne.io.read_raw_edf(file_path, verbose=40)
    info = edf_file.info
    file_channels = info['ch_names']
    
    i = file_list.index(file_path)

    # find which files do not have the clabel problem
    all_clabel_idx = [idx for idx, s in enumerate(file_list) if 'CLABEL' in s]
    all_file_idx = list(range(0,len(file_list)))
    not_clabel_idx = list(set(all_file_idx).symmetric_difference(set(all_clabel_idx)))
    not_clabel_idx.sort()
    not_clabel_idx = np.asarray(not_clabel_idx)
    
    # find the normal file that is closest in time to the clabel file
    ch_num_matches = False
    counter = 1
    try:
        while not(ch_num_matches):    
            print("Iteration " + str(counter) + " trying to fix C-Label File")
            counter = counter + 1
            closest_file_idx = np.argmin(not_clabel_idx - i)
            edf_normal_fname = file_list[not_clabel_idx[closest_file_idx]].replace("\\","/")
            edf_file_tmp = mne.io.read_raw_edf(edf_normal_fname, verbose=40)
            # if the num channels is different between the saved value and the num channels in the file
            if len(file_channels) != len(edf_file_tmp.info['ch_names']):
                np.delete(not_clabel_idx,(not_clabel_idx==1).nonzero()[0][0])
            else:
                ch_num_matches = True
            pat_ch_names = edf_file_tmp.info['ch_names']
            del edf_file_tmp
    except:
        print("problem finding labels for a c-label file")
    print("Clabel file fixed!")
    channel_names = pat_ch_names
    
    return edf_file, channel_names, info

def get_physical_and_digital_min_max(edf_fpath):
    edf_hdr = highlevel.read_edf_header(edf_fpath, read_annotations=False)
    
    physical_min = list()
    physical_max = list()
    digital_min = list()
    digital_max = list()
    
    for mm in range(len(edf_hdr['SignalHeaders'])):
        physical_min.append(edf_hdr['SignalHeaders'][mm]['physical_min'])
        physical_max.append(edf_hdr['SignalHeaders'][mm]['physical_max'])
        digital_min.append(edf_hdr['SignalHeaders'][mm]['digital_min'])
        digital_max.append(edf_hdr['SignalHeaders'][mm]['digital_max'])
    return physical_min, physical_max, digital_min, digital_max

def write_mat_files(start_datetime, stop_datetime, edf_file, info, channel_names, dst_path, save_fname, physical_min, physical_max, digital_min, digital_max):
    edf_file_starttime = info['meas_date'].replace(tzinfo=None)
    onset_time_diff_sec = (start_datetime - edf_file_starttime).seconds
    offset_time_diff_sec = (stop_datetime - edf_file_starttime).seconds
    data_path = dst_path + save_fname
    
    start_datetime_str = start_datetime.strftime("%Y-%m-%d %H:%M:%S")
    stop_datetime_str = stop_datetime.strftime("%Y-%m-%d %H:%M:%S")
    
    try:
        start_index = edf_file.time_as_index(onset_time_diff_sec)[0]
        stop_index = edf_file.time_as_index(offset_time_diff_sec)[0]
        
        # Check if the directory exists and if not create it  
        # if not os.path.exists(dst_path):
        #     os.mkdir(dst_path)
        
        # get the seizure data and save it
        temp_data = edf_file.get_data(start=start_index, stop=stop_index)
        num_zeros = np.sum(np.around(temp_data,6)==0)
        if (num_zeros / info['sfreq'] / temp_data.shape[0]) > 120:
            return 0
        else:
            save_dic = {"monopole":temp_data, "channels":channel_names, "sfreq":info['sfreq'], "digital_min":digital_min, "digital_max":digital_max,
                        "physical_min":physical_min, "physical_max":physical_max, "onset_datetime":start_datetime_str, "offset_datetime":stop_datetime_str}
            savemat(data_path, save_dic)
            return 1
    except:
        traceback.print_exc()
        print("problem when getting data and writing mat files")
        error_data_path = data_path.replace(".mat","_ERROR_writing_mat.txt")
        # with open(error_data_path, 'w') as fp:
        #     pass
        append_error_to_master_error(error_data_path)

def write_edf_data_from_times(pat_id, start_datetime, stop_datetime, dst_folder, save_fname, event_num, sz_type): 
    data_from_root = "Z:/000_Data/SEEG/SEEG_Entire_EMU_Downloads/data/"
    csv_list = glob.glob(data_from_root+pat_id+"/"+pat_id+"*.csv")
    # edf_list = csv_list.upper().replace("\\","/").replace(".CSV",".EDF")
    edf_list = [s.upper().replace("\\","/").replace(".CSV",".EDF") for s in csv_list]
    edf_fpath = ""
    for csv_fpath in csv_list:
        tmp_df = pd.read_csv(csv_fpath)
        tmp_edf_start = convert_csv_datetime_to_python_datetime(tmp_df.iloc[0]['startdate'])
        tmp_edf_stop = tmp_edf_start + datetime.timedelta(seconds=int(tmp_df.iloc[0]['Duration (sec)']))
        if (tmp_edf_start < start_datetime) and (tmp_edf_stop > stop_datetime):
            edf_fpath = csv_fpath.upper().replace(".CSV",".EDF")
            break
        
    if edf_fpath == "":
        if not os.path.exists(dst_folder):
            os.mkdir(dst_folder)
        save_fname = pat_id+"_"+event_num+"_"+sz_type+"_interictal_monopole_ERROR_no_file_for_seizure.txt"
        with open(dst_folder + save_fname, 'w') as fp:
            pass
        append_error_to_master_error(dst_folder + save_fname)
        return "None", 0
    
    edf_fpath = edf_fpath.upper().replace("\\","/").replace(".CSV",".EDF")
    
    if "CLABEL" in edf_fpath.replace("\\","/").split("/")[-1].split(".")[0].split("_")[-1]:
        print("Found C-Label File")
        edf_file, channel_names, info = fix_clabel_file(edf_fpath, edf_list)
    else:
        edf_file = mne.io.read_raw_edf(edf_fpath, verbose=40)
        info = edf_file.info
        channel_names = info['ch_names']
        
    physical_min, physical_max, digital_min, digital_max = get_physical_and_digital_min_max(edf_fpath)
    
    success_code = write_mat_files(start_datetime, stop_datetime, edf_file, info, channel_names, dst_folder, save_fname, physical_min, physical_max, digital_min, digital_max)
    
    return edf_file, success_code
    

def trim_zeros(periictal_data, fs):
    # inputs: 
    #   peri-ictal data as np
    #   pre_post_flag -> pre = 1, post = 2
    #   fs -> sampling frequency
    # outputs: 
    #   cleaned data
    #   success code for if we need to skip
    
    # duration of the window
    zero_window_dur_sec = 5
    zero_window_dur_samples = zero_window_dur_sec * fs
    
    zero_windows = range(0, np.shape(periictal_data)[1], int(zero_window_dur_samples))
    num_zero_windows = len(zero_windows)
    
    zero_window_flag = np.zeros((num_zero_windows,1))

    for ii in range(len(zero_windows)):
        if ii == (len(zero_windows) - 1):
            tmp_windowed_data = periictal_data[:,zero_windows[ii]:]
        else:
            tmp_windowed_data = periictal_data[:,zero_windows[ii]:zero_windows[ii+1]]
        tmp_summed_zeros = np.sum(np.around(tmp_windowed_data,decimals=6)==0.000000, axis=1) / zero_window_dur_samples
        zero_window_flag[ii] = np.sum(tmp_summed_zeros >= 1.0) > 0
    
    if np.sum(zero_window_flag) == 0:
        success_code = 1
        return success_code, periictal_data
    elif np.sum(zero_window_flag) > 0:
        print("There are " + str(np.sum(zero_window_flag)) + " total zeros :(" )
        return 0, periictal_data
    else:
        print("An error has occured in trim zeros")

def remove_windows_with_zeros(sorted_indexes, startTimes, interictal_timedelta, window_size, fname_list, pat_id):
    # inputs: 
    #   sorted_indexs
    #   startTimes
    #   interictal_timedelta -> how long the total window is
    #   window_size -> window size used to check for zeros
    #   fname_list -> list of fnames that correspond to the starttimes
    #   data_from_root
    #   pat_id
    # outputs:
    #   success_code
    #   best_index -> index (from sorted indexes) that is the best
    
    data_from_root = "Z:/000_Data/SEEG/SEEG_Entire_EMU_Downloads/data/"
    
    file_idx_without_zeros = -1
    
    # find which edf file needs to be read in
    past_fname = ""
    counter = 1
    # make it succeed twice in a row
    historical_success_code = 0
    for ii in sorted_indexes:
        print("Checking proposed interictal " + str(counter) + "/" + str(len(sorted_indexes)))
        counter = counter + 1
        curr_fname = fname_list[ii]
        if curr_fname != past_fname:
            # get rid of the old one if we're loading in a new one
            if 'edf_file' in locals():
                edf_file = None
            # load in the edf file
            past_fname = curr_fname
            edf_fpath = data_from_root + pat_id + "/" + curr_fname
            edf_file = mne.io.read_raw_edf(edf_fpath, verbose=40)
            info = edf_file.info
            edf_file_starttime = info['meas_date'].replace(tzinfo=None)
            fs = info['sfreq']
            
        onset_datetime = startTimes[ii]
        offset_datetime = startTimes[ii] + interictal_timedelta
        
        # get the index of the data to grab
        s_diff_start = (onset_datetime - edf_file_starttime).seconds
        ms_diff_start = (onset_datetime - edf_file_starttime).microseconds / 1000
        tot_s_diff_start = s_diff_start + ms_diff_start / 1000
        start_index = math.floor(tot_s_diff_start * fs)
        
        s_diff_stop = (offset_datetime - edf_file_starttime).seconds
        ms_diff_stop = (offset_datetime - edf_file_starttime).microseconds / 1000
        tot_s_diff_stop = s_diff_stop + ms_diff_stop / 1000
        stop_index = math.ceil(tot_s_diff_stop * fs)
        
        temp_data = edf_file.get_data(start=start_index, stop=stop_index)
        
        success_code, interictal_data = trim_zeros(temp_data, fs)
        
        # if there are zeros in the file, keep searching
        # if there aren't zeros in the file, we're done!
        if success_code == 1 and historical_success_code == 1:
            file_idx_without_zeros = ii
            break
        elif success_code == 1 and historical_success_code == 0:
            historical_success_code = 1
        elif success_code == 0 and historical_success_code == 1:
            historical_success_code = 0
            
        
    # check in case no clean file was found
    if file_idx_without_zeros == -1:
        return 0, file_idx_without_zeros
    else:
        return 1, file_idx_without_zeros
        

def write_edf_data_from_times_sleep_epochs(pat_id, start_datetime, stop_datetime, dst_folder, save_fname):
    data_from_root = "Z:/000_Data/SEEG/SEEG_Entire_EMU_Downloads/data/"
    csv_list = glob.glob(data_from_root+pat_id+"/"+pat_id+"*.csv")
    
    edf_list = [s.upper().replace("\\","/").replace(".CSV",".EDF") for s in csv_list]
    edf_fpath = ""
    for csv_fpath in csv_list:
        tmp_df = pd.read_csv(csv_fpath)
        tmp_edf_start = convert_csv_datetime_to_python_datetime(tmp_df.iloc[0]['startdate'])
        tmp_edf_stop = tmp_edf_start + datetime.timedelta(seconds=int(tmp_df.iloc[0]['Duration (sec)']))
        if (tmp_edf_start < start_datetime) and (tmp_edf_stop > stop_datetime):
            edf_fpath = csv_fpath.upper().replace(".CSV",".EDF")
            break
        
    if edf_fpath == "":
        if not os.path.exists(dst_folder):
            os.mkdir(dst_folder)
        save_fname = pat_id+"_"+"_interictal_monopole_ERROR_no_file_for_seizure.txt"
        with open(dst_folder + save_fname, 'w') as fp:
            pass
        append_error_to_master_error(dst_folder + save_fname)
        return "None", 0
    
    edf_fpath = edf_fpath.upper().replace("\\","/").replace(".CSV",".EDF")
    
    if "CLABEL" in edf_fpath.replace("\\","/").split("/")[-1].split(".")[0].split("_")[-1]:
        print("Found C-Label File")
        edf_file, channel_names, info = fix_clabel_file(edf_fpath, edf_list)
    else:
        edf_file = mne.io.read_raw_edf(edf_fpath, verbose=40)
        info = edf_file.info
        channel_names = info['ch_names']
        
    physical_min, physical_max, digital_min, digital_max = get_physical_and_digital_min_max(edf_fpath)
    
    success_code = write_mat_files(start_datetime, stop_datetime, edf_file, info, channel_names, dst_folder, save_fname, physical_min, physical_max, digital_min, digital_max)
    
    return edf_file, success_code



#%%

for n in range(len(pat_names)):
    # try:
    pat_id = pat_names[n]
    print("\npatient " + pat_id + ": " + str(n+1) + " out of " + str(len(pat_names)))

    ch_list_csv_fpath = glob.glob(data_from_root+pat_id+"/all_EDF_channel_list*.csv")
    ch_list_csv = pd.read_csv(ch_list_csv_fpath[0])
    # skip this patient if there are problems with the channels
    if ch_list_csv.isnull().values.any():
        print("Skipping " + pat_id + " because the channels have a problem")
        dst_err_folder = save_root + pat_id + "_ERROR_ch_names/"
        if not os.path.exists(dst_err_folder):
            os.mkdir(dst_err_folder)
        continue


    pat_time_df = all_time_df[all_time_df['Pat ID'] == pat_id]
    
    pat_time_df = pat_time_df.sort_values(by=['onset_datetime'])
    
    pat_file_df = pat_time_df[pat_time_df['Type'] == 'File']
    pat_file_df = pat_file_df.reset_index()
    
    pat_time_df = pat_time_df[(pat_time_df['Type'] == 'Seizure') | (pat_time_df['Type'] == 'Stimulation')]
    pat_time_df = pat_time_df.reset_index()
    

    
    idx_to_rem = pat_time_df[(pat_time_df['onset_datetime_buffer'] == 'None') | (pat_time_df['onset_datetime_buffer'] == 'Unknown') | (pat_time_df['onset_datetime_buffer']).isnull() ].index
    if idx_to_rem.shape[0] > 0:
        sorted_idx_to_rem = np.sort(np.asarray(idx_to_rem))[::-1]
        for k in range(len(sorted_idx_to_rem)):
            event_num = str(pat_time_df.iloc[sorted_idx_to_rem[k]]['Event Number (Matches EMU Final Report)'])
            sz_type = pat_time_df.iloc[sorted_idx_to_rem[k]]['Seizure Type (FAS; FIAS; FBTC; Non-electrographic; Subclinical; Unknown)']
            dst_folder = save_root + pat_id + "/"
            if not os.path.exists(dst_folder):
                os.mkdir(dst_folder)
            save_fname = pat_id+"_"+event_num+"_"+sz_type+"_interictal_monopole_ERROR_no_onset_time.txt"
            with open(dst_folder + save_fname, 'w') as fp:
                pass
            append_error_to_master_error(dst_folder + save_fname)
            pat_time_df = pat_time_df.drop(sorted_idx_to_rem[k])
    pat_time_df = pat_time_df.reset_index()
    
    sz_idxs = pat_time_df[pat_time_df['Type'] == 'Seizure'].index.values.astype(int)
    
    cols_to_loop_through = ['onset_datetime','offset_datetime','onset_datetime_buffer','offset_datetime_buffer']
        
    pat_file_df = change_str_to_datetime(pat_file_df, cols_to_loop_through)
    pat_time_df = change_str_to_datetime(pat_time_df, cols_to_loop_through)
    
    
    # read in the timeframes available for the sleep epochs
    clean_sleep_pat_df = clean_sleep_df[clean_sleep_df["PatID"] == pat_id]
    
    # loop through each sleep category
    for ii in range(len(sleep_categories)):
        tmp_sleep_cat_df = clean_sleep_pat_df[clean_sleep_pat_df["SleepCat"] == sleep_categories[ii]]
        
        # if no categories are available, skip
        if len(tmp_sleep_cat_df) == 0:
            continue
        
        # if the category is already completed, skip
        dst_folder = sleep_epoch_dir + pat_id + "/"
        if not os.path.isdir(dst_folder):
            os.mkdir(dst_folder)
        
        tmp_completed_fpaths = glob.glob(dst_folder + "*" + sleep_categories[ii] + "*_monopole.mat")
        if len(tmp_completed_fpaths) > 0:
            continue
        
        # sort the dataframe so that the highest certainty values are first
        tmp_sleep_cat_df = tmp_sleep_cat_df.sort_values(by="AvgCertainty", ascending=False, ignore_index=True)
        
        # try to extract the highest certainty stage
        for jj in range(len(tmp_sleep_cat_df)):
            onset_datetime = datetime.datetime.strptime(tmp_sleep_cat_df["OnsetDatetime"].iloc[jj], "%Y-%m-%d %H:%M:%S")
            total_duration = tmp_sleep_cat_df["Duration"].iloc[jj]
            onset_datetime = onset_datetime + datetime.timedelta(seconds=math.floor((total_duration - sleep_epoch_duration_s)/2))
            offset_datetime = onset_datetime + datetime.timedelta(seconds=sleep_epoch_duration_s)
            

            
            save_fname = pat_id + "_" + tmp_sleep_cat_df["SleepCat"].iloc[jj] + "_" + "POD" + str(tmp_sleep_cat_df["POD"].iloc[jj]) + "_monopole.mat"
            
            tmp_edf_file, success_code = write_edf_data_from_times_sleep_epochs(pat_id, onset_datetime, offset_datetime, dst_folder, save_fname)
            # if successfully written, go to next phase
            if success_code == 1:
                break
    # except:
    #     print('skipping')
            