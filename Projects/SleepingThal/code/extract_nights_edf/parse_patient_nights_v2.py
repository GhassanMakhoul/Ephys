# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 17:25:45 2023

@author: Derek
"""
#%%
import pandas as pd
import datetime
import mne
import os
import glob
from pyedflib import highlevel
import numpy as np
from math import ceil
from tqdm import tqdm

# ignore the dataframe warnings
pd.options.mode.chained_assignment = None  # default='warn'

root_data_folder = "Z:/000_Data/SEEG/SEEG_Sleep_Staging/data/nights_of_sleep/"
root_resting_folder = "Z:/000_Data/SEEG/SEEG_EyesClosed_RestingState/data/"
root_source_folder = "Z:/000_Data/SEEG/SEEG_Entire_EMU_Downloads/data/"

sleep_epoch_dir = "Z:/000_Data/SEEG/SEEG_Sleep_Staging/data/extracted_sleep_epochs_v2/"

log_fpath = "Z:/000_Data/SEEG/SEEG_Sleep_Staging/data/nights_of_sleep/log_file.txt"

conn_dir = "Z:/000_Data/SEEG/SEEG_Periictal/data/Connectivity/seizure_structs_pre10min_ictal_post10min/5sDur_1sStride/"
all_pat_folders = glob.glob(conn_dir+"*pat*")

pat_folders_completed = glob.glob(root_data_folder+"*pat*")
# pat_folders_completed = pat_folders_completed.replace("\\", "/")

pats_to_rem = list()
confident_sleep_categories = ["R", "W", "N2", "N3"]
pat_epoch_folders = glob.glob(sleep_epoch_dir + "*pat*")
for i in range(len(pat_epoch_folders)):
    tmp_epoch_folder = pat_epoch_folders[i].replace("\\", "/")
    tmp_file_types = confident_sleep_categories.copy()
    for tmp_monopole_file in glob.iglob(tmp_epoch_folder + "/*_monopole.mat"):
        tmp_monopole_file = tmp_monopole_file.replace("\\","/").split("/")[-1]
        tmp_sleep_type = tmp_monopole_file.split("_")[1]
        # could throw an error if N1 is found and it isn't in the list
        if tmp_sleep_type != "N1":
            tmp_file_types.remove(tmp_sleep_type)
    
    if len(tmp_file_types) == 0:
        pats_to_rem.append(pat_epoch_folders[i].replace("\\", "/").split("/")[-1])

# for i in range(len(pat_folders_completed)):
#     pats_to_rem.append(pat_folders_completed[i].replace("\\","/").split("/")[-1])

# pats_to_rem_extra = ["Epat20"]

pat_ids = list()

for i in range(len(all_pat_folders)):
    tmp_pat_name = all_pat_folders[i].replace("\\", "/").split("/")[-1]
    if tmp_pat_name not in pats_to_rem:
        pat_ids.append(tmp_pat_name)

#%% functions
def write_error_file(root_data_folder, pat_id, night_num, reason):
    path_to_file = root_data_folder + pat_id + "/" + pat_id + "_Night" + str(night_num) + "_" + reason + ".txt"
    with open(path_to_file, 'a') as f:
        f.write("PatientName: " + pat_id + reason)

def filter_edf_data(edf_file, fs):
    """FT Filt A, but in python

    Args:
        edf_file (_type_): mne python edf file
        fs (_type_): sampling frequency

    Returns:
        edf_file : filtered mne python edf file 
    """
    
    edf_file.load_data()
    
    # high pass filter
    iir_params = dict(order=3, ftype='butter', output='sos')
    iir_params = mne.filter.construct_iir_filter(iir_params=iir_params, f_pass=0.5, f_stop=None, sfreq=fs, btype='highpass', return_copy=False)
    edf_file = edf_file.filter(l_freq=None, h_freq=None, method="iir", iir_params=iir_params)

    # bandstop
    iir_params = dict(order=3, ftype='butter', output='sos')
    iir_params = mne.filter.construct_iir_filter(iir_params=iir_params, f_pass=[59, 61], f_stop=None, sfreq=fs, btype='bandstop', return_copy=False)
    edf_file = edf_file.filter(l_freq=None, h_freq=None, method="iir", iir_params=iir_params)

    # bandstop
    iir_params = dict(order=3, ftype='butter', output='sos')
    iir_params = mne.filter.construct_iir_filter(iir_params=iir_params, f_pass=[119, 121], f_stop=None, sfreq=fs, btype='bandstop', return_copy=False)
    edf_file = edf_file.filter(l_freq=None, h_freq=None, method="iir", iir_params=iir_params)

    # lowpass
    iir_params = dict(order=8, ftype='butter', output='sos')
    iir_params = mne.filter.construct_iir_filter(iir_params=iir_params, f_pass=150, f_stop=None, sfreq=fs, btype='lowpass', return_copy=False)
    edf_file = edf_file.filter(l_freq=None, h_freq=None, method="iir", iir_params=iir_params)
    
    return edf_file

def extract_night_from_edf_files(edf_fpath, start_datetime, stop_datetime, channels_used,
                                 save_folder, pat_id, POD, POD_segment, date_string):
    # check to see if it is a clabel file
    if "clabel" in edf_fpath:
        tmp_root_folder_idx = edf_fpath.rfind("/")
        tmp_root_folder = edf_fpath[:tmp_root_folder_idx]
        edf_files = glob.glob(tmp_root_folder + "/*.EDF")
        for kk in range(len(edf_files)):
            if "clabel" not in edf_files[kk]:
                file_with_correct_labels = edf_files[kk]
                break
        edf_file_with_correct_labels = mne.io.read_raw_edf(file_with_correct_labels)
        
        # Let's figure out which channels to drop
        # all channels from the edf file
        edf_ch_names = edf_file_with_correct_labels.ch_names
        og_edf_ch_names = edf_ch_names
        # remove whitespaces
        for kk in range(len(edf_ch_names)):
            edf_ch_names[kk] = edf_ch_names[kk].replace(" ", "")
        # channels from the channels used excel file
        bip_names = channels_used[0].to_list()
        first_ch, second_ch = list(zip(*(k.split(" - ") for k in bip_names)))
        first_ch = list(first_ch)
        second_ch = list(second_ch)
        all_used_channels = first_ch + second_ch
        # remove duplicates
        all_used_channels = [*set(all_used_channels)]
        # remove whitespace
        for kk in range(len(all_used_channels)):
            all_used_channels[kk] = all_used_channels[kk].replace(" ", "")
        
        # channels that will be removed
        ch_to_drop = list(set(edf_ch_names).difference(all_used_channels))
        
        # get the indexes of which ch need to be dropped
        ch_idx_to_drop = list()
        for kk in range(len(ch_to_drop)):
            ch_idx_to_drop.append( edf_ch_names.index(ch_to_drop[kk]) )
        
        # ch names to keep
        ch_names_to_keep = list()
        for kk in range(len(edf_ch_names)):
            if edf_ch_names[kk] in all_used_channels:
                ch_names_to_keep.append(og_edf_ch_names[kk])
        
        # load in the edf file
        edf_file = mne.io.read_raw_edf(edf_fpath)
        info = edf_file.info
        
        edf_ch_names = edf_file.ch_names
        # ch_to_drop = edf_ch_names[ch_idx_to_drop]
        ch_to_drop = [edf_ch_names[i] for i in ch_idx_to_drop]
        
        # save_fpath = save_fpath.replace(".edf", "_clabel.edf")
        clabel_string = "_clabel"
        
        del edf_file_with_correct_labels
    else:
        clabel_string = ""
        # load in the edf file
        edf_file = mne.io.read_raw_edf(edf_fpath)
        info = edf_file.info
        
        # Let's figure out which channels to drop
        # all channels from the edf file
        edf_ch_names = edf_file.ch_names
        # channels from the channels used excel file
        bip_names = channels_used[0].to_list()
        first_ch, second_ch = list(zip(*(k.split(" - ") for k in bip_names)))
        first_ch = list(first_ch)
        second_ch = list(second_ch)
        all_used_channels = first_ch + second_ch
        # remove duplicates
        all_used_channels = [*set(all_used_channels)]
        # channels that will be removed
        ch_to_drop = list(set(edf_ch_names).difference(all_used_channels))
    
    
    edf_file_starttime = info['meas_date'].replace(tzinfo=None)
    edf_file_endtime = edf_file_starttime + datetime.timedelta(seconds=edf_file._last_time)
    edf_file_duration_s = edf_file._last_time
           
    # filter the edf file
    fs = info["sfreq"]
    
    # get the start and end time
    starttime_in_seconds = (start_datetime - edf_file_starttime).seconds
    endtime_in_seconds = (stop_datetime - edf_file_starttime).seconds
    
    # list to store all of the save fpaths used
    save_fpath_list = list()
    
    # how long in s can the final file be given the size limits
    max_size_bytes = 40e9 # 40 gig in memory max
    
    # how long the saved segment can be
    max_save_duration_s = ((max_size_bytes / 8) / len(all_used_channels)) / fs
    
    # cut the file and save it
    segment_duration = (endtime_in_seconds - starttime_in_seconds)
    if segment_duration > max_save_duration_s:
        num_files_needed = ceil(segment_duration / max_save_duration_s)
        dur_of_each_file_s = segment_duration / num_files_needed
        for kk in range(num_files_needed):
            # save paths
            POD_segment = POD_segment+kk
            save_fname = pat_id + "_Night" + str(POD) + "-" + str(POD_segment) + "_" + date_string + "_filt" + clabel_string + ".edf"
            save_fpath = save_folder + save_fname
            
            save_fpath_list.append(save_fpath)
            
            # crop times
            crop_start_time_s = starttime_in_seconds + kk*dur_of_each_file_s
            crop_stop_time_s = starttime_in_seconds + (kk+1)*dur_of_each_file_s
            
            # cap the endtime at the end of the file
            if crop_stop_time_s >= edf_file_duration_s:
                crop_stop_time_s = edf_file_duration_s
            
            crop_start_datetime = start_datetime + datetime.timedelta(seconds=kk*dur_of_each_file_s)
            
            # load in the edf file
            edf_file = mne.io.read_raw_edf(edf_fpath)
            
            # remove unneeded channels
            edf_file = edf_file.drop_channels(ch_to_drop)
            
            # crop the file
            edf_file = edf_file.crop(tmin=crop_start_time_s, tmax=crop_stop_time_s)
            
            # set the new starttime of the edf file, which is the time it was cropped from
            start_datetime_utc = datetime.datetime(year=crop_start_datetime.year, 
                                            month=crop_start_datetime.month,
                                            day=crop_start_datetime.day, 
                                            hour=crop_start_datetime.hour,
                                            minute=crop_start_datetime.minute,
                                            second=crop_start_datetime.second,
                                            tzinfo=datetime.timezone.utc)
            edf_file.set_meas_date(start_datetime_utc)
            
            # filter the edf file
            edf_file = filter_edf_data(edf_file, fs)
            
            # save the file
            edf_file.export(save_fpath, verbose=True)
            
            # clean up the file for memory management
            edf_file.close()
            del edf_file
    else:
    
        save_fname = pat_id + "_Night" + str(POD) + "-" + str(POD_segment) + "_" + date_string + "_filt" + clabel_string + ".edf"
        save_fpath = save_folder + save_fname
        
        save_fpath_list.append(save_fpath)
    
        # drop unneeded channels
        edf_file = edf_file.drop_channels(ch_to_drop)
            
        # crop the file
        edf_file = edf_file.crop(tmin=starttime_in_seconds, tmax=endtime_in_seconds)
        
        # filter the file
        edf_file = filter_edf_data(edf_file, fs)
                
        # set the new starttime of the edf file, which is the time it was cropped from
        start_datetime_utc = datetime.datetime(year=start_datetime.year, 
                                        month=start_datetime.month,
                                        day=start_datetime.day, 
                                        hour=start_datetime.hour,
                                        minute=start_datetime.minute,
                                        second=start_datetime.second,
                                        tzinfo=datetime.timezone.utc)
        edf_file.set_meas_date(start_datetime_utc)

        # save the file
        edf_file.export(save_fpath, verbose=True)
        
        # clean up the file for memory management 
        edf_file.close()
        del edf_file
    
    # if it is a clabel file, then fix the channel names
    if clabel_string == "_clabel":
        for jj in range(len(save_fpath_list)):
            save_fpath = save_fpath_list[jj]
            clabel_file = highlevel.read_edf_header(save_fpath)
            clabel_ch = clabel_file["channels"]
            channel_rename_dict = {clabel_ch[i]: ch_names_to_keep[i] for i in range(len(clabel_ch))}
            
            renamed_save_fpath = save_fpath.replace("_clabel.edf", ".edf")
            highlevel.rename_channels(edf_file=save_fpath, mapping=channel_rename_dict, new_file=renamed_save_fpath)
        
    return POD_segment




#%%

# stopped at Spat05

# pat_order = [34, 17, 18, 24, 16, 19, 20, 21, 22, 23, 25, 26, 31, 32, 45, 2, 48, 54, 55, 38, 37, 12, 30, 7, 27, 15, 58, 39, 11, 51, 57, 1, 41, 42, 52, 53, 56, 44, 8, 46, 5, 49, 40, 13, 14, 28, 36, 47, 3, 4, 6, 43, 9, 10, 29, 0, 35, 50, 33]

# for p in range(len(pat_ids)):
for p in tqdm(range(len(pat_ids))):
    if pat_ids[p] == "Spat31":
        continue
        
    # try:
    pat_id = pat_ids[p]
    
    print(pat_id + " \n")
    
    root_patient_dest_folder = root_data_folder + pat_id + "/"
    if not os.path.isdir(root_patient_dest_folder):
        os.mkdir(root_patient_dest_folder)
    root_patient_source_folder = root_source_folder + pat_id + "/"
    
    event_csv = pd.read_csv("Z:/000_Data/SEEG/SEEG_Sleep_Staging/notes/all_time_data_01092023_112957.csv", delimiter="\t")
    pat_csv = event_csv[event_csv["Pat ID"] == pat_id]
    
    
    patient_resting_folder = glob.glob(root_resting_folder+pat_id+"*")
    patient_resting_folder = patient_resting_folder[0]
    patient_resting_folder = patient_resting_folder.replace("\\", "/")
    
    channels_used_fpath = patient_resting_folder + "/" + pat_id + "_Channels_Used.xlsx"
    channels_used = pd.read_excel(channels_used_fpath, header=None)
    
    # convert the pat_csv datetime strings to datetime
    for i in range(len(pat_csv)):
        try:
            if not isinstance(pat_csv["onset_datetime_buffer"].iloc[i], datetime.datetime) and not isinstance(pat_csv["onset_datetime_buffer"].iloc[i], float):
                if pat_csv["onset_datetime_buffer"].iloc[i] == 'None':
                    pat_csv["onset_datetime_buffer"].iloc[i] = None
                else:
                    try:
                        pat_csv["onset_datetime_buffer"].iloc[i] = datetime.datetime.strptime(pat_csv["onset_datetime_buffer"].iloc[i], "%Y-%m-%d %H:%M:%S.%f")
                    except:
                        pat_csv["onset_datetime_buffer"].iloc[i] = datetime.datetime.strptime(pat_csv["onset_datetime_buffer"].iloc[i], "%Y-%m-%d %H:%M:%S")
        except:
            print("skipping")
        
        try:
            if not isinstance(pat_csv["onset_datetime"].iloc[i], datetime.datetime) and not isinstance(pat_csv["onset_datetime"].iloc[i], float):
                if pat_csv["onset_datetime"].iloc[i] == 'None':
                    pat_csv["onset_datetime"].iloc[i] = None
                else:
                    try:
                        pat_csv["onset_datetime"].iloc[i] = datetime.datetime.strptime(pat_csv["onset_datetime"].iloc[i], "%Y-%m-%d %H:%M:%S.%f")
                    except:
                        pat_csv["onset_datetime"].iloc[i] = datetime.datetime.strptime(pat_csv["onset_datetime"].iloc[i], "%Y-%m-%d %H:%M:%S")
        except:
            print("skipping")
                    
        try:
            if not isinstance(pat_csv["offset_datetime_buffer"].iloc[i], datetime.datetime) and not isinstance(pat_csv["offset_datetime_buffer"].iloc[i], float):
                if pat_csv["offset_datetime_buffer"].iloc[i] == 'None':
                    pat_csv["offset_datetime_buffer"].iloc[i] = None
                else:
                    try:    
                        pat_csv["offset_datetime_buffer"].iloc[i] = datetime.datetime.strptime(pat_csv["offset_datetime_buffer"].iloc[i], "%Y-%m-%d %H:%M:%S.%f")
                    except:
                        pat_csv["offset_datetime_buffer"].iloc[i] = datetime.datetime.strptime(pat_csv["offset_datetime_buffer"].iloc[i], "%Y-%m-%d %H:%M:%S")
        except:
            print("skipping")
           
        try: 
            if not isinstance(pat_csv["offset_datetime"].iloc[i], datetime.datetime) and not isinstance(pat_csv["offset_datetime"].iloc[i], float):            
                if pat_csv["offset_datetime"].iloc[i] == 'None':
                    pat_csv["offset_datetime"].iloc[i] = None
                else:
                    try:
                        pat_csv["offset_datetime"].iloc[i] = datetime.datetime.strptime(pat_csv["offset_datetime"].iloc[i], "%Y-%m-%d %H:%M:%S.%f")
                    except:
                        pat_csv["offset_datetime"].iloc[i] = datetime.datetime.strptime(pat_csv["offset_datetime"].iloc[i], "%Y-%m-%d %H:%M:%S")
        except:
            print("skipping")
    
    # df with all files
    pat_file_csv = pat_csv[pat_csv["Type"] == "File"]
    
    # df with all seizures
    pat_sz_csv = pat_csv[pat_csv["Type"] == "Seizure"]
    
    # df with all stims
    pat_stim_csv = pat_csv[pat_csv["Type"] == "Stimulation"]
    
    # find the earliest night
    start_datetime = pat_file_csv["onset_datetime_buffer"].min()
    
    # find the latest night
    stop_datetime = pat_file_csv["offset_datetime_buffer"].max()
    
    
    # nights to loop through
    tot_nights = (stop_datetime - start_datetime).days

    # find out what nights were completed
    nightly_completed_files = glob.glob(root_patient_dest_folder + "*Night*")
    nights_completed = list()
    
    for i in range(tot_nights):
        
        # if the night was already completed, skip
        tmp_night_files = glob.glob(root_patient_dest_folder + "*Night" + str(i) + "*")
        if len(tmp_night_files) > 0:
            continue
        
        # 10pm the previous day
        start_of_night = datetime.datetime(year=start_datetime.year, 
                                        month=start_datetime.month,
                                        day=start_datetime.day, 
                                        hour=22,
                                        minute=0)
        start_of_night = start_of_night + datetime.timedelta(days=i)
        
        # 6am the next day
        end_of_night = datetime.datetime(year=start_datetime.year, 
                                        month=start_datetime.month,
                                        day=start_datetime.day, 
                                        hour=6,
                                        minute=0)
        end_of_night = end_of_night + datetime.timedelta(days=i+1)
        
        
        # get all the files for the night
        files_to_include = pat_file_csv[(pat_file_csv["onset_datetime_buffer"] < start_of_night) &
                                        (pat_file_csv["offset_datetime_buffer"] > end_of_night)]
        
        # if there isn't one file that encompasses the whole night, grab multiple
        if len(files_to_include) == 0:
            files_to_include = pat_file_csv[((pat_file_csv["onset_datetime_buffer"] < start_of_night) &
                                            (pat_file_csv["offset_datetime_buffer"] > start_of_night)) |
                                            ((pat_file_csv["onset_datetime_buffer"] > start_of_night) &
                                            (pat_file_csv["offset_datetime_buffer"] < end_of_night)) |
                                            ((pat_file_csv["onset_datetime_buffer"] < end_of_night) &
                                            (pat_file_csv["offset_datetime_buffer"] > end_of_night))]
            
        
        # if there are seizures during the night, trash it
        sz_df_during_night = pat_sz_csv[((pat_sz_csv["onset_datetime_buffer"] > start_of_night) &
                                        (pat_sz_csv["offset_datetime_buffer"] < end_of_night)) |
                                        ((pat_sz_csv["onset_datetime_buffer"] < start_of_night) & 
                                        (pat_sz_csv["offset_datetime_buffer"] > start_of_night)) |
                                        ((pat_sz_csv["onset_datetime_buffer"] < end_of_night) &
                                        (pat_sz_csv["offset_datetime_buffer"] > end_of_night))]
        if len(sz_df_during_night) > 0:
            write_error_file(root_data_folder, pat_id, i, "seizure")
            continue
        
        # if stim was completed during the night, trash it
        stim_df_during_night = pat_stim_csv[((pat_stim_csv["onset_datetime_buffer"] > start_of_night) &
                                            (pat_stim_csv["offset_datetime_buffer"] < end_of_night)) |
                                            ((pat_stim_csv["onset_datetime_buffer"] < start_of_night) & 
                                            (pat_stim_csv["offset_datetime_buffer"] > start_of_night)) |
                                            ((pat_stim_csv["onset_datetime_buffer"] < end_of_night) &
                                            (pat_stim_csv["offset_datetime_buffer"] > end_of_night))]
        if len(stim_df_during_night) > 0:
            write_error_file(root_data_folder, pat_id, i, "stim")
            continue
        
        
        # write the files 
        fnames_to_include = files_to_include["FileName"].to_list()
        
        POD_segment = 0
        
        for j in range(len(fnames_to_include)):
            # if there are mutiple files for the night, grab the onset time as the start time
            if files_to_include["onset_datetime_buffer"].iloc[j] < start_of_night:
                tmp_sleep_start_datetime = start_of_night
            else:
                tmp_sleep_start_datetime = files_to_include["onset_datetime_buffer"].iloc[j]
                
            # if there are multiple files for the night, grab the offset time as the end time
            if files_to_include["offset_datetime_buffer"].iloc[j] > end_of_night:
                tmp_sleep_stop_datetime = end_of_night
            else:
                tmp_sleep_stop_datetime = files_to_include["offset_datetime_buffer"].iloc[j]
            
            
            # edf file fpath to grab the night data from
            edf_source_fpath = root_patient_source_folder + fnames_to_include[j]
            
            POD_segment = extract_night_from_edf_files(edf_source_fpath, tmp_sleep_start_datetime,
                                        tmp_sleep_stop_datetime, channels_used,
                                        root_patient_dest_folder, pat_id, i, POD_segment,
                                        start_of_night.strftime("%m%d%Y"))
            
            POD_segment = POD_segment + 1
               
        break
    # except:
    #     with open(log_fpath, 'a') as f:
    #        f.write("ERROR WITH PatientName: " + pat_ids[p] + "\n")