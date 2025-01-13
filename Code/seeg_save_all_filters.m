%% Filters_all
function [] = seeg_save_all_filters(UNFILTERED_edf_path, ft_path, eeglab_path, out_dir)

% This script will take in a raw EDF and filter it in the the ways outlined
% below, and create a .set and .edf of each
% 1) EEG_lab pop_eegfiltnew (Hamming windowed sinc filter). 
%   Output name appended with:
%       1a: _FILTERED_WITH_FIR_filtfilt_hamming_windowed_sinc_EEGlab_1-59_61-119Hz.edf
%       1b: _FILTERED_WITH_FIR_filtfilt_hamming_windowed_sinc_EEGlab_1-59_61-119Hz.set
% 2) filtfilt 
%   Output name appended with:
%       2a: _FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119_121-179_181-239Hz.edf
%       2b: _FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119_121-179_181-239Hz.set
% 3) filtfilt 
%   Output name appended with:
%       3a: _FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119Hz.edf
%       3b: _FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119Hz.set
% 4) filtfilt 
%   Output name appended with:
%       4a: _FILTERED_WITH_IIR_filtfilt_butterworth_1-59_61-119Hz.edf
%       4b: _FILTERED_WITH_IIR_filtfilt_butterworth_1-59_61-119Hz.set
% 5) filtfilt 
%   Output name appended with:
%       5a: _FILTERED_WITH_IIR_filtfilt_butterworth_1-55_65-119Hz.edf
%       5b: _FILTERED_WITH_IIR_filtfilt_butterworth_1-55_65-119Hz.set


clear all
close all
restoredefaultpath

% Path to fieldtrip toolbox
% ft_path = 'C:\Users\johnsgw3\Documents\MATLAB\Toolboxes\fieldtrip-20190819';
addpath(ft_path);
ft_defaults

% EEGlab path
% eeglab_path = 'C:\Users\johnsgw3\Documents\MATLAB\Toolboxes\eeglab2019_1';
addpath(eeglab_path); 

% % Create Output directory
% % Get timestamp to not overwrite output directory 
% out_dir = 'Z:\000_Data\SEEG\data\raw+preprocessed_seeg\Epat27\Filtered_Chunks\Nine_Chunks_Bipole_UnusedChannelsDeleted_Filtered';
% 
% UNFILTERED_edf_path = 'Z:\000_Data\SEEG\data\raw+preprocessed_seeg\Epat27\Unfiltered_Chunks\Nine_Chunks_Bipole_UnusedChannelsDeleted_UnFiltered\Epat27_bipole_unfiltered_unusedDeleted_18m20m.edf';
%% 1) EEG lab - 1-119 Hz BP and 55-65 HZ Notch Hamming sinc FIR filter. 

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEGR1 = pop_biosig(UNFILTERED_edf_path,'importevent','off','importannot','off','blockepoch','off');

% Bandpass filter
BP_EEGR1 = pop_eegfiltnew(EEGR1, 1,119,[],0,[],0);
BP_EEGR1.setname = 'bp_filt';
BP_EEGR1 = eeg_checkset(BP_EEGR1);

% Notch filter
F_EEGR1 = pop_eegfiltnew(BP_EEGR1, 55,65,[],1,[],0);
F_EEGR1.setname = 'full_filt';
F_EEGR1 = eeg_checkset(F_EEGR1);

% Append the name
append_str_set = '_FILTERED_WITH_FIR_filtfilt_hamming_windowed_sinc_EEGlab_1-55_65-119Hz.set';
append_str_edf = '_FILTERED_WITH_FIR_filtfilt_hamming_windowed_sinc_EEGlab_1-55_65-119Hz.edf';
splits = split(UNFILTERED_edf_path,'\');
end_splits = split(splits(end),'.');
save_name_set = char(strcat(end_splits(1),append_str_set));
save_name_edf = char(strcat(end_splits(1),append_str_edf));
pop_saveset( F_EEGR1, 'filename',save_name_set,'filepath',char(out_dir));
pop_writeeeg(F_EEGR1, fullfile(out_dir,save_name_edf), 'TYPE','EDF');

% %Loading Subject from SET file using Fieldtrip
% cfg.dataset = char(fullfile(out_dir,save_name));
% cfg.continuous = 'yes';
% hdr = ft_read_header(char(fullfile(out_dir,save_name)));
% 
% data_ft = ft_preprocessing(cfg);
% dataRef = cell2mat(data_ft.trial);
% label = data_ft.label;
% sampFreq = hdr.Fs;

%Cleaning Up
cfg.method = 'channel';
cfg.channel = 'all';

% figure(1)
% plot(1:num_pts,dataRef(1,start_pt:start_pt+num_pts-1))
% legend('EDFBrowser','EEGLab')
% hold on

%% 2/A) Filtfilt IIR Butterworth

cfg = [];
cfg.dataset = char(UNFILTERED_edf_path);
cfg.continuous = 'yes';

% Initialize FieldTrip info
ft_defaults
hdr = ft_read_header(char(UNFILTERED_edf_path));
data = ft_preprocessing(cfg);
label = data.label;
dataRef = cell2mat(data.trial);
sampFreq = hdr.Fs; %get sampling rate from header
numChans = hdr.nChans;

restoredefaultpath % needed to use filtfilt


[b1,a1] = butter(3,0.01/(sampFreq/2),'high');
% figure(2)
% freqz(b1,a1,1e6,sampFreq)
f1_data = filtfilt(b1,a1,dataRef')';

[b2,a2] = butter(3,[59,61]/(sampFreq/2),'stop');
% figure(3)
% freqz(b2,a2,1e6,sampFreq)
f2_data = filtfilt(b2,a2,f1_data')';

[b3,a3] = butter(3,[119,121]/(sampFreq/2),'stop');
% figure(4)
% freqz(b3,a3,1e6,sampFreq)
f3_data = filtfilt(b3,a3,f2_data')';

[b4,a4] = butter(3,[179,181]/(sampFreq/2),'stop');
% figure(5)
% freqz(b4,a4,1e6,sampFreq)
f4_data = filtfilt(b4,a4,f3_data')';

[b5,a5] = butter(8,239/(sampFreq/2),'low');
% figure(6)
% freqz(b5,a5,1e6,sampFreq)
f5_data = filtfilt(b5,a5,f4_data')';

% figure(1)
% plot(1:num_pts,f5_data(1,start_pt:start_pt+num_pts-1))
% legend('EDFBrowser','EEGLab','filtfilt')
% hold on 
% 
% plot(1:num_pts,dataRef(1,start_pt:start_pt+num_pts-1))
% legend('EDFBrowser','EEGLab','filtfilt','raw')
% hold off

% Path to fieldtrip toolbox (TODO: TRY TO GET RID OF restoredefaultpath call)
addpath(ft_path);
ft_defaults

% Create save names
append_str_set = '_FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119_121-179_181-239Hz.set';
append_str_edf = '_FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119_121-179_181-239Hz.edf';
splits = split(UNFILTERED_edf_path,'\');
end_splits = split(splits(end),'.');
save_name_set = char(strcat(end_splits(1),append_str_set));
save_name_edf = char(strcat(end_splits(1),append_str_edf));

% Save as EDF
ft_write_data(fullfile(out_dir,save_name_edf),f5_data,'header',hdr,'dataformat','edf');

% Save as SET
% EEGlab path
addpath(eeglab_path); 

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEG_filtfilt = pop_biosig(fullfile(out_dir,save_name_edf),'importevent','off','importannot','off','blockepoch','off');
pop_saveset(EEG_filtfilt, 'filename',save_name_set,'filepath',char(out_dir));

%% 3/B) same as 2, but stop filters at 119 Hz

addpath(ft_path);
ft_defaults

cfg = [];
cfg.dataset = char(UNFILTERED_edf_path);
cfg.continuous = 'yes';

% Initialize FieldTrip info
ft_defaults
hdr = ft_read_header(char(UNFILTERED_edf_path));
data = ft_preprocessing(cfg);
label = data.label;
dataRef = cell2mat(data.trial);
sampFreq = hdr.Fs; %get sampling rate from header
numChans = hdr.nChans;

restoredefaultpath % needed to use filtfilt


[b1,a1] = butter(3,0.01/(sampFreq/2),'high');
% figure(2)
% freqz(b1,a1,1e6,sampFreq)
f1_data = filtfilt(b1,a1,dataRef')';

[b2,a2] = butter(3,[59,61]/(sampFreq/2),'stop');
% figure(3)
% freqz(b2,a2,1e6,sampFreq)
f2_data = filtfilt(b2,a2,f1_data')';

[b3,a3] = butter(8,119/(sampFreq/2),'low');
% figure(6)
% freqz(b5,a5,1e6,sampFreq)
f3_data = filtfilt(b3,a3,f2_data')';

% figure(1)
% plot(1:num_pts,f5_data(1,start_pt:start_pt+num_pts-1))
% legend('EDFBrowser','EEGLab','filtfilt')
% hold on 
% 
% plot(1:num_pts,dataRef(1,start_pt:start_pt+num_pts-1))
% legend('EDFBrowser','EEGLab','filtfilt','raw')
% hold off

% Path to fieldtrip toolbox (TODO: TRY TO GET RID OF restoredefaultpath call)
addpath(ft_path);
ft_defaults

% Create save names
append_str_set = '_FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119Hz.set';
append_str_edf = '_FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119Hz.edf';
splits = split(UNFILTERED_edf_path,'\');
end_splits = split(splits(end),'.');
save_name_set = char(strcat(end_splits(1),append_str_set));
save_name_edf = char(strcat(end_splits(1),append_str_edf));

% Save as EDF
ft_write_data(fullfile(out_dir,save_name_edf),f3_data,'header',hdr,'dataformat','edf');

% Save as SET
% EEGlab path
addpath(eeglab_path); 

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEG_filtfilt = pop_biosig(fullfile(out_dir,save_name_edf),'importevent','off','importannot','off','blockepoch','off');
pop_saveset(EEG_filtfilt, 'filename',save_name_set,'filepath',char(out_dir));


%% 4/C) Same as 3, but do 1 Hz filtering instead of dwon to 0.01 Hz
addpath(ft_path);
ft_defaults

cfg = [];
cfg.dataset = char(UNFILTERED_edf_path);
cfg.continuous = 'yes';

% Initialize FieldTrip info
ft_defaults
hdr = ft_read_header(char(UNFILTERED_edf_path));
data = ft_preprocessing(cfg);
label = data.label;
dataRef = cell2mat(data.trial);
sampFreq = hdr.Fs; %get sampling rate from header
numChans = hdr.nChans;

restoredefaultpath % needed to use filtfilt


[b1,a1] = butter(3,1/(sampFreq/2),'high');
% figure(2)
% freqz(b1,a1,1e6,sampFreq)
f1_data = filtfilt(b1,a1,dataRef')';

[b2,a2] = butter(3,[59,61]/(sampFreq/2),'stop');
% figure(3)
% freqz(b2,a2,1e6,sampFreq)
f2_data = filtfilt(b2,a2,f1_data')';

[b3,a3] = butter(8,119/(sampFreq/2),'low');
% figure(6)
% freqz(b5,a5,1e6,sampFreq)
f3_data = filtfilt(b3,a3,f2_data')';

% figure(1)
% plot(1:num_pts,f5_data(1,start_pt:start_pt+num_pts-1))
% legend('EDFBrowser','EEGLab','filtfilt')
% hold on 
% 
% plot(1:num_pts,dataRef(1,start_pt:start_pt+num_pts-1))
% legend('EDFBrowser','EEGLab','filtfilt','raw')
% hold off

% Path to fieldtrip toolbox (TODO: TRY TO GET RID OF restoredefaultpath call)
addpath(ft_path);
ft_defaults

% Create save names
append_str_set = '_FILTERED_WITH_IIR_filtfilt_butterworth_1-59_61-119Hz.set';
append_str_edf = '_FILTERED_WITH_IIR_filtfilt_butterworth_1-59_61-119Hz.edf';
splits = split(UNFILTERED_edf_path,'\');
end_splits = split(splits(end),'.');
save_name_set = char(strcat(end_splits(1),append_str_set));
save_name_edf = char(strcat(end_splits(1),append_str_edf));

% Save as EDF
ft_write_data(fullfile(out_dir,save_name_edf),f3_data,'header',hdr,'dataformat','edf');

% Save as SET
% EEGlab path
addpath(eeglab_path); 

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEG_filtfilt = pop_biosig(fullfile(out_dir,save_name_edf),'importevent','off','importannot','off','blockepoch','off');
pop_saveset(EEG_filtfilt, 'filename',save_name_set,'filepath',char(out_dir));

%% 5/D) Same as 4, but do 55-65 Notch instead of 59-61 Hz
addpath(ft_path);
ft_defaults

cfg = [];
cfg.dataset = char(UNFILTERED_edf_path);
cfg.continuous = 'yes';

% Initialize FieldTrip info
ft_defaults
hdr = ft_read_header(char(UNFILTERED_edf_path));
data = ft_preprocessing(cfg);
label = data.label;
dataRef = cell2mat(data.trial);
sampFreq = hdr.Fs; %get sampling rate from header
numChans = hdr.nChans;

restoredefaultpath % needed to use filtfilt

[b1,a1] = butter(3,1/(sampFreq/2),'high');
% figure(2)
% freqz(b1,a1,1e6,sampFreq)
f1_data = filtfilt(b1,a1,dataRef')';

[b2,a2] = butter(3,[55,65]/(sampFreq/2),'stop');
% figure(3)
% freqz(b2,a2,1e6,sampFreq)
f2_data = filtfilt(b2,a2,f1_data')';

[b3,a3] = butter(8,119/(sampFreq/2),'low');
% figure(6)
% freqz(b5,a5,1e6,sampFreq)
f3_data = filtfilt(b3,a3,f2_data')';

% figure(1)
% plot(1:num_pts,f5_data(1,start_pt:start_pt+num_pts-1))
% legend('EDFBrowser','EEGLab','filtfilt')
% hold on 
% 
% plot(1:num_pts,dataRef(1,start_pt:start_pt+num_pts-1))
% legend('EDFBrowser','EEGLab','filtfilt','raw')
% hold off

% Path to fieldtrip toolbox (TODO: TRY TO GET RID OF restoredefaultpath call)
addpath(ft_path);
ft_defaults

% Create save names
append_str_set = '_FILTERED_WITH_IIR_filtfilt_butterworth_1-55_65-119Hz.set';
append_str_edf = '_FILTERED_WITH_IIR_filtfilt_butterworth_1-55_65-119Hz.edf';
splits = split(UNFILTERED_edf_path,'\');
end_splits = split(splits(end),'.');
save_name_set = char(strcat(end_splits(1),append_str_set));
save_name_edf = char(strcat(end_splits(1),append_str_edf));

% Save as EDF
ft_write_data(fullfile(out_dir,save_name_edf),f3_data,'header',hdr,'dataformat','edf');

% Save as SET
% EEGlab path
addpath(eeglab_path); 

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEG_filtfilt = pop_biosig(fullfile(out_dir,save_name_edf),'importevent','off','importannot','off','blockepoch','off');
pop_saveset(EEG_filtfilt, 'filename',save_name_set,'filepath',char(out_dir));