%% Split into chunks and save filtered signals

% This script will take in a raw EDF and filter it in the the ways outlined
% below, and create a .set and .edf of each
% 1) EEG_lab pop_eegfiltnew (Hamming windowed sinc filter). 
%   Output name appended with:
%       1a: _FILTERED_WITH_FIR_filtfilt_hamming_windowed_sinc_EEGlab_1-59_61-119Hz.edf
%       1b: _FILTERED_WITH_FIR_filtfilt_hamming_windowed_sinc_EEGlab_1-59_61-119Hz.set
% A) filtfilt 
%   Output name appended with:
%       2a: _FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119_121-179_181-239Hz.edf
%       2b: _FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119_121-179_181-239Hz.set
% B) filtfilt 
%   Output name appended with:
%       3a: _FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119Hz.edf
%       3b: _FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119Hz.set
% C) filtfilt 
%   Output name appended with:
%       4a: _FILTERED_WITH_IIR_filtfilt_butterworth_1-59_61-119Hz.edf
%       4b: _FILTERED_WITH_IIR_filtfilt_butterworth_1-59_61-119Hz.set
% D) filtfilt 
%   Output name appended with:
%       5a: _FILTERED_WITH_IIR_filtfilt_butterworth_1-55_65-119Hz.edf
%       5b: _FILTERED_WITH_IIR_filtfilt_butterworth_1-55_65-119Hz.set
%
%
% Script will also split edf into chunks if desired
% NOTE: filtering is done before splitting into chunks to avoid boundary
% condition artifacts

clear all
close all
restoredefaultpath

patID = 'Spat60';
ernie_main_letter = 'Z';

% Path to fieldtrip toolbox
ft_path = [ernie_main_letter ':\shared_toolboxes\fieldtrip-20190819'];
addpath(ft_path);
ft_defaults

% EEGlab path
eeglab_path = [ernie_main_letter ':\shared_toolboxes\eeglab2021.1'];
addpath(eeglab_path); 

% If raw EDF is in 2 collections, run the first one, then run second one to
% tmp directory and rename
root_pat_folder = [ernie_main_letter ':\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\' patID ];
unfiltered_edf_path = [root_pat_folder '\Unfiltered_Chunks\First_Collection\Raw_Bipole_UnusedChannelsDeleted_Unfiltered\' patID '_Raw_bipole_unusedChannelsDeleted_Unfiltered.edf'];

% unfilt_chunk_dir us UNUSED if SPLIT_EDF is 0
unfilt_chunk_dir = [root_pat_folder '\Unfiltered_Chunks\First_Collection\All_2minChunks_Bipole_UnusedChannelsDeleted_UnFiltered'];

filt_chunk_dir = [root_pat_folder '\Filtered_Chunks\First_Collection\All_2minChunks_Bipole_UnusedChannelsDeleted_Filtered'];
tmp_fpath = [root_pat_folder '\Filtered_Chunks'];
if ~isfolder(tmp_fpath)
    mkdir(tmp_fpath);
end

tmp_fpath = [root_pat_folder '\Filtered_Chunks\First_Collection'];
if ~isfolder(tmp_fpath)
    mkdir(tmp_fpath);
end

tmp_fpath = [root_pat_folder '\Filtered_Chunks\First_Collection\All_2minChunks_Bipole_UnusedChannelsDeleted_Filtered'];
if ~isfolder(tmp_fpath)
    mkdir(tmp_fpath);
end

SPLIT_EDF = 1; %0 = leave edf whole, %1 = split into CHUNK_SIZE second chunks
CHUNK_SIZE = 120; % seconds

%% Run 

% Read in edf with FT
cfg = [];
cfg.dataset = char(unfiltered_edf_path);
cfg.continuous = 'yes';
hdr = ft_read_header(char(unfiltered_edf_path));
FT_data = ft_preprocessing(cfg);
label = FT_data.label;
FT_dataRef = cell2mat(FT_data.trial);
sampFreq = hdr.Fs; %get sampling rate from header
numChans = hdr.nChans;

% Read in edf with EEGLab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEGLab = pop_biosig(unfiltered_edf_path,'importevent','off','importannot','off','blockepoch','off');

% Filter with all possible filters 
EEGLab_filt = eeglab_filt(EEGLab);
restoredefaultpath % needed to run filtfilt
FT_data_filt_A = FT_filt_A(FT_dataRef,sampFreq);
FT_data_filt_B = FT_filt_B(FT_dataRef,sampFreq);
FT_data_filt_C = FT_filt_C(FT_dataRef,sampFreq);
FT_data_filt_D = FT_filt_D(FT_dataRef,sampFreq);

% Must add toolbox paths back and initialize again
addpath(ft_path);
ft_defaults
addpath(eeglab_path); 
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Prepare output names
eeglab_filt_append_str_set = '_FILTERED_WITH_FIR_filtfilt_hamming_windowed_sinc_EEGlab_1-55_65-119Hz.set';
eeglab_filt_append_str_edf = '_FILTERED_WITH_FIR_filtfilt_hamming_windowed_sinc_EEGlab_1-55_65-119Hz.edf';
FT_filt_A_append_str_set = '_FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119_121-179_181-239Hz.set';
FT_filt_A_append_str_edf = '_FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119_121-179_181-239Hz.edf';
FT_filt_B_append_str_set = '_FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119Hz.set';
FT_filt_B_append_str_edf = '_FILTERED_WITH_IIR_filtfilt_butterworth_0.01-59_61-119Hz.edf';
FT_filt_C_append_str_set = '_FILTERED_WITH_IIR_filtfilt_butterworth_1-59_61-119Hz.set';
FT_filt_C_append_str_edf = '_FILTERED_WITH_IIR_filtfilt_butterworth_1-59_61-119Hz.edf';
FT_filt_D_append_str_set = '_FILTERED_WITH_IIR_filtfilt_butterworth_1-55_65-119Hz.set';
FT_filt_D_append_str_edf = '_FILTERED_WITH_IIR_filtfilt_butterworth_1-55_65-119Hz.edf';

splits = split(unfiltered_edf_path,'\');
end_splits = split(splits(end),'.');

eeglab_filt_save_name_set = char(strcat(end_splits(1),eeglab_filt_append_str_set));
eeglab_filt_save_name_edf = char(strcat(end_splits(1),eeglab_filt_append_str_edf));
FT_filt_A_save_name_set = char(strcat(end_splits(1),FT_filt_A_append_str_set));
FT_filt_A_save_name_edf = char(strcat(end_splits(1),FT_filt_A_append_str_edf));
FT_filt_B_save_name_set = char(strcat(end_splits(1),FT_filt_B_append_str_set));
FT_filt_B_save_name_edf = char(strcat(end_splits(1),FT_filt_B_append_str_edf));
FT_filt_C_save_name_set = char(strcat(end_splits(1),FT_filt_C_append_str_set));
FT_filt_C_save_name_edf = char(strcat(end_splits(1),FT_filt_C_append_str_edf));
FT_filt_D_save_name_set = char(strcat(end_splits(1),FT_filt_D_append_str_set));
FT_filt_D_save_name_edf = char(strcat(end_splits(1),FT_filt_D_append_str_edf));


% Split into chunks if desired
if SPLIT_EDF == 1

    % Determine how many X second chunks can be made
    total_seconds = size(FT_dataRef,2)/sampFreq;
    num_chunks = floor(total_seconds/CHUNK_SIZE);
    
    % For each chunk
    for i = 1:num_chunks
        
        append_chunk_edf = sprintf('_chunk_%is_%is.edf',(i-1)*CHUNK_SIZE,i*CHUNK_SIZE);
        append_chunk_set = sprintf('_chunk_%is_%is.set',(i-1)*CHUNK_SIZE,i*CHUNK_SIZE);
        
        % Split the raw input data and save as unfiltered edf
        raw_split_data = FT_dataRef(:,((i-1)*CHUNK_SIZE*sampFreq + 1):((i-1)*CHUNK_SIZE*sampFreq + CHUNK_SIZE*sampFreq));
        split_hdr = hdr;
        split_hdr.nSamples = CHUNK_SIZE*sampFreq;
        split_save_name = char(strcat(end_splits(1),append_chunk_edf));
        ft_write_data(fullfile(unfilt_chunk_dir,split_save_name),-1*raw_split_data,'header',split_hdr,'dataformat','edf'); %Absolutley no idea how signal got flipped (added '-1*' to correct)
        
        % Split filtered data and save
        EEGLab_split_edf = EEGLab_filt;
        EEGLab_split_edf.pnts = CHUNK_SIZE*sampFreq;
        EEGLab_split_edf.xmax = CHUNK_SIZE;
        EEGLab_split_edf.times = EEGLab_filt.times(((i-1)*CHUNK_SIZE*sampFreq + 1):((i-1)*CHUNK_SIZE*sampFreq + CHUNK_SIZE*sampFreq));
        EEGLab_split_edf.data = -1*EEGLab_filt.data(:,((i-1)*CHUNK_SIZE*sampFreq + 1):((i-1)*CHUNK_SIZE*sampFreq + CHUNK_SIZE*sampFreq)); %Absolutley no idea how signal got flipped (added '-1*' to correct)
        
        split_data_A = FT_data_filt_A(:,((i-1)*CHUNK_SIZE*sampFreq + 1):((i-1)*CHUNK_SIZE*sampFreq + CHUNK_SIZE*sampFreq));
        split_data_B = FT_data_filt_B(:,((i-1)*CHUNK_SIZE*sampFreq + 1):((i-1)*CHUNK_SIZE*sampFreq + CHUNK_SIZE*sampFreq));
        split_data_C = FT_data_filt_C(:,((i-1)*CHUNK_SIZE*sampFreq + 1):((i-1)*CHUNK_SIZE*sampFreq + CHUNK_SIZE*sampFreq));
        split_data_D = FT_data_filt_D(:,((i-1)*CHUNK_SIZE*sampFreq + 1):((i-1)*CHUNK_SIZE*sampFreq + CHUNK_SIZE*sampFreq));
        
        % Save as filtered chunks
        % Save eeglab data
        s = split(eeglab_filt_save_name_edf,'.edf');
        split_eeglab_save_name_edf = char(strcat(s(1),append_chunk_edf))
        split_eeglab_save_name_set = char(strcat(s(1),append_chunk_set))
        pop_saveset(EEGLab_split_edf, 'filename',split_eeglab_save_name_set,'filepath',char(filt_chunk_dir));
        pop_writeeeg(EEGLab_split_edf, fullfile(filt_chunk_dir,split_eeglab_save_name_edf), 'TYPE','EDF');
        
        % Save FT data as SET and EDF
        % A)
        s = split(FT_filt_A_save_name_edf,'.edf');
        split_A_save_name_edf = char(strcat(s(1),append_chunk_edf));
        split_A_save_name_set = char(strcat(s(1),append_chunk_set));
        ft_write_data(fullfile(filt_chunk_dir,split_A_save_name_edf),split_data_A,'header',split_hdr,'dataformat','edf');
        EEG_filtfilt = pop_biosig(fullfile(filt_chunk_dir,split_A_save_name_edf),'importevent','off','importannot','off','blockepoch','off');
        pop_saveset(EEG_filtfilt, 'filename',split_A_save_name_set,'filepath',char(filt_chunk_dir));
        
        % B)
        s = split(FT_filt_B_save_name_edf,'.edf');
        split_B_save_name_edf = char(strcat(s(1),append_chunk_edf));
        split_B_save_name_set = char(strcat(s(1),append_chunk_set));
        ft_write_data(fullfile(filt_chunk_dir,split_B_save_name_edf),split_data_B,'header',split_hdr,'dataformat','edf');
        EEG_filtfilt = pop_biosig(fullfile(filt_chunk_dir,split_B_save_name_edf),'importevent','off','importannot','off','blockepoch','off');
        pop_saveset(EEG_filtfilt, 'filename',split_B_save_name_set,'filepath',char(filt_chunk_dir));
        
        % C)
        s = split(FT_filt_C_save_name_edf,'.edf');
        split_C_save_name_edf = char(strcat(s(1),append_chunk_edf));
        split_C_save_name_set = char(strcat(s(1),append_chunk_set));
        ft_write_data(fullfile(filt_chunk_dir,split_C_save_name_edf),split_data_C,'header',split_hdr,'dataformat','edf');
        EEG_filtfilt = pop_biosig(fullfile(filt_chunk_dir,split_C_save_name_edf),'importevent','off','importannot','off','blockepoch','off');
        pop_saveset(EEG_filtfilt, 'filename',split_C_save_name_set,'filepath',char(filt_chunk_dir));
        
        % D)
        s = split(FT_filt_D_save_name_edf,'.edf');
        split_D_save_name_edf = char(strcat(s(1),append_chunk_edf));
        split_D_save_name_set = char(strcat(s(1),append_chunk_set));
        ft_write_data(fullfile(filt_chunk_dir,split_D_save_name_edf),split_data_D,'header',split_hdr,'dataformat','edf');
        EEG_filtfilt = pop_biosig(fullfile(filt_chunk_dir,split_D_save_name_edf),'importevent','off','importannot','off','blockepoch','off');
        pop_saveset(EEG_filtfilt, 'filename',split_D_save_name_set,'filepath',char(filt_chunk_dir));
        
    end
% Else save whole thing as .edf and .set
else
    
    % Save eeglab data
    pop_saveset(EEGLab_filt, 'filename',eeglab_filt_save_name_set,'filepath',char(filt_chunk_dir));
    pop_writeeeg(EEGLab_filt, fullfile(filt_chunk_dir,eeglab_filt_save_name_edf), 'TYPE','EDF');
    
    % Save FT data as SET and EDF
    % A)
    ft_write_data(fullfile(filt_chunk_dir,FT_filt_A_save_name_edf),FT_data_filt_A,'header',hdr,'dataformat','edf');
    EEG_filtfilt = pop_biosig(fullfile(filt_chunk_dir,FT_filt_A_save_name_edf),'importevent','off','importannot','off','blockepoch','off');
    pop_saveset(EEG_filtfilt, 'filename',FT_filt_A_save_name_set,'filepath',char(filt_chunk_dir));
    
    % B)
    ft_write_data(fullfile(filt_chunk_dir,FT_filt_B_save_name_edf),FT_data_filt_B,'header',hdr,'dataformat','edf');
    EEG_filtfilt = pop_biosig(fullfile(filt_chunk_dir,FT_filt_B_save_name_edf),'importevent','off','importannot','off','blockepoch','off');
    pop_saveset(EEG_filtfilt, 'filename',FT_filt_B_save_name_set,'filepath',char(filt_chunk_dir));
    
    % C)
    ft_write_data(fullfile(filt_chunk_dir,FT_filt_C_save_name_edf),FT_data_filt_C,'header',hdr,'dataformat','edf');
    EEG_filtfilt = pop_biosig(fullfile(filt_chunk_dir,FT_filt_C_save_name_edf),'importevent','off','importannot','off','blockepoch','off');
    pop_saveset(EEG_filtfilt, 'filename',FT_filt_C_save_name_set,'filepath',char(filt_chunk_dir));
    
    % D)
    ft_write_data(fullfile(filt_chunk_dir,FT_filt_D_save_name_edf),FT_data_filt_D,'header',hdr,'dataformat','edf');
    EEG_filtfilt = pop_biosig(fullfile(filt_chunk_dir,FT_filt_D_save_name_edf),'importevent','off','importannot','off','blockepoch','off');
    pop_saveset(EEG_filtfilt, 'filename',FT_filt_D_save_name_set,'filepath',char(filt_chunk_dir));
    
end

close all

disp("Script Completed")