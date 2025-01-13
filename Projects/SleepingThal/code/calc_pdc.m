%% This script is now becoming the main operating script for our pipeline
%% With this script, we should be able to specify a directory of input files
%% and calculate PDC for them. 
%% Stated function: 
% EDF_chunk ->PDC from n_trials-> PDC        \
%                                             -> compare(PDC,Null_pdc) -> sig_pdc 
% FFT(EDF_chunk) -> random_phase -> Null_PDC /

%%%%%
addpath('/home/ghassan/Documents/Research/Ephys/Projects/SleepingThal/code/shared_toolboxes/fieldtrip-20190819/');
addpath("/mnt/ernie_main/000_Data/SEEG/SEEG_Sleep_Staging/data/connectivity/5sDur_1sStride_v2/");
addpath('/home/ghassan/Documents/Research/Ephys/Projects/SleepingThal/code/shared_toolboxes/')
filt_edf = load("/mnt/ernie_main/000_Data/SEEG/SEEG_Sleep_Staging/data/extracted_sleep_epochs_v2/Spat52/Spat52_N2_POD1.mat");
interictal_file = '/mnt/ernie_main/000_Data/SEEG/SEEG_EyesClosed_RestingState/data/Spat53/Filtered_Chunks/First_Collection/All_2minChunks_Bipole_UnusedChannelsDeleted_Filtered/Spat53_Raw_bipole_unusedChannelsDeleted_Unfiltered_FILTERED_WITH_FIR_filtfilt_hamming_windowed_sinc_EEGlab_1-55_65-119Hz_chunk_0s_120s.edf'; 
%fieldtrip expects chars use single ' quotes
SHUFFLE = 1;
subj ='Spat53'
% Define 5-second trials
cfg = [];
cfg.dataset = interictal_file;
cfg.trialfun = 'ft_trialfun_general'; % Use the general trial function
cfg.trialdef.triallength = 5; % Trial length in seconds
cfg.trialdef.ntrials = 24;
% cfg.trialdef.overlap=4; % No overlap between trials (optional)
cfg = ft_definetrial(cfg);
% Preprocess the data with the defined trials
data = ft_preprocessing(cfg);

nnodes = length(data.label);
ttimepoints = data.fsample*cfg.trialdef.triallength;
mshuff = cfg.trialdef.ntrials;
if SHUFFLE
    fname = sprintf('%s_SHUFFLED',subj);
    %shuffling on first trial only
    shuff_trial = zeros(ttimepoints,nnodes,mshuff);
    X = data.trial{1}';
    for n=1:nnodes
        shuff= phaseran(X(:,n),mshuff);
        shuff_trial(:,n,:) = [shuff; zeros(1,1,mshuff)];
    end

    % check if odd in future
    disp("Shuffled trials")
    for n=1:cfg.trialdef.ntrials
        data.trial{n} = shuff_trial(:,:,n)';
    end
else
    fname = subj;
end



disp(data)
% fs = filt_edf.sampling_freq;
% seeg_filt = filt_edf.filt_data;
% n_samps = size(seeg_filt,2);
% time = n_samps/fs;
% trial_len = 5;

% %% configure the raw SEEG trial
% data = [];
% [data.trial, data.time]  = sig_to_trials(seeg_filt, trial_len, 1, fs);

% data.fsample = [fs];
% data.label = cellstr(filt_edf.bip_montage_label);
%%
%get MVAR fit (assuming it's 5)
cfg = [];
cfg.order = 5;
cfg.method = 'bsmart';
%feeding in data from above simulation
mdata = ft_mvaranalysis(cfg, data);
mdata.coeffs;

%% computing the parametric spectral transfer matrix
cfg             = [];
cfg.method      = 'mvar';
mfreq           = ft_freqanalysis(cfg, mdata);
mfreq;
%%
figure
cfg           = [];
cfg.method    = 'pdc';
mpdc       = ft_connectivityanalysis(cfg, mfreq);

cfg           = [];
cfg.parameter = 'pdcspctrm';
cfg.zlim      = [0 1];

ft_connectivityplot(cfg, mpdc);



% Increase overall figure size
set(gcf, 'Position', [1, 1, 5000, 5000]); % [x, y, width, height]

% Get all subplot axes
axesHandles = findall(gcf, 'Type', 'axes');

% % Adjust each axis's position
% for i = 1:length(axesHandles)
%     % Example: Increase the size by scaling the position
%     pos = get(axesHandles(i), 'Position'); % Get current position [x, y, width, height]
%     pos(3) = pos(3) * 2; % Increase width
%     pos(4) = pos(4) * 2; % Increase height
%     set(axesHandles(i), 'Position', pos); % Apply new position
% end
saveas(gcf, sprintf('../viz/%s_interictal_pdc_23 trials.png', fname))
%% Summarize connectivity
pdc = mpdc;
theta = [4 8];
BANDS_TO_ANALYZE = theta;
idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
T = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
Tint = squeeze(max(pdc.pdcspctrm(:,:,idx),[],3));

q = size(T,1);                                %number of nodes
T(1:q+1:end)= nan;                          %replace diagonal with NaNs
Tint(1:q+1:end)= nan;                          %replace diagonal with NaNs


alpha = [8 12];
BANDS_TO_ANALYZE = alpha;
idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
A = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
q = size(A,1);                                %number of nodes
A(1:q+1:end)= nan;                          %replace diagonal with NaNs

beta = [13 30];
BANDS_TO_ANALYZE = beta;
idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
B = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
q = size(B,1);                                %number of nodes
B(1:q+1:end)= nan;                          %replace diagonal with NaNs

low_gamma = [31 80];
BANDS_TO_ANALYZE = low_gamma;
idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
G_low = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
q = size(G_low,1);                                %number of nodes
G_low(1:q+1:end)= nan;                          %replace diagonal with NaNs

high_gamma = [81 150];
BANDS_TO_ANALYZE = high_gamma;
idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
G_high = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
q = size(G_high,1);                                %number of nodes
G_high(1:q+1:end)= nan;                          %replace diagonal with NaNs
%%
figure
imagesc(T)
title("Theta Power")
colorbar
saveas(gcf, sprintf('../viz/%s_Theta_PDC_summary.png', fname))

figure
imagesc(G_high)
title("G high ")
colorbar

saveas(gcf, sprintf('../viz/%s_Gamma_high_PDC_summary.png',fname))

%%

function [trials,time] = sig_to_trials(sig, win_size, stride, fs)
%
%TODO preallocate
%if I zero padded then I'd have T/strid number of trials
n_samps = size(sig,2);
time = n_samps/fs;
num_trials = (time-stride-win_size)/stride;
s=1;
e=win_size*fs
trials = cell(1,num_trials);
time = cell(1,num_trials)
for i=1:num_trials
trials{i} = sig(:,s:e);
time{i} = 0:1/fs:win_size-1/fs;
s=s+stride*fs;
e=e+stride*fs;
end
end
