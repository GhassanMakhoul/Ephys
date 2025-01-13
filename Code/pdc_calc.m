function [D,T,A,B,G_low,G_high] = pdc_calc(filt_data,hdr,mvar_order,epoch_seg_length)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updated by DJD on 2/20/2024 to include delta. Old code commented out
% below.

% Calculate 2-minute values
ft_data.label = hdr.label';       % cell-array containing strings, Nchan*1
ft_data.fsample = hdr.frequency(1);   % sampling frequency in Hz, single number
ft_data.trial{1,1} = filt_data; % cell-array containing a data matrix for each trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
ft_data.time{1,1} = (0:size(filt_data,2)-1)/ft_data.fsample;      % cell-array containing a time axis for each trial (1*Ntrial), each time axis is a 1*Nsamples vector
ft_data.nSamples = size(filt_data,2);

cfg = [];
cfg.continuous   = 'yes';
[ft_data_preprocessed] = ft_preprocessing(cfg, ft_data);

cfg = [];
cfg.length = epoch_seg_length;
data_segmented = ft_redefinetrial(cfg, ft_data_preprocessed);

cfg = [];
cfg.order = mvar_order;
cfg.method = 'bsmart';
mdata = ft_mvaranalysis(cfg,data_segmented);

cfg = [];
cfg.method = 'mvar';
mfreq = ft_freqanalysis(cfg,mdata);

cfg = [];
cfg.method = 'pdc';
pdc = ft_connectivityanalysis(cfg,mfreq);


delta = [1 3];
BANDS_TO_ANALYZE = delta;
idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
D = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
q = size(D,1);                                %number of nodes
D(1:q+1:end)= nan;                          %replace diagonal with NaNs

theta = [4 8];
BANDS_TO_ANALYZE = theta;
idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
T = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
q = size(T,1);                                %number of nodes
T(1:q+1:end)= nan;                          %replace diagonal with NaNs

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






% function [T,A,B,G_low,G_high] = pdc_calc(filt_data,hdr,mvar_order,epoch_seg_length)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % updated by DJD on 2/20/2024 to include delta. Old code commented out
% % below.
% 
% % Calculate 2-minute values
% ft_data.label = hdr.label';       % cell-array containing strings, Nchan*1
% ft_data.fsample = hdr.frequency(1);   % sampling frequency in Hz, single number
% ft_data.trial{1,1} = filt_data; % cell-array containing a data matrix for each trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
% ft_data.time{1,1} = (0:size(filt_data,2)-1)/ft_data.fsample;      % cell-array containing a time axis for each trial (1*Ntrial), each time axis is a 1*Nsamples vector
% ft_data.nSamples = size(filt_data,2);
% 
% cfg = [];
% cfg.continuous   = 'yes';
% [ft_data_preprocessed] = ft_preprocessing(cfg, ft_data);
% 
% cfg = [];
% cfg.length = epoch_seg_length;
% data_segmented = ft_redefinetrial(cfg, ft_data_preprocessed);
% 
% cfg = [];
% cfg.order = mvar_order;
% cfg.method = 'bsmart';
% mdata = ft_mvaranalysis(cfg,data_segmented);
% 
% cfg = [];
% cfg.method = 'mvar';
% mfreq = ft_freqanalysis(cfg,mdata);
% 
% cfg = [];
% cfg.method = 'pdc';
% pdc = ft_connectivityanalysis(cfg,mfreq);
% 
% 
% theta = [4 8];
% BANDS_TO_ANALYZE = theta;
% idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
% T = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
% q = size(T,1);                                %number of nodes
% T(1:q+1:end)= nan;                          %replace diagonal with NaNs
% 
% alpha = [8 12];
% BANDS_TO_ANALYZE = alpha;
% idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
% A = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
% q = size(A,1);                                %number of nodes
% A(1:q+1:end)= nan;                          %replace diagonal with NaNs
% 
% beta = [13 30];
% BANDS_TO_ANALYZE = beta;
% idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
% B = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
% q = size(B,1);                                %number of nodes
% B(1:q+1:end)= nan;                          %replace diagonal with NaNs
% 
% low_gamma = [31 80];
% BANDS_TO_ANALYZE = low_gamma;
% idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
% G_low = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
% q = size(G_low,1);                                %number of nodes
% G_low(1:q+1:end)= nan;                          %replace diagonal with NaNs
% 
% high_gamma = [81 150];
% BANDS_TO_ANALYZE = high_gamma;
% idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
% G_high = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
% q = size(G_high,1);                                %number of nodes
% G_high(1:q+1:end)= nan;                          %replace diagonal with NaNs

