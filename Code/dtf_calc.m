function [T,A,B,G_low,G_high] = dtf_calc(filt_data,hdr,mvar_order,epoch_seg_length)

ft_data.label = hdr.label';       % cell-array containing strings, Nchan*1
ft_data.fsample = hdr.frequency(1);   % sampling frequency in Hz, single number
ft_data.trial{1,1} = filt_data; % cell-array containing a data matrix for each trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
ft_data.time{1,1} = (0:size(filt_data,2)-1)/ft_data.fsample;      % cell-array containing a time axis for each trial (1*Ntrial), each time axis is a 1*Nsamples vector
ft_data.nSamples = size(filt_data,2);

cfg = [];
cfg.continuous   = 'yes'; 
[ft_data] = ft_preprocessing(cfg, ft_data);

cfg = [];
cfg.length = epoch_seg_length;
data_segmented = ft_redefinetrial(cfg, ft_data);

cfg = [];
cfg.order = mvar_order;
cfg.method = 'bsmart';
mdata = ft_mvaranalysis(cfg,data_segmented);

cfg = [];
cfg.method = 'mvar';
mfreq = ft_freqanalysis(cfg,mdata);

cfg = [];
cfg.method = 'dtf';
dtf = ft_connectivityanalysis(cfg,mfreq);


theta = [4 8];
BANDS_TO_ANALYZE = theta; 
idx = find(dtf.freq >= BANDS_TO_ANALYZE(1) & dtf.freq <= BANDS_TO_ANALYZE(2));
T = squeeze(mean(dtf.dtfspctrm(:,:,idx),3));
q=size(T,1);                                %number of nodes
T(1:q+1:end)= 0;                          %replace diagonal with 0s

alpha = [8 12]; 
BANDS_TO_ANALYZE = alpha; 
idx = find(dtf.freq >= BANDS_TO_ANALYZE(1) & dtf.freq <= BANDS_TO_ANALYZE(2));
A = squeeze(mean(dtf.dtfspctrm(:,:,idx),3));
q=size(A,1);                                %number of nodes
A(1:q+1:end)= 0;                          %replace diagonal with 0s

beta = [13 30]; 
BANDS_TO_ANALYZE = beta; 
idx = find(dtf.freq >= BANDS_TO_ANALYZE(1) & dtf.freq <= BANDS_TO_ANALYZE(2));
B = squeeze(mean(dtf.dtfspctrm(:,:,idx),3));
q=size(B,1);                                %number of nodes
B(1:q+1:end)= 0;                          %replace diagonal with 0s

gamma_low = [31 80];
BANDS_TO_ANALYZE = gamma_low; 
idx = find(dtf.freq >= BANDS_TO_ANALYZE(1) & dtf.freq <= BANDS_TO_ANALYZE(2));
G_low = squeeze(mean(dtf.dtfspctrm(:,:,idx),3));
q=size(G_low,1);                                %number of nodes
G_low(1:q+1:end)= 0;                          %replace diagonal with 0s

gamma_high = [81 150];
BANDS_TO_ANALYZE = gamma_high; 
idx = find(dtf.freq >= BANDS_TO_ANALYZE(1) & dtf.freq <= BANDS_TO_ANALYZE(2));
G_high = squeeze(mean(dtf.dtfspctrm(:,:,idx),3));
q=size(G_high,1);                                %number of nodes
G_high(1:q+1:end)= 0;                          %replace diagonal with 0s


end