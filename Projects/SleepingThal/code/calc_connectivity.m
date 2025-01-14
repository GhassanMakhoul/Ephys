function [] = calc_connectivity(INP_F, OUT_F,trial_len, ntrials)
    %% This script is now becoming the main operating script for our pipeline
    %% With this script, we should be able to specify a directory of input files
    %% and calculate PDC for them. 
    %% Stated function: 
    % EDF_chunk ->PDC from n_trials-> PDC        \
    %                                             -> compare(PDC,Null_pdc) -> sig_pdc 
    % #TODO:FFT(EDF_chunk) -> random_phase -> Null_PDC /
    % Visualizations:
    %  1. input signal, should be saved to output path
    %  2. (N_nodex x N_node tiles with PDC vs freq plot for all areas
    %  3. PDC summarization
    %
    % Output structure fields
    %   patID
    %   sleep_stage : W, N2, N3, R (No N1)
    %   start_time : 
    %   end_time : -corresponds to end time for edf
    %   bip_labels_used : ? 
    %   region_names : dk_atlas regions
    %   soz_final : label of SOZ, PZ, NIZ, IZ etc
    %   raw_soz_string :
    %   window_dur_sec : trial_len
    %   window_stride_sec : stride_len
    %   fs: 
    %   PDC_all_windowed_verbose : [f_bands x i_strides x n_channels x n_channels]
    %   PDC_all_windowed : [6_CANONICAL_BANDS x i_strides x n_channels x n_channels]
    %
    % -----In Future update------
    %   May add: RMS, PWELCH, CROSS Spectral Connectivity?
    % NOTE: to play nice with python we will be using chars. 
    %%%%%
    addpath('/home/ghassan/Documents/Research/Ephys/Projects/SleepingThal/code/shared_toolboxes/fieldtrip-20190819/');
    addpath("/mnt/ernie_main/000_Data/SEEG/SEEG_Sleep_Staging/data/connectivity/5sDur_1sStride_v2/");
    addpath('/home/ghassan/Documents/Research/Ephys/Projects/SleepingThal/code/shared_toolboxes/')
    interictal_file = INP_F; 
    %fieldtrip expects chars use single ' quotes
    SHUFFLE = 1;
    subj ='Spat53'
    % Define 5-second trials
    cfg = [];
    cfg.dataset = interictal_file;
    cfg.trialfun = 'ft_trialfun_general'; % Use the general trial function
    cfg.trialdef.triallength = trial_len; % Trial length in seconds
    cfg.trialdef.ntrials = ntrials;
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
    exit; 
    %%
end