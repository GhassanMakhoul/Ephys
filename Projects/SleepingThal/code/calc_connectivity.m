function [network_state] = calc_connectivity(subj, inp_f, out_dir,metric, shuffle, trial_len, overlap, ntrials, summary,viz_path, soz_label_file, demographics_xls_file)
    
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
    % tst_file = '/mnt/ernie_main/000_Data/SEEG/SEEG_Sleep_Staging/data/extracted_sleep_epochs_v2//Spat34/Spat34_N2_POD1.mat'
    %%%%%
    %global vars
    addpath('/home/ghassan/Documents/Research/Ephys/Projects/SleepingThal/code/shared_toolboxes/edfRead')
    BANDS = ["delta"; "theta"; "alpha"; "beta"; "gamma_low";"gamma_high" ];
    disp(trial_len)
    %% input parsing to force numbers from bash input (important if using shyaml)
    if ischar(trial_len)
        %if true, then likely this file is called from bash 
        % numbers will likely be strings then
        trial_len = str2double(trial_len);
        ntrials = str2double(ntrials);
        shuffle = str2double(shuffle);
        overlap = str2double(overlap);
    end
    %% read in file and plot raw tracings, important for QA
    % [hdr, record] = edfread(inp_f); for now let's just read from .mat files and get out of EDF ASAP
    bip_data = load(inp_f);
    filt_data = bip_data.filt_data;
    [n_ch, n_samps] = size(filt_data);
    TT = array2timetable(filt_data', 'SampleRate', bip_data.sampling_freq, 'VariableNames', bip_data.bip_montage_label);
    %save edf visualization (good to sanity check for any artifacts) 
    f = figure('Position', [0,0,1000,10000]);          
    tiledlayout(n_ch,1);

    for i=1:n_ch             
        nexttile                                
        plot(TT.Time, filt_data(i,:))              
        ylabel(bip_data.bip_montage_label{i})                    
    end
    [~,fname,~] = fileparts(inp_f);

    saveas(f, sprintf("%s/%s_stacked_plot.png", "../viz/pdc_calculations",fname))          
    fprintf("Saved %s SEEG trace for chunk %s", subj, fname)
    ft_data_trials = mat_to_ft_trial(inp_f, 'trial_len', trial_len, 'overlap', overlap);
    network_state = struct;
    %% Calculate main connectivity metric of choice
    if strcmp(metric, 'pdc')
        %calc PDC connectivity
        mpdc = calc_pdc(ft_data_trials, shuffle, []);
        %plot full PDC connectivity
        cfg           = [];
        cfg.parameter = 'pdcspctrm';
        cfg.zlim      = [0 1];
        ft_connectivityplot(cfg, mpdc);
        % Increase overall figure size
    
        set(gcf, 'Position', [1, 1, 5000, 5000]); % [x, y, width, height]
        saveas(gcf, sprintf('%s/%s_Full_PDC.png',viz_path, subj))
        % Summarize connectivity
        %TODO modularize
        pdc = mpdc;
        [D, T, A, B, G_low, G_high] = summarize_pdc(pdc, summary);
        network_state.pdc = pdc;
        network_state.pdc_summarized = zeros(size(D));
        %
        banded_pdc = {D, T, A, B, G_low, G_high};
        for n=1:length(BANDS)
            % save out summary of bands connectivity
            conn_mat = banded_pdc{n};
            close all; figure
            imagesc(conn_mat)
            title(sprintf("%s Power", BANDS(n)))
            colorbar
            saveas(gcf, sprintf('%s/%s_%s_PDC_summary.png', viz_path , fname, BANDS(n)))
            
        end
    end


    %% generate output struct and add in obligatory measures
    % may want to consider making this modular and giving it its own script file
    % that way I can append to a struct if it exists
    % and I can also add more connectivity measures
    %save raw data
    network_state.raw_data = bip_data.filt_data;
    network_state.sampling_freq = bip_data.sampling_freq;
    % these conversiions to Chars are for python compatibility
    network_state.dk_atlas_label = convertStringsToChars(bip_data.bip_montage_region);
    network_state.bip_montage_region =convertStringsToChars(bip_data.bip_montage_label); 
    % save pdc
    % 
    network_state.pdc = pdc;
    % save wavelet decomp
    wvlt_cfg = [];
    wvlt_cfg.output = 'pow';
    wvlt_cfg.taper = 'hanning';
    wvlt_cfg.tampsofrq = 1;
    wvlt_cfg.method = 'mtmfft';
    %may change in future
    wvlt_cfg.foi = 1:128;
    disp("POW SPECTRUM")
    network_state.pow_hanning = calc_wavelet(ft_data_trials, wvlt_cfg,'trial_len',trial_len, overlap);
    
    
    % FOOOF 
    disp("FOOOF!")
    wvlt_cfg = [];
    wvlt_cfg.output = 'fooof_peaks';
    wvlt_cfg.method = 'mtmfft';
    wvlt_cfg.taper = 'hanning';
    network_state.pow_foof = calc_wavelet(ft_data_trials, wvlt_cfg);

    %%delta phase locking values 0.5-4hz
    disp("PhaseLocking Values!")
    wvlt = calc_wavelet(ft_data_trials, []);
    plv_cfg = [];
    plv_cfg.method = 'plv';
    plv_cfg.keeptrials = 'yes';
    network_state.plv = ft_connectivityanalysis(plv_cfg, wvlt);
    disp("DONE: with PLV")

    %%TODO - future analysis should check for thalamic placement and do 
    network_state.PAC = [];
    if any(contains(bip_data.bip_montage_region, 'thal','IgnoreCase', true))
        %mostly interested in spindle cross-freq-coupling
        % spindles oscillate at 11-15Hz and have a 
        % low frequency envelope with period 0.5s - 3 s (1/3Hz - 2Hz)
        disp("WE GOT SOME THALAMUS DATA!")
        cross_spect_cfg =[];
        cross_spect_cfg.freqlow_rng = 0.5:.5:4;
        cross_spect_cfg.freqhigh_rng = 11:0.5:15;
        cfg.chanlow_char = 'thal'; % thal channels
        cfg.chanhigh_char = 'ctx'; % ctx channels
        cfg.freqlow = [0.5 2];
        cfg.freqhigh = [11 14];
        cfg.method = 'pac';
        network_state.PAC = calc_crossspectral(ft_data_trials, cross_spect_cfg);
        %do PAC
    end
    %time series autocorr,
    %  supposed to be dynamic EI balance Nanda... Rubinov 2023: 
    % Time-resolved correlation of distributed brain activity tracks E-I balance and accounts for diverse scale-free phenomena
    % start with a lag of +1
    autocorr = calc_autocorr(inp_f, 1);
    network_state.autocorr = autocorr;

    %% annotate file with relevant metadata: SOZ, PZ, NIZ locations, also epilepsy type
    ch_table = readtable(soz_label_file);
    bool_ix =  contains(ch_table.Var1, subj);
    network_state.channel_designation = ch_table(bool_ix, [2,3]); % calling meta data getter with empty cell so that it will use default row values 
    network_state.demographics_data = get_pat_info(subj,demographics_xls_file, {});
    disp(network_state)
    out_f = sprintf("%s/%s_connectivity.mat", out_dir, fname);
    save(out_f, "network_state")
    exit; 
end