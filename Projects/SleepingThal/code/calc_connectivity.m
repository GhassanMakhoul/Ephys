function [] = calc_connectivity(subj, inp_f, out_dir,metric, shuffle, trial_len, ntrials, summary,viz_path)
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
    %global vars
    addpath('shared_toolboxes/edfRead')
    BANDS = ["delta"; "theta"; "alpha"; "beta"; "gamma_low";"gamma_high" ];


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
    disp(sprintf("Saved %s SEEG trace for chunk %s", subj, fname));
    connectivity = struct;
    if strcmp(metric, 'pdc')
        %calc PDC connectivity
        mpdc = calc_pdc(inp_f, shuffle, trial_len, ntrials);
        %plot full PDC connectivity
        cfg           = [];
        cfg.parameter = 'pdcspctrm';
        cfg.zlim      = [0 1];
        ft_connectivityplot(cfg, mpdc);
        % Increase overall figure size
    
        set(gcf, 'Position', [1, 1, 5000, 5000]); % [x, y, width, height]
        saveas(gcf, sprintf('%s/%s_Full_PDC.png','../viz/pdc_calculations', subj))
        % Summarize connectivity
        %TODO modularize
        pdc = mpdc;
        banded_pdc = summarize_pdc(pdc, summary);
        %
        for n=1:length(BANDS)
            % save out summary of bands connectivity
            conn_mat = banded_pdc(n);
            figure
            imagesc(conn_mat)
            title(sprintf("%s Power", BANDS(n)))
            colorbar
            saveas(gcf, sprintf('%s/%s_%s_PDC_summary.png', '../viz/pdc_calculations' , fname, band))
        end
    end

    %% save output struct 
    % may want to consider making this modular and giving it its own script file
    % that way I can append to a struct if it exists
    % and I can also add more connectivity measures
    % 
    save(sprintf("%s/%s_connectivity.mat", out_dir, fname), "connectivity")
    exit; 
    %%
end