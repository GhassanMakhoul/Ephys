function [mpdc] = calc_pdc(inp_f, shuffle, trial_len, ntrials)
    %% Stated function: 
    % EDF_chunk ->PDC from n_trials-> PDC        \
    %                                             -> compare(PDC,Null_pdc) -> sig_pdc 
    % FFT(EDF_chunk) -> random_phase -> Null_PDC /
    %%%%%
    %adding path to repo's fieldtrip
    tic
    addpath('/home/ghassan/Documents/Research/Ephys/Code/fieldtrip-master/');
    %fieldtrip expects chars use single ' quotes
    data = load(inp_f);
    % Define 5-second trials
    ft_data.label = cellstr(data.bip_montage_label);       % cell-array containing strings, Nchan*1
    ft_data.fsample = data.sampling_freq;   % sampling frequency in Hz, single number
    ft_data.trial{1,1} = data.filt_data; % cell-array containing a data matrix for each trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
    ft_data.time{1,1} = (0:size(data.filt_data,2)-1)/ft_data.fsample;      % cell-array containing a time axis for each trial (1*Ntrial), each time axis is a 1*Nsamples vector
    ft_data.nSamples = size(data.filt_data,2);
    
    cfg = [];
    cfg.continuous   = 'yes';
    [ft_data_preprocessed] = ft_preprocessing(cfg, ft_data);
    
    cfg = [];
    cfg.length = trial_len;
    data_segmented = ft_redefinetrial(cfg, ft_data_preprocessed);
 

    nnodes = length(ft_data.label);
    ttimepoints = ft_data.fsample*trial_len;
    % this should really get its own script too. 
    % TODO get to the bottom of this damn null modeling business
    disp(shuffle)
    if strcmp(shuffle, '1')
        mshuff = ntrials;
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
            ft_data.trial{n} = shuff_trial(:,:,n)';
        end
    end
    
    %get MVAR fit (assuming it's 5)
    cfg = [];
    cfg.order = 5;
    cfg.method = 'bsmart';
    %feeding in data from above simulation
    mdata = ft_mvaranalysis(cfg, data_segmented);
    mdata.coeffs;

    %% computing the parametric spectral transfer matrix
    cfg             = [];
    cfg.method      = 'mvar';
    mfreq           = ft_freqanalysis(cfg, mdata);

    %%
    figure
    cfg           = [];
    cfg.method    = 'pdc';
    mpdc       = ft_connectivityanalysis(cfg, mfreq);

    disp("called PDC")
    toc
    %TODO save out PDC files 
%%
end