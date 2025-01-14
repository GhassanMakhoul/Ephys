function [mpdc] = calc_pdc(inp_f, shuffle, trial_len, ntrials)
    %% Stated function: 
    % EDF_chunk ->PDC from n_trials-> PDC        \
    %                                             -> compare(PDC,Null_pdc) -> sig_pdc 
    % FFT(EDF_chunk) -> random_phase -> Null_PDC /
    %%%%%
    %adding path to repo's fieldtrip
    tic
    addpath('/home/ghassan/Documents/Research/Ephys/Projects/SleepingThal/code/shared_toolboxes/fieldtrip-20190819/');
    %fieldtrip expects chars use single ' quotes
    
    % Define 5-second trials
    cfg = [];
    cfg.dataset = inp_f;
    cfg.trialfun = 'ft_trialfun_general'; % Use the general trial function
    cfg.trialdef.triallength = trial_len; % Trial length in seconds
    cfg.trialdef.ntrials = ntrials;
    % cfg.trialdef.overlap=4; % No overlap between trials (optional)
    cfg = ft_definetrial(cfg);
    % Preprocess the data with the defined trials
    data = ft_preprocessing(cfg);

    nnodes = length(data.label);
    ttimepoints = data.fsample*cfg.trialdef.triallength;
    % this should really get its own script too. 
    % TODO get to the bottom of this damn null modeling business
    if shuffle
        mshuff = cfg.trialdef.ntrials;
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
    end
    disp(data)
    
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

    %%
    figure
    cfg           = [];
    cfg.method    = 'pdc';
    mpdc       = ft_connectivityanalysis(cfg, mfreq);

    disp("called PDC")
    toc
%%
end