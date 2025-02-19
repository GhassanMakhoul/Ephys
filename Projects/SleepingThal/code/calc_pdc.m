function [mpdc] = calc_pdc(ft_data, shuffle, cfg)
    %% Stated function: 
    % EDF_chunk ->PDC from n_trials-> PDC        \
    %                                             -> compare(PDC,Null_pdc) -> sig_pdc 
    % FFT(EDF_chunk) -> random_phase -> Null_PDC /
    %%%%%
    %adding path to repo's fieldtrip
    tic
    addpath('/home/ghassan/Documents/Research/Ephys/Code/fieldtrip-master/');
    %fieldtrip expects chars use single ' quotes

    if numel(cfg) == 0
        %get MVAR fit (assuming it's 5)
        cfg = [];
        cfg.order = 5; %TODO make modular
        cfg.method = 'bsmart';
    end
    
    % this should really get its own script too. 
    % TODO get to the bottom of this damn null modeling business
    disp(shuffle)
    if shuffle
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
    
    
    %feeding in data from above simulation
    mdata = ft_mvaranalysis(cfg, ft_data);
    mdata.coeffs;

    %% computing the parametric spectral transfer matrix
    cfg             = [];
    cfg.method      = 'mvar';
    cfg.keeptrials  = 'yes';
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