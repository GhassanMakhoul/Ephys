function [freq] = calc_wavelet(ft_data,wvlt_cfg, varargin)
    
    cfg=wvlt_cfg;
    %should have method
    %should have output
    %should have tapers and tsper smoothing if necessary
    if numel(cfg) == 0
        % default value picks up better low frequency 
        % often use for plv 
        cfg.method='mtmfft';
        cfg.output='fourier';
        cfg.taper='dpss';
        cfg.tapsmofrq = 1;
    end
    %these params are fixed
    cfg.trials='all';
    if any(contains(cfg.output, 'fourier', 'IgnoreCase', true))
        cfg.keeptrials = 'yes';
        cfg.keeptapers = 'yes';
    elseif any(contains(cfg.output, 'fooof','IgnoreCase', true))
        cfg.keeptrials = 'no';
        cfg.keeptapers = 'no';
    else
        cfg.keeptrials = 'yes';
        cfg.keeptapers = 'no';
    end
    cfg.toi = 'all';
    cfg.width = 7;
    cfg.gwitdh = 3;
    cfg.pad='nextpow2';
    [freq] = ft_freqanalysis(cfg, ft_data);
   
end