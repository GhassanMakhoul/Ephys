function [crss_spct] = calc_crossspectral(inp_f, cfg)


    addpath('/home/ghassan/Documents/Research/Ephys/Code/fieldtrip-master/');
    data = load(inp_f);
    bip_ch = data.bip_montage_label
    if strcmp(cfg.chanlow_char, 'all')
        cfg.freqlow = 'all';
    else
        low_ix = contains(regions, cfg.chanlow_char, 'IgnoreCase',true);
        cfg.freqlow = bip_ch(low_ix);
    end

    if strcmp(cfg.chanhigh_char, 'all')
        cfg.freqhigh = 'all';
    else
        hi_ix = contains(regions, cfg.chanhigh_char, "IgnoreCase", true);
        cfg.freqhigh = bip_ch(hi_ix)
    end
    %fieldtrip expects chars use single ' quotes
    % first get wvlts
    wvlt = calc_wavelet(inp_f, []);

    crss_spct = ft_crossfrequencyanalysis(cfg, wvlt, wvlt);    
end