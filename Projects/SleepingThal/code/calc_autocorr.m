function [autocorr] = calc_autocorr(inp_f, lag)

    %fieldtrip expects chars use single ' quotes
    data = load(inp_f);
    lfp_seeg =  data.filt_data; % 
    [n,t] = size(lfp_seeg);
    lfp_seeg = lfp_seeg - mean(lfp_seeg,2);
    lfp_seeg = lfp_seeg ./ std(lfp_seeg, [], 2);
    autocorr = zeros(n,t);
    for ii = 1:t-1
         inst_corr = corrcoef(lfp_seeg(:, ii), lfp_seeg(:, ii+lag));
         autocorr(:,ii) = inst_corr(1,2);
    end

end