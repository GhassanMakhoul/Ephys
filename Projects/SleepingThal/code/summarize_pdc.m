function [D,T,A,B,G_low,G_high] = summarize_pdc(pdc, method)
    %% summarizes PDC into canonical bands
    % uses summarization schema specifiec by METHOD
    %       if method is 'mean' uses the historical
    %           averaging technique to summarize data
    %       if method is 'peak' uses a peak finder to 
    %             specify areas of high coherence
    %           NOTE: this feature is under development and may
    %                   need some fine-tuning. Expect sparse connectivity
    if strcmp(method, 'mean') 


        
            delta = [1 3];
            BANDS_TO_ANALYZE = delta;
            idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
            D = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
            q = size(D,1);                                %number of nodes
            D(1:q+1:end)= nan;                          %replace diagonal with NaNs



           theta = [4 8];
            BANDS_TO_ANALYZE = theta;
            idx = find(pdc.freq >= BANDS_TO_ANALYZE(1) & pdc.freq <= BANDS_TO_ANALYZE(2));
            T = squeeze(mean(pdc.pdcspctrm(:,:,idx),3));
            q = size(T,1);                                %number of nodes
            T(1:q+1:end)= nan;                          %replace diagonal with NaNs

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
    elseif strcmp(method, "peak")
        disp("CURRENTLY NOT SUPPORTED")
    
    end

end