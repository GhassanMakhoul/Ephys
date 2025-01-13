function [T,A,B,G_low,G_high] = imcoh_calc(filt_data,hdr)

% Do the Rolston stuff
params.tapers = [3 5]; % [NW K] - time-bandwidth product, number of tapers)
params.Fs = hdr.frequency(1);
win = 1;
[Sc,Cmat,Ctot,Cvec,Cent,f] = CrossSpecMatc(filt_data',win,params);
%Sc = cross spectral matrix freq x chan x chan
%Cmat = coherence matrix freq x chan x chan
%Ctot = total coherence SV(1)^2/sum(SV^2)
%f = frequencies

Rtemp = abs(imag(Cmat));
%Rtemp = abs(real(Cmat));

BTA = [[4 8];[8 12];[13 30];[31 80];[81 150]];

% Fisher Z-Transformation
for b = 1:size(BTA,1)
    idx = find(f >= BTA(b,1) & f <= BTA(b,2));
    R{b} = squeeze(mean(Rtemp(idx,:,:),1));
    
    Z_temp{b} = atanh(R{b});
    idx = find(Z_temp{b} == inf);
    Z_temp{b}(idx) = nan;
    if ~isempty(idx)
        disp(['found nan in freq ' num2str(b)])
    end
end
T = Z_temp{1};
A = Z_temp{2};
B = Z_temp{3};
G_low = Z_temp{4};
G_high = Z_temp{5};

end