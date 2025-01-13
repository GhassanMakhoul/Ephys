function [FT_data_filt_C] = FT_filt_C(data, sampFreq)

[b1,a1] = butter(3,1/(sampFreq/2),'high');
% figure(2)
% freqz(b1,a1,1e6,sampFreq)
f1_data = filtfilt(b1,a1,data')';

[b2,a2] = butter(3,[59,61]/(sampFreq/2),'stop');
% figure(3)
% freqz(b2,a2,1e6,sampFreq)
f2_data = filtfilt(b2,a2,f1_data')';

[b3,a3] = butter(8,119/(sampFreq/2),'low');
% figure(6)
% freqz(b5,a5,1e6,sampFreq)
f3_data = filtfilt(b3,a3,f2_data')';

FT_data_filt_C = f3_data;

end

