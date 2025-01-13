function [FT_data_filt_A] = FT_filt_A(data, sampFreq)

[b1,a1] = butter(3,1/(sampFreq/2),'high');
% figure(2)
% freqz(b1,a1,1e6,sampFreq)
f1_data = filtfilt(b1,a1,data')';

[b2,a2] = butter(3,[59,61]/(sampFreq/2),'stop');
% figure(3)
% freqz(b2,a2,1e6,sampFreq)
f2_data = filtfilt(b2,a2,f1_data')';

[b3,a3] = butter(3,[119,121]/(sampFreq/2),'stop');
% figure(4)
% freqz(b3,a3,1e6,sampFreq)
f3_data = filtfilt(b3,a3,f2_data')';

[b4,a4] = butter(8,150/(sampFreq/2),'low');
% figure(6)
% freqz(b5,a5,1e6,sampFreq)
f4_data = filtfilt(b4,a4,f3_data')';

FT_data_filt_A = f4_data;

end

