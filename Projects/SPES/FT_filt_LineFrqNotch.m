function [data_filt] = FT_filt_LineFrqNotch(data, sampFreq)

[b1,a1] = butter(5,[57,63]/(sampFreq/2),'stop');
f1_data = filtfilt(b1,a1,data')';

[b2,a2] = butter(5,[117,123]/(sampFreq/2),'stop');
f2_data = filtfilt(b2,a2,f1_data')';

data_filt = f2_data;

end

