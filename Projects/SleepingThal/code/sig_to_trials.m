function [trials,time] = sig_to_trials(sig, win_size, stride, fs)
    %
    %TODO preallocate
    %if I zero padded then I'd have T/strid number of trials
    n_samps = size(sig,2);
    time = n_samps/fs;
    num_trials = (time-stride-win_size)/stride;
    s=1;
    e=win_size*fs
    trials = cell(1,num_trials);
    time = cell(1,num_trials)
    for i=1:num_trials
            trials{i} = sig(:,s:e);
            time{i} = 0:1/fs:win_size-1/fs;
            s=s+stride*fs;
            e=e+stride*fs;
    end
end