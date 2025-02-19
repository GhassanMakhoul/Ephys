function [ft_trial] = mat_to_ft_trial(inp_f, varargin)
    %%% This function takes an input .mat file which opens a struct of seeg 
    %%% data that has been preprocessed. Struct should have the folowing fields:
    %%%             sampling_freq, filt_data, bip_montage label, bip_montage_region (dk is  optional)
    %%% Function then defines a trial using fieldtrip's define trial functionality,
    %%% default triallength is 5 seconds with no overlap
    %% parse the input arguments
    addpath('/home/ghassan/Documents/Research/Ephys/Code/fieldtrip-master/');
    disp(inp_f)
    options = struct('trial_len',5, 'overlap',0);
    optionNames = fieldnames(options);
    %# count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
        error('EXAMPLE needs propertyName/propertyValue pairs')
    end

    %# count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
    error('EXAMPLE needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = lower(pair{1}); %# make case insensitive

        if any(strcmp(inpName,optionNames))
            %# overwrite options. If you want you can test for the right class here
            %# Also, if you find out that there is an option you keep getting wrong,
            %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
            options.(inpName) = pair{2};
        else
            error('%s is not a recognized parameter name',inpName)
        end
    end
    %% DONE parsing default arguments
    fprintf("Defining a trial with length %s and overlap %s", options.trial_len, options.overlap)
    addpath('/home/ghassan/Documents/Research/Ephys/Code/fieldtrip-master/');
    data = load(inp_f);

    ft_data.label = cellstr(data.bip_montage_label);
        
    % cell-array containing strings, Nchan*1

    ft_data.fsample = data.sampling_freq;   % sampling frequency in Hz, single number
    ft_data.trial{1,1} = data.filt_data; % cell-array containing a data matrix for each trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
    ft_data.time{1,1} = (0:size(data.filt_data,2)-1)/ft_data.fsample;      % cell-array containing a time axis for each trial (1*Ntrial), each time axis is a 1*Nsamples vector
    ft_data.nSamples = size(data.filt_data,2);

    cfg = [];
    cfg.continuous   = 'yes';
    [ft_data_preprocessed] = ft_preprocessing(cfg, ft_data);

    cfg = [];
    cfg.length = options.trial_len;
    cfg.overlap = options.overlap;
    [ft_trial] = ft_redefinetrial(cfg, ft_data_preprocessed);
end