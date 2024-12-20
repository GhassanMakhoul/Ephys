%% SPES Spectral Power Plots

% Make a struct of pats_spes with following fields:
% Pat ID Clean
% Stim labels
% Response labels (should be same as Channels_Used.xlsx)
% Cell of 5D matrices pre-train z-scored for each current used (stim x distThresh x time x response x freq) with PSD value as entries
% Order of 5D labels
% List of distance thresholds in 5D matrix
% List of time windows in 6D matrix
% List of freqs in 5D matrix
% 5D matrix z-scored acorss all 5 dimensions (this wrong to do?, maybe it hekps if some patients just are not large responders to stim...)
% SOZ stim labels
% SOZ response labels

%NOTE this code adds an additional agressive line frequency notch filter (57-63, 117-123 5th order BW)
% (should have already been filtered out, but was seeing artufacts on channels close to stim)

close all
clear all
restoredefaultpath
tic

% Start paralell pool (if not already running)
cores = feature('numcores');
pool = gcp('nocreate'); % If no pool, do not create new one.
if isempty(pool)
    pool = parpool('local',cores);
    disp(['Pool has been started with Num Workers ' num2str(pool.NumWorkers)]);
else
    poolsize = pool.NumWorkers;
    disp(['Pool is already running with Num Workers ' num2str(pool.NumWorkers)]);
end

only_run_file_check = 0; % If 1, no calculations will be made, but the code will run to check the presence of all necessary files for every patient.

source_dir = "X:\000_Data\SPES\data\preprocessed";
pulse_folder_subpath = "%s_CCEP_single_pulses";
pretrain_folder_subpath = "%s_CCEP_pre_train_baselines";
mat_keyword = ".mat";
soz_label_file = "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\labels\all_pats_bipole_ez_labels.csv";
dist_file = "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\labels\all_pats_euc_dist_list.csv";
metric_loc_in_struct = 3;


pat_dir_list = ["Epat26", "Epat27", "Epat30",...
    "Epat31", "Epat34", "Epat35", "Epat37",...
    "Epat38", "Epat39", "Spat30",...
    "Spat31", "Spat34", "Spat36", "Spat37",...
    "Epat43", "Spat41", "Spat42",...
    "Spat44", "Spat48", "Spat49", "Spat50",...
    "Spat52", "Spat53"];

% 0 = no mesial SOZs, 1 = pure mesial SOZs, 2 = mesial + other SOZs. 
pat_ep_type = [1, 1, 0, ...
    1, 1, 1, 0, ...
    1, 0, 2, ...
    1, 0, 2, 1,...
    0, 1, 0, ...
    0, 0, 0, 0,...
    0, 0];

pat_chUsed_files = [...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat26_dpat05\Epat26_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat27_dpat06\Epat27_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat30\Epat30_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat31_apat111\Epat31_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat34_dpat11\Epat34_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat35_apat115\Epat35_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat37_dpat13\Epat37_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat38\Epat38_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat39\Epat39_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat30\Spat30_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat31\Spat31_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat34\Spat34_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat36\Spat36_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat37\Spat37_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Epat43_apat117_vpat04\Epat43_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat41\Spat41_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat42\Spat42_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat44\Spat44_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat48\Spat48_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat49\Spat49_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat50\Spat50_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat52\Spat52_Channels_Used.xlsx";...
    "X:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data\Spat53\Spat53_Channels_Used.xlsx";...
    ];


% Include variable ability to exclude distance
% Starting at 5 mm avoids complicated logic of excluding channels that overlap
% with stim bipole pair (assuming all contact spacing is smaller than 5 mm)
dist_gap = 15; % creates a distance band within which edges will be included
dist_threshes = 5:dist_gap:80; % mm, nodes (i.e. response channels) outside this will be included

rbp_bands = [[4 7];[8 12];[13 30];[31 80];[81 150]];
rbp_wholeband = [1 150];

% mA choice (1 or 3 mA) - change both variables
mA_selection_keys = ["3mA"];


% Frequency vector of interst
psd_freqs = 1:150;
w_all = [300]; % Window lengths to calculate 

for w_iter = 1:length(w_all)
    w_length = w_all(w_iter);
    % X ms sliding window, start at 5 ms
    % w_length = 300; % ms
    w_stride = 10; % ms
    w_start = 5; % ms
    w_end = w_length + 6; % ms %#### CURRENTLY will only do one window
    
    end_reached = 0;
    windows = [[w_start w_start+w_length]];
    iter = 1;
    while ~end_reached
        iter = iter + 1;
        if windows(iter-1,2) + w_stride <= w_end
            windows = [windows;[windows(iter-1,1) + w_stride, windows(iter-1,2) + w_stride]];
        else
            end_reached = 1;
        end
    end
    
    
    
    
    %% Create pats_spes
    
    % Read in the distance table
    % Load in the bipole pair euclidean distance table
    Table = readtable(dist_file);
    Table_Cell = table2cell(Table);
    table_dist_patIDs = string(Table_Cell(:,1));
    table_dist_from_bip = string(Table_Cell(:,2));
    table_dist_from_bip = replace(table_dist_from_bip," ","");
    table_dist_from_bip = replace(table_dist_from_bip,"-","");
    table_dist_to_bip = string(Table_Cell(:,3));
    table_dist_to_bip = replace(table_dist_to_bip," ","");
    table_dist_to_bip = replace(table_dist_to_bip,"-","");
    table_dist_mm = str2double(string(Table_Cell(:,4)));
    
    % Read in the SOZ label table
    % Read in the SOZ table
    Table = readtable(soz_label_file);
    Table_Cell = table2cell(Table);
    table_soz_patIDs = string(Table_Cell(:,1));
    table_soz_bip = string(Table_Cell(:,2));
    table_soz_bip = replace(table_soz_bip," ","");
    table_soz_bip = replace(table_soz_bip,"-","");
    table_soz_desig = string(Table_Cell(:,3));
    
    
    % Iterate through patients
    for i = 1:length(pat_dir_list)% 1 to the amount of selected patients
        
        % Get the channels_used file read in for this patient
        resp_used = string(table2cell(readtable(pat_chUsed_files(i),'ReadVariableNames',false)));
        resp_used = replace(resp_used," ","");
        resp_used = replace(resp_used,"-","");
        num_response = length(resp_used);
        
        % Get the distance subset for this patient
        idxs = find(contains(table_dist_patIDs(:,1), pat_dir_list{i}));
        dist_subset_from_bip = table_dist_from_bip(idxs);
        dist_subset_to_bip = table_dist_to_bip(idxs);
        dist_subset_mm = table_dist_mm(idxs);
        
        
        % Get the SOZ subset for this patient
        idxs = find(contains(table_soz_patIDs(:,1), pat_dir_list{i}));
        soz_subset_bip = table_soz_bip(idxs);
        soz_subset_desig = table_soz_desig(idxs);
        
        
        % Get pulse directory contents
        curr_dir = fullfile(source_dir, pat_dir_list(i),sprintf(pulse_folder_subpath,pat_dir_list(i)));
        contents = dir(curr_dir);
        mat_idxs = find(contains({contents.name}, mat_keyword));
        if isempty(mat_idxs); error("%d .MATs found in pulse directory (expected > 0)\n%s",length(mat_idxs), curr_dir); end
        fprintf("%d MATs found in pulse directory: %s\n",length(mat_idxs), curr_dir)
        mat_names_pulse = string({contents(mat_idxs).name})';
        pulse_full_mat_paths = fullfile(curr_dir,mat_names_pulse)';
        
        % Get pretrain directory contents
        curr_dir = fullfile(source_dir, pat_dir_list(i),sprintf(pretrain_folder_subpath,pat_dir_list(i)));
        contents = dir(curr_dir);
        mat_idxs = find(contains({contents.name}, mat_keyword));
        if isempty(mat_idxs); error("%d .MATs found in pretrain directory (expected > 0)\n%s",length(mat_idxs), curr_dir); end
        fprintf("%d MATs found in pretrain directory: %s\n",length(mat_idxs), curr_dir)
        mat_names_pretrain = string({contents(mat_idxs).name})';
        pretrain_full_mat_paths = fullfile(curr_dir,mat_names_pretrain)';
        
        for mA_idx = 1:length(mA_selection_keys)
            
            % Initialize variables
            pat_spes = struct;
            
            % Save for future use
            pat_spes.patID_clean = pat_dir_list{i};
            pat_spes.currents = mA_selection_keys(mA_idx);
            pat_spes.time_win_ms = windows;
            pat_spes.dist_threshes_mm = dist_threshes;
            pat_spes.psd_freqs = psd_freqs;
            pat_spes.data_axis_labels_for_each_current = ["Stim"; "Dist Thresh"; "Time Window"; "Response Bipolar Pair"; "PSD Freq"];
            
            % Find the single pulse files for this mA & organize into groups (i.e.
            % stim labels)
            mat_names_at_current = mat_names_pulse(contains(mat_names_pulse,mA_selection_keys(mA_idx)));
            mat_paths_at_current = pulse_full_mat_paths(contains(mat_names_pulse,mA_selection_keys(mA_idx)));
            
            
            % Get the list of unique stim pairs
            splits_pulse = split(mat_names_at_current,"_");
            stim_bips_dups = splits_pulse(:,2);
            stims_unique = unique(stim_bips_dups);
            stims_unique = replace(stims_unique," ","");
            stims_unique = replace(stims_unique,"-","");
            
            % Find the stims that are actually used (i.e. in gray matter after
            % preprocessing)
            stim_idxs_used = [];
            for j = 1:length(stims_unique)
                idx = find(contains(resp_used,stims_unique(j)));
                if length(idx) ~= 1
                    % fprintf("Info: Found %d matches in 'Channels_Used.xlsx' for stim channel %s in %s\n",length(idx), stims_unique(j), pats_spes(i).patID_clean);
                else
                    stim_idxs_used = [stim_idxs_used;j];
                end
            end
            fprintf("Summary: Could not find %d out of %d stim channels for %s in Channels_Used.xlsx\n", length(stims_unique)-length(stim_idxs_used),length(stims_unique),pat_spes.patID_clean);
            stims_used = stims_unique(stim_idxs_used);
            num_stims_used = length(stims_used);
            
            
            % Find the stim and response SOZ labels
            soz_stim = nan(num_stims_used,1);
            soz_resp = nan(num_response,1);
            for s = 1:num_stims_used
                idx = [];
                idx = find(contains(soz_subset_bip,stims_used(s)));
                if length(idx) ~= 1; error("ERROR: Found %d SOZ designation matches for stim channel %s in pat %s",length(idx), stims_used(s), pat_spes.patID_clean); end
                soz_stim(s) = soz_subset_desig(idx);
            end
            for r = 1:num_response
                idx = [];
                idx = find(contains(soz_subset_bip,resp_used(r)));
                if length(idx) ~= 1; error("ERROR: Found %d SOZ designation matches for resp channel %s in pat %s",length(idx), resp_used(r), pat_spes.patID_clean); end
                soz_resp(r) = soz_subset_desig(idx);
            end
            
            % Store for future use
            pat_spes.stim_labels = stims_used;
            pat_spes.response_labels = resp_used;
            pat_spes.stim_soz = soz_stim;
            pat_spes.response_soz = soz_resp;
            pat_spes.rbp_bands = rbp_bands;
            
            % Initialize this patient's 5D matrix (stim x distThresh x time x response x freq)'
            data_PSD_SUB = nan(num_stims_used,length(dist_threshes),size(windows,1),num_response,length(psd_freqs),'single');
            data_PSD_Z = nan(num_stims_used,length(dist_threshes),size(windows,1),num_response,length(psd_freqs),'single');
            data_RBP_SUB = nan(num_stims_used,length(dist_threshes),size(windows,1),num_response,size(rbp_bands,1),'single');
            data_RBP_Z = nan(num_stims_used,length(dist_threshes),size(windows,1),num_response,size(rbp_bands,1),'single');
            data_PSD_SUB_AllDists = nan(num_stims_used,1,size(windows,1),num_response,length(psd_freqs),'single');
            data_PSD_Z_AllDists = nan(num_stims_used,1,size(windows,1),num_response,length(psd_freqs),'single');
            data_RBP_Z_AllDists = nan(num_stims_used,1,size(windows,1),num_response,size(rbp_bands,1),'single');
            data_RBP_SUB_AllDists = nan(num_stims_used,1,size(windows,1),num_response,size(rbp_bands,1),'single');
            
            % Initialize pretrain_PSD to store later if needed
            pretrain_PSD_mean = nan(num_stims_used,num_response,length(psd_freqs),'single');
            pretrain_PSD_std = nan(num_stims_used,num_response,length(psd_freqs),'single');
            
            pretrain_RBP_mean = nan(num_stims_used,num_response,size(rbp_bands,1),'single');
            pretrain_RBP_std = nan(num_stims_used,num_response,size(rbp_bands,1),'single');
            
            % Iterate through stims used and pull out pulse files
            pulse_paths = [];
            
%             % Create pool here to avoid memory leak problems with parfor
%             poolobj = parpool('local',num_pool_workers);
%             fprintf("Parallel pool with %d workers started\n",poolobj.NumWorkers)
            
            parfor s = 1:num_stims_used
                
                fprintf("w_length %d/%d (%dms): Pat %d/%d (%s): mA %d/%d (%s): stim %d/%d (%s)\n", ...
                    w_iter,length(w_all),w_all(w_iter), i,length(pat_dir_list),...
                    pat_spes.patID_clean, mA_idx, length(mA_selection_keys), mA_selection_keys(mA_idx), s, num_stims_used,stims_used(s))
                
                % Find the pulse full path names
                file_resp_bips = replace(splits_pulse(:,2)," ","");
                file_resp_bips = replace(file_resp_bips,"-","");
                pulse_paths = mat_paths_at_current(find(contains(file_resp_bips,stims_used(s))));
                
                % Find pre-train files for this stim train (& at this mA)
                splits_pretrain = split(mat_names_pretrain,"_");
                files_pretrain_bips = replace(splits_pretrain(:,2)," ","");
                files_pretrain_bips = replace(files_pretrain_bips,"-","");
                pretrain_path = pretrain_full_mat_paths(find(contains(files_pretrain_bips,stims_used(s)) & contains(splits_pretrain(:,5),mA_selection_keys(mA_idx))));
                
                % Check if more than one pretrain file was found
                if length(pretrain_path) ~= 1
                    error("ERROR: found %d pretrain .mats for %s in pat %s, expected exactly 1",length(pretrain_path),stims_used(s),pat_spes.patID_clean);
                end
                
                
                % Calculate pre-train values for the channels used (with same
                % window length and stride as used for post-stim)
                myVars = {'pre_train_1','pre_train_2','pre_train_3','fs','labels'};
                S = load(pretrain_path,myVars{:}); % There should be 3 pretrain segments (pre_train_1/2/3) in addition to 'fs' and 'labels'
                pre_train_1 = S.pre_train_1;
                pre_train_2 = S.pre_train_2;
                pre_train_3 = S.pre_train_3;
                fs = S.fs;
                labels = S.labels;
                
                labels_str = string(labels)';
                labels_str = replace(labels_str," ","");
                labels_str = replace(labels_str,"-","");
                
                % 'Make bipolar pairs of the pulse file lables (a "Naive"
                % scheme can be used with overlap between electrodes because
                % we will just search for the used names anyway
                % Also FILTER the data
                pretrain_bip_labels = strcat(labels_str(1:end-1), labels_str(2:end));
                pretrain_bip_data_1 = FT_filt_LineFrqNotch(pre_train_1(1:end-1,:) - pre_train_1(2:end,:),fs);
                pretrain_bip_data_2 = FT_filt_LineFrqNotch(pre_train_2(1:end-1,:) - pre_train_2(2:end,:),fs);
                pretrain_bip_data_3 = FT_filt_LineFrqNotch(pre_train_3(1:end-1,:) - pre_train_3(2:end,:),fs);
                
                pretrain_bip_data_1_respOrg = nan(num_response,size(pretrain_bip_data_1,2));
                pretrain_bip_data_2_respOrg = nan(num_response,size(pretrain_bip_data_2,2));
                pretrain_bip_data_3_respOrg = nan(num_response,size(pretrain_bip_data_3,2));
                
                % Pull out and organize only the data that corresponds to response channels
                for k = 1:num_response
                    idx = [];
                    idx = find(contains(pretrain_bip_labels,resp_used(k)));
                    if length(idx) ~= 1; error("ERROR: Found %d label matches for resp channel %s in pat %s",length(idx), resp_used(k), pat_spes.patID_clean); end
                    pretrain_bip_data_1_respOrg(k,:) = pretrain_bip_data_1(idx,:);
                    pretrain_bip_data_2_respOrg(k,:) = pretrain_bip_data_2(idx,:);
                    pretrain_bip_data_3_respOrg(k,:) = pretrain_bip_data_3(idx,:);
                end
                
                % Stop here if just checking for the presence of all the files
                % and correct label/SOZ matching for each patient
                if ~only_run_file_check
                    
                    % Assuming all three epochs are the same length, calculate
                    % windows for pretrain segments
                    pretrain_length_ms = (size(pretrain_bip_data_1,2) / fs) * 1000;
                    end_reached = 0;
                    pt_windows = [[0 w_length]];
                    iter = 1;
                    while ~end_reached
                        iter = iter + 1;
                        if pt_windows(iter-1,2) + w_stride <= pretrain_length_ms
                            pt_windows = [pt_windows;[pt_windows(iter-1,1) + w_stride, pt_windows(iter-1,2) + w_stride]];
                        else
                            end_reached = 1;
                        end
                    end
                    
                    % Calculate the average PSD for each response channel across the three pretrain epochs
                    pretrain_PSD_perW = nan(3 * size(pt_windows,1),num_response,length(psd_freqs));
                    
                    % In paralell we will be computing relative band power
                    % (RBP), so we need to calculate the pretrain RBP for
                    % the frequency bands of interest. Initialize pretrain
                    % vairable here.
                    pretrain_RBP_perW = nan(3 * size(pt_windows,1),num_response,size(rbp_bands,1));
                    
                    
                    for k = 1:size(pt_windows,1)
                        samp_start = round(pt_windows(k,1)/1000 * fs + 1);
                        samp_end = round(pt_windows(k,2)/1000 * fs);
                        
                        WINDOW = zeros(1,size(pretrain_bip_data_1_respOrg,2));
                        WINDOW(samp_start:samp_end) = 1;
                        
                        % Each column is a signal for periodogram (i.e. must flip)
                        pxx1 = pwelch(pretrain_bip_data_1_respOrg',WINDOW,[],psd_freqs,fs);
                        pxx2 = pwelch(pretrain_bip_data_2_respOrg',WINDOW,[],psd_freqs,fs);
                        pxx3 = pwelch(pretrain_bip_data_3_respOrg',WINDOW,[],psd_freqs,fs);
                        
                        % Assign PSDs to a single variable
                        pretrain_PSD_perW((k-1)*3+1,:,:) = pxx1';
                        pretrain_PSD_perW((k-1)*3+2,:,:) = pxx2';
                        pretrain_PSD_perW((k-1)*3+3,:,:) = pxx3';
                        
                        % Calculate relative band power (RBP) for this
                        % window
                        epoch_pt1 = pretrain_bip_data_1_respOrg(:,samp_start:samp_end);
                        epoch_pt2 = pretrain_bip_data_2_respOrg(:,samp_start:samp_end);
                        epoch_pt3 = pretrain_bip_data_3_respOrg(:,samp_start:samp_end);
                        
                        total_power_pt1 = bandpower(epoch_pt1',fs,rbp_wholeband)';
                        total_power_pt2 = bandpower(epoch_pt2',fs,rbp_wholeband)';
                        total_power_pt3 = bandpower(epoch_pt3',fs,rbp_wholeband)';
                        
                        % Iterate through the bands desired
                        for m = 1:size(rbp_bands,1)
                            pretrain_RBP_perW((k-1)*3+1,:,m) = bandpower(epoch_pt1',fs,rbp_bands(m,:))'./total_power_pt1;
                            pretrain_RBP_perW((k-1)*3+2,:,m) = bandpower(epoch_pt2',fs,rbp_bands(m,:))'./total_power_pt2;
                            pretrain_RBP_perW((k-1)*3+3,:,m) = bandpower(epoch_pt3',fs,rbp_bands(m,:))'./total_power_pt3;
                        end
                        
                    end
                    
                    % Average the pretrain metrics
                    pretrain_PSD_mean(s,:,:) = squeeze(nanmean(pretrain_PSD_perW,1));
                    pretrain_PSD_std(s,:,:) = squeeze(nanstd(pretrain_PSD_perW,[],1));
                    
                    pretrain_RBP_mean(s,:,:) = squeeze(nanmean(pretrain_RBP_perW,1));
                    pretrain_RBP_std(s,:,:) = squeeze(nanstd(pretrain_RBP_perW,[],1));
                    
                    % Get a subset of the subset from the distance table based
                    % on this stim bip
                    idxs = find(contains(dist_subset_from_bip,stims_used(s)));
                    dist_subset_to_subset_bip = dist_subset_to_bip(idxs);
                    dist_subset_mm_subset = dist_subset_mm(idxs);
                    
                    % Organize the distances based on the response channel
                    % order
                    dist_mm_ordered = nan(num_response,1);
                    for res = 1:num_response
                        idx = find(contains(dist_subset_to_subset_bip,resp_used(res)));
                        if length(idx)~=1; error("ERROR: %s not found in distance subset",resp_used(res)); end
                        dist_mm_ordered(res) = dist_subset_mm_subset(idx);
                    end
                    
                    % Initialize temp variable
                    dist_tw_SUB_pxx = nan(length(pulse_paths),length(dist_threshes),size(windows,1),num_response,length(psd_freqs));
                    dist_tw_ZSCORED_pxx = nan(length(pulse_paths),length(dist_threshes),size(windows,1),num_response,length(psd_freqs));
                    dist_tw_ZSCORED_RBP = nan(length(pulse_paths),length(dist_threshes),size(windows,1),num_response,size(rbp_bands,1));
                    dist_tw_SUB_RBP = nan(length(pulse_paths),length(dist_threshes),size(windows,1),num_response,size(rbp_bands,1));
                    dist_tw_SUB_pxx_AllDists = nan(length(pulse_paths),1,size(windows,1),num_response,length(psd_freqs));
                    dist_tw_ZSCORED_pxx_AllDists = nan(length(pulse_paths),1,size(windows,1),num_response,length(psd_freqs));
                    dist_tw_ZSCORED_RBP_AllDists = nan(length(pulse_paths),1,size(windows,1),num_response,size(rbp_bands,1));
                    dist_tw_SUB_RBP_AllDists = nan(length(pulse_paths),1,size(windows,1),num_response,size(rbp_bands,1));
                    
                    % Iterate through the pulses
                    for p = 1:length(pulse_paths)
                        
                        S = load(pulse_paths(p));
                        fs_pulse = S.fs;
                        if fs_pulse ~= fs; error("Sampling frequencies do not match"); end
                        labels_pulse = string(S.labels)';
                        %'
                        pulse_data = S.pulse;
                        labels_pulse = replace(labels_pulse," ","");
                        labels_pulse = replace(labels_pulse,"-","");
                        
                        % Make bipolar pairs of the pulse file lables (a "Naive"
                        % scheme can be used with overlap between electrodes because
                        % we will just search for the used names anyway)
                        % Also FILTER the data
                        pulse_bip_labels = strcat(labels_str(1:end-1), labels_str(2:end));
                        pulse_bip_data = FT_filt_LineFrqNotch(pulse_data(1:end-1,:) - pulse_data(2:end,:),fs);
                        pulse_bip_data_respOrg = nan(num_response,size(pulse_bip_data,2));
                        
                        % Pull out and organize only the data that corresponds to response channels
                        for k = 1:num_response
                            idx = [];
                            idx = find(contains(pulse_bip_labels,resp_used(k)));
                            if length(idx) ~= 1; error("ERROR: Found %d label matches for resp channel %s in pat %s",length(idx), resp_used(k), pat_spes.patID_clean); end
                            pulse_bip_data_respOrg(k,:) = pulse_bip_data(idx,:);
                        end
                        
                        
                        % Calculate PSD for chans_used_thresholded_data WITHIN
                        % the time window
                        for k = 1:size(windows,1)
                            samp_start = round(windows(k,1)/1000 * fs + 1);
                            samp_end = round(windows(k,2)/1000 * fs);
                            
                            WINDOW = zeros(1,size(pulse_bip_data_respOrg,2));
                            WINDOW(samp_start:samp_end) = 1;
                            
                            % Each column is a signal for periodogram (i.e. must flip)
                            pxx = pwelch(pulse_bip_data_respOrg',WINDOW,[],psd_freqs,fs)';
                            
                            % Calculate pulse RBP for desired bands
                            rbp = nan(num_response,size(rbp_bands,1));
                            rbp_pulse_wb = bandpower(pulse_bip_data_respOrg(:,samp_start:samp_end)',fs,rbp_wholeband)';
                            for m = 1:size(rbp_bands,1)
                                rbp(:,m) = bandpower(pulse_bip_data_respOrg(:,samp_start:samp_end)',fs,rbp_bands(m,:))'./rbp_pulse_wb;
                            end
                                
                            % Create an all distance metric (to be used for
                            % z-scoring)
                            dist_tw_SUB_pxx_AllDists(p,1,k,:,:) = pxx - squeeze(pretrain_PSD_mean(s,:,:));
                                
                            % Z-score each pulse's PSD to pre-train'
                            dist_tw_ZSCORED_pxx_AllDists(p,1,k,:,:) = (pxx - squeeze(pretrain_PSD_mean(s,:,:)))./squeeze(pretrain_PSD_std(s,:,:));
                            
                            % RBP for all distances
                            dist_tw_ZSCORED_RBP_AllDists(p,1,k,:,:) = (rbp - squeeze(pretrain_RBP_mean(s,:,:)))./squeeze(pretrain_RBP_std(s,:,:));
                            dist_tw_SUB_RBP_AllDists(p,1,k,:,:) = (rbp - squeeze(pretrain_RBP_mean(s,:,:)))./squeeze(pretrain_RBP_std(s,:,:));
                            
                            % Iterate through the distance thresholds
                            for d = 1:length(dist_threshes)
                                
                                % Re-assign the raw values so that
                                % previoulsy excluded values are returned
                                % into data
                                pxx_fresh = pxx;
                                rbp_fresh = rbp;
                                
                                % Find the response channels outside of the distance band
                                % Ensuring > 0 is always true will eliminate self channel
                                exclude_idxs = find(dist_mm_ordered <= dist_threshes(d) | dist_mm_ordered > dist_threshes(d) + dist_gap);
                                
                                % Only select the channels outside the distance
                                % threshold, NAN the other channel data
                                pxx_fresh(exclude_idxs,:) = nan;
                                rbp_fresh(exclude_idxs,:) = nan;
                                
                                % subtract the pretrain PSD from the
                                % post-stimulation PSD
                                dist_tw_SUB_pxx(p,d,k,:,:) = pxx_fresh - squeeze(pretrain_PSD_mean(s,:,:));
                                
                                % Z-score each pulses PSD to pre-train
                                dist_tw_ZSCORED_pxx(p,d,k,:,:) = (pxx_fresh - squeeze(pretrain_PSD_mean(s,:,:)))./squeeze(pretrain_PSD_std(s,:,:));
                                
                                % Subtract the pretrain RBP for each band
                                dist_tw_SUB_RBP(p,d,k,:,:) = (rbp_fresh - squeeze(pretrain_RBP_mean(s,:,:)));
                                
                                % Z-score the RBP for each band to pre-train
                                dist_tw_ZSCORED_RBP(p,d,k,:,:) = (rbp_fresh - squeeze(pretrain_RBP_mean(s,:,:)))./squeeze(pretrain_RBP_std(s,:,:));
                                
                            end
                        end
                    end
                    % Average the Z-scored PSD across the pulses for this stim pair (includes all distances and tws) and store in pats 5D matrix
                    % 5D matrix (stim x distThresh x time x response x freq)
                    data_PSD_SUB(s,:,:,:,:) = squeeze(nanmean(dist_tw_SUB_pxx,1));
                    data_PSD_Z(s,:,:,:,:) = squeeze(nanmean(dist_tw_ZSCORED_pxx,1));
                    data_PSD_SUB_AllDists(s,:,:,:,:) = squeeze(nanmean(dist_tw_SUB_pxx_AllDists,1));
                    data_PSD_Z_AllDists(s,:,:,:,:) = squeeze(nanmean(dist_tw_ZSCORED_pxx_AllDists,1));
                    
                    data_RBP_Z(s,:,:,:,:) = squeeze(nanmean(dist_tw_ZSCORED_RBP,1));
                    data_RBP_SUB(s,:,:,:,:) = squeeze(nanmean(dist_tw_SUB_RBP,1));
                    data_RBP_Z_AllDists(s,:,:,:,:) = squeeze(nanmean(dist_tw_ZSCORED_RBP_AllDists,1));
                    data_RBP_SUB_AllDists(s,:,:,:,:) = squeeze(nanmean(dist_tw_SUB_RBP_AllDists,1));
                    
                end
            end
            
            % Save to the output struct (must be outside parfor loop)
            pat_spes.mean_stim_pretrain_PSD = pretrain_PSD_mean;
            pat_spes.std_stim_pretrain_PSD = pretrain_PSD_std;
            pat_spes.data_SUB = data_PSD_SUB;
            pat_spes.data_Z = data_PSD_Z;
            pat_spes.data_RBP_Z = data_RBP_Z;
            pat_spes.data_RBP_SUB = data_RBP_SUB;
            pat_spes.data_RBP_Z_AllDists = data_RBP_Z_AllDists;
            pat_spes.data_RBP_SUB_AllDists = data_RBP_SUB_AllDists;
            pat_spes.data_SUB_AllDists = data_PSD_SUB_AllDists;
            pat_spes.data_Z_AllDists = data_PSD_Z_AllDists;
            
            
            % Save the struct, delete it, and shutdown pool (restarted on next loop) to avoid running out of memory
            save(sprintf('%s_%s_SPES_windowed_PSD_%d_length_%d_stride_DISTGAP%dmm.mat',pat_spes.patID_clean, mA_selection_keys(mA_idx),w_length,w_stride,dist_gap),'pat_spes','dist_gap')
            delete pat_spes
%             delete(poolobj)
            toc
        end
        toc
    end
    toc
end
toc