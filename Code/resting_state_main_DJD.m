
%% Connectivity calculation for 2-minute epochs (epochs already generated)

% Hard coded to do Theta[4-8 Hz], Alpha[8-12 Hz], Beta[13-30 Hz], Low
% Gamma[31-80 Hz], High Gamma[81-150 Hz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated on 2/20/2024 by DJD to include delta


close all
clear all
restoredefaultpath

source_dir = 'Z:\000_Data\SEEG\SEEG_EyesClosed_RestingState\data';
sub_folder_path = 'Unfiltered_Chunks\First_Collection\All_2minChunks_Bipole_UnusedChannelsDeleted_UnFiltered';
edf_keyword = 'Raw_bipole_unusedChannelsDeleted_Unfiltered';
soz_label_file = 'Z:\000_Data\SEEG\SEEG_EyesClosed_RestingState\labels\all_pats_bipole_soz_labels.csv';
metric_loc_in_struct = 3;

%Path to violin plot folder
addpath('Z:\shared_toolboxes\Violinplot-Matlab-master')
% Path to EDFread
addpath('Z:\shared_toolboxes\edfRead')
% Path to fieldtrip toolbox
addpath('Z:\shared_toolboxes\fieldtrip-20190819');
% Path to RolstonCode for Imaginary COherence
addpath(genpath('Z:\shared_toolboxes\RolstonCode'));

% Select which connectivty matrix to generate (and associated variables)

metric = 'PDC';
mvar_order = 10; % Will need to do 1 for short epochs
long_epoch_seg_length = 5; % seconds in which to define new trials

% metric = 'DTF';
% mvar_order = 10; % Will need to do 1 for short epochs
% long_epoch_seg_length = 5; % seconds in which to define new trials

% metric = 'GC';

% metric = 'PLI';

% metric = 'CFD'; ??? cross frequency one... https://www.medrxiv.org/content/10.1101/2021.12.30.21268524v1.full

% metric = 'ImCoh';

% Using a list of patients (full folder name listed under '000_data\SEEG\data, index into SEEG data folder and pull out the

% Select the patients we wish to get PDC data for

pat_dir_list = ["Epat26_dpat05", "Epat27_dpat06", "Epat30",...
    "Epat31_apat111", "Epat34_dpat11", "Epat35_apat115", "Epat37_dpat13",...
    "Epat38", "Epat39", "Spat28", "Spat30",...
    "Spat31", "Spat34", "Spat36", "Spat37",...
    "Epat43_apat117_vpat04", "Spat41", "Spat42", "Spat44",...
    "Spat39", "Spat40", "Epat02", "Epat03",...
    "Epat04_dpat07", "Epat05", "Epat06_apat102", "Epat08_dpat08",...
    "Epat09", "Epat10", "Epat11_apat113", "Epat13",...
    "Epat14", "Epat15", "Epat17_apat112", "Epat18_dpat03",...
    "Epat19_dpat04", "Epat20", "Epat21_pat107",...
    "Epat22", "Epat23_spat21_dpat02", "Epat24",...
    "Epat25_apat108_dpat10", "Epat28_apat106", "pat11", "pat33_dpat09",...
    "Spat12", "Spat13", "Spat14", "Spat17",...
    "Spat18", "Spat19", "Spat20", "Spat22",...
    "Spat23", "Spat24", "Spat26", "Spat29",...
    "Spat32", "Spat33", "Spat48", "Spat49",...
    "Spat50", "Spat51", "Spat52", "Spat53",...
    "Spat47", "Epat33", "Spat02", "Spat03",...
    "Spat05", "Spat06", "Spat07", "Spat08",...
    "Spat09", "Spat10", "Spat11", "Spat25",...,
    "Spat27", "Spat45", "Spat46", "Epat41"];



% name in the folder supercedes the name in excel


% Initialize the output
% Save into .mat
pats = struct;


% Iterate through the subject directories
for i = 1:length(pat_dir_list)% 1 to the amount of selected patients
    
    curr_dir = fullfile(source_dir, pat_dir_list(i),sub_folder_path);
    
    pats(i).subID = pat_dir_list(i);
    
    % Get directory contents
    contents = dir(curr_dir);
    edf_idxs = find(contains({contents.name}, edf_keyword));
    if isempty(edf_idxs); error("%d EDFs found in directory (expected > 0)\n%s",length(edf_idxs), curr_dir); end
    fprintf("\n%d EDFs found in directory:\n",length(edf_idxs))
    edf_names = string({contents(edf_idxs).name}')
    full_edf_paths = fullfile(curr_dir,edf_names);
    
    data = cell(length(edf_idxs),1);
    hdrs = cell(length(edf_idxs),1);
    
    % Import all of the EDFs
    for j = 1:length(edf_idxs)
        disp("Importing EDF...")
        [hdrs{j}, data{j}] = edfread(full_edf_paths(j));
        disp("EDF import complete")
    end
    
    %% Obtain connectivity metric for all 2-minute epochs for this subject
    
    % Initialize the 2-minute epoch results
    all_D = nan(length(edf_idxs),length(hdrs{1}.label),length(hdrs{1}.label));
    all_T = all_D;
    all_A = all_D;
    all_B = all_D;
    all_G_low = all_D;
    all_G_high = all_D;
    
    for j = 1:length(edf_idxs)
        
        fprintf("%d/%d Sub ID: %s, epoch %d\n",i,length(pat_dir_list),pat_dir_list(i),j)
        
        if j ==1
            % Save a copy of the labels for this patient
            pats(i).labels = string(hdrs{j}.label)';
        end
        
        % Filter the data
        data_curr = data{j};
        filt_data = FT_filt_A(data_curr,hdrs{j}.frequency(1));
        
        % Calculate the metric of choice - all hard coded for the following
        % bands: Theta[4-8 Hz], Alpha[8-12 Hz], Beta[13-30 Hz], Low Gamma[31-80 Hz], High Gamma[81-150 Hz]
        
        switch metric
            
            case 'PDC'
                [D,T,A,B,G_low,G_high] = pdc_calc(filt_data,hdrs{j},mvar_order,long_epoch_seg_length);
                
            case 'DTF'
                [T,A,B,G_low,G_high] = dtf_calc(filt_data,hdrs{j},mvar_order,long_epoch_seg_length);
                
            case 'ImCoh'
                [T,A,B,G_low,G_high] = imcoh_calc(filt_data,hdrs{j});
                
        end
        
        all_D(j,:,:) = D;
        all_T(j,:,:) = T;
        all_A(j,:,:) = A;
        all_B(j,:,:) = B;
        all_G_low(j,:,:) = G_low;
        all_G_high(j,:,:) = G_high;
        
    end
    
    % Average the PDC across 2-minute epochs
    avg_D = squeeze(nanmean(all_D,1));
    avg_T = squeeze(nanmean(all_T,1));
    avg_A = squeeze(nanmean(all_A,1));
    avg_B = squeeze(nanmean(all_B,1));
    avg_G_low = squeeze(nanmean(all_G_low,1));
    avg_G_high = squeeze(nanmean(all_G_high,1));
    
    % Assign to master patient variable
    pats(i).long.delta = avg_D;
    pats(i).long.theta = avg_T;
    pats(i).long.alpha = avg_A;
    pats(i).long.beta = avg_B;
    pats(i).long.gamma_low = avg_G_low;
    pats(i).long.gamma_high = avg_G_high;
    
    
end


%% Go through the patient connectivity matrices and assign SOZ labels

% Read in the SOZ table
Table = readtable(soz_label_file);
Table_Cell = table2cell(Table);

%Standardize Patient names based on the spreadsheet (i.e., E/Spat + a number)
for i = 1:length(pat_dir_list)
    if (pat_dir_list{i}(1) == 'p') || (pat_dir_list{i}(1) == 'P') % file is called pat
        idx_p = strfind(pat_dir_list{i}, '_');
        
        % Check if there is no second name (i.e. no "_" present)
        if isempty(idx_p)
            pat_dir_list{i} = pat_dir_list{i};
        else
            pat_dir_list{i} = pat_dir_list{i}(1:idx_p-1);
        end
        
    else % file is called epat, spat, dpat, etc.
        idx_np = strfind(pat_dir_list{i}, '_');
        
        % Check if there is no second name (i.e. no "_" present)
        if isempty(idx_np)
            pat_dir_list{i} = pat_dir_list{i};
        else
            pat_dir_list{i} =  pat_dir_list{i}(1:idx_np-1);
        end
    end
    
    % Save for future use
    pats(i).patID_clean = pat_dir_list{i};
end

%Standardize Bipole based on the spreadsheet (i.e., E/Spat + a number);
%Eliminate dashes and white spaces
for i = 1:length(Table_Cell)
    Table_Cell{i,2} = replace(Table_Cell{i,2},"-","");
    Table_Cell{i,2} = replace(Table_Cell{i,2}," ","");
end

% Make a new field for SOZ data
% Now, the patient names match up with the names in "Table"
% We can now find the SoZ designations and pair them to the data in pats
% Matching the SoZ designations for each patient

% Create a new field for the SOZ data in Pats
for FU = 1:size(pats,2)
    pats(FU).SOZ = [];
end

% Matching algorithm
for ii = 1:size(pat_dir_list,2) % 1: number of patients
    for jj = 1:size(Table_Cell,1) % 1: all contact pairs across all patients
        if Table_Cell(jj,1) == pat_dir_list(ii) % there is a match b/w the patient name in excel and the patient in Matlab...
            bipole_pair = Table_Cell(jj,2); % Then assign "bipole pair" for that patient to that category
            for kk = 1:size(pats(ii).labels,1) % Now, loop through each patient's number of contact pairs
                if pats(ii).labels(kk) == cell2mat(Table_Cell(jj,2)) % if the patient's Matlab contact pair is identical to the excel contact pair...
                    pats(ii).SOZ(kk) = cell2mat(Table_Cell(jj,3)); % Then, their SOZ designation will be assigned to the designation on the Excel sheet
                end
            end
        end
    end
    % Transpose the data so it is a column vector
    pats(ii).SOZ = [pats(ii).SOZ]';
end

%% Z-Scoring and Averaging Algorithm
% Make the struct to store the data
for n = 1:size(pat_dir_list,2) % loops through the number of patients
    pats(n).Avg_Soz = struct(); % each patient is assigned a struct for average SoZ
    pats(n).Avg_Soz.bands = []; % Empty array to store the values, for each band
end

% Initizalize ROC curve master variables
roc_bip_in = cell(1,6);
roc_bip_out = cell(1,6);
roc_soz_label = cell(1,6);

for a = 1:size(pat_dir_list,2) % loops through the number of patients
    
    fprintf("Sub ID: %s\n",pat_dir_list(a))
    
    fields = fieldnames(pats(a).long); % lists freqiency bands
    pats_cell = struct2cell(pats(a)); % re-computes for every patient!
    long_cell = struct2cell(pats_cell{metric_loc_in_struct});
    
    % Z-score the matrices to themselves
    for z = 1:size(long_cell,1) % 1 to 5 long of these for each band for each pt.
        mean_mat = nanmean(nanmean(long_cell{z}));
        std_mat = nanstd(nanstd(long_cell{z}));
        long_cell{z} = (long_cell{z}-mean_mat)/std_mat;
    end
    
    % Save the Z-scored matrices in pats struct
    pats(a).long_Z = long_cell;
    
    num_bands = length(fields);
    str = string(fields); % puts fields into a string
    
    % For this patient, enter each frequency band and do averaging
    
    for b = 1:num_bands % loops through each frequency band
        
        % Save the ROC values
        roc_bip_in{b} = [roc_bip_in{b}; nanmean(long_cell{b},1)'];
        roc_bip_out{b} = [roc_bip_out{b}; nanmean(long_cell{b},2)];
        roc_soz_label{b} = [roc_soz_label{b}; pats(a).SOZ];
        
        
        avg_soz_in = [];
        avg_soz_out = [];
        avg_nsoz_in = [];
        avg_nsoz_out = [];
        
        for c = 1:size(pats(a).SOZ,1) % loops through the number of bipole pairs for that patient
            if pats(a).SOZ(c,1) == 1 % if SOZ in this position
                avg_soz_in = [avg_soz_in, mean(long_cell{b}(:,c), 'omitnan')];
                avg_soz_out = [avg_soz_out, mean(long_cell{b}(c,:), 'omitnan')];
            elseif pats(a).SOZ(c,1) ~= 1 % if not an SOZ in this position
                avg_nsoz_in = [avg_nsoz_in, mean(long_cell{b}(:,c), 'omitnan')];
                avg_nsoz_out = [avg_nsoz_out, mean(long_cell{b}(c,:), 'omitnan')];
            end
        end
        
        % now, store the average SOZ and non SOZ in and out for each band
        % Row is a band; column is a different entry, as dictated below:
        
        % each row is a different band; each column is a different entry
        pats(a).Avg_Soz.bands(b,1) = mean(avg_soz_in);
        pats(a).Avg_Soz.bands(b,2) = mean(avg_soz_out);
        pats(a).Avg_Soz.bands(b,3) = mean(avg_nsoz_in);
        pats(a).Avg_Soz.bands(b,4) = mean(avg_nsoz_out);
    end
end

% save to struct
save('Z:\000_Data\SEEG\SEEG_EyesClosed_RestingState\code\connectivity_calculations\PDC_pats_with_delta.mat', 'pats');


%% Sort the Data Into its Own Struct

% results = struct;
% for i = 1:size(pat_dir_list,2)
%     results(i).pat_id = pat_dir_list(i);
%     results(i).freqbandnames = fields;
%     results(i).avg_PDC = pats(i).Avg_Soz.bands;
%     results(i).dataorder = {'SozIn','SozOut', 'nSozIn','nSozOut'};
% end

% %% Violin Plots
% 
% cum_band_data = struct;
% 
% for i = 1:size(long_cell,1) % loop through the bands
%     cum_band_data(i).SOZ_in = [];
%     cum_band_data(i).SOZ_out = [];
%     cum_band_data(i).nSOZ_in = [];
%     cum_band_data(i).nSOZ_out = [];
%     
%     for j = 1:size(pats,2) % loop through the patients
%         cum_band_data(i).SOZ_in = [cum_band_data(i).SOZ_in, pats(j).Avg_Soz.bands(i,1)];
%         cum_band_data(i).SOZ_out = [cum_band_data(i).SOZ_out, pats(j).Avg_Soz.bands(i,2)];
%         cum_band_data(i).nSOZ_in = [cum_band_data(i).nSOZ_in, pats(j).Avg_Soz.bands(i,3)];
%         cum_band_data(i).nSOZ_out = [cum_band_data(i).nSOZ_out,pats(j).Avg_Soz.bands(i,4)];
%     end
%     
% end
% 
% %% Cello Plots
% 
% for i = 1:size(long_cell,1)  % for each frequency band
%     violin_struct = struct;
%     figure;
%     clf
%     str = sprintf('Inward %s: %s Frequency Band', metric, results(1).freqbandnames{i,1});
%     %title(str);
%     violin_struct.SoZ = cum_band_data(i).SOZ_in;
%     violin_struct.NonSoZ = cum_band_data(i).nSOZ_in;
%     Inward_Viola = violinplot(violin_struct);
%     Inward_Viola(1,1).ViolinColor = [1 0 0];
%     Inward_Viola(1,2).ViolinColor = [65/256,105/256,225/256];
%     ylim([-2, 6])
%     ylabel(sprintf('Z-scored %s Inward Strength',metric))
%     ax = gca;
%     ax.XGrid = 'off';
%     ax.YGrid = 'on';
%     yticks([-2 0 2 4 6])
%     
%     violin_struct = struct;
%     figure;
%     clf
%     str = sprintf('Outward %s: %s Frequency Band', metric, results(1).freqbandnames{i,1});
%     %title(str);
%     violin_struct.SoZ = cum_band_data(i).SOZ_out;
%     violin_struct.NonSoZ = cum_band_data(i).nSOZ_out;
%     Outward_Viola = violinplot(violin_struct);
%     Outward_Viola(1,1).ViolinColor = [1 0 0];
%     Outward_Viola(1,2).ViolinColor = [65/256,105/256,225/256];
%     ylim([-2, 6])
%     ylabel(sprintf('Z-scored %s Outward Strength',metric))
%     ax = gca;
%     ax.XGrid = 'off';
%     ax.YGrid = 'on';
%     yticks([-2 0 2 4 6])
%     
%     
%     
%     % Stats
%     [H,P_in,CI,STATS] = ttest(cum_band_data(i).SOZ_in, cum_band_data(i).nSOZ_in);
%     [H,P_out,CI,STATS] = ttest(cum_band_data(i).SOZ_out, cum_band_data(i).nSOZ_out);
%     fprintf("%s, Inward diff (%0.3f Z) p=%0.3e, Outward diff (%0.3f Z): p=%0.3e\n",...
%         results(1).freqbandnames{i,1},...
%         nanmean(cum_band_data(i).SOZ_in) - nanmean(cum_band_data(i).nSOZ_in),...
%         P_in,...
%         nanmean(cum_band_data(i).SOZ_out) - nanmean(cum_band_data(i).nSOZ_out),...
%         P_out)
%     
%     
% end
% 
% %% Plot ROC Curves
% 
% for i = 1:length(fields)
%     
%     roc_labels = roc_soz_label{i} == 1;
%     
%     figure
%     [X,Y,T,AUC_in,OPTROCPT] = perfcurve(roc_labels,roc_bip_in{i},1);
%     plot(X,Y)
%     xlabel('False positive rate')
%     ylabel('True positive rate')
%     title(sprintf('ROC for %s Inward (SOZ vs Non-SOZ): AUC %0.3f',fields{i},AUC_in))
%     
%     % ROC curve Outward
%     figure
%     [X,Y,T,AUC_out,OPTROCPT] = perfcurve(roc_labels,roc_bip_out{i},0);
%     plot(X,Y)
%     xlabel('False positive rate')
%     ylabel('True positive rate')
%     title(sprintf('ROC for %s Outward (SOZ vs Non-SOZ): AUC %0.3f',fields{i}, AUC_out))
%     
%     
%     % ROC curve Inward - Outward
%     figure
%     [X,Y,T,AUC_inOut,OPTROCPT] = perfcurve(roc_labels,roc_bip_in{i} - roc_bip_out{i},1);
%     plot(X,Y)
%     xlabel('False positive rate')
%     ylabel('True positive rate')
%     title(sprintf('ROC for %s Inward-Outward (SOZ vs Non-SOZ): AUC %0.3f',fields{i},AUC_inOut))
%     
% end
% 
% 
