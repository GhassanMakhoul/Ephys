clear all
close all
clc

restoredefaultpath

addpath Z:\shared_toolboxes\SEEG_sleep_staging_DJD_edits\
addpath Z:\shared_toolboxes\MatlabProgressBar\
addpath Z:\shared_toolboxes\Derek_functions\

%
sleep_dir = "Z:\000_Data\SEEG\SEEG_Sleep_Staging\data\nights_of_sleep\";
all_sleep_csv_fpath = "Z:\000_Data\SEEG\SEEG_Sleep_Staging\notes\all_sleep_stages_06122023.csv";
log_fpath = "Z:\000_Data\SEEG\SEEG_Sleep_Staging\notes\sleep_processing_log_06132023.txt";

cohort_fpath = "Z:\000_Data\SEEG\SEEG_Sleep_Staging\notes\cohort6.xlsx";

pat_ids = readtable(cohort_fpath,'ReadVariableNames',false); 
pat_ids = string(pat_ids{:, 1});

% read in the csv file
opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Type", "PatID", "POD", "SleepCat", "OnsetDatetime", "OffsetDatetime", "Duration"];
opts.VariableTypes = ["string", "string", "double", "string", "datetime", "datetime", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Type", "PatID", "SleepCat"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Type", "PatID", "SleepCat"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "OnsetDatetime", "InputFormat", "yyyy-MM-dd HH:mm:ss");
opts = setvaropts(opts, "OffsetDatetime", "InputFormat", "yyyy-MM-dd HH:mm:ss");

% Import the data
all_sleep_data_extracted = readtable(all_sleep_csv_fpath, opts);

% [~, pat_ids] = get_files_in_dir(sleep_dir, "*pat*", 1);

% skipped Epat14 as the files were too big
% skipped Epat15 b/c they have diff channels
for pp = 1:length(pat_ids)
    pat_id = pat_ids(pp);

    % get the PODs
    [~, tmp_sleep_edf_files] = get_files_in_dir(sleep_dir+pat_id+"\","*_filt.edf",0);
    if length(tmp_sleep_edf_files) == 1 & tmp_sleep_edf_files(1) == ""
        continue;
    end

    per_patient_PODs = [""];
    for qq = 1:length(tmp_sleep_edf_files)
        tmp_split_edf_file = strsplit(tmp_sleep_edf_files(qq), "_");
        tmp_split_edf_file = tmp_split_edf_file(2);
        tmp_split_POD = strsplit(tmp_split_edf_file, "-");
        tmp_split_POD = strrep(tmp_split_POD(1), "Night", "");
        per_patient_PODs(qq) = tmp_split_POD;
    end

    per_patient_PODs = unique(per_patient_PODs);


    per_patient_csv = all_sleep_data_extracted(all_sleep_data_extracted{:, 2} == pat_id, :);

    for qq = 1:length(per_patient_PODs)

        try
            event_type = "Sleep";
           
            fprintf("%s %d/%d\n", pat_id, pp, length(pat_ids));

            POD_str = per_patient_PODs(qq);
            POD = str2num(POD_str);

            % check to see if this POD has already been completed
            per_pod_csv = per_patient_csv(per_patient_csv{:, 3} == POD, :);
            if size(per_pod_csv, 1) > 0
                continue;
            end
    
            % get all the files for a night of sleep
            sleep_edf_files = get_files_in_dir(sleep_dir+pat_id+"\","*Night"+POD_str+"*_filt.edf",0);
    
            if length(sleep_edf_files) == 0
                sleep_text_files = get_files_in_dir(sleep_dir+pat_id+"\","*.txt",0);
                if length(sleep_text_files) == 0
                    fprintf("No text files or edf files \n\n");
                else
                    fprintf("No edf files\n\n");
                end
    
            end
        
            % TODO -> for the entire analysis, I'll need to add a second loop for
            % POD
        
            % get which POD the file is
            
            
            % conver the file names into a cell array b/c that is what SleepSEEG
            % requires
            sleep_seeg_fpath_cell = {};
            for qq = 1:length(sleep_edf_files)
                sleep_seeg_fpath_cell{qq} = char(sleep_edf_files(qq));
            end
            
        
            % run the sleep staging
            [Summary,SleepStage]=SleepSEEG(sleep_seeg_fpath_cell, pat_id);
            
        
            % Save the raw data
            write_raw_sleep_stage_csv(Summary, all_sleep_csv_fpath, POD, pat_id, event_type);
            
            % compute the confident data
            sleep_stage_labels = ["R", "W", "N1", "N2", "N3"];
            
            sleep_stages = SleepStage;
            sleep_stages(sleep_stages(:, 4)<0.5, 3) = 0;
            
            datetimes = datetime(sleep_stages(:, 2), 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss');
            
            stages = sleep_stages(:, 3)';
            d = [true, diff(stages) ~= 0, true];
            n = diff(find(d));               % Number of repetitions
            
            % remove the last index of d as it runs over
            d = d(1:end-1);
            
            % what index does each run of sleep phase start
            idx_start_sleep_phases = find(d);
            idx_stop_sleep_phases = [idx_start_sleep_phases(2:end), length(d)];
            stop_start_idxs = [idx_start_sleep_phases', idx_stop_sleep_phases'];
            
            % start/stop datetimes
            start_datetimes = datetimes(idx_start_sleep_phases);
            stop_datetimes = datetimes(idx_stop_sleep_phases);
            start_stop_datetimes = [start_datetimes, stop_datetimes];
            
            % what is the sleep stage for each of the phases
            sleep_stage_per_phase = stages(d);
            
            % remove all that are less than 5 minutes long
            stop_start_idxs(n<10, :) = [];
            start_stop_datetimes(n<10, :) = [];
            sleep_stage_per_phase(n<10) = [];
            n(n<10) = [];
            
            % remove all that are not confident
            stop_start_idxs(sleep_stage_per_phase==0, :) = [];
            start_stop_datetimes(sleep_stage_per_phase==0, :) = [];
            n(sleep_stage_per_phase==0) = [];
            sleep_stage_per_phase(sleep_stage_per_phase==0) = [];
            avg_certainty = [];
            for ii = 1:size(stop_start_idxs, 1)
                avg_certainty(ii) = mean(sleep_stages(stop_start_idxs(ii, 1):stop_start_idxs(ii, 2), 4)); 
            end
        
            % there are no stages that fit the requested criteria
            if isempty(stop_start_idxs)
                fid = fopen(log_fpath, "a");
                fprintf(fid, "No cleaned sleep events for pat:%s POD:%d\n", pat_id, POD);
                fclose(fid);
            else
                sleep_table = table();
                sleep_table{:, "Type"} = repmat(["Sleep"], size(start_stop_datetimes, 1), 1);
                sleep_table{:, "PatID"} = repmat(pat_id, size(start_stop_datetimes, 1), 1);
                sleep_table{:, "POD"} = repmat(POD, size(start_stop_datetimes, 1), 1);
                sleep_table{:, "SleepCat"} = sleep_stage_labels(sleep_stage_per_phase)';
                sleep_table{:, "onset_datetime"} = start_stop_datetimes(:, 1);
                sleep_table{:, "offset_datetime"} = start_stop_datetimes(:, 2);
                sleep_table{:, "duration"} = n' * 30;
                sleep_table{:, "AvgCertainty"} = avg_certainty';
                
                % write the csv
                cleaned_sleep_csv_fpath = "Z:\000_Data\SEEG\SEEG_Sleep_Staging\notes\cleaned_sleep_times_06122023.csv";
                write_cleaned_csv(sleep_table, cleaned_sleep_csv_fpath)
            end
        catch
            fid = fopen(log_fpath, "a");
            fprintf(fid, "ERROR in pat:%s\n", pat_ids(pp));
            fclose(fid);
        end
    end
end

% send_gmail_derek("Finished Sleep Staging", "That's it");

%% functions

function [] = write_raw_sleep_stage_csv(Summary, save_fpath, POD, pat_id, event_type)
    for ii = 2:size(Summary, 1)
        if Summary{ii, 2} == ' '
            Summary{ii, 2} = Summary{ii-1, 2};
        end
    end
    
    %
    fid = fopen(save_fpath, "a");
    % fprintf(fid, "Type\tPatID\tPOD\tSleepCat\tOnsetDatetime\tOffsetDatetime\tDuration\n");
    

    
    for ii = 2:(size(Summary, 1) - 1)
        
        onset_datetime_string = string(datetime(strcat(Summary{ii, 2}, '--', Summary{ii, 3}), 'InputFormat', 'dd-MMM-yyyy--HH:mm:ss', 'Format', 'yyyy-MM-dd HH:mm:ss'));
        offset_datetime_string = string(datetime(strcat(Summary{ii+1, 2}, '--', Summary{ii+1, 3}), 'InputFormat', 'dd-MMM-yyyy--HH:mm:ss', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    
        sleep_category = string(Summary{ii, 4});
    
        sleep_duration = Summary{ii, 5} * 30;
    
        fprintf(fid, "%s\t%s\t%d\t%s\t%s\t%s\t%d\n", event_type, pat_id, POD, sleep_category, onset_datetime_string, offset_datetime_string, sleep_duration);
    end
    
    % last one
    sleep_duration = Summary{end, 5} * 30;
    
    onset_datetime = datetime(strcat(Summary{end, 2}, '--', Summary{end, 3}), 'InputFormat', 'dd-MMM-yyyy--HH:mm:ss', 'Format', 'yyyy-MM-dd HH:mm:ss');
    offset_datetime = onset_datetime + seconds(sleep_duration);
    
    onset_datetime_string = string(onset_datetime);
    offset_datetime_string = string(offset_datetime);
    
    sleep_category = string(Summary{end, 4});
    
    fprintf(fid, "%s\t%s\t%d\t%s\t%s\t%s\t%d\n", event_type, pat_id, POD, sleep_category, onset_datetime_string, offset_datetime_string, sleep_duration);
    
    fclose(fid);
end


function [] = write_cleaned_csv(sleep_table, save_fpath)
    fid = fopen(save_fpath, "a");

    for ii = 1:size(sleep_table, 1)
        event_type = sleep_table{ii, "Type"};
        pat_id = sleep_table{ii, "PatID"};
        POD = sleep_table{ii, "POD"};
        sleep_category = sleep_table{ii, "SleepCat"};
        onset_datetime_string = string(sleep_table{ii, "onset_datetime"});
        offset_datetime_string = string(sleep_table{ii, "offset_datetime"});
        sleep_duration = sleep_table{ii, "duration"};
        certainty = sleep_table{ii, "AvgCertainty"};
    
        fprintf(fid, "%s\t%s\t%d\t%s\t%s\t%s\t%d\t%.3f\n", event_type, pat_id, POD, sleep_category, onset_datetime_string, offset_datetime_string, sleep_duration, certainty);
    end

    fclose(fid);
end