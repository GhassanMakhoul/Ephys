function [pat_metadata] = get_pat_info(subj, demographics_xls_file,columns)
    %%% Reads in the NIH clinical_info sheet and generates a set of rows
    %%% to add as meta data for each SUBJ passed into the method
    %%% Data retrieved: 
    %               - unilateral_mTLE
    %               - Epilepsy Type
    %               - Gender, 
    %               - age at scan,
    %               - FAS_Freq, FIAS_Freq, FBTC_Freq,
    %               - Pathology, 
    %               - Engel 6 month Num, 1yr, 2yr, 3yr
    %%% 
    %%%
    %%%
    
    sheet = 'Patients'; %hardcoding based on NIH12_clinical_info
    % in the case where columns are not passed hardcoding the columns headers from the NIH_clinical info sheet
    % note that if you are running this with a different info sheet, these headers may not 
    % agreee with your
    if numel(columns) == 0
        columns = {'Unilateral mTLE? (1:Yes, 2:No)',...
         'Presumed Epilepsy Type (0=mTLE, 1=TLE, 2=PLE, 3=FLE, 4=Multifocal, 5=PVNH/Cav Malformation, 6=IGE, 7=PNES, 8=unknown, 9=insular, 10=other)',...
         'Gender  (1:female, 2:male)','Age at scan (or SEEG if Spat)','Age Onset', 'Duration, yrs', 'FAS_Freq_monthly',...
          'FIAS_Freq_monthly', 'FBTC_Freq_monthly', 'TotalSz_Freq_Monthly', ...
          'EngelOutcome6mo Num (0: Engel 1A, 1: Engel 1; 2: Engel 2; 3: Engel 3; 4 Engel 4)',...
          'EngelOutcome1yr Num (0: Engel 1A, 1: Engel 1; 2: Engel 2; 3: Engel 3; 4 Engel 4)',...
          'EngelOutcome2yr Num (0: Engel 1A, 1: Engel 1; 2: Engel 2; 3: Engel 3; 4 Engel 4)',...
          'EngelOutcome3yr Num (0: Engel 1A, 1: Engel 1; 2: Engel 2; 3: Engel 3; 4 Engel 4)'};
    end
    nih_table = readtable(demographics_xls_file, 'Sheet', sheet, 'TextType','string');

    og_table_headers = nih_table.Properties.VariableDescriptions;
    table_headers = nih_table.Properties.VariableNames;

    
    disp(table_headers{17})
    pat_ids = nih_table.Patient;
    matches = contains(pat_ids, subj, "IgnoreCase",true);
    pat_row = nih_table(matches, :);
    pat_metadata = struct;
    if height(pat_row) >1
        warning("Multiple matches found for patID: %s", subj)
    elseif height(pat_row) ==0
        alt_ids = nih_table.Patient_secondary_name;
        matches = contains(alt_ids, subj, "IgnoreCase",true);
        if height(matches) == 0
            warning("NO Match found for patient %s", subj)
        end
    end
    for i=1:length(columns)
        column_name = columns{i};
        bool_ix = contains(og_table_headers, column_name);
        if any(bool_ix)
            %going to pick out the nice header to use
            hdr = table_headers{bool_ix};
            pat_metadata.(hdr) = pat_row.(hdr);
        else
            warning('Column "%s"  does not exist in this verion of the NIH12_clinical info sheet', column_name)
        end
    end
end