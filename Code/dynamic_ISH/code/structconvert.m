function [] = structconvert(INP_F,OUT_F)
    %INP_F = '/mnt/ernie_main/PCA_project/data/81_structs_with_dkregion_and_soz/pat_structs_with_dkregion_and_soz.mat';
    %OUT_F = "/mnt/ernie_main/Ghassan/ephys/data/og_ISH.mat";
    
    STR_FIELDS = ["patID", "eventID", "sz_type", "awareness_label","bip_labels_used", "region_name"];
    
    pdc = load(INP_F);
    disp(["LOADED", INP_F])
    
    
    %for OG ISH cohort
    %n_entries = max(size(seizure.pats))
    %for fname = STR_FIELDS
    %    for i=1:n_entries
    %        seizure.pats(i).(fname) = convertStringsToChars(seizure.pats(i).(fname));
    %    end
    %end
    
    for field = STR_FIELDS
        pdc.seizure.(field) = convertStringsToChars(pdc.seizure.(field));
    end
    disp(["Saving to ", OUT_F])
    save(OUT_F, "pdc", '-v7.3');
    % exit;

end