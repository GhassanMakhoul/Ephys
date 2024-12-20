function [filt_data,bip_montage_label, bip_montage_region] = generate_bipole_montage_unix(monopole_data,channel_names,sfreq,pat_root_folder)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    pat_name = split(pat_root_folder,"/");
    pat_name = split(pat_name(end-1),"_");
    pat_name = pat_name(1);

    % read in the atlas assignment excel file
    csv_path = get_files_in_dir(pat_root_folder+"/", "*Assignments+EZdesignations.xlsx");
    
    if pat_name == "Epat09"
        csv_path = strcat(pat_root_folder, "/Epat09_COLOR_CODED_AtlasAssignments+EZdesignations_periictal_cleaned.xlsx");
    end

    atlas_assignments_csv = readcell(csv_path(1));
    
    bip_montage = [];
    bip_montage_label = string();
    bip_montage_region = string();

    channel_names = strrep(strrep(strrep(upper(channel_names),"ONE","1"),"TWO","2"),"THR","3");
    
    counter = 1;
    
    error_flag = 0; % if there is an error with the channels, then skip this file
    
    % loop through each region, skipping the header
    for jj = 2:size(atlas_assignments_csv,1)
    
        % loop through each bipole pair
        for ii = 4:size(atlas_assignments_csv,2)
            % stop if we've reached the end of the bipole pairs in a
            % region
            if ismissing(atlas_assignments_csv{jj,ii})
                break
            end
    
            contact_pair_together = string(atlas_assignments_csv{jj,ii});
            contact_pair = split(string(atlas_assignments_csv{jj,ii})," - ");

            contact_pair = strrep(strrep(strrep(upper(contact_pair),"ONE","1"),"TWO","2"),"THR","3");
            contact_pair = erase(erase(contact_pair,"-")," ");
    
            if (contact_pair(1)=="RTP9" | contact_pair(2)=="RTP9") & pat_name=="Epat09"
                continue
            end

            if isempty(find(strcmp(channel_names,contact_pair(1))))
                error_flag = -1;
                break
            end
            
            try
                contact_pair_idx(1) = find(strcmp(channel_names,contact_pair(1)));
                contact_pair_idx(2) = find(strcmp(channel_names,contact_pair(2)));
            catch
                fprintf("err");
            end
    
            bip_montage(counter,:) = monopole_data(contact_pair_idx(1),:) - monopole_data(contact_pair_idx(2),:);
            bip_montage_label(counter) = contact_pair_together;
            bip_montage_region(counter) = string(atlas_assignments_csv{jj,2});
    
            counter = counter + 1;
        
        end % end loop through each bipole pair
    
    end % end loop through each region 
    
    
    filt_data = bip_montage;



end
