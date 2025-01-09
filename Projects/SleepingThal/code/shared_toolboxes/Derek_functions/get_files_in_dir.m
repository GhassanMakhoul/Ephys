function [file_paths, file_names] = get_files_in_dir(search_dir, search_term, varargin)
%GET_FILES_IN_DIR Function to get the path to files in a directory as a
%string array that is filtered
% Inputs:
%   search_dir: string path to the directory
%   search_term: string to search for in the directory
% Outputs:
%   file_paths: string array (nx1) paths to files
%   file_names: string array (nx1) files names without path

    if ispc
        slh = '\';
    elseif isunix
        slh = '/';
    else
        error("Unable to detect OS type");
    end

    if ~isfolder(search_dir)
        error("Search directory isn't a folder")
    end

    % 0 = not dir
    % 1 = is dir
    is_dir_flag = 0;

    % directory or not flag
    if length(varargin)==1
        is_dir_flag = varargin{1};
    end

    tmp_char = char(search_dir);
    if tmp_char(end) ~= slh
        search_dir = strcat(search_dir, slh);
    end

    struct_files = dir(strcat(search_dir,search_term));
    struct_logical_idxs = [struct_files.isdir] == is_dir_flag;

    cell_files = struct2cell(struct_files(struct_logical_idxs));
    string_fnames = string([cell_files(1,:)])';
    file_paths = repmat(search_dir,length(string_fnames),1) + string_fnames;

    % skip over the . and .. 
    if is_dir_flag
        inds_to_keep = [];
        for ii = 1:length(file_paths)
            tmp_fname = strsplit(file_paths(ii), slh);
            tmp_fname = tmp_fname(end);
            if tmp_fname == "." 
                % do nothing
            elseif tmp_fname == ".."
                % do nothing
            else
                inds_to_keep(end+1) = ii;
            end
        end

        file_paths = file_paths(inds_to_keep);
    end

    file_names = [""];
    for ii = 1:length(file_paths)
        tmp_fname = strsplit(file_paths(ii), slh);
        file_names(ii,1) = tmp_fname(end);
    end

end

