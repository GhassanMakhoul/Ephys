function []= bipole_single_pulse(pat_id, output_dir)

addpath("/mnt/ernie_main/shared_toolboxes/Derek_functions");
pulse_root_dir = "/mnt/ernie_main/000_Data/SPES/data/preprocessed";
%% first get files 
inpdir = strcat(pulse_root_dir, "/", pat_id,"/",pat_id,"_CCEP_single_pulses" );
pulse_files = get_files_in_dir(inpdir, '*.mat',0);
root_dir = get_files_in_dir("/mnt/ernie_main/000_Data/SEEG/SEEG_EyesClosed_RestingState/data/",pat_id+"*", 1);
%% then run gen_bipole_montage_unix.m
for ii=1:numel(pulse_files)
    pulse_file = pulse_files(ii);
    tmp = split(pulse_file, "/");
    p_f_name = tmp(end);
    out_path = strcat(output_dir, "/", p_f_name );
  
    pulse = load(pulse_file);
    fs = pulse.fs;
    labels = string(pulse.labels);
    [filt_data, bip_montage_label, bip_montage_region] = generate_bipole_montage_unix(pulse.pulse, labels, fs, root_dir);
    labels = {};
    regions = {};
    for jj=1:numel(bip_montage_label)
        labels{jj,1} = char(bip_montage_label(jj));
        regions{jj,1} = char(bip_montage_label(jj));

    end
    save(out_path, 'filt_data','labels','regions', 'fs')
%% then save to output directory
end
exit;
end