#!bin/bash

INP_DIR=/mnt/ernie_main/000_Data/SEEG/SEEG_Periictal/data/Connectivity/seizure_structs_pre10min_ictal_post10min/5sDur_1sStride
OUT_DIR=/mnt/ernie_main/Ghassan/ephys/data/connectivity
for d in $INP_DIR/*pat*/; do
    subj_dir="$(basename $d)"
    echo $subj_dir
    ##Test to make sure directory exists before running
    if test ! -d $OUT_DIR/$subj_dir; then
        echo "DIR DNE: $OUT_DIR/$subj_dir"
        mkdir $OUT_DIR/$subj_dir
    fi
    #For each file (.mat) convert struct
    for f_path in $d/*.mat; do
        f_name="$(basename $f_path)"
        echo $f_name
        /usr/local/MATLAB/R2024a/bin/matlab -nodisplay -nosplash -nodesktop -r "structconvert('$f_path', '$OUT_DIR/$subj_dir/$f_name')"
    done
done