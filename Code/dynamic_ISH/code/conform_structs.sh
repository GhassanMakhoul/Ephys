#!bin/bash
> ../logs/conform.logs
INP_DIR=/mnt/ernie_main/000_Data/SEEG/SEEG_Periictal/data/Connectivity/seizure_structs_pre10min_ictal_post10min/5sDur_1sStride
OUT_DIR=/mnt/ernie_main/Ghassan/ephys/data/connectivity
SUB_LIST=../data/pat_list.txt
while read -r subj; do 
    echo $subj
    ##Test to make sure directory exists before running
    if test ! -d $OUT_DIR/$subj; then
        echo "DIR DNE: $OUT_DIR/$subj"
        mkdir $OUT_DIR/$subj
    fi
    #For each file (.mat) convert struct
    for f_path in $INP_DIR/$subj/*PDC.mat; do
        f_name="$(basename $f_path)"
        echo $f_name
       /usr/local/MATLAB/R2024a/bin/matlab -nodisplay -nosplash -nodesktop -r "structconvert('$f_path', '$OUT_DIR/$subj/$f_name')"
       echo Saved to: $OUT_DIR/$subj/$f_name
    done
done < "$SUB_LIST"