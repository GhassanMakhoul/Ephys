#!bin/bash
> ../logs/pdc.logs
#INP_DIR=/mnt/ernie_main/000_Data/SEEG/SEEG_Periictal/data/Connectivity/seizure_structs_pre10min_ictal_post10min/5sDur_1sStride
INP_DIR=/mnt/ernie_main/000_Data/SEEG/SEEG_Sleep_Staging/data/extracted_sleep_epochs_v2/
OUT_DIR=/mnt/ernie_main/Ghassan/ephys/data/sleep
CONFIG=config_sleepy.yml
SUB_LIST="/home/ghassan/Documents/Research/Ephys/Projects/SleepingThal/data/pat_list_tst.txt"


n_trials=$(cat $CONFIG | shyaml get-value connectivity.n_trials)
win_duration=$(cat $CONFIG | shyaml get-value connectivity.win_duration)
viz_path=$(cat $CONFIG | shyaml get-value connectivity.viz_path)
win_duration=$(cat $CONFIG | shyaml get-value connectivity.win_duration)
shuffle=$(cat $CONFIG | shyaml get-value connectivity.shuffle)
metric=$(cat $CONFIG | shyaml get-value connectivity.metric)
agg=$(cat $CONFIG | shyaml get-value connectivity.agg)

while read -r subj; do 
    echo $subj
    ##Test to make sure directory exists before running
    if test ! -d $OUT_DIR/$subj; then
        echo "DIR DNE: $OUT_DIR/$subj"
        mkdir $OUT_DIR/$subj
    fi
    #For each 
    for f_path in $INP_DIR/$subj/*[0-9].mat; do #ensures only using bipolar recordings
        f_name="$(basename $f_path)"
        echo $f_path
        echo $viz_path
        /usr/local/MATLAB/R2024a/bin/matlab -nodisplay -nosplash -nodesktop -r \
        "calc_connectivity('$subj', '$f_path','$OUT_DIR/$subj/','$metric','$shuffle', '$win_duration', '$n_trials', '$agg', '$viz_path')"
       echo OUT DIR: $OUT_DIR/$subj/$f_name
       echo 
       break
    done
done < "$SUB_LIST"