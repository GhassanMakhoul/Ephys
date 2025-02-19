#!bin/bash
> ../logs/pdc.logs
#INP_DIR=/mnt/ernie_main/000_Data/SEEG/SEEG_Periictal/data/Connectivity/seizure_structs_pre10min_ictal_post10min/5sDur_1sStride
INP_DIR=/mnt/ernie_main/000_Data/SEEG/SEEG_Sleep_Staging/data/extracted_sleep_epochs_v2/
OUT_DIR=/mnt/ernie_main/Ghassan/ephys/data/sleep
CONFIG=config_sleepy.yml
SUB_LIST="/home/ghassan/Documents/Research/Ephys/Projects/SleepingThal/data/pat_list_tst.txt"
SOZ_LABELS="/mnt/ernie_main/000_Data/SEEG/SEEG_EyesClosed_RestingState/labels/all_pats_bipole_soz_labels.csv"
NIH_SHEET='/mnt/ernie_main/Ghassan/NIH12_Clinical_Info_11_15.xlsx'

n_trials=$(cat $CONFIG | shyaml get-value connectivity.n_trials)
win_duration=$(cat $CONFIG | shyaml get-value connectivity.win_duration)
viz_path=$(cat $CONFIG | shyaml get-value connectivity.viz_path)
win_duration=$(cat $CONFIG | shyaml get-value connectivity.win_duration)
shuffle=$(cat $CONFIG | shyaml get-value connectivity.shuffle)
metric=$(cat $CONFIG | shyaml get-value connectivity.metric)
overlap=$(cat $CONFIG | shyaml get-value connectivity.overlap)
agg=$(cat $CONFIG | shyaml get-value connectivity.agg)
counter=0
while read -r subj; do 
    echo $subj
    ##Test to make sure output directory exists before running
    if test ! -d $OUT_DIR/$subj; then
        echo "DIR DNE: $OUT_DIR/$subj"
        mkdir $OUT_DIR/$subj
    fi

    if test ! -d $viz_path/$subj; then
        echo "DIR DNE: $viz_path/$subj"
        mkdir $viz_path/$subj
    fi
    #For each 
    for f_path in $INP_DIR/$subj/*[0-9].mat; do #ensures only using bipolar recordings
        f_name="$(basename $f_path)"
        echo $f_path
        echo $viz_path/$subj
        /usr/local/MATLAB/R2024a/bin/matlab -nodisplay -nosplash -nodesktop -r \
        "calc_connectivity('$subj', '$f_path','$OUT_DIR/$subj/','$metric','$shuffle', '$win_duration', '$overlap','$n_trials', '$agg', '$viz_path/$subj', '$SOZ_LABELS', '$NIH_SHEET')"
       echo OUT DIR: $OUT_DIR/$subj/$f_name
       ((counter++))
       echo $counter
       if [ $counter -eq 1 ]; then
            break
        fi
    done
done < "$SUB_LIST"