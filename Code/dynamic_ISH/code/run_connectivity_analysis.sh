#!bin bash

INP_DIR='/mnt/ernie_main/Ghassan/ephys/data/connectivity/'
OUT_DIR="/mnt/ernie_main/Price/ephys/data/periconnectivity/"
LOG_DIR="~/Documents/Research/Ephys/Code/dynamic_ISH/logs/"
CONFIG_FILE="config_connectivity.yml"
> ../logs/peri_conn.logs
for d in $INP_DIR*/; do 
    echo Running connectivity on $d 
    python connectivity_dynamics.py -d $d -p $OUT_DIR -l $LOG_DIR -c $CONFIG_FILE
done