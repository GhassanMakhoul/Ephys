#!bin bash

INP_DIR='/mnt/ernie_main/Ghassan/ephys/data/connectivity/'
OUT_DIR="/mnt/ernie_main/Ghassan/ephys/data/periconnectivity/"
LOG_DIR="~/Documents/Research/Ephys/Code/dynamic_ISH/logs/"

for d in $INP_DIR*/; do 
 python connectivity_dynamics.py -d $d -p $OUT_DIR -l $LOG_DIR
done