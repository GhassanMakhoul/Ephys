#!bin bash

INP_DIR='/mnt/ernie_main/Ghassan/ephys/data/connectivity/'
OUT_DIR="/mnt/ernie_main/Ghassan/ephys/data/ei_bal/"
LOG_DIR="~/Documents/Research/Ephys/Code/dynamic_ISH/logs/"
CONFIG_FILE="config_connectivity.yml"
SUBJ='all'

Help()
{
    #Display Help
    echo "This script calls the eibal.py python file to calculate channel level measures of peri-ictal power"
    echo "Syntax: scriptTemplate [-h|i|o|c|l|s]"
    echo "options:"
    echo "  -h     Print this Help."
    echo "  -i     input directory"
    echo "  -o     output directory"
    echo "  -c     config file (should be yml)"
    echo "  -l     log directory"
    echo "  -s     subject(s) to analyze" 
    echo "         NOTE: if passing more than one subject always include quotations around list of subjects"
    echo '         Example -s "Epat20 Spat5'
    echo "The following arguments will default to the following "
    echo "-i : Input Directory of connaectivity structs:"
    echo $INP_DIR
    echo "-o : Output directory"
    echo $OUT_DIR
    echo "-c : The config file"
    echo $CONFIG_FILE
    echo "-l :  the Log directory"
    echo $LOG_DIR
}
> ../logs/peri_eibal.logs

while getopts ":h:o:i:c:s:l:" option; do
    case $option in
        h)
            Help
            exit;;
        o)
            OUT_DIR=$OPTARG;;
        i)
            INP_DIR=$OPTARG;;
        c)
            CONFIG_FILE=$OPTARG;;
        s)
            SUBJ=$OPTARG;;
        l)
            LOG_DIR=$OPTARG;;
        \?)
            echo "ERROR: Invalid option"
            exit;;
    esac
done

if [ "$SUBJ" = "all" ]; then
    for d in $INP_DIR*/; do 
        echo Running Spectral Analysis on $d 
        python eibal.py -d $d -p $OUT_DIR -c $CONFIG_FILE
    done
    exit;
fi

read -r -a splitArray <<<"$SUBJ"

for a in "${splitArray[@]}"; do
    echo Running for $a
    python eibal.py -d $INP_DIR"$a"/ -p $OUT_DIR -c $CONFIG_FILE
done