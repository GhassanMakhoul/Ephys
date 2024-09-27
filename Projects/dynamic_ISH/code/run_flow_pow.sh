#!bin bash

POWDIR="/mnt/ernie_main/Ghassan/ephys/data/ei_bal/"
FLOWDIR="/mnt/ernie_main/Price/ephys/data/periconnectivity"
LOG_DIR="~/Documents/Research/Ephys/Code/dynamic_ISH/logs/"
CONFIG_FILE="config_connectivity.yml"
SUBJLIST=subjlist.csv
PLTNAME="default_plot.pdf"

Help()
{
    #Display Help
    echo "This script calls the eibal.py python file to plot level measures of peri-ictal flow vs power"
    echo "Syntax: scriptTemplate [-h|i|o|c|l|s]"
    echo "options:"
    echo "  -h     Print this Help."
    echo "  -fd    input directory for connectivity channel wise flow"
    echo "  -pd    input directory for channel wise power"
    echo "  -c     config file (should be yml)"
    echo "  -l     log directory"
    echo "  -s     subject(s) to analyze list" 
    echo "The following arguments will default to the following "
    echo "-pd : Input Directory of power structs:"
    echo $POWDIR
    echo "-fd : flow dir directory"
    echo $FLOWDIR
    echo "-c : The config file"
    echo $CONFIG_FILE
    echo "-s the location of the list of subjects to process"
    echo $SUBJLIST
    echo "-l :  the Log directory"
    echo $LOG_DIR
}
> ../logs/peri_flow_pow.logs

while  [ $# -gt 0 ] ; do  #getopts ":h:fd:pd:c:plt:-s:l:" option;
    case $1 in
        -h)
            Help
            exit;;
        -fd | --flowdir) FLOWDIR=$OPTARG;;
        -i) POWDIR="$2";;
        -c) CONFIG_FILE="$2";;
        -plt) PLTNAME="$2";;
        -s) SUBJLIST="$2";;
        -l) LOG_DIR="$2";;
        \?)
            echo "ERROR: Invalid option"
            exit;;
    esac
    shift
done

python plotpowerbal.py -s $SUBJLIST -f $FLOWDIR -p $POWDIR -t $PLTNAME -l $LOG_DIR -c $CONFIG_FILE