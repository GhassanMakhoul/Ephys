#!/bin/bash

while read p; do
  sub_id=$p
  echo "Running on " $sub_id
  outdir="/mnt/ernie_main/000_Data/SPES/data/preprocessed/"$sub_id"/"$sub_id"_CCEP_single_pulses_bipole/"
  mkdir -p $outdir
  /usr/local/MATLAB/R2023b/bin/matlab -nodisplay -nosplash -nodesktop -r "bipole_single_pulse('$sub_id','$outdir');"
  echo "SUCCESS! for " $sub_id
done < subj.txt


# change tst to subj.txt