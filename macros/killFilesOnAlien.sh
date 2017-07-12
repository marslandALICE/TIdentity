#!/bin/bash

#
#   Example usage: 
#  /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/killFilesOnAlien.sh /alice/cern.ch/user/m/marsland/EbyeIterPID_MC_FULL/LHC17c5b_pass5_lowIR/output/244917 EventStat_temp.root 1 
# 

#AnalysisResults.root
#EventStat_temp.root
#log_archive.zip
#root_archive.zip
#stderr
#stdout


alienpath=$1    # master job id to be taken from monalisa
filename=$2   
printmode=$3    # 1 to check the job list to be killed

if [ $printmode == 1 ]; then
    for a in $(alien_find $alienpath $filename | grep alice); do
      echo "alien_rm $a"
    done
elif [ $printmode == 0 ]; then
   for a in $(alien_find $alienpath $filename | grep alice); do
      alien_rm $a
   done
fi
  