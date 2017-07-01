#!/bin/bash

#
#   Example usage: 
#          /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/killAlienJobs.sh 0 marsland "WAITING" 1 
# 

masterjob=$1    # master job id to be taken from monalisa
username=$2   
jobmode=$3      # either RUNNING, WAITING, DONE, ASSIGED, EXPIRED, "" --> to kill all jobs
printmode=$4    # 1 to check the job list to be killed


if [ $printmode == 1 ]; then

  if [ $masterjob != 0 ]; then
    gbbox top -user $username -split $masterjob -all_status | grep "$jobmode" 
  else 
    gbbox top -user $username -all_status | grep "$jobmode" 
  fi
  
elif [ $printmode == 0 ]; then

  if [ $masterjob != 0 ]; then
    for jobid in $(gbbox top -user $username -split $masterjob -all_status | grep "$jobmode" | awk '{print $1}'); do 
      alien_kill $jobid; 
    done 
  else 
    for jobid in $(gbbox top -user $username -all_status | grep "$jobmode" | awk '{print $1}'); do 
      alien_kill $jobid; 
    done 
  fi
  
fi
  