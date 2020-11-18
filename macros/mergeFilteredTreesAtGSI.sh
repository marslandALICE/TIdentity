#!/bin/bash

###############################################################################
############################## Modification region ############################
localRunning=0   # 1 for testing locally, 0 for running on batch farm
batchCommand="/usr/bin/sbatch --get-user-env --time=8:00:00 --mem-per-cpu=6000 -J MCcl_Process  -o %j.%a.out -e  %j.%a.err"
###############################################################################
scriptDir=/lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts
filtTreeDirGSI=/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN2/LHC15o_pass1_NoSelection_06082019/mergedRuns/
runList=$filtTreeDirGSI/runs-2015-LHC15o-pass1.list
runListRCT=$filtTreeDirGSI/runs_15o_pass1_RCT.list
cdir=/u/marsland/PHD/macros/marsland_EbyeRatios
#cdir=$RUN_ON_GRID_DIR/Ebye/code

#
# // Post process the counting of chunks
# sort -nk6 Monalisa_15o_pass1.list > Monalisa_15o_pass1_sorted.list
# sort -nk5 chunkStat_2015o_pass1.list > chunkStat_2015o_pass1_sorted.list
# cat chunkStat_2015o_pass1_sorted.list | awk '{sum+=$5 ; print $0} END{print "sum=",sum}'
#

# nohup bash -c "source /u/marsland/PHD/macros/marsland_EbyeRatios/copyFromGRID.sh; copyFiles 000195721 /alice/data/2013/LHC13d/000195721/cpass1/ AliESDs_Barrel.root" &

# meld /u/marsland/PHD/macros/marsland_EbyeRatios/mergeFilteredTreesAtGSI.sh $RUN_ON_GRID_DIR/Ebye/code/mergeFilteredTreesAtGSI.sh

# split -a 4 -d -l 3000 results.list   --additional-suffix=.list Sub_
# split -a 4 -d -l 3000 tmp_hist.list  --additional-suffix=.list Sub_tmp_hist_
# awk '{$0="alien://"$0}1' Sub_tmp_0000.list > results_0000.list;
# awk '{$0="alien://"$0}1' Sub_tmp_hist_0000.list > results_hist_0000.list;

# tar -czvf AccScan.tar.gz AccScan    --> compress directory
# tar -xzvf AccScan.tar.gz            --> uncompress directory
# find . -type f | wc -l              --> count files

###########################################################################################################

# below is generic
###########################################################################################################
mergeFilteredTreesAtGSI()
{

  ###########################################################################################################
  # how to run:
  #
  if [ 1 -eq 0 ]; then
    cd /lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/ATO-465/data/LHC16g1/LHC16g1_pass1_140619
    cpEbyeMacros
    source /u/marsland/PHD/macros/marsland_EbyeRatios/mergeFilteredTreesAtGSI.sh;
    mergeFilteredTreesAtGSI 0   # 0: filtering; 1: filtering+mycode
  fi
  ###########################################################################################################

  trackTreesONOFF=$1
  if [ $trackTreesONOFF == 0 ]; then
    mergeDir=$filtTreeDirGSI/mergedFiltered
  elif [ $trackTreesONOFF == 1 ]; then
    mergeDir=$filtTreeDirGSI/mergeddEdx
  fi
  mkdir -p $mergeDir;  cd $mergeDir
  find $filtTreeDirGSI -iname "AnalysisResults*.root" | grep root > results_local_tree.list
  find $filtTreeDirGSI -iname "PtResHistograms*.root" | grep root > results_local_hist.list

  for run in $(cat $runList); do

    runName=000$run
    fileNameToMergeTree=results_local_tree_$runName.list;
    fileNameToMergeHist=results_local_hist_$runName.list;

    less $mergeDir/results_local_tree.list | grep $run > $fileNameToMergeTree;
    less $mergeDir/results_local_hist.list | grep $run > $fileNameToMergeHist;

    runDir=$mergeDir/$runName; mkdir -p $runDir;
    mv $fileNameToMergeTree $runDir
    mv $fileNameToMergeHist $runDir
    cd $runDir

    cp $cdir/mergeFilteredTreesAtGSI.sh .; chmod +x ./mergeFilteredTreesAtGSI.sh

    if [ $localRunning == 1 ]; then
      ./mergeFilteredTreesAtGSI.sh processMerge $runDir/$fileNameToMergeTree $runDir/$fileNameToMergeHist $trackTreesONOFF
    elif [ $localRunning == 0 ]; then
      eval $batchCommand ./mergeFilteredTreesAtGSI.sh processMerge $runDir/$fileNameToMergeTree $runDir/$fileNameToMergeHist $trackTreesONOFF
    fi

    cd $mergeDir

  done


}

###########################################################################################################
mergeFilteredTreesAndTracksAtGSI()
{

  ###########################################################################################################
  # how to run:
  #
  if [ 1 -eq 0 ]; then
    cpEbyeMacros
    source /u/marsland/PHD/macros/marsland_EbyeRatios/mergeFilteredTreesAtGSI.sh;
    mergeFilteredTreesAndTracksAtGSI 20
  fi
  ###########################################################################################################
  nSubFiles=$1
  fileNameToMergeHist="dummy.txt"

  mergeDir=$filtTreeDirGSI/mergeddEdx
  mkdir -p $mergeDir;  cd $mergeDir
  find $filtTreeDirGSI -iname "AnalysisResults.root" | grep root > results_local_tree.list


  for run in $(cat $runList); do

    runName=000$run
    fileNameRunWise=results_local_tree_$runName.list;

    less $mergeDir/results_local_tree.list | grep $run > $fileNameRunWise;

    runDir=$mergeDir/$runName; mkdir -p $runDir;
    mv $fileNameRunWise $runDir
    cd $runDir

    split -a 4 -d -l $nSubFiles $fileNameRunWise   --additional-suffix=.list Sub_$runName\_

    counter=0
    for isubFile in $( ls -lS Sub* | awk {'print$9'} ); do

      subDir=$runDir/Sub_$counter
      mkdir -p $subDir; cd $subDir; mv ../$isubFile .

      counter=$(($counter+1))
      cp $cdir/mergeFilteredTreesAtGSI.sh .; chmod +x ./mergeFilteredTreesAtGSI.sh

      if [ $localRunning == 1 ]; then
        ./mergeFilteredTreesAtGSI.sh processMerge $subDir/$isubFile "dummy.txt" 1
      elif [ $localRunning == 0 ]; then
        eval $batchCommand ./mergeFilteredTreesAtGSI.sh processMerge $subDir/$isubFile "dummy.txt" 1
      fi

      cd $runDir

    done

    cd $mergeDir
  done


}

###########################################################################################################
processMerge()
{

  fileNameToMergeTree=$1
  fileNameToMergeHist=$2
  trackTreesONOFF=$3
  ##
  if [ $trackTreesONOFF == 0 ]; then
    mkdir filteredTrees; cd filteredTrees;
    alihadd -i "V0s"    -i "highPt"          -i "dEdx"        -i "Laser" -i "MCEffTree" -i "CosmicPairs" \
            -i "events" -i "eventInfoTracks" -i "eventInfoV0" -s 2000000000 AnalysisResults_filteredTree.root @$fileNameToMergeTree
    alihadd AnalysisResults_filteredHist.root @$fileNameToMergeHist
    cd ..
  elif [ $trackTreesONOFF == 1 ]; then
    mkdir filteredTreesAndTracks; cd filteredTreesAndTracks;
    alihadd -i "mcMoms"  -i "mcGenMoms" -i "fullacc" -i "tracks"  -i "fArmPodTree" -i "eventVars" -i "dscaled" -i "mcFull" \
            -i "fTreeMC" -i "mcGen"     -i "V0s"     -i "highPt"  -i "dEdx"        -i "MCEffTree" -i "events" \
            -i "eventInfoTracks" -i "eventInfoV0" -s 2000000000 AnalysisResults_tree.root @$fileNameToMergeTree
    # alihadd -i "cleanHists" AnalysisResults_hist.root @$fileNameToMergeTree
    cd ..
  fi

}

###########################################################################################################
MakeEventTreesFromFiltered()
{

  ###########################################################################################################
  # how to run:
  #
  if [ 1 -eq 0 ]; then
    cpEbyeMacros
    source /u/marsland/PHD/macros/marsland_EbyeRatios/mergeFilteredTreesAtGSI.sh;

    MakeEventTreesFromFiltered 0   -100000 100000     0     0       0
    MakeEventTreesFromFiltered 0   -100000 1000       0.88  0       0

    MakeEventTreesFromFiltered 0    2000   100000     0.88  1       0
    MakeEventTreesFromFiltered 0    2000   100000     0.88  -1      0
    MakeEventTreesFromFiltered 0    6000   100000     0.88  1       0
    MakeEventTreesFromFiltered 0    6000   100000     0.88  -1      0

    MakeEventTreesFromFiltered 0    2000   8000       0.88  0       0

    hadd Debug.root     000*/Debug.root
    hadd HistsCorr.root 000*/HistsCorr.root
    hadd Hists.root     000*/Hists.root

    #RealData_FilterTreesMakeHists("ptot","p",0,file,0  ,2000,50000,    0.88, 1)   // full correction and selection
    # MakeEventTreesFromFiltered  period,  kPileUpLow,  kPileUpHigh,  kTimeSeriesEff,  kVzPileup
  fi
  ###########################################################################################################

  period=$1
  kPileUpLow=$2
  kPileUpHigh=$3
  kTimeSeriesEff=$4
  kVzPileup=$5
  kSystSetting=$6

  outbase=${PWD}

  AnaDir=$outbase/Settings_pileUp_$kPileUpLow\_$kPileUpHigh\_timeSeries_$kTimeSeriesEff\_vZPileUp_$kVzPileup\_SystSetting_$kSystSetting
  mkdir -p $AnaDir;  cd $AnaDir
  find $filtTreeDirGSI/mergeddEdx -iname "Analysis*tree*.root" | grep root > results_merged.list
  cp $scriptDir/RealData_FilterTreesMakeHists.C .
  cp $scriptDir/mergeFilteredTreesAtGSI.sh .

  for run in $(cat $runListRCT); do

    runName=000$run
    fileNameRunWise=results_merged_$runName.list;

    less $AnaDir/results_merged.list | grep $run > $fileNameRunWise;

    runDir=$AnaDir/$runName; mkdir -p $runDir;
    mv $fileNameRunWise $runDir
    cd $runDir

    NFiles=$(less $runDir/$fileNameRunWise | wc -l)
    ZBNFiles=$(($NFiles - 1))       # indexing of SLURM_ARRAY_TASK_ID starts from 0 that is why
    # submit jobs
    chmod +x $scriptDir/mergeFilteredTreesAtGSI.sh
    if [ $localRunning == 1 ]; then
      $scriptDir/mergeFilteredTreesAtGSI.sh ProcessEventTreesFromFiltered $runDir/$fileNameRunWise $period  $kPileUpLow $kPileUpHigh  $kTimeSeriesEff  $kVzPileup $kSystSetting
    elif [ $localRunning == 0 ]; then
      eval $batchCommand --array=0-$ZBNFiles $scriptDir/mergeFilteredTreesAtGSI.sh ProcessEventTreesFromFiltered  $runDir/$fileNameRunWise $period  $kPileUpLow $kPileUpHigh  $kTimeSeriesEff  $kVzPileup $kSystSetting
    fi

    cd $AnaDir
  done

  cd $outbase


}

###########################################################################################################
ProcessEventTreesFromFiltered()
{

  # inputs
  fileList=$1
  period=$2
  kPileUpLow=$3
  kPileUpHigh=$4
  kTimeSeriesEff=$5
  kVzPileup=$6
  kSystSetting=$7
  #
  if [ $localRunning == 1 ]; then
    index=0
  elif [ $localRunning == 0 ]; then
    index=$SLURM_ARRAY_TASK_ID
  fi
  counter=$(($index+1))
  line=`sed "${counter}q;d" $fileList`
  #
  # waiting job for merging
  if [ $index == 0 -a $localRunning == 0  ]; then
    rootFileHists=Hists.root
    rootFileHistsCorr=HistsCorr.root
    rootFileDebug=Debug.root
    eval $batchCommand -d afterany:$SLURM_ARRAY_JOB_ID $scriptDir/mergeFilteredTreesAtGSI.sh mergeSingleRunHists $rootFileHists $rootFileHistsCorr $rootFileDebug
  fi
  #
  # runthe main macro to create histogram
  mkdir file_$counter; cd file_$counter
  cp $scriptDir/mergeFilteredTreesAtGSI.sh $scriptDir/RealData_FilterTreesMakeHists.C $scriptDir/rootlogon.C .; chmod +x ./mergeFilteredTreesAtGSI.sh
  aliroot -q -b rootlogon.C RealData_FilterTreesMakeHists.C+\(\"ptot\",\"p\",$kSystSetting,\"$line\",$period,$kPileUpLow,$kPileUpHigh,$kTimeSeriesEff,$kVzPileup\)
  cd ..

}

###############################################################################
mergeSingleRunHists()
{
  rootFileHists=$1
  rootFileHistsCorr=$2
  rootFileDebug=$3
  hadd $rootFileDebug     `pwd`/file*/Debug_*.root
  hadd $rootFileHists     `pwd`/file*/Hists_*.root
  hadd $rootFileHistsCorr `pwd`/file*/HistsCorr_*.root
  mkdir Logs Files
  mv *.err *.out Logs;
  mv file_* Files;
  tar -czvf Files.tar.gz Files;
  tar -czvf Logs.tar.gz Logs;
  rm -rf Logs;
  rm -rf Files;
}

###########################################################################################################
copyFilteredTrees()
{
  ###############################################################################
  #  to run  -->
  if [ 1 -eq 0 ]; then
    ## run it in a screen session
    cd /lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/ATO-465/data/LHC16g1/LHC16g1_pass1_140619
    source /u/marsland/PHD/macros/marsland_EbyeRatios/mergeFilteredTreesAtGSI.sh;
    copyFilteredTrees /alice/cern.ch/user/p/pwg_pp/Skimmed_ESDs_16g1_pass1_140619/
  fi
  ###############################################################################
  alienBasePath=$1

  alien_find $alienBasePath *.root | grep root | grep -v EventStat > tmp.list
  awk '{$0="alien://"$0}1' tmp.list > results.list;
  alisync @results.list -o `pwd` -j 4 -t 3000
}

###########################################################################################################
compareLists()
{
  ###############################################################################
  #  to run  -->
  if [ 1 -eq 0 ]; then
    source $RUN_ON_GRID_DIR/Ebye/code/mergeFilteredTreesAtGSI.sh;
    compareLists
  fi
  ###############################################################################
  largeList=$RUN_ON_GRID_DIR/Ebye/lists/tiden_Filter_15o_pass1/chunkStat_2015o_pass1_sorted.list
  monalisaList=$RUN_ON_GRID_DIR/Ebye/lists/tiden_Filter_15o_pass1/Monalisa_15o_pass1_sorted.list

  matchRun=0
  countFiles=0
  counter=0
  fileName=runsSublist$counter\-2015-LHC15o-pass1.list
  fileNameAll=runsMonalisa-2015-LHC15o-pass1.list
  for run0 in $(less $monalisaList | awk {'print$1'}); do
    #
    matchRun=$(cat $largeList   | grep $run0 | awk {'print$1'})
    nEvents=$(cat $monalisaList | grep $run0 | awk {'print$6'})
    nFiles=$(cat $largeList     | grep $run0 | awk {'print$5'})

    countFiles=$(($countFiles+$nFiles))
    echo $countFiles  $nFiles

    if [ "$countFiles" -gt 150000 ]; then
      counter=$(($counter+1))
      countFiles=0
      fileName=runsSublist$counter\-2015-LHC15o-pass1.list
    fi

    less $largeList | grep $matchRun >> $fileName
    less $largeList | grep $matchRun >> tmp.list

    echo "crossCheck $matchRun = $run0   --> $nEvents   ---  $nFiles"

  done

  sort -nk5 tmp.list > $fileNameAll
  cat $fileNameAll | awk '{sum+=$5 ; print $0} END{print "sum=",sum}'
  rm tmp.list
  wc -l runs*

}

###########################################################################################################
RemoveZipFilesIn()
{
  ###############################################################################
  #  to run  -->
  if [ 1 -eq 0 ]; then
    source $RUN_ON_GRID_DIR/Ebye/code/mergeFilteredTreesAtGSI.sh;
    source /u/marsland/PHD/macros/marsland_EbyeRatios/mergeFilteredTreesAtGSI.sh;
    RemoveZipFilesIn /alice/cern.ch/user/p/pwg_pp/Skimmed_ESDs_dEdxPerformance

    # to check the size before and after
    source $AliPhysics_SRC/PWGPP/scripts/aliendu.sh
    aliendu /alice/cern.ch/user/p/pwg_pp/Skimmed_ESDs_dEdxPerformance  10000000000  1   20
  fi
  ###############################################################################

  fileDir=$1
  for zipFile in $(alien_find $fileDir/ "*.zip" | grep zip); do
    echo  "remove -->   $zipFile"
    alien_rm $zipFile
  done

}

###########################################################################################################
ResubmitJobsCopiedFromMonalisa()
{
  ###############################################################################
  # how to run:
  #
  if [ 1 -eq 0 ]; then
    cd /home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/test/tiden_filter_subsample0
    source /u/marsland/PHD/macros/marsland_EbyeRatios/mergeFilteredTreesAtGSI.sh;
    source $RUN_ON_GRID_DIR/Ebye/code/mergeFilteredTreesAtGSI.sh;
    ResubmitJobsCopiedFromMonalisa "TaskEbyeIterPIDMC"
    ResubmitJobsCopiedFromMonalisa "analysis"
  fi
  ###############################################################################
  taskStr=$1

  for job in $(less jobs.list | grep $taskStr | awk '{print $1}'); do
    for subjob in $(gbbox top -user pwg_pp -split $job  -all_status | grep -E 'ERROR|EXPIRED|ZOMBIE' | awk '{print $1}'); do
      echo $subjob >> subjobs.list
      alien_resubmit $subjob
    done
  done

}

###########################################################################################################
KillJobsCopiedFromMonalisa()
{
  ###############################################################################
  # how to run:
  #
  if [ 1 -eq 0 ]; then
    cd /home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/test/tiden_filter_subsample0
    source $RUN_ON_GRID_DIR/Ebye/code/mergeFilteredTreesAtGSI.sh;
    KillJobsCopiedFromMonalisa "TaskEbyeIterPIDMC"
    KillJobsCopiedFromMonalisa "analysis"
  fi
  ###############################################################################

  taskStr=$1

  for job in $(less jobs.list | grep  $taskStr | awk '{print $1}'); do
    alien_kill $job
  done


}

###########################################################################################################
CountFilesMC()
{

  ###############################################################################
  ###############################################################################
  #  to run  -->
  if [ 1 -eq 0 ]; then
    source $RUN_ON_GRID_DIR/Ebye/code/mergeFilteredTreesAtGSI.sh;
    CountFilesMC "/alice/sim/2017/LHC17c5b" "AliESDs.root"
  fi
  ###############################################################################
  ###############################################################################

  alienpath=$1
  fileType=$2

  counter=0
  for rundir in $(alien_ls $alienpath); do

    nFiles=$(alien_find $alienpath/$rundir $fileType | wc -l)
    echo $nFiles "        $rundir        $alienpath/$rundir"     >> chunkStat.list
    counter=$(($counter+$nFiles))

  done

  echo " ---------------------------------------------"
  echo " total number of files $counter for $alienpath"
  echo " ---------------------------------------------"

}

###########################################################################################################
CountFilesRealData()
{

  ###############################################################################
  ###############################################################################
  #  to run  -->
  if [ 1 -eq 0 ]; then
    source $RUN_ON_GRID_DIR/Ebye/code/mergeFilteredTreesAtGSI.sh;
    CountFilesRealData "/alice/data/2015/LHC15o" "pass1" "AliESDs.root"

    alien_ls /alice/data/2015/LHC15o/000*/pass5_lowIR/*/AliESDs.root | wc -l
    alien_ls /alice/data/2015/LHC15o/000*/pass1/*/AliESDs.root | wc -l
  fi
  ###############################################################################
  ###############################################################################

  alienpath=$1
  pass=$2
  fileType=$3

  rm chunkStat.list
  touch chunkStat.list

  nRuns=$(alien_ls $alienpath | grep 000 | wc -l)
  echo "Total number of files  =  $nRuns"

  counter=0
  for rundir in $(alien_ls $alienpath); do

    nFiles=$(alien_find $alienpath/$rundir/$pass $fileType | grep root | wc -l)
    echo " $counter        $rundir        $alienpath/$rundir   nFiles =  $nFiles "
    if [ "$nFiles" -gt "1" ]; then
      echo "being written to file   $rundir        $alienpath/$rundir   nFiles =  $nFiles "
      echo "        $rundir        $alienpath/$rundir   nFiles =  $nFiles "     >> chunkStat.list
    fi
    counter=$(($counter+$nFiles))

  done

  echo " ---------------------------------------------"
  echo " total number of files $counter for $alienpath"
  echo " ---------------------------------------------"

}

###########################################################################################################
CopyAnaResultsFromAlienDir()
{

  ###############################################################################
  ###############################################################################
  #  to run  -->
  if [ 1 -eq 0 ]; then
    source $RUN_ON_GRID_DIR/Ebye/code/mergeFilteredTreesAtGSI.sh;
    CopyAnaResultsFromAlienDir "/alice/cern.ch/user/p/pwg_pp/Skimmed_ESDs_dEdxPerformance_test/LHC16g1_pass1_20190809_051/" "AnalysisResults.root" "g1_1"
  fi
  ###############################################################################
  ###############################################################################

  alienpath=$1
  fileName=$2
  suffix=$3

  alien_find $alienpath $fileName | grep root | awk '{$0="alien://"$0}1' > results_$suffix\.list
  alisync @results_$suffix\.list -o `pwd` -j 4 -t 3000
  find ${PWD}/ -iname "*.root" > local_$suffix\.list

}

###############################################################################
main()
{
  eval "$@"
}
main "$@"
