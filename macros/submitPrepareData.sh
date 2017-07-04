#!/bin/bash

#
#  Hists     -> job submit to produce histograms of different eta and cent 
#  Cleans    -> job submit to produce clean Samples of different eta and cent 
#  DirCheck  -> to check if directory exist
#  FileCheck -> to check if file exist
#  processCleans -> run aliroot for a single eta-cent conf for clean samples
#  processHists  -> run aliroot for a single eta-cent conf for clean samples
#  SentAllJos    -> sent all jobs to to farm (~10000)
#  mergeData     -> merge All data fir different eta-cent conf
#

# some global variables
DATE=`date +%Y%m%d_%H_%M`
DATE=`date +%m%d_%H_%M`
# qsubCommand="qsub -cwd -V -l h_rt=24:0:0,h_rss=4G -P alice -b y -r y -o out.log -e err.log"
# qsubCommand="qsub -cwd -V -l h_rt=24:0:0,h_rss=6G -b y -r y -o out.log -e err.log"
qsubCommand="/usr/bin/sbatch --get-user-env --time=8:00:00 --mem-per-cpu=6000  -o out.log -e err.log"
histName="cleanHists"
shDir=/lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts
allDataDir=/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data

############################## Modification region ############################
###############################################################################

localRunning=0   # 1 for testing locally, 0 for running on batch farm
testScripts=0    # 1 for just to see if macro works do not process any macro. 0 to run the macros either local on grid

################### real data setting ####################
mainDir=$allDataDir/PbPb/Real/RUN1/10h_pass2_cRows_80_10EtaBin_mombin20MeV
mcFile=$allDataDir/PbPb/MC/RUN1/11a10abis_MCclosure_cRows_80_10EtaBin_mombin20MeV/mergedPeriods/AnalysisResults.root
hisDir=$mainDir/Hists
dataDir=$mainDir/mergedPeriods               
histFile=$dataDir/AnalysisResults_Hists.root
idenFile=$dataDir/AnalysisResults.list
cleanFile=$dataDir/CleanSampleTree.root

SourceOfHists=2    # 1 for datatree, 0 for identree, 2 for the ThnSparses are used for the hist creation
nSubSample=21
treeFileName="fIdenTree"
# etaDownArray=(-0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7)
# etaUpArray=(  -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1  0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
# etaDownArray=(-0.8 -0.75 -0.7 -0.65 -0.6 -0.55 -0.5 -0.45 -0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75)
# etaUpArray=(-0.75 -0.7 -0.65 -0.6 -0.55 -0.5 -0.45 -0.4 -0.35 -0.3 -0.25 -0.2 -0.15 -0.1 -0.05 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8)
etaDownArray=(-1.0 -0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8)
etaUpArray=(  -0.8 -0.6 -0.4 -0.2  0.0 0.2 0.4 0.6 0.8 1.0)
centDownArray=(0.0 5.0  10.0 20.0 30.0 40.0 50.0 60.0 70.0)
centUpArray=(  5.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0)
# centDownArray=(0.0 2.5 5.0 7.5   10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0 75.0)
# centUpArray=(  2.5 5.0 7.5 10.0  15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0 75.0 80.0)
fitSkewnessArr=(0.7 0.5 0.5)
fitKurtosisArrPi=(1.)
fitKurtosisArrKaPr=(1.04 1. 0.96)
signArr=(-1 0 1)
##########################################################

################### MC data setting ####################
# SourceOfHists=2    # 1 for datatree, 0 for identree, 2 for the ThnSparses are used for the hist creation
# fastGen=0          # 1 for only mcGen, 0 for both MCgen and MCrec
# mainDir=$allDataDir/PbPb/MC/RUN1/11a10abis_MCclosure_cRows_80_10EtaBin_mombin20MeV
# EffCheckDir=$mainDir/EfficiencyCheck  
# hisDir=$mainDir/Hists
# dataDir=$mainDir/mergedPeriods               
# histFile=$dataDir/AnalysisResults_Hists.root
# mcFile=$dataDir/marsland_MCmoments.root
# idenFile=$dataDir/AnalysisResults.list
# cleanFile=$dataDir/CleanSampleTree.root

# nSubSample=16
# treeFileName="fIdenTreeMC"
# signArr=(-1 0 1)
# etaDownArray=(-1.0 -0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8)
# etaUpArray=(  -0.8 -0.6 -0.4 -0.2  0.0 0.2 0.4 0.6 0.8 1.0)
# centDownArray=(0.0 5.0  10.0 20.0 30.0 40.0 50.0 60.0 70.0)
# centUpArray=(  5.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0)
# 
# # etaDownArrayEffScan=(-0.8)
# # etaUpArrayEffScan=(   0.8)
# # pDownArrayEffScan=(0.2)
# # pUpArrayEffScan=(  1.5)
# # etaDownArrayEffScan=(-0.5 -0.8)
# # etaUpArrayEffScan=(   0.5  0.8)
# # pDownArrayEffScan=(0.2 0.2 0.3 0.3)
# # pUpArrayEffScan=(  1.5 1.8 1.5 1.8)
# etaDownArrayEffScan=(-0.8)
# etaUpArrayEffScan=(   0.8)
# pDownArrayEffScan=(0.2 0.2 0.6 0.6)
# pUpArrayEffScan=(  1.5 1.8 1.5 1.8)
#########################################################

########################## end of Modification region #########################
###############################################################################

echo " ============= All functions are sourced, have FUN !!!! ============= "
(source $AliPhysics_SRC/PWGPP/scripts/utilities.sh; hostInfo) 1>hostInfo.log

SentHists()
{
###############################################################################
###############################################################################
#  to run  -->   
if [ 1 -eq 0 ]; then
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_Reference/Hists
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_Reference/Hists
  cpEbyeMacros
  source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitPrepareData.sh 
  SentHists > out.log 2>&1 & 
fi
###############################################################################
###############################################################################
outbase=`pwd`

# define eta values
counter=0
for (( icent = 0; icent < ${#centDownArray[@]}; icent++ )); do
  for (( ieta = 0; ieta < ${#etaDownArray[@]}; ieta++ )); do
  
   centLow=${centDownArray[$icent]}
   centUp=${centUpArray[$icent]}
   etaLow=${etaDownArray[$ieta]}
   etaUp=${etaUpArray[$ieta]}
   counter=$(($counter+1))
   
   echo $counter " centrality interval =   " $centLow - $centUp "       eta Window =   " $etaLow - $etaUp
   cp $shDir/rootlogon.C $shDir/CreateAllTIDENHists.C $shDir/submitPrepareData.sh  .
   chmod +x submitPrepareData.sh
   [[ $testScripts == 1 ]] && continue
   
   sleep 0.5
   # submit jobs
   if [ $localRunning == 1 ]; then
      ./submitPrepareData.sh processHists $etaLow $etaUp $centLow $centUp
   elif [ $localRunning == 0 ]; then
      eval $qsubCommand ./submitPrepareData.sh processHists $etaLow $etaUp $centLow $centUp
   fi
   
  done 
done 

cd $outbase

}
###############################################################################
processHists()
{

# inputs
etaLow=$1
etaUp=$2
centLow=$3
centUp=$4

# Aliroot
# source  /lustre/nyx/alice/users/marsland/bin/.tpcdevEnv

# run aliroot for each file
if [ $SourceOfHists == 1 ]; then
   echo  "make hists using dataTree  -->  $etaLow $etaUp $centLow $centUp "
   aliroot -q -b $shDir/rootlogon.C $shDir/CreateAllTIDENHists.C\(0\,0\,0\,\"\",\"\",\"$dataFile\",\"$mcFile\",\"$histFile\",$etaLow,$etaUp,$centLow,$centUp\) 
elif [ $SourceOfHists == 0 ]; then
   echo  "make hists using Identree  -->  $etaLow $etaUp $centLow $centUp "
   aliroot -q -b $shDir/rootlogon.C $shDir/CreateAllTIDENHists.C\(3\,0\,0\,\"\",\"$idenFile\",\"$dataFile\",\"$mcFile\",\"$histFile\",$etaLow,$etaUp,$centLow,$centUp\) 
elif [ $SourceOfHists == 2 ]; then
   echo  "make hists using Identree  -->  $etaLow $etaUp $centLow $centUp "
   aliroot -q -b $shDir/rootlogon.C $shDir/CreateAllTIDENHists.C\(4\,0\,0\,\"\",\"\",\"$cleanFile\",\"$mcFile\",\"$histFile\",$etaLow,$etaUp,$centLow,$centUp\) 
fi

}
###############################################################################
SentMainSubSample()
{
###############################################################################
###############################################################################
#  to run  -->   
if [ 1 -eq 0 ]; then
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_Reference/SubSamples
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_Reference/SubSamples
  cpEbyeMacros
  source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitPrepareData.sh 
  SentMainSubSample > out.log 2>&1 & 
fi
###############################################################################
###############################################################################
outbase=`pwd`

# loop over centralities
for (( icent = 0; icent < ${#centDownArray[@]}; icent++ )); do

   centDir=$outbase/cent_${centDownArray[$icent]}
   DirCheckDelete $centDir; mkdir $centDir; cd $centDir
   echo "centrality bin =   " $icent 
   cp $shDir/rootlogon.C $shDir/CreateAllTIDENHists.C $shDir/submitPrepareData.sh .
   chmod +x submitPrepareData.sh
   [[ $testScripts == 1 ]] && continue
   
   if [ $localRunning == 1 ]; then
       ./submitPrepareData.sh processSubSamples 1 $icent 0 $idenFile  
   elif [ $localRunning == 0 ]; then
       eval $qsubCommand ./submitPrepareData.sh processSubSamples 1 $icent 0 $idenFile
   fi
   
done 
cd $outbase

}
###############################################################################
###############################################################################
SentSubSamples()
{
###############################################################################
###############################################################################
#  to run  -->   
if [ 1 -eq 0 ]; then
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_Reference/SubSamples
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_Reference/SubSamples
  cpEbyeMacros
  source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitPrepareData.sh 
  SentSubSamples > out.log 2>&1 & 
fi
###############################################################################
###############################################################################
outbase=`pwd`

# loop over centralities
for (( icent = 0; icent < ${#centDownArray[@]}; icent++ )); do

   centDir=$outbase/cent_${centDownArray[$icent]}
   cd $centDir
   rootFilePath=$(ls $(pwd)/*ss0.root)
   echo "root file to be processed -->  "$rootFilePath
   for (( isamp = 1; isamp < $nSubSample; isamp++ )); do
   
     echo "centrality bin =   " $icent "   sample no =  " $isamp
     cp $shDir/rootlogon.C $shDir/CreateAllTIDENHists.C $shDir/submitPrepareData.sh .
     chmod +x submitPrepareData.sh
     [[ $testScripts == 1 ]] && continue

     if [ $localRunning == 1 ]; then
        ./submitPrepareData.sh processSubSamples 2 $icent $isamp $rootFilePath
     elif [ $localRunning == 0 ]; then
        eval $qsubCommand ./submitPrepareData.sh processSubSamples 2 $icent $isamp $rootFilePath
     fi
     
   done
done 
cd $outbase

}
###############################################################################
SentMCEffCheck()
{
###############################################################################
###############################################################################
#  to run  -->   
if [ 1 -eq 0 ]; then
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_HIJING_TightCuts/EfficiencyCheck
  cpEbyeMacros
  source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitPrepareData.sh 
  SentMCEffCheck > out.log 2>&1 & 
fi
###############################################################################
###############################################################################
outbase=`pwd`
mcEffCheckDir=$outbase/MC_EffCheck_$DATE
DirCheckDelete $mcEffCheckDir; mkdir $mcEffCheckDir; cd $mcEffCheckDir
# loop over centralities
for (( idatatype = 0; idatatype < 2; idatatype++ )); do
   for (( ieta = 0; ieta < ${#etaDownArrayEffScan[@]}; ieta++ )); do
      for (( imom = 0; imom < ${#pDownArrayEffScan[@]}; imom++ )); do
        for (( isamp = 0; isamp < $nSubSample; isamp++ )); do
          [[ $fastGen == 1 && $idatatype == 0 ]] && continue
          momDown=${pDownArrayEffScan[$imom]}
          momUp=${pUpArrayEffScan[$imom]}
          etaDown=${etaDownArrayEffScan[$ieta]}
          etaUp=${etaUpArrayEffScan[$ieta]}
          echo "Momentum Range =   " $momDown - $momUp "   eta Range =  " $etaDown - $etaUp "   data type =  " $idatatype

          # create a dir and make production inside
          mometaDir=$mcEffCheckDir/Mom_$isamp\_$idatatype\_$momDown\_$momUp\_Eta_$etaDown\_$etaUp
          mkdir $mometaDir; cd $mometaDir

          cp $shDir/rootlogon.C $shDir/TIdentitySubSampleMC.C $shDir/submitPrepareData.sh .
          chmod +x submitPrepareData.sh
          [[ $testScripts == 1 ]] && continue
     
          # call the submitter
          if [ $localRunning == 1 ]; then
             ./submitPrepareData.sh processEffCheck $etaDown $etaUp $momDown $momUp $idatatype $isamp
          elif [ $localRunning == 0 ]; then
             eval $qsubCommand ./submitPrepareData.sh processEffCheck $etaDown $etaUp $momDown $momUp $idatatype $isamp 
          fi
        done
      done 
   done
done 
cd $outbase

}
###############################################################################
processEffCheck()
{

# inputs
etaDown=$1
etaUp=$2
momDown=$3
momUp=$4
datatype=$5
sampNo=$6

# run aliroot for each file
aliroot -q -l -b rootlogon.C TIdentitySubSampleMC.C\(\"$mcFile\"\,$sampNo\,$datatype\,$momDown\,$momUp\,$etaDown\,$etaUp\) 

}
###############################################################################
processSubSamples()
{

# inputs
idenOpt=$1
centBin=$2
subSample=$3
rootFile=$4

# Aliroot
# source  /lustre/nyx/alice/users/marsland/bin/.tpcdevEnv

# run aliroot for each file
aliroot -q -l -b $shDir/rootlogon.C $shDir/CreateAllTIDENHists.C\($idenOpt,$centBin\,$subSample,\"$treeFileName\",\"$rootFile\",\"\",\"\",\"\",0\,0\,0\,0\) 

}
###############################################################################
FitScan()
{
###############################################################################
###############################################################################
if [ 1 -eq 0 ]; then
   cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_mombin20MeV/Fits
   cpEbyeMacros
   source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitPrepareData.sh 
   FitScan 150 0.2 3.2 6 > outslice.log 2>&1 &
fi
###############################################################################
###############################################################################

# inputs
nSlice=$1
ptMin=$2
ptMax=$3
maxIter=$4

# base directory
outbase=`pwd`

# create the fit results directory
fitScanDir=$outbase/FitScan_$DATE
DirCheckDelete $fitScanDir; mkdir $fitScanDir; cd $fitScanDir

# define eta values only for skewness scan
# for (( i = 0; i < ${#fitSkewnessArr[@]}; i++ )); do
#   for (( j = 0; j < ${#fitSkewnessArr[@]}; j++ )); do
#     for (( k = 0; k < ${#fitSkewnessArr[@]}; k++ )); do
#        piSkew=${fitSkewnessArr[$i]}
#        kaSkew=${fitSkewnessArr[$j]}
#        prSkew=${fitSkewnessArr[$k]}
#        $shDir/submitPrepareData.sh Slices $nSlice $ptMin $ptMax $maxIter $piSkew $kaSkew $prSkew
#      done 
#   done 
# done 
    
for (( i = 0; i < ${#fitKurtosisArrPi[@]}; i++ )); do
  for (( j = 0; j < ${#fitKurtosisArrKaPr[@]}; j++ )); do
    for (( k = 0; k < ${#fitKurtosisArrKaPr[@]}; k++ )); do

      piKurt=${fitKurtosisArrPi[$i]}
      kaKurt=${fitKurtosisArrKaPr[$j]}
      prKurt=${fitKurtosisArrKaPr[$k]}

      piSkew=${fitSkewnessArr[0]}
      kaSkew=${fitSkewnessArr[1]}
      prSkew=${fitSkewnessArr[2]}

      # submit jobs
      $shDir/submitPrepareData.sh Slices $nSlice $ptMin $ptMax $maxIter $piSkew $kaSkew $prSkew $piKurt $kaKurt $prKurt

    done 
  done 
done 

}
###############################################################################
Slices()
{

###############################################################################
###############################################################################
if [ 1 -eq 0 ]; then
   cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_Reference/Fits
   cpEbyeMacros
   source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitPrepareData.sh 
   Slices 300 0.2 3.2 7 0.8 0.5 0.5 > outslice.log 2>&1 &
fi
###############################################################################
###############################################################################

# inputs
nSlice=$1
ptMin=$2
ptMax=$3
maxIter=$4
piSkew=$5
kaSkew=$6
prSkew=$7
piKurt=$8
kaKurt=$9
prKurt=${10}

# base directory
outbase=`pwd`

# create the fit results directory
fitResDir=$outbase/FitResults_piS$piSkew\_kaS$kaSkew\_prS$prSkew\_piK$piKurt\_kaK$kaKurt\_prK$prKurt\_$DATE
DirCheckDelete $fitResDir; mkdir $fitResDir; cd $fitResDir

# loop over files and apply the iterative fitting
for dataFile in $(ls -lS $hisDir | awk {'print$9'}); do  

   # from the file name retrieve the directory name "Hists_PbPb_eta_0.4_0.6_cent_20_25.root --> PbPb_eta_0.4_0.6_cent_10_15"
   A=$dataFile;  
   [[ $A == *"Hists_"* ]] || continue
   resultDirName="${A:6:${#A}-11}"
   echo $resultDirName
   cd $fitResDir
   for (( isign = 0; isign < ${#signArr[@]}; isign++ )); do
     sign=${signArr[$isign]}
     fitEtaCentDir=$fitResDir/$resultDirName\_sign$sign; DirCheckDelete $fitEtaCentDir;  mkdir $fitEtaCentDir;  cd $fitEtaCentDir
   
     fileName=$hisDir/$dataFile
     cp $shDir/rootlogon.C $shDir/PIDIterativeFitting*.C $shDir/submitPrepareData.sh .
     chmod +x submitPrepareData.sh
     [[ $testScripts == 1 ]] && continue 
   
     # Submit the jobs
     echo " =============== Real data fits are being performed =============== "
     if [ $localRunning == 1 ]; then
       ./submitPrepareData.sh processRealSlices $sign $fileName $nSlice $ptMin $ptMax $maxIter $piSkew $kaSkew $prSkew $piKurt $kaKurt $prKurt
     elif [ $localRunning == 0 ]; then
       eval $qsubCommand ./submitPrepareData.sh processRealSlices $sign $fileName $nSlice $ptMin $ptMax $maxIter $piSkew $kaSkew $prSkew $piKurt $kaKurt $prKurt
     fi
   done 

done

cd $outbase

}
###############################################################################
ParameterCheck()
{

###############################################################################
###############################################################################
if [ 1 -eq 0 ]; then
   cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_mombin20MeV/Fits/FitScan_20150111_02_12
   cpEbyeMacros
   source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitPrepareData.sh 
   ParameterCheck
fi
###############################################################################
###############################################################################

# base directory
outbase=`pwd`

# loop over files and apply the iterative fitting
for dirName in $(ls -lS $outbase | awk {'print$9'}); do  

# from the file name retrieve the directory name "Hists_PbPb_eta_0.4_0.6_cent_20_25.root --> PbPb_eta_0.4_0.6_cent_10_15"
   A=$dirName; cd $dirName
   cp $shDir/rootlogon.C $shDir/CreateAllParamPlots.C $shDir/submitPrepareData.sh .
   chmod +x submitPrepareData.sh
   echo $A
   
   [[ $testScripts == 1 ]] && continue
   
# Submit the jobs
   if [ $localRunning == 1 ]; then
     ./submitPrepareData.sh processParameterCheck 
   elif [ $localRunning == 0 ]; then
     eval $qsubCommand ./submitPrepareData.sh processParameterCheck
   fi
   
   cd $outbase
done

}
###############################################################################
processMCSlices()
{

# inputs
rootFileAllHists=$1
nSlice=$2
ptMin=$3
ptMax=$4
maxIter=$5

# Aliroot
# source  /lustre/nyx/alice/users/marsland/bin/.tpcdevEnv

# run aliroot for each file
aliroot -q -b rootlogon.C PIDIterativeFittingMC.C+\(\"$rootFileAllHists\",$nSlice\,$ptMin\,$ptMax\,$maxIter\)  2> err.log

}
###############################################################################
processRealSlices()
{

# inputs
sign=$1
rootFileAllHists=$2
nSlice=$3
ptMin=$4
ptMax=$5
maxIter=$6
piSkew=$7
kaSkew=$8
prSkew=$9
piKurt=${10}
kaKurt=${11}
prKurt=${12}

# Aliroot
# source  /lustre/nyx/alice/users/marsland/bin/.tpcdevEnv

# run aliroot for each file
aliroot -q -b rootlogon.C PIDIterativeFitting.C+\($sign,\"$rootFileAllHists\",$nSlice\,$ptMin\,$ptMax\,$maxIter\,$piSkew\,$kaSkew\,$prSkew,$piKurt\,$kaKurt\,$prKurt\)  2> err.log

}
###############################################################################
processParameterCheck()
{

treeFileName=AllPIDtrees.root
paramFileName=ParamTree.root

# check the file if exist delete
FileCheckDelete $treeFileName

# merge trees
hadd AllPIDtrees.root */Tree*.root

# produce paramtree
aliroot -q -l -b rootlogon.C CreateAllParamPlots.C+\(0\,\"AllPIDtrees.root\",\"\"\)  2> err.log

# run aliroot for each file
aliroot -q -l -b rootlogon.C CreateAllParamPlots.C+\(1\,\"\",\"ParamTree.root\"\)  2> err.log

}
SentMergeHists()
{
###############################################################################
###############################################################################
#  to run  -->   
if [ 1 -eq 0 ]; then
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics_cRows_80_16EtaBin_mombin20MeV_ITSOFF/mergedPeriods
  cpEbyeMacros
  source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitPrepareData.sh 
  SentMergeHists > out.log 2>&1 & 
fi
###############################################################################
###############################################################################
outbase=`pwd`

# define eta values
counter=0
for ifileName in $(less $idenFile ); do  
   
   counter=$(($counter+1))
   outName=$(echo "00$counter\_AnalysisResults_Hists.root")
   echo  $ifileName $outName 
   
   cp $shDir/submitPrepareData.sh  .; chmod +x submitPrepareData.sh
   [[ $testScripts == 1 ]] && continue
   
   # submit jobs
   if [ $localRunning == 1 ]; then
      ./submitPrepareData.sh mergeHists $ifileName $histName $outName 
   elif [ $localRunning == 0 ]; then
      eval $qsubCommand ./submitPrepareData.sh mergeHists $ifileName $histName $outName 
   fi
   
done 

cd $outbase

}
###############################################################################
mergeHists()
{

# inputs
fileName=$1
histName=$2
outName=$3


# run aliroot for each file
aliroot -q -b $shDir/mergeHists.C\(\"$fileName\",\"$histName\",\"$outName\"\) 

}

###############################################################################
DirCheckDelete()
{

dirName=$1
if [ -d "$dirName" ]; then
echo "AAAA  " $dirName "  it exist already lets delete it  "
rm -rf $dirName 
fi

}
###############################################################################
DirCheckCreate()
{

dirName=$1
if [ ! -d "$dirName" ]; then
echo " ========================================================================  " $dirName " create a new one  "
mkdir $dirName 
fi

}
###############################################################################
FileCheckDelete()
{

file=$1
if [ -f "$file" ]; then
echo " ========================================================================  " $file "  It will be deleted  "
rm -rf $file 
fi

}
###############################################################################
FileCheckExit()
{

file=$1
if [ -f "$file" ]; then
echo " ========================================================================  " $file "  It will be deleted  "
return 1 
fi

}
###############################################################################
main()
{
eval "$@"
}

main $@





