#!/bin/bash

#  DirCheck  -> to check if directory exist
#  FileCheck -> to check if file exist
#

# some global variables
# DATE=`date +%Y_%m_%d_%H_%M_%S`
# dataType  : MC(fTree):0, Real(fTree):1
alirootVersion=ALICE/vAN-20141214-prod
DATE=`date +%Y_%m_%d_%H_%M`
qsubCommand="qsub -cwd -V -l h_rt=24:0:0,h_rss=6G -b y -r y -o out.log -e err.log"
# qsubCommand="qsub -cwd -V -l h_rt=24:0:0,h_rss=4G -P alice -b y -r y -o out.log -e err.log"
IdBase=/lustre/nyx/alice/users/marsland/pFluct/files/analysis/TIdentity/example
shDir=/lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts
chmod +x $shDir/submitIdentity.sh

############################## Modification region ############################
###############################################################################

localRunning=0       # 1 for testing locally, 0 for running on batch farm
maxSubSampleIndex=1 # if 1 only run over full stat do not process subsamples (default 26)

############ real data setting ############
dataSetLocation=PbPb/Real/Systematics_cRows_80_16EtaBin_mombin20MeV
paramTree=$dataSetLocation/ParamTrees/FitResults_piS0.7_kaS0.5_prS0.5_pikaprKauto_kaK1.04_FitErrors.root
dataTreeDir=$dataSetLocation/SubSamples
iterArray=(6)
pDownArray=(0.2 0.2 0.3 0.3)
pUpArray=(  1.5 1.8 1.5 1.8)
etaDownArray=(-0.8 -0.5)
etaUpArray=(   0.8  0.5)
# pDownArray=(0.2)
# pUpArray=(  1.5)
# etaDownArray=(-0.1 -0.15 -0.2 -0.25 -0.3 -0.35 -0.4 -0.45 -0.5 -0.55 -0.6 -0.65 -0.7 -0.75 -0.8)
# etaUpArray=(   0.1  0.15  0.2  0.25  0.3  0.35  0.4  0.45  0.5  0.55  0.6  0.65  0.7  0.75  0.8)
# etaDownArray=(-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8)
# etaUpArray=(   0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8)
centDownArray=(0.0 5.0  10.0 20.0 30.0 40.0 50.0 60.0 70.0)
centUpArray=(  5.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0)
# centDownArray=(0.0 2.5 5.0  7.5  10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0 75.0)
# centUpArray=(  2.5 5.0 7.5  10.0 15.0 20.0 25.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 65.0 70.0 75.0 80.0) 
############################################

############## MC data setting #############
# dataSetLocation=PbPb/MC/MC_Closure_Test_HIJING
# paramTree=$dataSetLocation/ParamTrees/FitResults_piS0.7_kaS0.5_prS0.5_pikaprKauto.root
# dataTreeDir=$dataSetLocation/SubSamples
# iterArray=(3)
# pDownArray=(0.2)
# pUpArray=(  1.5)
# # pDownArray=(0.2 0.2 0.3 0.3)
# # pUpArray=(  1.5 1.8 1.5 1.8)
# etaDownArray=(-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8)
# etaUpArray=(   0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8)
# # etaDownArray=(-0.5 -0.8)
# # etaUpArray=(   0.5  0.8)
# centDownArray=(0.0 5.0  10.0 20.0 30.0 40.0 50.0 60.0 70.0)
# centUpArray=(  5.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0)
############################################

########################## end of Modification region #########################
###############################################################################

outbase=/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data

echo " ============= All functions are sourced, have FUN !!!! ============= "

## make the executable "testiden"
if [ 1 -eq 0 ]; then
  cpEbyeMacros
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/TIdentity
  make clean; make
  cpEbyeMacros
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/TIdentity/example
  make clean; make 
fi
###############################################################################
############################  MAIN FUNCTIONS ##################################
###############################################################################

IdenFullScan()
{
###############################################################################
###############################################################################
# Send all jobs to farm;  loop over all centrality bins
#  to run  -->   
if [ 1 -eq 0 ]; then
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_cRows/TIdenResults
  cpEbyeMacros
  source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitIdentity.sh 
  IdenFullScan iter  > out.log 2>&1 &  
fi
###############################################################################
###############################################################################

#inputs
loopName=$1

curDir=`pwd`
etaMomScanDir=$curDir/TIdenResults_EtaMomScan_$DATE
DirCheckDelete $etaMomScanDir; mkdir $etaMomScanDir; cd $etaMomScanDir
cp $outbase/$paramTree $etaMomScanDir

# Momentum Scan
for (( ieta = 0; ieta < ${#etaDownArray[@]}; ieta++ )); do
  for (( ip = 0; ip < ${#pDownArray[@]}; ip++ )); do
    
     etaDown=${etaDownArray[$ieta]}
     etaUp=${etaUpArray[ieta]}
     pDown=${pDownArray[$ip]}
     pUp=${pUpArray[$ip]}
      
      # create a dir for each scenario
     etaMomDir=$etaMomScanDir/Eta_$etaDown\_$etaUp\_Mom_$pDown\_$pUp
     if [ -d "$etaMomDir" ]; then
       continue
     fi
     DirCheckDelete $etaMomDir; mkdir $etaMomDir; cd $etaMomDir
     echo $loopName $etaDown $etaUp $pDown $pUp
     IdenIterScan $loopName $etaDown $etaUp $pDown $pUp
      
  done
done

cd $curDir

}
###############################################################################
IdenIterScan()
{
###############################################################################
###############################################################################
# Send all jobs to farm;  loop over all centrality bins
#  to run  -->   
if [ 1 -eq 0 ]; then
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/TIdenResults_PbPb/MC/DataTrees_IdMCtracks_OerderedEvents
  cpEbyeMacros
  source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitIdentity.sh 
  IdenIterScan iter 0 0.2 0.3 0.35  > out.log 2>&1 &  
fi
###############################################################################
###############################################################################

#inputs
loopName=$1
etaDown=$2
etaUp=$3
pDown=$4
pUp=$5

curDir=`pwd`
for (( iter = 0; iter < ${#iterArray[@]}; iter++ )); do
   iterIndex=${iterArray[$iter]}
   echo $iterIndex $loopName $etaDown $etaUp $pDown $pUp
   cp $shDir/submitIdentity.sh  . 
   chmod +x submitIdentity.sh
   ./submitIdentity.sh SentIdenJobs $iterIndex $loopName $etaDown $etaUp $pDown $pUp   
done
cd $curDir

}
###############################################################################
SentIdenJobs()
{
###############################################################################
###############################################################################
# Send all jobs to farm;  loop over all centrality bins
#  to run  -->   
if [ 1 -eq 0 ]; then
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/TIdenResults_PbPb/MC/DataTrees
  cpEbyeMacros
  source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitIdentity.sh 
  SentIdenJobs 1 iter 0 0.2 0.3 0.35           > out.log 2>&1 &  
fi
###############################################################################
###############################################################################

iterIndex=$1           # iter number of the fits 
loopName=$2
etaDown=$3
etaUp=$4
pDown=$5
pUp=$6

curDir=`pwd`

# directory for hists and samples 
IdResDir=$curDir/TIdenResults_$loopName\_$iterIndex
DirCheckDelete $IdResDir
mkdir $IdResDir
cp $outbase/$paramTree $IdResDir

##===========================================================
# loop over centrality bins
for icentDown in ${centDownArray[@]}
   do
   
   ## if PP is processed there is only one cent bin
   if [ ! -d "$outbase/$dataTreeDir/cent_$icentDown" ]; then
   continue
   fi
   
   DirCheckDelete $IdResDir/cent_$icentDown
   mkdir $IdResDir/cent_$icentDown
   
   # get into dataTreeDir and read the file name
   cd $outbase/$dataTreeDir/cent_$icentDown
   
   fileIndex=-1
   for subSampleFile in $(ls -lS *.root | awk {'print$9'}) 
      do
      fileIndex=$(($fileIndex+1))
      [[ $fileIndex -ge $maxSubSampleIndex ]] && continue
      #   from the file name retrieve the directory name "aaa.root --> aaa" and create it
      A=$subSampleFile
      subSampleDir="${A:0:${#A}-5}"     # 0 character from front and 5 from end
      
      DirCheckDelete $IdResDir/cent_$icentDown/$subSampleDir
      mkdir $IdResDir/cent_$icentDown/$subSampleDir
      cd $IdResDir/cent_$icentDown/$subSampleDir
   
      # run the real macro       
      echo "running ===   " testId $dataTreeDir/cent_$icentDown/$subSampleFile $icentDown $fileIndex $iterIndex $etaDown $etaUp $pDown $pUp
      cp $shDir/submitIdentity.sh . ; chmod +x submitIdentity.sh
      ./submitIdentity.sh testId $dataTreeDir/cent_$icentDown/$subSampleFile $icentDown $fileIndex $iterIndex $etaDown $etaUp $pDown $pUp
#       testId pp/Real/pp_cRows_80_16EtaBin_mombin20MeV_LooseCuts/SubSamples/cent_-105/SubSample_cent-105_ss1.root -105 1 6 0 0.2 0.2 0.3 
#       testId pp/Real/pp_cRows_80_16EtaBin_mombin20MeV_LooseCuts/SubSamples/cent_-105/SubSample_cent-105_ss0.root -105 0 6 -0.8 0.8 0.4 0.6
#       testId pp/Real/pp_cRows_80_16EtaBin_mombin20MeV_LooseCuts/SubSamples/cent_-105/SubSample_cent-105_ss6.root -105 4 6 -0.8 0.8 0.4 0.6
   done

done

cd $curDir

}
###############################################################################
testId()
{

###############################################################################
###############################################################################
# to run  -->   
if [ 1 -eq 0 ]; then   
   cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/test
   cpEbyeMacros
   source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitIdentity.sh 
   testId PbPb/Real/PbPb_cRows80/SubSamples/cent_0/SubSample_cent0_ss1.root 0 1 1 0 0.2 0.2 2   
   testId PbPb/MC/MC_PbPb_cRows/SubSamples/cent_0/SubSample_cent0_ss1.root 0 1 1 0 0.2 0.2 2   
   testId PbPb/MC/MC_PbPb_Reference/SubSamples/cent_0/SubSample_cent0_ss1.root 0 1 1 0 0.2 0.3 0.35   
fi
###############################################################################
###############################################################################

# inputs
dataTree=$1
cent=$2
subSampleIndex=$3
iterIndex=$4
etaDown=$5
etaUp=$6
pDown=$7
pUp=$8

# echo "$dataTree -- $paramTree -- $cent -- $subSampleIndex -- $iterIndex -- $etaDown -- $etaUp -- $pDown -- $pUp"

if [ $localRunning == 1 ]; then
   $IdBase/testIden $dataTree $paramTree $cent $subSampleIndex $iterIndex $etaDown $etaUp $pDown $pUp
elif [ $localRunning == 0 ]; then
   eval $qsubCommand $IdBase/testIden $dataTree $paramTree $cent $subSampleIndex $iterIndex $etaDown $etaUp $pDown $pUp
fi

}
###############################################################################
MakeFinalPlots()
{
###############################################################################
###############################################################################
#  to run  -->   
if [ 1 -eq 0 ]; then
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/TIdenResults/TIdenResults_EtaMomScan_2015_03_14_13_10
  cpEbyeMacros
  source /lustre/nyx/alice/users/marsland/pFluct/files/analysis/scripts/submitIdentity.sh 
  MakeFinalPlots 6 > out.log 2>&1 & 
fi
###############################################################################
###############################################################################
outbase=`pwd`
fitIterIndex=$1

# loop over centralities
for name in $(ls -lS $outbase | awk {'print$9'}); do 

   substr="root"
   [[ "$name" == *"$substr"* ]] && continue
   macroDir=/u/marsland/PHD/macros/marsland_EbyeRatios
   dirName=$outbase/$name
   cd $dirName; echo $dirName
   rm TIdenMoments.root
   hadd TIdenMoments.root TI*/cen*/*/TImoment*.root 
   aliroot -q -b rootlogon.C $macroDir/TIdentityStatisticalError.C+\(\"TIdenMoments.root\",$fitIterIndex\) 
   
done 
cd $outbase

}
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





