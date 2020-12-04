class  AliAnalysisManager;
class  AliAnalysisAlien;
//
//
//
// modes: "test" to run over a small set of files (requires alien connection but everything stored locally),
//        "full" to run over everything on the grid,
//        "terminate" to merge results after "full"
//
/*

//
// valgrindOption --> 0 --> Normal, 1--> valgrind, 2--> callgrind, 3-->Massif
// mode           --> "test" or "full" or "terminate"
// localOrGrid    --> 0 --> Use only one run and run locally, 1 --> run for all runs on the grid
// list           --> defines the data set according to year, period and pass --> "test-2015-LHC15o-pass5_lowIR.list"
// fname          --> output directory name in my home folder in alien
// isMC           --> 0:Data , 1:MCfull,  2:fastMC
// setType        --> setting type encoded in Config file
// lhcPeriod      --> 1 for RUN1, and 2 for RUN2
//

Notes:
runs-2017-LHC17i2-pass1.list  --> MultSelectionTask should be disabled

How to run:

cd $RUN_ON_GRID_DIR/Ebye/test/testMC
rm *.xml Task* *.root stderr out* Ali*
cp $RUN_ON_GRID_DIR/Ebye/code/*.*  .; rm *.so *.d


// run Massif: 0 --> Normal, 1--> valgrind, 2--> callgrind, 3-->Massif
// ----------------------------------------------------------------------------------------------------------------------------------------
// tests real data
aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runsTest-2015-LHC15o-pass1.list"                           ,"Skimmed_ESDs_dEdxPerformance_test",0,0,2,"vAN-20190820-1")'
aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/tiden_Filter_15o_pass1/runsSublist9-2015-LHC15o-pass1.list","Skimmed_ESDs_dEdxPerformance"     ,0,0,2,"vAN-20190820-1")'

aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/tiden_Filter_15o_pass1/runsBigList1-2015-LHC15o-pass1.list","Skimmed_ESDs_dEdxPerformance"     ,0,0,2,"vAN-20190820-1")'
aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/tiden_Filter_15o_pass1/runsBigList2-2015-LHC15o-pass1.list","Skimmed_ESDs_dEdxPerformance"     ,0,0,2,"vAN-20190820-1")'

//
// for testing only
aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runsTest-2015-LHC15o-pass1.list","Skimmed_ESDs_dEdxPerformance_test",0,0,2,"vAN-20190820-1")'
aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runsMonalisa-2018-LHC18q-pass1.list","Skimmed_ESDs_dEdxPerformance_test",0,0,3,"vAN-20190820-1")'
aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runsMonalisa-2018-LHC18r-pass1.list","Skimmed_ESDs_dEdxPerformance_test",0,0,3,"vAN-20190820-1")'
//
//
aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/ESD_Filtering/lists/runsTest-2016-LHC16g1-pass1.list","Skimmed_ESDs_dEdxPerformance_test",1,0,2,"vAN-20190820-1")'
aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/ESD_Filtering/lists/runs-2016-LHC16g1-pass1.list","Skimmed_ESDs_dEdxPerformance",1,0,2,"vAN-20190820-1")'

aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runs5Large-2018-LHC18q-pass1.list","Skimmed_ESDs_dEdxPerformance"     ,0,0,3,"vAN-20190820-1")'
aliroot -b -q 'runTiden_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runs5Large-2018-LHC18r-pass1.list","Skimmed_ESDs_dEdxPerformance"     ,0,0,3,"vAN-20190820-1")'



Problems faced:
saving --> delete the output dir on alien /alice/cern.ch/user/m/marsland/EbyeIterPID/LHC15o_pass5_lowIR

###
// command to check valgrind output
less xxx.txt | grep TIdentity | awk {'print$5'} | sort | uniq -u

###
// count chunks
source $RUN_ON_GRID_DIR/Ebye/code/helpers/alihadd_GRIDoutput.sh;
CountFilesRealData "/alice/data/2015/LHC15o" "pass1" "AliESDs.root"

###
// Post process the counting of chunks
sort -nk5 chunkStat_2015o_pass1.list > chunkStat_2015o_pass1_sorted.list
cat chunkStat_2015o_pass1_sorted.list | awk '{sum+=$5 ; print $0} END{print "sum=",sum}'

#### Copy and merge afterwards:
Copy:

cd /lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/ATO-465/data/LHC16g1/LHC16g1_pass1_140619
source /u/marsland/PHD/macros/marsland_EbyeRatios/mergeFilteredTreesAtGSI.sh;
copyFilteredTrees /alice/cern.ch/user/p/pwg_pp/Skimmed_ESDs_16g1_pass1_140619/

Merge:

filtTreeDirGSI=/lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/ATO-465/data/LHC16g1c/LHC16g1c_pass1_140519
cd $filtTreeDirGSI
source /u/marsland/PHD/macros/marsland_EbyeRatios/mergeFilteredTreesAtGSI.sh;
mergeFilteredTreesAtGSI



*/

Bool_t fUseMultSelection = kTRUE;
const Int_t nTestFiles = 1;
Int_t nChunksPerJob = 10;
Bool_t fRunLocalFiles = kFALSE;
TString aliPhysicsTag = "";

TString fValgrind  = "/usr/bin/valgrind --leak-check=full --leak-resolution=high --num-callers=40 --error-limit=no --show-reachable=yes  --log-file=xxx.txt --suppressions=$ROOTSYS/etc/valgrind-root.supp  -v ";
TString fCallgrind = "/usr/bin/valgrind --tool=callgrind --log-file=cpu.txt   --num-callers=40 -v  --trace-children=yes ";
TString fMassif    = "/usr/bin/valgrind --tool=massif ";

void runTiden_FilteredRaw(Int_t valgrindOption = 0, TString mode="test",Int_t localOrGrid=0, TString passStr="1", TString list = "$RUN_ON_GRID_DIR/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list", TString fname="EbyeIterPID", Int_t isMC=0, Int_t setType=0, Int_t lhcPeriod=2, TString physicsTagForFullTest="vAN-20190619-1")
{
  //
  // valgrindOption --> 0 --> Normal, 1--> valgrind, 2--> callgrind, 3-->Massif
  // mode           --> "test" or "full" or "terminate"
  // localOrGrid    --> 0 --> Use only one run and run locally, 1 --> run for all runs on the grid
  // list           --> defines the data set according to year, period and pass --> "test-2015-LHC15o-pass5_lowIR.list"
  // fname          --> output directory name in my home folder in alien
  // isMC           --> 0:Data , 1:MCfull,  2:fastMC
  // setType        --> setting type encoded in Config file
  // lhcPeriod      --> 1 for RUN1, and 2 for RUN2
  //
  if (isMC==0) nChunksPerJob=30;

  aliPhysicsTag=physicsTagForFullTest;
  AliLog::SetGlobalDebugLevel(0);
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  //
  // Create and configure the alien handler plugin
  //
  AliAnalysisGrid *alienHandler;
  if (!fRunLocalFiles){
    alienHandler = CreateAlienHandler(valgrindOption,mode,localOrGrid,list,fname,isMC);
    if (!alienHandler) return;
    // Connect plug-in to the analysis manager
    mgr->SetGridHandler(alienHandler);
    //   AliAODInputHandler* aodH = new AliAODInputHandler();
    //   AliESDInputHandler* aodH = new AliESDInputHandler();
    //   mgr->SetInputEventHandler(aodH)
  }
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  // Add handlers
  AliVEventHandler* handler=0x0;
  if (isMC>0){
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C"));
    handler = AddMCHandler(kFALSE);
    ((AliMCEventHandler*)handler)->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C"));
    handler = AddESDHandler();
    ((AliESDInputHandler*)handler)->SetReadFriends(kFALSE);
    ((AliESDInputHandler*)handler)->SetNeedField();
  } else
  {
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C"));
    handler = AddESDHandler();
  }
  mgr->SetInputEventHandler(handler);
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  if (!(isMC==2 || isMC==4)){
    // Add Additional tasks
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"));
    AddTaskPhysicsSelection(isMC);
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
    AliAnalysisTaskPIDResponse *taskPID = NULL;
    // if (isMC==0) taskPID=AddTaskPIDResponse(kFALSE,kTRUE,kTRUE,1);
    if (isMC==0) taskPID=AddTaskPIDResponse(kFALSE,kTRUE,kFALSE,passStr);
    if (isMC==1 || isMC==3) taskPID=AddTaskPIDResponse(kTRUE,kTRUE,kTRUE,passStr);
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C"));
    AliCentralitySelectionTask *taskCentrality=AddTaskCentrality();
    if(fUseMultSelection){
      gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
      AddTaskMultSelection();
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskConfigOCDB.C");
      AddTaskConfigOCDB("raw://");
    }
  }
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  // to run a local task not in AliPhysics (there will be conflicts if a task in AliPhysics has the same name)
  //   gROOT->LoadMacro("AliAnalysisTaskNetLambdaMC.cxx+g");
  //   gROOT->LoadMacro("AddTaskNetLambdaMC.C");
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  // to run a task in AliPhysics
  //   gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/EBYE/IdentityMethodEbyeFluctuations/macros/AddTask_marsland_EbyeIterPID.C");
  //   AliAnalysisTask *ana = AddTask_marsland_EbyeIterPID(kTRUE,"Config_marsland_EbyeIterPID.C",setType);
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  // to run locally
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  //
  // Filtered tree
  gROOT->LoadMacro("./AddTaskFilteredTreeLocal.C");
  AliAnalysisTask *ana = AddTaskFilteredTree("",isMC);
  //
  // My task
  gROOT->LoadMacro("$RUN_ON_GRID_DIR/Ebye/code/AliAnalysisTaskTIdentityPID.cxx+g");
  gROOT->LoadMacro("$RUN_ON_GRID_DIR/Ebye/code/AddTask_marsland_TIdentityPID.C");
  AliAnalysisTask *ana = AddTask_marsland_TIdentityPID(kFALSE,"Config_marsland_TIdentityPID.C",setType,lhcPeriod);
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if (!fRunLocalFiles) {
    // Start analysis in grid.
    mgr->StartAnalysis("grid");   // to set the number of events --> mgr->StartAnalysis("grid",nEvents);
  } else {
    // to run over files stored locally, uncomment this section,
    // and comment out the above lines related to alienHandler and StartAnalysis("grid")
    TChain *chain = new TChain("esdTree");
    chain->AddFile("$RUN_ON_GRID_DIR/Ebye/test/testReal/testFiles_00246272/AliESDs_0.root");
    chain->AddFile("$RUN_ON_GRID_DIR/Ebye/test/testReal/testFiles_00246272/AliESDs_1.root");
    chain->AddFile("$RUN_ON_GRID_DIR/Ebye/test/testReal/testFiles_00246272/AliESDs_2.root");
    chain->Print();
    mgr->StartAnalysis("local",chain);
  }

}

// ----------------------------------------------------------------------------------------------------------------------
AliAnalysisGrid* CreateAlienHandler(Int_t valgrindOption = 0,TString mode="test",Int_t localOrGrid=0,TString list = "$RUN_ON_GRID_DIR/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list",TString fname="testName",Int_t isMC=0)
{

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  switch (valgrindOption) {
    case 1: plugin->SetExecutableCommand(Form("%s aliroot -q -b",fValgrind.Data())); break;
    case 2: plugin->SetExecutableCommand(Form("%s aliroot -q -b",fCallgrind.Data())); break;
    case 3: plugin->SetExecutableCommand(Form("%s aliroot -q -b",fMassif.Data())); break;
    default: plugin->SetExecutableCommand("aliroot -q -b"); break;
  }
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(nTestFiles);
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion(aliPhysicsTag); // change to something up-to-date vAN-20170717-1 vAN-20180403-1
  plugin->SetNrunsPerMaster(1);
  plugin->SetSplitMaxInputFileNumber(nChunksPerJob); // 3 in the LEGO trains
  if (!isMC) plugin->SetRunPrefix("000");  // fort the data
  //
  // -----------------------------------------------------------------------------------------
  // ------------------------- Read runs from the list----------------------------------------
  // -----------------------------------------------------------------------------------------
  //
  // Read runs from the list
  gSystem->Exec(Form("grep -v run  %s | awk '{print $1}' > tmp.list",list.Data()));
  FILE *fileTmp=fopen("tmp.list","r");
  Int_t nRuns=0, run=0, runPrev=0;
  std::cout << " ----------- list of good runs ----------- " << std::endl;
  while (!feof(fileTmp)) {
    fscanf(fileTmp,"%d",&run);
    if (runPrev != run) runPrev=run;
    else continue;
    plugin->AddRunNumber(run);
    std::cout << nRuns+1 << "  " << run << "  is included in the processing "<< std::endl;
    nRuns++;
    if ( localOrGrid==0 && nRuns==1 ) {
      std::cout << " in test mode process only one run " << std::endl;
      break;
    }
  }
  gSystem->Exec("rm tmp.list");
  TObjArray *objArr1 = list.Tokenize("/");
  TString fileName = ((objArr1->At((Int_t)objArr1->GetLast()))->GetName());
  TObjArray *objArr2  = fileName.Tokenize("-");
  Int_t year         = atoi((objArr2->At(1))->GetName());
  TString yearStr    = (objArr2->At(1))->GetName();
  TString period     = ((objArr2->At(2))->GetName());
  TString passLong   = ((objArr2->At(3))->GetName());
  TObjArray *objArr3  = passLong.Tokenize(".");
  TString pass = ((objArr3->At(0))->GetName());
  std::cout << " --------------------------------------- " << std::endl;
  std::cout << " year = " << year << "   period = " << period << "   pass = " << pass << std::endl;
  std::cout << " --------------------------------------- " << std::endl;
  //
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  //
  // Set filenames, input and output directories on alien
  TDatime date;
  plugin->SetGridWorkingDir(Form("%s/%s_%s_%d_%d%d/",fname.Data(),period.Data(),pass.Data(),date.GetDate(),date.GetHour(),date.GetMinute()));
  if (isMC==0) {       // data
    std::cout << " Data SOURCE = REAL DATA " << std::endl;
    plugin->SetGridDataDir(Form("/alice/data/%d/%s/",year,period.Data())); // /alice/data/2015/LHC15o/000246858/pass1/15000246858039.402/AliESDs.root
    plugin->SetDataPattern(Form("/%s/*/AliESDs.root",pass.Data()));
  }
  else if (isMC==1) {  // RUN2 full MC gen+rec
    std::cout << " Data SOURCE = RUN2 full MC gen+rec " << std::endl;
    plugin->SetGridDataDir(Form("/alice/sim/%d/%s/",year,period.Data()));
    plugin->SetDataPattern("/*/AliESDs.root");
  }
  else if (isMC==2) {  // RUN2 fast MC gen
    std::cout << " Data SOURCE = RUN2 fast MC gen " << std::endl;
    plugin->SetAdditionalLibs("pythia6 Tree Geom VMC Physics Minuit Gui Minuit2 STEERBase ESD OADB ANALYSIS ANALYSISalice CDB STEER CORRFW EMCALUtils EMCALrec VZERObase VZEROrec");
    plugin->SetAdditionalRootLibs("libVMC.so libPhysics.so libTree.so libMinuit.so libProof.so libSTEERBase.so libESD.so libAOD.so");
    plugin->SetMCLoop(kTRUE);
    plugin->SetUseMCchain();
    plugin->SetNMCjobs(1000);
    plugin->SetNMCevents(100);
    plugin->SetGridDataDir(Form("/alice/sim/%d/%s/",year,period.Data()));
    plugin->SetDataPattern("/*/galice.root");
    //       plugin->SetDataPattern("/*/root_archive.zip#galice.root");
  }
  else if (isMC==3) {  // RUN1 full MC gen+rec HIJING
    std::cout << " Data SOURCE = RUN1 full MC gen+rec HIJING " << std::endl;
    plugin->SetGridDataDir(Form("/alice/sim/%s/",period.Data()));
    plugin->SetDataPattern("/*/AliESDs.root");
  }
  else if (isMC==4) {  // RUN1 fast MC gen
    std::cout << " Data SOURCE = RUN1 fast MC gen " << std::endl;
    plugin->SetAdditionalLibs("pythia6 Tree Geom VMC Physics Minuit Gui Minuit2 STEERBase ESD OADB ANALYSIS ANALYSISalice CDB STEER CORRFW EMCALUtils EMCALrec VZERObase VZEROrec");
    plugin->SetAdditionalRootLibs("libVMC.so libPhysics.so libTree.so libMinuit.so libProof.so libSTEERBase.so libESD.so libAOD.so");
    plugin->SetMCLoop(kTRUE);
    plugin->SetUseMCchain();
    plugin->SetNMCjobs(1000);
    plugin->SetNMCevents(100);
    //       plugin->SetSplitMode(Form("production:1-%d", 100));
    plugin->SetGridDataDir(Form("/alice/sim/%d/%s/",year,period.Data()));
    plugin->SetDataPattern("/*/galice.root");
  }

  plugin->SetAnalysisMacro(Form("TaskEbyeIterPIDMC_%s_%s.C",period.Data(),pass.Data()));
  plugin->SetExecutable(Form("TaskEbyeIterPIDMC_%s_%s.sh",period.Data(),pass.Data()));
  plugin->SetJDLName(Form("TaskEbyeIterPIDMC_%s_%s.jdl",period.Data(),pass.Data()));
  //
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // include additional libs
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  // to run locally
  plugin->SetAnalysisSource("AliAnalysisTaskTIdentityPID.cxx");
  plugin->SetAdditionalLibs("AliAnalysisTaskTIdentityPID.cxx AliAnalysisTaskTIdentityPID.h");
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  plugin->SetGridOutputDir(yearStr); // In this case will be $HOME/work/output

  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  //plugin->SetDefaultOutputs(0);
  //plugin->SetOutputFiles("AnalysisResults.root");
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  //plugin->SetMaxMergeFiles(40);
  plugin->SetMaxMergeStages(4);

  plugin->SetTTL(56399);
  // plugin->SetTTL(86399);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);
  // Optionally modify split mode (default 'se')
  plugin->SetSplitMode("se");
  //plugin->SetSplitMaxInputFileNumber();
  plugin->SetKeepLogs(kFALSE);
  plugin->SetOutputToRunNo(1);
  // my settings
  // if (localOrGrid==0) plugin->SetKeepLogs(kTRUE); // keep the log files
  plugin->SetKeepLogs(kTRUE); // keep the log files

  return plugin;
}
