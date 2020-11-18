class  AliAnalysisManager;
class  AliAnalysisAlien;

/*

JIRA: PWGPP-340, ATO-465

--> create a list from monalisa list format
less runs-2010-LHC10h-pass2_monalisa.list   | sed -e $'s/,\ /\\\n/g'  >> runs-2010-LHC10h-pass2.list
less runs-2011-LHC11h_2-pass2_monalisa.list | sed -e $'s/,\ /\\\n/g'  >> runs-2011-LHC11h_2-pass2.list
less runs-2015-LHC15o-pass1_monalisa.list   | sed -e $'s/,\ /\\\n/g'  >> runs-2015-LHC15o-pass1.list

How to run
cp $RUN_ON_GRID_DIR/Ebye/code/*.*  . ; rm *.d *.so


// 10h
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"2","$RUN_ON_GRID_DIR/Ebye/lists/runs-2010-LHC10h-pass2.list","Skimmed_ESDs_10h_pass2_140519",0, "vAN-20190902-1")'
// 11h
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"2","$RUN_ON_GRID_DIR/Ebye/lists/runs-2011-LHC11h_2-pass2.list","Skimmed_ESDs_11h_pass2_140519",0, "vAN-20190902-1")'
// 15o
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runsTest-2015-LHC15o-pass1.list","Skimmed_ESDs_dEdxPerformance_test",0, "vAN-20190902-1")'
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runs-2015-LHC15o-pass1.list","Skimmed_ESDs_dEdxPerformance",0, "vAN-20190902-1")'
// 18q
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runsMonalisa-2018-LHC18q-pass1.list","Skimmed_ESDs",0, "vAN-20190902-1")'
// 18r
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runsMonalisa-2018-LHC18r-pass1.list","Skimmed_ESDs",0, "vAN-20190902-1")'
// 16g1a
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runs-2016-LHC16g1-pass1.list","Skimmed_ESDs_16g1_pass1_140619",1,"vAN-20190902-1")'
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runs-2016-LHC16g1a-pass1.list","Skimmed_ESDs_16g1a_pass1_140619",1,"vAN-20190902-1")'
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runs-2016-LHC16g1b-pass1.list","Skimmed_ESDs_16g1b_pass1_140619",1,"vAN-20190902-1")'
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runs-2016-LHC16g1c-pass1.list","Skimmed_ESDs_16g1c_pass1_140619",1,"vAN-20190902-1")'
//
//
//
// LHC18l8a -->  	Pb-Pb, 5.02 TeV - General-purpose Monte Carlo production anchored to LHC18q/r, minimum bias, ALIROOT-8157
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runsMonalisa5Runs-2018-LHC18l8a-pass1.list","Skimmed_ESDs",1,"vAN-20190902-1")'
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runsMonalisa5Runs-2018-LHC18q-pass1.list","Skimmed_ESDs",0, "vAN-20190902-1")'
//
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runs17Large-2016-LHC16g1-pass1.list","Skimmed_ESDs",1,"vAN-20190902-1")'
aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"1","$RUN_ON_GRID_DIR/Ebye/lists/runs5Runs-2015-LHC15o-pass1.list","Skimmed_ESDs",0, "vAN-20190902-1")'


// ------------------------------------- HOW to USE --------------------------------------------------------------------

Input Parameters:
  valgrindOption --> 0 --> Normal, 1--> valgrind, 2--> callgrind, 3-->Massif
  mode           --> "test" to run over a small set of files (requires alien connection but everything stored locally),
                     "full" to run over everything on the grid,
                     "terminate" to merge results after "full"
  localOrGrid    --> 0 ; Use only one run and run locally,
                     1 ; run for all runs on the grid
  list           --> defines the data set according to year, period and pass --> "test-2015-LHC15o-pass5_lowIR.list"
  fname          --> output directory name in my home folder in alien
  isMC           --> 0:Data , 1:MCfull

Usage:
  TEST MODE Example --> to run over "nTestFiles = 5" locally
      aliroot -b -q 'runGrid_FilteredRaw.C(0,"test",0,"<path to run list>","<outputDir name on alien>",0)'

  FULL MODE example --> run over full run list on the grid  with "nChunksPerJob = 150" and "aliPhysicsTag = "vAN-20170610-1""
      aliroot -b -q 'runGrid_FilteredRaw.C(0,"full",1,"<path to run list>","<outputDir name on alien>",0)'

// ----------------------------------------------------------------------------------------------------------------------

*/

//
// ----------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------
//
TString aliPhysicsTag = "vAN-20190902-1"; //  	vAN-20180828-1  vAN-20181119-1  vAN-20190105_ROOT6-1
Int_t nTestFiles = 1;
Int_t nChunksPerJob = 100;
TString fValgrind  = "/usr/bin/valgrind --leak-check=full --leak-resolution=high --num-callers=40 --error-limit=no --show-reachable=yes  --log-file=xxx.txt --suppressions=$ROOTSYS/etc/valgrind-root.supp  -v ";
TString fCallgrind = "/usr/bin/valgrind --tool=callgrind --log-file=cpu.txt   --num-callers=40 -v  --trace-children=yes ";
TString fMassif    = "/usr/bin/valgrind --tool=massif ";
//
// ----------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------
//
// ----------------------------------------------------------------------------------------------------------------------
void runGrid_FilteredRaw(Int_t valgrindOption = 0, TString mode="test", Int_t localOrGrid=0, TString passStr="1", TString list = "$RUN_ON_GRID_DIR/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list", TString fname="ESD_Filtering", Int_t isMC=kFALSE, TString physicsTagForFullTest="vAN-20190902-1")
{
  //
  // valgrindOption --> 0 --> Normal, 1--> valgrind, 2--> callgrind, 3-->Massif
  // mode           --> "test" or "full" or "terminate"
  // localOrGrid    --> 0 --> Use only one run and run locally, 1 --> run for all runs on the grid
  // list           --> defines the data set according to year, period and pass --> "test-2015-LHC15o-pass5_lowIR.list"
  // fname          --> output directory name in my home folder in alien
  // isMC           --> 0:Data , 1:MCfull
  //

  aliPhysicsTag=physicsTagForFullTest;

  // load libraries
  gSystem->Load("libCore.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGTools.so");
  gSystem->Load("libPWGCFCorrelationsBase.so");
  gSystem->Load("libPWGCFCorrelationsDPhi.so");
  AliLog::SetGlobalDebugLevel(5);
  //
  // Create and configure the alien handler plugin
  //
  AliAnalysisGrid *alienHandler = CreateAlienHandler(valgrindOption,mode,localOrGrid,list,fname,isMC);
  if (!alienHandler) return;
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  // Add handlers
  AliVEventHandler* handler=0x0;
  if (isMC>0){
    nChunksPerJob = nChunksPerJob*2;
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C"));
    handler = AddMCHandler(kFALSE);
    ((AliMCEventHandler*)handler)->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C"));
    handler = AddESDHandler();
    ((AliESDInputHandler*)handler)->SetReadFriends(kFALSE);
    ((AliESDInputHandler*)handler)->SetNeedField();
  } else {
    gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C"));
    handler = AddESDHandler();
  }
  mgr->SetInputEventHandler(handler);
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  // Add Additional tasks
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(isMC);
  //
  AliAnalysisTaskPIDResponse *taskPID=NULL;
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  if (isMC==0) taskPID=AddTaskPIDResponse(kFALSE,kTRUE,kFALSE,passStr);
  if (isMC>0)  taskPID=AddTaskPIDResponse(kTRUE,kTRUE,kTRUE,passStr);
  //
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality=AddTaskCentrality();
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AddTaskMultSelection();
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/AddTaskConfigOCDB.C");
  AddTaskConfigOCDB("raw://");
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  // to run a task in AliPhysics
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  //     gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/macros/AddTaskFilteredTree.C");
  gROOT->LoadMacro("./AddTaskFilteredTreeLocal.C");
  AliAnalysisTask *ana = AddTaskFilteredTree("",isMC);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");   // to set the number of events --> mgr->StartAnalysis("grid",nEvents);

}

// ----------------------------------------------------------------------------------------------------------------------
AliAnalysisGrid* CreateAlienHandler(Int_t valgrindOption = 0,TString mode="test",Int_t localOrGrid=0,TString list = "$RUN_ON_GRID_DIR/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list",TString fname="testName",Int_t isMC=0){

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
  plugin->SetAliPhysicsVersion(aliPhysicsTag); // change to something up-to-date vAN-20170717-1
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
    plugin->SetGridDataDir(Form("/alice/data/%d/%s/",year,period.Data())); // /alice/data/2015/LHC15o/000246858/pass1/15000246858039.402/AliESDs.root
    plugin->SetDataPattern(Form("/%s/*/AliESDs.root",pass.Data()));
  }
  else if (isMC==1) {  // RUN2 full MC gen+rec
    plugin->SetGridDataDir(Form("/alice/sim/%d/%s/",year,period.Data()));
    plugin->SetDataPattern("/*/AliESDs.root");
  }
  else if (isMC==3) {  // RUN1 full MC gen+rec
    plugin->SetGridDataDir(Form("/alice/sim/%s/",period.Data()));
    plugin->SetDataPattern("/*/AliESDs.root");
  }

  TString extraLibs;
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  plugin->SetGridOutputDir(yearStr); // In this case will be $HOME/work/output

  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  //   plugin->SetDefaultOutputs(0);
  //   plugin->SetOutputFiles("FilterEvents_Trees.root");
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  //plugin->SetMaxMergeFiles(40);
  plugin->SetMaxMergeStages(3);

  plugin->SetTTL(86399);
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
  if (localOrGrid==0) plugin->SetKeepLogs(kTRUE); // keep the log files

  return plugin;
}
