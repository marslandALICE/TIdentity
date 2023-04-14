R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/train/AddESDHandler.C>
#include <ANALYSIS/macros/train/AddAODHandler.C>
#include <ANALYSIS/macros/train/AddMCHandler.C>
#include <ANALYSIS/macros/AddTaskPIDResponse.C>

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <OADB/macros/AddTaskPhysicsSelection.C>
#include <OADB/macros/AddTaskCentrality.C>
#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>
#include <PWGPP/TPC/macros/AddTaskConfigOCDB.C>
#include "PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C"
//include "PWGJE/EMCALJetTasks/macros/AddTask_siweyhmi_JetHadro.C" //CHANGE use this when built into AliPhysics

R__ADD_INCLUDE_PATH($PWD)
#include "AddTask_siweyhmi_JetHadro.C" //CHANGE use this when not built into AliPhysics
#include <AddTaskFilteredTreeLocal.C>

#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"

AliAnalysisGrid* CreateAlienHandler(Int_t, TString, Int_t, TString, TString, Int_t);
class  AliAnalysisManager;
class  AliAnalysisAlien;
class  AliAnalysisTaskRho;

/*

Example usage:

cd /home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/test/root6_based/4thMoment_29092021
aliroot -b -q 'runGrid_JetHadro.C(0,0,"test",0,"3","$RUN_ON_GRID_DIR/Ebye/lists/runsONERUN-2020-LHC20e3a-pass3.list","PWGPP695_MC_remapping",1,65,2018,"18q",3,"vAN-20210925_ROOT6-1")'
aliroot -b -q 'runGrid_JetHadro.C(0,0,"test",0,"3","$RUN_ON_GRID_DIR/Ebye/lists/runs-2020-LHC20e3a-pass3.list","PWGPP695_MC_remapping",1,65,2018,"18q",3,"vAN-20210925_ROOT6-1")'

valgrindOption --> 0 --> Normal, 1--> valgrind, 2--> callgrind, 3-->Massif
modes          --> "test" --> to run over a small set of files (requires alien connection but everything stored locally), "full" --> to run over everything on the grid, "terminate" --> to merge results after "full"
localOrGrid    --> 0 --> Use only one run and run locally, 1 --> run for all runs on the grid
list           --> defines the data set according to year, period and pass --> "test-2015-LHC15o-pass5_lowIR.list"
fname          --> output directory name in my home folder in alien
isMC           --> 0:Data , 1:MCfull,  2:fastMC
setType        --> setting type encoded in Config file
lhcYear        --> year

*/

Bool_t fAddTIdentityTask = kTRUE;
Bool_t fAddJetFinderTask = kTRUE;
Bool_t fAddFilteredTrees = kTRUE;
Bool_t fUseMultSelection = kTRUE;
const Int_t nTestFiles = 10;
const Int_t nChunksPerJob = 20;
TString dataBaseDir = "/Users/sierraweyhmiller/Sierra/Caines_Research_Yr2/JetHadro/data";
TString aliPhysicsTag = "vAN-20201124-1"; //  	vAN-20180828-1  vAN-20181119-1  vAN-20190105_ROOT6-1
//
// debugging options
TString fValgrind  = "/usr/bin/valgrind --leak-check=full --leak-resolution=high --num-callers=40 --error-limit=no --show-reachable=yes  --log-file=xxx.txt --suppressions=$ROOTSYS/etc/valgrind-root.supp  -v ";
TString fCallgrind = "/usr/bin/valgrind --tool=callgrind --log-file=cpu.txt   --num-callers=40 -v  --trace-children=yes ";
TString fMassif    = "/usr/bin/valgrind --tool=massif ";

Bool_t fDoAOD = kFALSE;

void runGrid_JetHadro(Bool_t fRunLocalFiles = kTRUE, Int_t valgrindOption = 0, TString mode="test",Int_t localOrGrid=0, TString passStr="1", TString list = "", TString fname="EbyeIterPID", Int_t isMC=0, Int_t setType=3, Int_t lhcYear=2015, TString periodName="15o", Int_t passIndex=2, TString physicsTagForFullTest="vAN-20201124-1")
{

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
  }
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  // Add handlers
  AliVEventHandler* handler=0x0;
  if (isMC>0){
    AliMCEventHandler* mcHandler = AddMCHandler(kFALSE);
    mcHandler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mcHandler);
    if (fDoAOD) {
      // AliAODInputHandler* aodHandler = AddAODHandler();
      AliAODInputHandler* aodHandler = new AliAODInputHandler();
      // aodHandler->SetReadFriends(kFALSE);
      aodHandler->SetNeedField();
      mgr->SetInputEventHandler(aodHandler);
    } else {
      // AliESDInputHandler* esdHandler = AddESDHandler();
      AliESDInputHandler* esdHandler = new AliESDInputHandler();
      esdHandler->SetReadFriends(kFALSE);
      esdHandler->SetNeedField();
      mgr->SetInputEventHandler(esdHandler);
    }
  } else {
    if (fDoAOD) {
      AliAODInputHandler* aodHandler = AddAODHandler();
      aodHandler->SetNeedField();
      mgr->SetInputEventHandler(aodHandler);
    } else {
      AliESDInputHandler* esdHandler = AddESDHandler();
      esdHandler->SetReadFriends(kFALSE);
      esdHandler->SetNeedField();
      mgr->SetInputEventHandler(esdHandler);
    }
  }
  // mgr->SetInputEventHandler(handler);
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  if (!(isMC==2 || isMC==4)){
    // Add Additional tasks
    AddTaskPhysicsSelection(isMC);
    AliAnalysisTaskPIDResponse *taskPID = NULL;
    // if (isMC==0) taskPID=AddTaskPIDResponse(kFALSE,kTRUE,kTRUE,1);
    if (lhcYear<2016){
      if (isMC==0) taskPID=AddTaskPIDResponse(kFALSE,kTRUE,kFALSE,passStr);
      if (isMC==1 || isMC==3 || isMC==5 ) taskPID=AddTaskPIDResponse(kTRUE,kTRUE,kTRUE,passStr);
    } else {
      // if (isMC==0) taskPID=AddTaskPIDResponse(kFALSE,kTRUE,kFALSE,passStr,kFALSE,"TPC-OADB:COMMON/PID/data/TPCPIDResponseOADB_pileupCorr.root;TPC-Maps:$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCetaMaps_pileupCorr.root" );
      if (isMC==0) taskPID=AddTaskPIDResponse(kFALSE,kTRUE,kFALSE,passStr,kFALSE);
      // if (isMC==1 || isMC==3 || isMC==5 ) taskPID=AddTaskPIDResponse(kTRUE,kTRUE,kTRUE,passStr,kFALSE,"TPC-OADB:COMMON/PID/data/TPCPIDResponseOADB_pileupCorr.root;TPC-Maps:$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCetaMaps_pileupCorr.root" );
      if (isMC==1 || isMC==3 || isMC==5 ) taskPID=AddTaskPIDResponse(kTRUE,kTRUE,kTRUE,passStr,kFALSE);
    }
    AliCentralitySelectionTask *taskCentrality=AddTaskCentrality(kTRUE, fDoAOD);
    if(fUseMultSelection){
      AliMultSelectionTask* multTask = AddTaskMultSelection();
      if (fAddFilteredTrees)
      {
        AddTaskConfigOCDB("raw://");
      }
      std::cout << "period name = " << periodName << std::endl;
      if(periodName.Contains("15o")) multTask->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
      //
    }
  }
  //
  // Filtered tree
  if (fAddFilteredTrees && isMC==0) {
    AliAnalysisTask *ana = AddTaskFilteredTreeLocal("",isMC);
  }
  //
  // My task --> has to be compiled here instead of including
  if (fAddTIdentityTask)
  {
    gROOT->LoadMacro("AliAnalysisJetHadro.cxx++g");
    AliAnalysisTask *ana = AddTask_siweyhmi_JetHadro(kTRUE, "Config_siweyhmi_JetHadro.C", setType, lhcYear, periodName, passIndex);
  }
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  TChain *chain;
  if (!fRunLocalFiles) {
    //
    // Start analysis on grid.
    mgr->StartAnalysis("grid");   // to set the number of events --> mgr->StartAnalysis("grid",nEvents);
  } else {
    //
    // real data local files
    if (isMC==0)
    {
      chain = new TChain("esdTree");
      // TString localFiles[] =
      // {
      //   "/alice/data/2017/LHC17p/000282008/pass1_FAST/17000282008038.211/AliESDs.root"
      // };
      TString localFiles[] =
      {
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.401/AliESDs.root", //5.02 TeV PbPb data //410
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.402/AliESDs.root", //5.02 TeV PbPb data //410
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.403/AliESDs.root", //5.02 TeV PbPb data //410
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.404/AliESDs.root", //5.02 TeV PbPb data //410
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.405/AliESDs.root", //5.02 TeV PbPb data //410
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.406/AliESDs.root", //5.02 TeV PbPb data //410
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.407/AliESDs.root", //5.02 TeV PbPb data //410
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.408/AliESDs.root", //5.02 TeV PbPb data //410
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.409/AliESDs.root", //5.02 TeV PbPb data //410
        "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.410/AliESDs.root" //5.02 TeV PbPb data
        //"/alice/data/2017/LHC17p/000282008/pass1_FAST/17000282008038.211/AliESDs.root", //5.02 TeV pp data
        //"/alice/data/2017/LHC17q/000282365/pass1_FAST/17000282365039.911/AliESDs.root" //5.02 TeV pp data
      };
      for (int ifile =0; ifile<nTestFiles; ifile++) chain->AddFile(dataBaseDir+localFiles[ifile]);
    }
    // Full MC local files
    if (isMC==1)
    {
      chain = new TChain("esdTree");
      TString localFiles[] =
      {
        "/alice/sim/2020/LHC20e3a/297379/001/AliESDs.root",
        "/alice/sim/2020/LHC20e3a/297379/002/AliESDs.root"
      };
      for (int ifile =0; ifile<nTestFiles; ifile++) chain->AddFile(dataBaseDir+localFiles[ifile]);
    }
    // Generated level MC local files
    if (isMC==2)
    {
      chain = new TChain("TE");
      TString localFiles[] =
      {
        // "/media/marsland/Samsung_T5/workdir/ThirdMoment_paper/data/alice/sim/2022/LHC22d1d/244917/001/galice.root" // EPOS
        "/media/marsland/Samsung_T5/workdir/ThirdMoment_paper/data/alice/sim/2022/LHC22d1a/297595/001/galice.root" // HIJING
      };
      for (int ifile =0; ifile<nTestFiles; ifile++) chain->AddFile(dataBaseDir+localFiles[ifile]);
    }
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
  plugin->SetRunMode(mode.Data()); //"test" or "full"
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
  //alienHandler->AddRunNumber(296378);
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
    // /alice/sim/2022/LHC22d1d/244917/001
    plugin->SetAdditionalLibs("pythia6 Tree Geom VMC Physics Minuit Gui Minuit2 STEERBase ESD OADB ANALYSIS ANALYSISalice CDB STEER CORRFW EMCALUtils EMCALrec VZERObase VZEROrec");
    plugin->SetAdditionalRootLibs("libVMC.so libPhysics.so libTree.so libMinuit.so libProof.so libSTEERBase.so libESD.so libAOD.so");
    plugin->SetMCLoop(kTRUE);
    plugin->SetUseMCchain();
    plugin->SetNMCjobs(1000);
    plugin->SetNMCevents(100);
    plugin->SetGridDataDir(Form("/alice/sim/%d/%s/",year,period.Data()));
    plugin->SetDataPattern("/*/galice.root");
    //       plugin->SetDataPattern("/*/root_archive.zip#galice.root");
    plugin->SetTreeName("TE");
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
    plugin->SetTreeName("TE");
  }
  else if (isMC==5) {  // RUN1 full MC gen+rec AMPT
    std::cout << " Data SOURCE = RUN1 full MC gen+rec AMPT " << std::endl;
    plugin->SetGridDataDir(Form("/alice/sim/%d/%s/",year,period.Data()));
    plugin->SetDataPattern("/*/AliESDs.root");
  } else {
    std::cout << " Unknown data source: isMC = " << isMC << std::endl;
  }

  plugin->SetAnalysisMacro(Form("JetHadro_%s_%s.C",period.Data(),pass.Data()));
  plugin->SetExecutable(Form("JetHadro_%s_%s.sh",period.Data(),pass.Data()));
  plugin->SetJDLName(Form("JetHadro_%s_%s.jdl",period.Data(),pass.Data()));
  //
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // include additional libs
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  // to run locally
  plugin->SetAnalysisSource("AliAnalysisJetHadro.cxx");
  plugin->SetAdditionalLibs("AliAnalysisJetHadro.cxx AliAnalysisJetHadro.h");
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
  // if (localOrGrid==0) plugin->SetKeepLogs(kTRUE); // keep the log files
  plugin->SetKeepLogs(kTRUE); // keep the log files

  return plugin;
}
