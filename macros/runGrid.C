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

R__ADD_INCLUDE_PATH($PWD)
#include <AddTask_marsland_TIdentityPID.C>
#include <AddTaskFilteredTreeLocal.C>

#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"

#include <string>
#include <fstream>

using namespace std;

AliAnalysisGrid* CreateAlienHandler(Int_t, TString, Int_t, TString, TString, Int_t);
class  AliAnalysisManager;
class  AliAnalysisAlien;

vector<TString> readFilelist(TString filelist);

/*

Example usage: 

cd /home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/test/root6_based/4thMoment_29092021
aliroot -b -q 'runGrid.C(0,0,"test",0,"3","$RUN_ON_GRID_DIR/Ebye/lists/runsONERUN-2020-LHC20e3a-pass3.list","PWGPP695_MC_remapping",1,65,2018,"18q",3,"vAN-20210925_ROOT6-1")'
aliroot -b -q 'runGrid.C(0,0,"test",0,"3","$RUN_ON_GRID_DIR/Ebye/lists/runs-2020-LHC20e3a-pass3.list","PWGPP695_MC_remapping",1,65,2018,"18q",3,"vAN-20210925_ROOT6-1")'

valgrindOption --> 0 --> Normal, 1--> valgrind, 2--> callgrind, 3-->Massif
modes          --> "test"      --> to run over a small set of files (requires alien connection but everything stored locally), 
                   "full"      --> to run over everything on the grid, 
                   "terminate" --> to merge results after "full"
localOrGrid    --> 0 --> Use only one run and run locally, 1 --> run for all runs on the grid
list           --> defines the data set according to year, period and pass --> "test-2015-LHC15o-pass5_lowIR.list"
fname          --> output directory name in my home folder in alien
isMC           --> 0:Data , 1:MCfull,  2:fastMC
setType        --> setting type encoded in Config file
lhcYear        --> year

*/

Bool_t fAddFilteredTrees = kFALSE;
Bool_t fUseMultSelection = kTRUE;
Int_t nTestFiles = 1;
const Int_t nChunksPerJob = 10;
// TString dataBaseDir = "/eos/user/m/marsland/data";
// TString dataBaseDir = "/media/marsland/Samsung_T5/data";
// TString dataBaseDir = "../data";
TString dataBaseDir = "/home/ceres/fokin/work/jalien/jalien-cache/LFN/";
TString aliPhysicsTag = "vAN-20201124-1"; //  	vAN-20180828-1  vAN-20181119-1  vAN-20190105_ROOT6-1
//
// debugging options
TString fValgrind  = "/usr/bin/valgrind --leak-check=full --leak-resolution=high --num-callers=40 --error-limit=no --show-reachable=yes  --log-file=xxx.txt --suppressions=$ROOTSYS/etc/valgrind-root.supp  -v ";
TString fCallgrind = "/usr/bin/valgrind --tool=callgrind --log-file=cpu.txt   --num-callers=40 -v  --trace-children=yes ";
TString fMassif    = "/usr/bin/valgrind --tool=massif ";

Bool_t fDoAOD = kFALSE;

void runGrid(Bool_t fRunLocalFiles = kTRUE, Int_t valgrindOption = 0, TString mode="test",Int_t localOrGrid=0, TString passStr="1", TString list = "", TString fname="EbyeIterPID", Int_t isMC=0, Int_t setType=3, Int_t lhcYear=2015, TString periodName="15o", Int_t passIndex=2, TString physicsTagForFullTest="vAN-20201124-1", TString localFilesList="files.txt")
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
  if (fAddFilteredTrees) {
    AliAnalysisTask *ana = AddTaskFilteredTreeLocal("",isMC);
  }
  //
  // My task --> has to be compiled here instead of including
  gROOT->LoadMacro("AliAnalysisTaskTIdentityPID.cxx++g");
  AliAnalysisTask *ana = AddTask_marsland_TIdentityPID(kFALSE,"Config_marsland_TIdentityPID.C",setType,lhcYear,periodName,passIndex);
  //
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  //
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  TChain *chain;
  if (!fRunLocalFiles) {
    // Start analysis in grid.
    mgr->StartAnalysis("grid");   // to set the number of events --> mgr->StartAnalysis("grid",nEvents);
  } else {
    // to run over files stored locally, uncomment this section,
    // and comment out the above lines related to alienHandler and StartAnalysis("grid")
    TChain *chain = new TChain("esdTree");
    TString localFiles_LHC18q[] = {
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622033.515/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622038.309/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622025.131/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622027.508/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622028.427/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622028.121/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622039.612/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622029.529/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622030.123/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622032.508/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.205/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622037.314/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622024.431/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622037.511/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622032.113/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622028.521/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622039.100/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622030.214/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622030.427/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622022.514/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622033.606/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.403/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622020.429/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622027.408/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622028.411/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622023.226/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622038.118/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622033.107/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622025.505/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622033.425/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622020.500/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.402/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622026.522/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.315/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622038.108/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622035.325/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622026.131/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622026.616/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622021.201/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622022.200/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622038.419/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622026.529/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622026.323/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622029.427/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622025.323/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622026.221/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622023.330/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622022.206/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622039.230/AliESDs.root",
      "/alice/data/2018/LHC18q/000296622/pass3/18000296622033.206/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025052.200/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025052.130/AliESDs.root"
    };
    TString localFiles_LHC15o_pass2[] = {
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540021.114/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540026.106/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540029.112/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540022.104/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540033.103/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540034.114/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540034.103/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540036.200/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540024.100/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540033.109/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540025.209/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540019.113/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540036.110/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540026.107/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540024.109/AliESDs.root",
      // "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540034.107/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540020.105/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540026.109/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540033.100/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540034.110/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540030.105/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540025.102/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540022.112/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540020.106/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540025.112/AliESDs.root",
      // "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540032.107/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540023.102/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540032.112/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540032.106/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540022.206/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540027.100/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540021.104/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540029.102/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540027.113/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540020.100/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540039.108/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540038.104/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540038.103/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540023.111/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540039.110/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540035.114/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540022.107/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540022.100/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540022.114/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540028.101/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540024.108/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540024.106/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540033.101/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540023.100/AliESDs.root",
      "../jalien-cache/LFN/alice/data/2015/LHC15o/000246540/pass2/15000246540027.106/AliESDs.root"
    };
    TString localFiles_LHC15o_pass1[] = 
    {
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540021.114/AliESDs.root",
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540025.209/AliESDs.root",
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540020.105/AliESDs.root",
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540032.112/AliESDs.root",
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540021.104/AliESDs.root",
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540038.104/AliESDs.root",
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540035.114/AliESDs.root",
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540022.100/AliESDs.root",
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540022.114/AliESDs.root",
      "/alice/data/2015/LHC15o/000246540/pass1/15000246540023.100/AliESDs.root"
    };
    TString localFiles_LHC10h[] = 
    {
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025052.200/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025052.130/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025030.90/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025062.80/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025032.30/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025055.50/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025016.60/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025028.70/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025013.40/AliESDs.root",
      "/alice/data/2010/LHC10h/000139025/pass3/10000139025028.200/AliESDs.root"
    };
    TString* localFiles = localFiles_LHC18q;

    if (localFilesList == "") {
      for (int ifile =0; ifile<nTestFiles; ifile++) chain->AddFile(dataBaseDir+localFiles_LHC10h[ifile]);
    } else {
      vector<TString> localFilesTxt = readFilelist(localFilesList);
      for (TString localFile : localFilesTxt) chain->AddFile(localFile);
    }
    chain->Print();
    mgr->StartAnalysis("local", chain);
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

vector<TString> readFilelist(TString filelist) {
    string line;
    ifstream inputfile;
    inputfile.open(filelist);

    vector<TString> ret;

    while (getline(inputfile, line)) {
        ret.push_back(line);
    }

    inputfile.close();

    return ret;
}