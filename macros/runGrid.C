class  AliAnalysisManager;
class  AliAnalysisAlien;
//
// modes: "test" to run over a small set of files (requires alien connection but everything stored locally),
//        "full" to run over everything on the grid,
//        "terminate" to merge results after "full"
/*
 
 Estimate the number of input files 
 
    source /u/marsland/PHD/macros/marsland_EbyeRatios/alihadd_GRIDoutput.sh; 
    CountFiles "/alice/sim/2017/LHC17c5b" "AliESDs.root"   --> 252442 
    CountFiles "/alice/sim/2016/LHC16k3b" "AliESDs.root"   --> 90082
    CountFiles "/alice/sim/2016/LHC16k3b2" "AliESDs.root"  --> 92380
    CountFiles "/alice/sim/2016/LHC16g1" "AliESDs.root"    --> 321832

 How to run:
 
    cd /home/marsland/Desktop/RUN_ON_GRID/Ebye/massif
    rm *.xml Task* *.root stderr out* Ali*
    cp /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/*.*  .

 // run Massif: 0 --> Normal, 1--> valgrind, 2--> callgrind, 3-->Massif 
 aliroot -b -q 'runGrid.C(3,"test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2015-LHC15o-pass5_lowIR.list","TPC_dEdx_Info",0,21,2)' &> stdout &

 // ----------------------------------------------------------------------------------------------------------------------------------------
 // tests real data 
 aliroot -b -q 'runGrid.C(0,"test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runForMCcomp-2015-LHC15o-pass1.list","TPC_dEdx_Info",0,21,2)' 
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // tests real data 
 aliroot -b -q 'runGrid.C(0,"test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2015-LHC15o-pass5_lowIR.list","TPC_dEdx_Info",0,21,2)' &> stdout &
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // tests RUN2 MC full
 aliroot -b -q 'runGrid.C(0,"test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass3.list","TPC_dEdx_Info",1,110,2)' 
 aliroot -b -q 'runGrid.C(0,"test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b2-pass3.list","TPC_dEdx_Info",1,110,2)' 
 aliroot -b -q 'runGrid.C(0,"test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2017-LHC17c5b-pass5_lowIR.list","TPC_dEdx_Info",1,110,2)' 
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // tests RUN2 MC Fast
 aliroot -b -q 'runGrid.C(0,"test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass3.list","TPC_dEdx_Info",2,111,2)' 
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // tests RUN1 MC full
 aliroot -b -q 'runGrid.C(0,"test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2010-LHC11a10a_bis-pass1.list","TPC_dEdx_Info",3,110,1)' 
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // tests RUN1 MC fast
 aliroot -b -q 'runGrid.C(0,"test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2013-LHC13f3b-pass1.list","TPC_dEdx_Info",4,111,1)' 

 
 
 LHC16k3b     -->  Pb-Pb, 5.02 TeV - HIJING with pileup anchored to 15o (pass3)   (one run -->  246751)
 LHC16k3b2    -->  Pb-Pb, 5.02 TeV - HIJING no pile-up (ion tail/xtalk) anchored to 15o (pass3)  (one run -->  246751)
 LHC17c5b     -->  Pb-Pb, 5.02 TeV - HIJING anchored to Low IR runs LHC15o (pass5)  (same list as --> runs-2015-LHC15o-pass5_lowIR.list)
 
 LHC16g1      -->  Pb-Pb, 5.02 TeV - HIJING min bias - anchored to LHC15o (pass1 highIR)
 
 LHC16k3a2    -->  p-p, 5.02 TeV - Pythia6 no pile-up (ion tail/xtalk) anchored to 15n (pass2)
 LHC16h8a     -->  p-p, 5.02 TeV - Pythia8_Monash2013 anchored to LHC15n
 LHC16h8b     -->  p-p, 5.02 TeV - Pythia6_Perugia2011 anchored to LHC15n

 Problems faced:
 saving --> delete the output dir on alien /alice/cern.ch/user/m/marsland/EbyeIterPID/LHC15o_pass5_lowIR
 
 */

const Int_t nTestFiles = 2;
const Int_t nChunksPerJob = 70;

TString fValgrind  = "/usr/bin/valgrind --leak-check=full --leak-resolution=high --num-callers=40 --error-limit=no --show-reachable=yes  --log-file=xxx.txt --suppressions=$ROOTSYS/etc/valgrind-root.supp  -v ";
TString fCallgrind = "/usr/bin/valgrind --tool=callgrind --log-file=cpu.txt   --num-callers=40 -v  --trace-children=yes ";
TString fMassif    = "/usr/bin/valgrind --tool=massif ";
 
void runGrid(Int_t valgrindOption = 0, TString mode="test",Int_t localOrGrid=0, TString list = "/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list", TString fname="EbyeIterPID", Int_t isMC=0, Int_t setType=109, Int_t lhcPeriod=1) 
{
    //
    // valgrindOption --> 0 --> Normal, 1--> valgrind, 2--> callgrind, 3-->Massif 
    // mode           --> "test" or "full" or "terminate"
    // list           --> defines the data set according to year, period and pass --> "test-2015-LHC15o-pass5_lowIR.list"
    // fname          --> output directory name in my home folder in alien
    // isMC           --> 0:Data , 1:MCfull,  2:fastMC
    // setType        --> setting type encoded in Config file
    // lhcPeriod      --> 1 for RUN1, and 2 for RUN2
    //
    
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
  //   AliAODInputHandler* aodH = new AliAODInputHandler();
  //   AliESDInputHandler* aodH = new AliESDInputHandler();
  //   mgr->SetInputEventHandler(aodH);
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
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
  } else {
      gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C"));
      handler = AddESDHandler();
  }
  mgr->SetInputEventHandler(handler);
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  if (!(isMC==2 || isMC==4)){
      // Add Additional tasks
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      AddTaskPhysicsSelection(isMC);
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
      AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(kTRUE,kTRUE,kTRUE);
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
      AliCentralitySelectionTask *taskCentrality=AddTaskCentrality();
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
      AddTaskMultSelection();
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
  gROOT->LoadMacro("/home/marsland/Desktop/RUN_ON_GRID/Ebye/code/AliAnalysisTaskTIdentityPID.cxx+g"); 
  gROOT->LoadMacro("/home/marsland/Desktop/RUN_ON_GRID/Ebye/code/AddTask_marsland_TIdentityPID.C");  
  AliAnalysisTask *ana = AddTask_marsland_TIdentityPID(kFALSE,"Config_marsland_TIdentityPID.C",setType,lhcPeriod);
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------

  if (!mgr->InitAnalysis())
    return;
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");   // to set the number of events --> mgr->StartAnalysis("grid",nEvents);
  
  // to run over files stored locally, uncomment this section,
  // and comment out the above lines related to alienHandler and StartAnalysis("grid")
  /*
    TChain *chain = new TChain("aodTree");
    chain->AddFile("./traintest/AliAOD_0100.root");
    chain->AddFile("./traintest/AliAOD_0101.root");
    chain->AddFile("./traintest/AliAOD_0102.root");
    chain->AddFile("./traintest/AliAOD_0103.root");
    chain->AddFile("./traintest/AliAOD_0104.root");
    chain->AddFile("./traintest/AliAOD_0105.root");
    chain->AddFile("./traintest/AliAOD_0106.root");
    chain->AddFile("./traintest/AliAOD_0107.root");
    chain->AddFile("./traintest/AliAOD_0108.root");
    chain->AddFile("./traintest/AliAOD_0109.root");
    chain->Print();
    mgr->StartAnalysis("local",chain);
  */
  
}

AliAnalysisGrid* CreateAlienHandler(Int_t valgrindOption = 0,TString mode="test",Int_t localOrGrid=0,TString list = "/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list",TString fname="testName",Int_t isMC=0){
    
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
  plugin->SetAliPhysicsVersion("vAN-20170507-1"); // change to something up-to-date vAN-20170717-1
  //   plugin->SetAliPhysicsVersion("vAN-20170717-1"); // change to something up-to-date vAN-20170717-1   vAN-20170627-1
  plugin->SetNrunsPerMaster(1);
  plugin->SetSplitMaxInputFileNumber(nChunksPerJob); // 3 in the LEGO trains
  
  if (!isMC) plugin->SetRunPrefix("000");  // fort the data
  // -----------------------------------------------------------------------------------------
  // ------------------------- Read runs from the list----------------------------------------
  // -----------------------------------------------------------------------------------------
  // Read runs from the list
  gSystem->Exec(Form("grep -v run  %s | awk '{print $1}' > tmp.list",list.Data()));
  gSystem->Exec(Form("more  %s",list.Data()));
  gSystem->Exec("more tmp.list");
  FILE *fileTmp=fopen("tmp.list","r");
  Int_t nRuns=0, run=0, runPrev=0;
  cout << " ----------- list of good runs ----------- " << endl;
  while (!feof(fileTmp)) {
      fscanf(fileTmp,"%d",&run);
      if (runPrev != run) runPrev=run;
      else continue;
      plugin->AddRunNumber(run);
      cout << run << endl;
      nRuns++;
      if ( localOrGrid==0 && nRuns==1 ) {
          cout << " in test mode process only one run " << endl; 
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
  cout << " --------------------------------------- " << endl;
  cout << " year = " << year << "   period = " << period << "   pass = " << pass << endl;  
  cout << " --------------------------------------- " << endl;
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // Set filenames, input and output directories on alien
  plugin->SetGridWorkingDir(Form("%s/%s_%s/",fname.Data(),period.Data(),pass.Data()));
  if (isMC==0) {       // data
      plugin->SetGridDataDir(Form("/alice/data/%d/%s/",year,period.Data())); // /alice/data/2015/LHC15o/000246858/pass1/15000246858039.402/AliESDs.root
      plugin->SetDataPattern(Form("/%s/*/AliESDs.root",pass.Data()));
  }
  else if (isMC==1) {  // RUN2 full MC gen+rec
      plugin->SetGridDataDir(Form("/alice/sim/%d/%s/",year,period.Data()));
      plugin->SetDataPattern("/*/AliESDs.root");
  }
  else if (isMC==2) {  // RUN2 fast MC gen
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
  else if (isMC==3) {  // RUN1 full MC gen+rec
      plugin->SetGridDataDir(Form("/alice/sim/%s/",period.Data()));
      plugin->SetDataPattern("/*/AliESDs.root");
  }
  else if (isMC==4) {  // RUN1 fast MC gen
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
  //plugin->SetOutputFiles("AnalysisResults.root.root");
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
  if (localOrGrid==0) plugin->SetKeepLogs(kTRUE); // keep the log files

  return plugin;
}
