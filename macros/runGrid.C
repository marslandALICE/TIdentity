class  AliAnalysisManager;
class  AliAnalysisAlien;

//
// modes: "test" to run over a small set of files (requires alien connection but everything stored locally),
//        "full" to run over everything on the grid,
//        "terminate" to merge results after "full"
/*
 
 How to run:
 
 cd /home/marsland/Desktop/RUN_ON_GRID/Ebye/work_dir
 rm *.xml Task* *.root stderr out* Ali*
 cp /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/*.* .
 
 // tests
 aliroot -b -q 'runGrid.C("test","/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list","EbyeIterPID",0,22,2)' &> stdout &
 aliroot -b -q 'runGrid.C("full","/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list","EbyeIterPID",0,22,2)' &> stdout &
 aliroot -b -q 'runGrid.C("terminate","/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list","EbyeIterPID",0,22,2)' 

 // full stat submission
 aliroot -b -q 'runGrid.C("full","/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2015-LHC15o-pass5_lowIR.list","EbyeIterPID",0,22,2)'
 aliroot -b -q 'runGrid.C("terminate","/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2015-LHC15o-pass5_lowIR.list","EbyeIterPID",0,22,2)'

 */

void runGrid(TString mode="test",TString list = "/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list", TString fname="EbyeIterPID", Int_t isMC=0, Int_t setType=109, Int_t lhcPeriod=1) 
{
    //
    // mode      --> "test" or "full" or "terminate"
    // list      --> defines the data set according to year, period and pass --> "test-2015-LHC15o-pass5_lowIR.list"
    // fname     --> output directory name in my home folder in alien
    // isMC      --> 0:Data , 1:MCfull,  2:fastMC
    // setType   --> setting type encoded in Config file
    // lhcPeriod --> 1 for RUN1, and 2 for RUN2
    //
    
  AliLog::SetGlobalDebugLevel(5);

  //__________________________________________________________________________
  // Use AliRoot includes to compile our task
  //gROOT->ProcessLine(".include $ALICE_PHYSICS/include $ALICE_ROOT/include");

  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  // Create and configure the alien handler plugin
  AliAnalysisGrid *alienHandler = CreateAlienHandler(mode,list,fname,isMC);  
  if (!alienHandler) return;
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  //   AliAODInputHandler* aodH = new AliAODInputHandler();
  AliESDInputHandler* aodH = new AliESDInputHandler();
  mgr->SetInputEventHandler(aodH);

  // Add PIDResponse task
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(isMC);
  
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(kTRUE,kTRUE,kTRUE);
    
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality=AddTaskCentrality();
  
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AddTaskMultSelection();
  
  // ----------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------
  // Add handlers
  if (isMC==2){
      gROOT->LoadMacro(gSystem->ExpandPathName("$AliRoot_SRC/ANALYSIS/macros/train/AddMCHandler.C"));
      AliVEventHandler* handler = AddMCHandler(kFALSE);  
  }
  //   gROOT->LoadMacro(gSystem->ExpandPathName("$AliRoot_SRC/ANALYSIS/macros/train/AddESDHandler.C"));
  //   AliVEventHandler* handler = AddESDHandler();
  //   ((AliESDInputHandler*)handler)->SetReadFriends(kFALSE);
  //   ((AliESDInputHandler*)handler)->SetNeedField();
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
  /*TChain *chain = new TChain("aodTree");
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
  mgr->StartAnalysis("local",chain);*/
  
}

AliAnalysisGrid* CreateAlienHandler(TString mode="test",TString list = "/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list",TString fname="testName",Bool_t isMC=kFALSE){
    
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(2);
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20170507-1"); // change to something up-to-date
  plugin->SetNrunsPerMaster(1);
  plugin->SetSplitMaxInputFileNumber(30); // 3 in the LEGO trains
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
      //       if ( mode.Contains("test") && nRuns==1 ) {cout << " in test mode process only one run " << endl; break;}
  }
  gSystem->Exec("rm tmp.list");
  TObjArray *objArr1 = list.Tokenize("/");
  TString fileName = ((objArr1->At((Int_t)objArr1->GetLast()))->GetName());
  TObjArray *objArr2  = fileName.Tokenize("-");
  Int_t year         = atoi((objArr2->At(1))->GetName());
  TString period     = ((objArr2->At(2))->GetName());
  TString passLong   = ((objArr2->At(3))->GetName());
  TObjArray *objArr3  = passLong.Tokenize(".");
  TString pass = ((objArr3->At(0))->GetName());
  TDatime tt;
  Int_t date   = tt.GetDate();
  Int_t hour   = tt.GetHour();
  Int_t minute = tt.GetMinute();
  cout << " --------------------------------------- " << endl;
  cout << " year = " << year << "   period = " << period << "   pass = " << pass << endl;  
  cout << " --------------------------------------- " << endl;
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // Set filenames, input and output directories on alien
  //   plugin->SetGridWorkingDir(Form("%s/%s_%s_%d_%d_%d/",fname.Data(),period.Data(),pass.Data(),date,hour,minute));
  plugin->SetGridWorkingDir(Form("%s/%s_%s/",fname.Data(),period.Data(),pass.Data()));
  // set input path "isMC = Data:0 MCfull:1  MCfast:2" 
  if (isMC==0) {
      plugin->SetDataPattern(Form("/%s/*/AliESDs.root",pass.Data()));
      plugin->SetGridDataDir(Form("/alice/data/%d/%s/",year,period.Data()));
  }
  else if (isMC==1) {
      plugin->SetDataPattern("/*/AliESDs.root");
      plugin->SetGridDataDir(Form("/alice/sim/%s/",period.Data()));
  }
  else if (isMC==2) {
      plugin->SetDataPattern("/*/galice.root");
      plugin->SetGridDataDir(Form("/alice/sim/%s/",period.Data()));
  }
  plugin->SetAnalysisMacro(Form("TaskEbyeIterPIDMC_%s_%s.C",period.Data(),pass.Data()));
  plugin->SetExecutable(Form("TaskEbyeIterPIDMC_%s_%s.sh",period.Data(),pass.Data()));
  plugin->SetJDLName(Form("TaskEbyeIterPIDMC_%s_%s.jdl",period.Data(),pass.Data()));
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

  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  //plugin->SetDefaultOutputs(0);
  //plugin->SetOutputFiles("AnalysisResults.root.root");
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  //plugin->SetMaxMergeFiles(40);
  plugin->SetMaxMergeStages(1);
  
  plugin->SetTTL(86399);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  //plugin->SetSplitMaxInputFileNumber();
  plugin->SetSplitMode("se");
  
  // my settings
  plugin->SetOutputToRunNo(1);
  plugin->SetKeepLogs(kTRUE); // keep the log files

  return plugin;
}
