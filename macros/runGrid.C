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
 
    cd /home/marsland/Desktop/RUN_ON_GRID/Ebye/testMC
    rm *.xml Task* *.root stderr out* Ali*
    cp /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/*.*  .
 
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------
 
 // tests real data 
 aliroot -b -q 'runGrid.C("test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2015-LHC15o-pass5_lowIR.list","EbyeIterPID_data",0,22,2)' &> stdout &
 aliroot -b -q 'runGrid.C("full",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2015-LHC15o-pass5_lowIR.list","EbyeIterPID_data",0,22,2)' &> stdout &
 
 aliroot -b -q 'runGrid.C("full",1,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2015-LHC15o-pass5_lowIR.list","EbyeIterPID_data",0,22,2)' 
 aliroot -b -q 'runGrid.C("terminate",1,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2015-LHC15o-pass5_lowIR.list","EbyeIterPID_data",0,22,2)' 

 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------

 // tests MC full
 aliroot -b -q 'runGrid.C("test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass1.list","EbyeIterPID_MC_FULL",1,110,2)' &> stdout &
 aliroot -b -q 'runGrid.C("full",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass1.list","EbyeIterPID_MC_FULL",1,110,2)' &> stdout &
 aliroot -b -q 'runGrid.C("terminate",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass1.list","EbyeIterPID_MC_FULL",1,110,2)' 

 // full stat MC full
 aliroot -b -q 'runGrid.C("full",1,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2017-LHC17c5b-pass5_lowIR.list","EbyeIterPID_MC_FULL",1,110,2)'
 aliroot -b -q 'runGrid.C("full",1,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass3.list","EbyeIterPID_MC_FULL",1,110,2)'
 aliroot -b -q 'runGrid.C("full",1,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b2-pass3.list","EbyeIterPID_MC_FULL",1,110,2)'
 aliroot -b -q 'runGrid.C("full",1,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16g1-pass1.list","EbyeIterPID_MC_FULL",1,110,2)'

 
 
 runs-2016-LHC16k3b-pass1.list
 runs-2016-LHC16k3b2-pass1.list
 runs-2017-LHC17c5b-pass5_lowIR.list
 runs-2016-LHC16g1-pass1.list

 
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------

 // tests MC Fast
 aliroot -b -q 'runGrid.C("test",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass1.list","EbyeIterPID_MC_FAST",2,111,2)' &> stdout &
 aliroot -b -q 'runGrid.C("full",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass1.list","EbyeIterPID_MC_FAST",2,111,2)' &> stdout &
 aliroot -b -q 'runGrid.C("terminate",0,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass1.list","EbyeIterPID_MC_FAST",2,111,2)' 

 // full stat MC full
 aliroot -b -q 'runGrid.C("full",1,"/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/runs-2016-LHC16k3b-pass1.list","EbyeIterPID_MC_FAST",2,111,2)'

 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------
 // ----------------------------------------------------------------------------------------------------------------------------------------

 
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

void runGrid(TString mode="test",Int_t localOrGrid=0, TString list = "/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list", TString fname="EbyeIterPID", Int_t isMC=0, Int_t setType=109, Int_t lhcPeriod=1) 
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
  AliAnalysisGrid *alienHandler = CreateAlienHandler(mode,localOrGrid,list,fname,isMC);  
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
  if (isMC==2 || isMC==1){
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

AliAnalysisGrid* CreateAlienHandler(TString mode="test",Int_t localOrGrid=0,TString list = "/home/marsland/Desktop/RUN_ON_GRID/Ebye/lists/test-2015-LHC15o-pass5_lowIR.list",TString fname="testName",Bool_t isMC=kFALSE){
    
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(2);
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20170507-1"); // change to something up-to-date
  plugin->SetNrunsPerMaster(1);
  plugin->SetSplitMaxInputFileNumber(130); // 3 in the LEGO trains
  
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
      if ( localOrGrid=0 && nRuns==1 ) {
          cout << " in test mode process only one run " << endl; 
          break;
    }
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
  cout << " --------------------------------------- " << endl;
  cout << " year = " << year << "   period = " << period << "   pass = " << pass << endl;  
  cout << " --------------------------------------- " << endl;
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // -----------------------------------------------------------------------------------------
  // Set filenames, input and output directories on alien
  plugin->SetGridWorkingDir(Form("%s/%s_%s/",fname.Data(),period.Data(),pass.Data()));
  if (isMC==0) {       // data
      plugin->SetGridDataDir(Form("/alice/data/%d/%s/",year,period.Data()));
      plugin->SetDataPattern(Form("/%s/*/AliESDs.root",pass.Data()));
  }
  else if (isMC==1) {  // full MC gen+rec
      plugin->SetGridDataDir(Form("/alice/sim/%d/%s/",year,period.Data()));
      plugin->SetDataPattern("/*/AliESDs.root");
  }
  else if (isMC==2) {  // fast MC gen
      plugin->SetGridDataDir(Form("/alice/sim/%d/%s/",year,period.Data()));
      plugin->SetDataPattern("/*/galice.root");
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
  plugin->SetMaxMergeStages(4);
  
  plugin->SetTTL(86399);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  //plugin->SetSplitMaxInputFileNumber();
  plugin->SetSplitMode("se");
  plugin->SetKeepLogs(kFALSE); 
  plugin->SetOutputToRunNo(1);
  // my settings
  if (localOrGrid==0) plugin->SetKeepLogs(kTRUE); // keep the log files

  return plugin;
}
