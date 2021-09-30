AliESDtrackCuts* CreateCuts(Bool_t=kTRUE, Bool_t=kTRUE);
#include "AliAnalysisTaskFilteredTree.h"
#include "TRandom.h"
#include "AliFilteredTreeEventCuts.h"
#include "AliFilteredTreeAcceptanceCuts.h"
#include "TSystem.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"


AliAnalysisTask* AddTaskFilteredTreeLocal(TString outputFile="", Int_t isMC=0)
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGUDbase");
  gSystem->Load("libTPCcalib");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWGLFspectra");


  gRandom->SetSeed(0);

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTaskFilteredTree", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kDebug);
  mgr->SetDebugLevel(1);



  //
  // Create physics trigger selection class
  //
  //AliPhysicsSelection *physTrigSel =  new AliPhysicsSelection();

  //
  // Create event cuts
  //
  Float_t zvWindow = 30. ;

  AliFilteredTreeEventCuts *evtCuts = new AliFilteredTreeEventCuts("AliFilteredTreeEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kFALSE);
  //evtCuts->SetTriggerRequired(kTRUE);

  //
  // Create geom. acceptance cuts
  //
  Float_t etaWindow = 1.0 ;
  Float_t ptMin = 0.15 ;

  AliFilteredTreeAcceptanceCuts *accCuts = new AliFilteredTreeAcceptanceCuts("AliFilteredTreeAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);
  accCuts->SetMaxDCAr(3.0);
  accCuts->SetMaxDCAz(30.0);

  //
  // Create standard esd track cuts
  //

  AliESDtrackCuts* esdTrackCuts = CreateCuts();
  if (!esdTrackCuts) {
    printf("ERROR: esdTrackCuts could not be created\n");
    return 0;
  } else {
    esdTrackCuts->SetHistogramsOn(kTRUE);
  }

  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //
  // Create task
  //
  AliAnalysisTaskFilteredTree *task = new AliAnalysisTaskFilteredTree("AliAnalysisTaskFilteredTree");
  //task->SetUseMCInfo(hasMC);
  // Double_t dScaleFactorPt = (isMC==0) ? 5000 : 100;
  // Double_t dScaleFaorV0 = (isMC==0) ? 250 : 10;
  Double_t dScaleFactorPt = (isMC==0) ? 20000 : 100;
  Double_t dScaleFactorV0 = (isMC==0) ? 1000 : 10;
  std::cout <<  " Pt downscaling = " << dScaleFactorPt << std::endl;
  std::cout <<  " V0 downscaling = " << dScaleFactorV0 << std::endl;
  task->SetLowPtTrackDownscaligF(dScaleFactorPt);
  task->SetLowPtV0DownscaligF(dScaleFactorV0);
  //task->SetLowPtTrackDownscaligF(1.e5);   // default
  //task->SetLowPtV0DownscaligF(2.e3);      // default
  task->SetProcessAll(kTRUE);
  task->SetProcessCosmics(kTRUE);
  task->SetUseESDfriends(kTRUE);
  //task->SetProcessAll(kFALSE);
  //task->SetFillTrees(kFALSE); // only histograms are filled

  // trigger
  //task->SelectCollisionCandidates(AliVEvent::kMB);
  //task->SelectCollisionCandidates(AliVEvent::kINT7);

  //
  // set analysis options from the Helper here !!!
  //
  // AlidNdPtHelper::OutputObject outputObject = AlidNdPtHelper::kCutAnalysisPbPb;
  // AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kAllPart ;
  //AlidNdPtHelper::AnalysisMode analysisMode = AlidNdPtHelper::kTPCITS;

  task->SetEventCuts(evtCuts);
  task->SetAcceptanceCuts(accCuts);
  task->SetTrackCuts(esdTrackCuts);
  task->SetAnalysisMode(AliAnalysisTaskFilteredTree::kTPCITSAnalysisMode);
  task->SetCentralityEstimator("V0M");

  // Add task
  mgr->AddTask(task);

  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  if (outputFile.IsNull())
    outputFile=Form("%s", AliAnalysisManager::GetCommonFileName());

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("filtered1", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());
  mgr->ConnectOutput(task, 1, coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("filtered2", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());
  mgr->ConnectOutput(task, 2, coutput2);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("filtered3", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());
  mgr->ConnectOutput(task, 3, coutput3);
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("filtered4", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());
  mgr->ConnectOutput(task, 4, coutput4);
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("filtered5", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());
  mgr->ConnectOutput(task, 5, coutput5);
  AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("filtered6", TTree::Class(), AliAnalysisManager::kOutputContainer, outputFile.Data());
  mgr->ConnectOutput(task, 6, coutput6);


 // store histograms in the separate file
  TString outputFileHisto = "PtResHistograms.root";

  AliAnalysisDataContainer *coutput7 = mgr->CreateContainer("histo7", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileHisto.Data());
  mgr->ConnectOutput(task, 7, coutput7);


  return task;
}

AliESDtrackCuts* CreateCuts(Bool_t fieldOn = kTRUE, Bool_t hists = kTRUE)
{
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
  if(!esdTrackCuts) return 0;

  if (hists)
    esdTrackCuts->DefineHistograms(1);

  Double_t cov1, cov2, cov3, cov4, cov5;
  Double_t nSigma;
  Double_t maxDCAtoVertex, maxDCAtoVertexXY, maxDCAtoVertexZ;
  Double_t minNClustersTPC;
  Double_t maxChi2PerClusterTPC;
  Double_t minPt, maxPt;

  //
  // TPC
  //
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  //
  // ITS
  //
  //esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  //

  TString tag = "TPC+ITS refit and KinkRejection required - for cut studies";

  // cuts for data without field
  if (!fieldOn)
  {
    cov5 = 1e10;
    tag += " without field";
  }

  Printf("Created track cuts for: %s", tag.Data());

  return esdTrackCuts;
}
