//
// Configuration file for the AliAnalysisJetHadro.cxx class
//
void SetDefaultsJets(AliAnalysisJetHadro *defaultTask, Int_t year, TString periodName, Int_t passIndex);
//
//
AliAnalysisJetHadro* Config_siweyhmi_JetHadro(Bool_t getFromAlien, Int_t settingType, Int_t year, TString periodName, Int_t passIndex, TString combinedName) {
  //
  //
  std::cout << " Info::siweyhmi: ===== In the Config --> Running with year = ";
  std::cout << year << " --- period name = " << periodName << " --- pass = " << passIndex << " --- settingType = " << settingType << std::endl;
  //
  // Configuration for the jet task
  //   AliEmcalJetTask* AddTaskEmcalJet(
  //   const char *nTracks                        = "usedefault",
  //   const char *nClusters                      = "usedefault",
  //   const AliJetContainer::EJetAlgo_t jetAlgo  = AliJetContainer::antikt_algorithm,
  //   const Double_t radius                      = 0.4,
  //   const AliJetContainer::EJetType_t jetType  = AliJetContainer::kFullJet,
  //   const Double_t minTrPt                     = 0.15,
  //   const Double_t minClPt                     = 0.30,
  //   const Double_t ghostArea                   = 0.005,
  //   const AliJetContainer::ERecoScheme_t reco  = AliJetContainer::pt_scheme,
  //   const char *tag                            = "Jet",
  //   const Double_t minJetPt                    = 0.,
  //   const Bool_t lockTask                      = kTRUE,
  //   const Bool_t bFillGhosts                   = kFALSE,
  //   const char *suffix                         = ""
  // )
  //
  // Create specific track selection
  UInt_t kPhysSel = AliVEvent::kCentral;
  if (year==2010) kPhysSel = AliVEvent::kMB;   // select minimum bias events for LHC10h
  if (year==2015) kPhysSel = AliVEvent::kINT7; // select minimum bias events for LHC15o
  // if (year==2018) task->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral) ;

  //Rho task
  TString sRhoName = "Rho";
  float emcjetradius = 0.4;
  float emcbgjetradius = 0.2;
  float emcGhostArea = 0.005;
  float fTrackPt = 0.15;
  float pT_sub_min = 40.0;    // can be 40, 60, 80
  double particleEtaCut = 0.9;
  //
  // track cuts for the container
  AliESDtrackCuts *fESDtrackCuts = new AliESDtrackCuts;
  fESDtrackCuts->SetEtaRange(-1*particleEtaCut,particleEtaCut);
  fESDtrackCuts->SetPtRange(fTrackPt,1000.);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts->SetMaxChi2PerClusterITS(36);
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);    // ?? FROM MARIAN
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetMinNCrossedRowsTPC(70);
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0208+0.04/pt^1.01");
  fESDtrackCuts->SetMaxDCAToVertexXY(2.4);   // hybrid cuts  TODO
  fESDtrackCuts->SetMaxDCAToVertexZ(3.2);    // hybrid cuts  TODO
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts->SetDCAToVertex2D(kTRUE);  // fESDtrackCuts->SetDCAToVertex2D(kFALSE);    TODO
  fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  if ( (year==2015 && passIndex==2) || (year==2018 && passIndex==3) ){
    fESDtrackCuts->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
  }

  //Background jet task
  AliEmcalJetTask *pKtChJetTask = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "usedefault", AliJetContainer::kt_algorithm, emcbgjetradius, AliJetContainer::kChargedJet, fTrackPt, 0, emcGhostArea, AliJetContainer::E_scheme, "Jet", 0., kFALSE, kFALSE);
  if (year<2018) pKtChJetTask->SelectCollisionCandidates(kPhysSel);
  pKtChJetTask->SetMinJetPt(fTrackPt);
  pKtChJetTask->SetUseNewCentralityEstimation(kTRUE);
  if (year == 2017) {pKtChJetTask->SetForceBeamType(AliAnalysisTaskEmcal::kpp);}
  if (year == 2018) {pKtChJetTask->SetForceBeamType(AliAnalysisTaskEmcal::kAA);}
  // pKtChJetTask->SetNCentBins(5);

  //Rho task
  AliAnalysisTaskRho *RhoTask = AliAnalysisTaskRho::AddTaskRhoNew("usedefault", "", sRhoName, emcbgjetradius, AliEmcalJet::kTPCfid, AliJetContainer::kChargedJet, kTRUE, AliJetContainer::E_scheme);
  if (year<2018) RhoTask->SelectCollisionCandidates(kPhysSel);
  RhoTask->SetExcludeLeadJets(2);
  RhoTask->SetUseNewCentralityEstimation(kTRUE);
  if (year == 2017) {RhoTask->SetForceBeamType(AliAnalysisTaskEmcal::kpp);}
  if (year == 2018) {RhoTask->SetForceBeamType(AliAnalysisTaskEmcal::kAA);}
  // RhoTask->SetNCentBins(5);

  AliTrackContainer * trackContbg = pKtChJetTask->GetTrackContainer("Tracks");
  auto trackcontbg = pKtChJetTask->AddTrackContainer("Tracks");
  trackcontbg->SetTrackFilterType(AliEmcalTrackSelection::ETrackFilterType_t::kCustomTrackFilter);
  trackcontbg->AddTrackCuts(new PWG::EMCAL::AliEmcalESDtrackCutsWrapper("ESDcuts", fESDtrackCuts));
  trackcontbg->SetName("bgdetTracks");

  RhoTask->RemoveJetContainer(0);
  RhoTask->AddJetContainer("Jet_KTChargedR020_Tracks_pT0150_E_scheme");
  AliJetContainer* bgjetCont = RhoTask->GetJetContainer("Jet_KTChargedR020_Tracks_pT0150_E_scheme");
  bgjetCont->SetJetPtCut(fTrackPt);
  bgjetCont->SetJetEtaLimits(-0.9+emcbgjetradius, 0.9-emcbgjetradius);
  bgjetCont->PrintCuts();

  // Signal jet task
  AliEmcalJetTask *pChJet02Task = AliEmcalJetTask::AddTaskEmcalJet("Tracks", "", AliJetContainer::antikt_algorithm, emcjetradius, AliJetContainer::kChargedJet, fTrackPt, 0, emcGhostArea, AliJetContainer::E_scheme, "SigJet", fTrackPt, kFALSE, kFALSE);
  if (year<2018) pChJet02Task->SelectCollisionCandidates(kPhysSel);

  // Create the track container and apply track cuts
  AliTrackContainer * trackCont = pChJet02Task->GetTrackContainer("Tracks");
  auto trackcont = pChJet02Task->AddTrackContainer("Tracks");
  trackcont->SetTrackFilterType(AliEmcalTrackSelection::ETrackFilterType_t::kCustomTrackFilter);
  trackcont->AddTrackCuts(new PWG::EMCAL::AliEmcalESDtrackCutsWrapper("ESDcuts", fESDtrackCuts));
  trackcont->SetName("detTracks");
  //
  // Configuration for the main analysis task
  AliAnalysisJetHadro *task = new AliAnalysisJetHadro(combinedName);
  SetDefaultsJets(task,year,periodName,passIndex);
  if (year<2018) task->SelectCollisionCandidates(kPhysSel);
  AliJetContainer* jetCont02 = task->AddJetContainer("SigJet_AKTChargedR040_Tracks_pT0150_E_scheme"); // has to modified accordingly with AliEmcalJetTask
  jetCont02->SetJetPtCut(fTrackPt);
  jetCont02->SetJetEtaLimits(-0.9+emcjetradius, 0.9-emcjetradius);

  jetCont02->SetJetRadius(emcjetradius);
  jetCont02->SetName("detJets");
  if (year != 2017) {
    jetCont02->SetRhoName(sRhoName);
    TString fRhoName = jetCont02->GetRhoName();
    std::cout << " Rho Name is " << fRhoName << endl;
  }

  jetCont02->PrintCuts();
  //
  // Other Specific settings

  switch (settingType) {
    //
    // ====================================================================================
    // ============================= Real Data Settings ===================================
    // ====================================================================================
    //
    case 0:{
      std::cout << " SETTING TYPE = " << settingType << " Info::siweyhmi: run over real data " << std::endl;
      task->SetUseCouts(kFALSE);
      task->SetDefaultTrackCuts(kTRUE);
      task->SetDefaultEventCuts(kTRUE);
      task->SetRunOnGrid(kTRUE);
      task->SetRunNumberForExpecteds(0);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      //
      task->SetLeadingJetCut(5);
      task->SetFillJetsBG(kFALSE);
      task->SetFillFastJet(kTRUE);
      //Set these in the wagon configuration CHANGE
      task->SetJetMinPtSub(pT_sub_min);
      task->SetJetAreaCut(0.0); // sets so it saves 1 out of every n events inclusive = 400. Jets = 40 w/ 40 GeV min jet requirement
      //
    }
    break;
    case 1:{
      std::cout << " SETTING TYPE = " << settingType << " Info::siweyhmi: run over MC data " << std::endl;
      task->SetIsMCtrue(kTRUE);
      task->SetUseCouts(kTRUE);
      task->SetMCTrackOriginType(0);   // 0:full scan, 1: prim
      task->SetUsePtCut(1); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      // acceptance
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 2;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = { 0.6, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 2.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);

      //
      task->SetLeadingJetCut(5);
      task->SetFillJetsBG(kFALSE);
      task->SetFillFastJet(kTRUE);

      //Set these in the wagon configuration CHANGE
      task->SetJetMinPtSub(pT_sub_min);
      task->SetPercentageOfEvents(0); //sets so it saves 1 out of every n events inclusive = 400. Jets = 40 w/ 40 GeV min jet requirement
      //
    }

  }

  // Finally initialize the task after physics selection
  task->Initialize();
  return task;

}
// ____________________________________________________________________________________________
void SetDefaultsJets(AliAnalysisJetHadro *defaultTask, Int_t year, TString periodName, Int_t passIndex)
{

  // Setters for the eta momentum dEdx and centrality bins
  std::cout << " Info::siweyhmi: ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " Info::siweyhmi: ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " Info::siweyhmi: ------------------- Set default settings for the task object ------------------------ " << std::endl;
  std::cout << " Info::siweyhmi: ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " Info::siweyhmi: ------------------------------------------------------------------------------------- " << std::endl;

  defaultTask->SetNSettings(22);
  defaultTask->SetCorrectForMissCl(0);
  defaultTask->SetYear(year);
  defaultTask->SetPeriodName(periodName);
  defaultTask->SetPassIndex(passIndex);
  defaultTask->SetSampleDeDxUpperEdge(400.);
  defaultTask->SetDeDxBinWidth(2.5);
  defaultTask->SetDeDxLowerEdge(20.);
  defaultTask->SetDeDxUpperEdge(2020.);
  defaultTask->SetNEtabins(16);
  defaultTask->SetEtaLowerEdge(-0.8);
  defaultTask->SetEtaUpperEdge( 0.8);
  defaultTask->SetNMomBins(250);
  defaultTask->SetMomLowerEdge(0.2);
  defaultTask->SetMomUpperEdge(5.2);
  defaultTask->SetNGenprotonBins(100);
  defaultTask->SetPercentageOfEvents(0);
  defaultTask->SetLeadingJetCut(2);

  // DEFAULT SETTINGS
  const Int_t tmpCentbins  = 14;
  const Int_t tmpEtaBinsMC = 3;
  const Int_t tmpMomBinsMC = 4;
  Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80,85,90,95,100};
  Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.5,-0.8,-1.};
  Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.5, 0.8, 1.};
  Float_t tmppDownArr[tmpMomBinsMC] = { 0.2, 0.6, 0.2, 0.6};
  Float_t tmppUpArr[tmpMomBinsMC]   = { 1.5, 1.5, 1.8, 1.8};

  // centrality binning and Eta Momentum Scans for MC
  defaultTask->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
  defaultTask->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
  defaultTask->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);

  // Boolians which are by default === ON ===
  defaultTask->SetRunOnGrid(kFALSE);
  defaultTask->SetIsMCtrue(kFALSE);
  defaultTask->SetIncludeITScuts(kTRUE);
  defaultTask->SetFillJetsBG(kFALSE);
  defaultTask->SetFillFastJet(kFALSE);
  defaultTask->SetUsePtCut(1);
  defaultTask->SetMCTrackOriginType(1);   // 0:full scan, 1: prim
  defaultTask->SetRapidityType(0);      // 0:pseudorapidity, 1: rapidity
  defaultTask->SetEffMatrixNSigmasTOF(2.5);

  // Extra Boolians which are by default === OFF ===
  defaultTask->SetDeDxCheck(kFALSE);
  defaultTask->SetFillOnlyHists(kFALSE);
  defaultTask->SetRunFastSimulation(kFALSE);
  defaultTask->SetFillHigherMomentsMCclosure(kFALSE);
  defaultTask->SetIncludeTOF(kFALSE);
  defaultTask->SetUseCouts(kFALSE);
  defaultTask->SetFillEventInfo(kFALSE);
  defaultTask->SetDefaultTrackCuts(kTRUE);
  defaultTask->SetDefaultEventCuts(kFALSE);

  // Setters for the systematic uncertainty checks
  defaultTask->SetSystCentEstimator(0);

}
