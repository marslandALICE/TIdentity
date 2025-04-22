
void SetDefaults(AliAnalysisTaskTIdentityPID *defaultTask, Int_t year, TString periodName, Int_t passIndex);
TTree *GetLookUpTable(Bool_t runOnGrid, Int_t index);
vector<THnF*> GetEffMatrixObjects(Bool_t runOnGrid, TString filename);
vector<TH2F*> GetQvecCorrObjects(Bool_t runOnGrid, TString filename);
//
//
AliAnalysisTaskTIdentityPID* Config_marsland_TIdentityPID(Bool_t getFromAlien, Int_t settingType, Int_t year, TString periodName, Int_t passIndex, Int_t lookUpTableIndex, TString combinedName) {
  //
  // Configuration file for the AliAnalysisTaskTIdentityPID.cxx class
  //
  AliAnalysisTaskTIdentityPID *task = new AliAnalysisTaskTIdentityPID(combinedName);
  SetDefaults(task,year,periodName,passIndex);
  std::cout << " Info::Config::marsland: ===== Physics selection ? = " << std::endl;
  if (year==2010) task->SelectCollisionCandidates(AliVEvent::kMB);   // select minimum bias events for LHC10h
  if (year==2015) task->SelectCollisionCandidates(AliVEvent::kINT7); // select minimum bias events for LHC15o
  if (year==2018) task->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral);
  //
  // Get the lookup table
  TTree *lookUpTree=NULL;
  if (lookUpTableIndex>0 && settingType>199) lookUpTree = GetLookUpTable(getFromAlien,lookUpTableIndex);

  std::cout << " Info::marsland: ===== In the Config --> Running with year = ";
  std::cout << year << " --- period name = " << periodName << " --- pass = " << passIndex << " --- lookUpTableIndex = " << lookUpTableIndex << " --- settingType = " << settingType << std::endl;
  // Other Specific settings

  switch (settingType) {
    //
    // ====================================================================================
    // ============================= Real Data Settings ===================================
    // ====================================================================================
    //
    case 0:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Real data " << std::endl;
      // Real data settings
      if( (passIndex==3 && periodName.Contains("18")) || (passIndex==2 && periodName.Contains("15")) ) {
        task->SetDefaultEventCuts(kTRUE);
      }
      task->SetDownsampleTrees(kTRUE);
      task->SetFillQvectorHists(kFALSE);
      task->SetApplyQVectorCorr(kFALSE);
      task->SetCollisionType(0); // 0 for PbPb, 1 for pp
      task->SetTaskSelection(2); // 0; both jet+net-p, 1: only jet, 2: only net-p, 3: only cutbased
      task->SetDownscalingFactor(0.001);
      task->SetFillJetsBG(2); // 0: no PID info for const + no BG, 1: fill BG tree, 2: no BG tree + PID info
      task->SetUseCouts(kFALSE);
      task->SetFillDebug(kFALSE);
      //
      task->SetFillEventInfo(kTRUE);
      task->SetFillDistributions(kTRUE);
      task->SetFillArmPodTree(kTRUE);
      task->SetV0InvMassHists(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      if (periodName.Contains("18q")) task->SetRunNumberForExpecteds(296433);
      else if (periodName.Contains("18r")) task->SetRunNumberForExpecteds(296749);
      else if (periodName.Contains("15o")) task->SetRunNumberForExpecteds(246087);
      //
      // acceptance & settings
      task->SetNSettings(1);
      task->SetSettings({0});
      //
      // task->SetNSettings(13);
      // task->SetSettings({0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 16, 17});
      //
      const Int_t tmpCentbins  = 14;
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 13;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = {0.3, 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.6, 1.5, 1.5, 1.5, 3.0, 0.4};
      Float_t tmppUpArr[tmpMomBinsMC]   = {1.5, 2.0, 3.0, 5.0, 1.5, 2.0, 3.0, 5.0, 3.0, 4.0, 5.0, 5.0, 5.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);

      std::vector<Double_t> effMatrixMomBins = {
        0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38,
        0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58,
        0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78,
        0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98,
        1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.18,
        1.20, 1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 1.38,
        1.40, 1.44, 1.48, 1.52, 1.56, 1.60, 1.64, 1.68, 1.72, 1.76,
        1.80, 1.84, 1.88, 1.92, 1.96, 2.00, 2.04, 2.08, 2.12, 2.16,
        2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00,
        4.20, 4.40, 4.60, 4.80, 5.00, 5.20
      };
      std::vector<Double_t> effMatrixCentBins = {0,5,10,20,30,40,50,60,70,80,90};
      task->SetEffMatrixMomBins(effMatrixMomBins);
      task->SetEffMatrixCentBins(effMatrixCentBins);
      task->SetNSigmaTPC({3.0, 3.5, 3.5});
      task->SetNSigmaTOF({-2.5, -3.0, -3.5}, {2.5, 2.5, 2.5});
      //
      // Read eff matrix
      // vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kFALSE, "AnalysisResults_hists.root");
      vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kTRUE, "AnalysisResults_hists.root");
      task->SetEffMatrixObjects(effMatrixObjects[0], effMatrixObjects[1], effMatrixObjects[2], effMatrixObjects[3]);

      vector<TH2F*> qvecCorrObjects = GetQvecCorrObjects(kTRUE, "AnalysisResults_qvecs.root");
      task->SetQvecCorrObjects(qvecCorrObjects[0], qvecCorrObjects[1], qvecCorrObjects[2], qvecCorrObjects[3]);

    }
    break;
    case 1:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Real data; only fill Qvectors " << std::endl;
      // Real data settings
      if( (passIndex==3 && periodName.Contains("18")) || (passIndex==2 && periodName.Contains("15")) ) {
        task->SetDefaultEventCuts(kTRUE);
      }
      task->SetFillQvectorHists(kTRUE);
      task->SetApplyQVectorCorr(kFALSE);
      task->SetCollisionType(0); // 0 for PbPb, 1 for pp
      task->SetTaskSelection(2); // 0; both jet+net-p, 1: only jet, 2: only net-p
      task->SetDownscalingFactor(0.001);
      task->SetFillJetsBG(2); // 0: no PID info for const + no BG, 1: fill BG tree, 2: no BG tree + PID info
      task->SetUseCouts(kFALSE);
      task->SetFillDebug(kFALSE);
      task->SetFillEventInfo(kTRUE);
      //
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      if (periodName.Contains("18q")) task->SetRunNumberForExpecteds(296433);
      else if (periodName.Contains("18r")) task->SetRunNumberForExpecteds(296749);
      else if (periodName.Contains("15o")) task->SetRunNumberForExpecteds(246087);
      //
      // acceptance & settings
      task->SetNSettings(6);
      task->SetSettings({0, 1, 2, 3, 4, 5});
      //
      const Int_t tmpCentbins  = 14;
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 8;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = {0.4, 0.4, 0.4, 0.6, 0.6, 1.5, 1.5, 1.5};
      Float_t tmppUpArr[tmpMomBinsMC]   = {1.0, 1.5, 2.0, 1.5, 2.0, 3.0, 4.0, 5.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);

      std::vector<Double_t> effMatrixMomBins = {
        0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38,
        0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58,
        0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78,
        0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98,
        1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.18,
        1.20, 1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 1.38,
        1.40, 1.44, 1.48, 1.52, 1.56, 1.60, 1.64, 1.68, 1.72, 1.76,
        1.80, 1.84, 1.88, 1.92, 1.96, 2.00, 2.04, 2.08, 2.12, 2.16,
        2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00,
        4.20, 4.40, 4.60, 4.80, 5.00, 5.20
      };
      std::vector<Double_t> effMatrixCentBins = {0,5,10,20,30,40,50,60,70,80,90};
      task->SetEffMatrixMomBins(effMatrixMomBins);
      task->SetEffMatrixCentBins(effMatrixCentBins);
      task->SetNSigmaTPC({3.0, 3.5, 3.5});
      task->SetNSigmaTOF({-2.5, -3.0, -3.5}, {2.5, 2.5, 2.5});
      //
      // Read eff matrix
      // vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kFALSE, "AnalysisResults_hists.root");
      vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kTRUE, "AnalysisResults_hists.root");
      task->SetEffMatrixObjects(effMatrixObjects[0], effMatrixObjects[1], effMatrixObjects[2], effMatrixObjects[3]);

    }
    break;
    case 5:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Real data " << std::endl;
      // Real data settings
      if( (passIndex==3 && periodName.Contains("18")) || (passIndex==2 && periodName.Contains("15")) ) {
        task->SetDefaultEventCuts(kTRUE);
      }
      task->SetFillQvectorHists(kFALSE);
      task->SetApplyQVectorCorr(kTRUE);
      task->SetCollisionType(0); // 0 for PbPb, 1 for pp
      task->SetTaskSelection(3); // 0; both jet+net-p, 1: only jet, 2: only net-p, 3: only cutbased
      task->SetDownscalingFactor(0.001);
      task->SetFillJetsBG(2); // 0: no PID info for const + no BG, 1: fill BG tree, 2: no BG tree + PID info
      task->SetUseCouts(kFALSE);
      task->SetFillDebug(kFALSE);
      task->SetFillArmPodTree(kFALSE);
      //
      task->SetFillEventInfo(kTRUE);
      task->SetFillDistributions(kTRUE);

      task->SetV0InvMassHists(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      if (periodName.Contains("18q")) task->SetRunNumberForExpecteds(296433);
      else if (periodName.Contains("18r")) task->SetRunNumberForExpecteds(296749);
      else if (periodName.Contains("15o")) task->SetRunNumberForExpecteds(246087);
      //
      // acceptance & settings
      // task->SetNSettings(6);
      // task->SetSettings({0, 5, 7, 11, 16, 17});
      //
      task->SetNSettings(13);
      task->SetSettings({0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 16, 17});
      //
      const Int_t tmpCentbins  = 14;
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 11;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = {0.3, 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.6, 1.5, 1.5, 1.5};
      Float_t tmppUpArr[tmpMomBinsMC]   = {1.5, 2.0, 3.0, 5.0, 1.5, 2.0, 3.0, 5.0, 3.0, 4.0, 5.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);

      std::vector<Double_t> effMatrixMomBins = {
        0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38,
        0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58,
        0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78,
        0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98,
        1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.18,
        1.20, 1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 1.38,
        1.40, 1.44, 1.48, 1.52, 1.56, 1.60, 1.64, 1.68, 1.72, 1.76,
        1.80, 1.84, 1.88, 1.92, 1.96, 2.00, 2.04, 2.08, 2.12, 2.16,
        2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00,
        4.20, 4.40, 4.60, 4.80, 5.00, 5.20
      };
      std::vector<Double_t> effMatrixCentBins = {0,5,10,20,30,40,50,60,70,80,90};
      task->SetEffMatrixMomBins(effMatrixMomBins);
      task->SetEffMatrixCentBins(effMatrixCentBins);
      task->SetNSigmaTPC({3.0, 3.5, 3.5});
      task->SetNSigmaTOF({-2.5, -3.0, -3.5}, {2.5, 2.5, 2.5});
      //
      // Read eff matrix
      // vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kFALSE, "AnalysisResults_hists.root");
      vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kTRUE, "AnalysisResults_hists.root");
      task->SetEffMatrixObjects(effMatrixObjects[0], effMatrixObjects[1], effMatrixObjects[2], effMatrixObjects[3]);

      vector<TH2F*> qvecCorrObjects = GetQvecCorrObjects(kTRUE, "AnalysisResults_qvecs.root");
      task->SetQvecCorrObjects(qvecCorrObjects[0], qvecCorrObjects[1], qvecCorrObjects[2], qvecCorrObjects[3]);

    }
    break;
    case 11:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Real data " << std::endl;
      // Real data settings
      if( (passIndex==3 && periodName.Contains("18")) || (passIndex==2 && periodName.Contains("15")) ) {
        task->SetDefaultEventCuts(kTRUE);
      }
      // task->SetNSettings(4);
      // task->SetSettings({0, 1, 16, 17});
      task->SetCollisionType(1); // 0 for PbPb, 1 for pp
      task->SetTaskSelection(0); // 0; both jet+net-p, 1: only jet, 2: only net-p
      task->SetNSettings(1);
      task->SetSettings({0});
      task->SetFillJetsBG(2);
      task->SetUseCouts(kFALSE);
      //
      task->SetFillEventInfo(kTRUE);
      task->SetFillDistributions(kTRUE);
      task->SetFillArmPodTree(kTRUE);
      task->SetV0InvMassHists(kTRUE);
      const Int_t tmpCentbins = 14;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;

      if (periodName.Contains("18q"))
        task->SetRunNumberForExpecteds(296433);
      else if (periodName.Contains("18r"))
        task->SetRunNumberForExpecteds(296749);
      else if (periodName.Contains("15o"))
        task->SetRunNumberForExpecteds(246087);
    }
    break;
    //
    // ====================================================================================
    // =================================== Fill only EffMatrix ============================
    // ====================================================================================
    //
    case 40:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Fill only EffMatrix " << std::endl;
      //
      // Main command of this case
      task->SetFillQvectorHists(kTRUE);
      task->SetFillEffMatrix(kTRUE); // It fills only eff histograms
      task->SetFillTracksMCgen(kFALSE);
      task->SetFillEventInfo(kTRUE);
      task->SetUseCouts(kTRUE);
      task->SetCollisionType(0); // 0 for PbPb, 1 for pp
      task->SetIsMCtrue(kTRUE);
      task->SetFillDebug(kFALSE);
      //
      std::cout << "period and pass = " << periodName << "    " << passIndex << std::endl;
      if( (passIndex==3) || (passIndex==2) ) {
        task->SetDefaultEventCuts(kTRUE);
        std::cout << " special settings for 18q pass3 and 15o pass2 " << std::endl;
      }
      task->SetMCTrackOriginType(0);   // 0:full scan, 1: prim
      task->SetRapidityType(0);      // 0: pseudorapidity, 1: rapidity
      task->SetUsePtCut(0);          // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      task->SetIncludeITScuts(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      //
      // acceptance & settings
      task->SetNSettings(13);
      task->SetSettings({0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 16, 17});
      //
      const Int_t tmpCentbins  = 14;
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 15;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = {0.3, 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.6, 1.5, 1.5, 1.5, 3.0, 0.4, 0.4,  0.4 };
      Float_t tmppUpArr[tmpMomBinsMC]   = {1.5, 2.0, 3.0, 5.0, 1.5, 2.0, 3.0, 5.0, 3.0, 4.0, 5.0, 5.0, 5.0, 10.0, 20.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);

      std::vector<Double_t> effMatrixMomBins = {
        0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38,
        0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58,
        0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78,
        0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98,
        1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.18,
        1.20, 1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 1.38,
        1.40, 1.44, 1.48, 1.52, 1.56, 1.60, 1.64, 1.68, 1.72, 1.76,
        1.80, 1.84, 1.88, 1.92, 1.96, 2.00, 2.04, 2.08, 2.12, 2.16,
        2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00,
        4.20, 4.40, 4.60, 4.80, 5.00, 5.20
      };
      std::vector<Double_t> effMatrixCentBins = {0,5,10,20,30,40,50,60,70,80,90};
      task->SetEffMatrixMomBins(effMatrixMomBins);
      task->SetEffMatrixCentBins(effMatrixCentBins);
      task->SetNSigmaTPC({3.0, 3.5, 3.5});
      task->SetNSigmaTOF({-2.5, -3.0, -3.5}, {2.5, 2.5, 2.5});

      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      //
      // baryons to be included for netbaryon analysis --> light and strange baryons
      // {p,n,delta++,delta+,delta0,delta-,Lambda,}
      const Int_t tmpNbaryons = 18;
      Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112,2224,2214,2114,1114,3122,3222,3212,3112,3224,3214,3114,3322,3312,3324,3314,3334};
      task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);

      vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kTRUE, "AnalysisResults_hists.root");
      task->SetEffMatrixObjects(effMatrixObjects[0], effMatrixObjects[1], effMatrixObjects[2], effMatrixObjects[3]);

    }
    break;
    //
    // ====================================================================================
    // =================================== FULL MC  =======================================
    // ====================================================================================
    //
    case 50:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Full MC PbPb " << std::endl;
      //
      // Main task of this case
      task->SetFillEffMatrix(kFALSE);
      task->SetCollisionType(0); // 0 for PbPb, 1 for pp
      task->SetDownscalingFactor(0.001);
      task->SetFillTracksMCgen(kTRUE);
      task->SetFillJetsBG(2);
      task->SetTaskSelection(2); // 0; both jet+net-p, 1: only jet, 2: only net-p
      task->SetIsMCtrue(kTRUE);
      task->SetSisterCheck(0);
      task->SetEtaUpperEdge(0.8);
      //
      // mostly for debugging
      task->SetFillDebug(kFALSE);
      task->SetApplyQVectorCorr(kFALSE);
      task->SetFillQvectorHists(kTRUE);
      task->SetUseCouts(kTRUE);
      task->SetFillTracksMCgen(kTRUE);
      task->SetFillArmPodTree(kTRUE);
      task->SetV0InvMassHists(kTRUE);
      //
      std::cout << "period and pass = " << periodName << "    " << passIndex << std::endl;
      if( (passIndex==3) || (passIndex==2) ) {
        task->SetDefaultEventCuts(kTRUE);
        std::cout << " special settings for 18q pass3 and 15o pass2 " << std::endl;
      }
      task->SetFillResonances(kTRUE);
      task->SetMCTrackOriginType(0);   // 0:full scan, 1: prim
      task->SetCorrectForMissCl(0);
      task->SetRapidityType(0);      // 0: pseudorapidity, 1: rapidity
      task->SetUsePtCut(0);          // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      task->SetFillTreeMC(kTRUE);
      task->SetFillEventInfo(kTRUE);
      task->SetIncludeITScuts(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      //
      // acceptance & settings
      task->SetNSettings(4);
      task->SetSettings({0, 2, 4, 5});
      //
      const Int_t tmpCentbins  = 14;
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 8;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = {0.4, 0.4, 0.4, 0.6, 0.6, 1.5, 1.5, 1.5};
      Float_t tmppUpArr[tmpMomBinsMC]   = {1.0, 1.5, 2.0, 1.5, 2.0, 3.0, 4.0, 5.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);

      std::vector<Double_t> effMatrixMomBins = {
        0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38,
        0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58,
        0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78,
        0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98,
        1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.18,
        1.20, 1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 1.38,
        1.40, 1.44, 1.48, 1.52, 1.56, 1.60, 1.64, 1.68, 1.72, 1.76,
        1.80, 1.84, 1.88, 1.92, 1.96, 2.00, 2.04, 2.08, 2.12, 2.16,
        2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00,
        4.20, 4.40, 4.60, 4.80, 5.00, 5.20
      };
      std::vector<Double_t> effMatrixCentBins = {0,5,10,20,30,40,50,60,70,80,90};
      task->SetEffMatrixMomBins(effMatrixMomBins);
      task->SetEffMatrixCentBins(effMatrixCentBins);
      task->SetNSigmaTPC({3.0, 3.5, 3.5});
      task->SetNSigmaTOF({-2.5, -3.0, -3.5}, {2.5, 2.5, 2.5});
      //
      // resonances to exclude e.g. {"p","n","delta++","delta+","delta0","delta-","Lambda"}
      const Int_t tmpNresonances = 4;
      TString tmpResArr[tmpNresonances] = {"rho","phi","omega","eta"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      //
      // baryons to be included for netbaryon analysis --> light and strange baryons
      const Int_t tmpNbaryons = 18; //  {p,n,delta++,delta+,delta0,delta-,lambda,sigmas, xis, omega-}
      Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112,2224,2214,2114,1114,3122,3222,3212,3112,3224,3214,3114,3322,3312,3324,3314,3334};
      task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);

      // vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kFALSE, "grid_output/AnalysisResults.root");
      vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kTRUE, "AnalysisResults_hists.root");
      task->SetEffMatrixObjects(effMatrixObjects[0], effMatrixObjects[1], effMatrixObjects[2], effMatrixObjects[3]);

      vector<TH2F*> qvecCorrObjects = GetQvecCorrObjects(kTRUE, "AnalysisResults_qvecs.root");
      task->SetQvecCorrObjects(qvecCorrObjects[0], qvecCorrObjects[1], qvecCorrObjects[2], qvecCorrObjects[3]);

    }
    break;
    case 51:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Full MC pp " << std::endl;
      //
      // Main task of this case
      task->SetCollisionType(1); // 0 for PbPb, 1 for pp
      //
      task->SetFillJetsBG(2);
      task->SetTaskSelection(0); // 0; both jet+net-p, 1: only jet, 2: only net-p
      task->SetIsMCtrue(kTRUE);
      task->SetUseCouts(kFALSE);
      task->SetSisterCheck(0);
      //
      std::cout << "period and pass = " << periodName << "    " << passIndex << std::endl;
      if( (passIndex==3) || (passIndex==2) ) {
        task->SetDefaultEventCuts(kTRUE);
        std::cout << " special settings for 18q pass3 and 15o pass2 " << std::endl;
      }
      task->SetFillArmPodTree(kFALSE);
      task->SetFillResonances(kTRUE);
      task->SetMCTrackOriginType(0);   // 0:full scan, 1: prim
      task->SetCorrectForMissCl(0);
      task->SetRapidityType(0);      // 0: pseudorapidity, 1: rapidity
      task->SetUsePtCut(0);          // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      task->SetFillTreeMC(kTRUE);
      task->SetFillEventInfo(kTRUE);
      task->SetIncludeITScuts(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      //
      // acceptance & settings
      // task->SetNSettings(4);
      // task->SetSettings({0, 1, 16, 17});
      task->SetNSettings(1);
      task->SetSettings({0});
      //
      const Int_t tmpCentbins  = 14;
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 15;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = {0.3, 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.6, 1.5, 1.5, 1.5, 3.0, 0.4, 0.4,  0.4 };
      Float_t tmppUpArr[tmpMomBinsMC]   = {1.5, 2.0, 3.0, 5.0, 1.5, 2.0, 3.0, 5.0, 3.0, 4.0, 5.0, 5.0, 5.0, 10.0, 20.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);

      std::vector<Double_t> effMatrixMomBins = {
        0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38,
        0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58,
        0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78,
        0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98,
        1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.18,
        1.20, 1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 1.38,
        1.40, 1.44, 1.48, 1.52, 1.56, 1.60, 1.64, 1.68, 1.72, 1.76,
        1.80, 1.84, 1.88, 1.92, 1.96, 2.00, 2.04, 2.08, 2.12, 2.16,
        2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00,
        4.20, 4.40, 4.60, 4.80, 5.00, 5.20
      };
      std::vector<Double_t> effMatrixCentBins = {0,5,10,20,30,40,50,60,70,80,90};
      task->SetEffMatrixMomBins(effMatrixMomBins);
      task->SetEffMatrixCentBins(effMatrixCentBins);
      task->SetNSigmaTPC({3.0, 3.5, 3.5});
      task->SetNSigmaTOF({-2.5, -3.0, -3.5}, {2.5, 2.5, 2.5});

      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      //
      // baryons to be included for netbaryon analysis --> light and strange baryons
      // {p,n,delta++,delta+,delta0,delta-,Lambda,}
      const Int_t tmpNbaryons = 18;
      Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112,2224,2214,2114,1114,3122,3222,3212,3112,3224,3214,3114,3322,3312,3324,3314,3334};
      task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);

      vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kTRUE, "AnalysisResults_hists.root");
      task->SetEffMatrixObjects(effMatrixObjects[0], effMatrixObjects[1], effMatrixObjects[2], effMatrixObjects[3]);

    }
    break;
    //
    //
    //
    //
    //
    case 60:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Full MC PbPb " << std::endl;
      //
      // Main task of this case
      task->SetUseAODsForMC(kTRUE);
      task->SetFillEffMatrix(kFALSE);
      task->SetCollisionType(0); // 0 for PbPb, 1 for pp
      task->SetDownscalingFactor(0.001);
      task->SetFillTracksMCgen(kTRUE);
      task->SetFillJetsBG(2);
      task->SetTaskSelection(2); // 0; both jet+net-p, 1: only jet, 2: only net-p
      task->SetIsMCtrue(kTRUE);
      task->SetSisterCheck(0);
      task->SetEtaUpperEdge(0.8);
      //
      // mostly for debugging
      task->SetFillDebug(kTRUE);
      task->SetApplyQVectorCorr(kFALSE);
      task->SetFillQvectorHists(kTRUE);
      task->SetUseCouts(kTRUE);
      task->SetFillTracksMCgen(kTRUE);
      task->SetFillArmPodTree(kTRUE);
      task->SetV0InvMassHists(kTRUE);
      //
      std::cout << "period and pass = " << periodName << "    " << passIndex << std::endl;
      if( (passIndex==3) || (passIndex==2) ) {
        task->SetDefaultEventCuts(kTRUE);
        std::cout << " special settings for 18q pass3 and 15o pass2 " << std::endl;
      }
      task->SetFillResonances(kTRUE);
      task->SetMCTrackOriginType(0);   // 0:full scan, 1: prim
      task->SetCorrectForMissCl(0);
      task->SetRapidityType(0);      // 0: pseudorapidity, 1: rapidity
      task->SetUsePtCut(0);          // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      task->SetFillTreeMC(kTRUE);
      task->SetFillEventInfo(kTRUE);
      task->SetIncludeITScuts(kTRUE);
      task->fEventCuts.fUseVariablesCorrelationCuts = true;
      //
      // acceptance & settings
      task->SetNSettings(4);
      task->SetSettings({0, 2, 4, 5});
      //
      const Int_t tmpCentbins  = 14;
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 8;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = {0.4, 0.4, 0.4, 0.6, 0.6, 1.5, 1.5, 1.5};
      Float_t tmppUpArr[tmpMomBinsMC]   = {1.0, 1.5, 2.0, 1.5, 2.0, 3.0, 4.0, 5.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);

      std::vector<Double_t> effMatrixMomBins = {
        0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38,
        0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58,
        0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78,
        0.80, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98,
        1.00, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.18,
        1.20, 1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 1.38,
        1.40, 1.44, 1.48, 1.52, 1.56, 1.60, 1.64, 1.68, 1.72, 1.76,
        1.80, 1.84, 1.88, 1.92, 1.96, 2.00, 2.04, 2.08, 2.12, 2.16,
        2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 4.00,
        4.20, 4.40, 4.60, 4.80, 5.00, 5.20
      };
      std::vector<Double_t> effMatrixCentBins = {0,5,10,20,30,40,50,60,70,80,90};
      task->SetEffMatrixMomBins(effMatrixMomBins);
      task->SetEffMatrixCentBins(effMatrixCentBins);
      task->SetNSigmaTPC({3.0, 3.5, 3.5});
      task->SetNSigmaTOF({-2.5, -3.0, -3.5}, {2.5, 2.5, 2.5});
      //
      // resonances to exclude e.g. {"p","n","delta++","delta+","delta0","delta-","Lambda"}
      const Int_t tmpNresonances = 4;
      TString tmpResArr[tmpNresonances] = {"rho","phi","omega","eta"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      //
      // baryons to be included for netbaryon analysis --> light and strange baryons
      const Int_t tmpNbaryons = 18; //  {p,n,delta++,delta+,delta0,delta-,lambda,sigmas, xis, omega-}
      Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112,2224,2214,2114,1114,3122,3222,3212,3112,3224,3214,3114,3322,3312,3324,3314,3334};
      task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);

      // vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kFALSE, "grid_output/AnalysisResults.root");
      vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kTRUE, "AnalysisResults_hists.root");
      task->SetEffMatrixObjects(effMatrixObjects[0], effMatrixObjects[1], effMatrixObjects[2], effMatrixObjects[3]);
  
      vector<TH2F*> qvecCorrObjects = GetQvecCorrObjects(kTRUE, "AnalysisResults_qvecs.root");
      task->SetQvecCorrObjects(qvecCorrObjects[0], qvecCorrObjects[1], qvecCorrObjects[2], qvecCorrObjects[3]);

    }
    break;
    //
    // ====================================================================================
    // ================================ FastGen  Settings =================================
    // ====================================================================================
    //
    case 200:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Fast Gen volume fluctuations " << std::endl;
      task->SetTaskSelection(2); // 0; both jet+net-p, 1: only jet, 2: only net-p
      task->SetFillResonances(kTRUE);
      task->SetCollisionType(0); // 0 for PbPb, 1 for pp
      task->SetRunFastSimulation(kTRUE);
      task->SetFillTracksMCgen(kTRUE);
      task->SetFillDebug(kFALSE);
      task->SetDownscalingFactor(0.001);
      task->SetUseCouts(kTRUE);
      task->SetIsMCtrue(kTRUE);
      task->SetUsePtCut(0); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      task->SetSisterCheck(0);
      task->SetRapidityType(0);      // 0: pseudorapidity, 1: rapidity
      //
      // Main task of this case
      // const Int_t tmpCentbins  = 14;
      // const Int_t tmpEtaBinsMC = 8;
      // const Int_t tmpMomBinsMC = 15;
      // Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      // Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      // Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      // Float_t tmppDownArr[tmpMomBinsMC] = {0.3, 0.3, 0.3, 0.3, 0.6, 0.6, 0.6, 0.6, 1.5, 1.5, 1.5, 3.0, 0.4, 0.4,  0.4 };
      // Float_t tmppUpArr[tmpMomBinsMC]   = {1.5, 2.0, 3.0, 5.0, 1.5, 2.0, 3.0, 5.0, 3.0, 4.0, 5.0, 5.0, 5.0, 10.0, 20.0};
      const Int_t tmpCentbins  = 14;
      const Int_t tmpEtaBinsMC = 8;
      const Int_t tmpMomBinsMC = 6;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,65,70,75,80,85,90};
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
      Float_t tmppDownArr[tmpMomBinsMC] = {0.4, 0.4, 0.4, 0.6, 0.6, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = {1.0, 1.5, 2.0, 1.0, 1.5, 2.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      //
      // resonances to exclude
      const Int_t tmpNresonances = 4;
      TString tmpResArr[tmpNresonances] = {"rho","phi","omega","eta"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      //
      // baryons to be included for netbaryon analysis --> light and strange baryons
      const Int_t tmpNbaryons = 18; //  {p,n,delta++,delta+,delta0,delta-,lambda,sigmas, xis, omega-}
      Int_t tmpBaryonArr[tmpNbaryons] = {2212,2112,2224,2214,2114,1114,3122,3222,3212,3112,3224,3214,3114,3322,3312,3324,3314,3334};
      task->SetMCBaryonArray(tmpNbaryons,tmpBaryonArr);

      // vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kFALSE, "AnalysisResults_hists.root");
      vector<THnF*> effMatrixObjects = GetEffMatrixObjects(kTRUE, "AnalysisResults_hists.root");
      task->SetEffMatrixObjects(effMatrixObjects[0], effMatrixObjects[1], effMatrixObjects[2], effMatrixObjects[3]);

    }
    break;
    case 201:{
      std::cout << " SETTING TYPE = " << settingType << " Info::marsland: Fast Gen for scoping document " << std::endl;
      task->SetFillJetsBG(2);
      task->SetFillResonances(kTRUE);
      task->SetCollisionType(0); // 0 for PbPb, 1 for pp
      task->SetRunOnGrid(kFALSE);
      task->SetRunFastSimulation(kTRUE);
      task->SetDownscalingFactor(0.001);
      task->SetUseCouts(kTRUE);
      task->SetIsMCtrue(kTRUE);
      task->SetUsePtCut(1); // 0: tpc momcut, 1: vertex momcut, 2: pT cut
      // eta bin scan
      const Int_t tmpEtaBinsMC = 15;
      Float_t tmpetaDownArr[tmpEtaBinsMC] = {-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7,-0.8,-1.,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0};
      Float_t tmpetaUpArr[tmpEtaBinsMC]   = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1., 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
      task->SetMCEtaScanArray(tmpEtaBinsMC, tmpetaDownArr, tmpetaUpArr);
      // mom bin scan
      const Int_t tmpMomBinsMC = 3;
      Float_t tmppDownArr[tmpMomBinsMC] = {0.6,1.5, 0.6};
      Float_t tmppUpArr[tmpMomBinsMC]   = {1.5,10.0,10.0};
      task->SetMCMomScanArray(tmpMomBinsMC, tmppDownArr,   tmppUpArr);
      // cent bins
      const Int_t tmpCentbins  = 12;
      Float_t tmpfxCentBins[tmpCentbins] = {0,5,10,20,30,40,50,60,70,80,90,100};
      task->SetCentralityBinning(tmpCentbins,tmpfxCentBins);
      // resonances to exclude
      const Int_t tmpNresonances = 1;
      TString tmpResArr[tmpNresonances] = {"xxx"};
      task->SetMCResonanceArray(tmpNresonances,tmpResArr);
      // const Int_t tmpNresonances = 5;
      // TString tmpResArr[tmpNresonances] = {"rho","phi","omega","eta","Delta"};
      // task->SetMCResonanceArray(tmpNresonances,tmpResArr);
    }
    break;

  }

  // Finally initialize the task after physics selection
  task->Initialize();
  return task;

}
// ____________________________________________________________________________________________
void SetDefaults(AliAnalysisTaskTIdentityPID *defaultTask, Int_t year, TString periodName, Int_t passIndex)
{

  // Setters for the eta momentum dEdx and centrality bins
  std::cout << " Info::marsland: ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " Info::marsland: ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " Info::marsland: ------------------- Set default settings for the task object ------------------------ " << std::endl;
  std::cout << " Info::marsland: ------------------------------------------------------------------------------------- " << std::endl;
  std::cout << " Info::marsland: ------------------------------------------------------------------------------------- " << std::endl;

  defaultTask->SetSisterCheck(0);
  defaultTask->SetNSettings(4);
  defaultTask->SetSettings({0, 1, 16, 17});
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
  defaultTask->SetDownscalingFactor(0.001);
  defaultTask->SetFillQvectorHists(kFALSE);
  defaultTask->SetApplyQVectorCorr(kFALSE);
  defaultTask->SetDownsampleTrees(kFALSE);
  defaultTask->SetUseAODsForMC(kFALSE);


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

  std::vector<Double_t> effMatrixMomBins = {0.2};
  std::vector<Double_t> effMatrixCentBins = {0.};
  std::vector<Double_t> effMatrixEtaBins = {-0.8};

  for (Int_t iMom = 1; iMom <= 50; iMom++) {
    effMatrixMomBins.push_back(effMatrixMomBins.back() + (5.2 - 0.2) / 50.);
  }
  for (Int_t iCent = 1; iCent <= 10; iCent++) {
    effMatrixCentBins.push_back(effMatrixCentBins.back() + 10.);
  }
  for (Int_t iEta = 1; iEta <= 16; iEta++) {
    effMatrixEtaBins.push_back(effMatrixEtaBins.back() + 0.1);
  }

  defaultTask->SetEffMatrixMomBins(effMatrixMomBins);
  defaultTask->SetEffMatrixCentBins(effMatrixCentBins);
  defaultTask->SetEffMatrixEtaBins(effMatrixEtaBins);
  defaultTask->SetNSigmaTOF({-2.5, -3.0, -3.5}, {2.5, 2.5, 2.5});
  defaultTask->SetTOFMomCut(0.7);

  // Boolians which are by default === ON ===
  defaultTask->SetRunOnGrid(kTRUE);
  defaultTask->SetIsMCtrue(kFALSE);
  defaultTask->SetIncludeITScuts(kTRUE);
  defaultTask->SetFillArmPodTree(kFALSE);
  defaultTask->SetUsePtCut(1);
  defaultTask->SetMCTrackOriginType(1);   // 0:full scan, 1: prim
  defaultTask->SetRapidityType(0);      // 0:pseudorapidity, 1: rapidity

  // Extra Boolians which are by default === OFF ===
  defaultTask->SetDeDxCheck(kFALSE);
  defaultTask->SetFillEffMatrix(kFALSE);
  defaultTask->SetFillDebug(kFALSE);
  defaultTask->SetFillTracksMCgen(kFALSE);
  defaultTask->SetRunFastSimulation(kFALSE);
  defaultTask->SetFillHigherMomentsMCclosure(kFALSE);
  defaultTask->SetIncludeTOF(kFALSE);
  defaultTask->SetUseCouts(kFALSE);
  defaultTask->SetWeakAndMaterial(kFALSE);
  defaultTask->SetFillEventInfo(kFALSE);
  defaultTask->SetFillTreeMC(kFALSE);
  defaultTask->SetFillNudynFastGen(kFALSE);
  defaultTask->SetDefaultEventCuts(kFALSE);
  defaultTask->SetFillJetsBG(0);


  // Setters for the systematic uncertainty checks
  defaultTask->SetSystCentEstimator(0);

}
// ____________________________________________________________________________________________
TTree *GetLookUpTable(Bool_t runOnGrid, Int_t index)
{

  std::cout << " Info::marsland: Copy LookUp table from alien " << std::endl;
  // Define the lookup table file name
  TString fileName  ="";
  if (index==1) fileName="MomentsTree_AccCan_HIJING.root";
  if (index==2) fileName="MomentsTree_AccCan_LHC13f3a.root";
  if (index==3) fileName="MomentsTree_AccCan_LHC13f3b.root";
  if (index==4) fileName="MomentsTree_AccCan_LHC13f3c.root";
  if (index==5) fileName="MomentsTree_AccCan_HIJING_FullAcc.root";
  if (index==6) fileName="MomentsTree_AccCan_LHC13f3a_FullAcc.root";
  if (index==7) fileName="MomentsTree_AccCan_LHC13f3b_FullAcc.root";
  if (index==8) fileName="MomentsTree_AccCan_LHC13f3c_FullAcc.root";
  if (index==9) fileName="LookupTable_HighOrderCumulants_10h_HIJING.root";
  if (index==10) fileName="LookupTable_res_0_cent_HighOrderCumulants_10h_HIJING.root";
  if (index==11) fileName="LookupTable_res_0_centimp_HighOrderCumulants_10h_HIJING.root";

  //
  TTree *tree=NULL;
  TFile *fInputLookUp=NULL;
  TString lookUpDir ="";
  TString lookUpPath="";
  //
  // connect to alien for the lookup table
  std::cout << " Info::marsland: Connecting to GRID for the lookup table " << std::endl;

  if (runOnGrid){
    TGrid * alien = TGrid::Connect("alien://",0,0,"t");
    lookUpDir = "alien:///alice/cern.ch/user/m/marsland/PWGCF/EBYE/IdentityMethodEbyeFluctuations/macros";
    gSystem->Exec(Form("alien_cp %s/%s .",lookUpDir.Data(),fileName.Data()));
    lookUpPath = Form("%s/%s",gSystem->pwd(),fileName.Data());
  } else {
    lookUpDir = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/LookUpTables";
    lookUpPath = Form("%s/%s",lookUpDir.Data(),fileName.Data());
  }

  //
  // retrieve the ttree
  std::cout << " Info::marsland: LookUp table used is = " << lookUpPath << std::endl;
  fInputLookUp = new TFile(lookUpPath);
  tree = (TTree*)fInputLookUp->Get("mcMoments");
  //
  // return the lookup tree for further processing
  if(tree) return tree;
  else { std::cout << " Error::marsland: There is no lookUp table" << std::endl; return 0;}

}

vector<THnF*> GetEffMatrixObjects(Bool_t runOnGrid, TString filename) {

  TString lookUpPath = "";
  if (gSystem->AccessPathName(filename)){
    if (runOnGrid) {
      TGrid* alien = TGrid::Connect("alien://",0,0,"t");
      TString lookUpDir = "alien:///alice/cern.ch/user/i/ifokin/ebye/lookuptables/";
      printf("%s/%s\n",lookUpDir.Data(),filename.Data());
      gSystem->Exec(Form("alien_cp %s/%s file://.",lookUpDir.Data(),filename.Data()));
      lookUpPath = Form("%s/%s",gSystem->pwd(),filename.Data());
    } else {
      TString lookUpDir = "/home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/code";
      lookUpPath = Form("%s/%s",lookUpDir.Data(),filename.Data());
    }
  } else {
    lookUpPath = filename;
  }

  printf(" Info::ilya: Efficiency histograms used from %s\n", lookUpPath.Data());

  TFile* inputFile = new TFile(lookUpPath, "read");
  TList* cleanHists = (TList*) inputFile->Get("cleanHists");

  THnF* effMatrixPosRec  = (THnF*) cleanHists->FindObject("fHistPosEffMatrixScanRec");
  THnF* effMatrixNegRec  = (THnF*) cleanHists->FindObject("fHistNegEffMatrixScanRec");
  THnF* effMatrixPosGen  = (THnF*) cleanHists->FindObject("fHistPosEffMatrixScanGen");
  THnF* effMatrixNegGen  = (THnF*) cleanHists->FindObject("fHistNegEffMatrixScanGen");
  //
  return {effMatrixPosGen, effMatrixNegGen, effMatrixPosRec, effMatrixNegRec};
}

vector<TH2F*> GetQvecCorrObjects(Bool_t runOnGrid, TString filename) {

  TString lookUpPath = "";
  if (gSystem->AccessPathName(filename)){
    if (runOnGrid) {
      TGrid* alien = TGrid::Connect("alien://",0,0,"t");
      TString lookUpDir = "alien:///alice/cern.ch/user/i/ifokin/ebye/lookuptables/";
      printf("%s/%s\n",lookUpDir.Data(),filename.Data());
      gSystem->Exec(Form("alien_cp %s/%s file://.",lookUpDir.Data(),filename.Data()));
      lookUpPath = Form("%s/%s",gSystem->pwd(),filename.Data());
    } else {
      TString lookUpDir = "/home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/code";
      lookUpPath = Form("%s/%s",lookUpDir.Data(),filename.Data());
    }
  } else {
    lookUpPath = filename;
  }

  printf(" Info::ilya: Qvector correction histograms used from %s\n", lookUpPath.Data());

  TFile* inputFile = new TFile(lookUpPath, "read");
  TList* cleanHists = (TList*) inputFile->Get("cleanHists");

  TH2F* hist_EP_2_Qx_Qy_pos  = (TH2F*) cleanHists->FindObject("fHist_EP_2_Qx_Qy_pos");
  TH2F* hist_EP_2_Qx_Qy_neg  = (TH2F*) cleanHists->FindObject("fHist_EP_2_Qx_Qy_neg");
  TH2F* hist_EP_3_Qx_Qy_pos  = (TH2F*) cleanHists->FindObject("fHist_EP_3_Qx_Qy_pos");
  TH2F* hist_EP_3_Qx_Qy_neg  = (TH2F*) cleanHists->FindObject("fHist_EP_3_Qx_Qy_neg");

  //
  return {hist_EP_2_Qx_Qy_pos, hist_EP_2_Qx_Qy_neg, hist_EP_3_Qx_Qy_pos, hist_EP_3_Qx_Qy_neg};
}
