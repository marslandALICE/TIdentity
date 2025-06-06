#ifndef ALIANALYSISEBYERATIOS_H
#define ALIANALYSISEBYERATIOS_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// This analysis extracts pT-spectra of charged kaons, protons, and pions.  //
// It is based on particles identifation via the dE/dx signal of the TPC.   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TH1;
class THn;
class TH1F;
class TH2F;
class TH3D;
class TList;
class TTree;
class TObjArray;
class AliESDEvent;
class AliESDpid;
class AliESDtools;
class AliESDtrack;
class AliAODTrack;
class AliESDtrackCuts;
class AliHeader;
class AliMCParticle;
class AliPIDCombined;
class AliPIDResponse;
class AliStack;


#include "AliAnalysisTaskSE.h"
#include "AliPIDCombined.h"
#include "AliTPCdEdxInfo.h"
#include "AliESDv0KineCuts.h"
#include "AliESDv0.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "THn.h"
#include "TVectorF.h"
#include "TCutG.h"
#include "TTreeStream.h"
#include "AliESDv0Cuts.h"
#include "AliEventCuts.h"
#include "TRandom3.h"

#include <vector>
#include <utility>
#include <algorithm>

// class AliAnalysisTaskPIDetaTreeElectrons : public AliAnalysisTaskPIDV0base {
class AliAnalysisTaskTIdentityPID : public AliAnalysisTaskSE {
public:

  AliEventCuts fEventCuts;     /// Event cuts

  // ---------------------------------------------------------------------------------
  //                           Constructor and Destructor
  // ---------------------------------------------------------------------------------

  AliAnalysisTaskTIdentityPID(const char *name);
  AliAnalysisTaskTIdentityPID();
  virtual ~AliAnalysisTaskTIdentityPID();

  enum momentType {kPi=0,kKa=1,kPr=2,kPiPi=3,kKaKa=4,kPrPr=5,kPiKa=6,kPiPr=7,kKaPr=8,kLa=9,kLaLa=10,kCh=11,kChCh=12,kBa=13,kBaBa=14};

  enum momentTypeUnlike {
    kPiPosPiNeg=0,
    kKaPosKaNeg=1,
    kPrPosPrNeg=2,
    kPiPosKaNeg=3,
    kPiPosPrNeg=4,
    kKaPosPiNeg=5,
    kKaPosPrNeg=6,
    kPrPosPiNeg=7,
    kPrPosKaNeg=8,
    kLaPosLaNeg=9,
    kChPosChNeg=10,
    kBaPosBaNeg=11,
  };

  enum trackCutBit {
    kNCrossedRowsTPC70=0,
    kNCrossedRowsTPC80=1,
    kNCrossedRowsTPC90=2,
    kMaxChi2PerClusterTPCSmall=3,
    kMaxChi2PerClusterTPC=4,
    kMaxChi2PerClusterTPCLarge=5,
    kMaxDCAToVertexXYPtDep=6,
    kMaxDCAToVertexXYPtDepLarge=7,
    kVertexZSmall=8,
    kVertexZ=9,
    kEventVertexZ=10,
    kEventVertexZLarge=11,
    kTPCSignalNSmall=12,
    kTPCSignalN=13,
    kTPCSignalNLarge=14,
    kCleanPrTOF=15,
    kCleanKaTOF=16,
    kCleanKaTOFTRD=17,
    kITSPixelCut=18,
    kCleanPiTOF=19,
    kCleanDeTOF=20,
    kPileup=21,
    kPileupLoose=22,
    kSharedCls=23,
    kSharedClsLoose=24,
    kFindableCls=25,
    kFindableClsTight=26,
    kFindableClsLoose=27,
    kBFieldPos=28,
    kBFieldNeg=29,
    kNSigmaTOFLoose=30,
    kNSigmaTOFLoose2=31
  };

  enum cutSettings {
    kCutReference=0,
    kCutCrossedRowsTPC70=1,
    kCutCrossedRowsTPC90=2,
    kCutMaxChi2PerClusterTPCSmall=3,
    kCutMaxChi2PerClusterTPCLarge=4,
    kCutMaxDCAToVertexXYPtDepLarge=5,
    kCutVertexZSmall=6,
    kCutEventVertexZLarge=7,
    kCutSharedCls=8,
    kCutFindableClsTight=9,
    kCutFindableClsLoose=10,
    kCutPileupLoose=11,
    kCutBFieldPos=12,
    kCutBFieldNeg=13,
    kCutTPCSignalNSmall=14,
    kCutTPCSignalNLarge=15,
    kCutITSPixel=16,
    kCutNSigmaTOFLoose=17,
    kCutNSigmaTOFLoose2=18,
    kCutSettingsCount = 19
  };

  enum centEst {
    kV0M=0,
    kCL0=1,
    kCL1=2,
  };

  enum kPDGpart{
    kPDGel=11,
    kPDGpi=211,
    kPDGka=321,
    kPDGpr=2212,
    kPDGde=1000010020,
    kPDGmu=13,
    kPDGla=3122,
    kPDGxi=3312,
    kPDGd0=421,
    kPDGph=333,
  };

  enum kNetMoments{
    kA=0,
    kB=1,
    kAA=2,
    kBB=3,
    kAB=4,
    kAAA=5,
    kBBB=6,
    kAAB=7,
    kBBA=8,
    kABBB=9,
    kAABB=10,
    kAAAB=11,
    kAAAA=12,
    kBBBB=13,
    kAmB=14
  };

  /*
  kV0M=0,           // Centrality from V0A+V0C
  kCL0=1,           // Centrality from Clusters in layer 0
  kCL1=2,           // Centrality from Clusters in layer 1
  kTRK=3,           // Centrality from tracks
  kTKL=4,           // Centrality from tracklets
  kV0MvsFMD=5,      // Centrality from V0 vs FMD
  kTKLvsV0M=6,      // Centrality from tracklets vs V0
  kZEMvsZDC=7,      // Centrality from ZEM vs ZDC
  kV0A=8,           // Centrality from V0A
  kV0C=9,           // Centrality from V0C
  kZNA=10,          // Centrality from ZNA
  kZNC=11,          // Centrality from ZNC
  kZPA=12,          // Centrality from ZPA
  kZPC=13,          // Centrality from ZPC
  kCND=14,          // Centrality from tracks (candle condition)
  kFMD=15,          // Centrality from FMD
  kNPA=16,          // Centrality from Npart (MC)
  kV0A0=17,         // Centrality from V0A0
  kV0A123=18,       // Centrality from V0A123
  kV0A23=19,        // Centrality from V0A23
  kV0C01=20,        // Centrality from V0C01
  kV0S=21,          // Centrality from V0S
  kV0MEq=22,        // Centrality from V0A+V0C equalized channel
  kV0AEq=23,        // Centrality from V0A equalized channel
  kV0CEq=24,        // Centrality from V0C equalized channel
  kSPDClusters=25,  // Centrality from SPD Clusters
  kSPDTracklets=26, // Centrality from SPD Tracklets
  */
  //
  // ---------------------------------------------------------------------------------
  //                                    Methods
  // ---------------------------------------------------------------------------------
  //
  virtual void   UserCreateOutputObjects();            // create output objects
  virtual void   UserExec(Option_t *option);           // run over event-by-event and fill output objects
  virtual void   Terminate(Option_t *);                // run only once and terminate

  // ---------------------------------------------------------------------------------
  //                                    Settings
  // ---------------------------------------------------------------------------------

  void   Initialize();
  void   PrintNumInBinary(UInt_t num);
  void   SetESDtrackCuts(AliESDtrackCuts * trackCuts)                 {fESDtrackCuts        = trackCuts;};
  void   SetIsMCtrue(Bool_t isMCdata = kTRUE)                         {fMCtrue              = isMCdata;};
  void   SetYear(const Int_t ifYear = 0)                              {fYear                = ifYear;}
  void   SetPeriodName(const TString ifPeriodName = "")               {fPeriodName          = ifPeriodName;}
  void   SetPassIndex(const Int_t ifPassIndex = 0)                    {fPassIndex           = ifPassIndex;}
  // Some boolian settings
  void   SetRunOnGrid(const Bool_t ifRunOnGrid = kTRUE)               {fRunOnGrid           = ifRunOnGrid;}
  void   SetIncludeITScuts(const Bool_t ifITSCuts = kTRUE)            {fIncludeITS          = ifITSCuts;}
  void   SetFillArmPodTree(const Bool_t ifArmpodTree = kTRUE)         {fFillArmPodTree      = ifArmpodTree;}
  void   SetDeDxCheck(const Bool_t ifDeDxCheck = kFALSE)              {fDEdxCheck           = ifDeDxCheck;}
  void   SetFillEffMatrix(const Bool_t ifFillEffMatrix = kFALSE)      {fFillEffMatrix       = ifFillEffMatrix;}
  void   SetFillDebug(const Bool_t ifFillDebug = kFALSE)              {fFillDebug           = ifFillDebug;}
  void   SetFillTracksMCgen(const Bool_t ifFillTracksMCgen = kFALSE)  {fFillTracksMCgen     = ifFillTracksMCgen;}
  void   SetFillEffLookUpTable(const Bool_t ifEffLookUpTable = kFALSE){fFillEffLookUpTable  = ifEffLookUpTable;}
  void   SetFillHigherMomentsMCclosure(const Bool_t ifHigherMomentsMCclosure = kFALSE){fFillHigherMomentsMCclosure  = ifHigherMomentsMCclosure;}
  void   SetRunFastSimulation(const Bool_t ifFastSimul = kFALSE)      {fRunFastSimulation   = ifFastSimul;}
  void   SetRunFastHighMomentCal(const Bool_t ifFastHighMom = kFALSE) {fRunFastHighMomentCal= ifFastHighMom;}
  void   SetRunCutBasedMethod(const Bool_t ifCutBased = kFALSE)       {fRunCutBasedMethod   = ifCutBased;}
  void   SetFillDistributions(const Bool_t ifGenDistributions = kFALSE) {fFillDistributions= ifGenDistributions;}
  void   SetFillTreeMC(const Bool_t ifTreeMC = kFALSE)                {fFillTreeMC= ifTreeMC;}

  void   SetDefaultEventCuts(const Bool_t ifDefaultEventCuts = kFALSE){fDefaultEventCuts= ifDefaultEventCuts;}
  void   SetFillNudynFastGen(const Bool_t ifNudynFastGen = kFALSE)    {fFillNudynFastGen= ifNudynFastGen;}
  void   SetCorrectForMissCl(const Int_t ifCorrectForMissCl = kFALSE) {fCorrectForMissCl= ifCorrectForMissCl;}
  void   SetUsePtCut(const Int_t ifUsePtCut = 1)                      {fUsePtCut            = ifUsePtCut;}
  void   SetMCTrackOriginType(const Int_t ifTrackOriginOnlyPrimary = 0) {fTrackOriginOnlyPrimary     = ifTrackOriginOnlyPrimary;}
  void   SetRapidityType(const Int_t ifRapidityType = 0)              {fRapidityType        = ifRapidityType;}
  void   SetSisterCheck(const Int_t ifSisterCheck = 0)                {fSisterCheck         = ifSisterCheck;}
  void   SetIncludeTOF(const Bool_t ifIncludeTOF = kFALSE)            {fIncludeTOF          = ifIncludeTOF;}
  void   SetUseCouts(const Bool_t ifUseCouts = kFALSE)                {fUseCouts            = ifUseCouts;}
  void   SetWeakAndMaterial(const Bool_t ifWeakAndMaterial = kFALSE)  {fWeakAndMaterial     = ifWeakAndMaterial;}
  void   SetFillEventInfo(const Bool_t ifEventInfo = kFALSE)          {fEventInfo           = ifEventInfo;}
  void   SetDownscalingFactor(const Float_t nDownscalingFactor = 0)   {fDownscalingFactor   = nDownscalingFactor;}
  void   SetV0InvMassHists(const Bool_t ifV0InvMassHists = kFALSE)    {fV0InvMassHists      = ifV0InvMassHists;}
  void   SetRunNumberForExpecteds(const Int_t ifRunNumberForExpecteds = 0)    {fRunNumberForExpecteds = ifRunNumberForExpecteds;}
  void   SetFillJetsBG(const Int_t ifFillJetsBG = kTRUE)              {fFillJetsBG          = ifFillJetsBG;}
  void   SetTaskSelection(const Int_t ifTaskSelection = kTRUE)        {fTaskSelection       = ifTaskSelection;}
  void   SetFillResonances(const Bool_t ifFillResonances = kFALSE)    {fFillResonances      = ifFillResonances;}
  void   SetFillQvectorHists(const Bool_t ifFillQvectorHists = kFALSE){fFillQvectorHists    = ifFillQvectorHists;}
  void   SetApplyQVectorCorr(const Bool_t ifApplyQVectorCorr = kFALSE){fApplyQVectorCorr    = ifApplyQVectorCorr;}
  void   SetDownsampleTrees(const Bool_t ifDownsampleTrees = kFALSE)  {fDownsampleTrees     = ifDownsampleTrees;}
  void   SetUseAODsForMC(const Bool_t ifUseAODsForMC = kFALSE)  {fUseAODsForMC     = ifUseAODsForMC;}

  void   SetSettings(const std::vector<Int_t> ifSystSettings) {
    fSystSettings = ifSystSettings;
    fNSettings = fSystSettings.size();
  }
  void   SetNSettings(const Int_t nSettings = 22) {
    std::vector<Int_t> tempSettings(nSettings);
    for (Int_t i = 0; i < nSettings; i++)
      tempSettings[i] = i;
    SetSettings(tempSettings);
  }

  //
  Bool_t GetRunOnGrid() const { return fRunOnGrid; }

  // Setters for the systematic uncertainty checks
  void   SetSystCentEstimator(const Int_t systCentEstimator = 0)  {fSystCentEstimatetor = systCentEstimator;}

  // Setters for the eta momentum dEdx and centrality bins
  void   SetSampleDeDxUpperEdge(const Float_t dEdxCleanUp = 200.) {fDEdxCleanUp         = dEdxCleanUp;}
  void   SetDeDxBinWidth(const Float_t dEdxBinWidth = 2.5)        {fDEdxBinWidth        = dEdxBinWidth;}
  void   SetDeDxLowerEdge(const Float_t dEdxLowerEdge = 20.)      {fDEdxDown            = dEdxLowerEdge;}
  void   SetDeDxUpperEdge(const Float_t dEdxUpperEdge = 1020.)    {fDEdxUp              = dEdxUpperEdge;}

  void   SetEtaLowerEdge(const Float_t etaLowerEdge = -0.8)       {fEtaMin              = etaLowerEdge;}
  void   SetEtaUpperEdge(const Float_t etaUpperEdge = 0.8)        {fEtaMax              = etaUpperEdge;}
  void   SetNEtabins(const Int_t nEtaBins = 20)                   {fNEtaBins            = nEtaBins;}
  void   SetMomLowerEdge(const Float_t momLowerEdge = 0.)         {fMomDown             = momLowerEdge;}
  void   SetMomUpperEdge(const Float_t momUpperEdge = 12.)        {fMomUp               = momUpperEdge;}
  void   SetNMomBins(const Int_t nMombins = 600)                  {fNMomBins            = nMombins;}
  void   SetCollisionType(Int_t collisionType = 0)          {fCollisionType       = collisionType;}

  void   SetEffMatrixMomBins(const std::vector<Double_t> nEffMatrixMomBins) {fEffMatrixMomBins = nEffMatrixMomBins;}
  void   SetEffMatrixCentBins(const std::vector<Double_t> nEffMatrixCentBins) {fEffMatrixCentBins = nEffMatrixCentBins;}
  void   SetEffMatrixEtaBins(const std::vector<Double_t> nEffMatrixEtaBins) {fEffMatrixEtaBins = nEffMatrixEtaBins;}

  void   SetNSigmaTPC(const std::vector<Double_t>& vecNSigmaTPC)  { fNSigmaTPC = vecNSigmaTPC; }
  void   SetNSigmaTPC(const Double_t nSigmaTPC)                   {
    std::vector<Double_t> vecNSigmaTPC(3, nSigmaTPC);
    SetNSigmaTPC(vecNSigmaTPC);
  }
  void   SetNSigmaTOF(const std::vector<Double_t>& vecNSigmaTOFDown, const std::vector<Double_t>& vecNSigmaTOFUp) {
    fNSigmaTOFDown = vecNSigmaTOFDown;
    fNSigmaTOFUp = vecNSigmaTOFUp;
  }
  void   SetTOFMomCut(const Double_t tofMomCut)                   { fTOFMomCut = tofMomCut; }

  std::vector<Double_t> GetNSigmas(const Int_t setting) {
    std::vector<Double_t> ret(4);
    switch (setting)
    {
      case kCutNSigmaTOFLoose:
      ret = {-fNSigmaTPC[1], fNSigmaTPC[1], fNSigmaTOFDown[1], fNSigmaTOFUp[1]};
      break;
      case kCutNSigmaTOFLoose2:
      ret = {-fNSigmaTPC[2], fNSigmaTPC[2], fNSigmaTOFDown[2], fNSigmaTOFUp[2]};
      break;
      default:
      ret = {-fNSigmaTPC[0], fNSigmaTPC[0], fNSigmaTOFDown[0], fNSigmaTOFUp[0]};
      break;
    }
    return ret;
  }


  // Set the binning of centrality
  void SetCentralityBinning(const Int_t tmpCentbins, Float_t tmpfxCentBins[])
  {
    // Create the histograms to be used in the binning of eta, cent and momentum
    std::cout << " Info::marsland: !!!!!! Centrality binning is being set !!!!!!! " << std::endl;
    fHistCent =  new TH1F("fHistCent","Centrality Bins",tmpCentbins-1    ,tmpfxCentBins );
    fHistPhi  =  new TH1F("fHistPhi" ,"Phi Bins"       ,36               ,-TMath::Pi(), TMath::Pi());
    // ==========================================
    // prepare real data centrality bins
    fNCentbinsData = tmpCentbins;
    fNCentBinsMC   = tmpCentbins-1;
    fxCentBins.resize(fNCentbinsData);
    for (Int_t i=0; i<fNCentbinsData; i++) fxCentBins[i] =  tmpfxCentBins[i];
    fcentDownArr.resize(fNCentBinsMC);
    fcentUpArr.resize(fNCentBinsMC);
    for (Int_t i=0; i<fNCentbinsData-1; i++) fcentDownArr[i] =  tmpfxCentBins[i];
    for (Int_t i=1; i<fNCentbinsData; i++)   fcentUpArr[i-1] =  tmpfxCentBins[i];
  }

  void SetMCEtaScanArray(const Int_t tmpEtaBinsMC, Float_t tmpetaDownArr[], Float_t tmpetaUpArr[])
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetMCEtaScanArray is being set !!!!!!! " << std::endl;
    fNEtaWinBinsMC = tmpEtaBinsMC;
    fetaDownArr.resize(fNEtaWinBinsMC);
    fetaUpArr.resize(fNEtaWinBinsMC);
    for (Int_t i=0; i<fNEtaWinBinsMC; i++) {
      fetaDownArr[i] =  tmpetaDownArr[i];
      fetaUpArr[i]   =  tmpetaUpArr[i];
    }
  }

  void SetMCResonanceArray(const Int_t tmpNRes, TString tmpResArr[])
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetMCResonanceArray is being set !!!!!!! " << std::endl;
    fNResBins = tmpNRes;
    fResonances.resize(fNResBins);
    for (Int_t i=0; i<fNResBins; i++) fResonances[i] = tmpResArr[i];

  }

  void SetMCBaryonArray(const Int_t tmpNBar, Int_t tmpBarArr[])
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetMCBaryonArray is being set !!!!!!! " << std::endl;
    fNBarBins = tmpNBar;
    fBaryons.resize(fNBarBins);
    for (Int_t i=0; i<fNBarBins; i++) fBaryons[i] = tmpBarArr[i];
  }

  void SetMCMomScanArray(const Int_t tmpMomBinsMC, Float_t tmppDownArr[], Float_t tmppUpArr[])
  {
    // set MC momentum values to scan
    std::cout << " Info::marsland: !!!!!! SetMCMomScanArray is being set !!!!!!! " << std::endl;
    fNMomBinsMC = tmpMomBinsMC;
    fpDownArr.resize(fNMomBinsMC);
    fpUpArr.resize(fNMomBinsMC);
    for (Int_t i=0; i<fNMomBinsMC; i++) {
      fpDownArr[i] =  tmppDownArr[i];
      fpUpArr[i]   =  tmppUpArr[i];
    }
  }

  void SetLookUpTableFirstMoments(TTree *lookUpTree, Int_t partType, Float_t pArr[],Float_t centArr[],Float_t etaArr[],const Int_t tmpMomBinsMC, const Int_t tmpCentbins, const Int_t tmpEtaBinsMC)
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetLookUpTableFirstMoments is being set !!!!!!!   " << std::endl;
    //
    // fill arrays from lookup table
    TH1D *h=NULL, *h1=NULL;
    for (Int_t imom=0; imom<tmpMomBinsMC; imom++){
      for (Int_t icent=0; icent<tmpCentbins; icent++){
        for (Int_t ieta=0; ieta<tmpEtaBinsMC; ieta++){
          //
          // with resonances
          lookUpTree->Draw(Form("momentPos.fElements[%d]-momentNeg.fElements[%d]",partType,partType),Form("abs(etaUp-%f)<0.01&&abs(pDown-%f)<0.01&&abs(centDown-%f)<0.01",etaArr[ieta],pArr[imom],centArr[icent]),"goff");
          h= (TH1D*)lookUpTree->GetHistogram()->Clone(); h-> SetName("Res");
          if (partType==0)  fNetPiFirstMoments[0][imom][icent][ieta] = h->GetMean();
          if (partType==1)  fNetKaFirstMoments[0][imom][icent][ieta] = h->GetMean();
          if (partType==2)  fNetPrFirstMoments[0][imom][icent][ieta] = h->GetMean();
          if (partType==9)  fNetLaFirstMoments[0][imom][icent][ieta] = h->GetMean();
          if (partType==11) fNetChFirstMoments[0][imom][icent][ieta] = h->GetMean();
          delete h;
          //
          // without resonances
          lookUpTree->Draw(Form("noResmomentPos.fElements[%d]-noResmomentNeg.fElements[%d]",partType,partType),Form("abs(etaUp-%f)<0.01&&abs(pDown-%f)<0.01&&abs(centDown-%f)<0.01",etaArr[ieta],pArr[imom],centArr[icent]),"goff");
          h1= (TH1D*)lookUpTree->GetHistogram()->Clone(); h1-> SetName("noRes");
          if (partType==0)  fNetPiFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          if (partType==1)  fNetKaFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          if (partType==2)  fNetPrFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          if (partType==9)  fNetLaFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          if (partType==11) fNetChFirstMoments[1][imom][icent][ieta] = h1->GetMean();
          delete h1;

        }
      }
    }

  }

  void SetLookUpTableEfficiencyCorrection(TTree *lookUpTree, Int_t partType, Float_t pDownArr[], Float_t pUpArr[], Float_t centArr[],Float_t etaArr[],const Int_t tmpMomBinsMC, const Int_t tmpCentbins, const Int_t tmpEtaBinsMC)
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetLookUpTableEfficiencyCorrection is being set !!!!!!!   " << std::endl;
    //
    // fill arrays from lookup table
    TH1D *hGenNet=NULL,   *hRecNet=NULL;
    TH1D *hGenCross=NULL, *hRecCross=NULL;
    TH1D *hGenPos=NULL,   *hRecPos=NULL;
    TH1D *hGenNeg=NULL,   *hRecNeg=NULL;
    for (Int_t imom=0; imom<tmpMomBinsMC; imom++){
      for (Int_t ieta=0; ieta<tmpEtaBinsMC; ieta++){
        for (Int_t icent=0; icent<tmpCentbins; icent++){
          //
          // ----------------------------
          // NET Particle cumulants
          // ----------------------------
          //
          TString etaMomCentCut = Form("abs(etaUp-%f)<0.01 && abs(pDown-%f)<0.01 && abs(pUp-%f)<0.01 && abs(centDown-%f)<0.01",etaArr[ieta],pDownArr[imom],pUpArr[imom],centArr[icent]);
          // generated net particles
          TString genStr = Form("momentPosGen.fElements[%d]-momentNegGen.fElements[%d]",partType,partType);
          lookUpTree->Draw(genStr,etaMomCentCut,"goff");
          hGenNet = (TH1D*)lookUpTree->GetHistogram()->Clone(); hGenNet-> SetName("netProtonGen");
          if (partType==0)  fNetPiFirstMomentsGen[imom][icent][ieta] = hGenNet->GetMean();
          if (partType==1)  fNetKaFirstMomentsGen[imom][icent][ieta] = hGenNet->GetMean();
          if (partType==2)  fNetPrFirstMomentsGen[imom][icent][ieta] = hGenNet->GetMean();
          //
          // reconstructed net particles
          TString recStr = Form("momentPosRec.fElements[%d]-momentNegRec.fElements[%d]",partType,partType);
          lookUpTree->Draw(recStr,etaMomCentCut,"goff");
          hRecNet= (TH1D*)lookUpTree->GetHistogram()->Clone(); hRecNet-> SetName("netProtonRec");
          if (partType==0)  fNetPiFirstMomentsRec[imom][icent][ieta] = hRecNet->GetMean();
          if (partType==1)  fNetKaFirstMomentsRec[imom][icent][ieta] = hRecNet->GetMean();
          if (partType==2)  fNetPrFirstMomentsRec[imom][icent][ieta] = hRecNet->GetMean();
          //
          // ----------------------------
          // Cross cumulants
          // ----------------------------
          //
          // generated net particles
          TString momentCrossGenStr = Form("momentCrossGen.fElements[%d]",partType);
          lookUpTree->Draw(momentCrossGenStr,etaMomCentCut,"goff");
          hGenCross = (TH1D*)lookUpTree->GetHistogram()->Clone(); hGenCross-> SetName("crossProtonGen");
          if (partType==0)  fCrossPiFirstMomentsGen[imom][icent][ieta] = hGenCross->GetMean();
          if (partType==1)  fCrossKaFirstMomentsGen[imom][icent][ieta] = hGenCross->GetMean();
          if (partType==2)  fCrossPrFirstMomentsGen[imom][icent][ieta] = hGenCross->GetMean();
          //
          // reconstructed net particles
          TString momentCrossRecStr = Form("momentCrossRec.fElements[%d]",partType);
          lookUpTree->Draw(momentCrossRecStr,etaMomCentCut,"goff");
          hRecCross= (TH1D*)lookUpTree->GetHistogram()->Clone(); hRecCross-> SetName("crossProtonRec");
          if (partType==0)  fCrossPiFirstMomentsRec[imom][icent][ieta] = hRecCross->GetMean();
          if (partType==1)  fCrossKaFirstMomentsRec[imom][icent][ieta] = hRecCross->GetMean();
          if (partType==2)  fCrossPrFirstMomentsRec[imom][icent][ieta] = hRecCross->GetMean();
          //
          // ----------------------------
          // Single  Particle cumulants
          // ----------------------------
          //
          // generated particles
          TString momentPosGenStr = Form("momentPosGen.fElements[%d]",partType);
          lookUpTree->Draw(momentPosGenStr,etaMomCentCut,"goff");
          hGenPos = (TH1D*)lookUpTree->GetHistogram()->Clone(); hGenPos-> SetName("protonPosGen");
          if (partType==0)  fPiFirstMomentsGen[0][imom][icent][ieta] = hGenPos->GetMean();
          if (partType==1)  fKaFirstMomentsGen[0][imom][icent][ieta] = hGenPos->GetMean();
          if (partType==2)  fPrFirstMomentsGen[0][imom][icent][ieta] = hGenPos->GetMean();
          //
          TString momentNegGenStr = Form("momentNegGen.fElements[%d]",partType);
          lookUpTree->Draw(momentNegGenStr,etaMomCentCut,"goff");
          hGenNeg = (TH1D*)lookUpTree->GetHistogram()->Clone(); hGenNeg-> SetName("protonNegGen");
          if (partType==0)  fPiFirstMomentsGen[1][imom][icent][ieta] = hGenNeg->GetMean();
          if (partType==1)  fKaFirstMomentsGen[1][imom][icent][ieta] = hGenNeg->GetMean();
          if (partType==2)  fPrFirstMomentsGen[1][imom][icent][ieta] = hGenNeg->GetMean();
          //
          // reconstruced particles
          TString momentPosRecStr = Form("momentPosRec.fElements[%d]",partType);
          lookUpTree->Draw(momentPosRecStr,etaMomCentCut,"goff");
          hRecPos = (TH1D*)lookUpTree->GetHistogram()->Clone(); hRecPos-> SetName("protonPosRec");
          if (partType==0)  fPiFirstMomentsRec[0][imom][icent][ieta] = hRecPos->GetMean();
          if (partType==1)  fKaFirstMomentsRec[0][imom][icent][ieta] = hRecPos->GetMean();
          if (partType==2)  fPrFirstMomentsRec[0][imom][icent][ieta] = hRecPos->GetMean();
          //
          // reconstruced negative particles
          TString momentNegRecStr = Form("momentNegRec.fElements[%d]",partType);
          lookUpTree->Draw(momentNegRecStr,etaMomCentCut,"goff");
          hRecNeg = (TH1D*)lookUpTree->GetHistogram()->Clone(); hRecNeg-> SetName("protonNegRec");
          if (partType==0)  fPiFirstMomentsRec[1][imom][icent][ieta] = hRecNeg->GetMean();
          if (partType==1)  fKaFirstMomentsRec[1][imom][icent][ieta] = hRecNeg->GetMean();
          if (partType==2)  fPrFirstMomentsRec[1][imom][icent][ieta] = hRecNeg->GetMean();
          //
          // delete pointers to tmp histograms
          delete hGenNet;
          delete hRecNet;
          delete hGenPos;
          delete hGenNeg;
          delete hRecPos;
          delete hRecNeg;
          delete hGenCross;
          delete hRecCross;


          //

        }
      }
    }

  }

  void SetEffMatrixObjects(THnF* effMatrixGenPos, THnF* effMatrixGenNeg, THnF* effMatrixRecPos, THnF* effMatrixRecNeg)
  {
    // set MC eta values to scan
    std::cout << " Info::ilya: !!!!!! SetEffMatrixObjects is being set !!!!!!!   " << std::endl;
    //
    fEffMatrixGenPos = (THnF*) effMatrixGenPos->Clone();
    fEffMatrixGenNeg = (THnF*) effMatrixGenNeg->Clone();
    fEffMatrixRecPos = (THnF*) effMatrixRecPos->Clone();
    fEffMatrixRecNeg = (THnF*) effMatrixRecNeg->Clone();
  }

  void SetQvecCorrObjects(TH2F* histEP2QxQypos, TH2F* histEP2QxQyneg, TH2F* histEP3QxQypos, TH2F* histEP3QxQyneg)
  {
    // set MC eta values to scan
    std::cout << " Info::ilya: !!!!!! SetQvecCorrObjects is being set !!!!!!!   " << std::endl;
    //
    fHistEP2QxQypos = (TH2F*) histEP2QxQypos->Clone();
    fHistEP2QxQyneg = (TH2F*) histEP2QxQyneg->Clone();
    fHistEP3QxQypos = (TH2F*) histEP3QxQypos->Clone();
    fHistEP3QxQyneg = (TH2F*) histEP3QxQyneg->Clone();
  }

  void SetLookUpTable_MissCl(TClonesArray *lookUpArray)
  {
    // set MC eta values to scan
    std::cout << " Info::marsland: !!!!!! SetLookUpTable_MissCl is being set !!!!!!!   " << std::endl;
    //
    const Int_t nParticles=4;
    const Int_t nMultBins=16;
    for (Int_t ipart=0; ipart<nParticles; ipart++){
      for (Int_t icent=0; icent<nMultBins; icent++){
        TString objname = Form("hDiffMissCl_part_%d_cent_%d",ipart,icent);
        fH2MissCl[ipart][icent] = *((TH2F*)(lookUpArray->FindObject(objname))->Clone());
        std::cout << fH2MissCl[ipart][icent].GetName() << std::endl;
      }
    }

  }

  Bool_t IsFromPileup(Int_t label) {
    Bool_t isTPCPileup = kFALSE;
    Bool_t isITSPileup = kFALSE;

    if (fCollisionType == 0) {
      isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(label,fMCEvent);
      isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
    }

    return isTPCPileup || isITSPileup;
  }

private:

  AliAnalysisTaskTIdentityPID(const AliAnalysisTaskTIdentityPID&);
  AliAnalysisTaskTIdentityPID& operator=(const AliAnalysisTaskTIdentityPID&);

  // ---------------------------------------------------------------------------------
  //                                   Functions
  // ---------------------------------------------------------------------------------

  void FillTPCdEdxReal();                   // Main function to fill all info + TIden
  void CalculateMoments_CutBasedMethod();
  void FillMCFull_NetParticles();
  void FillTreeMC();
  void FillTreeMCAOD();
  void FastGen_NetParticles();              // Run over galice.root for Fastgen higher moments
  void MCclosureHigherMoments();   // Calculate higher moments for REC and GEN
  void FillEffMatrix();            // Prepare efficiency matrix
  void FillCleanSamples();                    // Fill Clean Pions
  void SelectCleanSamplesFromV0s(AliESDv0 *v0, AliESDtrack *track0, AliESDtrack *track1);
  void SetSpecialV0Cuts(AliESDv0KineCuts* cuts);
  void BinLogAxis(TH1 *h);
  void CalculateEventInfo();
  void DumpDownScaledTree();
  void GetExpecteds(AliESDtrack *track, Double_t closestPar[3]);
  void CreateEventInfoTree();
  //
  Int_t CountEmptyEvents(Int_t counterBin, Int_t setting);
  Int_t CacheTPCEventInformation();
  UInt_t SetCutBitsAndSomeTrackVariables(AliESDtrack *track, Int_t particleType);
  Bool_t CheckIfFromResonance(AliMCParticle *trackMCgen);
  Bool_t CheckIfBaryonAOD(AliAODMCParticle *trackMCgen);
  Bool_t CheckIfBaryon(AliMCParticle *trackMCgen);
  Bool_t CheckIfFromAnyResonance(AliMCParticle *trackMCgen, Float_t etaLow, Float_t etaUp, Float_t pDown, Float_t pUp);
  Bool_t ApplyDCAcutIfNoITSPixel(AliESDtrack *track);
  Bool_t ApplyDCAcutIfNoITSPixelAOD(AliAODTrack *track);
  Bool_t GetSystematicClassIndex(UInt_t cut,Int_t syst);
  Bool_t CheckPsiPair(const AliESDv0* v0);
  Double_t GetTrackEfficiency(const Int_t& part, const Double_t& ptot, const Double_t& eta, const Int_t& setting, const Int_t& sign, const Int_t& pidCutType=0);


  //
  // Jet Functions
  void FindJetsFJ();
  void FindJetsFJGen();
  Int_t MakeEventPlane(Int_t doEP_Psi2, Int_t doEP_Psi3, Int_t setting);
  void GetFlatenicityMC();
  void GetFlatenicity();
  Double_t ComputeSpherocity(Int_t setting);
  void FillEventInfoMCAOD();
  void FillEventInfoMC();
  Float_t RelativePhi(Float_t mphi, Float_t vphi);

  //
  // ---------------------------------------------------------------------------------
  //                                   Members
  // ---------------------------------------------------------------------------------
  //
  AliPIDResponse   * fPIDResponse;            //! PID response object
  AliESDEvent      * fESD;                    //! ESD object
  AliAODEvent      * fAOD;                    //! AOD object
  TList            * fListHist;               //! list for histograms
  AliESDtrackCuts  * fESDtrackCuts;           //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit96;     //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit128;    //! basic cut variables
  AliESDtrackCuts  * fESDtrackCuts_Bit768;    //! basic cut variables
  AliESDtrackCuts  * fESDtrackCutsLoose;      //! basic cut variables for debugging
  AliESDv0Cuts     * fESDtrackCutsV0;         //! basic cut variables for V0
  AliESDtrackCuts  * fESDtrackCutsCleanSamp;  //! basic cut variables for clean pion and electron form V0s
  AliPIDCombined   * fPIDCombined;            //! combined PID object
  AliTPCdEdxInfo   * fTPCdEdxInfo;            //! detailed dEdx info
  AliStack         * fMCStack;                //! stack object to get Mc info
  AliESDv0KineCuts * fV0OpenCuts;             // v0 strong filter for tagged V0s
  AliESDv0KineCuts * fV0StrongCuts;           // v0 strong filter for tagged V0s
  AliAnalysisCuts  * fK0sPionCuts;            // filter for pions from K0s
  AliAnalysisCuts  * fLambdaProtonCuts;       // filter for protons from Lambda
  AliAnalysisCuts  * fLambdaPionCuts;         // filter for pions from Lambda
  AliAnalysisCuts  * fGammaElectronCuts;      // filter for electrons from gamma conversions
  const AliESDVertex * fVertex;               // primary vertex
  AliESDtools      * fESDtool;                 // tools to calculate derived variables from the ESD

  TTree            * fArmPodTree;             // Tree for clean pion and proton selection
  TTreeSRedirector * fTreeSRedirector;        /// temp tree to dump output
  TTree            * fTreeMomentsMCfull;      // tree for reconstructed moments
  TTree            * fTreeTracksMCgen;        // tree for reconstructed moments
  TTree            * fTreeDebug3;             // tree for dnch/deta calculation
  TTree            * fTreeTracksMCrec;        // tree for mc samples
  TTree            * fTreeDebug4;             // tree to check dEdx performance for a small data sample
  TTree            * fTreeTracks;             // tree to save all variables for control plots
  TTree            * fTreeDebug;              // tree for debugging
  TTree            * fTreeResonances;         // tree with full acceptance filled with MC
  TTree            * fTreeMomentsMCgen;       // tree with higher moment calculations
  TTree            * fTreeEvents;
  TTree            * fTreeEventsMC;
  TTree            * fTreeDScaled;
  TTree            * fTreeDebug2;
  TTree            * fTreeExpecteds;          // tree with expected dE/dx
  TTree            * fTreeCutBased;           // tree with moments from cut based method
  TTree            * fTreejetsFJ;             // tree for fastjet signal jets
  TTree            * fTreejetsFJBG;           // tree for fastjet background jets
  TTree            * fTreejetsFJconst;          // tree for fastjet signal jet constituents
  TTree            * fTreejetsFJGen;
  TTree            * fTreejetsFJBGGen;
  TTree            * fTreejetsFJconstGen;
  TRandom3         fRandom;


  TString            fPeriodName;
  Int_t              fYear;
  Int_t              fPassIndex;
  UInt_t             fPileUpBit;
  TH1F             * fHistCent;               // helper histogram for TIdentity tree
  TH1F             * fHistPhi;
  TH2F             * fHistRapDistFullAccPr;
  TH2F             * fHistRapDistFullAccAPr;
  TH1F             * fHistNhard;
  TH1F             * fHistNproj;
  TH1F             * fHistNtarget;
  TH1F             * fHistNN;
  TH1F             * fHistNNW;
  TH1F             * fHistNWN;
  TH1F             * fHistNWNW;
  TH1F             * fHistNpionsInTPC;
  TH1F             * fHistNkaonsInTPC;
  TH1F             * fHistNprotonsInTPC;
  TH1F             * fHistNchargedInTPC;
  TH1F             * fHistNallInTPC;
  TH1F             * fHistNpionsInV0M;
  TH1F             * fHistNkaonsInV0M;
  TH1F             * fHistNprotonsInV0M;
  TH1F             * fHistNchargedInV0M;
  TH1F             * fHistNallInV0M;

  TH2F             * fHist_EP_2_Qx_Qy_pos;
  TH2F             * fHist_EP_2_Qx_Qy_neg;
  TH1F             * fHist_EP_2_Psi_pos;
  TH1F             * fHist_EP_2_Psi_neg;
  TH1F             * fHist_EP_2_Psi;
  TH2F             * fHist_EP_3_Qx_Qy_pos;
  TH2F             * fHist_EP_3_Qx_Qy_neg;
  TH1F             * fHist_EP_3_Psi_pos;
  TH1F             * fHist_EP_3_Psi_neg;
  TH1F             * fHist_EP_3_Psi;


  TH1F             * fHistInvK0s;             // helper histogram for TIdentity tree
  TH1F             * fHistInvLambda;          // helper histogram for TIdentity tree
  TH1F             * fHistInvAntiLambda;      // helper histogram for TIdentity tree
  TH1F             * fHistInvPhoton;          // helper histogram for TIdentity tree
  //
  TH1F             * fHistPhiTPCcounterA;     // helper histogram for TIdentity tree
  TH1F             * fHistPhiTPCcounterC;     // helper histogram for TIdentity tree
  TH1F             * fHistPhiTPCcounterAITS;  // helper histogram for TIdentity tree
  TH1F             * fHistPhiTPCcounterCITS;  // helper histogram for TIdentity tree
  TH1F             * fHistPhiITScounterA;     // helper histogram for TIdentity tree
  TH1F             * fHistPhiITScounterC;     // helper histogram for TIdentity tree

  TH2F              fH2MissCl[4][16];          // histogram which hold all PIDresponse info

  TString           fChunkName;

  UInt_t            fTrackCutBits;           // integer which hold all cut variations as bits
  Int_t             fSystClass;
  Double_t          fEtaMin;
  Double_t          fEtaMax;
  Int_t             fNEtaBins;
  Float_t           fDownscalingFactor;     // when only a fDownscalingFactor is enough
  Int_t             fEventIDinFile;         // event id in file
  Int_t             fChunkIDinJob;          // event id in file

  Bool_t            fRunOnGrid;              // flag if real data or MC is processed
  Bool_t            fMCtrue;                 // flag if real data or MC is processed
  Bool_t            fEventInfo;              // flag if event info and downscaled track tree is filled
  Bool_t            fWeakAndMaterial;        // flag for the Weak and Material analysis
  Bool_t            fFillEffMatrix;          // flag for efficiency matrix filling
  Bool_t            fFillDebug;              // flag for efficiency matrix filling
  Bool_t            fDEdxCheck;              // flag to check only the dEdx performance
  Bool_t            fIncludeITS;             // decide whether to use ITS or not
  Bool_t            fFillTracksMCgen;          //
  Bool_t            fFillEffLookUpTable;     //
  Bool_t            fFillHigherMomentsMCclosure;
  Bool_t            fFillArmPodTree;         // switch whether to fill clean sample tree
  Bool_t            fRunFastSimulation;      // when running over galice.root do not fill other objects
  Bool_t            fRunFastHighMomentCal;   // when running over galice.root do not fill other objects
  Bool_t            fRunCutBasedMethod;      // moments from cut based method as cross check
  Bool_t            fFillDistributions;   // when running over galice.root do not fill other objects
  Bool_t            fFillTreeMC;
  Bool_t            fDefaultEventCuts;
  Bool_t            fFillNudynFastGen;
  Bool_t            fFillResonances;
  Int_t             fCollisionType;
  Int_t             fCorrectForMissCl;       // 0; defaults crows, 1; ncls used wo correction, 2; ncls used with correction
  Int_t             fUsePtCut;
  Int_t             fTrackOriginOnlyPrimary;
  Int_t             fRapidityType;
  Int_t             fSisterCheck;           // 0: reject the mother anyways, 1: if both girls are in acceptance rejet mother

  Bool_t            fIncludeTOF;             // Include TOF information to investigate the efficiency loss effects on observable
  Bool_t            fUseCouts;               // for debugging
  Bool_t            fV0InvMassHists;         // V0 invariant mass for QA
  Int_t             fRunNumberForExpecteds;  // Run number in which to fill the expecteds tree
  Bool_t            fFillExpecteds;
  Bool_t            fFillQvectorHists;
  Bool_t            fApplyQVectorCorr;
  Bool_t            fDefaultCuts;
  Bool_t            fDownsampleTrees;        // downsample the trees to save disk space
  Bool_t            fUseAODsForMC;

  Int_t             fNSettings;
  std::vector<Int_t> fSystSettings;
  Int_t             fNMomBins;               // number of mombins --> for 20MeV slice 150 and 10MeV 300
  Float_t           fMomDown;                // bottom limit for the momentum range (default 0.2)
  Float_t           fMomUp;                  // uppper limit for the momentum range (default 3.2)
  Float_t           fDEdxBinWidth;           // bin width for the dEdx histograms (default 2.5)
  Float_t           fDEdxUp;                 // bottom limit for dEdx histogram (default 20)
  Float_t           fDEdxDown;               // upper limit for dEdx histogram (default 1020)
  Float_t           fDEdxCleanUp;            // upper limit for dEdx histogram of clean kaons and electrons (default 140)

  Float_t           fArmPodTPCSignal;
  Float_t           fArmPodptot;
  Float_t           fArmPodEta;
  Float_t           fArmPodCentrality;
  Float_t           fQt;
  Float_t           fAlfa;
  Float_t           fCosPA;

  Float_t           fNSigmasElTPC;           // TPC N sigma for Electron
  Float_t           fNSigmasPiTPC;           // TPC N sigma for Pion
  Float_t           fNSigmasKaTPC;           // TPC N sigma for Kaon
  Float_t           fNSigmasPrTPC;           // TPC N sigma for Proton
  Float_t           fNSigmasDeTPC;           // TPC N sigma for Deuteron

  Float_t           fDEdxEl;                 // Expected Electron dEdx
  Float_t           fDEdxKa;                 // Expected Kaon dEdx
  Float_t           fDEdxPi;                 // Expected Pion dEdx
  Float_t           fDEdxPr;                 // Expected Proton dEdx
  Float_t           fDEdxDe;                 // Expected Deuteron dEdx

  Float_t           fSigmaEl;                // Expected Electron sigma
  Float_t           fSigmaKa;                // Expected Kaon sigma
  Float_t           fSigmaPi;                // Expected Pion sigma
  Float_t           fSigmaPr;                // Expected Proton sigma
  Float_t           fSigmaDe;                // Expected Deuteron sigma

  Float_t           fNSigmasElTOF;           // TOF N sigma for Electron
  Float_t           fNSigmasPiTOF;           // TOF N sigma for Pion
  Float_t           fNSigmasKaTOF;           // TOF N sigma for Kaon
  Float_t           fNSigmasPrTOF;           // TOF N sigma for Proton
  Float_t           fNSigmasDeTOF;           // TOF N sigma for Deuteron

  Float_t           fTOFSignalEl;            // Expected Electron TOF signal
  Float_t           fTOFSignalKa;            // Expected Kaon TOF signal
  Float_t           fTOFSignalPi;            // Expected Pion TOF signal
  Float_t           fTOFSignalPr;            // Expected Proton TOF signal
  Float_t           fTOFSignalDe;            // Expected Deuteron TOF signal

  Float_t           fTOFSigmaEl;             // Expected Electron TOF sigma
  Float_t           fTOFSigmaKa;             // Expected Kaon TOF sigma
  Float_t           fTOFSigmaPi;             // Expected Pion TOF sigma
  Float_t           fTOFSigmaPr;             // Expected Proton TOF sigma
  Float_t           fTOFSigmaDe;             // Expected Deuteron TOF sigma

  Float_t           fNSigmasPrITS;           // ITS N sigma for Proton
  Float_t           fNSigmasKaITS;           // ITS N sigma for Kaon

  Float_t           fTPCSignalMC;
  Float_t           fPtotMC;
  Float_t           fPtotMCtruth;
  Float_t           fPtMC;
  Float_t           fEtaMC;
  Int_t             fSignMC;

  Float_t           fPxMC;                     // x component of momentum
  Float_t           fPyMC;                     // y component of momentum
  Float_t           fPzMC;                     // z component of momentum

  Float_t           fElMC;
  Float_t           fPiMC;
  Float_t           fKaMC;
  Float_t           fPrMC;
  Float_t           fDeMC;
  Float_t           fMuMC;
  Float_t           fLaMC;

  Double_t          fMCImpactParameter;
  Int_t             fNHardScatters;            // Number of hard scatterings
  Int_t             fNProjectileParticipants;  // Number of projectiles participants
  Int_t             fNTargetParticipants;      // Number of target participants
  Int_t             fNNColl;                   // Number of N-N collisions
  Int_t             fNNwColl;                  // Number of N-Nwounded collisions
  Int_t             fNwNColl;                  // Number of Nwounded-N collisons
  Int_t             fNwNwColl;                 // Number of Nwounded-Nwounded collisions


  Float_t           fElMCgen;
  Float_t           fXiMCgen;
  Float_t           fLaMCgen;
  Float_t           fPiMCgen;
  Float_t           fKaMCgen;
  Float_t           fPrMCgen;
  Float_t           fDeMCgen;
  Float_t           fMuMCgen;
  Float_t           fBaMCgen;
  Float_t           fPhMCgen;


  Float_t           fPx;                     // x component of momentum
  Float_t           fPy;                     // y component of momentum
  Float_t           fPz;                     // z component of momentum
  Float_t           fPtot;                   // TPC momentum
  Float_t           fPVertex;                // TPC momentum
  Float_t           fPt;                     // Transverse momentum
  Float_t           fY;                      // rapidity

  Int_t              fMultiplicity;           // Multiplicity in case of PbPb
  Int_t              fMultiplicityMC;
  Float_t            fCentrality;             // centrality information
  Float_t            fCentImpBin;
  Double_t           fVz;                     // Vertex position
  ULong64_t          fEventGID;               // global Event Id
  Int_t              fEventGIDMC;             // global MC event id
  Int_t              fEventCountInFile;       // event count per job
  Int_t              fEvent;                  // Event counter for Christian
  Int_t              fEventMC;                // Event id for MC data
  Int_t              fEventMCgen;             // Event id for MC generated
  std::vector<std::vector<Int_t>> fTrackCounter;      // counter for tracks in different acceptances

  Float_t            fTPCSignal;              // Measured dE/dx
  Float_t            fEta;                    // pseudo rapidity
  Float_t            fNContributors;          // Ntracks
  Float_t            fTheta;                  // theta
  Float_t            fPhi;                    // azimut angle
  Int_t              fSign;                   // sign of the particle
  Int_t              fTPCShared;              // number of shared clusters
  Int_t              fTPCFindable;            // number of findable clusters
  Int_t              fNcl;                    // number of points used for dEdx
  Int_t              fNclCorr;                // number of points used for dEdx
  Float_t            fITSSignal;              // ITS signal

  Int_t              fNResBins;
  Int_t              fNBarBins;
  Int_t              fNEtaWinBinsMC;
  Int_t              fNMomBinsMC;
  Int_t              fNCentBinsMC;

  std::vector<Double_t> fEffMatrixMomBins;
  std::vector<Double_t> fEffMatrixCentBins;
  std::vector<Double_t> fEffMatrixEtaBins;

  std::vector<Double_t> fNSigmaTPC;              // n sigma TPC for cut based method
  std::vector<Double_t> fNSigmaTOFDown;
  std::vector<Double_t> fNSigmaTOFUp;
  Double_t              fTOFMomCut;

  Int_t              fNResModeMC;
  Int_t              fNCentbinsData;
  Float_t            fMissingCl;
  Int_t              fTPCMult;
  Int_t              fEventMult;
  Double_t           fTimeStamp;
  Float_t            fIntRate;
  Int_t              fRunNo;
  Float_t            fBField;
  TString            fBeamType;
  Bool_t             fIsMCPileup;
  Int_t              fMCGeneratorIndex;
  //
  // Jet variables
  Float_t            fLeadingJetCut;
  Double_t           fJetPt;
  Double_t           fJetEta;
  Double_t           fJetPhi;
  Float_t            fjetRhoVal;
  Float_t            fRhoFJ;
  Bool_t             fhasAcceptedFJjet;
  Bool_t             fhasRealFJjet;
  Int_t              fFillJetsBG;          // switch whether to fill jetsEMC constituent tree
  Int_t              fTaskSelection;       // switch whether to fill jetsEMC constituent tree
  TH1F              *fJetHistptSub;        // control histogram for rho subtrated jet pt

  //
  // Event shape
  Double_t fEP_2_Qx_neg;
  Double_t fEP_2_Qx_pos;
  Double_t fEP_2_Qy_neg;
  Double_t fEP_2_Qy_pos;
  Double_t fEP_2_Psi_pos;
  Double_t fEP_2_Psi_neg;
  Double_t fEP_2_Psi;
  Double_t fEP_ntracks_neg;
  Double_t fEP_ntracks_pos;
  //
  Double_t fEP_3_Qx_neg;
  Double_t fEP_3_Qx_pos;
  Double_t fEP_3_Qy_neg;
  Double_t fEP_3_Qy_pos;
  Double_t fEP_3_Psi_pos;
  Double_t fEP_3_Psi_neg;
  Double_t fEP_3_Psi;
  //
  Double_t fFlatenicity;
  Double_t fFlatenicityScaled;
  Double_t fSpherocity;

  //
  // Cut variables
  Double_t fTrackProbElTPC;
  Double_t fTrackProbPiTPC;
  Double_t fTrackProbKaTPC;
  Double_t fTrackProbPrTPC;
  Bool_t   fTrackProbDeTPC;
  Double_t fTrackProbElTOF;
  Double_t fTrackProbPiTOF;
  Double_t fTrackProbKaTOF;
  Double_t fTrackProbPrTOF;
  Bool_t   fTrackProbDeTOF;

  Float_t fTrackTPCCrossedRows;
  Float_t fTrackChi2TPC;
  Float_t fTrackChi2TPCcorr;
  Float_t fTrackDCAxy;
  Float_t fTrackDCAz;
  Float_t fTrackLengthInActiveZone;
  Float_t fTrackTPCSignalN;
  Bool_t  fTrackIsFirstITSlayer;
  Bool_t  fTrackIsSecondITSlayer;
  Bool_t  fTrackNewITScut;
  Bool_t  fTrackRequireITSRefit;

  // Additional cuts from marian
  Bool_t             fIsITSpixel01;           // if track has hits in innermost 2 pixels of ITS
  Int_t              fNITSclusters;           // number of ITS clusters
  Float_t            fPrimRestriction;        // prim vertex cut recommended by marian
  Float_t            fTPCvZ;                  // TPC vertex
  Float_t            fSPDvZ;                  // SPD vertex

  //   CleanSample cuts
  Bool_t             fCleanPionsFromK0;
  Bool_t             fCleanPion0FromK0;
  Bool_t             fCleanPion1FromK0;
  Bool_t             fCleanPion0FromLambda;
  Bool_t             fCleanPion1FromLambda;
  Bool_t             fCleanProton0FromLambda;
  Bool_t             fCleanProton1FromLambda;
  Bool_t             fHasTrack0FirstITSlayer;
  Bool_t             fHasTrack1FirstITSlayer;
  Bool_t             fHasV0FirstITSlayer;

  //  Variables for systematic uncertainty checks
  //  B field configurations -->  use default settings and analyse the following set of runs
  //  ------------------------------------------------
  //  Field (++)  --> run interval is [137161, 138275]
  //  Field (--)  --> run interval is [138364, 139510]
  //  ------------------------------------------------
  Int_t              fSystCentEstimatetor;   // 0 --> "V0M"   ||| -1 -->  "TRK" ||| +1 --> "CL1"
  Float_t            fNetPiFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetKaFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetPrFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetLaFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetChFirstMoments[2][4][10][8];    //[fNResModeMC][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]

  Float_t            fNetPiFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetKaFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetPrFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetPiFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetKaFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fNetPrFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]

  Float_t            fCrossPiFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossKaFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossPrFirstMomentsRec[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossPiFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossKaFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fCrossPrFirstMomentsGen[4][10][8];    //[fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]

  Float_t            fPiFirstMomentsGen[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fKaFirstMomentsGen[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fPrFirstMomentsGen[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fPiFirstMomentsRec[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fKaFirstMomentsRec[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]
  Float_t            fPrFirstMomentsRec[2][4][10][8];    //[2][fNMomBinsMC][fNCentBinsMC][fNEtaWinBinsMC]


  std::vector<float>  fetaDownArr;
  std::vector<float>  fetaUpArr;
  std::vector<float>  fcentDownArr;
  std::vector<float>  fcentUpArr;
  std::vector<float>  fpDownArr;
  std::vector<float>  fpUpArr;
  std::vector<float>  fxCentBins;
  std::vector<std::string> fResonances;
  std::vector<int>    fBaryons;
  std::vector<float>  fCounterEtaBins;
  std::vector<float>  fCounterMomBins;

  //
  // control and QA histograms
  //
  THnF             * fHistPosEffMatrixRec;       // histogram efficiency matrix --> reconstructed traks
  THnF             * fHistNegEffMatrixRec;       // histogram efficiency matrix --> generated traks
  THnF             * fHistPosEffMatrixGen;       // histogram efficiency matrix --> reconstructed pions
  THnF             * fHistNegEffMatrixGen;       // histogram efficiency matrix --> generated pions
  THnF             * fHistPosEffMatrixScanRec;       // histogram efficiency matrix --> reconstructed traks
  THnF             * fHistNegEffMatrixScanRec;       // histogram efficiency matrix --> generated traks
  THnF             * fHistPosEffMatrixScanGen;       // histogram efficiency matrix --> reconstructed pions
  THnF             * fHistNegEffMatrixScanGen;       // histogram efficiency matrix --> generated pions
  TH1F             * fHistEmptyEvent;         // control histogram for empty event
  TH1F             * fHistCentrality;         // control histogram for centrality
  TH1F             * fHistCentralityImpPar;         // control histogram for centrality
  TH1F             * fHistImpParam;           // control histogram for impact parameter
  TH1F             * fHistVertex;             // control histogram for vertexZ
  TH1D             * fHistReturns;               // control histogram for reasons of returns
  TH1D             * fHistPileup;                // control histogram for pileup cut
  TH3D             * fHistPileup2D;                // control histogram for pileup cut
  TH2F             * fHistArmPod;             // control histogram for Armanteros Podolanski plot
  THnF             * fEffMatrixGenPos;           // histogram efficiency matrix read from file
  THnF             * fEffMatrixGenNeg;           // histogram efficiency matrix read from file
  THnF             * fEffMatrixRecPos;           // histogram efficiency matrix read from file
  THnF             * fEffMatrixRecNeg;           // histogram efficiency matrix read from file
  TH2F             * fHistEP2QxQypos;           // histogram efficiency matrix read from file
  TH2F             * fHistEP2QxQyneg;           // histogram efficiency matrix read from file
  TH2F             * fHistEP3QxQypos;           // histogram efficiency matrix read from file
  TH2F             * fHistEP3QxQyneg;           // histogram efficiency matrix read from file

  std::vector<std::vector<std::vector<std::vector<TH2F*>>>> fEffMatrixProjections;  // container for efficiency matrix projections
  //
  // Counters for Marian
  //
  TVectorF         * fEventInfo_PhiTPCdcarA;  // track counter
  TVectorF         * fEventInfo_PhiTPCdcarC; // dedx info counter
  TVectorF         * fEventInfo_CacheTrackCounters;  // track counter
  TVectorF         * fEventInfo_CacheTrackdEdxRatio; // dedx info counter
  TVectorF         * fEventInfo_CacheTrackNcl;       // ncl counter
  TVectorF         * fEventInfo_CacheTrackChi2;      // chi2 counter
  TVectorF         * fEventInfo_CacheTrackMatchEff;  // matchEff counter
  TVectorF         * fEventInfo_CentralityEstimates;
  TGraph           * fEventInfo_LumiGraph;           // grap for the interaction rate info for a run
  TH1F             * fEventInfo_HisTPCVertexA;
  TH1F             * fEventInfo_HisTPCVertexC;
  TH1F             * fEventInfo_HisTPCVertexACut;
  TH1F             * fEventInfo_HisTPCVertexCCut;
  TH1F             * fEventInfo_HisTPCVertex;
  TVectorF         * fEventInfo_CacheTrackTPCCountersZ; // track counter with DCA z cut
  static const char*  fEventInfo_centEstStr[];              //!centrality types

  AliEventCuts* fPileUpTightnessCut4;
  AliEventCuts* fPileUpTightnessCut3;
  AliEventCuts* fPileUpTightnessCut2;
  AliEventCuts* fPileUpTightnessCut1;

  ClassDef(AliAnalysisTaskTIdentityPID, 8);

};

#endif
