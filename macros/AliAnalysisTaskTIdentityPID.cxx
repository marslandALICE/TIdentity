/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
// Analysis for net particle fluctuations and jet-substructure studies   //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGrid.h"
#include "TSystem.h"
#include "TCutG.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH3.h"
#include "THn.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TMatrixF.h"
#include "TVectorF.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "AliPDG.h"
#include "AliMathBase.h"
#include "AliESDFMD.h"
#include "AliFMDFloatMap.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliTPCdEdxInfo.h"
#include "AliESDv0KineCuts.h"
#include "AliKFVertex.h"
#include "AliLumiTools.h"
#include "AliKFParticle.h"
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEposEventHeader.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDv0Cuts.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliCentrality.h"
#include "AliESDUtils.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliTPCParam.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliKFParticle.h"
#include "AliAnalysisTaskFilteredTree.h"
#include "AliAnalysisTaskTIdentityPID.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliRunLoader.h"
#include "AliEventCuts.h"
#include "AliAnalysisUtils.h"
#include "AliESDtools.h"
#include "AliFJWrapper.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <bitset>
using namespace std;
using std::cout;
using std::setw;



ClassImp(AliAnalysisTaskTIdentityPID)

const char* AliAnalysisTaskTIdentityPID::fEventInfo_centEstStr[] = {"V0M","CL0","CL1"};
// {"V0M","V0A","V0C","FMD","TRK","TKL","CL0","CL1","V0MvsFMD","ZNA","TKLvsV0M","ZEMvsZDC","V0A123","V0A0","V0S", "MB", "Ref", "V0av"};
// const char* AliAnalysisTaskTIdentityPID::fEventInfo_centEstStr[] = {"V0M","CL0","CL1","TRK","TKL","V0MvsFMD","TKLvsV0M","ZEMvsZDC","kV0A"
// ,"V0C","ZNA","ZNC","ZPA","ZPC","CND","FMD","NPA","V0A0","V0A123","V0A23","V0C01","V0S","V0MEq","V0AEq","V0CEq","SPDClusters","SPDTracklets"};

#define USE_STREAMER 1


// -----------------------------------------------------------------------
//                            Constructor and Destructor
// -----------------------------------------------------------------------
//________________________________________________________________________
AliAnalysisTaskTIdentityPID::AliAnalysisTaskTIdentityPID()
: AliAnalysisTaskSE("TaskEbyeRatios"), fEventCuts(0), fPIDResponse(0),fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_Bit96(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCutsLoose(0),
fESDtrackCutsV0(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fV0OpenCuts(0x0),
fV0StrongCuts(0x0),
fK0sPionCuts(0x0),
fLambdaProtonCuts(0x0),
fLambdaPionCuts(0x0),
fGammaElectronCuts(0x0),
fVertex(0x0),
fESDtool(nullptr),
fArmPodTree(0x0),
fTreeSRedirector(0x0),
fTreeMomentsMCfull(0x0),
fTreeTracksMCgen(0x0),
fTreeDebug3(0x0),
fTreeTracksMCrec(0x0),
fTreeDebug4(0x0),
fTreeTracks(0x0),
fTreeDebug(0x0),
fTreeResonances(0x0),
fTreeMomentsMCgen(0x0),
fTreeEvents(0x0),
fTreeEventsMC(0x0),
fTreeDScaled(0x0),
fTreeDebug2(0x0),
fTreeExpecteds(0x0),
fTreeCutBased(0x0),
fTreejetsFJ(0x0),
fTreejetsFJBG(0x0),
fTreejetsFJconst(0x0),
fTreejetsFJGen(0x0),
fTreejetsFJBGGen(0x0),
fTreejetsFJconstGen(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPileUpBit(0),
fHistCent(0),
fHistPhi(0),
fHistGenMult(0),
fHistRapDistFullAccPr(0),
fHistRapDistFullAccAPr(0),
fHistInvK0s(0),
fHistInvLambda(0),
fHistInvAntiLambda(0),
fHistInvPhoton(0),
fHistPhiTPCcounterA(0),
fHistPhiTPCcounterC(0),
fHistPhiTPCcounterAITS(0),
fHistPhiTPCcounterCITS(0),
fHistPhiITScounterA(0),
fHistPhiITScounterC(0),
fH2MissCl(),
fChunkName(""),
fTrackCutBits(0),
fSystClass(0),
fEtaMin(0),
fEtaMax(0),
fNEtaBins(0),
fDownscalingFactor(0),
fEventIDinFile(0),
fChunkIDinJob(0),
fRunOnGrid(kFALSE),
fMCtrue(kFALSE),
fEventInfo(kFALSE),
fWeakAndMaterial(kFALSE),
fFillEffMatrix(kFALSE),
fFillDebug(kFALSE),
fDEdxCheck(kFALSE),
fIncludeITS(kTRUE),
fFillTracksMCgen(kFALSE),
fFillEffLookUpTable(kFALSE),
fFillHigherMomentsMCclosure(kFALSE),
fFillArmPodTree(kTRUE),
fRunFastSimulation(kFALSE),
fRunFastHighMomentCal(kFALSE),
fRunCutBasedMethod(kFALSE),
fFillDistributions(kFALSE),
fFillTreeMC(kFALSE),
fDefaultTrackCuts(kFALSE),
fDefaultEventCuts(kFALSE),
fFillNudynFastGen(kFALSE),
fFillResonances(kFALSE),
fCollisionType(0),
fCorrectForMissCl(0),
fUsePtCut(1),
fTrackOriginOnlyPrimary(0),
fRapidityType(0),
fSisterCheck(0),
fIncludeTOF(kFALSE),
fUseCouts(kFALSE),
fV0InvMassHists(kFALSE),
fRunNumberForExpecteds(0),
fFillExpecteds(kFALSE),
fDefaultCuts(kFALSE),
fNSettings(18),
fSystSettings(0),
fNMomBins(0),
fMomDown(0),
fMomUp(0),
fDEdxBinWidth(0),
fDEdxUp(0),
fDEdxDown(0),
fDEdxCleanUp(0),
fArmPodTPCSignal(0),
fArmPodptot(0),
fArmPodEta(0),
fArmPodCentrality(0),
fQt(0),
fAlfa(0),
fCosPA(0),
fNSigmasElTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtotMCtruth(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fPxMC(0),
fPyMC(0),
fPzMC(0),
fElMC(0),
fPiMC(0),
fKaMC(0),
fPrMC(0),
fDeMC(0),
fMuMC(0),
fLaMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fElMCgen(0),
fXiMCgen(0),
fPiMCgen(0),
fKaMCgen(0),
fPrMCgen(0),
fDeMCgen(0),
fMuMCgen(0),
fLaMCgen(0),
fBaMCgen(0),
fPhMCgen(0),
fPx(0),
fPy(0),
fPz(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fMultiplicity(0),
fMultiplicityMC(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventGID(0),
fEventGIDMC(0),
fEventCountInFile(0),
fEvent(0),
fEventMC(0),
fEventMCgen(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fTheta(0),
fPhi(0),
fSign(0),
fTPCShared(0),
fTPCFindable(0),
fNcl(0),
fNclCorr(0),
fNResBins(0),
fNBarBins(0),
fNEtaWinBinsMC(-100),
fNMomBinsMC(-100),
fNCentBinsMC(-100),
fGenprotonBins(-100),
fEffMatrixMomBins(0),
fEffMatrixCentBins(0),
fEffMatrixEtaBins(0),
fNSigmaTPC(0),
fNSigmaTOFDown(0),
fNSigmaTOFUp(0),
fTOFMomCut(0),
fNResModeMC(2),
fNCentbinsData(14),
fMissingCl(0.),
fTPCMult(0),
fEventMult(0),
fTimeStamp(0),
fIntRate(0),
fRunNo(0),
fBField(0),
fBeamType(0),
fIsMCPileup(0),
fMCGeneratorIndex(0),
fLeadingJetCut(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
fRhoFJ(0),
fhasAcceptedFJjet(0),
fhasRealFJjet(0),
fFillJetsBG(0),
fTaskSelection(0),
fJetHistptSub(0),
fEP_2_Qx_neg(0),
fEP_2_Qx_pos(0),
fEP_2_Qy_neg(0),
fEP_2_Qy_pos(0),
fEP_2_Psi_pos(0),
fEP_2_Psi_neg(0),
fEP_2_Psi(0),
fEP_ntracks_neg(0),
fEP_ntracks_pos(0),
fEP_3_Qx_neg(0),
fEP_3_Qx_pos(0),
fEP_3_Qy_neg(0),
fEP_3_Qy_pos(0),
fEP_3_Psi_pos(0),
fEP_3_Psi_neg(0),
fEP_3_Psi(0),
fFlatenicity(0),
fFlatenicityScaled(0),
fSpherocity(0),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fTrackTPCCrossedRows(0),
fTrackChi2TPC(0),
fTrackChi2TPCcorr(0),
fTrackDCAxy(0),
fTrackDCAz(0),
fTrackLengthInActiveZone(0),
fTrackTPCSignalN(0),
fTrackIsFirstITSlayer(0),
fTrackIsSecondITSlayer(0),
fTrackNewITScut(0),
fTrackRequireITSRefit(0),
fIsITSpixel01(0),
fNITSclusters(0),
fPrimRestriction(0),
fTPCvZ(0),
fSPDvZ(0),
fCleanPionsFromK0(0),
fCleanPion0FromK0(0),
fCleanPion1FromK0(0),
fCleanPion0FromLambda(0),
fCleanPion1FromLambda(0),
fCleanProton0FromLambda(0),
fCleanProton1FromLambda(0),
fHasTrack0FirstITSlayer(0),
fHasTrack1FirstITSlayer(0),
fHasV0FirstITSlayer(0),
fSystCentEstimatetor(0),
fetaDownArr(),
fetaUpArr(),
fcentDownArr(),
fcentUpArr(),
fpDownArr(),
fpUpArr(),
fxCentBins(),
fResonances(),
fBaryons(),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistPosEffMatrixScanRec(0),
fHistNegEffMatrixScanRec(0),
fHistPosEffMatrixScanGen(0),
fHistNegEffMatrixScanGen(0),
fEffMatrixProjections(0),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistCentralityImpPar(0),
fHistImpParam(0),
fHistVertex(0),
fHistReturns(0),
fHistPileup(0),
fHistPileup2D(0),
fHistArmPod(0),
fEffMatrixGenPos(0),
fEffMatrixGenNeg(0),
fEffMatrixRecPos(0),
fEffMatrixRecNeg(0),
fEventInfo_PhiTPCdcarA(0),
fEventInfo_PhiTPCdcarC(0),
fEventInfo_CacheTrackCounters(0),
fEventInfo_CacheTrackdEdxRatio(0),
fEventInfo_CacheTrackNcl(0),
fEventInfo_CacheTrackChi2(0),
fEventInfo_CacheTrackMatchEff(0),
fEventInfo_CentralityEstimates(0),
fEventInfo_LumiGraph(0),
fEventInfo_HisTPCVertexA(0),
fEventInfo_HisTPCVertexC(0),
fEventInfo_HisTPCVertexACut(0),
fEventInfo_HisTPCVertexCCut(0),
fEventInfo_HisTPCVertex(0),
fEventInfo_CacheTrackTPCCountersZ(0),
fPileUpTightnessCut1(0),
fPileUpTightnessCut2(0),
fPileUpTightnessCut3(0),
fPileUpTightnessCut4(0)
{
  // default Constructor
  /* fast compilation test
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/AliAnalysisTaskTIdentityPID.cxx++
  .L /lustre/nyx/alice/users/marsland/train/trunk/marsland_EbyeRatios/AliAnalysisTaskTIdentityPID.cxx++
  */
}

//________________________________________________________________________
AliAnalysisTaskTIdentityPID::AliAnalysisTaskTIdentityPID(const char *name)
: AliAnalysisTaskSE(name), fEventCuts(0), fPIDResponse(0), fESD(0), fListHist(0),
fESDtrackCuts(0),
fESDtrackCuts_Bit96(0),
fESDtrackCuts_Bit128(0),
fESDtrackCuts_Bit768(0),
fESDtrackCutsLoose(0),
fESDtrackCutsV0(0),
fESDtrackCutsCleanSamp(0),
fPIDCombined(0x0),
fTPCdEdxInfo(0x0),
fMCStack(0x0),
fV0OpenCuts(0x0),
fV0StrongCuts(0x0),
fK0sPionCuts(0x0),
fLambdaProtonCuts(0x0),
fLambdaPionCuts(0x0),
fGammaElectronCuts(0x0),
fVertex(0x0),
fESDtool(nullptr),
fArmPodTree(0x0),
fTreeSRedirector(0x0),
fTreeMomentsMCfull(0x0),
fTreeTracksMCgen(0x0),
fTreeDebug3(0x0),
fTreeTracksMCrec(0x0),
fTreeDebug4(0x0),
fTreeTracks(0x0),
fTreeDebug(0x0),
fTreeResonances(0x0),
fTreeMomentsMCgen(0x0),
fTreeEvents(0x0),
fTreeEventsMC(0x0),
fTreeDScaled(0x0),
fTreeDebug2(0x0),
fTreeExpecteds(0x0),
fTreeCutBased(0x0),
fTreejetsFJ(0x0),
fTreejetsFJBG(0x0),
fTreejetsFJconst(0x0),
fTreejetsFJGen(0x0),
fTreejetsFJBGGen(0x0),
fTreejetsFJconstGen(0x0),
fRandom(0),
fPeriodName(""),
fYear(0),
fPassIndex(0),
fPileUpBit(0),
fHistCent(0),
fHistPhi(0),
fHistGenMult(0),
fHistRapDistFullAccPr(0),
fHistRapDistFullAccAPr(0),
fHistInvK0s(0),
fHistInvLambda(0),
fHistInvAntiLambda(0),
fHistInvPhoton(0),
fHistPhiTPCcounterA(0),
fHistPhiTPCcounterC(0),
fHistPhiTPCcounterAITS(0),
fHistPhiTPCcounterCITS(0),
fHistPhiITScounterA(0),
fHistPhiITScounterC(0),
fH2MissCl(),
fChunkName(""),
fTrackCutBits(0),
fSystClass(0),
fEtaMin(0),
fEtaMax(0),
fNEtaBins(0),
fDownscalingFactor(0),
fEventIDinFile(0),
fChunkIDinJob(0),
fRunOnGrid(kFALSE),
fMCtrue(kFALSE),
fEventInfo(kFALSE),
fWeakAndMaterial(kFALSE),
fFillEffMatrix(kFALSE),
fFillDebug(kFALSE),
fDEdxCheck(kFALSE),
fIncludeITS(kTRUE),
fFillTracksMCgen(kFALSE),
fFillEffLookUpTable(kFALSE),
fFillHigherMomentsMCclosure(kFALSE),
fFillArmPodTree(kTRUE),
fRunFastSimulation(kFALSE),
fRunFastHighMomentCal(kFALSE),
fRunCutBasedMethod(kFALSE),
fFillDistributions(kFALSE),
fFillTreeMC(kFALSE),
fDefaultTrackCuts(kFALSE),
fDefaultEventCuts(kFALSE),
fFillNudynFastGen(kFALSE),
fFillResonances(kFALSE),
fCollisionType(0),
fCorrectForMissCl(0),
fUsePtCut(1),
fTrackOriginOnlyPrimary(0),
fRapidityType(0),
fSisterCheck(0),
fIncludeTOF(kFALSE),
fUseCouts(kFALSE),
fV0InvMassHists(kFALSE),
fRunNumberForExpecteds(0),
fFillExpecteds(kFALSE),
fDefaultCuts(kFALSE),
fNSettings(18),
fSystSettings(0),
fNMomBins(0),
fMomDown(0),
fMomUp(0),
fDEdxBinWidth(0),
fDEdxUp(0),
fDEdxDown(0),
fDEdxCleanUp(0),
fArmPodTPCSignal(0),
fArmPodptot(0),
fArmPodEta(0),
fArmPodCentrality(0),
fQt(0),
fAlfa(0),
fCosPA(0),
fNSigmasElTOF(0),
fNSigmasPiTOF(0),
fNSigmasKaTOF(0),
fNSigmasPrTOF(0),
fNSigmasDeTOF(0),
fDEdxEl(0),
fDEdxKa(0),
fDEdxPi(0),
fDEdxPr(0),
fDEdxDe(0),
fSigmaEl(0),
fSigmaKa(0),
fSigmaPi(0),
fSigmaPr(0),
fSigmaDe(0),
fNSigmasElTPC(0),
fNSigmasPiTPC(0),
fNSigmasKaTPC(0),
fNSigmasPrTPC(0),
fNSigmasDeTPC(0),
fTPCSignalMC(0),
fPtotMC(0),
fPtotMCtruth(0),
fPtMC(0),
fEtaMC(0),
fSignMC(0),
fPxMC(0),
fPyMC(0),
fPzMC(0),
fElMC(0),
fPiMC(0),
fKaMC(0),
fPrMC(0),
fDeMC(0),
fMuMC(0),
fLaMC(0),
fMCImpactParameter(0),
fNHardScatters(0),
fNProjectileParticipants(0),
fNTargetParticipants(0),
fNNColl(0),
fNNwColl(0),
fNwNColl(0),
fNwNwColl(0),
fElMCgen(0),
fXiMCgen(0),
fPiMCgen(0),
fKaMCgen(0),
fPrMCgen(0),
fDeMCgen(0),
fMuMCgen(0),
fLaMCgen(0),
fBaMCgen(0),
fPhMCgen(0),
fPx(0),
fPy(0),
fPz(0),
fPtot(0),
fPVertex(0),
fPt(0),
fY(0),
fMultiplicity(0),
fMultiplicityMC(0),
fCentrality(0),
fCentImpBin(0),
fVz(0),
fEventGID(0),
fEventGIDMC(0),
fEventCountInFile(0),
fEvent(0),
fEventMC(0),
fEventMCgen(0),
fTPCSignal(0),
fEta(0),
fNContributors(0),
fTheta(0),
fPhi(0),
fSign(0),
fTPCShared(0),
fTPCFindable(0),
fNcl(0),
fNclCorr(0),
fNResBins(0),
fNBarBins(0),
fNEtaWinBinsMC(-100),
fNMomBinsMC(-100),
fNCentBinsMC(-100),
fGenprotonBins(-100),
fEffMatrixMomBins(0),
fEffMatrixCentBins(0),
fEffMatrixEtaBins(0),
fNSigmaTPC(0),
fNSigmaTOFDown(0),
fNSigmaTOFUp(0),
fTOFMomCut(0),
fNResModeMC(2),
fNCentbinsData(14),
fMissingCl(0.),
fTPCMult(0),
fEventMult(0),
fTimeStamp(0),
fIntRate(0),
fRunNo(0),
fBField(0),
fBeamType(0),
fIsMCPileup(0),
fMCGeneratorIndex(0),
fLeadingJetCut(0),
fJetPt(0),
fJetEta(0),
fJetPhi(0),
fjetRhoVal(0),
fRhoFJ(0),
fhasAcceptedFJjet(0),
fhasRealFJjet(0),
fFillJetsBG(0),
fTaskSelection(0),
fJetHistptSub(0),
fEP_2_Qx_neg(0),
fEP_2_Qx_pos(0),
fEP_2_Qy_neg(0),
fEP_2_Qy_pos(0),
fEP_2_Psi_pos(0),
fEP_2_Psi_neg(0),
fEP_2_Psi(0),
fEP_ntracks_neg(0),
fEP_ntracks_pos(0),
fEP_3_Qx_neg(0),
fEP_3_Qx_pos(0),
fEP_3_Qy_neg(0),
fEP_3_Qy_pos(0),
fEP_3_Psi_pos(0),
fEP_3_Psi_neg(0),
fEP_3_Psi(0),
fFlatenicity(0),
fFlatenicityScaled(0),
fSpherocity(0),
fTrackProbElTPC(0),
fTrackProbPiTPC(0),
fTrackProbKaTPC(0),
fTrackProbPrTPC(0),
fTrackProbDeTPC(0),
fTrackProbElTOF(0),
fTrackProbPiTOF(0),
fTrackProbKaTOF(0),
fTrackProbPrTOF(0),
fTrackProbDeTOF(0),
fTrackTPCCrossedRows(0),
fTrackChi2TPC(0),
fTrackChi2TPCcorr(0),
fTrackDCAxy(0),
fTrackDCAz(0),
fTrackLengthInActiveZone(0),
fTrackTPCSignalN(0),
fTrackIsFirstITSlayer(0),
fTrackIsSecondITSlayer(0),
fTrackNewITScut(0),
fTrackRequireITSRefit(0),
fIsITSpixel01(0),
fNITSclusters(0),
fPrimRestriction(0),
fTPCvZ(0),
fSPDvZ(0),
fCleanPionsFromK0(0),
fCleanPion0FromK0(0),
fCleanPion1FromK0(0),
fCleanPion0FromLambda(0),
fCleanPion1FromLambda(0),
fCleanProton0FromLambda(0),
fCleanProton1FromLambda(0),
fHasTrack0FirstITSlayer(0),
fHasTrack1FirstITSlayer(0),
fHasV0FirstITSlayer(0),
fSystCentEstimatetor(0),
fetaDownArr(),
fetaUpArr(),
fcentDownArr(),
fcentUpArr(),
fpDownArr(),
fpUpArr(),
fxCentBins(),
fResonances(),
fBaryons(),
fHistPosEffMatrixRec(0),
fHistNegEffMatrixRec(0),
fHistPosEffMatrixGen(0),
fHistNegEffMatrixGen(0),
fHistPosEffMatrixScanRec(0),
fHistNegEffMatrixScanRec(0),
fHistPosEffMatrixScanGen(0),
fHistNegEffMatrixScanGen(0),
fEffMatrixProjections(0),
fHistEmptyEvent(0),
fHistCentrality(0),
fHistCentralityImpPar(0),
fHistImpParam(0),
fHistVertex(0),
fHistReturns(0),
fHistPileup(0),
fHistPileup2D(0),
fHistArmPod(0),
fEffMatrixGenPos(0),
fEffMatrixGenNeg(0),
fEffMatrixRecPos(0),
fEffMatrixRecNeg(0),
fEventInfo_PhiTPCdcarA(0),
fEventInfo_PhiTPCdcarC(0),
fEventInfo_CacheTrackCounters(0),
fEventInfo_CacheTrackdEdxRatio(0),
fEventInfo_CacheTrackNcl(0),
fEventInfo_CacheTrackChi2(0),
fEventInfo_CacheTrackMatchEff(0),
fEventInfo_CentralityEstimates(0),
fEventInfo_LumiGraph(0),
fEventInfo_HisTPCVertexA(0),
fEventInfo_HisTPCVertexC(0),
fEventInfo_HisTPCVertexACut(0),
fEventInfo_HisTPCVertexCCut(0),
fEventInfo_HisTPCVertex(0),
fEventInfo_CacheTrackTPCCountersZ(0),
fPileUpTightnessCut1(0),
fPileUpTightnessCut2(0),
fPileUpTightnessCut3(0),
fPileUpTightnessCut4(0)
{
  //
  //         standard constructur which should be used
  //
  /* fast compilation test
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I/lustre/nyx/alice/users/marsland/alicehub/sw/ubuntu1404_x86-64/AliPhysics/master-1/include");
  gSystem->AddIncludePath("-I/lustre/nyx/alice/users/marsland/alicehub/AliPhysics/PWGPP");
  .L /home/marsland/Desktop/RUN_ON_GRID/Ebye/code/AliAnalysisTaskTIdentityPID.cxx++
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/AliAnalysisTaskTIdentityPID.cxx++
  .L /lustre/nyx/alice/users/marsland/train/trunk/marsland_EbyeRatios/AliAnalysisTaskTIdentityPID.cxx++
  */
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:***************** CONSTRUCTOR CALLED: AliAnalysisTaskTIdentityPID  *****************"<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  std::cout << " Info::marsland:===================================================================================="<< std::endl;
  // ==========================================
  //
  // ==========================================
  // Initialize arrays
  for (Int_t ires=0; ires<2; ires++){
    for (Int_t imom=0; imom<4; imom++){
      for (Int_t icent=0; icent<10; icent++){
        for (Int_t ieta=0; ieta<8; ieta++){
          fNetPiFirstMoments[ires][imom][icent][ieta]=0.;
          fNetKaFirstMoments[ires][imom][icent][ieta]=0.;
          fNetPrFirstMoments[ires][imom][icent][ieta]=0.;
          fNetLaFirstMoments[ires][imom][icent][ieta]=0.;
          fNetChFirstMoments[ires][imom][icent][ieta]=0.;
        }
      }
    }
  }
  // Initialize arrays
  for (Int_t imom=0; imom<4; imom++){
    for (Int_t icent=0; icent<10; icent++){
      for (Int_t ieta=0; ieta<8; ieta++){
        fNetPiFirstMomentsGen[imom][icent][ieta]=0.;
        fNetKaFirstMomentsGen[imom][icent][ieta]=0.;
        fNetPrFirstMomentsGen[imom][icent][ieta]=0.;
        fNetPiFirstMomentsRec[imom][icent][ieta]=0.;
        fNetKaFirstMomentsRec[imom][icent][ieta]=0.;
        fNetPrFirstMomentsRec[imom][icent][ieta]=0.;
      }
    }
  }
  for (Int_t imom=0; imom<4; imom++){
    for (Int_t icent=0; icent<10; icent++){
      for (Int_t ieta=0; ieta<8; ieta++){
        fCrossPiFirstMomentsGen[imom][icent][ieta]=0.;
        fCrossKaFirstMomentsGen[imom][icent][ieta]=0.;
        fCrossPrFirstMomentsGen[imom][icent][ieta]=0.;
        fCrossPiFirstMomentsRec[imom][icent][ieta]=0.;
        fCrossKaFirstMomentsRec[imom][icent][ieta]=0.;
        fCrossPrFirstMomentsRec[imom][icent][ieta]=0.;
      }
    }
  }
  // Initialize arrays
  for (Int_t isign=0; isign<2; isign++){
    for (Int_t imom=0; imom<4; imom++){
      for (Int_t icent=0; icent<10; icent++){
        for (Int_t ieta=0; ieta<8; ieta++){
          fPiFirstMomentsGen[isign][imom][icent][ieta]=0.;
          fKaFirstMomentsGen[isign][imom][icent][ieta]=0.;
          fPrFirstMomentsGen[isign][imom][icent][ieta]=0.;
          fPiFirstMomentsRec[isign][imom][icent][ieta]=0.;
          fKaFirstMomentsRec[isign][imom][icent][ieta]=0.;
          fPrFirstMomentsRec[isign][imom][icent][ieta]=0.;
        }
      }
    }
  }

  // Initialize miaaing cluter map
  // const Int_t nParticles=4;
  // const Int_t nCentBins=9;
  // for (Int_t ipart=0; ipart<nParticles; ipart++){
  //   for (Int_t icent=0; icent<nCentBins; icent++){
  //     fH2MissCl[ipart][icent] = 0;
  //   }
  // }
  //
  // ==========================================
  //
  // Define outputs
  DefineOutput(1,  TList::Class());
  DefineOutput(2,  TTree::Class());
  DefineOutput(3,  TTree::Class());
  DefineOutput(4,  TTree::Class());
  DefineOutput(5,  TTree::Class());
  DefineOutput(6,  TTree::Class());
  DefineOutput(7,  TTree::Class());
  DefineOutput(8,  TTree::Class());
  DefineOutput(9,  TTree::Class());
  DefineOutput(10, TTree::Class());
  DefineOutput(11, TTree::Class());
  DefineOutput(12, TTree::Class());
  DefineOutput(13, TTree::Class());
  DefineOutput(14, TTree::Class());
  DefineOutput(15, TTree::Class());
  DefineOutput(16, TTree::Class());
  DefineOutput(17, TTree::Class());
  DefineOutput(18, TTree::Class());
  DefineOutput(19, TTree::Class());
  DefineOutput(20, TTree::Class());
  DefineOutput(21, TTree::Class());
  DefineOutput(22, TTree::Class());
  DefineOutput(23, TTree::Class());

  // ==========================================

}
//________________________________________________________________________
AliAnalysisTaskTIdentityPID::~AliAnalysisTaskTIdentityPID()
{

  //
  // Destructor
  //
  std::cout << " Info::marsland: ===== In the Destructor ===== " << std::endl;
  if (fHistPosEffMatrixRec) delete fHistPosEffMatrixRec;
  if (fHistNegEffMatrixRec) delete fHistNegEffMatrixRec;
  if (fHistPosEffMatrixGen) delete fHistPosEffMatrixGen;
  if (fHistNegEffMatrixGen) delete fHistNegEffMatrixGen;
  if (fHistPosEffMatrixScanRec) delete fHistPosEffMatrixScanRec;
  if (fHistNegEffMatrixScanRec) delete fHistNegEffMatrixScanRec;
  if (fHistPosEffMatrixScanGen) delete fHistPosEffMatrixScanGen;
  if (fHistNegEffMatrixScanGen) delete fHistNegEffMatrixScanGen;
  if (fHistEmptyEvent)      delete fHistEmptyEvent;
  if (fHistCentrality)      delete fHistCentrality;
  if (fHistCentralityImpPar)delete fHistCentralityImpPar;
  if (fHistImpParam)        delete fHistImpParam;
  if (fHistVertex)          delete fHistVertex;
  if (fHistArmPod)          delete fHistArmPod;
  if (fHistCent)            delete fHistCent;
  if (fHistPhi)             delete fHistPhi;
  if (fHistGenMult)         delete fHistGenMult;
  if (fHistRapDistFullAccPr)   delete fHistRapDistFullAccPr;
  if (fHistRapDistFullAccAPr)  delete fHistRapDistFullAccAPr;
  if (fHistInvK0s)          delete fHistInvK0s;
  if (fHistInvLambda)       delete fHistInvLambda;
  if (fHistInvAntiLambda)   delete fHistInvAntiLambda;
  if (fHistInvPhoton)       delete fHistInvPhoton;
  //
  // Marians histograms
  if (fHistPhiTPCcounterA)  delete fHistPhiTPCcounterA;
  if (fHistPhiTPCcounterC)  delete fHistPhiTPCcounterC;
  if (fHistPhiTPCcounterAITS)  delete fHistPhiTPCcounterAITS;
  if (fHistPhiTPCcounterCITS)  delete fHistPhiTPCcounterCITS;
  if (fHistPhiITScounterA)  delete fHistPhiITScounterA;
  if (fHistPhiITScounterC)  delete fHistPhiITScounterC;
  if (fEventInfo_HisTPCVertexA)  delete fEventInfo_HisTPCVertexA;
  if (fEventInfo_HisTPCVertexC)  delete fEventInfo_HisTPCVertexC;
  if (fEventInfo_HisTPCVertex)  delete fEventInfo_HisTPCVertex;
  if (fEventInfo_HisTPCVertexACut)  delete fEventInfo_HisTPCVertexACut;
  if (fEventInfo_HisTPCVertexCCut)  delete fEventInfo_HisTPCVertexCCut;
  if (fEventInfo_PhiTPCdcarA)  delete fEventInfo_PhiTPCdcarA;
  if (fEventInfo_PhiTPCdcarC)  delete fEventInfo_PhiTPCdcarC;
  if (fEventInfo_CacheTrackCounters)  delete fEventInfo_CacheTrackCounters;
  if (fEventInfo_CacheTrackdEdxRatio)  delete fEventInfo_CacheTrackdEdxRatio;
  if (fEventInfo_CacheTrackNcl)  delete fEventInfo_CacheTrackNcl;
  if (fEventInfo_CacheTrackChi2)  delete fEventInfo_CacheTrackChi2;
  if (fEventInfo_CacheTrackMatchEff)  delete fEventInfo_CacheTrackMatchEff;
  if (fEventInfo_CentralityEstimates)  delete fEventInfo_CentralityEstimates;
  if (fEventInfo_CacheTrackTPCCountersZ)  delete fEventInfo_CacheTrackTPCCountersZ;
  if (fPIDCombined) delete fPIDCombined;
  if (fESDtrackCuts)          delete fESDtrackCuts;
  if (fESDtrackCuts_Bit96)    delete fESDtrackCuts_Bit96;
  if (fESDtrackCuts_Bit128)   delete fESDtrackCuts_Bit128;
  if (fESDtrackCuts_Bit768)   delete fESDtrackCuts_Bit768;
  if (fESDtrackCutsLoose)     delete fESDtrackCutsLoose;
  if (fESDtrackCutsV0)        delete fESDtrackCutsV0;
  if (fTreeSRedirector)       delete fTreeSRedirector;
  if (fESDtrackCutsCleanSamp) delete fESDtrackCutsCleanSamp;
  if (fPileUpTightnessCut4) delete fPileUpTightnessCut4;
  if (fPileUpTightnessCut3) delete fPileUpTightnessCut3;
  if (fPileUpTightnessCut2) delete fPileUpTightnessCut2;
  if (fPileUpTightnessCut1) delete fPileUpTightnessCut1;




}
//
// ---------------------------------------------------------------------------------
//                                     Functions
// ---------------------------------------------------------------------------------
//
void AliAnalysisTaskTIdentityPID::Initialize()
{
  //
  // updating parameters in case of changes (standard cuts and the eta window)
  //
  std::cout << " Info::marsland: ===== In the Initialize ===== " << std::endl;
  if (fRunFastSimulation)    { std::cout << " Info::marsland: !!! We are running fast simulation return !!! " << std::endl; return; }
  if (fRunFastHighMomentCal) { std::cout << " Info::marsland: !!! We are running fast high moment calculation return !!! " << std::endl; return; }
  //
  // fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  //
  // ------------------------------------------------
  //
  // tight DCA cut used by Emil
  fESDtrackCuts_Bit96     = new AliESDtrackCuts("fESDtrackCuts_Bit96","");
  fESDtrackCuts_Bit96->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit96->SetPtRange(0.,100000.);
  fESDtrackCuts_Bit96->SetMinNClustersTPC(70); // ???? should be --> fESDtrackCuts_Bit96->SetMinNCrossedRowsTPC(70);
  fESDtrackCuts_Bit96->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts_Bit96->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts_Bit96->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts_Bit96->SetRequireITSRefit(kTRUE);
  fESDtrackCuts_Bit96->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fESDtrackCuts_Bit96->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
  fESDtrackCuts_Bit96->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts_Bit96->SetMaxDCAToVertexZ(2);
  fESDtrackCuts_Bit96->SetDCAToVertex2D(kFALSE);
  fESDtrackCuts_Bit96->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts_Bit96->SetMaxChi2PerClusterITS(36);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit96->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit96->SetMaxChi2PerClusterTPC(4);
  }
  //
  // TPC only tracks
  fESDtrackCuts_Bit128 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fESDtrackCuts_Bit128->SetName("Bit128");
  fESDtrackCuts_Bit128->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit128->SetPtRange(0.,100000.);
  fESDtrackCuts_Bit128->SetMinNClustersTPC(70);
  fESDtrackCuts_Bit128->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts_Bit128->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts_Bit128->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts_Bit128->SetDCAToVertex2D(kTRUE);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit128->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit128->SetMaxChi2PerClusterTPC(4);
  }
  //
  // Hybrid tracks
  fESDtrackCuts_Bit768 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fESDtrackCuts_Bit768->SetName("Bit768");
  fESDtrackCuts_Bit768->SetEtaRange(-100.,100.);
  fESDtrackCuts_Bit768->SetPtRange(0.,100000.);
  fESDtrackCuts_Bit768->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts_Bit768->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts_Bit768->SetDCAToVertex2D(kTRUE);
  fESDtrackCuts_Bit768->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts_Bit768->SetMaxFractionSharedTPCClusters(0.4);
  fESDtrackCuts_Bit768->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  fESDtrackCuts_Bit768->SetRequireITSRefit(kTRUE);
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts_Bit768->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts_Bit768->SetMaxChi2PerClusterTPC(4);
  }
  //
  // for the systematic check fill all tracks and tag them with cutbit but for MC do not
  //
  fESDtrackCuts = new AliESDtrackCuts("esdTrackCuts","");
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  fESDtrackCuts->SetPtRange(0.15,1000.);
  fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fESDtrackCuts->SetMaxChi2PerClusterITS(36);
  fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);    // ?? FROM MARIAN
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetRequireITSRefit(kTRUE);
  fESDtrackCuts->SetMinNCrossedRowsTPC(70);
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0208+0.04/pt^1.01");
  fESDtrackCuts->SetMaxDCAToVertexXY(2.4);
  fESDtrackCuts->SetMaxDCAToVertexZ(3.2);
  fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDtrackCuts->SetDCAToVertex2D(kTRUE);  // fESDtrackCuts->SetDCAToVertex2D(kFALSE);    TODO
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    fESDtrackCuts->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if (fIncludeITS) {
    // require ITS pixels  -->  Reason for the empty events and structure in phi
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  }
  //
  // very loose cuts --> cuts will be tightened using the bipmap
  fESDtrackCutsLoose = new AliESDtrackCuts("esdTrackCutsLoose","");
  fESDtrackCutsLoose->SetEtaRange(-0.9,0.9);
  fESDtrackCutsLoose->SetPtRange(0.1,100000.);
  fESDtrackCutsLoose->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsLoose->SetMinNClustersTPC(50);
  fESDtrackCutsLoose->SetMinNCrossedRowsTPC(50);
  fESDtrackCutsLoose->SetMaxDCAToVertexXY(10);
  fESDtrackCutsLoose->SetMaxDCAToVertexZ(10);
  //
  // track cuts to be used for v0s
  fESDtrackCutsCleanSamp = new AliESDtrackCuts("AliESDtrackCutsV0","");
  fESDtrackCutsCleanSamp -> SetEtaRange(-1.5,1.5);
  fESDtrackCutsCleanSamp -> SetPtRange(0.,1e10);
  fESDtrackCutsCleanSamp -> SetMinNCrossedRowsTPC(80);
  fESDtrackCutsCleanSamp -> SetRequireTPCRefit(kTRUE);
  fESDtrackCutsCleanSamp -> SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsCleanSamp -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fESDtrackCutsCleanSamp -> SetMaxChi2PerClusterITS(36);
  fESDtrackCutsCleanSamp -> SetMaxFractionSharedTPCClusters(0.4);
  // ------------------------------------------------
  // V0 selection
  // ------------------------------------------------
  fESDtrackCutsV0 = new AliESDv0Cuts("AliESDCutsV0","");
  fESDtrackCutsV0 ->SetMaxDcaV0Daughters(1.0);
  // ------------------------------------------------
  //
  // Special selection for clean samples
  fV0OpenCuts   = new AliESDv0KineCuts();
  fV0StrongCuts = new AliESDv0KineCuts();
  SetSpecialV0Cuts(fV0OpenCuts);
  SetSpecialV0Cuts(fV0StrongCuts);

  fPileUpTightnessCut4 = new AliEventCuts(kFALSE);
  fPileUpTightnessCut3 = new AliEventCuts(kFALSE);
  fPileUpTightnessCut2 = new AliEventCuts(kFALSE);
  fPileUpTightnessCut1 = new AliEventCuts(kFALSE);

  fPileUpTightnessCut4->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,4);
  fPileUpTightnessCut3->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,3);
  fPileUpTightnessCut2->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,2);
  fPileUpTightnessCut1->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1);
  //
  //
  std::cout << " Info::marsland: ===================================================== " << std::endl;
  std::cout << " Info::marsland: =============== Summary of Track Cuts =============== " << std::endl;
  std::cout << " Info::marsland: ===================================================== " << std::endl;
}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::UserCreateOutputObjects()
{
  //
  // Create output histograms, trees of the analysis (called once)
  //
  Initialize();
  std::cout << " Info::marsland: ===== In the UserCreateOutputObjects ===== " << std::endl;
  // ------------  setup PIDCombined  ---------------
  fPIDCombined = new AliPIDCombined("pidCombined","");
  //
  // **********************   Input handler to get the PID object *********************
  if (!(fRunFastSimulation || fRunFastHighMomentCal)) {
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    if (!inputHandler)
    AliFatal("Input handler needed");
    else {
      fPIDResponse = inputHandler->GetPIDResponse();       // PID response object
      if (!fPIDResponse) std::cout << " Info::marsland: ======= PIDResponse object was not created ====== " << std::endl;
    }
  }
  //
  // ************************************************************************
  //   OpenFile output --> one can open several files
  // ************************************************************************
  //
  OpenFile(1);
  fTreeSRedirector = new TTreeSRedirector();
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  //
  if (fDefaultEventCuts) fEventCuts.AddQAplotsToList(fListHist); /// fList is your output TList
  //
  // ************************************************************************
  //   Efficiency matrix histograms
  // ************************************************************************
  //
  const Int_t ndim=5;
  Int_t nbins0[ndim] ={
    3,
    (Int_t) fEffMatrixCentBins.size() - 1,
    (Int_t) fEffMatrixMomBins.size() - 1,
    (Int_t) fEffMatrixEtaBins.size() - 1,
    50
  };

  std::vector<Double_t> particlesBins = {0., 1., 2., 3.};
  std::vector<Double_t> centBins = fEffMatrixCentBins;
  std::vector<Double_t> momBins = fEffMatrixMomBins;
  std::vector<Double_t> etaBins = fEffMatrixEtaBins;
  std::vector<Double_t> phiBins = {0.};
  for (Int_t i = 1; i <= 50; i++) {
    phiBins.push_back(i * 6.25 / 50.);
  }

  std::vector<std::vector<Double_t>> effMatrixBins = {
    particlesBins,
    centBins,
    momBins,
    etaBins,
    phiBins
  };

  fHistPosEffMatrixRec  =new THnF("fHistPosEffMatrixRec", "fHistPosEffMatrixRec", ndim, nbins0, effMatrixBins);
  fHistNegEffMatrixRec  =new THnF("fHistNegEffMatrixRec", "fHistNegEffMatrixRec", ndim, nbins0, effMatrixBins);
  fHistPosEffMatrixGen  =new THnF("fHistPosEffMatrixGen", "fHistPosEffMatrixGen", ndim, nbins0, effMatrixBins);
  fHistNegEffMatrixGen  =new THnF("fHistNegEffMatrixGen", "fHistNegEffMatrixGen", ndim, nbins0, effMatrixBins);
  TString axisNameEff[ndim]  = {"particle"      ,"Centrality"     ,"momentum"      ,"eta"  ,"phi"};
  TString axisTitleEff[ndim] = {"particle type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta" ,"#phi"};
  for (Int_t iEff=0; iEff<ndim;iEff++){
    fHistPosEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
    fHistNegEffMatrixRec->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixRec->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
    fHistPosEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistPosEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
    fHistNegEffMatrixGen->GetAxis(iEff)->SetName(axisNameEff[iEff]);  fHistNegEffMatrixGen->GetAxis(iEff)->SetTitle(axisTitleEff[iEff]);
  }
  fListHist->Add(fHistPosEffMatrixRec);
  fListHist->Add(fHistNegEffMatrixRec);
  fListHist->Add(fHistPosEffMatrixGen);
  fListHist->Add(fHistNegEffMatrixGen);
  //
  //
  // const Int_t lastSystSetting = *max_element(fSystSettings.begin(), fSystSettings.end());
  const Int_t ndimScan=7;
  Int_t nbinsScan[ndimScan] ={
    3,
    2,
    kCutSettingsCount,
    3,
    (Int_t) fEffMatrixCentBins.size() - 1,
    (Int_t) fEffMatrixMomBins.size() - 1,
    (Int_t) fEffMatrixEtaBins.size() - 1
  };

  std::vector<Double_t> detectorBins = {0., 1., 2., 3.};
  std::vector<Double_t> originBins = {0., 1., 2.};
  std::vector<Double_t> settingsBins = {0.};
  for (Int_t i = 1; i <= kCutSettingsCount; i++) {
    settingsBins.push_back(i);
  }

  std::vector<std::vector<Double_t>> effMatrixScanBins = {
    detectorBins,
    originBins,
    settingsBins,
    particlesBins,
    centBins,
    momBins,
    etaBins
  };

  fHistPosEffMatrixScanRec  =new THnF("fHistPosEffMatrixScanRec", "fHistPosEffMatrixScanRec", ndimScan, nbinsScan, effMatrixScanBins);
  fHistNegEffMatrixScanRec  =new THnF("fHistNegEffMatrixScanRec", "fHistNegEffMatrixScanRec", ndimScan, nbinsScan, effMatrixScanBins);
  fHistPosEffMatrixScanGen  =new THnF("fHistPosEffMatrixScanGen", "fHistPosEffMatrixScanGen", ndimScan, nbinsScan, effMatrixScanBins);
  fHistNegEffMatrixScanGen  =new THnF("fHistNegEffMatrixScanGen", "fHistNegEffMatrixScanGen", ndimScan, nbinsScan, effMatrixScanBins);

  TString axisNameEffScan[ndimScan]  = {"detector", "origin",    "systematic", "particle"      ,"Centrality"     ,"momentum"                ,"eta"  };
  TString axisTitleEffScan[ndimScan] = {"detector", "isprimary", "setting",    "particle type" ,"Centrality (%)" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta"};
  for (Int_t iEff=0; iEff<ndimScan;iEff++){
    fHistPosEffMatrixScanRec->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistPosEffMatrixScanRec->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
    fHistNegEffMatrixScanRec->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistNegEffMatrixScanRec->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
    fHistPosEffMatrixScanGen->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistPosEffMatrixScanGen->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
    fHistNegEffMatrixScanGen->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistNegEffMatrixScanGen->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
  }
  fListHist->Add(fHistPosEffMatrixScanRec);
  fListHist->Add(fHistNegEffMatrixScanRec);
  fListHist->Add(fHistPosEffMatrixScanGen);
  fListHist->Add(fHistNegEffMatrixScanGen);

  // array to store efficiency matrix projections
  // dimensions: [setting][cent][eta][isTOF][sign]
  const size_t etaDim = fEffMatrixRecPos->GetAxis(6)->GetNbins();
  const size_t centDim = fEffMatrixRecPos->GetAxis(4)->GetNbins();
  const size_t settingsDim = kCutSettingsCount;
  const size_t signDim = 2;

  fEffMatrixProjections = vector<vector<vector<vector<TH2F*>>>>(centDim);

  // resize the vectors to the appropriate dimensions
  for (size_t iCent = 0; iCent < centDim; iCent++) {
    fEffMatrixProjections[iCent] = vector<vector<vector<TH2F*>>>(settingsDim);
    for (size_t iSetting = 0; iSetting < settingsDim; iSetting++) {
      fEffMatrixProjections[iCent][iSetting] = vector<vector<TH2F*>>(signDim);
      for (size_t iSign = 0; iSign < signDim; iSign++) {
        fEffMatrixProjections[iCent][iSetting][iSign] = vector<TH2F*>(3);
      }
    }
  }
  //
  // ************************************************************************
  //   Event histograms
  // ************************************************************************
  //
  Int_t nSafeBins = (fRunFastSimulation) ? 1200 : 5; // TODO
  fHistEmptyEvent        = new TH1F("hEmptyEvent",           "control histogram to count empty events"    , 10,  0., 10.);
  fHistCentrality        = new TH1F("hCentrality",           "control histogram for centrality"           , 10,  0., 100.);
  fHistCentralityImpPar  = new TH1F("hCentralityImpPar",     "control histogram for centrality imppar"    , 10,  0., 100.);
  fHistImpParam          = new TH1F("hImpParam",             "control histogram for impact parameter"     , 200, 0., 20.);
  fHistVertex            = new TH1F("hVertex",               "control histogram for vertex Z position"    , 200, -20., 20.);
  fHistReturns           = new TH1D("hReturns",              "control histogram for reasons of returns"   , 7, 0, 7);
  fHistPileup            = new TH1D("hPileup",               "control histogram for pileup cut"           , 5, 0, 5);
  fHistPileup2D          = new TH3D("hPileup2D",             "control histogram for pileup cut"           , 100, 0, 50000, 180, 0, 4500000, 5, 0, 5);

  vector<TString> returnLabels = {"All", "Vz>12", "Vertex ok", "Cent > 90", "Empty", "Finished"};
  for (size_t i = 0; i < returnLabels.size(); i++) {
    fHistReturns->GetXaxis()->SetBinLabel(i + 1, returnLabels[i]);
  }

  fHistGenMult           = new TH1F("hGenPrMult",            "generated protons"                          , fGenprotonBins,0., 200.);
  fHistRapDistFullAccPr  = new TH2F("hRapDistFullAccPr",     "rapidity dist of protons in full acceptance"    , nSafeBins, -15, 15., 10, 0., 100.);
  fHistRapDistFullAccAPr = new TH2F("hRapDistFullAccApr",    "rapidity dist of antiprotons in full acceptance", nSafeBins, -15, 15., 10, 0., 100.);
  fJetHistptSub          = new TH1F("hjetptsub",             "control histogram for rho subtracted jetpt" , 1200, -200., 400.);

  fListHist->Add(fHistEmptyEvent);
  fListHist->Add(fHistCentrality);
  if (fMCtrue) {
    fListHist->Add(fHistCentralityImpPar);
    fListHist->Add(fHistImpParam);
    fListHist->Add(fHistGenMult);
    fListHist->Add(fHistRapDistFullAccPr);
    fListHist->Add(fHistRapDistFullAccAPr);
  }
  fListHist->Add(fHistVertex);
  fListHist->Add(fHistReturns);
  fListHist->Add(fHistPileup);
  fListHist->Add(fHistPileup2D);
  fEventInfo_CentralityEstimates  = new TVectorF(3);
  for (Int_t i=0;i<3;i++) (*fEventInfo_CentralityEstimates)[i]=-10.;
  //
  // ************************************************************************
  //   Marians counters
  // ************************************************************************
  //
  if (fEventInfo)
  {
    // vectors
    fEventInfo_PhiTPCdcarA         = new TVectorF(36);
    fEventInfo_PhiTPCdcarC         = new TVectorF(36);
    fEventInfo_CacheTrackCounters  = new TVectorF(20);
    fEventInfo_CacheTrackdEdxRatio = new TVectorF(30);
    fEventInfo_CacheTrackNcl       = new TVectorF(20);
    fEventInfo_CacheTrackChi2      = new TVectorF(20);
    fEventInfo_CacheTrackMatchEff  = new TVectorF(20);
    fEventInfo_CacheTrackTPCCountersZ = new TVectorF(8);
    for (Int_t i=0;i<8;i++) (*fEventInfo_CacheTrackTPCCountersZ)[i]=0.;
    for (Int_t i=0;i<36;i++){
      (*fEventInfo_PhiTPCdcarA)[i]=0.;
      (*fEventInfo_PhiTPCdcarC)[i]=0.;
    }
    for (Int_t i=0;i<20;i++){
      (*fEventInfo_CacheTrackCounters)[i]=0.;
      (*fEventInfo_CacheTrackdEdxRatio)[i]=0.;
      (*fEventInfo_CacheTrackNcl)[i]=0.;
      (*fEventInfo_CacheTrackChi2)[i]=0.;
      (*fEventInfo_CacheTrackMatchEff)[i]=0.;
    }
    // Hists
    fHistPhiTPCcounterA    = new TH1F("hPhiTPCcounterC",       "control histogram to count tracks on the A side in phi ", 36, 0.,18.);
    fHistPhiTPCcounterC    = new TH1F("hPhiTPCcounterA",       "control histogram to count tracks on the C side in phi ", 36, 0.,18.);
    fHistPhiTPCcounterAITS = new TH1F("hPhiTPCcounterAITS",    "control histogram to count tracks on the A side in phi ", 36, 0.,18.);
    fHistPhiTPCcounterCITS = new TH1F("hPhiTPCcounterCITS",    "control histogram to count tracks on the C side in phi ", 36, 0.,18.);
    fHistPhiITScounterA    = new TH1F("hPhiITScounterA",       "control histogram to count tracks on the A side in phi ", 36, 0.,18.);
    fHistPhiITScounterC    = new TH1F("hPhiITScounterC",       "control histogram to count tracks on the C side in phi ", 36, 0.,18.);
    fEventInfo_HisTPCVertexA = new TH1F("hisTPCZA", "hisTPCZA", 1000, -250, 250);
    fEventInfo_HisTPCVertexC = new TH1F("hisTPCZC", "hisTPCZC", 1000, -250, 250);
    fEventInfo_HisTPCVertex = new TH1F("hisTPCZ", "hisTPCZ", 1000, -250, 250);
    fEventInfo_HisTPCVertexACut = new TH1F("hisTPCZACut", "hisTPCZACut", 1000, -250, 250);
    fEventInfo_HisTPCVertexCCut = new TH1F("hisTPCZCCut", "hisTPCZCCut", 1000, -250, 250);
    fEventInfo_HisTPCVertex->SetLineColor(1);
    fEventInfo_HisTPCVertexA->SetLineColor(2);
    fEventInfo_HisTPCVertexC->SetLineColor(4);
    fEventInfo_HisTPCVertexACut->SetLineColor(3);
    fEventInfo_HisTPCVertexCCut->SetLineColor(6);
  }
  //
  // ************************************************************************
  //   Clean sample helper histograms
  // ************************************************************************
  //
  if (fV0InvMassHists && fFillArmPodTree)
  {
    fHistArmPod            = new TH2F("hArmPod",           "Armenteros-Podolanski plot"                     , 100,-1.,1., 110,0.,0.22);
    fHistInvK0s            = new TH1F("fHistInvK0s",       "control histogram for K0s invariant mass"       , 1000, 0.3,  0.70);
    fHistInvLambda         = new TH1F("fHistInvLambda",    "control histogram for lambda invariant mass"    , 1000, 1.07, 1.16);
    fHistInvAntiLambda     = new TH1F("fHistInvAntiLambda","control histogram for antilambda invariant mass", 1000, 1.07, 1.16);
    fHistInvPhoton         = new TH1F("fHistInvPhoton",    "control histogram for photon invariant mass"    , 1000, 0.,   0.05);
    fListHist->Add(fHistInvK0s);
    fListHist->Add(fHistInvLambda);
    fListHist->Add(fHistInvAntiLambda);
    fListHist->Add(fHistInvPhoton);
    fListHist->Add(fHistArmPod);
  }
  //
  // ************************************************************************
  //   Trees
  // ************************************************************************
  //
  fArmPodTree         = ((*fTreeSRedirector)<<"cleanSamp").GetTree();      // 2
  fTreeEvents         = ((*fTreeSRedirector)<<"eventInfo").GetTree();      // 3
  fTreeEventsMC       = ((*fTreeSRedirector)<<"eventInfoMC").GetTree();    // 4
  fTreeTracks         = ((*fTreeSRedirector)<<"tracks").GetTree();         // 5
  fTreeTracksMCgen    = ((*fTreeSRedirector)<<"tracksMCgen").GetTree();    // 6
  fTreeTracksMCrec    = ((*fTreeSRedirector)<<"tracksMCrec").GetTree();    // 7
  fTreeDScaled        = ((*fTreeSRedirector)<<"tracksdscaled").GetTree();  // 8
  fTreeMomentsMCfull  = ((*fTreeSRedirector)<<"momentsMCrec").GetTree();   // 9
  fTreeMomentsMCgen   = ((*fTreeSRedirector)<<"momentsMCgen").GetTree();   // 10
  fTreeDebug          = ((*fTreeSRedirector)<<"debug").GetTree();          // 11
  fTreeDebug2         = ((*fTreeSRedirector)<<"debug2").GetTree();         // 12
  fTreeDebug3         = ((*fTreeSRedirector)<<"debug3").GetTree();         // 13
  fTreeDebug4         = ((*fTreeSRedirector)<<"debug4").GetTree();         // 14
  fTreeResonances     = ((*fTreeSRedirector)<<"resonances").GetTree();     // 15
  fTreeExpecteds      = ((*fTreeSRedirector)<<"expecteds").GetTree();      // 16
  fTreeCutBased       = ((*fTreeSRedirector)<<"cutBased").GetTree();       // 17
  fTreejetsFJ         = ((*fTreeSRedirector)<<"jetsFJ").GetTree();         // 18
  fTreejetsFJBG       = ((*fTreeSRedirector)<<"jetsFJBG").GetTree();       // 19
  fTreejetsFJconst    = ((*fTreeSRedirector)<<"jetsFJconst").GetTree();    // 20
  fTreejetsFJGen      = ((*fTreeSRedirector)<<"jetsFJGen").GetTree();      // 21
  fTreejetsFJBGGen    = ((*fTreeSRedirector)<<"jetsFJBGGen").GetTree();    // 22
  fTreejetsFJconstGen = ((*fTreeSRedirector)<<"jetsFJconstGen").GetTree(); // 23
  //
  // ************************************************************************
  //   Send output objects to container (remember to also add the outputs using DefineOutput above and in the AddTask)
  // ************************************************************************
  //
  PostData(1,  fListHist);
  PostData(2,  fArmPodTree);
  PostData(3,  fTreeEvents);
  PostData(4,  fTreeEventsMC);
  PostData(5,  fTreeTracks);
  PostData(6,  fTreeTracksMCgen);
  PostData(7,  fTreeTracksMCrec);
  PostData(8,  fTreeDScaled);
  PostData(9,  fTreeMomentsMCfull);
  PostData(10, fTreeMomentsMCgen);
  PostData(11, fTreeDebug);
  PostData(12, fTreeDebug2);
  PostData(13, fTreeDebug3);
  PostData(14, fTreeDebug4);
  PostData(15, fTreeResonances);
  PostData(16, fTreeExpecteds);
  PostData(17, fTreeCutBased);
  PostData(18, fTreejetsFJ);
  PostData(19, fTreejetsFJBG);
  PostData(20, fTreejetsFJconst);
  PostData(21, fTreejetsFJGen);
  PostData(22, fTreejetsFJBGGen);
  PostData(23, fTreejetsFJconstGen);

  fEventCuts.SetManualMode();

  std::cout << " Info::marsland: ===== Out of UserCreateOutputObjects ===== " << std::endl;

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::UserExec(Option_t *)
{
  //
  // main event loop
  //
  if (fUseCouts) std::cout << "event = " << fEventCountInFile << " Info::marsland: ===== In the UserExec ===== " << std::endl;
  //
  // Check Monte Carlo information and other access first:
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) fMCtrue = kFALSE;
  //
  //  get the filename
  TTree *chain = (TChain*)GetInputData(0);
  if(!chain) { Printf(" Error::marsland: Could not receive input chain"); return; }
  TString tmpChunkname = fChunkName;
  TObjString fileName(chain->GetCurrentFile()->GetName());
  fChunkName = (TString)fileName.GetString();
  if (tmpChunkname != fChunkName) {
    fEventCountInFile = 0;
    fChunkIDinJob++;
    std::cout << " chunk ID = " << fChunkIDinJob << " Info::marsland: ===== Current chunk name is ===== " << fChunkName << std::endl;
  }
  fEventCountInFile++;
  //
  // ======================================================================
  // ========================== See if MC or Real =========================
  // ======================================================================
  //
  if (eventHandler) fMCEvent = eventHandler->MCEvent();
  AliGenEventHeader* genHeader = 0x0;
  TString genheaderName;
  if (fMCEvent){
    genHeader = fMCEvent->GenEventHeader();
    genheaderName = genHeader->GetName();
    if(!genHeader){ Printf(" Error::marsland: Event generator header not available!!!\n"); return; }
  }
  //
  // If ESDs exist get some event variables
  //
  fCentrality = -5;
  fCentImpBin =-10.;
  AliCentrality    *esdCentrality = 0x0;
  AliMultSelection *MultSelection = 0x0;
  AliMultSelectionTask *MultSelectionTask = 0x0;
  ULong64_t gid=0;
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (fESD)
  {
    //
    // Init magnetic filed for golden chi2 cut
    fESD->InitMagneticField();
    //
    // event selection
    fPileUpBit=0;
    if (fDefaultEventCuts && fESD){

      if ( (fPassIndex==3 || fPassIndex==2) && fYear>2013){
        //
        // pileup bit: 0bxxxx, where the leftmost bit is the tightest and the rightmost is the loosest cut
        // OOB pileup cut (for Pb-Pb) based on ITS and TPC clusters: 0-> no cut; 1-> default cut (remove all OOB pileup); 2-> looser cut; 3-> even more looser cut; 4-> very loose cut
        vector<AliEventCuts*> pileUpCuts = {fPileUpTightnessCut1, fPileUpTightnessCut2, fPileUpTightnessCut3, fPileUpTightnessCut4};
        fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,0); // do not apply any pile cut

        Double_t nTPCClusters = fESD->GetNumberOfTPCClusters();
        Double_t nITSClusters = 0;
        AliVMultiplicity *multiObj = fESD->GetMultiplicity();
        for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);

        if (fEventCuts.AcceptEvent(fESD)) {
          fHistPileup->Fill(0);
          fHistPileup2D->Fill(nITSClusters, nTPCClusters, 0.);
        }

        for (size_t iPileupCut = 0; iPileupCut < pileUpCuts.size(); iPileupCut++) {
          if (pileUpCuts[iPileupCut]->AcceptEvent(fESD)) {
            fPileUpBit |= 1 << iPileupCut;
            Double_t binIndex = pileUpCuts.size() - iPileupCut;
            fHistPileup->Fill(binIndex);
            fHistPileup2D->Fill(nITSClusters, nTPCClusters, binIndex);
          }
        }
        // if (!fEventCuts.AcceptEvent(fESD)) {cout<< "pileup event " << endl; fHistEmptyEvent->Fill(0); return;} // TODO
      }
    }
    //
    // calculate TPC mult
    fTPCMult = 0;
    for (Int_t itrack=0;itrack<fESD->GetNumberOfTracks();++itrack){
      AliESDtrack *track = fESD->GetTrack(itrack);
      if (track->IsOn(AliESDtrack::kTPCin)) fTPCMult++;
    }
    //
    //
    esdCentrality = fESD->GetCentrality();
    MultSelection = (AliMultSelection*) fESD-> FindListObject("MultSelection");
    fRunNo = fESD->GetRunNumber();
    //
    // if the run number is specified, fill expecteds only for that run. otherwise fill for any run
    if (fRunNumberForExpecteds > 0) fFillExpecteds = (fRunNo == fRunNumberForExpecteds);
    else fFillExpecteds = kTRUE;
    //
    static Int_t timeStampCache = -1;
    if (!fMCtrue) {
      fTimeStamp = fESD->GetTimeStampCTPBCCorr();
      const char *ocdb;
      if(!fRunOnGrid) ocdb = Form("local:///cvmfs/alice.cern.ch/calibration/data/%d/OCDB",fYear);
      else ocdb = "raw://";
      //
      // retrieve interaction rate
      if (timeStampCache!=Int_t(fTimeStamp)) timeStampCache=Int_t(fTimeStamp);
      if (!gGrid && timeStampCache>0) {
        AliInfo("Trying to connect to AliEn ...");
        TGrid::Connect("alien://");
      } else {
        fEventInfo_LumiGraph = (TGraph*)AliLumiTools::GetLumiFromCTP(fRunNo,ocdb);
        fIntRate   = fEventInfo_LumiGraph->Eval(fTimeStamp); delete fEventInfo_LumiGraph;
      }
    }
    fEventMult = fESD->GetNumberOfTracks();
    fBField    = fESD->GetMagneticField();
    fEvent     = fESD->GetEventNumberInFile();
    fBeamType  = fESD->GetBeamType();
    //
    // Global id for the event --> which is made with Hashing
    //
    ULong64_t orbitID      = (ULong64_t)fESD->GetOrbitNumber();
    ULong64_t bunchCrossID = (ULong64_t)fESD->GetBunchCrossNumber();
    ULong64_t periodID     = (ULong64_t)fESD->GetPeriodNumber();
    gid = ((periodID << 36) | (orbitID << 12) | bunchCrossID);
    fEventGID = gid;    // uniqe event id for real data
    // fTimeStamp = fESD->GetTimeStamp();
    // fEventGID  = TMath::Abs(Int_t(TString::Hash(&gid,sizeof(Int_t))));    // uniqe event id for real data
  }
  //
  // Get rid of "E-AliESDpid::GetTPCsignalTunedOnData: Tune On Data requested, but MC event not set. Call SetCurrentMCEvent before!" errors
  if (!fPIDResponse && !(fRunFastSimulation || fRunFastHighMomentCal)) fPIDResponse = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  if (fMCEvent && fPIDResponse) fPIDResponse->SetCurrentMCEvent(fMCEvent);
  //
  if(fMCtrue)
  {
    //
    // Add different generator particles to PDG Data Base to avoid problems when reading MC generator particles
    AliPDG::AddParticlesToPdgDataBase();
    //
    // ========================== MC =========================
    //
    fMCStack = fMCEvent->Stack();
    if (!fMCStack) { Printf(" Error::marsland: No MC stack available !!!\n"); return;}
    //
    if (MultSelection) {
      fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
      if (fUseCouts)  std::cout << " Info::marsland: Centrality is taken from MultSelection " << fCentrality << std::endl;
    } else if (esdCentrality) {
      fCentrality = esdCentrality->GetCentralityPercentile("V0M");
      if (fUseCouts)  std::cout << " Info::marsland: Centrality is taken from esdCentrality " << fCentrality << std::endl;
    }
    //
    // AliCollisionGeometry* colGeometry = dynamic_cast<AliCollisionGeometry*>(genHeader);
    // std::cout << " aaaa  " <<  ((AliGenEposEventHeader*) genHeader)->ImpactParameter() << std::endl;
    // std::cout << " bbbb  " <<  colGeometry->ImpactParameter() << std::endl;
    //
    // OLD --> impact parameters to use: 0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5
    // OLD --> corresponding Centrality:  0     5    10    20    30     40     50    60      70    80
    // Double_t impParArr[10] = {0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5};
    //
    // source for the impact parameters --> https://arxiv.org/pdf/1710.07098.pdf
    Double_t impParArr[12] = {0.0, 3.49, 4.93, 6.98, 8.55, 9.87, 11.0, 12.1, 13.1, 14.0, 14.9, 20.0};
    AliGenHijingEventHeader *lHIJINGHeader = 0x0;  // event header for HIJING
    AliGenHepMCEventHeader *lHepMCHeader = 0x0;    // event header for EPOS
    if ( fMCEvent ){
      //
      // If EPOS based MC
      if ( genheaderName.Contains("EPOSLHC") ){
        if (genHeader->InheritsFrom(AliGenHepMCEventHeader::Class())) lHepMCHeader = (AliGenHepMCEventHeader*)genHeader;
        if (lHepMCHeader ){
          fNHardScatters = lHepMCHeader->Ncoll_hard(); // Number of hard scatterings
          fNProjectileParticipants = lHepMCHeader->Npart_proj(); // Number of projectile participants
          fNTargetParticipants     = lHepMCHeader->Npart_targ(); // Number of target participants
          fNNColl   = lHepMCHeader->Ncoll(); // Number of NN (nucleon-nucleon) collisions
          fNNwColl  = lHepMCHeader->N_Nwounded_collisions(); // Number of N-Nwounded collisions
          fNwNColl  = lHepMCHeader->Nwounded_N_collisions(); // Number of Nwounded-N collisons
          fNwNwColl = lHepMCHeader->Nwounded_Nwounded_collisions();// Number of Nwounded-Nwounded collisions
          fMCImpactParameter = lHepMCHeader->impact_parameter();
          if (fUseCouts)  std::cout << " Info::marsland: EPOS: Centrality is taken from ImpactParameter = " << fMCImpactParameter << "  "  << ((AliGenEposEventHeader*) genHeader)->GetName() << std::endl;
        }
        //
        // If HIJING based MC
      } else if (!TMath::IsNaN(((AliGenHijingEventHeader*) genHeader)->ImpactParameter()) ){
        lHIJINGHeader = (AliGenHijingEventHeader*) genHeader;
        if (lHIJINGHeader ){
          fNHardScatters = lHIJINGHeader->HardScatters();
          fNProjectileParticipants = lHIJINGHeader->ProjectileParticipants();
          fNTargetParticipants     = lHIJINGHeader->TargetParticipants();
          fNNColl   = lHIJINGHeader->NN();
          fNNwColl  = lHIJINGHeader->NNw();
          fNwNColl  = lHIJINGHeader->NwN();
          fNwNwColl = lHIJINGHeader->NwNw();
          fMCImpactParameter = lHIJINGHeader->ImpactParameter();
          if (fUseCouts)  std::cout << " Info::marsland: HIJING: Centrality is taken from ImpactParameter = " << fMCImpactParameter << "  "  << ((AliGenHijingEventHeader*) genHeader)->GetName() << std::endl;
        }
      }
      fHistImpParam->Fill(fMCImpactParameter);
      if (fMCImpactParameter>=impParArr[0]  && fMCImpactParameter<impParArr[1])  fCentImpBin=2.5;
      if (fMCImpactParameter>=impParArr[1]  && fMCImpactParameter<impParArr[2])  fCentImpBin=7.5;
      if (fMCImpactParameter>=impParArr[2]  && fMCImpactParameter<impParArr[3])  fCentImpBin=15.;
      if (fMCImpactParameter>=impParArr[3]  && fMCImpactParameter<impParArr[4])  fCentImpBin=25.;
      if (fMCImpactParameter>=impParArr[4]  && fMCImpactParameter<impParArr[5])  fCentImpBin=35.;
      if (fMCImpactParameter>=impParArr[5]  && fMCImpactParameter<impParArr[6])  fCentImpBin=45.;
      if (fMCImpactParameter>=impParArr[6]  && fMCImpactParameter<impParArr[7])  fCentImpBin=55.;
      if (fMCImpactParameter>=impParArr[7]  && fMCImpactParameter<impParArr[8])  fCentImpBin=65.;
      if (fMCImpactParameter>=impParArr[8]  && fMCImpactParameter<impParArr[9])  fCentImpBin=75.;
      if (fMCImpactParameter>=impParArr[9]  && fMCImpactParameter<impParArr[10]) fCentImpBin=85.;
      if (fMCImpactParameter>=impParArr[10] && fMCImpactParameter<impParArr[11]) fCentImpBin=95.;
      if (fMCImpactParameter<=impParArr[0]  && fMCImpactParameter>impParArr[11]) fCentImpBin=-10.;
    }
    if (fCentrality<0) fCentrality=fCentImpBin;
    //
    // Use file name in Hashing to create unique event ID
    fEventGIDMC  = TMath::Abs(Int_t(TString::Hash(&fEventCountInFile,sizeof(Int_t))));    // uniqe event id for real data
    // fEventGIDMC += TMath::Abs(Int_t(fCentrality)+fEventCountInFile+(1000*TMath::Abs(fMCImpactParameter)));  // ????
    fEventGID    = fEventGIDMC;
    if (fUseCouts) {
      std::cout << " Info::marsland: ========================================================================================== " << std::endl;
      std::cout << " Info::marsland: " << fEventCountInFile << " ----- eventIDMC = " << fEventGIDMC << "   " << fChunkName << std::endl;
      std::cout << " Info::marsland: Centrality = " << fCentrality << " ----- Impact Param = " << fMCImpactParameter << " fCentralityImp = " << fCentImpBin << std::endl;
      std::cout << " Info::marsland: ========================================================================================== " << std::endl;
    }
    //
  }
  // the event exists, we start counting from here
  fHistReturns->Fill(0);
  //
  if (!(fRunFastSimulation || fRunFastHighMomentCal))
  {
    //
    // ========================== Real =========================
    //
    if (!fESD)          { Printf(" Error::marsland: fESD not available"); return; }
    if (!fESDtrackCuts) { Printf(" Error::marsland: fESDtrackCuts not available"); return; }
    //
    // ------------------------------------------------
    // ------- monitor vertex position =---------------
    // ------------------------------------------------
    //
    Bool_t isVertexOk = kTRUE;
    fVertex = fESD->GetPrimaryVertexTracks();
    const AliESDVertex *vertexSPD = fESD->GetPrimaryVertexTracks();
    const AliESDVertex *vertexTPC = fESD->GetPrimaryVertexTracks();
    if( fVertex->GetNContributors()<1) isVertexOk = kFALSE;
    if( fVertex->GetNContributors()>1) {
      vertexSPD = fESD->GetPrimaryVertexSPD();    // SPD vertex
      vertexTPC = fESD->GetPrimaryVertexTPC();    // TPC vertex
      fTPCvZ = vertexTPC->GetZ();
      fSPDvZ = vertexSPD->GetZ();
      fVz    = fVertex->GetZ();
      TString vertexType = fVertex->GetTitle();    // ??? Put condition Abs(vertex-vertexTPC) as a bool_t into ttree
      if ( vertexType.Contains("vertexer: Z") && (fVertex->GetDispersion() > 0.04 || fVertex->GetZRes() > 0.25) ) isVertexOk = kFALSE; // TODO
    }
    fMultiplicity    = fVertex->GetNContributors();    // fMultiplicity = fESD -> GetNumberOfTracks();
    fNContributors   = fVertex->GetNContributors();
    fMultiplicityMC  = fMultiplicity;
    //
    // ------------------------------------------------
    // ------- event vertex cut along Z ---------------
    // ------------------------------------------------
    //
    // if (fMCtrue && TMath::Abs(fVz) > 15) return;   // For MC put fixed cut
    // if (fDefaultTrackCuts && (TMath::Abs(fVz)>7 || TMath::Abs(fVz)<0.15) ) return;
    // else if (TMath::Abs(fVz)>15) return;
    if (TMath::Abs(fVz)>12) {
      fHistReturns->Fill(1);
      return;
    }
    //
    if (fVertex && isVertexOk) {
      fHistVertex->Fill(fVz);
    } else {
      fHistReturns->Fill(2);
      return;
    }
    //
    // ------------------------------------------------
    // ---------- Centrality definition ---------------
    // ------------------------------------------------
    //
    Int_t nEst = sizeof(fEventInfo_centEstStr)/sizeof(char*);
    fEventInfo_CentralityEstimates->Zero();  // matchEff counter
    if (fBeamType.CompareTo("A-A") == 0) { // PbPb
      if (MultSelection) {
        if(fSystCentEstimatetor == -1) fCentrality = MultSelection->GetMultiplicityPercentile("TRK");
        if(fSystCentEstimatetor ==  0) fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
        if(fSystCentEstimatetor ==  1) fCentrality = MultSelection->GetMultiplicityPercentile("CL1");
        for (Int_t i=0;i<nEst;i++) (*fEventInfo_CentralityEstimates)[i]=MultSelection->GetMultiplicityPercentile(fEventInfo_centEstStr[i]);
      } else if (esdCentrality) {
        if(fSystCentEstimatetor == -1) fCentrality = esdCentrality->GetCentralityPercentile("TRK");
        if(fSystCentEstimatetor ==  0) fCentrality = esdCentrality->GetCentralityPercentile("V0M");
        if(fSystCentEstimatetor ==  1) fCentrality = esdCentrality->GetCentralityPercentile("CL1");
        for (Int_t i=0;i<nEst;i++) (*fEventInfo_CentralityEstimates)[i]=esdCentrality->GetCentralityPercentile(fEventInfo_centEstStr[i]);
      } else {
        std::cout << " Info::marsland: Error: There is no cent info " << std::endl;
      }
    }
    //
    if (fUseCouts) {
      std::cout << " Info::marsland: =============================================================================================== " << std::endl;
      std::cout << " Info::marsland: Event counter = " << fEventCountInFile << " - cent =  " << fCentrality << " = gid = " << gid << " = fEventGID = " << fEventGID << " file: " << fChunkName << "   beamtype = " << fBeamType << std::endl;
      std::cout << " Info::marsland: =============================================================================================== " << std::endl;
    }
  }
  if (fCentrality > 90.) {
    fHistReturns->Fill(3);
    return;
  }
  fHistCentrality->Fill(fCentrality);  // count events after physics and vertex selection
  //
  // if (!fESDtool) {
  //   fESDtool = new AliESDtools();
  //   fESDtool->SetStreamer(fTreeSRedirector);
  // }
  // fESDtool->Init(NULL,fESD);
  // fESDtool->CalculateEventVariables();
  // fESDtool->SetMCEvent(fMCEvent);
  // fESDtool->DumpEventVariables();
  //
  // ======================================================================
  //   Filling part
  // ======================================================================
  //
  // Real Data Analysis
  //
  if (!fMCtrue && fESD){
    if (CountEmptyEvents(2,0)<1) {
      fHistReturns->Fill(4);
      return;
    }
    if (fUseCouts)  std::cout << " Info::marsland: (Real Data Analysis) Filling = " << fEventCountInFile << std::endl;
    CalculateEventInfo(); CreateEventInfoTree();
    if (fFillArmPodTree) FillCleanSamples();
    if (fTaskSelection==0) {FillTPCdEdxReal(); CalculateMoments_CutBasedMethod(); FindJetsFJ();} // bot jet and net-p
    if (fTaskSelection==1) {FindJetsFJ();} // only jet
    if (fTaskSelection==2) {FillTPCdEdxReal(); CalculateMoments_CutBasedMethod();} // only net-p
    if (fUseCouts)  std::cout << " Info::marsland: (Real Data Analysis) End of Filling part = " << fEventCountInFile << std::endl;

    fHistReturns->Fill(5);
    return;
  }
  //
  // Full MC analysis
  //
  if (fMCtrue && fESD && !fRunFastSimulation){
    if (CountEmptyEvents(2,0)<1) {
      fHistReturns->Fill(4);
      return;
    }
    if (fFillEffMatrix) {
      if (fUseCouts)  std::cout << " Info::marsland: Fill only EffMatrix, event and downscaled trees = " << fEventCountInFile << std::endl;
      FillEffMatrix();
      fHistReturns->Fill(5);
      return;
    } else {
      if (fUseCouts) std::cout << " Info::marsland: (full MC analysis) End of Filling part = " << fEventCountInFile << std::endl;
      FillEffMatrix();
      CalculateEventInfo(); 
      CreateEventInfoTree();
      FillEventInfoMC(); 
      FillTreeMC();
      if (fFillArmPodTree && fFillDebug) FillCleanSamples();
      if (fTaskSelection==0) {FillMCFull_NetParticles(); CalculateMoments_CutBasedMethod(); FindJetsFJ(); FindJetsFJGen();} // bot jet and net-p
      if (fTaskSelection==1) {FindJetsFJ(); FindJetsFJGen();} // only jet
      if (fTaskSelection==2) {FillMCFull_NetParticles(); CalculateMoments_CutBasedMethod();} // only net-p
      fHistReturns->Fill(5);
      return;
    }
  }
  //
  if (fMCtrue && fFillHigherMomentsMCclosure){
    if (CountEmptyEvents(2,0)<1) {
      fHistReturns->Fill(4);
      return;
    }
    if (fUseCouts)  std::cout << " Info::marsland: (MCclosure Higher Moments) End of Filling part = " << fEventCountInFile << std::endl;
    MCclosureHigherMoments();
    fHistReturns->Fill(5);
    return;
  }
  //
  // FastGen analysis
  //
  if (fRunFastSimulation) {
    if (CountEmptyEvents(2,-1)<1) {
      fHistReturns->Fill(4);
      return;
    }
    if (fUseCouts)  std::cout << " Info::marsland: FastGen Filling = " << fEventCountInFile << std::endl;
    FillEventInfoMC();
    if (fTaskSelection==0) {FindJetsFJGen(); FastGen_NetParticles();} // bot jet and net-p
    if (fTaskSelection==1) {FindJetsFJGen(); } // only jet
    if (fTaskSelection==2) {FastGen_NetParticles();} // only net-p
    fHistReturns->Fill(5);
    return;
  }


}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::FillTPCdEdxReal()
{
  //
  // Fill dEdx information for the TPC and also clean kaon and protons
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillTPCdEdxReal ===== " << std::endl;
  // --------------------------------------------------------------
  // Get the event
  AliVEvent *event=InputEvent();
  //
  // --------------------------------------------------------------
  //  Main track loop
  // --------------------------------------------------------------
  //
  Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
  const AliMultiplicity *multObj = fESD->GetMultiplicity();
  Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();
  for (Int_t itrack=0;itrack<event->GetNumberOfTracks();++itrack) {   // Track loop

    fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
    fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);
    Int_t label = track->GetLabel();
    if (!track) continue;
    //
    // --------------------------------------------------------------
    //      Get relevant track info and set cut bits
    // --------------------------------------------------------------
    //
    Bool_t fBit96_base   = fESDtrackCuts_Bit96->AcceptTrack(track);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(track);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(track);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(track);
    if (!track->GetInnerParam()) continue;               // Ask if track is in the TPC
    if (!fESDtrackCutsLoose->AcceptTrack(track))  continue;    // Loose cuts
    if (!(track->GetTPCsignalN()>0)) continue;
    //
    // Get the track variables
    Double_t closestPar[3];
    GetExpecteds(track,closestPar);
    SetCutBitsAndSomeTrackVariables(track,0);
    Int_t tpcNcls = track->GetTPCncls();
    //
    // --------------------------------------------------------------
    //  Some print out
    // --------------------------------------------------------------
    //
    // Tree for the all cut variables
    if (fUseCouts && fEventCountInFile==5) {
      std::cout << " Info::marsland: CutBinMap --> " <<fTrackTPCCrossedRows << " " << fTrackChi2TPC << " " <<  fTrackNewITScut  << std::endl;
      PrintNumInBinary(fTrackCutBits);
    }
    //
    // --------------------------------------------------------------
    //   Fill the trees
    // --------------------------------------------------------------
    //
    Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
    Int_t nITSClusters = 0;
    AliVMultiplicity *multiObj = fESD->GetMultiplicity();
    for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
    //
    // different dca cuts
    // TMath::Abs(fTrackDCAxy)< 0.3
    Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPt,1.1);    // 10h tuned loose cut
    Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPt,1.01);    // 10h tuned loose cut
    Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPt,1.01);  // 10h tuned loose cut

    UShort_t tpcFindableCls = track->GetTPCNclsF();
    UShort_t tpcSharedCls = track->GetTPCnclsS();
    //
    Double_t tofSignalTunedOnData = track->GetTOFsignalTunedOnData();
    Double_t length = track->GetIntegratedLength();
    Double_t tofSignal = track->GetTOFsignal();
    Double_t beta = -.05;
    if((length > 0) && (tofSignal > 0)) beta = length / 2.99792458e-2 / tofSignal;
    //
    // Get systematic settings
    TVectorF settings(kCutSettingsCount);
    for (Int_t i = 0; i < kCutSettingsCount;i++) {
      settings[i] = (Float_t) GetSystematicClassIndex(fTrackCutBits, i);
    }

    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"tracks"<<
    //
    "gid=" << fEventGID << // global event ID
    "cutBit=" << fTrackCutBits << // Systematic Cuts
    "dEdx=" << fTPCSignal <<  // dEdx of the track
    "beta=" << beta << // TOF beta  --> beta = GetIntegratedLength / 2.99792458e-2 / tofSignal;
    "sign=" << fSign <<  // charge
    "ptot=" << fPtot <<  // TPC momentum
    "p=" << fPVertex <<  // momentum at vertex
    "pT=" << fPt <<  // transverse momentum
    "eta=" << fEta <<  // eta
    "phi=" << fPhi <<  // phi
    "cent=" << fCentrality; // centrality from V0
    if (fFillDebug){
      (*fTreeSRedirector)<<"tracks"<<
      "settings.=" << &settings << // Systematic settings
      "run=" << fRunNo << // run Number
      "bField=" << fBField << // magnetic filed
      "intrate=" << fIntRate <<  // interaction rate
      "pileupbit=" << fPileUpBit << // flag for pileup selection
      "primMult=" << fNContributors << //  #prim tracks
      "tpcClMult=" << tpcClusterMultiplicity << // TPC cluster multiplicity
      "tpcmult=" << fTPCMult << //  TPC track multiplicity
      "tpcclmult=" << nTPCClusters << // TPC cluster multiplicity
      "itsmult=" << itsNumberOfTracklets << // ITS multiplicity
      "itsclmult=" << nITSClusters << // ITS cluster multiplicity
      "ncltpc=" << fNcl <<  // number of clusters
      //
      // track level
      "label=" << label <<
      "eventtime=" << fTimeStamp << // event timeStamp
      "defCut=" << fDefaultCuts << // default cuts tuned by hand
      "bit96=" << fBit96_base << // tight cuts of 2011 tuned data
      "bit128=" << fBit128 << // TPC only tracks cuts
      "bit768=" << fBit768 << // Hybrid track cuts
      "pixCut=" << ifDCAcutIfNoITSPixel << // cut: apply a DCAcut If No ITS Pixel
      "dcabase=" << dcaBaseCut << // dcaxy base cut
      "dca10h=" << dca10h << // dcaxy cut tuned to 10h
      "dca11h=" << dca11h << // dcaxy cut tuned to 11h
      "tpcFindableCls=" << tpcFindableCls << // number of findable clusters
      "tpcSharedCls=" << tpcSharedCls << // number of shared clusters
      "tpcSignalN=" << fTrackTPCSignalN << // number of cl used in dEdx
      "lengthInActiveZone=" << fTrackLengthInActiveZone << // track length in active zone
      "dcaxy=" << fTrackDCAxy << // dcaxy
      "dcaz=" << fTrackDCAz << // dcaz
      "cRows=" << fTrackTPCCrossedRows << // crossed Rows in TPC
      "chi2tpc=" << fTrackChi2TPC << // TPC chi2
      "missCl=" << fMissingCl << // fraction of missing clusters
      //
      "eltpcpid=" << fNSigmasElTPC << // nsigma TPC for electrons
      "pitpcpid=" << fNSigmasPiTPC << // nsigma TPC for pions
      "katpcpid=" << fNSigmasKaTPC << // nsigma TPC for kaons
      "prtpcpid=" << fNSigmasPrTPC << // nsigma TPC for protons
      "detpcpid=" << fNSigmasPrTOF << // nsigma TPC for deuterons
      //
      "eltofpid=" << fNSigmasElTOF << // nsigma TOF for electrons
      "pitofpid=" << fNSigmasPiTOF << // nsigma TOF for pions
      "katofpid=" << fNSigmasKaTOF << // nsigma TOF for kaons
      "prtofpid=" << fNSigmasPrTOF << // nsigma TOF for protons
      "detofpid=" << fNSigmasDeTOF << // nsigma TOF for deuterons
      //
      "tofSignal=" << tofSignal << // tof signal --> GetTOFsignal
      "tofSignalTOD=" << tofSignalTunedOnData << // tof signal tuneondata --> GetTOFsignalTunedOnData
      "prtofpid=" << fNSigmasPrTOF ;
    }
    (*fTreeSRedirector)<<"tracks"<<"\n";

    //
    // --------------------------------------------------------------
    //  Fill the THnSparseF for the Expected values form PID response
    // --------------------------------------------------------------
    //
    // define acceptance of interest
    Bool_t etaAcc  = (fEta >=fEtaMin        && fEta<=fEtaMax);
    Bool_t momAcc  = (fPVertex>=fMomDown    && fPVertex<=fMomUp);
    Bool_t dEdxAcc = (fTPCSignal>=fDEdxDown && fTPCSignal<=fDEdxUp);
    Bool_t fAcceptance = (etaAcc && momAcc && dEdxAcc);
    Bool_t nSigmasElTPCCut = (TMath::Abs(fNSigmasElTPC)<2);
    Bool_t nSigmasPiTPCCut = (TMath::Abs(fNSigmasPiTPC)<2);
    Bool_t nSigmasKaTPCCut = (TMath::Abs(fNSigmasKaTPC)<2);
    Bool_t nSigmasPrTPCCut = (TMath::Abs(fNSigmasPrTPC)<2);
    Bool_t nSigmasDeTPCCut = (TMath::Abs(fNSigmasDeTPC)<2);
    Bool_t nSigmaTPCall = (nSigmasElTPCCut || nSigmasPiTPCCut || nSigmasKaTPCCut || nSigmasPrTPCCut || nSigmasDeTPCCut);
    Bool_t ndEdxTPCall  = (fDEdxEl>20 || fDEdxPi>20 || fDEdxKa>20 || fDEdxPr>20 || fDEdxDe>20);

    if(fAcceptance && !fMCtrue){
      if (fFillExpecteds && fEvent < 5 && (fNSigmasPiTPC >= 3 || (fNSigmasPiTPC < 3 && fRandom.Rndm() < 0.001))) {
        Double_t sign = static_cast<Double_t>(fSign);
        if(!fTreeSRedirector) return;
        (*fTreeSRedirector)<<"expecteds"<<
        "cent=" << fCentrality   <<
        "sign=" << sign <<
        "ptot=" << fPtot <<
        "eta=" << fEta <<
        "phi=" << fPhi <<
        "dEdxEl=" << fDEdxEl <<
        "dEdxPi=" << fDEdxPi <<
        "dEdxKa=" << fDEdxKa <<
        "dEdxPr=" << fDEdxPr <<
        "dEdxDe=" << fDEdxDe <<
        "sigmaEl=" << fSigmaEl <<
        "sigmaPi=" << fSigmaPi <<
        "sigmaKa=" << fSigmaKa <<
        "sigmaPr=" << fSigmaPr <<
        "sigmaDe=" << fSigmaDe <<
        "dEdx=" << fTPCSignal <<
        //
        "\n";
      }
    }
  } // end of track loop
}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::CalculateMoments_CutBasedMethod()
{

  //
  // Calculate moments of net-protons with TPC and TOF dE/dx cuts
  //
  // Assign subsample index
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  // Int_t runNumber = (Int_t)fESD->GetRunNumber();
  if (fUseCouts) std::cout << " Info::marsland: ===== In the CalculateMoments_CutBasedMethod ===== " << std::endl;
  //
  // check if the event is full
  Int_t nStackTracks = fESD->GetNumberOfTracks();
  if (nStackTracks <= 1) return;

  // get centrality index
  Int_t centIndex = -1;
  for (Int_t iCent = 0; iCent < fNCentBinsMC; iCent++) {
    if (fCentrality >= fxCentBins[iCent] && fCentrality < fxCentBins[iCent+1]) {
      centIndex = iCent;
      break;
    }
  }
  if (fCollisionType == 1) centIndex = 0;

  const size_t settingsDim = fSystSettings.size();
  const size_t etaDim      = fetaDownArr.size();
  const size_t momentumDim = fpDownArr.size();
  const size_t signDim     = 2;

  const size_t nMoments = 14 + 1; // moments to fourth order + p-pbar

  size_t orderR = 4;
  size_t orderS = 4;

  UInt_t recProtonCounter[settingsDim][etaDim][momentumDim][signDim];

  // store qs up to the given order
  Double_t arrQ[settingsDim][etaDim][momentumDim][orderR][orderS];

  for (size_t iSetting = 0; iSetting < settingsDim; iSetting++) {
    for (size_t iEta = 0; iEta < etaDim; iEta++) {
      for (size_t iMomentum = 0; iMomentum < momentumDim; iMomentum++) {
        for (size_t iSign = 0; iSign < signDim; iSign++) {
          recProtonCounter[iSetting][iEta][iMomentum][iSign] = 0;
        }
        for (size_t iOrderR = 0; iOrderR < orderR; iOrderR++) {
          for (size_t iOrderS = 0; iOrderS < orderS; iOrderS++) {
            arrQ[iSetting][iEta][iMomentum][iOrderR][iOrderS] = 0.;
          }
        }
      }
    }
  }

  //
  Bool_t bCutReference           = (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15);
  Bool_t bEventVertexZLarge      = (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1);
  UInt_t counterTracksRec = 0;
  //
  // -----------------------------------------------------------------------------------------
  // ----------------------------  Real Data with PID cuts  ----------------------------------
  // -----------------------------------------------------------------------------------------
  //
  for (Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++) {
    // track loop
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *trackReal = fESD->GetTrack(irectrack);
    if (trackReal==NULL) continue;
    // apply detector cuts
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
    if (!trackReal-> GetInnerParam()) continue;
    if (!fESDtrackCutsLoose->AcceptTrack(trackReal))  continue;    // Loose Cuts
    if (!(trackReal->GetTPCsignalN()>0)) continue;
    if (!ifDCAcutIfNoITSPixel) continue;
    //
    // Get the cut bit information apply track cuts
    SetCutBitsAndSomeTrackVariables(trackReal,0);
    //
    // acceptance cuts
    Double_t ptotCut = 0.;
    if (fUsePtCut == 0) ptotCut = fPtot;
    if (fUsePtCut == 1) ptotCut = fPVertex;
    if (fUsePtCut == 2) ptotCut = fPt;

    //
    // Acceptance scan
    for (size_t ieta = 0; ieta < etaDim; ieta++) {
      Bool_t etaAcc  = (fEta >= fetaDownArr[ieta] && fEta < fetaUpArr[ieta]);
      Bool_t etaAccMaxWindow = (fEta >= fetaDownArr[0]  && fEta <= fetaUpArr[etaDim - 1]);
      if (!etaAcc) continue;

      for (size_t imom = 0; imom < momentumDim; imom++) {
        Bool_t momAcc  = (ptotCut >= fpDownArr[imom]  && ptotCut < fpUpArr[imom]);
        Bool_t momAccMaxWindow = (ptotCut >= fpDownArr[0] && ptotCut <= fpUpArr[momentumDim - 1]);
        //
        // count first moments for given Centrality and momentum window
        if (etaAccMaxWindow && momAccMaxWindow) counterTracksRec++;

        if (!momAcc) continue;
        for (size_t iset = 0; iset < settingsDim; iset++) {
          Int_t setting = fSystSettings[iset];
          // event Vz cuts
          if (setting == kCutEventVertexZLarge && !bEventVertexZLarge) continue;
          else if (setting != kCutEventVertexZLarge && !bCutReference) continue;

          if (GetSystematicClassIndex(fTrackCutBits,setting)) {

            const Int_t signIndex = (fSign < 0); // +1 -> 0, -1 -> 1

            Double_t eff = 1e-5;

            // check pid
            std::vector<Double_t> nSigmaLimits = GetNSigmas(setting);
            Double_t nSigmaTPCDown = nSigmaLimits[0];
            Double_t nSigmaTPCUp   = nSigmaLimits[1];
            Double_t nSigmaTOFDown = nSigmaLimits[2];
            Double_t nSigmaTOFUp   = nSigmaLimits[3];

            Bool_t prTPC = (fNSigmasPrTPC > nSigmaTPCDown && fNSigmasPrTPC <= nSigmaTPCUp);
            Bool_t prTOF = (fNSigmasPrTOF > nSigmaTOFDown && fNSigmasPrTOF <= nSigmaTOFUp);

            Bool_t passesPIDCut = (ptotCut < fTOFMomCut && prTPC) || (ptotCut >= fTOFMomCut && prTOF && prTPC);

            if (passesPIDCut) {
              recProtonCounter[iset][ieta][imom][signIndex]++;
              // sum the qs
              if (fEffMatrixGenPos) {
                eff = GetTrackEfficiency(2, ptotCut, fEta, setting, fSign, 2);
              }
              for (size_t iOrderR = 0; iOrderR < orderR; iOrderR++) {
                for (size_t iOrderS = 0; iOrderS < orderS; iOrderS++) {
                  Int_t r = iOrderR + 1;
                  Int_t s = iOrderS + 1;
                  if (eff > 0) arrQ[iset][ieta][imom][iOrderR][iOrderS] += pow(fSign, r) / pow(eff, s);
                }
              }
            }

            if (fFillDebug){
              (*fTreeSRedirector)<<"debug"<<
              "eff=" << eff <<
              "nSigmaTPCDown=" << nSigmaTPCDown <<
              "nSigmaTPCUp=" << nSigmaTPCUp <<
              "nSigmaTOFDown=" << nSigmaTOFDown <<
              "prTPC=" << prTPC <<
              "prTOF=" << prTOF <<
              "passesPIDCut=" << passesPIDCut <<
              "ptotCut=" << ptotCut <<
              "p=" << fPVertex <<  // momentum at vertex
              "pT=" << fPt <<  // transverse momentum
              "fEta=" << fEta <<
              "setting=" << setting <<
              "fSign=" << fSign <<
              "momAcc=" << momAcc <<
              "etaAcc=" << etaAcc <<
              "dEdx=" << fTPCSignal <<  // dEdx of the track
              "cent=" << fCentrality << // centrality from V0
              "\n";
            }
          }
        } // ======= end of settings loop =======
      } // ======= end of momentum loop =======
    } // ======= end of eta loop =======
  } // ======= end of track loop =======

  for (size_t ieta = 0; ieta < etaDim; ieta++) {
    for (size_t imom = 0; imom < momentumDim; imom++) {
      for (size_t iset = 0; iset < settingsDim; iset++) {
        Int_t setting = fSystSettings[iset];

        TVectorF fMomNetPrRec(nMoments);

        Int_t recPos = recProtonCounter[iset][ieta][imom][0];
        Int_t recNeg = recProtonCounter[iset][ieta][imom][1];

        fMomNetPrRec[kA]    = recPos;
        fMomNetPrRec[kB]    = recNeg;
        fMomNetPrRec[kAA]   = recPos * recPos;
        fMomNetPrRec[kBB]   = recNeg * recNeg;
        fMomNetPrRec[kAB]   = recPos * recNeg;
        fMomNetPrRec[kAAA]  = recPos * recPos * recPos;
        fMomNetPrRec[kBBB]  = recNeg * recNeg * recNeg;
        fMomNetPrRec[kAAB]  = recPos * recPos * recNeg;
        fMomNetPrRec[kBBA]  = recNeg * recNeg * recPos;
        fMomNetPrRec[kABBB] = recPos * recNeg * recNeg * recNeg;
        fMomNetPrRec[kAABB] = recPos * recPos * recNeg * recNeg;
        fMomNetPrRec[kAAAB] = recPos * recPos * recPos * recNeg;
        fMomNetPrRec[kAAAA] = recPos * recPos * recPos * recPos;
        fMomNetPrRec[kBBBB] = recNeg * recNeg * recNeg * recNeg;
        fMomNetPrRec[kAmB]  = recPos - recNeg;

        TMatrixF qMatrix(orderR, orderS);

        for (size_t iOrderR = 0; iOrderR < orderR; ++iOrderR) {
          for (size_t iOrderS = 0; iOrderS < orderS; ++iOrderS) {
            qMatrix(iOrderR, iOrderS) = arrQ[iset][ieta][imom][iOrderR][iOrderS];
          }
        }

        // fill tree which contains moments
        if(!fTreeSRedirector) return;
        if (counterTracksRec > 0) {
          (*fTreeSRedirector) << "cutBased" <<
          "gid=" << fEventGID << // global event ID
          "syst=" << setting << // systematic setting index
          "isample=" << sampleNo << // sample id for subsample method
          "cent=" << fCentrality << // centrality from V0
          "vZ=" << fVz << // event vertex z
          //
          "pDown=" << fpDownArr[imom] << // lower edge of momentum bin
          "pUp=" << fpUpArr[imom] << // upper edge of momentum bin
          "etaDown=" << fetaDownArr[ieta] << // lower edge of eta bin
          "etaUp=" << fetaUpArr[ieta] << // upper edge of eta bin
          //
          "netPrMomRec.=" << &fMomNetPrRec << // moments up to 4th order with cut based method
          "qMatrix.=" << &qMatrix << // matrix with qs for efficiency correction
          "\n";
        }
      } // ======= end of settings loop =======
    } // ======= end of momentum loop =======
  } // ======= end of eta loop =======
}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::FillMCFull_NetParticles()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  // Assign subsample index
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  // Int_t runNumber = (Int_t)fESD->GetRunNumber();
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillMCFull_NetParticles ===== " << std::endl;
  //
  // check if the MC event is full
  Int_t nStackTracks = fMCEvent->GetNumberOfTracks();
  if (nStackTracks>1) fHistCentralityImpPar->Fill(fCentImpBin);
  else return;
  //
  // Count charged primary particles in gen level
  Int_t nChGen = 0;
  for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
  {
    AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
    if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
    Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
    Int_t sign = (pdg>0) ? 1 : -1;
    Double_t ptotMCgen = trackMCgen->P();
    Double_t etaMCgen = trackMCgen->Eta();
    Bool_t etaAcc = (etaMCgen>=-0.8 && etaMCgen<=0.8);
    Bool_t momAcc  = (ptotMCgen>=0.15  && ptotMCgen<100.);
    if (TMath::Abs(sign)==1) nChGen++;
  }
  //
  // Get V0 amplitudes
  Int_t multV0A = 0, multV0C = 0;
  for (Int_t i=0;i<32;i++) { multV0A += fESD->GetVZEROData()-> GetMultiplicityV0A(i); }
  for (Int_t i=0;i<32;i++) { multV0C += fESD->GetVZEROData()-> GetMultiplicityV0C(i); }
  //
  // ======================================================================
  // --------------   MC information with ideal PID   ---------------------
  // ======================================================================
  //
  const Int_t nParticles  = 3;
  const Int_t nMoments    = 14;
  const Int_t nPileUpSettings = 2;
  Int_t nOriginType = (fTrackOriginOnlyPrimary!=1) ? 4 : 1;
  //
  // counters with resonances         // counters without resonances
  TVectorF genPos(nParticles);        TVectorF nRgenPos(nParticles);
  TVectorF genNeg(nParticles);        TVectorF nRgenNeg(nParticles);
  TVectorF fMomNetPiGen(nMoments);    TVectorF fNRMomNetPiGen(nMoments);
  TVectorF fMomNetKaGen(nMoments);    TVectorF fNRMomNetKaGen(nMoments);
  TVectorF fMomNetPrGen(nMoments);    TVectorF fNRMomNetPrGen(nMoments);
  TVectorF recPos(nParticles);        TVectorF nRrecPos(nParticles);
  TVectorF recNeg(nParticles);        TVectorF nRrecNeg(nParticles);
  TVectorF fMomNetPiRec(nMoments);    TVectorF fNRMomNetPiRec(nMoments);
  TVectorF fMomNetKaRec(nMoments);    TVectorF fNRMomNetKaRec(nMoments);
  TVectorF fMomNetPrRec(nMoments);    TVectorF fNRMomNetPrRec(nMoments);

  const size_t etaDim = fetaDownArr.size();
  const size_t momDim = fNMomBinsMC;
  Bool_t isInBinRec[etaDim][momDim];
  Bool_t isInBinGen[etaDim][momDim];
  Int_t arrGenPos[etaDim][momDim][nParticles];
  Int_t arrGenNeg[etaDim][momDim][nParticles];
  Int_t arrRecPos[etaDim][momDim][nParticles];
  Int_t arrRecNeg[etaDim][momDim][nParticles];
  Int_t arrNrGenPos[etaDim][momDim][nParticles];
  Int_t arrNrGenNeg[etaDim][momDim][nParticles];
  Int_t arrNrRecPos[etaDim][momDim][nParticles];
  Int_t arrNrRecNeg[etaDim][momDim][nParticles];
  Int_t arrNtracksRec[etaDim][momDim];
  Int_t arrNtracksGen[etaDim][momDim];
  Int_t arrNtracksTPC[etaDim][momDim];
  Int_t arrNtracksITS[etaDim][momDim];
  //
  Bool_t bCutReference           = (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15);
  Bool_t bEventVertexZLarge      = (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1);
  //
  // setting scan
  for (size_t iset=0; iset<fSystSettings.size(); iset++) {
    Int_t setting = fSystSettings[iset];
    //
    // event Vz cuts
    if (setting == kCutEventVertexZLarge && !bEventVertexZLarge) continue;
    else if (setting != kCutEventVertexZLarge && !bCutReference) continue;
    //
    for (Int_t iorig=0; iorig<nOriginType; iorig++) {
      for (Int_t ipileup = 0; ipileup < nPileUpSettings; ipileup++) {
        // initialize counters
        for (size_t iEta = 0; iEta < etaDim; iEta++) {
          for (size_t iMom = 0; iMom < momDim; iMom++) {
            isInBinRec[iEta][iMom] = kFALSE;
            isInBinGen[iEta][iMom] = kFALSE;
            arrNtracksRec[iEta][iMom] = 0;
            arrNtracksGen[iEta][iMom] = 0;
            arrNtracksTPC[iEta][iMom] = 0;
            arrNtracksITS[iEta][iMom] = 0;
            for (Int_t iPart = 0; iPart < nParticles; iPart++) {
              arrGenPos[iEta][iMom][iPart] = 0;
              arrGenNeg[iEta][iMom][iPart] = 0;
              arrRecPos[iEta][iMom][iPart] = 0;
              arrRecNeg[iEta][iMom][iPart] = 0;
              arrNrGenPos[iEta][iMom][iPart] = 0;
              arrNrGenNeg[iEta][iMom][iPart] = 0;
              arrNrRecPos[iEta][iMom][iPart] = 0;
              arrNrRecNeg[iEta][iMom][iPart] = 0;
            }
          }
        }
        // -----------------------------------------------------------------------------------------
        // ----------------------------   MC generated pure MC particles  --------------------------
        // -----------------------------------------------------------------------------------------
        //
        for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++) {
          fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.;
          AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
          if (!trackMCgen) continue;
          //
          // Select real trigger event and reject other pile up vertices
          if (ipileup==0 && IsFromPileup(iTrack)) continue;
          //
          // kinematic variables
          Double_t ptotMCgen = 0.;
          if(fUsePtCut==0) ptotMCgen = trackMCgen->P();
          if(fUsePtCut==1) ptotMCgen = trackMCgen->P();
          if(fUsePtCut==2) ptotMCgen = trackMCgen->Pt();
          Float_t phiMCGen  = trackMCgen->Phi();
          Double_t etaMCgen = trackMCgen->Eta();
          Bool_t etaAccLargeWindow = (etaMCgen>=-0.8 && etaMCgen<=0.8);
          Bool_t momAccLargeWindow = (ptotMCgen>=0.15  && ptotMCgen<100.);
          if ( !(etaAccLargeWindow && etaAccLargeWindow) ) continue;
          //
          fMCGeneratorIndex = trackMCgen->GetGeneratorIndex();
          //
          TString genname = "";
          fMCEvent->GetCocktailGenerator(trackMCgen->GetLabel(), genname);
          Bool_t isHijing = genname.Contains("Hijing", TString::kIgnoreCase);
          Bool_t isPythia = genname.Contains("Pythia", TString::kIgnoreCase);
          //
          // check the origin of the track
          Bool_t bPrim     = fMCStack->IsPhysicalPrimary(iTrack);
          Bool_t bMaterial = fMCStack->IsSecondaryFromMaterial(iTrack);
          Bool_t bWeak     = fMCStack->IsSecondaryFromWeakDecay(iTrack);
          Bool_t bAcceptOrigin = kFALSE;
          if (iorig==0) bAcceptOrigin = bPrim;
          if (iorig==1) bAcceptOrigin = (bPrim || bWeak);
          if (iorig==2) bAcceptOrigin = (bPrim || bMaterial);
          if (iorig==3) bAcceptOrigin = (bPrim || bMaterial || bWeak);
          if (!bAcceptOrigin) continue;
          //
          // select particle of interest
          Int_t iPart = -10;
          Int_t pdg = trackMCgen->Particle()->GetPdgCode();
          if (TMath::Abs(pdg) == kPDGpi) {iPart = 0; fPiMCgen = iPart;} // select pi+
          if (TMath::Abs(pdg) == kPDGka) {iPart = 1; fKaMCgen = iPart;} // select ka+
          if (TMath::Abs(pdg) == kPDGpr) {iPart = 2; fPrMCgen = iPart;} // select pr+
          if (iPart == -10) continue; // perfect PID cut
          //
          for (size_t iEta = 0; iEta < etaDim; iEta++) {
            for (size_t iMom = 0; iMom < momDim; iMom++) {
              Bool_t etaAcc = (etaMCgen>=fetaDownArr[iEta] && etaMCgen<=fetaUpArr[iEta]);
              Bool_t momAcc = (ptotMCgen>=fpDownArr[iMom]  && ptotMCgen<fpUpArr[iMom]);
              isInBinGen[iEta][iMom] = (etaAcc && momAcc);
            }
          }
          //
          // fill the moments
          if (bPrim && !fIsMCPileup){
            for (size_t iEta = 0; iEta < etaDim; iEta++) {
              for (size_t iMom = 0; iMom < momDim; iMom++) {
                if (isInBinGen[iEta][iMom]) {
                  arrNtracksGen[iEta][iMom]++;

                  if ( fPiMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][kPi]++;
                  if ( fKaMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][kKa]++;
                  if ( fPrMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][kPr]++;
                  //
                  if ( fPiMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][kPi]++;
                  if ( fKaMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][kKa]++;
                  if ( fPrMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][kPr]++;
                  //
                  // reject resonances
                  Bool_t acceptRes = CheckIfFromAnyResonance(trackMCgen,fetaDownArr[iEta],fetaUpArr[iEta],fpDownArr[iMom],fpUpArr[iMom]);
                  if ( acceptRes ) {
                    if ( fPiMCgen>-1 && pdg<0) arrNrGenNeg[iEta][iMom][kPi]++;
                    if ( fKaMCgen>-1 && pdg<0) arrNrGenNeg[iEta][iMom][kKa]++;
                    if ( fPrMCgen>-1 && pdg<0) arrNrGenNeg[iEta][iMom][kPr]++;
                    //
                    if ( fPiMCgen>-1 && pdg>0) arrNrGenPos[iEta][iMom][kPi]++;
                    if ( fKaMCgen>-1 && pdg>0) arrNrGenPos[iEta][iMom][kKa]++;
                    if ( fPrMCgen>-1 && pdg>0) arrNrGenPos[iEta][iMom][kPr]++;
                  }
                }
              }
            }
          }
        } // ======= end of track loop for generated particles =======
        //
        // -----------------------------------------------------------------------------------------
        // ----------------------------   reconstructed MC particles  ------------------------------
        // -----------------------------------------------------------------------------------------
        //
        for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++) {
          //
          // track loop
          //
          // initialize the dummy particle id
          fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.;
          // Esd track
          fTrackCutBits=0;  // reset the bits for the next track
          AliESDtrack *trackReal = fESD->GetTrack(irectrack);
          if (trackReal == NULL) continue;
          Int_t lab = TMath::Abs(trackReal->GetLabel());           // avoid from negatif labels, they include some garbage
          AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab);
          // select pile up
          if (ipileup==0 && IsFromPileup(lab)) continue;
          //
          // check the origin of the track
          fMCGeneratorIndex = trackMCgen->GetGeneratorIndex();
          //
          TString genname = "";
          fMCEvent->GetCocktailGenerator(trackMCgen->GetLabel(), genname);
          Bool_t isHijing = genname.Contains("Hijing", TString::kIgnoreCase);
          Bool_t isPythia = genname.Contains("Pythia", TString::kIgnoreCase);
          //
          Bool_t bPrim     = fMCStack->IsPhysicalPrimary(lab);
          Bool_t bMaterial = fMCStack->IsSecondaryFromMaterial(lab);
          Bool_t bWeak     = fMCStack->IsSecondaryFromWeakDecay(lab);
          Bool_t bAcceptOrigin = kFALSE;
          if (iorig == 0) bAcceptOrigin = bPrim;
          if (iorig == 1) bAcceptOrigin = (bPrim || bWeak);
          if (iorig == 2) bAcceptOrigin = (bPrim || bMaterial);
          if (iorig == 3) bAcceptOrigin = (bPrim || bMaterial || bWeak);
          if (!bAcceptOrigin) continue;
          //
          // Identify particle wrt pdg code
          Int_t pdg = trackMCgen->Particle()->GetPdgCode();
          Int_t iPart = -10;
          if (TMath::Abs(pdg) == kPDGpi) { iPart = 0; fPiMC = trackReal->GetTPCsignal(); } // select pi+
          if (TMath::Abs(pdg) == kPDGka) { iPart = 1; fKaMC = trackReal->GetTPCsignal(); } // select ka+
          if (TMath::Abs(pdg) == kPDGpr) { iPart = 2; fPrMC = trackReal->GetTPCsignal(); } // select pr+
          if (iPart == -10) continue; // perfect PID cut
          //
          // apply detector cuts
          Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
          if (!trackReal-> GetInnerParam()) continue;
          if (!fESDtrackCutsLoose->AcceptTrack(trackReal)) continue;
          if (!(trackReal->GetTPCsignalN()>0)) continue;
          if (!ifDCAcutIfNoITSPixel) continue;
          //
          // kinematic variables
          Double_t ptotMCrec = 0.;
          if(fUsePtCut==0) ptotMCrec = trackReal->GetInnerParam()->GetP();
          if(fUsePtCut==1) ptotMCrec = trackReal->P();
          if(fUsePtCut==2) ptotMCrec = trackReal->Pt();
          Float_t phiMCRec  = trackReal->Phi();
          Double_t etaMCrec = trackReal->Eta();
          //
          for (size_t iEta = 0; iEta < etaDim; iEta++) {
            for (size_t iMom = 0; iMom < momDim; iMom++) {
              Bool_t etaAcc = (etaMCrec>=fetaDownArr[iEta] && etaMCrec<=fetaUpArr[iEta]);
              Bool_t momAcc  = (ptotMCrec>=fpDownArr[iMom]  && ptotMCrec<fpUpArr[iMom]);
              isInBinRec[iEta][iMom] = (etaAcc && momAcc);
            }
          }
          //
          // Get the cut bit information apply track cuts
          SetCutBitsAndSomeTrackVariables(trackReal,iPart);
          if (!GetSystematicClassIndex(fTrackCutBits,setting)) continue;
          //
          if (setting == 0 && iorig == 0 && ipileup == 0) {
            Double_t closestPar[3];
            GetExpecteds(trackReal,closestPar);
          }
          //
          // count first moments for given eta and momentum window
          for (size_t iEta = 0; iEta < etaDim; iEta++) {
            for (size_t iMom = 0; iMom < momDim; iMom++) {
              if (isInBinRec[iEta][iMom]) {
                if (trackReal->IsOn(AliESDtrack::kTPCrefit)) arrNtracksTPC[iEta][iMom]++;
                if (trackReal->IsOn(AliESDtrack::kITSrefit)) arrNtracksITS[iEta][iMom]++;
                arrNtracksRec[iEta][iMom]++;
                if ( fPiMC>-1 && pdg<0) arrRecNeg[kPi][iEta][iMom]++;
                if ( fKaMC>-1 && pdg<0) arrRecNeg[kKa][iEta][iMom]++;
                if ( fPrMC>-1 && pdg<0) arrRecNeg[kPr][iEta][iMom]++;
                //
                if ( fPiMC>-1 && pdg>0) arrRecPos[kPi][iEta][iMom]++;
                if ( fKaMC>-1 && pdg>0) arrRecPos[kKa][iEta][iMom]++;
                if ( fPrMC>-1 && pdg>0) arrRecPos[kPr][iEta][iMom]++;
                //
                // count first moments for given Centrality and momentum window without resonances
                Bool_t acceptRes = CheckIfFromAnyResonance(trackMCgen,fetaDownArr[iEta],fetaUpArr[iEta],fpDownArr[iMom],fpUpArr[iMom]);
                if ( acceptRes ) {
                  if ( fPiMC>-1 && pdg<0) arrNrRecNeg[kPi][iEta][iMom]++;
                  if ( fKaMC>-1 && pdg<0) arrNrRecNeg[kKa][iEta][iMom]++;
                  if ( fPrMC>-1 && pdg<0) arrNrRecNeg[kPr][iEta][iMom]++;
                  //
                  if ( fPiMC>-1 && pdg>0) arrNrRecPos[kPi][iEta][iMom]++;
                  if ( fKaMC>-1 && pdg>0) arrNrRecPos[kKa][iEta][iMom]++;
                  if ( fPrMC>-1 && pdg>0) arrNrRecPos[kPr][iEta][iMom]++;
                }
              }
            }
          }
        } // ======= end of track loop =======
        //
        // -----------------------------------------------------------------------------------------
        // --------------------   Calculation of moments on the event level  -----------------------
        // -----------------------------------------------------------------------------------------
        //
        for (size_t iEta = 0; iEta < etaDim; iEta++) {
          for (size_t iMom = 0; iMom < momDim; iMom++) {
            for(Int_t i = 0;i < nMoments; i++) {
              fMomNetPiGen[i]=0.;       fNRMomNetPiGen[i]=0.;
              fMomNetKaGen[i]=0.;       fNRMomNetKaGen[i]=0.;
              fMomNetPrGen[i]=0.;       fNRMomNetPrGen[i]=0.;
              fMomNetPiRec[i]=0.;       fNRMomNetPiRec[i]=0.;
              fMomNetKaRec[i]=0.;       fNRMomNetKaRec[i]=0.;
              fMomNetPrRec[i]=0.;       fNRMomNetPrRec[i]=0.;
            }
            for (Int_t iPart = 0; iPart < nParticles; iPart++) {
              genPos[iPart] = arrGenPos[iEta][iMom][iPart];
              genNeg[iPart] = arrGenNeg[iEta][iMom][iPart];
              recPos[iPart] = arrRecPos[iEta][iMom][iPart];
              recNeg[iPart] = arrRecNeg[iEta][iMom][iPart];
              nRgenPos[iPart] = arrNrGenPos[iEta][iMom][iPart];
              nRgenNeg[iPart] = arrNrGenNeg[iEta][iMom][iPart];
              nRrecPos[iPart] = arrNrRecPos[iEta][iMom][iPart];
              nRrecNeg[iPart] = arrNrRecNeg[iEta][iMom][iPart];
            }
            // ************************************************************************
            //   Moments with resonances
            // ************************************************************************
            //
            // Generated with reosnances                                                    // Reconstructed with reosnances
            // Net Pions
            fMomNetPiGen[kA]    = genPos[kPi];                                              fMomNetPiRec[kA]    = recPos[kPi];
            fMomNetPiGen[kB]    = genNeg[kPi];                                              fMomNetPiRec[kB]    = recNeg[kPi];
            fMomNetPiGen[kAA]   = genPos[kPi]*genPos[kPi];                                  fMomNetPiRec[kAA]   = recPos[kPi]*recPos[kPi];
            fMomNetPiGen[kBB]   = genNeg[kPi]*genNeg[kPi];                                  fMomNetPiRec[kBB]   = recNeg[kPi]*recNeg[kPi];
            fMomNetPiGen[kAB]   = genPos[kPi]*genNeg[kPi];                                  fMomNetPiRec[kAB]   = recPos[kPi]*recNeg[kPi];
            fMomNetPiGen[kAAA]  = genPos[kPi]*genPos[kPi]*genPos[kPi];                      fMomNetPiRec[kAAA]  = recPos[kPi]*recPos[kPi]*recPos[kPi];
            fMomNetPiGen[kBBB]  = genNeg[kPi]*genNeg[kPi]*genNeg[kPi];                      fMomNetPiRec[kBBB]  = recNeg[kPi]*recNeg[kPi]*recNeg[kPi];
            fMomNetPiGen[kAAB]  = genPos[kPi]*genPos[kPi]*genNeg[kPi];                      fMomNetPiRec[kAAB]  = recPos[kPi]*recPos[kPi]*recNeg[kPi];
            fMomNetPiGen[kBBA]  = genNeg[kPi]*genNeg[kPi]*genPos[kPi];                      fMomNetPiRec[kBBA]  = recNeg[kPi]*recNeg[kPi]*recPos[kPi];
            fMomNetPiGen[kABBB] = genPos[kPi]*genNeg[kPi]*genNeg[kPi]*genNeg[kPi];          fMomNetPiRec[kABBB] = recPos[kPi]*recNeg[kPi]*recNeg[kPi]*recNeg[kPi];
            fMomNetPiGen[kAABB] = genPos[kPi]*genPos[kPi]*genNeg[kPi]*genNeg[kPi];          fMomNetPiRec[kAABB] = recPos[kPi]*recPos[kPi]*recNeg[kPi]*recNeg[kPi];
            fMomNetPiGen[kAAAB] = genPos[kPi]*genPos[kPi]*genPos[kPi]*genNeg[kPi];          fMomNetPiRec[kAAAB] = recPos[kPi]*recPos[kPi]*recPos[kPi]*recNeg[kPi];
            fMomNetPiGen[kAAAA] = genPos[kPi]*genPos[kPi]*genPos[kPi]*genPos[kPi];          fMomNetPiRec[kAAAA] = recPos[kPi]*recPos[kPi]*recPos[kPi]*recPos[kPi];
            fMomNetPiGen[kBBBB] = genNeg[kPi]*genNeg[kPi]*genNeg[kPi]*genNeg[kPi];          fMomNetPiRec[kBBBB] = recNeg[kPi]*recNeg[kPi]*recNeg[kPi]*recNeg[kPi];

            //
            // Net Kaons
            fMomNetKaGen[kA]    = genPos[kKa];                                              fMomNetKaRec[kA]    = recPos[kKa];
            fMomNetKaGen[kB]    = genNeg[kKa];                                              fMomNetKaRec[kB]    = recNeg[kKa];
            fMomNetKaGen[kAA]   = genPos[kKa]*genPos[kKa];                                  fMomNetKaRec[kAA]   = recPos[kKa]*recPos[kKa];
            fMomNetKaGen[kBB]   = genNeg[kKa]*genNeg[kKa];                                  fMomNetKaRec[kBB]   = recNeg[kKa]*recNeg[kKa];
            fMomNetKaGen[kAB]   = genPos[kKa]*genNeg[kKa];                                  fMomNetKaRec[kAB]   = recPos[kKa]*recNeg[kKa];
            fMomNetKaGen[kAAA]  = genPos[kKa]*genPos[kKa]*genPos[kKa];                      fMomNetKaRec[kAAA]  = recPos[kKa]*recPos[kKa]*recPos[kKa];
            fMomNetKaGen[kBBB]  = genNeg[kKa]*genNeg[kKa]*genNeg[kKa];                      fMomNetKaRec[kBBB]  = recNeg[kKa]*recNeg[kKa]*recNeg[kKa];
            fMomNetKaGen[kAAB]  = genPos[kKa]*genPos[kKa]*genNeg[kKa];                      fMomNetKaRec[kAAB]  = recPos[kKa]*recPos[kKa]*recNeg[kKa];
            fMomNetKaGen[kBBA]  = genNeg[kKa]*genNeg[kKa]*genPos[kKa];                      fMomNetKaRec[kBBA]  = recNeg[kKa]*recNeg[kKa]*recPos[kKa];
            fMomNetKaGen[kABBB] = genPos[kKa]*genNeg[kKa]*genNeg[kKa]*genNeg[kKa];          fMomNetKaRec[kABBB] = recPos[kKa]*recNeg[kKa]*recNeg[kKa]*recNeg[kKa];
            fMomNetKaGen[kAABB] = genPos[kKa]*genPos[kKa]*genNeg[kKa]*genNeg[kKa];          fMomNetKaRec[kAABB] = recPos[kKa]*recPos[kKa]*recNeg[kKa]*recNeg[kKa];
            fMomNetKaGen[kAAAB] = genPos[kKa]*genPos[kKa]*genPos[kKa]*genNeg[kKa];          fMomNetKaRec[kAAAB] = recPos[kKa]*recPos[kKa]*recPos[kKa]*recNeg[kKa];
            fMomNetKaGen[kAAAA] = genPos[kKa]*genPos[kKa]*genPos[kKa]*genPos[kKa];          fMomNetKaRec[kAAAA] = recPos[kKa]*recPos[kKa]*recPos[kKa]*recPos[kKa];
            fMomNetKaGen[kBBBB] = genNeg[kKa]*genNeg[kKa]*genNeg[kKa]*genNeg[kKa];          fMomNetKaRec[kBBBB] = recNeg[kKa]*recNeg[kKa]*recNeg[kKa]*recNeg[kKa];

            //
            // Net Protons
            fMomNetPrGen[kA]    = genPos[kPr];                                              fMomNetPrRec[kA]    = recPos[kPr];
            fMomNetPrGen[kB]    = genNeg[kPr];                                              fMomNetPrRec[kB]    = recNeg[kPr];
            fMomNetPrGen[kAA]   = genPos[kPr]*genPos[kPr];                                  fMomNetPrRec[kAA]   = recPos[kPr]*recPos[kPr];
            fMomNetPrGen[kBB]   = genNeg[kPr]*genNeg[kPr];                                  fMomNetPrRec[kBB]   = recNeg[kPr]*recNeg[kPr];
            fMomNetPrGen[kAB]   = genPos[kPr]*genNeg[kPr];                                  fMomNetPrRec[kAB]   = recPos[kPr]*recNeg[kPr];
            fMomNetPrGen[kAAA]  = genPos[kPr]*genPos[kPr]*genPos[kPr];                      fMomNetPrRec[kAAA]  = recPos[kPr]*recPos[kPr]*recPos[kPr];
            fMomNetPrGen[kBBB]  = genNeg[kPr]*genNeg[kPr]*genNeg[kPr];                      fMomNetPrRec[kBBB]  = recNeg[kPr]*recNeg[kPr]*recNeg[kPr];
            fMomNetPrGen[kAAB]  = genPos[kPr]*genPos[kPr]*genNeg[kPr];                      fMomNetPrRec[kAAB]  = recPos[kPr]*recPos[kPr]*recNeg[kPr];
            fMomNetPrGen[kBBA]  = genNeg[kPr]*genNeg[kPr]*genPos[kPr];                      fMomNetPrRec[kBBA]  = recNeg[kPr]*recNeg[kPr]*recPos[kPr];
            fMomNetPrGen[kABBB] = genPos[kPr]*genNeg[kPr]*genNeg[kPr]*genNeg[kPr];          fMomNetPrRec[kABBB] = recPos[kPr]*recNeg[kPr]*recNeg[kPr]*recNeg[kPr];
            fMomNetPrGen[kAABB] = genPos[kPr]*genPos[kPr]*genNeg[kPr]*genNeg[kPr];          fMomNetPrRec[kAABB] = recPos[kPr]*recPos[kPr]*recNeg[kPr]*recNeg[kPr];
            fMomNetPrGen[kAAAB] = genPos[kPr]*genPos[kPr]*genPos[kPr]*genNeg[kPr];          fMomNetPrRec[kAAAB] = recPos[kPr]*recPos[kPr]*recPos[kPr]*recNeg[kPr];
            fMomNetPrGen[kAAAA] = genPos[kPr]*genPos[kPr]*genPos[kPr]*genPos[kPr];          fMomNetPrRec[kAAAA] = recPos[kPr]*recPos[kPr]*recPos[kPr]*recPos[kPr];
            fMomNetPrGen[kBBBB] = genNeg[kPr]*genNeg[kPr]*genNeg[kPr]*genNeg[kPr];          fMomNetPrRec[kBBBB] = recNeg[kPr]*recNeg[kPr]*recNeg[kPr]*recNeg[kPr];
            //
            // ************************************************************************
            //   Moments without resonances
            // ************************************************************************
            //
            // Generated with reosnances
            //
            // Net Pions
            fNRMomNetPiGen[kA]    = nRgenPos[kPi];                                                    fNRMomNetPiRec[kA]   = nRrecPos[kPi];
            fNRMomNetPiGen[kB]    = nRgenNeg[kPi];                                                    fNRMomNetPiRec[kB]   = nRrecNeg[kPi];
            fNRMomNetPiGen[kAA]   = nRgenPos[kPi]*nRgenPos[kPi];                                      fNRMomNetPiRec[kAA]  = nRrecPos[kPi]*nRrecPos[kPi];
            fNRMomNetPiGen[kBB]   = nRgenNeg[kPi]*nRgenNeg[kPi];                                      fNRMomNetPiRec[kBB]  = nRrecNeg[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAB]   = nRgenPos[kPi]*nRgenNeg[kPi];                                      fNRMomNetPiRec[kAB]  = nRrecPos[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAAA]  = nRgenPos[kPi]*nRgenPos[kPi]*nRgenPos[kPi];                        fNRMomNetPiRec[kAAA] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecPos[kPi];
            fNRMomNetPiGen[kBBB]  = nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi];                        fNRMomNetPiRec[kBBB] = nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAAB]  = nRgenPos[kPi]*nRgenPos[kPi]*nRgenNeg[kPi];                        fNRMomNetPiRec[kAAB] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kBBA]  = nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenPos[kPi];                        fNRMomNetPiRec[kBBA] = nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecPos[kPi];
            fNRMomNetPiGen[kABBB] = nRgenPos[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi];          fNRMomNetPiRec[kABBB] = nRrecPos[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAABB] = nRgenPos[kPi]*nRgenPos[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi];          fNRMomNetPiRec[kAABB] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAAAB] = nRgenPos[kPi]*nRgenPos[kPi]*nRgenPos[kPi]*nRgenNeg[kPi];          fNRMomNetPiRec[kAAAB] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecPos[kPi]*nRrecNeg[kPi];
            fNRMomNetPiGen[kAAAA] = nRgenPos[kPi]*nRgenPos[kPi]*nRgenPos[kPi]*nRgenPos[kPi];          fNRMomNetPiRec[kAAAA] = nRrecPos[kPi]*nRrecPos[kPi]*nRrecPos[kPi]*nRrecPos[kPi];
            fNRMomNetPiGen[kBBBB] = nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi]*nRgenNeg[kPi];          fNRMomNetPiRec[kBBBB] = nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi]*nRrecNeg[kPi];

            //
            // Net Kaons
            fNRMomNetKaGen[kA]    = nRgenPos[kKa];                                                    fNRMomNetKaRec[kA]   = nRrecPos[kKa];
            fNRMomNetKaGen[kB]    = nRgenNeg[kKa];                                                    fNRMomNetKaRec[kB]   = nRrecNeg[kKa];
            fNRMomNetKaGen[kAA]   = nRgenPos[kKa]*nRgenPos[kKa];                                      fNRMomNetKaRec[kAA]  = nRrecPos[kKa]*nRrecPos[kKa];
            fNRMomNetKaGen[kBB]   = nRgenNeg[kKa]*nRgenNeg[kKa];                                      fNRMomNetKaRec[kBB]  = nRrecNeg[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAB]   = nRgenPos[kKa]*nRgenNeg[kKa];                                      fNRMomNetKaRec[kAB]  = nRrecPos[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAAA]  = nRgenPos[kKa]*nRgenPos[kKa]*nRgenPos[kKa];                        fNRMomNetKaRec[kAAA] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecPos[kKa];
            fNRMomNetKaGen[kBBB]  = nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa];                        fNRMomNetKaRec[kBBB] = nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAAB]  = nRgenPos[kKa]*nRgenPos[kKa]*nRgenNeg[kKa];                        fNRMomNetKaRec[kAAB] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kBBA]  = nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenPos[kKa];                        fNRMomNetKaRec[kBBA] = nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecPos[kKa];
            fNRMomNetKaGen[kABBB] = nRgenPos[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa];          fNRMomNetKaRec[kABBB] = nRrecPos[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAABB] = nRgenPos[kKa]*nRgenPos[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa];          fNRMomNetKaRec[kAABB] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAAAB] = nRgenPos[kKa]*nRgenPos[kKa]*nRgenPos[kKa]*nRgenNeg[kKa];          fNRMomNetKaRec[kAAAB] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecPos[kKa]*nRrecNeg[kKa];
            fNRMomNetKaGen[kAAAA] = nRgenPos[kKa]*nRgenPos[kKa]*nRgenPos[kKa]*nRgenPos[kKa];          fNRMomNetKaRec[kAAAA] = nRrecPos[kKa]*nRrecPos[kKa]*nRrecPos[kKa]*nRrecPos[kKa];
            fNRMomNetKaGen[kBBBB] = nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa]*nRgenNeg[kKa];          fNRMomNetKaRec[kBBBB] = nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa]*nRrecNeg[kKa];

            //
            // Net Protons
            fNRMomNetPrGen[kA]    = nRgenPos[kPr];                                                    fNRMomNetPrRec[kA]   = nRrecPos[kPr];
            fNRMomNetPrGen[kB]    = nRgenNeg[kPr];                                                    fNRMomNetPrRec[kB]   = nRrecNeg[kPr];
            fNRMomNetPrGen[kAA]   = nRgenPos[kPr]*nRgenPos[kPr];                                      fNRMomNetPrRec[kAA]  = nRrecPos[kPr]*nRrecPos[kPr];
            fNRMomNetPrGen[kBB]   = nRgenNeg[kPr]*nRgenNeg[kPr];                                      fNRMomNetPrRec[kBB]  = nRrecNeg[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAB]   = nRgenPos[kPr]*nRgenNeg[kPr];                                      fNRMomNetPrRec[kAB]  = nRrecPos[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAAA]  = nRgenPos[kPr]*nRgenPos[kPr]*nRgenPos[kPr];                        fNRMomNetPrRec[kAAA] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecPos[kPr];
            fNRMomNetPrGen[kBBB]  = nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr];                        fNRMomNetPrRec[kBBB] = nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAAB]  = nRgenPos[kPr]*nRgenPos[kPr]*nRgenNeg[kPr];                        fNRMomNetPrRec[kAAB] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kBBA]  = nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenPos[kPr];                        fNRMomNetPrRec[kBBA] = nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecPos[kPr];
            fNRMomNetPrGen[kABBB] = nRgenPos[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr];          fNRMomNetPrRec[kABBB] = nRrecPos[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAABB] = nRgenPos[kPr]*nRgenPos[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr];          fNRMomNetPrRec[kAABB] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAAAB] = nRgenPos[kPr]*nRgenPos[kPr]*nRgenPos[kPr]*nRgenNeg[kPr];          fNRMomNetPrRec[kAAAB] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecPos[kPr]*nRrecNeg[kPr];
            fNRMomNetPrGen[kAAAA] = nRgenPos[kPr]*nRgenPos[kPr]*nRgenPos[kPr]*nRgenPos[kPr];          fNRMomNetPrRec[kAAAA] = nRrecPos[kPr]*nRrecPos[kPr]*nRrecPos[kPr]*nRrecPos[kPr];
            fNRMomNetPrGen[kBBBB] = nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr]*nRgenNeg[kPr];          fNRMomNetPrRec[kBBBB] = nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr]*nRrecNeg[kPr];
            //
            // fill tree which contains moments
            if(!fTreeSRedirector) return;
            Int_t tempNtracksTPC = arrNtracksTPC[iEta][iMom];
            Int_t tempNtracksITS = arrNtracksITS[iEta][iMom];
            (*fTreeSRedirector)<<"momentsMCrec"<<
            "gid=" << fEventGID <<
            "binMult=" << arrNtracksRec[iEta][iMom] <<
            "isample=" << sampleNo << // sample id for subsample method
            "vZ=" << fVz << // event vertex z
            "centimp=" << fCentImpBin << // centraltiy from impact parameter
            "impPar=" << fMCImpactParameter << // impact parameter taken from MC event header
            //
            "nhard=" << fNHardScatters << // Number of hard scatterings
            "nproj=" << fNProjectileParticipants << // Number of projectiles participants
            "ntarget=" << fNTargetParticipants << // Number of target participants
            "nn=" << fNNColl << // Number of N-N collisions
            "nnw=" << fNNwColl << // Number of N-Nwounded collisions
            "nwn=" << fNwNColl << // Number of Nwounded-N collisons
            "nwnw=" << fNwNwColl << // Number of Nwounded-Nwounded collisions
            "nch=" << nChGen << // Number of charged particles in 4pi
            //
            "syst=" << setting << // systematic setting index
            "orig=" << iorig << // origin type primary, or several combinations
            "pileup=" << ipileup << // pileup or not
            "pDown=" << fpDownArr[iMom] << // lower edge of momentum bin
            "pUp=" << fpUpArr[iMom] << // upper edge of momentum bin
            "etaDown=" << fetaDownArr[iEta] << // lower edge of eta bin
            "etaUp=" << fetaUpArr[iEta] << // upper edge of eta bin
            "cent=" << fCentrality << // centrality from V0
            //
            "netPiMomGen.=" << &fMomNetPiGen << // momnets up to 3rd order for (net)pions on generated level with resonances
            "netKaMomGen.=" << &fMomNetKaGen << // momnets up to 3rd order for (net)kaons on generated level with resonances
            "netPrMomGen.=" << &fMomNetPrGen << // momnets up to 3rd order for (net)protons on generated level with resonances
            //
            "netPiMomRec.=" << &fMomNetPiRec << // momnets up to 3rd order for (net)pions on reconstruced level with resonances
            "netKaMomRec.=" << &fMomNetKaRec << // momnets up to 3rd order for (net)kaons on reconstruced level with resonances
            "netPrMomRec.=" << &fMomNetPrRec << // momnets up to 3rd order for (net)protons on reconstruced level with resonances
            //
            "nRnetPiMomGen.=" << &fNRMomNetPiGen << // momnets up to 3rd order for (net)pions on generated level without resonances
            "nRnetKaMomGen.=" << &fNRMomNetKaGen << // momnets up to 3rd order for (net)kaons on generated level without resonances
            "nRnetPrMomGen.=" << &fNRMomNetPrGen << // momnets up to 3rd order for (net)protons on generated level without resonances
            //
            "nRnetPiMomRec.=" << &fNRMomNetPiRec << // momnets up to 3rd order for (net)pions on reconstruced level without resonances
            "nRnetKaMomRec.=" << &fNRMomNetKaRec << // momnets up to 3rd order for (net)kaons on reconstruced level without resonances
            "nRnetPrMomRec.=" << &fNRMomNetPrRec << // momnets up to 3rd order for (net)protons on reconstruced level without resonances
            "\n";

          } // momemtum loop
        } // eta loop
      } // pileup loop
    } // track origin loop
  } // settings loop
}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::FillEventInfoMC()
{

  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillEventInfoMC ===== " << std::endl;
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  const Int_t nMultType = 5;
  TVectorF fMultTPC(nMultType);
  TVectorF fMultV0M(nMultType);
  for(Int_t i=0;i<nMultType; i++){ fMultTPC[i]=0.; fMultV0M[i]=0.; }
  Int_t nStackTracks = fMCEvent->GetNumberOfTracks();

  //
  // Acceptance counter
  AliMCParticle *trackMCgen;
  for (Int_t iTrack = 0; iTrack < nStackTracks; iTrack++)
  {
    trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
    if (!trackMCgen) continue;
    if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
    Int_t sign = trackMCgen->Charge();
    Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
    Int_t absPDG = TMath::Abs(pdg);
    Int_t iPart = -10;
    if (absPDG == kPDGpi) {iPart = 0;} // select pi+
    if (absPDG == kPDGka) {iPart = 1;} // select ka+
    if (absPDG == kPDGpr) {iPart = 2;} // select pr+
    if (absPDG == kPDGde) {iPart = 3;} // select de
    if (absPDG == kPDGel) {iPart = 4;} // select el
    if (absPDG == kPDGxi) {iPart = 5;} // select ksi
    if (absPDG == kPDGla) {iPart = 6;} // select La
    if (absPDG == kPDGph) {iPart = 7;} // select phi
    Float_t etaMCgen  = trackMCgen->Eta();
    //
    // multiplicity counters for V0 and TPC acceptances
    Bool_t bV0Macc = (etaMCgen < 5.1 && etaMCgen > 2.8) || (etaMCgen < -1.7 && etaMCgen > -3.7);
    Bool_t bTPCacc = (etaMCgen < 0.8 && etaMCgen > -0.8);
    //
    // for the V0M acceptance all particles
    if (bV0Macc) {
      if (iPart==0) fMultV0M[iPart]++;
      if (iPart==1) fMultV0M[iPart]++;
      if (iPart==2) fMultV0M[iPart]++;
      if (TMath::Abs(sign)>1) fMultV0M[3]++;
      if (TMath::Abs(sign)<1) fMultV0M[4]++;
    }
    //
    // for the TPC acc only charged particles
    if (bTPCacc) {
      if (iPart==0) fMultTPC[iPart]++;
      if (iPart==1) fMultTPC[iPart]++;
      if (iPart==2) fMultTPC[iPart]++;
      if (TMath::Abs(sign)>1) fMultTPC[3]++;
      if (TMath::Abs(sign)<1) fMultTPC[4]++;
    }
  }
  //
  // Fill track kinematics
  if (fMultTPC[0]>1 &&  fMultTPC[1]>1 && fMultTPC[2]>1){ // Crucial condition --> some chunks and events in EPOS are buggy
    AliMCParticle *trackMCmother;
    AliMCParticle *trackMCGmother;
    for (Int_t iTrack = 0; iTrack < nStackTracks; iTrack++)
    {
      trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
      if (!trackMCgen) continue;
      if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
      Int_t sign = trackMCgen->Charge();
      Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
      TObjString parName(trackMCgen->Particle()->GetName());
      Int_t absPDG = TMath::Abs(pdg);
      Int_t iPart = -10;
      if (absPDG == kPDGpi) {iPart = 0;} // select pi+
      if (absPDG == kPDGka) {iPart = 1;} // select ka+
      if (absPDG == kPDGpr) {iPart = 2;} // select pr+
      if (absPDG == kPDGde) {iPart = 3;} // select de
      if (absPDG == kPDGel) {iPart = 4;} // select el
      if (absPDG == kPDGxi) {iPart = 5;} // select ksi
      if (absPDG == kPDGla) {iPart = 6;} // select La
      if (absPDG == kPDGph) {iPart = 7;} // select phi
      Float_t ptotMCgen = trackMCgen->P();
      Float_t pxMCgen = trackMCgen->Px();
      Float_t pyMCgen = trackMCgen->Py();
      Float_t pzMCgen = trackMCgen->Pz();
      Float_t pTMCgen   = trackMCgen->Pt();
      Float_t phiMCGen  = trackMCgen->Phi();
      Float_t etaMCgen  = trackMCgen->Eta();
      Float_t rapMCgen  = trackMCgen->Y();
      //
      // get the pdg info for mother and daugh ter
      TObjString momName="xxx", momGName="ggg";
      Int_t labMom = 0, labGMom = 0, pdgMom = 0, pdgGMom = 0;
      labMom = trackMCgen->Particle()->GetFirstMother();
      if ((labMom>=0) && (labMom < nStackTracks)){
        pdgMom  = fMCStack->Particle(labMom)->GetPdgCode();
        momName.SetString(fMCStack->Particle(labMom)->GetName());
        trackMCmother = (AliMCParticle *)fMCEvent->GetTrack(labMom);
        //
        labGMom = trackMCmother->GetMother();
        if ((labGMom>=0) && (labGMom < nStackTracks)){
          pdgGMom = fMCStack->Particle(labGMom)->GetPdgCode();
          momGName.SetString(fMCStack->Particle(labGMom)->GetName());
        }
      }
      //
      // dump all tracks
      Bool_t bV0Macc = (etaMCgen < 5.1 && etaMCgen > 2.8) || (etaMCgen < -1.7 && etaMCgen > -3.7);
      Bool_t bTPCacc = (etaMCgen < 0.8 && etaMCgen > -0.8);
      // if (fFillTracksMCgen && (bTPCacc || bV0Macc) ){
      if (fFillTracksMCgen && bTPCacc && TMath::Abs(sign)>1 ){
        Int_t detectorAcc = (bV0Macc) ? 0 : 1;
        Bool_t ifBar = CheckIfBaryon(trackMCgen);
        Bool_t acceptRes = CheckIfFromResonance(trackMCgen);
        Bool_t acceptAnyRes = CheckIfFromAnyResonance(trackMCgen,-0.8,0.8,0.2,10.);
        (*fTreeSRedirector)<<"tracksMCgen"<<
        "gid=" << fEventGID << // global event ID
        "part=" << iPart << // particle index --> pi, ka, pr, el, de, ksi, la, phi
        "ifBar=" << ifBar <<
        "acceptRes=" << acceptRes <<
        "acceptAnyRes=" << acceptAnyRes <<
        "sign=" << sign << // sign
        "p=" << ptotMCgen << // vertex momentum
        "px=" << pxMCgen << // vertex momentum
        "py=" << pyMCgen << // vertex momentum
        "pz=" << pzMCgen << // vertex momentum
        "eta=" << etaMCgen << // mc pseudorapidity
        "rap=" << rapMCgen << // mc rapidity
        "phi=" << phiMCGen << // mc phi
        "cent=" << fCentrality; // Centrality
        if (fFillDebug){
          (*fTreeSRedirector)<<"tracksMCgen"<<
          "pdg=" << pdg << // pdg of prim particle
          "pdgMom=" << pdgMom << // pdg of mother
          "pdgGMom=" << pdgGMom << // pdg of grandmother
          "parName.=" << &parName << // name of particle
          "momName.=" << &momName << // name of mother
          "momGName.=" << &momGName ; // name of grandmother
        }
        (*fTreeSRedirector)<<"tracksMCgen"<<"\n";
      }
    }
    //
    // event shape calculations
    if (fCollisionType == 1) { // for pp collisions
      GetFlatenicityMC();
      fSpherocity = ComputeSpherocity(-1);
    }
    MakeEventPlane(1, 1, -1); // for PbPb
    (*fTreeSRedirector)<<"eventInfoMC"<<
    "gid=" << fEventGID << // global event ID
    // "event=" << fEventCountInFile <<
    // "chunk=" << fChunkIDinJob <<
    "run=" << fRunNo << // run Number
    "ntracks_neg=" << fEP_ntracks_neg << // total number of negative tracks
    "ntracks_pos=" << fEP_ntracks_pos << // total number of positive tracks
    "momtype=" << fUsePtCut << // momentum type pt, ptot_tpc, ptot_vertex
    "isample=" << sampleNo << // sample id for subsample method
    "cent=" << fCentrality << // centrality from V0
    "centimp=" << fCentImpBin << // centraltiy from impact parameter
    // event shape observables
    "spher=" << fSpherocity <<
    "flat=" << fFlatenicity <<
    "flatS=" << fFlatenicityScaled <<
    "Qx2_neg=" << fEP_2_Qx_neg <<
    "Qx2_pos=" << fEP_2_Qx_pos <<
    "Qy2_neg=" << fEP_2_Qy_neg <<
    "Qy2_pos=" << fEP_2_Qy_pos <<
    "psi2_pos=" << fEP_2_Psi_pos <<
    "psi2_neg=" << fEP_2_Psi_neg <<
    "psi2=" << fEP_2_Psi <<
    "Qx3_neg=" << fEP_3_Qx_neg <<
    "Qx3_pos=" << fEP_3_Qx_pos <<
    "Qy3_neg=" << fEP_3_Qy_neg <<
    "Qy3_pos=" << fEP_3_Qy_pos <<
    "psi3_pos=" << fEP_3_Psi_pos <<
    "psi3_neg=" << fEP_3_Psi_neg <<
    "psi3=" << fEP_3_Psi <<
    // wounded nucleon info
    "impPar=" << fMCImpactParameter << // impact parameter taken from MC event header
    "nhard=" << fNHardScatters << // Number of hard scatterings
    "nproj=" << fNProjectileParticipants << // Number of projectiles participants
    "ntarget=" << fNTargetParticipants << // Number of target participants
    "nn=" << fNNColl << // Number of N-N collisions
    "nnw=" << fNNwColl << // Number of N-Nwounded collisions
    "nwn=" << fNwNColl << // Number of Nwounded-N collisons
    "nwnw=" << fNwNwColl << // Number of Nwounded-Nwounded collisions
    // multilicity estimators
    "multTPC.=" << &fMultTPC << // momnets up to 4th order for (net)pions on gen level
    "multV0M.=" << &fMultV0M << // momnets up to 4th order for (net)kaons on gen level
    "\n";
  }


}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::FastGen_NetParticles()
{
  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  // Assign subsample index
  fEventIDinFile = fMCEvent -> GetEventNumberInFile();
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FastGen_NetParticles ===== event = " << fEventCountInFile << std::endl;
  //
  // ======================================================================
  // --------------   MC information with ideal PID   ---------------------
  // ======================================================================
  //
  Int_t nStackTracks = fMCEvent->GetNumberOfTracks();
  if (nStackTracks>1) fHistCentralityImpPar->Fill(fCentImpBin);
  else return;
  //
  // keep onepermile of the events as debug
  Bool_t eventToDebug = kFALSE;
  if ( fRandom.Rndm() < fDownscalingFactor ) eventToDebug = kTRUE;
  //
  const Int_t nParticles = 8;
  const Int_t nMoments   = 14;
  TVectorF genPos(nParticles);
  TVectorF genNeg(nParticles);
  TVectorF fMomNetPiGen(nMoments);
  TVectorF fMomNetKaGen(nMoments);
  TVectorF fMomNetPrGen(nMoments);
  TVectorF fMomNetDeGen(nMoments);
  TVectorF fMomNetElGen(nMoments);
  TVectorF fMomNetXiGen(nMoments);
  TVectorF fMomNetLaGen(nMoments);
  TVectorF fMomNetPhGen(nMoments);
  //
  const size_t etaDim = fetaDownArr.size();
  const size_t momDim = fNMomBinsMC;
  Bool_t isInBinGen[etaDim][momDim];
  Int_t arrGenPos[etaDim][momDim][nParticles];
  Int_t arrGenNeg[etaDim][momDim][nParticles];

  for (size_t iEta = 0; iEta < etaDim; iEta++) {
    for (size_t iMom = 0; iMom < momDim; iMom++) {
      isInBinGen[iEta][iMom] = kFALSE;
      for (Int_t iPart = 0; iPart < nParticles; iPart++) {
        arrGenPos[iEta][iMom][iPart] = 0;
        arrGenNeg[iEta][iMom][iPart] = 0;
      }
    }
  }
  //
  // ----------------------------   MC generated pure MC particles  --------------------------
  AliMCParticle *trackMCgen;
  AliMCParticle *trackMCmother;
  AliMCParticle *trackMCGmother;
  for (Int_t iTrack = 0; iTrack < nStackTracks; iTrack++)
  {
    //
    // initialize the dummy particle id
    fXiMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.; fLaMCgen = -100.; fPhMCgen = -100.; fDeMCgen = -100.; fElMCgen = -100.;
    trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
    if (!trackMCgen) continue;
    if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
    //
    // select particle of interest
    Int_t sign = trackMCgen->Charge();
    Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
    Float_t ptotMCgen = trackMCgen->P();
    Float_t pTMCgen   = trackMCgen->Pt();
    Float_t phiMCGen  = trackMCgen->Phi();
    Float_t etaMCgen  = trackMCgen->Eta();
    Float_t rapMCgen  = trackMCgen->Y();
    Float_t energyMCgen  = trackMCgen->E();
    UInt_t pcode = trackMCgen->MCStatusCode();
    Int_t absPDG = TMath::Abs(pdg);
    Int_t iPart = -10;
    if (absPDG == kPDGpi) {iPart = 0; fPiMCgen = iPart;} // select pi+
    if (absPDG == kPDGka) {iPart = 1; fKaMCgen = iPart;} // select ka+
    if (absPDG == kPDGpr) {iPart = 2; fPrMCgen = iPart;} // select pr+
    if (absPDG == kPDGde) {iPart = 3; fDeMCgen = iPart;} // select de+
    if (absPDG == kPDGel) {iPart = 4; fElMCgen = iPart;} // select el+
    if (absPDG == kPDGxi) {iPart = 5; fXiMCgen = iPart;} // select ksi
    if (absPDG == kPDGla) {iPart = 6; fLaMCgen = iPart;} // select La
    if (absPDG == kPDGph) {iPart = 7; fPhMCgen = iPart;} // select phi
    if (iPart == -10) continue;
    Bool_t parInterest = (fPiMCgen>-1 || fKaMCgen>-1 || fPrMCgen>-1 || fXiMCgen>-1 || fLaMCgen>-1 || fPhMCgen>-1 || fDeMCgen>-1 || fElMCgen>-1) ? kTRUE : kFALSE;
    //
    // Fil fastgen histograms for protons
    if (fPrMCgen>0)
    {
      if (pdg>0) fHistRapDistFullAccPr->Fill(etaMCgen, fCentrality, 1.);
      if (pdg<0) fHistRapDistFullAccAPr->Fill(etaMCgen, fCentrality, 1.);
    }
    //
    // dump resonance and gen distributions
    if(fFillResonances)
    {
      //
      // gen kinematics
      Int_t labMom = 0, labGMom = 0, pdgMom = 0, pdgGMom = 0;
      TObjString momName="xxx";
      TObjString momGName="ggg";
      TObjString parName(trackMCgen->Particle()->GetName());
      //
      // get the pdg info for mother and daughter
      labMom = trackMCgen->Particle()->GetFirstMother();
      if ((labMom>=0) && (labMom < nStackTracks)){
        pdgMom  = fMCStack->Particle(labMom)->GetPdgCode();
        momName.SetString(fMCStack->Particle(labMom)->GetName());
        trackMCmother = (AliMCParticle *)fMCEvent->GetTrack(labMom);
        //
        labGMom = trackMCmother->GetMother();
        if ((labGMom>=0) && (labGMom < nStackTracks)){
          pdgGMom = fMCStack->Particle(labGMom)->GetPdgCode();
          momGName.SetString(fMCStack->Particle(labGMom)->GetName());
        }
      }
      //
      Bool_t etaAccMaxWindow = (etaMCgen>=-0.8  && etaMCgen<=0.8);
      Bool_t momAccMaxWindow = (ptotMCgen>=0.15 && ptotMCgen<=100.);
      Bool_t ifBar = CheckIfBaryon(trackMCgen);
      Bool_t acceptRes = CheckIfFromResonance(trackMCgen);
      Bool_t acceptAnyRes = CheckIfFromAnyResonance(trackMCgen,-0.8,0.8,0.2,10.);
      //
      // Fill resonance tree for full acceptance
      if(!fTreeSRedirector) return;
      if ( eventToDebug && (etaAccMaxWindow && momAccMaxWindow) )  {
        (*fTreeSRedirector)<<"resonances"<<
        "gid=" << fEventGID <<
        "event=" << fEventCountInFile <<
        "chunk=" << fChunkIDinJob <<
        "ifBar=" << ifBar <<
        "acceptRes=" << acceptRes <<
        "acceptAnyRes=" << acceptAnyRes <<
        "pcode=" << pcode <<
        "part=" << iPart <<
        "parInterest=" << parInterest << // only pi, ka, and proton
        "sign=" << sign << // sign
        "E=" << energyMCgen << // energy
        "p=" << ptotMCgen << // vertex momentum
        "pT=" << pTMCgen << // transverse momentum
        "eta=" << etaMCgen << // mc eta
        "rap=" << rapMCgen << // mc eta
        "phi=" << phiMCGen << // mc eta
        "cent=" << fCentrality << // cent bin
        "imp=" << fMCImpactParameter << // cent bin
        "pdg=" << pdg << // pdg of prim particle
        "pdgMom=" << pdgMom << // pdg of mother
        "pdgGMom=" << pdgGMom << // pdg of mother
        "parName.=" << &parName << // full path - file name with ESD
        "momName.=" << &momName << // full path - file name with ESD
        "momGName.=" << &momGName << // full path - file name with ESD
        "\n";
      }
    }
    //
    // Acceptance selection
    if(fUsePtCut==0) ptotMCgen = trackMCgen->P();
    if(fUsePtCut==1) ptotMCgen = trackMCgen->P();
    if(fUsePtCut==2) ptotMCgen = trackMCgen->Pt();
    for (size_t iEta = 0; iEta < etaDim; iEta++) {
      for (size_t iMom = 0; iMom < momDim; iMom++) {
        Bool_t etaAcc = (etaMCgen>=fetaDownArr[iEta] && etaMCgen<=fetaUpArr[iEta]);
        Bool_t momAcc  = (ptotMCgen>=fpDownArr[iMom]  && ptotMCgen<fpUpArr[iMom]);
        isInBinGen[iEta][iMom] = (etaAcc && momAcc);
      }
    }
    for (size_t iEta = 0; iEta < etaDim; iEta++) {
      for (size_t iMom = 0; iMom < momDim; iMom++) {
        if (isInBinGen[iEta][iMom]) {
          if ( fPiMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][0]++;
          if ( fKaMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][1]++;
          if ( fPrMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][2]++;
          if ( fDeMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][3]++;
          if ( fElMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][4]++;
          if ( fXiMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][5]++;
          if ( fLaMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][6]++;
          if ( fPhMCgen>-1 && pdg<0) arrGenNeg[iEta][iMom][7]++;
          //
          if ( fPiMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][0]++;
          if ( fKaMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][1]++;
          if ( fPrMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][2]++;
          if ( fDeMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][3]++;
          if ( fElMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][4]++;
          if ( fXiMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][5]++;
          if ( fLaMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][6]++;
          if ( fPhMCgen>-1 && pdg>0) arrGenPos[iEta][iMom][7]++;
        }
      }
    }
  }
  for (Int_t iEta=0; iEta<fNEtaWinBinsMC; iEta++){
    for (Int_t iMom=0; iMom<fNMomBinsMC; iMom++){
      for(Int_t i = 0;i < nMoments; i++) {
        fMomNetPiGen[i]=0.;
        fMomNetKaGen[i]=0.;
        fMomNetPrGen[i]=0.;
        fMomNetDeGen[i]=0.;
        fMomNetElGen[i]=0.;
        fMomNetXiGen[i]=0.;
        fMomNetLaGen[i]=0.;
        fMomNetPhGen[i]=0.;
      }
      for (Int_t iPart = 0; iPart < nParticles; iPart++) {
        genPos[iPart] = arrGenPos[iEta][iMom][iPart];
        genNeg[iPart] = arrGenNeg[iEta][iMom][iPart];
      }
      //
      // -----------------------------------------------------------------------------------------
      // --------------------   Calculation of moments on the event level  -----------------------
      // -----------------------------------------------------------------------------------------
      //
      // Net Pions
      fMomNetPiGen[kA]    = genPos[kPi];
      fMomNetPiGen[kB]    = genNeg[kPi];
      fMomNetPiGen[kAA]   = genPos[kPi]*genPos[kPi];
      fMomNetPiGen[kBB]   = genNeg[kPi]*genNeg[kPi];
      fMomNetPiGen[kAB]   = genPos[kPi]*genNeg[kPi];
      fMomNetPiGen[kAAA]  = genPos[kPi]*genPos[kPi]*genPos[kPi];
      fMomNetPiGen[kBBB]  = genNeg[kPi]*genNeg[kPi]*genNeg[kPi];
      fMomNetPiGen[kAAB]  = genPos[kPi]*genPos[kPi]*genNeg[kPi];
      fMomNetPiGen[kBBA]  = genNeg[kPi]*genNeg[kPi]*genPos[kPi];
      fMomNetPiGen[kABBB] = genPos[kPi]*genNeg[kPi]*genNeg[kPi]*genNeg[kPi];
      fMomNetPiGen[kAABB] = genPos[kPi]*genPos[kPi]*genNeg[kPi]*genNeg[kPi];
      fMomNetPiGen[kAAAB] = genPos[kPi]*genPos[kPi]*genPos[kPi]*genNeg[kPi];
      fMomNetPiGen[kAAAA] = genPos[kPi]*genPos[kPi]*genPos[kPi]*genPos[kPi];
      fMomNetPiGen[kBBBB] = genNeg[kPi]*genNeg[kPi]*genNeg[kPi]*genNeg[kPi];
      //
      // Net Kaons
      fMomNetKaGen[kA]    = genPos[kKa];
      fMomNetKaGen[kB]    = genNeg[kKa];
      fMomNetKaGen[kAA]   = genPos[kKa]*genPos[kKa];
      fMomNetKaGen[kBB]   = genNeg[kKa]*genNeg[kKa];
      fMomNetKaGen[kAB]   = genPos[kKa]*genNeg[kKa];
      fMomNetKaGen[kAAA]  = genPos[kKa]*genPos[kKa]*genPos[kKa];
      fMomNetKaGen[kBBB]  = genNeg[kKa]*genNeg[kKa]*genNeg[kKa];
      fMomNetKaGen[kAAB]  = genPos[kKa]*genPos[kKa]*genNeg[kKa];
      fMomNetKaGen[kBBA]  = genNeg[kKa]*genNeg[kKa]*genPos[kKa];
      fMomNetKaGen[kABBB] = genPos[kKa]*genNeg[kKa]*genNeg[kKa]*genNeg[kKa];
      fMomNetKaGen[kAABB] = genPos[kKa]*genPos[kKa]*genNeg[kKa]*genNeg[kKa];
      fMomNetKaGen[kAAAB] = genPos[kKa]*genPos[kKa]*genPos[kKa]*genNeg[kKa];
      fMomNetKaGen[kAAAA] = genPos[kKa]*genPos[kKa]*genPos[kKa]*genPos[kKa];
      fMomNetKaGen[kBBBB] = genNeg[kKa]*genNeg[kKa]*genNeg[kKa]*genNeg[kKa];
      //
      // Net Protons
      fMomNetPrGen[kA]    = genPos[kPr];
      fMomNetPrGen[kB]    = genNeg[kPr];
      fMomNetPrGen[kAA]   = genPos[kPr]*genPos[kPr];
      fMomNetPrGen[kBB]   = genNeg[kPr]*genNeg[kPr];
      fMomNetPrGen[kAB]   = genPos[kPr]*genNeg[kPr];
      fMomNetPrGen[kAAA]  = genPos[kPr]*genPos[kPr]*genPos[kPr];
      fMomNetPrGen[kBBB]  = genNeg[kPr]*genNeg[kPr]*genNeg[kPr];
      fMomNetPrGen[kAAB]  = genPos[kPr]*genPos[kPr]*genNeg[kPr];
      fMomNetPrGen[kBBA]  = genNeg[kPr]*genNeg[kPr]*genPos[kPr];
      fMomNetPrGen[kABBB] = genPos[kPr]*genNeg[kPr]*genNeg[kPr]*genNeg[kPr];
      fMomNetPrGen[kAABB] = genPos[kPr]*genPos[kPr]*genNeg[kPr]*genNeg[kPr];
      fMomNetPrGen[kAAAB] = genPos[kPr]*genPos[kPr]*genPos[kPr]*genNeg[kPr];
      fMomNetPrGen[kAAAA] = genPos[kPr]*genPos[kPr]*genPos[kPr]*genPos[kPr];
      fMomNetPrGen[kBBBB] = genNeg[kPr]*genNeg[kPr]*genNeg[kPr]*genNeg[kPr];
      //
      // Net De
      fMomNetXiGen[kA]    = genPos[3];
      fMomNetXiGen[kB]    = genNeg[3];
      fMomNetXiGen[kAA]   = genPos[3]*genPos[3];
      fMomNetXiGen[kBB]   = genNeg[3]*genNeg[3];
      fMomNetXiGen[kAB]   = genPos[3]*genNeg[3];
      fMomNetXiGen[kAAA]  = genPos[3]*genPos[3]*genPos[3];
      fMomNetXiGen[kBBB]  = genNeg[3]*genNeg[3]*genNeg[3];
      fMomNetXiGen[kAAB]  = genPos[3]*genPos[3]*genNeg[3];
      fMomNetXiGen[kBBA]  = genNeg[3]*genNeg[3]*genPos[3];
      fMomNetXiGen[kABBB] = genPos[3]*genNeg[3]*genNeg[3]*genNeg[3];
      fMomNetXiGen[kAABB] = genPos[3]*genPos[3]*genNeg[3]*genNeg[3];
      fMomNetXiGen[kAAAB] = genPos[3]*genPos[3]*genPos[3]*genNeg[3];
      fMomNetXiGen[kAAAA] = genPos[3]*genPos[3]*genPos[3]*genPos[3];
      fMomNetXiGen[kBBBB] = genNeg[3]*genNeg[3]*genNeg[3]*genNeg[3];
      //
      // Net El
      fMomNetXiGen[kA]    = genPos[4];
      fMomNetXiGen[kB]    = genNeg[4];
      fMomNetXiGen[kAA]   = genPos[4]*genPos[4];
      fMomNetXiGen[kBB]   = genNeg[4]*genNeg[4];
      fMomNetXiGen[kAB]   = genPos[4]*genNeg[4];
      fMomNetXiGen[kAAA]  = genPos[4]*genPos[4]*genPos[4];
      fMomNetXiGen[kBBB]  = genNeg[4]*genNeg[4]*genNeg[4];
      fMomNetXiGen[kAAB]  = genPos[4]*genPos[4]*genNeg[4];
      fMomNetXiGen[kBBA]  = genNeg[4]*genNeg[4]*genPos[4];
      fMomNetXiGen[kABBB] = genPos[4]*genNeg[4]*genNeg[4]*genNeg[4];
      fMomNetXiGen[kAABB] = genPos[4]*genPos[4]*genNeg[4]*genNeg[4];
      fMomNetXiGen[kAAAB] = genPos[4]*genPos[4]*genPos[4]*genNeg[4];
      fMomNetXiGen[kAAAA] = genPos[4]*genPos[4]*genPos[4]*genPos[4];
      fMomNetXiGen[kBBBB] = genNeg[4]*genNeg[4]*genNeg[4]*genNeg[4];
      //
      // Net Xi
      fMomNetXiGen[kA]    = genPos[5];
      fMomNetXiGen[kB]    = genNeg[5];
      fMomNetXiGen[kAA]   = genPos[5]*genPos[5];
      fMomNetXiGen[kBB]   = genNeg[5]*genNeg[5];
      fMomNetXiGen[kAB]   = genPos[5]*genNeg[5];
      fMomNetXiGen[kAAA]  = genPos[5]*genPos[5]*genPos[5];
      fMomNetXiGen[kBBB]  = genNeg[5]*genNeg[5]*genNeg[5];
      fMomNetXiGen[kAAB]  = genPos[5]*genPos[5]*genNeg[5];
      fMomNetXiGen[kBBA]  = genNeg[5]*genNeg[5]*genPos[5];
      fMomNetXiGen[kABBB] = genPos[5]*genNeg[5]*genNeg[5]*genNeg[5];
      fMomNetXiGen[kAABB] = genPos[5]*genPos[5]*genNeg[5]*genNeg[5];
      fMomNetXiGen[kAAAB] = genPos[5]*genPos[5]*genPos[5]*genNeg[5];
      fMomNetXiGen[kAAAA] = genPos[5]*genPos[5]*genPos[5]*genPos[5];
      fMomNetXiGen[kBBBB] = genNeg[5]*genNeg[5]*genNeg[5]*genNeg[5];
      //
      // Net La
      fMomNetLaGen[kA]    = genPos[6];
      fMomNetLaGen[kB]    = genNeg[6];
      fMomNetLaGen[kAA]   = genPos[6]*genPos[6];
      fMomNetLaGen[kBB]   = genNeg[6]*genNeg[6];
      fMomNetLaGen[kAB]   = genPos[6]*genNeg[6];
      fMomNetLaGen[kAAA]  = genPos[6]*genPos[6]*genPos[6];
      fMomNetLaGen[kBBB]  = genNeg[6]*genNeg[6]*genNeg[6];
      fMomNetLaGen[kAAB]  = genPos[6]*genPos[6]*genNeg[6];
      fMomNetLaGen[kBBA]  = genNeg[6]*genNeg[6]*genPos[6];
      fMomNetLaGen[kABBB] = genPos[6]*genNeg[6]*genNeg[6]*genNeg[6];
      fMomNetLaGen[kAABB] = genPos[6]*genPos[6]*genNeg[6]*genNeg[6];
      fMomNetLaGen[kAAAB] = genPos[6]*genPos[6]*genPos[6]*genNeg[6];
      fMomNetLaGen[kAAAA] = genPos[6]*genPos[6]*genPos[6]*genPos[6];
      fMomNetLaGen[kBBBB] = genNeg[6]*genNeg[6]*genNeg[6]*genNeg[6];
      //
      // Net Phi
      fMomNetPhGen[kA]    = genPos[7];
      fMomNetPhGen[kB]    = genNeg[7];
      fMomNetPhGen[kAA]   = genPos[7]*genPos[7];
      fMomNetPhGen[kBB]   = genNeg[7]*genNeg[7];
      fMomNetPhGen[kAB]   = genPos[7]*genNeg[7];
      fMomNetPhGen[kAAA]  = genPos[7]*genPos[7]*genPos[7];
      fMomNetPhGen[kBBB]  = genNeg[7]*genNeg[7]*genNeg[7];
      fMomNetPhGen[kAAB]  = genPos[7]*genPos[7]*genNeg[7];
      fMomNetPhGen[kBBA]  = genNeg[7]*genNeg[7]*genPos[7];
      fMomNetPhGen[kABBB] = genPos[7]*genNeg[7]*genNeg[7]*genNeg[7];
      fMomNetPhGen[kAABB] = genPos[7]*genPos[7]*genNeg[7]*genNeg[7];
      fMomNetPhGen[kAAAB] = genPos[7]*genPos[7]*genPos[7]*genNeg[7];
      fMomNetPhGen[kAAAA] = genPos[7]*genPos[7]*genPos[7]*genPos[7];
      fMomNetPhGen[kBBBB] = genNeg[7]*genNeg[7]*genNeg[7]*genNeg[7];
      //
      // fill tree which contains moments
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"momentsMCgen"<<
      "gid=" << fEventGID <<
      "momtype=" << fUsePtCut <<
      "isample=" << sampleNo << // sample id for subsample method
      "cent=" << fCentrality << // centrality from V0
      "impPar=" << fMCImpactParameter << // impact parameter taken from MC event header
      //
      "nhard=" << fNHardScatters << // Number of hard scatterings
      "nproj=" << fNProjectileParticipants << // Number of projectiles participants
      "ntarget=" << fNTargetParticipants << // Number of target participants
      "nn=" << fNNColl << // Number of N-N collisions
      "nnw=" << fNNwColl << // Number of N-Nwounded collisions
      "nwn=" << fNwNColl << // Number of Nwounded-N collisons
      "nwnw=" << fNwNwColl << // Number of Nwounded-Nwounded collisions
      //
      "pDown=" << fpDownArr[iMom] << // lower edge of momentum bin
      "pUp=" << fpUpArr[iMom] << // upper edge of momentum bin
      "etaDown=" << fetaDownArr[iEta] << // lower edge of eta bin
      "etaUp=" << fetaUpArr[iEta] << // upper edge of eta bin
      //
      "netPiMomGen.=" << &fMomNetPiGen << // momnets up to 4th order for (net)pions on gen level
      "netKaMomGen.=" << &fMomNetKaGen << // momnets up to 4th order for (net)kaons on gen level
      "netPrMomGen.=" << &fMomNetPrGen << // momnets up to 4th order for (net)protons on gen level
      "netDeMomGen.=" << &fMomNetDeGen << // momnets up to 4th order for (net)protons on gen level
      "netElMomGen.=" << &fMomNetElGen << // momnets up to 4th order for (net)protons on gen level
      "netXiMomGen.=" << &fMomNetXiGen << // momnets up to 4th order for (net)xi on gen level
      "netLaMomGen.=" << &fMomNetLaGen << // momnets up to 4th order for (net)La on gen level
      "netPhMomGen.=" << &fMomNetPhGen << // momnets up to 4th order for (net)xi on gen level
      "\n";

    } // ======= end of momentum loop =======
  } // ======= end of eta loop =======

  if (fUseCouts) std::cout << " Info::marsland: ===== Out of the FastGen_NetParticles ===== event = " << fEventCountInFile << std::endl;

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::FillTreeMC()
{

  Int_t trackOrigin = -10;
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillTreeMC ===== " << std::endl;
  //
  // ======================================================================
  // ------   reconstructed MC particles with dEdx information-------------
  // ======================================================================
  //
  Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
  const AliMultiplicity *multObj = fESD->GetMultiplicity();
  Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();
  AliMCParticle *trackMCgen;
  AliMCParticle *trackMCmother;
  AliMCParticle *trackMCGmother;
  Int_t nStackTracks = fMCEvent->GetNumberOfTracks();
  for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
  {
    //
    // Esd track
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *trackReal = fESD->GetTrack(irectrack);
    if (trackReal==NULL) continue;
    Int_t lab = TMath::Abs(trackReal->GetLabel());
    //
    // Track cuts from dtector
    Bool_t fBit96_base = fESDtrackCuts_Bit96->AcceptTrack(trackReal);
    Bool_t fBit128     = fESDtrackCuts_Bit128->AcceptTrack(trackReal);
    Bool_t fBit768     = fESDtrackCuts_Bit768->AcceptTrack(trackReal);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
    if (!trackReal->GetInnerParam()) continue;     // TODO        // Ask if track is in the TPC
    if (!fESDtrackCutsLoose->AcceptTrack(trackReal))  continue;    // TODO
    if (!(trackReal->GetTPCsignalN()>0)) continue; // TODO
    //
    // Get generated track info
    trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab);
    Int_t pdg = trackMCgen->Particle()->GetPdgCode();
    TObjString parName(trackMCgen->Particle()->GetName());
    Float_t ptotMCgen = trackMCgen->P();
    Float_t pTMCgen   = trackMCgen->Pt();
    Float_t phiMCGen  = trackMCgen->Phi();
    Float_t etaMCgen  = trackMCgen->Eta();
    Float_t rapMCgen  = trackMCgen->Y();
    TObjString momName="xxx", momGName="ggg";
    Int_t labMom = 0, labGMom = 0, pdgMom = 0, pdgGMom = 0;
    labMom = trackMCgen->Particle()->GetFirstMother();
    if ((labMom>=0) && (labMom < nStackTracks)){
      pdgMom  = fMCStack->Particle(labMom)->GetPdgCode();
      momName.SetString(fMCStack->Particle(labMom)->GetName());
      trackMCmother = (AliMCParticle *)fMCEvent->GetTrack(labMom);
      //
      labGMom = trackMCmother->GetMother();
      if ((labGMom>=0) && (labGMom < nStackTracks)){
        pdgGMom = fMCStack->Particle(labGMom)->GetPdgCode();
        momGName.SetString(fMCStack->Particle(labGMom)->GetName());
      }
    }
    //
    // check the origin of the track
    TString genname = "";
    fMCEvent->GetCocktailGenerator(trackMCgen->GetLabel(), genname);
    Bool_t isHijing = genname.Contains("Hijing", TString::kIgnoreCase);
    Bool_t isPythia = genname.Contains("Pythia", TString::kIgnoreCase);
    trackOrigin = -10;
    if (fMCStack->IsPhysicalPrimary(lab))        trackOrigin = 0;
    if (fMCStack->IsSecondaryFromMaterial(lab))  trackOrigin = 1;
    if (fMCStack->IsSecondaryFromWeakDecay(lab)) trackOrigin = 2;
    if (!isHijing && !isPythia)                  trackOrigin = 3;
    if (trackOrigin<-1) continue; // TODO
    //
    // get pileup info
    Bool_t isTPCPileup=kFALSE, isITSPileup=kFALSE;
    if (fCollisionType==0){
      isTPCPileup = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(lab,fMCEvent);
      isITSPileup = AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent, "Hijing");
    }
    //
    // match the track with mc track
    Int_t iPart = -10;
    if (TMath::Abs(pdg) == kPDGpi) { iPart = 0; } // select pi
    if (TMath::Abs(pdg) == kPDGka) { iPart = 1; } // select ka
    if (TMath::Abs(pdg) == kPDGpr) { iPart = 2; } // select pr
    if (TMath::Abs(pdg) == kPDGde) { iPart = 3; } // select de
    if (TMath::Abs(pdg) == kPDGel) { iPart = 4; } // select el
    //
    Double_t closestPar[3];
    GetExpecteds(trackReal,closestPar);
    SetCutBitsAndSomeTrackVariables(trackReal,iPart);
    //
    if (trackReal-> GetInnerParam()){
      fPtotMC       = trackReal->GetInnerParam()->GetP();
      fTPCSignalMC  = trackReal->GetTPCsignal();
    }
    fEtaMC        = trackReal->Eta();
    fPtMC         = trackReal->Pt();
    fSignMC       = trackReal->GetSign();
    Float_t pMC   = trackReal->P();
    Float_t fPhiMC= trackReal->Phi();
    Int_t nTPCClusters = fESD->GetNumberOfTPCClusters();
    Int_t nITSClusters = 0;
    AliVMultiplicity *multiObj = fESD->GetMultiplicity();
    for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
    Bool_t dca11h     = TMath::Abs(fTrackDCAxy)<0.0105+0.0350/TMath::Power(fPtMC,1.1);    // 10h tuned loose cut
    Bool_t dca10h     = TMath::Abs(fTrackDCAxy)<0.0182+0.0350/TMath::Power(fPtMC,1.01);    // 10h tuned loose cut
    Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPtMC,1.01);  // 10h tuned loose cut
    UShort_t tpcFindableCls = trackReal->GetTPCNclsF();
    UShort_t tpcSharedCls = trackReal->GetTPCnclsS();
    Float_t dca[2], covar[3];
    trackReal->GetImpactParameters(dca, covar);
    Double_t tofSignalTunedOnData = trackReal->GetTOFsignalTunedOnData();
    Double_t length = trackReal->GetIntegratedLength();
    Double_t tofSignal = trackReal->GetTOFsignal();
    Double_t beta = -.05;
    if((length > 0) && (tofSignal > 0)) beta = length / 2.99792458e-2 / tofSignal;
    //
    // Get systematic settings
    TVectorF settings(kCutSettingsCount);
    for (Int_t i = 0; i < kCutSettingsCount;i++) {
      settings[i] = (Float_t) GetSystematicClassIndex(fTrackCutBits, i);
    }
    //
    // get efficiency per track
    Double_t effNoCut   = GetTrackEfficiency(2, fPtotMC, fEta, 0, fSign);
    Double_t effTOF     = GetTrackEfficiency(2, fPtotMC, fEta, 0, fSign, 1);
    Double_t effTPCTOF  = GetTrackEfficiency(2, fPtotMC, fEta, 0, fSign, 2);
    //
    Bool_t acceptRes    = CheckIfFromResonance(trackMCgen);
    Bool_t acceptAnyRes = CheckIfFromAnyResonance(trackMCgen,-0.8,0.8,0.2,10.);
    //
    // Fill MC closure tree
    if(!fTreeSRedirector) return;
    (*fTreeSRedirector)<<"tracksMCrec"<<
    "gid=" << fEventGID <<  //  global event ID
    "orig=" << trackOrigin <<   // origin of the track
    "part=" << iPart <<
    "acceptRes=" << acceptRes <<
    "acceptAnyRes=" << acceptAnyRes <<
    "dEdx=" << fTPCSignalMC <<    // dEdx of mc track
    "beta=" << beta <<
    "cutBit=" << fTrackCutBits <<  //  Systematic Cuts
    "sign=" << fSignMC <<         // sign
    "ptot=" << fPtotMC <<         // tpc momentum
    "p=" << pMC <<             // vertex momentum
    "pT=" << fPtMC <<           // transverse momentum
    "eta=" << fEtaMC <<          // mc eta
    "phi=" << fPhiMC <<          // mc eta
    "cent=" << fCentrality <<     // Centrality
    "tpcpileup=" << isTPCPileup <<
    "itspileup=" << isITSPileup <<
    "nsigmatofel=" << fNSigmasElTOF <<  // interaction rate
    "nsigmatofpi=" << fNSigmasPiTOF <<  // interaction rate
    "nsigmatofka=" << fNSigmasKaTOF <<  // interaction rate
    "nsigmatofpr=" << fNSigmasPrTOF <<  // interaction rate
    "nsigmatpcel=" << fNSigmasElTPC <<  // interaction rate
    "nsigmatpcpi=" << fNSigmasPiTPC <<  // interaction rate
    "nsigmatpcka=" << fNSigmasKaTPC <<  // interaction rate
    "nsigmatpcpr=" << fNSigmasPrTPC <<  // interaction rate
    "effnocut=" << effNoCut <<
    "efftof=" << effTOF <<
    "efftpctof=" << effTPCTOF ;
    if (fFillDebug){
      (*fTreeSRedirector)<<"tracksMCrec"<<
      "run=" << fRunNo <<  // run Number
      "settings.=" << &settings <<  //  Systematic settings
      "bField=" << fBField <<  // run Number
      "pileupbit=" << fPileUpBit <<
      "generatorIndex=" << fMCGeneratorIndex <<  // generator index
      "vZ=" << fVz <<
      "centimp=" << fCentImpBin <<
      "impPar=" << fMCImpactParameter <<      // impact parameter taken from MC event header
      "primmult=" << fNContributors <<  //  #prim tracks
      "ncltpc=" << fNcl <<  //  centrality
      "ncltpccorr="<< fNclCorr <<  //  centrality
      "tpcmult=" << fTPCMult <<  //  TPC multiplicity
      "itsmult=" << itsNumberOfTracklets <<
      "itsclmult=" << nITSClusters <<    // ITS multiplicity
      "tpcclmult=" << nTPCClusters <<    // ITS multiplicity
      //
      // track level info
      "defCut=" << fDefaultCuts <<  // default cut
      "bit96=" << fBit96_base <<  // run Number
      "bit128=" << fBit128 <<  // run Number
      "bit768=" << fBit768 <<  // run Number
      "pixCut=" << ifDCAcutIfNoITSPixel <<  // run Number
      "tpcFindableCls=" << tpcFindableCls << // number of findable clusters
      "tpcSharedCls=" << tpcSharedCls << // number of shared clusters
      "tpcSignalN=" << fTrackTPCSignalN            <<  //  number of cl used in dEdx
      "lengthInActiveZone=" << fTrackLengthInActiveZone <<  //  track length in active zone
      "tofSignal=" << tofSignal <<
      "tofSignalTOD=" << tofSignalTunedOnData <<
      "dcaxy=" << fTrackDCAxy <<
      "dcaz=" << fTrackDCAz <<
      "cRows=" << fTrackTPCCrossedRows  <<
      "chi2tpc=" << fTrackChi2TPC <<
      "missCl=" << fMissingCl <<
      "chi2tpccorr=" << fTrackChi2TPCcorr <<
      "dcabase=" << dcaBaseCut <<  //  TPC multiplicity
      "dca10h=" << dca10h <<  //  TPC multiplicity
      "dca11h=" << dca11h <<  //  TPC multiplicity
      "fCdd=" << covar[0] <<
      "fCdz=" << covar[1] <<
      "fCzz=" << covar[2] <<
      //
      // gen level info
      "pgen=" << ptotMCgen <<             // vertex momentum
      "pTgen=" << pTMCgen <<           // transverse momentum
      "etagen=" << etaMCgen <<          // mc eta
      "phigen=" << phiMCGen <<          // mc eta
      "pdg=" << pdg <<        // pdg of prim particle
      "pdgMom=" << pdgMom <<        // pdg of mother
      "pdgGMom=" << pdgGMom <<          // pdg of mother
      "parName.=" << &parName <<          // full path - file name with ESD
      "momName.=" << &momName <<          // full path - file name with ESD
      "momGName.=" << &momGName ;          // full path - file name with ESD
    }
    (*fTreeSRedirector)<<"tracksMCrec"<<"\n";

  } // ======= end of track loop for MC dEdx filling =======


  if (fUseCouts) std::cout << " Info::marsland: ===== Out of FillTreeMC ===== " << std::endl;
}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::MCclosureHigherMoments()
{

  //
  // Fill dEdx information for the TPC and also the clean kaon and protons
  //
  Int_t sampleNo = 0;
  Int_t nSubSample = 20;
  sampleNo = Int_t(fEventGID)%nSubSample;
  fEventIDinFile = fMCEvent -> GetEventNumberInFile();
  if (fUseCouts) {
    std::cout << " Info::marsland: ===== In the MCclosureHigherMoments ===== sampleNo = " << sampleNo;
    std::cout << " gid = " << fEventGID << " evtNuminFile " << fEventIDinFile << std::endl;
  }
  //
  // ======================================================================
  // ======================================================================
  //
  // ========= Efficiency Check Eta momentum and Centrality scan ==========
  const Int_t nHighMoments = 5;
  const Int_t nMoments = 3;
  for (Int_t ieta=0; ieta<fNEtaWinBinsMC; ieta++){
    for (Int_t imom=0; imom<fNMomBinsMC; imom++){
      for (Int_t icent=0; icent<fNCentbinsData-1; icent++){

        Double_t centBin = (fcentDownArr[icent]+fcentUpArr[icent])/2.;
        Float_t nTracksgen=0, nTracksrec=0, trCountgen=0, trCountrec=0;
        // vectors to hold moments
        TVectorF netPiGen(nHighMoments);   TVectorF piPosGen(nHighMoments);   TVectorF piNegGen(nHighMoments);
        TVectorF netKaGen(nHighMoments);   TVectorF kaPosGen(nHighMoments);   TVectorF kaNegGen(nHighMoments);
        TVectorF netPrGen(nHighMoments);   TVectorF prPosGen(nHighMoments);   TVectorF prNegGen(nHighMoments);
        TVectorF netPiRec(nHighMoments);   TVectorF piPosRec(nHighMoments);   TVectorF piNegRec(nHighMoments);
        TVectorF netKaRec(nHighMoments);   TVectorF kaPosRec(nHighMoments);   TVectorF kaNegRec(nHighMoments);
        TVectorF netPrRec(nHighMoments);   TVectorF prPosRec(nHighMoments);   TVectorF prNegRec(nHighMoments);
        for(Int_t i=0;i<nHighMoments; i++){
          netPiGen[i]=0.;  piPosGen[i]=0.;  piNegGen[i]=0.;
          netKaGen[i]=0.;  kaPosGen[i]=0.;  kaNegGen[i]=0.;
          netPrGen[i]=0.;  prPosGen[i]=0.;  prNegGen[i]=0.;
          netPiRec[i]=0.;  piPosRec[i]=0.;  piNegRec[i]=0.;
          netKaRec[i]=0.;  kaPosRec[i]=0.;  kaNegRec[i]=0.;
          netPrRec[i]=0.;  prPosRec[i]=0.;  prNegRec[i]=0.;
        }
        // vectors to hold moments
        TVectorF genMomentsPos(nMoments);     TVectorF recMomentsPos(nMoments);
        TVectorF genMomentsNeg(nMoments);     TVectorF recMomentsNeg(nMoments);
        TVectorF genMomentsCross(nMoments);   TVectorF recMomentsCross(nMoments);
        for(Int_t i=0;i<nMoments; i++){
          genMomentsPos[i]=0.;    recMomentsPos[i]=0.;
          genMomentsNeg[i]=0.;    recMomentsNeg[i]=0.;
          genMomentsCross[i]=0.;  recMomentsCross[i]=0.;
        }
        //
        // ************************************************************************
        //   Constructed track counters
        // ************************************************************************
        //
        for(Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++)
        { // track loop
          //
          // initialize the dummy particle id
          Int_t piMCrec =-100., kaMCrec =-100., prMCrec =-100.;
          //
          // Esd track
          fTrackCutBits=0;  // reset the bits for the next track
          AliESDtrack *trackReal = fESD->GetTrack(iTrack);
          if (trackReal==NULL) continue;
          // Get generated track info
          Int_t lab = TMath::Abs(trackReal->GetLabel());
          if (!fMCStack->IsPhysicalPrimary(lab))continue;
          TParticle *trackMC  = fMCStack->Particle(lab);
          Int_t pdg = trackMC->GetPdgCode();
          //
          // acceptance cuts
          Double_t ptotMCrec = trackReal->P();
          Double_t etaMCrec  = trackReal->Eta();
          if (etaMCrec<fetaDownArr[ieta] || etaMCrec>fetaUpArr[ieta]) continue;
          if (ptotMCrec<fpDownArr[imom]  || ptotMCrec>fpUpArr[imom]) continue;
          //
          // detector cuts
          Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
          if (!trackReal -> GetInnerParam()) continue;
          if (!fESDtrackCuts->AcceptTrack(trackReal)) continue;  // real track cuts
          if (!ifDCAcutIfNoITSPixel) continue;  // TODO
          //
          Int_t iPart = -10;
          if (TMath::Abs(pdg) == kPDGpi) {iPart = 1; piMCrec = iPart;} // select pi+
          if (TMath::Abs(pdg) == kPDGka) {iPart = 2; kaMCrec = iPart;} // select ka+
          if (TMath::Abs(pdg) == kPDGpr) {iPart = 3; prMCrec = iPart;} // select pr+
          if (iPart == -10) continue;
          // count first moments
          if ( (fCentrality>=fcentDownArr[icent]) && (fCentrality<fcentUpArr[icent]) )
          {
            nTracksrec++;
            // Count identified particles
            if ( piMCrec>-1 || kaMCrec>-1 || prMCrec>-1) trCountrec++;
            //
            if ( piMCrec>-1 && pdg<0) recMomentsNeg[kPi]++;
            if ( kaMCrec>-1 && pdg<0) recMomentsNeg[kKa]++;
            if ( prMCrec>-1 && pdg<0) recMomentsNeg[kPr]++;
            //
            if ( piMCrec>-1 && pdg>0) recMomentsPos[kPi]++;
            if ( kaMCrec>-1 && pdg>0) recMomentsPos[kKa]++;
            if ( prMCrec>-1 && pdg>0) recMomentsPos[kPr]++;
          }

        } // ======= end of rec track loop =======
        //
        // ************************************************************************
        //   Generated track counters
        // ************************************************************************
        //
        for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
        { // track loop
          //
          // Select real trigger event and reject other pile up vertices
          if (IsFromPileup(iTrack)) continue;
          //
          // initialize the dummy particle id
          Int_t piMCgen =-100., kaMCgen =-100., prMCgen =-100.;
          AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
          Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
          //
          // apply primary track and acceptance cuts
          Double_t ptotMCgen = trackMCgen->P();
          Double_t etaMCgen = (fRapidityType==0) ? trackMCgen->Eta() :  trackMCgen->Y();
          if (etaMCgen<fetaDownArr[ieta]  || etaMCgen>fetaUpArr[ieta]) continue;
          if (ptotMCgen<fpDownArr[imom]   || ptotMCgen>fpUpArr[imom]) continue;
          if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
          //
          // select particle of interest
          Int_t iPart = -10;
          if (TMath::Abs(pdg) == kPDGpi) {iPart = 1; piMCgen = iPart;} // select pi+
          if (TMath::Abs(pdg) == kPDGka) {iPart = 2; kaMCgen = iPart;} // select ka+
          if (TMath::Abs(pdg) == kPDGpr) {iPart = 3; prMCgen = iPart;} // select pr+
          if (iPart == -10) continue;
          // count first moments
          if ( (fCentrality>=fcentDownArr[icent]) && (fCentrality<fcentUpArr[icent]) )
          {
            nTracksgen++;
            // Count identified particles
            if ( piMCgen>-1 || kaMCgen>-1 || prMCgen>-1) trCountgen++;
            //
            if ( piMCgen>-1 && pdg<0) genMomentsNeg[kPi]++;
            if ( kaMCgen>-1 && pdg<0) genMomentsNeg[kKa]++;
            if ( prMCgen>-1 && pdg<0) genMomentsNeg[kPr]++;
            //
            if ( piMCgen>-1 && pdg>0) genMomentsPos[kPi]++;
            if ( kaMCgen>-1 && pdg>0) genMomentsPos[kKa]++;
            if ( prMCgen>-1 && pdg>0) genMomentsPos[kPr]++;
          }
        } // ======= end of track loop =======
        //
        // net particle cumulants  ---- GEN                                // net particle cumulants  ---- REC
        netPiGen[0]=(genMomentsPos[kPi]-genMomentsNeg[kPi]);               netPiRec[0]=(recMomentsPos[kPi]-recMomentsNeg[kPi]);
        netKaGen[0]=(genMomentsPos[kKa]-genMomentsNeg[kKa]);               netKaRec[0]=(recMomentsPos[kKa]-recMomentsNeg[kKa]);
        netPrGen[0]=(genMomentsPos[kPr]-genMomentsNeg[kPr]);               netPrRec[0]=(recMomentsPos[kPr]-recMomentsNeg[kPr]);
        netPiGen[1]=netPiGen[0]-fNetPiFirstMomentsGen[imom][icent][ieta];  netPiRec[1]=netPiRec[0]-fNetPiFirstMomentsRec[imom][icent][ieta];
        netKaGen[1]=netKaGen[0]-fNetKaFirstMomentsGen[imom][icent][ieta];  netKaRec[1]=netKaRec[0]-fNetKaFirstMomentsRec[imom][icent][ieta];
        netPrGen[1]=netPrGen[0]-fNetPrFirstMomentsGen[imom][icent][ieta];  netPrRec[1]=netPrRec[0]-fNetPrFirstMomentsRec[imom][icent][ieta];
        netPiGen[2]=netPiGen[1]*netPiGen[1];                               netPiRec[2]=netPiRec[1]*netPiRec[1];
        netKaGen[2]=netKaGen[1]*netKaGen[1];                               netKaRec[2]=netKaRec[1]*netKaRec[1];
        netPrGen[2]=netPrGen[1]*netPrGen[1];                               netPrRec[2]=netPrRec[1]*netPrRec[1];
        netPiGen[3]=netPiGen[1]*netPiGen[1]*netPiGen[1];                   netPiRec[3]=netPiRec[1]*netPiRec[1]*netPiRec[1];
        netKaGen[3]=netKaGen[1]*netKaGen[1]*netKaGen[1];                   netKaRec[3]=netKaRec[1]*netKaRec[1]*netKaRec[1];
        netPrGen[3]=netPrGen[1]*netPrGen[1]*netPrRec[1];                   netPrRec[3]=netPrRec[1]*netPrRec[1]*netPrRec[1];
        netPiGen[4]=netPiGen[1]*netPiGen[1]*netPiGen[1]*netPiGen[1];       netPiRec[4]=netPiRec[1]*netPiRec[1]*netPiRec[1]*netPiRec[1];
        netKaGen[4]=netKaGen[1]*netKaGen[1]*netKaGen[1]*netKaGen[1];       netKaRec[4]=netKaRec[1]*netKaRec[1]*netKaRec[1]*netKaRec[1];
        netPrGen[4]=netPrGen[1]*netPrGen[1]*netPrGen[1]*netPrGen[1];       netPrRec[4]=netPrRec[1]*netPrRec[1]*netPrRec[1]*netPrRec[1];
        //
        // particle cumulants  ---- GEN                                    // particle cumulants  ---- REC
        piPosGen[0]=genMomentsPos[kPi];                                    piPosRec[0]=recMomentsPos[kPi];
        kaPosGen[0]=genMomentsPos[kKa];                                    kaPosRec[0]=recMomentsPos[kKa];
        prPosGen[0]=genMomentsPos[kPr];                                    prPosRec[0]=recMomentsPos[kPr];
        piPosGen[1]=piPosGen[0]-fPiFirstMomentsGen[0][imom][icent][ieta];  piPosRec[1]=piPosRec[0]-fPiFirstMomentsRec[0][imom][icent][ieta];
        kaPosGen[1]=kaPosGen[0]-fKaFirstMomentsGen[0][imom][icent][ieta];  kaPosRec[1]=kaPosRec[0]-fKaFirstMomentsRec[0][imom][icent][ieta];
        prPosGen[1]=prPosGen[0]-fPrFirstMomentsGen[0][imom][icent][ieta];  prPosRec[1]=prPosRec[0]-fPrFirstMomentsRec[0][imom][icent][ieta];
        piPosGen[2]=piPosGen[1]*piPosGen[1];                               piPosRec[2]=piPosRec[1]*piPosRec[1];
        kaPosGen[2]=kaPosGen[1]*kaPosGen[1];                               kaPosRec[2]=kaPosRec[1]*kaPosRec[1];
        prPosGen[2]=prPosGen[1]*prPosGen[1];                               prPosRec[2]=prPosRec[1]*prPosRec[1];
        piPosGen[3]=piPosGen[1]*piPosGen[1]*piPosGen[1];                   piPosRec[3]=piPosRec[1]*piPosRec[1]*piPosRec[1];
        kaPosGen[3]=kaPosGen[1]*kaPosGen[1]*kaPosGen[1];                   kaPosRec[3]=kaPosRec[1]*kaPosRec[1]*kaPosRec[1];
        prPosGen[3]=prPosGen[1]*prPosGen[1]*prPosGen[1];                   prPosRec[3]=prPosRec[1]*prPosRec[1]*prPosRec[1];
        //
        // Anti particle cumulants  ---- GEN                               // Anti particle cumulants  ---- REC
        piNegGen[0]=genMomentsNeg[kPi];                                    piNegRec[0]=recMomentsNeg[kPi];
        kaNegGen[0]=genMomentsNeg[kKa];                                    kaNegRec[0]=recMomentsNeg[kKa];
        prNegGen[0]=genMomentsNeg[kPr];                                    prNegRec[0]=recMomentsNeg[kPr];
        piNegGen[1]=piNegGen[0]-fPiFirstMomentsGen[1][imom][icent][ieta];  piNegRec[1]=piNegRec[0]-fPiFirstMomentsRec[1][imom][icent][ieta];
        kaNegGen[1]=kaNegGen[0]-fKaFirstMomentsGen[1][imom][icent][ieta];  kaNegRec[1]=kaNegRec[0]-fKaFirstMomentsRec[1][imom][icent][ieta];
        prNegGen[1]=prNegGen[0]-fPrFirstMomentsGen[1][imom][icent][ieta];  prNegRec[1]=prNegRec[0]-fPrFirstMomentsRec[1][imom][icent][ieta];
        piNegGen[2]=piNegGen[1]*piNegGen[1];                               piNegRec[2]=piNegRec[1]*piNegRec[1];
        kaNegGen[2]=kaNegGen[1]*kaNegGen[1];                               kaNegRec[2]=kaNegRec[1]*kaNegRec[1];
        prNegGen[2]=prNegGen[1]*prNegGen[1];                               prNegRec[2]=prNegRec[1]*prNegRec[1];
        piNegGen[3]=piNegGen[1]*piNegGen[1]*piNegGen[1];                   piNegRec[3]=piNegRec[1]*piNegRec[1]*piNegRec[1];
        kaNegGen[3]=kaNegGen[1]*kaNegGen[1]*kaNegGen[1];                   kaNegRec[3]=kaNegRec[1]*kaNegRec[1]*kaNegRec[1];
        prNegGen[3]=prNegGen[1]*prNegGen[1]*prNegGen[1];                   prNegRec[3]=prNegRec[1]*prNegRec[1]*prNegRec[1];
        //
        // fill tree which contains moments
        if(!fTreeSRedirector) return;
        (*fTreeSRedirector)<<"debug"<<
        "gid=" << fEventGID <<
        "isample=" << sampleNo << // sample id for subsample method
        "impPar=" << fMCImpactParameter << // impact parameter taken from MC event header
        "centDown=" << fcentDownArr[icent] << // lower edge of cent bin
        "centUp=" << fcentUpArr[icent] << // upper edge of cent bin
        "centBin=" << centBin << // cent bin
        "pDown=" << fpDownArr[imom] << // lower edge of momentum bin
        "pUp=" << fpUpArr[imom] << // upper edge of momentum bin
        "etaDown=" << fetaDownArr[ieta] << // lower edge of eta bin
        "etaUp=" << fetaUpArr[ieta] << // upper edge of eta bin
        //
        "trCountGen=" << nTracksgen << // number of identified tracks within the given cent and mom range
        "trCountRec=" << nTracksrec << // number of identified tracks within the given cent and mom range
        //
        "momPosGen.=" << &genMomentsPos << // second moment of positive particles
        "momNegGen.=" << &genMomentsNeg << // second moment of negative particles
        "netPiGen.=" << &netPiGen << // second moments for particle+antiparticle
        "netKaGen.=" << &netKaGen << // second moment of positive particles
        "netPrGen.=" << &netPrGen << // second moment of negative particles
        //
        "piPosGen.=" << &piPosGen << // second moments for particle+antiparticle
        "kaPosGen.=" << &kaPosGen << // second moment of positive particles
        "prPosGen.=" << &prPosGen << // second moment of negative particles
        "piNegGen.=" << &piNegGen << // second moments for particle+antiparticle
        "kaNegGen.=" << &kaNegGen << // second moment of positive particles
        "prNegGen.=" << &prNegGen << // second moment of negative particles
        //
        "momPosRec.=" << &recMomentsPos << // second moment of positive particles
        "momNegRec.=" << &recMomentsNeg << // second moment of negative particles
        "netPiRec.=" << &netPiRec << // second moments for particle+antiparticle
        "netKaRec.=" << &netKaRec << // second moment of positive particles
        "netPrRec.=" << &netPrRec << // second moment of negative particles
        //
        "piPosRec.=" << &piPosRec << // second moments for particle+antiparticle
        "kaPosRec.=" << &kaPosRec << // second moment of positive particles
        "prPosRec.=" << &prPosRec << // second moment of negative particles
        "piNegRec.=" << &piNegRec << // second moments for particle+antiparticle
        "kaNegRec.=" << &kaNegRec << // second moment of positive particles
        "prNegRec.=" << &prNegRec << // second moment of negative particles
        "\n";

      } // ======= end of Centrality loop =======
    }// ======= end of momentum loop =======
  } // ======= end of eta loop =======
  //
  // ======================================================================
  //
  // ======================================================================

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::FillEffMatrix()
{

  //
  // Fill Efficiency Matrix histograms for rec and gen
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillEffMatrix ===== " << std::endl;
  //
  // check if the MC event is full
  Int_t nStackTracks = fMCEvent->GetNumberOfTracks();
  if (nStackTracks>1) fHistCentralityImpPar->Fill(fCentImpBin);
  else return;
  //
  const Int_t nParticles  = 3;
  Int_t nOriginType = (fTrackOriginOnlyPrimary!=1) ? 4 : 1;
  //
  Bool_t bCutReference      = (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15);
  Bool_t bEventVertexZLarge = (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1);
  //
  // setting scan
  for (size_t iset=0; iset<kCutSettingsCount; iset++){
    Int_t setting = iset;
    // event Vz cuts
    if (setting == kCutEventVertexZLarge && !bEventVertexZLarge) continue;
    else if (setting != kCutEventVertexZLarge && !bCutReference) continue;
    //
    for (Int_t iorig=0; iorig<nOriginType; iorig++){
      //
      // Initialize counters
      Int_t nTracksrec=0, nTracksgen=0, nTracksTPC = 0, nTracksITS = 0;
      //
      // -----------------------------------------------------------------------------------------
      // ----------------------------   MC generated pure MC particles  --------------------------
      // -----------------------------------------------------------------------------------------
      //
      for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
      {
        //
        // initialize the dummy particle id
        fElMCgen =-100.; fPiMCgen =-100.; fKaMCgen =-100.; fPrMCgen =-100.;
        AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
        if (!trackMCgen) continue;
        if (IsFromPileup(iTrack)) continue; // pile up selection
        //
        // check the origin of the track
        Bool_t bPrim     = fMCStack->IsPhysicalPrimary(iTrack);
        Bool_t bMaterial = fMCStack->IsSecondaryFromMaterial(iTrack);
        Bool_t bWeak     = fMCStack->IsSecondaryFromWeakDecay(iTrack);
        Bool_t bAcceptOrigin = kFALSE;
        if (iorig==0) bAcceptOrigin = bPrim;
        if (iorig==1) bAcceptOrigin = (bPrim || bWeak);
        if (iorig==2) bAcceptOrigin = (bPrim || bMaterial);
        if (iorig==3) bAcceptOrigin = (bPrim || bMaterial || bWeak);
        if (!bAcceptOrigin) continue;
        //
        // select particle of interest
        Int_t iPart = -10;
        Int_t pdg = trackMCgen->Particle()->GetPdgCode();
        if (TMath::Abs(pdg) == kPDGpi) {iPart = 0; fPiMCgen = iPart;} // select pi+
        if (TMath::Abs(pdg) == kPDGka) {iPart = 1; fKaMCgen = iPart;} // select ka+
        if (TMath::Abs(pdg) == kPDGpr) {iPart = 2; fPrMCgen = iPart;} // select pr+
        if (iPart == -10) continue; // perfect PID cut
        //
        Double_t ptotMCgen = 0.;
        if(fUsePtCut==0) ptotMCgen = trackMCgen->P();
        if(fUsePtCut==1) ptotMCgen = trackMCgen->P();
        if(fUsePtCut==2) ptotMCgen = trackMCgen->Pt();
        Float_t phiMCGen  = trackMCgen->Phi();
        Double_t etaMCgen = trackMCgen->Eta();
        Bool_t etaAccMaxWindow = (etaMCgen>=fEtaMin   && etaMCgen<=fEtaMax);
        Bool_t momAccMaxWindow = (ptotMCgen>=fMomDown && ptotMCgen<=fMomUp);
        //
        // Fill eff Matrix
        if (etaAccMaxWindow && momAccMaxWindow){
          //
          // Eff matrix phi etc.
          if (setting==kCutReference && bPrim){
            Double_t xxxGen[5]={Float_t(iPart),fCentrality,ptotMCgen,etaMCgen,phiMCGen};
            if (pdg>0) fHistPosEffMatrixGen->Fill(xxxGen);
            if (pdg<0) fHistNegEffMatrixGen->Fill(xxxGen);
          }
          //
          // systematic scan
          if (bPrim && (iorig == 0 || iorig == 3)) {
            // if ((iorig == 0 || iorig == 3)) {
            Double_t xxxGenSystScan[7]={0.,static_cast<Double_t>(iorig > 0),Float_t(setting),Float_t(iPart),fCentrality,ptotMCgen,etaMCgen};
            if (pdg>0) fHistPosEffMatrixScanGen->Fill(xxxGenSystScan);
            if (pdg<0) fHistNegEffMatrixScanGen->Fill(xxxGenSystScan);
            // generated does not know about TOF
            Double_t xxxGenSystScanTOF[7]={1.,static_cast<Double_t>(iorig > 0),Float_t(setting),Float_t(iPart),fCentrality,ptotMCgen,etaMCgen};
            if (pdg>0) fHistPosEffMatrixScanGen->Fill(xxxGenSystScanTOF);
            if (pdg<0) fHistNegEffMatrixScanGen->Fill(xxxGenSystScanTOF);
            Double_t xxxGenSystScanTPCTOF[7]={2.,static_cast<Double_t>(iorig > 0),Float_t(setting),Float_t(iPart),fCentrality,ptotMCgen,etaMCgen};
            if (pdg>0) fHistPosEffMatrixScanGen->Fill(xxxGenSystScanTPCTOF);
            if (pdg<0) fHistNegEffMatrixScanGen->Fill(xxxGenSystScanTPCTOF);
          }
        }
      } // ======= end of track loop for generated particles =======
      //
      // -----------------------------------------------------------------------------------------
      // ----------------------------   reconstructed MC particles  ------------------------------
      // -----------------------------------------------------------------------------------------
      //
      for(Int_t irectrack = 0; irectrack < fESD->GetNumberOfTracks(); irectrack++)
      {
        //
        // initialize the dummy particle id
        fElMC =-100.; fPiMC =-100.; fKaMC =-100.; fPrMC =-100.;
        fTrackCutBits=0;  // reset the bits for the next track
        AliESDtrack *trackReal = fESD->GetTrack(irectrack);
        if (trackReal==NULL) continue;
        Int_t lab = TMath::Abs(trackReal->GetLabel());  // avoid from negatif labels, they include some garbage
        if (IsFromPileup(lab)) continue;  // pile up selection
        //
        // check the origin of the track
        AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab); // TParticle *trackMC  = fMCStack->Particle(lab);
        TString genname = "";
        fMCEvent->GetCocktailGenerator(trackMCgen->GetLabel(), genname);
        Bool_t isHijing = genname.Contains("Hijing", TString::kIgnoreCase);
        Bool_t isPythia = genname.Contains("Pythia", TString::kIgnoreCase);
        //
        // check the origin of the track
        Bool_t bPrim     = fMCStack->IsPhysicalPrimary(lab);
        Bool_t bMaterial = fMCStack->IsSecondaryFromMaterial(lab);
        Bool_t bWeak     = fMCStack->IsSecondaryFromWeakDecay(lab);
        Bool_t bAcceptOrigin = kFALSE;
        if (iorig==0) bAcceptOrigin = bPrim;
        if (iorig==1) bAcceptOrigin = (bPrim || bWeak);
        if (iorig==2) bAcceptOrigin = (bPrim || bMaterial);
        if (iorig==3) bAcceptOrigin = (bPrim || bMaterial || bWeak);
        if (!bAcceptOrigin) continue;
        //
        // Identify particle wrt pdg code
        Int_t pdg = trackMCgen->Particle()->GetPdgCode();                     // Int_t pdg = trackMC->GetPdgCode();
        //
        Int_t iPart = -10;
        if (TMath::Abs(pdg) == kPDGpi) { iPart = 0; fPiMC = trackReal->GetTPCsignal(); } // select pi+
        if (TMath::Abs(pdg) == kPDGka) { iPart = 1; fKaMC = trackReal->GetTPCsignal(); } // select ka+
        if (TMath::Abs(pdg) == kPDGpr) { iPart = 2; fPrMC = trackReal->GetTPCsignal(); } // select pr+
        if (iPart == -10) continue; // perfect PID cut
        //
        // apply detector cuts
        Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(trackReal);
        if (!trackReal-> GetInnerParam()) continue;
        if (!fESDtrackCutsLoose->AcceptTrack(trackReal))  continue;    // Loose Cuts
        if (!(trackReal->GetTPCsignalN()>0)) continue;
        if (!ifDCAcutIfNoITSPixel) continue;
        //
        // acceptance cuts
        Double_t ptotMCrec = 0.;
        if(fUsePtCut==0) ptotMCrec = trackReal->GetInnerParam()->GetP();
        if(fUsePtCut==1) ptotMCrec = trackReal->P();
        if(fUsePtCut==2) ptotMCrec = trackReal->Pt();
        Float_t phiMCRec  = trackReal->Phi();
        Double_t etaMCrec = trackReal->Eta();
        Bool_t etaAccMaxWindow = (etaMCrec>=fEtaMin   && etaMCrec<=fEtaMax);
        Bool_t momAccMaxWindow = (ptotMCrec>=fMomDown && ptotMCrec<=fMomUp);
        //
        // Get the cut bit information apply track cuts
        Double_t closestPar[3];
        GetExpecteds(trackReal,closestPar);
        SetCutBitsAndSomeTrackVariables(trackReal,iPart);
        if (!GetSystematicClassIndex(fTrackCutBits,setting)) continue;
        //
        // Fill Efficiency Matrices
        if (etaAccMaxWindow && momAccMaxWindow){
          //
          // Eff matrix with phi etc. for the default setting
          if (setting==kCutReference && bPrim){
            Double_t xxxRec[5]={Float_t(iPart),fCentrality,ptotMCrec,etaMCrec,phiMCRec};
            if (pdg>0) fHistPosEffMatrixRec->Fill(xxxRec);
            if (pdg<0) fHistNegEffMatrixRec->Fill(xxxRec);
          }
          //
          if ((iorig == 0 || iorig == 3)) {
            // TPC eff matrix for all settings
            Double_t xxxRecSystScan[7]={0.,static_cast<Double_t>(iorig > 0),Float_t(setting),Float_t(iPart),fCentrality,ptotMCrec,etaMCrec};
            if (pdg>0) fHistPosEffMatrixScanRec->Fill(xxxRecSystScan);
            if (pdg<0) fHistNegEffMatrixScanRec->Fill(xxxRecSystScan);
            //
            // TOF+TPC eff matrix for all settings
            vector<Double_t> nSigmaLimits = GetNSigmas(setting);
            Double_t nSigmaTPCDown = nSigmaLimits[0];
            Double_t nSigmaTPCUp = nSigmaLimits[1];
            Double_t nSigmaTOFDown = nSigmaLimits[2];
            Double_t nSigmaTOFUp = nSigmaLimits[3];

            Bool_t piTPC = (fNSigmasPiTPC > nSigmaTPCDown && fNSigmasPiTPC <= nSigmaTPCUp);
            Bool_t kaTPC = (fNSigmasKaTPC > nSigmaTPCDown && fNSigmasKaTPC <= nSigmaTPCUp);
            Bool_t prTPC = (fNSigmasPrTPC > nSigmaTPCDown && fNSigmasPrTPC <= nSigmaTPCUp);

            Bool_t piTOF = (fNSigmasPiTOF > nSigmaTOFDown && fNSigmasPiTOF <= nSigmaTOFUp);
            Bool_t kaTOF = (fNSigmasKaTOF > nSigmaTOFDown && fNSigmasKaTOF <= nSigmaTOFUp);
            Bool_t prTOF = (fNSigmasPrTOF > nSigmaTOFDown && fNSigmasPrTOF <= nSigmaTOFUp);

            if ( (piTOF && iPart==0) ||  (kaTOF && iPart==1) || (prTOF && iPart==2) ) {
              Double_t xxxRecTOF[7]={1.,static_cast<Double_t>(iorig > 0),Float_t(setting),Float_t(iPart),fCentrality,ptotMCrec,etaMCrec};
              if (pdg>0) fHistPosEffMatrixScanRec->Fill(xxxRecTOF);
              if (pdg<0) fHistNegEffMatrixScanRec->Fill(xxxRecTOF);
            }

            if ( (piTOF && piTPC && iPart==0) ||  (ptotMCrec >= 0.35 && kaTOF && kaTPC && iPart==1) || (((ptotMCrec < fTOFMomCut && prTPC) || (ptotMCrec >= fTOFMomCut && prTOF && prTPC)) && iPart==2) ) {
              Double_t xxxRecTPCTOF[7]={2.,static_cast<Double_t>(iorig > 0),Float_t(setting),Float_t(iPart),fCentrality,ptotMCrec,etaMCrec};
              if (pdg>0) fHistPosEffMatrixScanRec->Fill(xxxRecTPCTOF);
              if (pdg<0) fHistNegEffMatrixScanRec->Fill(xxxRecTPCTOF);
            }
          }
        }
      } // ======= end of track loop =======
    } // track origin loop
  } // settings loop
}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::FillCleanSamples()
{

  // Fill Clean Pions from K0s
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FillCleanSamples ===== " << std::endl;

  if (fPIDResponse) {
    fPIDResponse->GetTPCResponse().SetBetheBlochParameters(1.28778e+00/50., 3.13539e+01, TMath::Exp(-3.16327e+01), 1.87901e+00, 6.41583e+00);
  }
  AliKFParticle::SetField(fESD->GetMagneticField());
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
  Double_t mm[3] = {0,0,0};
  const Double_t cProtonMass  =TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  const Double_t cPionMass    =TDatabasePDG::Instance()->GetParticle(211)->Mass();
  const Double_t cElectronMass=TDatabasePDG::Instance()->GetParticle(11)->Mass();
  //
  // Selection From Ionut
  const AliESDVertex *primaryVertex = fESD->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*primaryVertex);
  if(fV0OpenCuts) {
    fV0OpenCuts->SetEvent(fESD);
    fV0OpenCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  if(fV0StrongCuts) {
    fV0StrongCuts->SetEvent(fESD);
    fV0StrongCuts->SetPrimaryVertex(&primaryVertexKF);
  }
  //
  //
  TObjArray* listCrossV0 = fESDtrackCutsV0->GetAcceptedV0s(fESD);
  Int_t nGoodV0s         = listCrossV0->GetEntries();
  delete listCrossV0;
  //
  // Loop over V0s
  Int_t pdgV0=0; Int_t pdgP=0; Int_t pdgN=0;
  for(Int_t iV0MI = 0; iV0MI < nGoodV0s; iV0MI++) {
    //
    UInt_t v0purity = 0;
    AliESDv0 * fV0s = fESD->GetV0(iV0MI);
    Int_t lOnFlyStatus = 0;
    lOnFlyStatus = fV0s->GetOnFlyStatus();
    if (!lOnFlyStatus) {fTrackCutBits=0; continue;}
    //
    AliESDtrack* trackPosTest = fESD->GetTrack(fV0s->GetPindex());
    AliESDtrack* trackNegTest = fESD->GetTrack(fV0s->GetNindex());
    //
    // ----------------------------------------------------------------------------------------------------------
    //  Selections from ionuts
    // ----------------------------------------------------------------------------------------------------------
    //
    //
    // protect against floating point exception in asin(deltat/chipair) calculation in AliESDv0KineCuts::PsiPair
    if (!CheckPsiPair(fV0s)) {
      continue;
    }

    if(trackPosTest->GetSign() == trackNegTest->GetSign()) {fTrackCutBits=0; continue;}
    Bool_t v0ChargesAreCorrect = (trackPosTest->GetSign()==+1 ? kTRUE : kFALSE);
    trackPosTest = (!v0ChargesAreCorrect ? fESD->GetTrack(fV0s->GetNindex()) : trackPosTest);
    trackNegTest = (!v0ChargesAreCorrect ? fESD->GetTrack(fV0s->GetPindex()) : trackNegTest);
    //
    Bool_t goodK0s = kTRUE, goodLambda = kTRUE, goodALambda = kTRUE, goodGamma = kTRUE;
    if(fV0OpenCuts) {
      goodK0s = kFALSE, goodLambda = kFALSE, goodALambda = kFALSE, goodGamma = kFALSE;
      pdgV0=0; pdgP=0; pdgN=0;
      Bool_t processV0 = fV0OpenCuts->ProcessV0(fV0s, pdgV0, pdgP, pdgN);
      if (processV0 && TMath::Abs(pdgV0)==310 &&  TMath::Abs(pdgP)==kPDGpi && TMath::Abs(pdgN)==kPDGpi) {
        goodK0s = kTRUE;
        if(fK0sPionCuts && (!fK0sPionCuts->IsSelected(trackPosTest) || !fK0sPionCuts->IsSelected(trackNegTest))) goodK0s = kFALSE;
      }
      if (processV0 && pdgV0== kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) {
        goodLambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(trackPosTest)) goodLambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(trackNegTest)) goodLambda = kFALSE;
      }
      if (processV0 && pdgV0==-kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) {
        goodALambda = kTRUE;
        if(fLambdaProtonCuts && !fLambdaProtonCuts->IsSelected(trackNegTest)) goodALambda = kFALSE;
        if(fLambdaPionCuts && !fLambdaPionCuts->IsSelected(trackPosTest)) goodALambda = kFALSE;
      }
      if (processV0 && TMath::Abs(pdgV0)==22  &&  TMath::Abs(pdgP)==kPDGel && TMath::Abs(pdgN)==kPDGel) {
        goodGamma = kTRUE;
        if(fGammaElectronCuts && (!fGammaElectronCuts->IsSelected(trackPosTest) || !fGammaElectronCuts->IsSelected(trackNegTest))) goodGamma = kFALSE;
      }
    }
    //
    Bool_t veryGoodK0s = kFALSE, veryGoodLambda = kFALSE, veryGoodALambda = kFALSE, veryGoodGamma = kFALSE;
    if(fV0StrongCuts && (goodK0s || goodLambda || goodALambda || goodGamma)) {
      pdgV0=0; pdgP=0; pdgN=0;
      Bool_t processV0 = fV0StrongCuts->ProcessV0(fV0s, pdgV0, pdgP, pdgN);
      if (processV0 && goodK0s     && TMath::Abs(pdgV0)==310 &&  TMath::Abs(pdgP)==kPDGpi && TMath::Abs(pdgN)==kPDGpi) veryGoodK0s = kTRUE;
      if (processV0 && goodLambda  && pdgV0== kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) veryGoodLambda = kTRUE;
      if (processV0 && goodALambda && pdgV0==-kPDGla         && (TMath::Abs(pdgP)==kPDGpi || TMath::Abs(pdgP)==kPDGpr) && (TMath::Abs(pdgN)==kPDGpi || TMath::Abs(pdgN)==kPDGpr)) veryGoodALambda = kTRUE;
      if (processV0 && goodGamma   && TMath::Abs(pdgV0)==22  &&  TMath::Abs(pdgP)==kPDGel && TMath::Abs(pdgN)==kPDGel) veryGoodGamma = kTRUE;
    }
    //
    // assign V0 quality selection
    if (goodK0s)                           (v0purity |= 1 << 0);
    if (goodLambda || goodALambda)         (v0purity |= 1 << 1);
    if (goodGamma)                         (v0purity |= 1 << 2);
    if (veryGoodK0s)                       (v0purity |= 1 << 3);
    if (veryGoodLambda || veryGoodALambda) (v0purity |= 1 << 4);
    if (veryGoodGamma)                     (v0purity |= 1 << 5);
    //
    // ----------------------------------------------------------------------------------------------------------
    //  My cuts
    // ----------------------------------------------------------------------------------------------------------
    //
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackPosTest)) {fTrackCutBits=0; continue;} // To FIX
    if (!fESDtrackCutsCleanSamp->AcceptTrack(trackNegTest)) {fTrackCutBits=0; continue;} // To FIX
    if (!trackPosTest->GetInnerParam()) {fTrackCutBits=0; continue;}
    if (!trackNegTest->GetInnerParam()) {fTrackCutBits=0; continue;}

    if( trackPosTest->GetSign() >0 && trackNegTest->GetSign() <0){
      fV0s->GetNPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0s->GetPPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }

    if( trackPosTest->GetSign() <0 && trackNegTest->GetSign() >0){
      fV0s->GetPPxPyPz(mn[0],mn[1],mn[2]); //reconstructed cartesian momentum components of negative daughter
      fV0s->GetNPxPyPz(mp[0],mp[1],mp[2]); //reconstructed cartesian momentum components of positive daughter
    }

    fV0s->GetPxPyPz(mm[0],mm[1],mm[2]); //reconstructed cartesian momentum components of mother

    TVector3 vecN(mn[0],mn[1],mn[2]);
    TVector3 vecP(mp[0],mp[1],mp[2]);
    TVector3 vecM(mm[0],mm[1],mm[2]);

    if ((vecP.Mag() * vecM.Mag())<0.00001) {fTrackCutBits=0; continue;}
    if ((vecN.Mag() * vecM.Mag())<0.00001) {fTrackCutBits=0; continue;}

    if (abs((vecP * vecM)/(vecP.Mag() * vecM.Mag())) > 1) continue;
    if (abs((vecN * vecM)/(vecN.Mag() * vecM.Mag())) > 1) continue;

    Double_t thetaP  = acos((vecP * vecM)/(vecP.Mag() * vecM.Mag()));
    Double_t thetaN  = acos((vecN * vecM)/(vecN.Mag() * vecM.Mag()));
    if ( ((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN)) <0.00001) {fTrackCutBits=0; continue;}
    fCosPA = fV0s->GetV0CosineOfPointingAngle();
    fAlfa = ((vecP.Mag())*cos(thetaP)-(vecN.Mag())*cos(thetaN))/((vecP.Mag())*cos(thetaP)+(vecN.Mag())*cos(thetaN));
    fQt   = vecP.Mag()*sin(thetaP);
    if (fV0InvMassHists) fHistArmPod->Fill(fAlfa,fQt);
    // fV0s->ChangeMassHypothesis(22);   // ????
    // fV0s->ChangeMassHypothesis(310); // ????
    //
    // main armentoros podolanki cuts
    if (TMath::Abs(fAlfa)>0.9) {fTrackCutBits=0; continue;}
    if (fQt >0.22) {fTrackCutBits=0; continue;}
    if (fQt >0.02 && fQt<0.12 && TMath::Abs(fAlfa)<0.4) {fTrackCutBits=0; continue;}
    SelectCleanSamplesFromV0s(fV0s,trackPosTest,trackNegTest);
    //
    TLorentzVector posE, negE, photon, posP, negP, posPi, negPi, lambda, antiLambda, kaon, posProton, k0sProton;
    negE.SetXYZM(mn[0],mn[1],mn[2],cElectronMass);
    posE.SetXYZM(mp[0],mp[1],mp[2],cElectronMass);
    negPi.SetXYZM(mn[0],mn[1],mn[2],cPionMass);
    posPi.SetXYZM(mp[0],mp[1],mp[2],cPionMass);
    negP.SetXYZM(mn[0],mn[1],mn[2],cProtonMass);
    posP.SetXYZM(mp[0],mp[1],mp[2],cProtonMass);
    kaon=posPi+negPi;
    photon=posE+negE;
    lambda=posP+negPi;
    antiLambda=posPi+negP;
    //
    //
    Float_t ptotForBetaGamma0 = trackPosTest->GetInnerParam()->GetP();
    Float_t ptotForBetaGamma1 = trackNegTest->GetInnerParam()->GetP();
    Float_t ptotForBetaGammaThr = 0.2;
    Double_t posNTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kPion));
    Double_t negNTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kPion));
    Double_t posNTPCSigmaPr = (ptotForBetaGamma0>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kProton)) : 0.;
    Double_t negNTPCSigmaPr = (ptotForBetaGamma1>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kProton)) : 0.;
    Double_t posNTPCSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackPosTest, AliPID::kElectron));
    Double_t negNTPCSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackNegTest, AliPID::kElectron));
    //
    // --------------------------------------------------------------
    //  Invariant mass cuts
    // --------------------------------------------------------------
    //
    Double_t k0sMass = kaon.M();
    Double_t lambdaMass = lambda.M();
    Double_t antiLambdaMass = antiLambda.M();
    Double_t photonMass = photon.M();
    //
    Bool_t isK0sMass        = (k0sMass>0.475 && k0sMass<0.52); // ( kaon.M()>0.490 && kaon.M()<0.504 );
    Bool_t isLambdaMass     = (lambdaMass>1.109 && lambdaMass<1.122); //( lambda.M()>1.113 && lambda.M()<1.118 );
    Bool_t isAntiLambdaMass = (antiLambdaMass>1.109 && antiLambdaMass<1.12); // ( antiLambda.M()>1.113 && antiLambda.M()<1.118 );
    Bool_t isPhotonMass     = (photonMass<0.005); // (photon.M()<0.005);
    Double_t oneLegSigma = 3.5;
    if (fQt<0.02 && TMath::Abs(fAlfa)<0.5){
      // Cuts concerns Electrons
      //
      // Apply one leg cut for electrons
      if (!(negNTPCSigmaEl<2. || posNTPCSigmaEl<2.)) {fTrackCutBits=0; continue;}
      //
      if (fV0InvMassHists) fHistInvPhoton->Fill(photonMass);
      if (isK0sMass) {fTrackCutBits=0; continue;}
      if (isLambdaMass) {fTrackCutBits=0; continue;}
      if (isAntiLambdaMass) {fTrackCutBits=0; continue;}
      if (!isPhotonMass) {fTrackCutBits=0; continue;}
    } else {
      if (fV0InvMassHists) {
        fHistInvK0s->Fill(k0sMass);
        fHistInvLambda->Fill(lambdaMass);
        fHistInvAntiLambda->Fill(antiLambdaMass);
      }
      if ( !(isK0sMass || isLambdaMass || isAntiLambdaMass) ) {fTrackCutBits=0; continue;} // Apply inv mass cut
      if (fQt>0.11 && (!(negNTPCSigmaPi < oneLegSigma || posNTPCSigmaPi < oneLegSigma))) {fTrackCutBits=0; continue;} // Apply one leg cut for K0s
    }
    //
    // Set the variables to be filled in the tree
    Bool_t selectPosLegs = (posNTPCSigmaPi<oneLegSigma || posNTPCSigmaPr<oneLegSigma || posNTPCSigmaEl<oneLegSigma); // Apply one leg cut for lambda and gamma
    Bool_t selectNegLegs = (negNTPCSigmaPi<oneLegSigma || negNTPCSigmaPr<oneLegSigma || negNTPCSigmaEl<oneLegSigma); // Apply one leg cut for lambda and gamma
    //
    // positive leg
    UInt_t cutBit0 = 0, cutBit1=0, fTrackCutBits=0;
    cutBit0 = SetCutBitsAndSomeTrackVariables(trackPosTest,0);
    Float_t dEdx0    = trackPosTest->GetTPCsignal();
    Float_t itsdEdx0 = trackPosTest->GetITSsignal();
    Float_t ptot0  = trackPosTest->GetInnerParam()->GetP();
    Float_t p0     = trackPosTest->P();
    Float_t pT0    = trackPosTest->Pt();
    Float_t eta0   = trackPosTest->Eta();
    Int_t sign0    = trackPosTest->Charge();
    Float_t phi0   = trackPosTest->Phi();
    Float_t nSigmasPiTOF0 = fPIDResponse->NumberOfSigmasTOF(trackPosTest, AliPID::kPion,  fPIDResponse->GetTOFResponse().GetTimeZero());
    Float_t nSigmasPrTOF0 = fPIDResponse->NumberOfSigmasTOF(trackPosTest, AliPID::kProton,fPIDResponse->GetTOFResponse().GetTimeZero());
    //
    // negative leg
    cutBit1 = SetCutBitsAndSomeTrackVariables(trackNegTest,0);
    Float_t dEdx1  = trackNegTest->GetTPCsignal();
    Float_t itsdEdx1 = trackNegTest->GetITSsignal();
    Float_t ptot1  = trackNegTest->GetInnerParam()->GetP();
    Float_t p1     = trackNegTest->P();
    Float_t pT1    = trackNegTest->Pt();
    Float_t eta1   = trackNegTest->Eta();
    Int_t sign1    = trackNegTest->Charge();
    Float_t phi1   = trackNegTest->Phi();
    Float_t nSigmasPiTOF1 = fPIDResponse->NumberOfSigmasTOF(trackNegTest, AliPID::kPion,  fPIDResponse->GetTOFResponse().GetTimeZero());
    Float_t nSigmasPrTOF1 = fPIDResponse->NumberOfSigmasTOF(trackNegTest, AliPID::kProton,fPIDResponse->GetTOFResponse().GetTimeZero());
    //
    // --------------------------------------------------------------
    //  Fill Clean Samples tree
    // --------------------------------------------------------------
    //
    if (selectNegLegs || selectPosLegs)
    {
      if (fFillArmPodTree)
      {
        if(!fTreeSRedirector) return;
        (*fTreeSRedirector)<<"cleanSamp"<<
        "gid=" << fEventGID << // global event ID
        "eventtime=" << fTimeStamp <<
        "intrate=" << fIntRate <<  // interaction rate
        "piFromK0=" << fCleanPionsFromK0 << // K0s cut for pions
        "v0haspixel=" << fHasV0FirstITSlayer  << // ITS pixel cout
        "purity=" << v0purity <<
        "lambdaMass=" << lambdaMass << // lambda mass
        "antiLambdaMass=" << antiLambdaMass <<  // anti lambda mass
        "k0sMass=" << k0sMass << // k0s mass
        "photonMass=" << photonMass <<  // photon mass
        "cosPA=" << fCosPA << // cosine of pointing angle
        "qt=" << fQt << // qT
        "alfa=" << fAlfa <<  // alpha
        "cent=" << fCentrality << // centrality
        //
        "cutBit0=" << cutBit0 << // cut bits
        "itsdEdx0=" << itsdEdx0 << // TPC dEdx
        "dEdx0=" << dEdx0 << // TPC dEdx
        "sign0=" << sign0 <<
        "ptot0=" << ptot0 << // momentum
        "p0=" << p0 <<
        "pT0=" << pT0 <<
        "eta0=" << eta0 <<  // eta
        "phi0=" << phi0 <<  // eta
        "nSigmasPiTOF0=" << nSigmasPiTOF0 << // TOF nsigma cut for pions
        "nSigmasPrTOF0=" << nSigmasPrTOF0 << // TOF nsigma cut for protons
        //
        "cutBit1=" << cutBit1 <<  // cut bits
        "dEdx1=" << dEdx1 <<  // TPC dEdx
        "itsdEdx1=" << itsdEdx1 <<  // TPC dEdx
        "sign1=" << sign1 <<
        "ptot1=" << ptot1 << // momentum
        "p1=" << p1 <<
        "pT1=" << pT1 <<
        "eta1=" << eta1 << // eta
        "phi1=" << phi1 << // eta
        "nSigmasPiTOF1=" << nSigmasPiTOF1 << // TOF nsigma cut for pions
        "nSigmasPrTOF1=" << nSigmasPrTOF1 << // TOF nsigma cut for protons
        //
        "\n";

      }
    }
    cutBit0 = 0, cutBit1=0; fTrackCutBits=0;

  } // end of V0 loop

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::GetExpecteds(AliESDtrack *track, Double_t closestPar[3])
{

  //
  // bettaGamma is not well deifned below bg=0.01 --> below 200MeV protons and deuterons
  Double_t ptotForBetaGamma = track->GetInnerParam()->GetP();
  Double_t ptotForBetaGammaThr = 0.2;
  // if (ptotForBetaGamma<ptotForBetaGammaThr) return;
  //
  // --------------------------------------------------------------
  //  Calculates expected sigma and dEdx for a given track and returns colesest expected particle and its index
  // --------------------------------------------------------------
  //
  fNSigmasElTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  fNSigmasPiTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fNSigmasKaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //
  fNSigmasElTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron, fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion,     fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon,     fPIDResponse->GetTOFResponse().GetTimeZero());
  //
  if (ptotForBetaGamma>ptotForBetaGammaThr){
    fNSigmasPrTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    fNSigmasDeTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    fNSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton,   fPIDResponse->GetTOFResponse().GetTimeZero());
    fNSigmasDeTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron, fPIDResponse->GetTOFResponse().GetTimeZero());
  }

  //
  //
  Int_t nSigmaTmp = (fEventInfo) ? 10000 : 2;
  //
  // Electron Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasElTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxEl  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaEl = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxEl  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kElectron);
    fSigmaEl = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kElectron);
  }
  //
  // Pion Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasPiTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxPi  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion,     AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaPi = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kPion,     AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxPi  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kPion);
    fSigmaPi = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kPion);
  }
  //
  // Kaon Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasKaTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxKa  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kKaon,  AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaKa = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kKaon,   AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxKa  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kKaon);
    fSigmaKa = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kKaon);
  }
  //
  // Proton Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasPrTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxPr  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaPr = fPIDResponse->GetTPCResponse().GetExpectedSigma(track,  AliPID::kProton, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxPr  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kProton);
    fSigmaPr = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kProton);
  }
  //
  // Deuteron Expected mean and sigma within 2nsigmaTPC
  if (TMath::Abs(fNSigmasDeTPC)<nSigmaTmp && ptotForBetaGamma>ptotForBetaGammaThr) {
    // fDEdxDe  = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kDeuteron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    // fSigmaDe = fPIDResponse->GetTPCResponse().GetExpectedSigma(track,  AliPID::kDeuteron, AliTPCPIDResponse::kdEdxDefault,fPIDResponse->UseTPCEtaCorrection(),fPIDResponse->UseTPCMultiplicityCorrection(),fPIDResponse->UseTPCPileupCorrection());
    fDEdxDe  = fPIDResponse->GetExpectedSignal(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
    fSigmaDe = fPIDResponse->GetExpectedSigma(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
  }
  //
  // --------------------------------------------------------------
  //  Find Closest dEdx its corresponding particle and mass
  // --------------------------------------------------------------
  //
  Double_t values[] = {fDEdxEl, fDEdxPi, fDEdxKa, fDEdxPr, fDEdxDe};
  Double_t tpcdEdx = track->GetTPCsignal();
  Double_t smallestDiff = TMath::Abs(tpcdEdx - values[0]);
  Int_t closestIndex = 0;
  for (Int_t i = 0; i < 5; i++) {
    Double_t currentDiff = TMath::Abs(tpcdEdx - values[i]);
    if (currentDiff < smallestDiff) {
      smallestDiff = currentDiff;
      closestIndex = i;
    }
  }
  //
  // TF1 f1("f1","AliExternalTrackParam::BetheBlochAleph(x/0.1)",0.1,5)
  Double_t partMass = 0.;
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  if (closestIndex == 0 ) partMass = pdg->GetParticle(kPDGel)->Mass();   // GetParticle("e+")
  if (closestIndex == 1 ) partMass = pdg->GetParticle(kPDGpi)->Mass();  // GetParticle("pi+")
  if (closestIndex == 2 ) partMass = pdg->GetParticle(kPDGka)->Mass();  // GetParticle("K+")
  if (closestIndex == 3 ) partMass = pdg->GetParticle(kPDGpr)->Mass(); // GetParticle("proton")
  if (closestIndex == 4 ) partMass = 2.01410178;                     // pdg->GetParticle(1000010020)->Mass();
  //
  closestPar[0]=values[closestIndex];
  closestPar[1]=closestIndex;
  closestPar[2]=partMass;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskTIdentityPID::CheckIfFromResonance(AliMCParticle *trackMCgen)
{
  //
  // default is accept resonances
  Bool_t acceptRes = kTRUE;
  //
  TObjString momName="xxx";
  Int_t labMom = trackMCgen->Particle()->GetFirstMother();
  Int_t pdgMom = 0;
  if ( (labMom>=0) && (labMom < fMCEvent->GetNumberOfTracks()) ){
    pdgMom  = fMCStack->Particle(labMom)->GetPdgCode();
    momName.SetString(fMCStack->Particle(labMom)->GetName());
  }
  //
  //Check if the particle is in the black list of resonances
  for (Int_t ires=0;ires<fNResBins;ires++){
    if ( fResonances[ires]=="xxx" ){
      if (!(momName.GetString().Contains(fResonances[ires]))) {acceptRes=kFALSE; break;}
    } else {
      if (momName.GetString().Contains(fResonances[ires])) {acceptRes=kFALSE; break;}
    }
  }
  return acceptRes;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskTIdentityPID::CheckIfBaryon(AliMCParticle *trackMCgen)
{
  //
  // default is accept baryons
  Bool_t ifBar = kFALSE;
  Int_t pdg = trackMCgen->Particle()->GetPdgCode();
  //
  //Check if the particle is in the baryon list
  for (Int_t ibar=0;ibar<fNBarBins;ibar++){
    if ( fBaryons[ibar] == pdg ) {
      ifBar = kTRUE;
      break;
    }
  }
  return ifBar;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskTIdentityPID::CheckIfFromAnyResonance(AliMCParticle *trackMCgen, Float_t etaLow, Float_t etaUp, Float_t pDown, Float_t pUp)
{

  Bool_t sisterInAcceptance = kTRUE;
  Bool_t motherInAcceptance = kTRUE;

  Int_t labMom = trackMCgen->GetMother();
  Int_t pdgMom = 0;
  if ((labMom>=0) && (labMom < fMCEvent->GetNumberOfTracks())){
    //
    // Check if the mother is also in the acceptance
    AliMCParticle *momTrack = (AliMCParticle *)fMCEvent->GetTrack(labMom);
    pdgMom = momTrack->Particle()->GetPdgCode();
    if ( pdgMom!=0 ) {
      Int_t    motherSign = momTrack->Charge();
      Double_t motherEta  = (fRapidityType==0) ? momTrack->Eta() :  momTrack->Y();
      Double_t motherMom  = momTrack->P();
      Bool_t etaAcc  = (motherEta>=etaLow && motherEta<etaUp);
      Bool_t momAcc  = (motherMom>=pDown  && motherMom<pUp  );
      if ( !(etaAcc && momAcc) ) motherInAcceptance=kFALSE;
    }
    //
    // Check if the sister is also in the acceptance
    Int_t labSister = momTrack->GetDaughterLast();
    if ((labSister>=0) && (labSister < fMCEvent->GetNumberOfTracks())){
      AliMCParticle *sisterTrack = (AliMCParticle *)fMCEvent->GetTrack(labSister);
      Int_t    sisterSign = sisterTrack->Charge();
      Double_t sisterEta  = (fRapidityType==0) ? sisterTrack->Eta() :  sisterTrack->Y();
      Double_t sisterMom  = sisterTrack->P();
      Bool_t etaAcc  = (sisterEta>=etaLow && sisterEta<etaUp);
      Bool_t momAcc  = (sisterMom>=pDown  && sisterMom<pUp  );
      if ( !(sisterSign!=0 && etaAcc && momAcc) ) sisterInAcceptance=kFALSE;
    }
  }
  // default is accept resonances
  Bool_t acceptRes = kTRUE;
  if ( pdgMom!=0 &&  fSisterCheck==0 ) acceptRes = kFALSE; // in anycase if mother exist reject particle
  if ( pdgMom!=0 &&  motherInAcceptance &&  sisterInAcceptance && fSisterCheck==1 ) acceptRes = kFALSE; // if sister and mother are in acc reject particle
  if ( pdgMom!=0 &&  motherInAcceptance && !sisterInAcceptance && fSisterCheck==2 ) acceptRes = kFALSE; // if sister and mother are in acc reject particle
  if ( pdgMom!=0 &&  motherInAcceptance && !sisterInAcceptance && fSisterCheck==3 ) acceptRes = kTRUE; // if sister and mother are in acc accept particle
  //
  if ( pdgMom!=0 && !motherInAcceptance && !sisterInAcceptance && fSisterCheck==4 ) acceptRes = kFALSE; // if sister and mother are in acc accept particle
  if ( pdgMom!=0 && !motherInAcceptance &&  sisterInAcceptance && fSisterCheck==5 ) acceptRes = kFALSE; // if sister and mother are in acc accept particle
  if ( pdgMom!=0 && !motherInAcceptance &&  fSisterCheck==6 ) acceptRes = kFALSE;    // if sister and mother are in acc accept particle
  return acceptRes;

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::SelectCleanSamplesFromV0s(AliESDv0 *v0, AliESDtrack *track0, AliESDtrack *track1)
{
  //
  // SetAliases and Metadata for the V0 trees
  //
  AliKFParticle kfparticle; //
  AliAnalysisTaskFilteredTree filteredV0;
  filteredV0.GetKFParticle(v0,fESD,kfparticle);
  //
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  Double_t massLambda = pdg->GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg->GetParticle("K0")->Mass();
  const Double_t massProton  =pdg->GetParticle(kPDGpr)->Mass();
  const Double_t massPion    =pdg->GetParticle(kPDGpi)->Mass();
  const Double_t livetimeK0=2.684341668932;  // livetime in cm (surpisely missing info in PDG - see root forum)
  const Double_t livetimeLambda=7.8875395;  // livetime in cm (missing info in PDG - see root forum)
  fHasTrack0FirstITSlayer = track0->HasPointOnITSLayer(0);
  fHasTrack1FirstITSlayer = track1->HasPointOnITSLayer(0);
  fHasV0FirstITSlayer = (fHasTrack0FirstITSlayer||fHasTrack1FirstITSlayer);
  //
  Double_t v0Rr = v0->GetRr();  // rec position of the vertex CKBrev
  Double_t v0P  = v0->P();      // TMath::Sqrt(Px()*Px()+Py()*Py()+Pz()*Pz())
  Double_t livetimeLikeK0 = TMath::Exp(-v0Rr/(TMath::Sqrt((v0P/massK0)*(v0P/massK0)+1)*livetimeK0));
  Double_t livetimeLikeLambda = TMath::Exp(-v0Rr/(TMath::Sqrt((v0P/massLambda)*(v0P/massLambda)+1)*livetimeLambda));
  Double_t livetimeLikeGamma = v0Rr/80.;
  // Double_t livetimeLikeBkg   = v0Rr/80.;

  // delta of mass
  Double_t K0Delta = v0->GetEffMass(2,2)-massK0;        //   tree->SetAlias("K0Delta","(v0.GetEffMass(2,2)-massK0)");
  Double_t LDelta  = v0->GetEffMass(4,2)-massLambda;    //   tree->SetAlias("LDelta","(v0.GetEffMass(4,2)-massLambda)");
  Double_t ALDelta = v0->GetEffMass(2,4)-massLambda;    //   tree->SetAlias("ALDelta","(v0.GetEffMass(2,4)-massLambda)");
  Double_t EDelta  = v0->GetEffMass(0,0);               //   tree->SetAlias("EDelta","(v0.GetEffMass(0,0))");

  // pull of the mass
  if (v0->GetKFInfo(2,2,1)==0. || v0->GetKFInfo(4,2,1)==0. || v0->GetKFInfo(2,4,1)==0. || v0->GetKFInfo(0,0,1)==0.) return;
  Double_t K0Pull = (v0->GetEffMass(2,2)-massK0)/v0->GetKFInfo(2,2,1);        //   tree->SetAlias("K0Pull","(v0.GetEffMass(2,2)-massK0)/v0.GetKFInfo(2,2,1)");
  Double_t LPull  = (v0->GetEffMass(4,2)-massLambda)/v0->GetKFInfo(4,2,1);    //   tree->SetAlias("LPull","(v0.GetEffMass(4,2)-massLambda)/v0.GetKFInfo(4,2,1)");
  Double_t ALPull = (v0->GetEffMass(2,4)-massLambda)/v0->GetKFInfo(2,4,1);    //   tree->SetAlias("ALPull","(v0.GetEffMass(2,4)-massLambda)/v0.GetKFInfo(2,4,1)");
  Double_t EPull  = EDelta/v0->GetKFInfo(0,0,1);                              //   tree->SetAlias("EPull","EDelta/v0.GetKFInfo(0,0,1)");

  Double_t K0Like0 = TMath::Exp(-K0Pull*K0Pull)*livetimeLikeK0;
  Double_t LLike0  = TMath::Exp(-LPull*LPull)*livetimeLikeLambda;
  Double_t ALLike0 = TMath::Exp(-ALPull*ALPull)*livetimeLikeLambda;
  Double_t ELike0  = TMath::Exp(-abs(EPull)*0.2)*livetimeLikeGamma;
  Double_t V0Like  = TMath::Exp(-TMath::ACos(v0->GetV0CosineOfPointingAngle())*v0Rr/0.36)*TMath::Exp(-TMath::Sqrt(kfparticle.GetChi2())/0.5);
  //
  Int_t ntracks = fESD->GetNumberOfTracks();
  Double_t BkgLike = 0.000005*ntracks;    // backround coeefecint  to be fitted - depends on other cuts
  Double_t LLike = LLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);
  Double_t ALLike = ALLike0/(K0Like0+LLike0+ALLike0+ELike0+BkgLike);
  //
  Double_t tr0NTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track0, AliPID::kPion));
  Double_t tr1NTPCSigmaPi = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion));
  Float_t ptotForBetaGamma0 = track0->GetInnerParam()->GetP();
  Float_t ptotForBetaGamma1 = track1->GetInnerParam()->GetP();
  Float_t ptotForBetaGammaThr = 0.2;
  Double_t tr0NTPCSigmaPr = (ptotForBetaGamma0>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track0, AliPID::kProton)) : 0.;
  Double_t tr1NTPCSigmaPr = (ptotForBetaGamma1>ptotForBetaGammaThr) ? TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kProton)) : 0.;
  //
  fCleanPion0FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPi<2 && (v0->PtArmV0()>0.06);
  fCleanPion1FromK0 = (K0Like0>0.05) && (K0Like0>(LLike0+ALLike0+ELike0)*3) && (TMath::Abs(K0Delta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPi<2 && (v0->PtArmV0()>0.06);
  fCleanPion0FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPr<2;
  fCleanPion1FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPr<2;
  fCleanProton0FromLambda = (LLike>0.05) && (LLike0>(K0Like0+ALLike0+ELike0)*3) && (TMath::Abs(LDelta)<0.006) && (V0Like>0.1) && tr1NTPCSigmaPi<2;
  fCleanProton1FromLambda = (ALLike>0.05) && (ALLike0>(K0Like0+LLike0+ELike0)*3) && (TMath::Abs(ALDelta)<0.006) && (V0Like>0.1) && tr0NTPCSigmaPi<2;
  fCleanPionsFromK0 =  (fCleanPion0FromK0 || fCleanPion1FromK0);

}
//________________________________________________________________________
Bool_t AliAnalysisTaskTIdentityPID::ApplyDCAcutIfNoITSPixel(AliESDtrack *track)
{

  Float_t p[2],cov[3];
  track->GetImpactParameters(p,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  Bool_t isFirstITSlayer  = track->HasPointOnITSLayer(0);
  Bool_t isSecondITSlayer = track->HasPointOnITSLayer(1);

  fIsITSpixel01 = (isFirstITSlayer || isSecondITSlayer);
  fNITSclusters = track->GetNumberOfITSClusters();

  if (!cov[0] || !cov[2]) {
    return kFALSE;
  } else {
    fPrimRestriction = TMath::Sqrt((p[0]*p[0])/cov[0] + (p[1]*p[1])/cov[2]);
    return (fPrimRestriction<2 && fNITSclusters>2) || (fPrimRestriction<5 && fIsITSpixel01);
  }

}
//________________________________________________________________________
UInt_t AliAnalysisTaskTIdentityPID::SetCutBitsAndSomeTrackVariables(AliESDtrack *track, Int_t particleType)
{
  //
  // Set some track variables
  //
  //
  // --------------------------------------------------------------
  //  calculate some variables by hand
  // --------------------------------------------------------------
  //
  // Double_t p[3];
  // track->GetPxPyPz(p);
  // Double_t momentum = TMath::Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
  // Double_t pt       = TMath::Sqrt(p[0]*p[0] + p[1]*p[1]);
  // Double_t mass     = track->GetMass();  // assumed to be pion mass
  // Double_t energy   = TMath::Sqrt(mass*mass + momentum*momentum);
  // Float_t eta = -100.;
  // Float_t rap   = -100.;
  // if((momentum != TMath::Abs(p[2]))&&(momentum != 0)) eta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));
  // if((energy != TMath::Abs(p[2]))&&(energy != 0))     rap = 0.5*TMath::Log((energy + p[2])/(energy - p[2]));
  //
  // --------------------------------------------------------------
  //  some extra getters
  // --------------------------------------------------------------
  //
  fDefaultCuts = fESDtrackCuts->AcceptTrack(track);
  Double_t goldenChi2 = 0.;
  fPVertex = track->P();
  fTheta=track->Theta();
  fSign= track->GetSign();
  fPx  =track->Px();
  fPy  =track->Py();
  fPz  =track->Pz();
  fPt  =track->Pt();
  fY   =track->Y();
  fPhi =track->Phi()-TMath::Pi();
  fEta = track->Eta();
  Float_t pv[2],cov[3];
  track->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  fTrackDCAxy = pv[0];
  fTrackDCAz  = pv[1];

  //
  // TPC related quantities
  Bool_t cleanDeTPC = kFALSE;
  if (track->GetInnerParam()){
    fPtot      = track->GetInnerParam()->GetP();
    fTPCSignal = track->GetTPCsignal();
    fTrackTPCSignalN     = track->GetTPCsignalN();
    fTrackTPCCrossedRows = Float_t(track->GetTPCCrossedRows());
    fTPCShared = track->GetTPCnclsS();
    fTPCFindable = track->GetTPCNclsF();
    fMissingCl = track->GetTPCClusterInfo(3,0,0,159);
    goldenChi2 = track->GetChi2TPCConstrainedVsGlobal(fVertex);
    fTrackLengthInActiveZone = track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0);
    //
    Float_t ptotForBetaGamma = track->GetInnerParam()->GetP();
    Float_t ptotForBetaGammaThr = 0.2;
    Double_t nSigmasDeTPC = (ptotForBetaGamma>ptotForBetaGammaThr) ? fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron) : 0.;
    cleanDeTPC = ((TMath::Abs(nSigmasDeTPC)<=2.));
    fNcl       = track->GetTPCncls();
    // fTrackChi2TPC  = (fNcl>0) ? TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/fNcl)) : -1;  // ???
    fTrackChi2TPC  = (fNcl>0) ? TMath::Abs(track->GetTPCchi2()/fNcl) : -1;  // ???
    //
    // correct for the missing clusters
    fNclCorr = fNcl;
    fTrackChi2TPCcorr = fTrackChi2TPC;
    if (fCorrectForMissCl==2){
      if ( (fPtot>0.2 && fPtot<3.2) && TMath::Abs(fEta)<0.8 && fCentrality<5 && particleType<4 ){
        Float_t missCl_binContent = 0.;
        Int_t multBin = 0;
        if (fTPCMult>=6000  && fTPCMult<6500)  multBin=0;
        if (fTPCMult>=6500  && fTPCMult<7000)  multBin=1;
        if (fTPCMult>=7000  && fTPCMult<7500)  multBin=2;
        if (fTPCMult>=7500  && fTPCMult<8000)  multBin=3;
        if (fTPCMult>=8000  && fTPCMult<8500)  multBin=4;
        if (fTPCMult>=8500  && fTPCMult<9000)  multBin=5;
        if (fTPCMult>=9000  && fTPCMult<9500)  multBin=6;
        if (fTPCMult>=9500  && fTPCMult<10000) multBin=7;
        if (fTPCMult>=10000 && fTPCMult<10500) multBin=8;
        if (fTPCMult>=10500 && fTPCMult<11000) multBin=9;
        if (fTPCMult>=11000 && fTPCMult<11500) multBin=10;
        if (fTPCMult>=11500 && fTPCMult<12000) multBin=11;
        if (fTPCMult>=12000 && fTPCMult<12500) multBin=12;
        if (fTPCMult>=12500 && fTPCMult<13000) multBin=13;
        if (fTPCMult>=13000 && fTPCMult<13500) multBin=14;
        if (fTPCMult>=13500 && fTPCMult<14000) multBin=15;
        missCl_binContent = fH2MissCl[particleType][multBin].GetBinContent(fH2MissCl[particleType][multBin].FindBin(fPtot,fEta));
        if(missCl_binContent<30 && missCl_binContent>0) {
          fNclCorr = fNcl-(missCl_binContent/100.)*fNcl;
          fTrackChi2TPCcorr  = (fNclCorr>0) ? TMath::Abs(track->GetTPCchi2()/fNclCorr) : -1;  // ???
        }
      }
    }
    //
    // --------------------------------------------------------------
    //      Bayesian PID part
    // --------------------------------------------------------------
    //
    if (ptotForBetaGamma>ptotForBetaGammaThr) {
      fPIDCombined->SetDefaultTPCPriors();
      Double_t probTPC[AliPID::kSPECIES]={0.};
      Double_t probTOF[AliPID::kSPECIES]={0.};
      // Get TPC probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPC);
      fTrackProbPiTPC = probTPC[AliPID::kPion];
      fTrackProbKaTPC = probTPC[AliPID::kKaon];
      fTrackProbPrTPC = probTPC[AliPID::kProton];
      // fTrackProbDeTPC = probTPC[AliPID::kDeuteron];
      // Get TOF probabilities
      fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);
      fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTOF);
      fTrackProbPiTOF = probTOF[AliPID::kPion];
      fTrackProbKaTOF = probTOF[AliPID::kKaon];
      fTrackProbPrTOF = probTOF[AliPID::kProton];
    }
  }
  //
  // its restrictions
  Bool_t isFirstITSlayer  = track->HasPointOnITSLayer(0);
  Bool_t isSecondITSlayer = track->HasPointOnITSLayer(1);
  fIsITSpixel01    = (isFirstITSlayer || isSecondITSlayer);
  fNITSclusters    = track->GetNumberOfITSClusters();
  if (cov[0]>0 && cov[1]>0){
    fPrimRestriction = TMath::Sqrt((pv[0]*pv[0])/cov[0] + (pv[1]*pv[1])/cov[2]);
  }
  //
  fTrackRequireITSRefit  = track->IsOn(AliESDtrack::kITSrefit); // track->IsOn(AliESDtrack::kTPCrefit);
  fTrackIsFirstITSlayer  = track->HasPointOnITSLayer(0);
  fTrackIsSecondITSlayer = track->HasPointOnITSLayer(1);
  fTrackNewITScut        = ApplyDCAcutIfNoITSPixel(track);
  //
  Double_t nclsTRD     = (Float_t)track->GetTRDncls();
  Double_t TOFSignalDx = track->GetTOFsignalDx();
  Double_t TOFSignalDz = track->GetTOFsignalDz();
  //
  //
  fNSigmasPiTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion,    fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasKaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon,    fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasPrTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton,  fPIDResponse->GetTOFResponse().GetTimeZero());
  fNSigmasDeTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron,fPIDResponse->GetTOFResponse().GetTimeZero());
  //
  vector<Double_t> tempNSigma = GetNSigmas(kCutReference);
  Double_t nSigmaTOFDown = tempNSigma[2];
  Double_t nSigmaTOFUp   = tempNSigma[3];
  //
  Bool_t cleanPiTOF = nSigmaTOFDown < fNSigmasPiTOF && fNSigmasPiTOF <= nSigmaTOFUp;
  Bool_t cleanKaTOF = nSigmaTOFDown < fNSigmasKaTOF && fNSigmasKaTOF <= nSigmaTOFUp;
  Bool_t cleanPrTOF = nSigmaTOFDown < fNSigmasPrTOF && fNSigmasPrTOF <= nSigmaTOFUp;
  Bool_t cleanDeTOF = nSigmaTOFDown < fNSigmasDeTOF && fNSigmasDeTOF <= nSigmaTOFUp;
  //
  Bool_t cleanKaTOFTRD = ((TMath::Abs(fNSigmasKaTOF)<=1.2) && TOFSignalDz<1. && TOFSignalDx<1. && nclsTRD>100);
  //
  tempNSigma = GetNSigmas(kCutNSigmaTOFLoose);
  Double_t nSigmaTPCDownLoose = tempNSigma[0];
  Double_t nSigmaTPCUpLoose   = tempNSigma[1];
  Double_t nSigmaTOFDownLoose = tempNSigma[2];
  Double_t nSigmaTOFUpLoose   = tempNSigma[3];
  Bool_t prTPCLoose = nSigmaTPCDownLoose < fNSigmasPrTPC && fNSigmasPrTPC <= nSigmaTPCUpLoose;
  Bool_t prTOFLoose = nSigmaTOFDownLoose < fNSigmasPrTOF && fNSigmasPrTOF <= nSigmaTOFUpLoose;
  Bool_t prPIDLoose = (fPtot < fTOFMomCut && prTPCLoose) || (fPtot >= fTOFMomCut && prTOFLoose && prTPCLoose);
  //
  tempNSigma = GetNSigmas(kCutNSigmaTOFLoose2);
  Double_t nSigmaTPCDownLoose2 = tempNSigma[0];
  Double_t nSigmaTPCUpLoose2   = tempNSigma[1];
  Double_t nSigmaTOFDownLoose2 = tempNSigma[2];
  Double_t nSigmaTOFUpLoose2   = tempNSigma[3];
  Bool_t prTPCLoose2 = nSigmaTPCDownLoose2 < fNSigmasPrTPC && fNSigmasPrTPC <= nSigmaTPCUpLoose2;
  Bool_t prTOFLoose2 = nSigmaTOFDownLoose2 < fNSigmasPrTOF && fNSigmasPrTOF <= nSigmaTOFUpLoose2;
  Bool_t prPIDLoose2 = (fPtot < fTOFMomCut && prTPCLoose2) || (fPtot >= fTOFMomCut && prTOFLoose2 && prTPCLoose2);
  //
  Bool_t dcaBaseCut = TMath::Abs(fTrackDCAxy)<0.0208+0.04/TMath::Power(fPt,1.01);
  Bool_t dcaLoose   = TMath::Abs(fTrackDCAxy)<0.4;  // 10h tuned loose cut
  // pileup
  fIsMCPileup = kFALSE;
  if (fMCtrue) {
    fIsMCPileup = IsFromPileup(track->GetLabel());
  }
  //
  // Systematic settings
  fTrackCutBits=0;
  //
  // Crossed rows
  if (fCorrectForMissCl==1){
    if (fNcl>=70) (fTrackCutBits |= 1 << kNCrossedRowsTPC70);
    if (fNcl>=80) (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (fNcl>=90) (fTrackCutBits |= 1 << kNCrossedRowsTPC90);
  } else if (fCorrectForMissCl==2){
    if (fNclCorr>=70) (fTrackCutBits |= 1 << kNCrossedRowsTPC70);
    if (fNclCorr>=80) (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (fNclCorr>=90) (fTrackCutBits |= 1 << kNCrossedRowsTPC90);
    // cout <<  "ncls  = " << fNclCorr <<  " --- " << fNcl << endl;
  } else {
    if (fTrackTPCCrossedRows>=70) (fTrackCutBits |= 1 << kNCrossedRowsTPC70);
    if (fTrackTPCCrossedRows>=80) (fTrackCutBits |= 1 << kNCrossedRowsTPC80);
    if (fTrackTPCCrossedRows>=90) (fTrackCutBits |= 1 << kNCrossedRowsTPC90);
  }
  //
  // Special treatment of the 2018 pass3 and 2015 pass2 data
  // Chi2 TPC
  if ( (fYear==2015&&fPassIndex==2) || (fYear==2018&&fPassIndex==3) ){
    if (fTrackChi2TPC<2.2) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);
    if (fTrackChi2TPC<2.5) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);
    if (fTrackChi2TPC<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);
  } else {
    //
    // correction for missing clusters
    if (fCorrectForMissCl==2){
      if (fTrackChi2TPCcorr<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);   // ????
      if (fTrackChi2TPCcorr<4.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);        // ????
      if (fTrackChi2TPCcorr<5.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);   // ????
    } else {
      if (fTrackChi2TPC<3.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCSmall);   // ????
      if (fTrackChi2TPC<4.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPC);        // ????
      if (fTrackChi2TPC<5.0) (fTrackCutBits |= 1 << kMaxChi2PerClusterTPCLarge);   // ????
    }
  }
  //
  // Shared TPC clusters
  Bool_t sharedCls = kFALSE;
  Bool_t sharedClsLoose = kFALSE;
  if (fTrackTPCCrossedRows > 0 && fNcl > 0) {
    sharedCls = (fTPCShared / fTrackTPCCrossedRows < 0.25) && (fTPCShared / static_cast<Float_t>(fNcl) < 0.3);
    sharedClsLoose = fTPCShared / fTrackTPCCrossedRows < 0.25;
  }
  if (sharedCls) (fTrackCutBits |= 1 << kSharedCls);
  if (sharedClsLoose) (fTrackCutBits |= 1 << kSharedClsLoose);
  //
  // Found TPC clusters
  if (fTPCFindable > 0) {
    if (fTrackTPCCrossedRows / fTPCFindable > 0.80) (fTrackCutBits |= 1 << kFindableCls);
    if (fTrackTPCCrossedRows / fTPCFindable > 0.85) (fTrackCutBits |= 1 << kFindableClsTight);
    if (fTrackTPCCrossedRows / fTPCFindable > 0.75) (fTrackCutBits |= 1 << kFindableClsLoose);
  }
  //
  // DCAxy
  if (dcaBaseCut) (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDep);
  if (dcaLoose)   (fTrackCutBits |= 1 << kMaxDCAToVertexXYPtDepLarge);
  //
  // DCAz
  if (TMath::Abs(fTrackDCAz)<0.15) (fTrackCutBits |= 1 << kVertexZSmall);
  if (TMath::Abs(fTrackDCAz)<1.00) (fTrackCutBits |= 1 << kVertexZ);
  //
  // Event vertex z
  if (TMath::Abs(fVz)<7 && TMath::Abs(fVz)>0.15) (fTrackCutBits |= 1 << kEventVertexZ);
  if (TMath::Abs(fVz)<8 && TMath::Abs(fVz)>0.1 ) (fTrackCutBits |= 1 << kEventVertexZLarge);
  //
  // NCl in dEdx calculation
  if (fTrackTPCSignalN>=60) (fTrackCutBits |= 1 << kTPCSignalNSmall);
  if (fTrackTPCSignalN>=70) (fTrackCutBits |= 1 << kTPCSignalN);
  if (fTrackTPCSignalN>=80) (fTrackCutBits |= 1 << kTPCSignalNLarge);
  //
  // pile-up
  if (!fMCtrue) { // real data
    if (fPileUpBit & 1 << 0) (fTrackCutBits |= 1 << kPileup);
    if (fPileUpBit & 1 << 1) (fTrackCutBits |= 1 << kPileupLoose);
  } else {
    if (fCollisionType == 0) { // PbPb
      if (!fIsMCPileup) (fTrackCutBits |= 1 << kPileup);
      fTrackCutBits |= 1 << kPileupLoose; // fill for all events if no pileup rejection
    } else if (fCollisionType == 1) { // pp
      fTrackCutBits |= 1 << kPileup; // no pileup in pp mc, so we always set this bit
    }
  }
  //
  // B field polarity
  if (fBField > 0) (fTrackCutBits |= 1 << kBFieldPos);
  if (fBField < 0) (fTrackCutBits |= 1 << kBFieldNeg);
  //
  // --------------------------------------------------------------------
  //                    Clean sample selections
  // --------------------------------------------------------------------
  //
  // variable nsigma TOF protons and kaons for amplitude estimation
  if (cleanPiTOF) (fTrackCutBits |= 1 << kCleanPiTOF);
  if (cleanPrTOF) (fTrackCutBits |= 1 << kCleanPrTOF);
  if (cleanKaTOF) (fTrackCutBits |= 1 << kCleanKaTOF);
  //
  // Clean Kaons protons and deuterons
  if (cleanKaTOFTRD)        (fTrackCutBits |= 1 << kCleanKaTOFTRD);
  if (fTrackProbKaTOF>=0.8) (fTrackCutBits |= 1 << kTrackProbKaTOF);
  // if (fTrackProbPrTOF>=0.8) (fTrackCutBits |= 1 << kTrackProbPrTOF);
  if (cleanDeTOF && cleanDeTPC) (fTrackCutBits |= 1 << kCleanDeTOF);
  //
  if (prPIDLoose)  (fTrackCutBits |= 1 << kNSigmaTOFLoose);
  if (prPIDLoose2) (fTrackCutBits |= 1 << kNSigmaTOFLoose2);
  //
  return fTrackCutBits;

}
//________________________________________________________________________
Bool_t AliAnalysisTaskTIdentityPID::GetSystematicClassIndex(UInt_t cut,Int_t syst)
{
  /*
  syst:
  0 -->  Reference
  1 -->  CRows70
  2 -->  CRows90
  3 -->  ActiveZone
  4 -->  Chi2TPCSmall
  5 -->  Chi2TPCLarge
  6 -->  kMaxDCAToVertexXYPtDepLarge
  7 -->  kVertexZSmall
  8 -->  kEventVertexZLarge
  9 -->  kSharedCls
  10 -->  kFindableClsTight
  11 -->  kFindableClsLoose
  12 -->  kPileupLoose
  13 -->  kBFieldPos
  14 -->  kBFieldNeg
  15 -->  kTPCSignalNSmall
  16 -->  kTPCSignalNLarge
  */

  std::vector<Int_t> fCutArr;

  switch(syst) {

    case kCutReference:   // 0 -->  Reference
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutCrossedRowsTPC70:  // 1 -->  kNCrossedRowsTPC70
    {
      fCutArr = {kNCrossedRowsTPC70,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutCrossedRowsTPC90:  // 2 -->  kNCrossedRowsTPC90
    {
      fCutArr = {kNCrossedRowsTPC90,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCSmall:   // 3 -->  kMaxChi2PerClusterTPCSmall
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCSmall, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCLarge:   // 4 -->  kMaxChi2PerClusterTPCLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCLarge, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutMaxDCAToVertexXYPtDepLarge:   // 5 -->  kMaxDCAToVertexXYPtDepLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDepLarge, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutVertexZSmall:   // 6 -->  kVertexZSmall
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZSmall, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutEventVertexZLarge:  // 7 -->  kEventVertexZLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZLarge, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutSharedCls:   // 8 -->  kSharedClsLoose
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedClsLoose, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutFindableClsTight:   // 9 -->  kFindableClsTight
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableClsTight,kTPCSignalN};
    }
    break;
    //
    case kCutFindableClsLoose:   // 10 -->  kFindableClsLoose
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableClsLoose,kTPCSignalN};
    }
    break;
    //
    case kCutPileupLoose:   // 11 -->  kPileupLoose
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileupLoose, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutBFieldPos:   // 12 -->  kBFieldPos
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN,kBFieldPos};
    }
    break;
    //
    case kCutBFieldNeg:   // 13 --> kBFieldNeg
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN,kBFieldNeg};
    }
    break;
    //
    case kCutTPCSignalNSmall:   // 14 --> kTPCSignalNSmall
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalNSmall};
    }
    break;
    //
    case kCutTPCSignalNLarge:   // 15 --> kTPCSignalNLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalNLarge};
    }
    break;
    //
    case kCutNSigmaTOFLoose:   // 16 --> kNSigmaTOFLoose
    {
      // default cuts, different pid cuts are applied in settings or in post processing where necessary
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutNSigmaTOFLoose2:   // 17 --> kNSigmaTOFLoose2
    {
      // default cuts, different pid cuts are applied in settings or in post processing where necessary
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls, kTPCSignalN};
    }
    break;
    //
    default:
    {
      fCutArr = {};
    }

  }
  //
  //  Apply conditions
  for (UInt_t i=0;i<fCutArr.size();i++){
    if( ((cut >> fCutArr[i]) & 1) == 0 )
    {
      return kFALSE;
    }
  }
  return kTRUE;

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::SetSpecialV0Cuts(AliESDv0KineCuts* cuts)
{

  cuts->SetMode(0, 0);    // cuts->SetMode(mode, type); mode 0: purely kinematical selection   type 0: pp 1:PbPb not yet ready
  //
  // leg cuts
  cuts->SetNTPCclusters(50);
  cuts->SetTPCrefit(kTRUE);
  cuts->SetTPCchi2perCls(4.0);
  cuts->SetTPCclusterratio(0.6);
  cuts->SetNoKinks(kTRUE);
  //
  // gamma cuts
  cuts->SetGammaCutChi2NDF(10.0);
  Float_t cosPoint[2] = {0.0, 0.02};
  cuts->SetGammaCutCosPoint(cosPoint);
  Float_t cutDCA[2] = {0.0, 0.25};
  cuts->SetGammaCutDCA(cutDCA);
  Float_t vtxR[2] = {3.0, 90.0};
  cuts->SetGammaCutVertexR(vtxR);
  Float_t psiPairCut[2]={0.0,0.05};
  cuts->SetGammaCutPsiPair(psiPairCut);
  cuts->SetGammaCutInvMass(0.05);
  // K0s cuts
  cuts->SetK0CutChi2NDF(10.0);
  Float_t cosPointK0s[2] = {0.0, 0.02};
  cuts->SetK0CutCosPoint(cosPointK0s);
  Float_t cutDCAK0s[2] = {0.0, 0.2};
  cuts->SetK0CutDCA(cutDCAK0s);
  Float_t vtxRK0s[2] = {2.0, 30.0};
  cuts->SetK0CutVertexR(vtxRK0s);
  Float_t k0sInvMass[2] = {0.486, 0.508};
  cuts->SetK0CutInvMass(k0sInvMass);
  // Lambda and anti-Lambda cuts
  cuts->SetLambdaCutChi2NDF(10.0);
  Float_t cosPointLambda[2] = {0.0, 0.02};
  cuts->SetLambdaCutCosPoint(cosPointLambda);
  Float_t cutDCALambda[2] = {0.0, 0.2};
  cuts->SetLambdaCutDCA(cutDCALambda);
  Float_t vtxRLambda[2] = {2.0, 40.0};
  cuts->SetLambdaCutVertexR(vtxRLambda);
  Float_t lambdaInvMass[2] = {1.11, 1.12};
  cuts->SetLambdaCutInvMass(lambdaInvMass);

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::BinLogAxis(TH1 *h)
{
  //
  // Method for the correct logarithmic binning of histograms
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the BinLogAxis ===== " << std::endl;
  TAxis *axis       = h->GetXaxis();
  Int_t bins        = axis->GetNbins();

  Double_t from     = axis->GetXmin();
  Double_t to       = axis->GetXmax();
  std::vector<double>  newBins;
  newBins.resize(bins + 1);

  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);

  for (Int_t i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins.data());

}
//________________________________________________________________________
Int_t AliAnalysisTaskTIdentityPID::CountEmptyEvents(Int_t counterBin, Int_t setting)
{

  //
  // Detect Empty Events
  Int_t trackcounter = 0;
  if (setting == -1){ // for MC gen level track counting
    if (fUseCouts) std::cout << " Info::marsland: ===== In the CountEmptyEvents from MC gen info ===== " << std::endl;
    for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
    {
      AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
      if (!trackMCgen) continue;
      if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
      Float_t pTMCgen   = trackMCgen->Pt();
      Float_t etaMCgen  = trackMCgen->Eta();
      Bool_t momAcceptance = (pTMCgen > 0.15 && pTMCgen < 10.);
      Bool_t etaAcceptance = (etaMCgen > -0.8 && etaMCgen < 0.8);
      if (momAcceptance && etaAcceptance) trackcounter++;
    }
  } else { // for data event plane
    if (fUseCouts) std::cout << " Info::marsland: ===== In the CountEmptyEvents from ESD tracks ===== " << std::endl;
    for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
      AliESDtrack* track = fESD->GetTrack(iTrack);
      if (!fESDtrackCuts->AcceptTrack(track)) continue;
      if (!track->GetInnerParam()) continue;
      if (!(track->GetTPCsignalN()>0)) continue;
      SetCutBitsAndSomeTrackVariables(track,0);
      Bool_t momAcceptance = (fPt > 0.15 && fPt < 10.);
      Bool_t etaAcceptance = (fEta > -0.8 && fEta < 0.8);
      if (momAcceptance && etaAcceptance) trackcounter++;
    }
  }
  //
  // check if the event is empty
  if (trackcounter<1) {
    fHistEmptyEvent->Fill(counterBin);
    std::cout << " Info::marsland: AAAAAA Empty event in " << fChunkName << std::endl;
  }
  return trackcounter;

}
//
void AliAnalysisTaskTIdentityPID::PrintNumInBinary(UInt_t num)
{
  TString bin="";
  Int_t numberOfBits = sizeof(UInt_t)*8;
  for (Int_t i=numberOfBits-1; i>=0; i--) {
    Bool_t isBitSet = (num & (1<<i));
    if (isBitSet) {
      bin+="1";
    } else {
      bin+="0";
    }
  }
  std::cout << "Info::marsland: fTrackCutBits = " << bin << std::endl;
}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  std::cout << " Info::marsland: ===== In the Terminate ===== " << std::endl;

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::CreateEventInfoTree()
{

  if (fUseCouts) std::cout << " Info::marsland: ===== In the CreateEventInfoTree ===== " << std::endl;
  //
  MakeEventPlane(1, 1, 0);
  if (fCollisionType == 1) {
    GetFlatenicity();
    fSpherocity = (fMCEvent) ? ComputeSpherocity(-1): ComputeSpherocity(0);
  }

  Int_t tpcClusterMultiplicity   = fESD->GetNumberOfTPCClusters();
  const AliMultiplicity *multObj = fESD->GetMultiplicity();
  Int_t itsNumberOfTracklets   = multObj->GetNumberOfTracklets();

  TVectorF phiCountA(36);
  TVectorF phiCountC(36);
  TVectorF phiCountAITS(36);
  TVectorF phiCountCITS(36);
  TVectorF phiCountAITSonly(36);
  TVectorF phiCountCITSonly(36);
  TVectorF tzeroMult(24);  for (Int_t i=1;i<24;i++) tzeroMult[i] = 0.;
  TVectorF vzeroMult(64);  for (Int_t i=1;i<64;i++) vzeroMult[i] = 0.;
  TVectorF itsClustersPerLayer(6); for (Int_t i=1;i<6;i++) itsClustersPerLayer[i] = 0.;
  //
  for (Int_t i=1;i<37;i++){
    phiCountA[i-1] = fHistPhiTPCcounterA->GetBinContent(i);
    phiCountC[i-1] = fHistPhiTPCcounterC->GetBinContent(i);
    phiCountAITS[i-1] = fHistPhiTPCcounterAITS->GetBinContent(i);
    phiCountCITS[i-1] = fHistPhiTPCcounterCITS->GetBinContent(i);
    phiCountAITSonly[i-1] = fHistPhiITScounterA->GetBinContent(i);
    phiCountCITSonly[i-1] = fHistPhiITScounterC->GetBinContent(i);
  }
  //
  // Additional counters for ITS TPC V0 and T0
  AliESDFMD* esdFMD = fESD->GetFMDData();
  AliFMDFloatMap fmdMult = esdFMD->MultiplicityMap();
  const AliESDTZERO *esdTzero = fESD->GetESDTZERO();
  const Double32_t *t0amp=esdTzero->GetT0amplitude();
  ULong64_t triggerMask = fESD->GetTriggerMask();
  Short_t   eventMult = fESD->GetNumberOfTracks();
  Int_t eventMultESD = fESD->GetNumberOfESDTracks();
  Int_t tpcTrackBeforeClean =fESD->GetNTPCTrackBeforeClean();

  for (Int_t i=0;i<24;i++) { tzeroMult[i] = t0amp[i]; }
  for (Int_t i=0;i<64;i++) { vzeroMult[i] = fESD->GetVZEROData()-> GetMultiplicity(i); }
  for (Int_t i=0;i<6;i++)  { itsClustersPerLayer[i] = multObj->GetNumberOfITSClusters(i); }

  if(!fTreeSRedirector) return;
  DumpDownScaledTree();
  (*fTreeSRedirector)<<"eventInfo"<<
  "gid=" << fEventGID << // global event ID
  "run=" << fRunNo <<  // run Number
  "spher=" << fSpherocity <<
  "flat=" << fFlatenicity <<
  "flatS=" << fFlatenicityScaled <<
  "Qx2_neg=" << fEP_2_Qx_neg <<
  "Qx2_pos=" << fEP_2_Qx_pos <<
  "Qy2_neg=" << fEP_2_Qy_neg <<
  "Qy2_pos=" << fEP_2_Qy_pos <<
  "psi2_pos=" << fEP_2_Psi_pos <<
  "psi2_neg=" << fEP_2_Psi_neg <<
  "psi2=" << fEP_2_Psi <<
  "Qx3_neg=" << fEP_3_Qx_neg <<
  "Qx3_pos=" << fEP_3_Qx_pos <<
  "Qy3_neg=" << fEP_3_Qy_neg <<
  "Qy3_pos=" << fEP_3_Qy_pos <<
  "psi3_pos=" << fEP_3_Psi_pos <<
  "psi3_neg=" << fEP_3_Psi_neg <<
  "psi3=" << fEP_3_Psi <<
  "ntracks_neg=" << fEP_ntracks_neg <<
  "ntracks_pos=" << fEP_ntracks_pos <<
  "pileupbit=" << fPileUpBit <<
  "bField=" << fBField << // run Number
  "timestamp=" << fTimeStamp << // timestamp
  "triggerMask=" << triggerMask << //trigger mask
  "intrate=" << fIntRate << // interaction rate
  "centV0M=" << fCentrality << // centrality
  "centrality.=" << fEventInfo_CentralityEstimates << // track counter
  "vZ=" << fVz <<  // vertex Z
  "tpcvz=" << fTPCvZ <<
  "spdvz=" << fSPDvZ <<
  "tpcMult=" << fTPCMult << // TPC multiplicity
  "eventMult=" << fEventMult << // event multiplicity
  "eventMultESD=" << eventMultESD << // event multiplicity ESD
  "primMult=" << fNContributors << // #prim tracks
  "tpcClusterMult=" << tpcClusterMultiplicity << // tpc cluster multiplicity
  "tpcTrackBeforeClean=" << tpcTrackBeforeClean << // tpc track before cleaning
  "itsTracklets=" << itsNumberOfTracklets << // number of ITS tracklets
  // "fmdMult.=" << &fmdMult << // T0 multiplicity
  "tZeroMult.=" << &tzeroMult << // T0 multiplicity
  "vZeroMult.=" << &vzeroMult <<  // V0 multiplicity
  "itsClustersPerLayer.=" << &itsClustersPerLayer << // its clusters per layer
  "trackCounters.=" << fEventInfo_CacheTrackCounters << // track counter
  "trackdEdxRatio.=" << fEventInfo_CacheTrackdEdxRatio << // dEdx conter
  "trackNcl.=" << fEventInfo_CacheTrackNcl << // nCluster counter
  "trackChi2.=" << fEventInfo_CacheTrackChi2 << // Chi2 counter
  "trackMatchEff.=" << fEventInfo_CacheTrackMatchEff << // Chi2 counter
  "trackTPCCountersZ.=" << fEventInfo_CacheTrackTPCCountersZ << // Chi2 counter
  "hisTPCVertexA.=" << fEventInfo_HisTPCVertexA << // Chi2 counter
  "hisTPCVertexC.=" << fEventInfo_HisTPCVertexC << // Chi2 counter
  "hisTPCVertex.=" << fEventInfo_HisTPCVertex << // Chi2 counter
  "hisTPCVertexACut.=" << fEventInfo_HisTPCVertexACut << // Chi2 counter
  "hisTPCVertexCCut.=" << fEventInfo_HisTPCVertexCCut << // Chi2 counter
  "phiTPCdcarA.=" << fEventInfo_PhiTPCdcarA << // track counter
  "phiTPCdcarC.=" << fEventInfo_PhiTPCdcarC << // dEdx conter
  "phiCountA.=" << &phiCountA << // TPC track count on A side
  "phiCountC.=" << &phiCountC << // TPC track count on C side
  "phiCountAITS.=" << &phiCountAITS << // track count fitted ITS on A side
  "phiCountCITS.=" << &phiCountCITS <<  // track count fitted ITS on C side
  "phiCountAITSOnly.=" << &phiCountAITSonly << // track count only ITS on A side
  "phiCountCITSOnly.=" << &phiCountCITSonly << // track count only ITS on C side
  "\n";

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::CalculateEventInfo()
{

  if (fUseCouts) std::cout << " Info::marsland: ===== In the CalculateEventInfo ===== " << std::endl;

  AliVEvent *event=InputEvent();
  CacheTPCEventInformation();
  //
  //
  const Int_t kNclTPCcut=60;
  const Int_t kDCACut=5;  // 5 cm primary cut
  const Int_t kMindEdxClustersRegion=15;
  const Float_t kTglCut=1.5;
  const Float_t kPtCut=0.100;
  const Float_t kDCAtpcNULL = -10000;
  // const Float_t kNTrackletCut=1.5;
  //
  fEventInfo_PhiTPCdcarA->Zero();
  fEventInfo_PhiTPCdcarC->Zero();
  fEventInfo_CacheTrackCounters->Zero();   // track counter
  fEventInfo_CacheTrackdEdxRatio->Zero(); // dedx info counter
  fEventInfo_CacheTrackNcl->Zero();       // ncl counter
  fEventInfo_CacheTrackChi2->Zero();      // chi2 counter
  fEventInfo_CacheTrackMatchEff->Zero();  // matchEff counter
  //
  if (fHistPhiTPCcounterA)    fHistPhiTPCcounterA->Reset();
  if (fHistPhiTPCcounterC)    fHistPhiTPCcounterC->Reset();
  if (fHistPhiTPCcounterAITS) fHistPhiTPCcounterAITS->Reset();
  if (fHistPhiTPCcounterCITS) fHistPhiTPCcounterCITS->Reset();
  if (fHistPhiITScounterA)    fHistPhiITScounterA->Reset();
  if (fHistPhiITScounterC)    fHistPhiITScounterC->Reset();
  //
  //
  Int_t nNumberOfTracks = event->GetNumberOfTracks();
  Float_t tpcDCAarrPhiA[36][nNumberOfTracks];
  Float_t tpcDCAarrPhiC[36][nNumberOfTracks];
  for (Int_t i=0;i<36;i++){
    for (Int_t j=0;j<nNumberOfTracks;j++){
      tpcDCAarrPhiA[i][j]=kDCAtpcNULL;
      tpcDCAarrPhiC[i][j]=kDCAtpcNULL;
    }
  }
  //
  // --------------------------------------------------------------
  //      Track LOOP
  // --------------------------------------------------------------
  //
  AliTPCdEdxInfo tpcdEdxInfo;
  for (Int_t itrack=0;itrack<nNumberOfTracks;++itrack)
  {

    //
    Double_t eta=-100., phiTPC=0.,sectorNumber=0., tpcdEdx=0., ptotTPC=0.;
    //
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (track == NULL) continue;
    //
    Double_t tgl        = track->Pz()/track->Pt();
    Double_t phiGlobal  = track->Phi()-TMath::Pi(); // ????
    // Int_t sign          = track->GetSign();
    Double_t phi        = track->GetParameterAtRadius(85,5,7);
    Double_t sectorNumbertmp = (9*phi/TMath::Pi()+18*(phi<0));
    eta = track->Eta();
    if (TMath::Abs(eta)>fEtaMax) continue;
    Bool_t isOnITS = track->IsOn(AliESDtrack::kITSrefit);
    Bool_t isOnTPC = track->IsOn(AliESDtrack::kTPCrefit);
    // Bool_t isOnTRD = track->IsOn(AliESDtrack::kTRDrefit);
    //
    // --------------------------------------------------------------
    //      TPC track information
    // --------------------------------------------------------------
    //
    if (track->GetInnerParam()) {
      tpcdEdx = track->GetTPCsignal();
      ptotTPC = track->GetInnerParam()->GetP();
      phiTPC  = track->GetInnerParam()->GetParameterAtRadius(85,5,7);
      sectorNumber = (9*phiTPC/TMath::Pi()+18*(phiTPC<0));
    }
    //
    // --------------------------------------------------------------
    //      Count only ITS tracks
    // --------------------------------------------------------------
    //
    if ( isOnITS && !isOnTPC ) {
      if (TMath::Abs(phi)>1e-10){
        if (tgl>0) fHistPhiITScounterA->Fill(sectorNumbertmp);
        if (tgl<0) fHistPhiITScounterC->Fill(sectorNumbertmp);
      }
    }
    //
    if (!track->GetInnerParam()) continue;  // ????
    if (track->IsOn(AliVTrack::kTPCout)==kFALSE)  continue;  // ????
    (*fEventInfo_CacheTrackCounters)[4]++;      // all TPC track with out flag
    // TPC track counters with DCAZ
    for (Int_t izCut=1; izCut<4; izCut++){
      Float_t impactParam[2];
      track->GetImpactParameters(impactParam[0],impactParam[1]);
      if (TMath::Abs(impactParam[0])>kDCACut) continue;
      if (TMath::Abs(track->GetInnerParam()->GetParameter()[1])<10.*(izCut+1.)) (*fEventInfo_CacheTrackTPCCountersZ)[izCut]++;
      if (TMath::Abs(impactParam[1])<10.*(izCut+1.)) (*fEventInfo_CacheTrackTPCCountersZ)[izCut+4]++;
    }
    //
    //
    Float_t dcaRPhi, dcaZ;
    track->GetImpactParameters(dcaRPhi, dcaZ);
    Int_t nclTPC    = track->GetTPCncls(); if (nclTPC<1) nclTPC=-1;
    Int_t nclITS    = track->GetITSNcls(); if (nclITS<1) nclITS=-1;
    Int_t nclTRD    = track->GetTRDncls(); if (nclTRD<1) nclTRD=-1;
    // Int_t nclTOF    = track->IsOn(AliVTrack::kTOFout);
    Double_t chi2TPC = TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/nclTPC));
    Double_t chi2ITS = TMath::Sqrt(TMath::Abs(track->GetITSchi2()));
    Double_t chi2TRD = TMath::Sqrt(TMath::Abs(track->GetTRDchi2()));
    // Double_t qP      = track->Charge()/track->P();
    // Double_t ptot0   = track->GetP();
    //
    // --------------------------------------------------------------
    //      Some track selections
    // --------------------------------------------------------------
    //
    if (nclTPC<kNclTPCcut) continue;
    if (TMath::Abs(tgl)>kTglCut) continue;
    if (track->Pt()<kPtCut) continue;
    if (TMath::Abs(dcaRPhi)>kDCACut || TMath::Abs(dcaZ)>kDCACut) continue;
    // if ( !( isOnITS||isOnTRD ) ) continue;   // ????
    //
    // --------------------------------------------------------------
    //      Fill TPC dca information for a given phi bin for each track
    // --------------------------------------------------------------
    //
    Int_t phiBin = fHistPhi->FindBin(phi)-1;
    Float_t pTPC[2],covTPC[3];          // p[0]=fdTPC; p[1]=fzTPC; cov[0]=fCddTPC; cov[1]=fCdzTPC; cov[2]=fCzzTPC;
    track->GetImpactParametersTPC(pTPC,covTPC);
    if (tgl>0) tpcDCAarrPhiA[phiBin][itrack]=pTPC[0];
    if (tgl<0) tpcDCAarrPhiC[phiBin][itrack]=pTPC[0];
    //
    // --------------------------------------------------------------
    //      TPC Phi counter
    // --------------------------------------------------------------
    //
    (*fEventInfo_CacheTrackCounters)[5]++;
    if (TMath::Abs(phiTPC)>1e-10){
      if (tgl>0) fHistPhiTPCcounterA->Fill(sectorNumber);
      if (tgl<0) fHistPhiTPCcounterC->Fill(sectorNumber);
      if(isOnITS){
        if (tgl>0) fHistPhiTPCcounterAITS->Fill(sectorNumber);
        if (tgl<0) fHistPhiTPCcounterCITS->Fill(sectorNumber);
      }
    }
    //
    // --------------------------------------------------------------
    //      track counter after pile up  ????
    // --------------------------------------------------------------
    //
    if (TMath::Min(chi2TPC,100.)<0) continue;
    (*fEventInfo_CacheTrackCounters)[1]++;
    Bool_t itsOK=track->IsOn(AliVTrack::kITSout) && nclITS>2  && chi2ITS>0;
    Bool_t trdOK=track->IsOn(AliVTrack::kTRDout) && nclTRD>35 && chi2TRD>0;
    Bool_t tofOK=track->IsOn(AliVTrack::kTOFout);
    //
    // --------------------------------------------------------------
    //      number of clusters cut
    // --------------------------------------------------------------
    //
    (*fEventInfo_CacheTrackNcl)[4]+=track->GetTPCncls(0, 63);
    (*fEventInfo_CacheTrackNcl)[5]+=track->GetTPCncls(64, 127);
    (*fEventInfo_CacheTrackNcl)[6]+=track->GetTPCncls(128, 159);
    (*fEventInfo_CacheTrackNcl)[1] += nclTPC;
    (*fEventInfo_CacheTrackChi2)[1]+= (chi2TPC>0) ? TMath::Sqrt(chi2TPC) : 2;   // sometimes negative chi2?

    if (itsOK && track->GetTPCdEdxInfo(tpcdEdxInfo)){

      Bool_t isOK=(tpcdEdxInfo.GetNumberOfCrossedRows(0)>kMindEdxClustersRegion);
      isOK&=(tpcdEdxInfo.GetNumberOfCrossedRows(1)>kMindEdxClustersRegion);
      isOK&=(tpcdEdxInfo.GetNumberOfCrossedRows(2)>kMindEdxClustersRegion);
      isOK&=((tpcdEdxInfo.GetSignalMax(0)>0) && (tpcdEdxInfo.GetSignalMax(1)>0) && (tpcdEdxInfo.GetSignalMax(2)>0));
      isOK&=((tpcdEdxInfo.GetSignalTot(0)>0) && (tpcdEdxInfo.GetSignalTot(1)>0) && (tpcdEdxInfo.GetSignalTot(2)>0));
      isOK&=(itsOK||trdOK);      // stronger pile-up cut requiring ITS or TRD

      if (isOK) {
        (*fEventInfo_CacheTrackCounters)[6]+=1;         // Counter with accepted TPC dEdx info
        (*fEventInfo_CacheTrackdEdxRatio)[0]+=TMath::Log(tpcdEdxInfo.GetSignalMax(3));
        (*fEventInfo_CacheTrackdEdxRatio)[1]+=TMath::Log(tpcdEdxInfo.GetSignalTot(3));
        (*fEventInfo_CacheTrackdEdxRatio)[2]+=TMath::Log(tpcdEdxInfo.GetSignalMax(0)/tpcdEdxInfo.GetSignalTot(0));
        (*fEventInfo_CacheTrackdEdxRatio)[3]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalTot(1));
        (*fEventInfo_CacheTrackdEdxRatio)[4]+=TMath::Log(tpcdEdxInfo.GetSignalMax(2)/tpcdEdxInfo.GetSignalTot(2));
        (*fEventInfo_CacheTrackdEdxRatio)[5]+=TMath::Log(tpcdEdxInfo.GetSignalMax(3)/tpcdEdxInfo.GetSignalTot(3));
        (*fEventInfo_CacheTrackdEdxRatio)[6]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalMax(0));
        (*fEventInfo_CacheTrackdEdxRatio)[7]+=TMath::Log(tpcdEdxInfo.GetSignalMax(1)/tpcdEdxInfo.GetSignalMax(2));
        (*fEventInfo_CacheTrackdEdxRatio)[8]+=TMath::Log(tpcdEdxInfo.GetSignalTot(1)/tpcdEdxInfo.GetSignalTot(0));
        (*fEventInfo_CacheTrackdEdxRatio)[9]+=TMath::Log(tpcdEdxInfo.GetSignalTot(1)/tpcdEdxInfo.GetSignalTot(2));
        //
        // --------------------------------------------------------------
        //      dEdx counter wrt splines and Bethe bloch
        // --------------------------------------------------------------
        //
        Double_t closestPar[3];    // closestPar[0] --> closest spline, Int_t(closestPar[1]) --> particle index,  closestPar[2] --> corresponding particle mass
        GetExpecteds(track,closestPar);
        if (closestPar[0]>0) (*fEventInfo_CacheTrackdEdxRatio)[10]+=TMath::Log(tpcdEdx/closestPar[0]);   // TODO
        (*fEventInfo_CacheTrackdEdxRatio)[11]+=TMath::Log((tpcdEdxInfo.GetSignalMax(0)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[12]+=TMath::Log((tpcdEdxInfo.GetSignalMax(1)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[13]+=TMath::Log((tpcdEdxInfo.GetSignalMax(2)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[14]+=TMath::Log((tpcdEdxInfo.GetSignalMax(3)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[15]+=TMath::Log((tpcdEdxInfo.GetSignalTot(0)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[16]+=TMath::Log((tpcdEdxInfo.GetSignalTot(1)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[17]+=TMath::Log((tpcdEdxInfo.GetSignalTot(2)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[18]+=TMath::Log((tpcdEdxInfo.GetSignalTot(3)/50.)/AliExternalTrackParam::BetheBlochAleph(ptotTPC/closestPar[2]));
        (*fEventInfo_CacheTrackdEdxRatio)[19]+=TMath::Log((tpcdEdxInfo.GetSignalMax(0)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[20]+=TMath::Log((tpcdEdxInfo.GetSignalMax(1)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[21]+=TMath::Log((tpcdEdxInfo.GetSignalMax(2)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[22]+=TMath::Log((tpcdEdxInfo.GetSignalMax(3)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[23]+=TMath::Log((tpcdEdxInfo.GetSignalTot(0)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[24]+=TMath::Log((tpcdEdxInfo.GetSignalTot(1)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[25]+=TMath::Log((tpcdEdxInfo.GetSignalTot(2)/50.));
        (*fEventInfo_CacheTrackdEdxRatio)[26]+=TMath::Log((tpcdEdxInfo.GetSignalTot(3)/50.));

      }
    }

    if (itsOK) {  // ITS
      (*fEventInfo_CacheTrackCounters)[0]++;
      (*fEventInfo_CacheTrackNcl)[0] += nclITS;
      (*fEventInfo_CacheTrackChi2)[0] += TMath::Min(TMath::Sqrt(chi2ITS),10.); // cutoff chi2 10
      (*fEventInfo_CacheTrackMatchEff)[2]+=trdOK;
      (*fEventInfo_CacheTrackMatchEff)[3]+=tofOK;
      (*fEventInfo_CacheTrackChi2)[4]+= (fTrackChi2TPC>0) ? TMath::Sqrt(fTrackChi2TPC) : 2; // TPC chi2 in case prolongation to ITS
      // long tracks properties
      if (nclITS>4){
        (*fEventInfo_CacheTrackCounters)[7]++;
        (*fEventInfo_CacheTrackNcl)[7] += nclITS;
        (*fEventInfo_CacheTrackChi2)[7]+=TMath::Min(TMath::Sqrt(chi2ITS),10.);
      }
    }
    if (trdOK) {// TRD    ///TODO - why chi2TRD could be smaller than 0?
      (*fEventInfo_CacheTrackCounters)[2]++;
      (*fEventInfo_CacheTrackNcl)[2] += nclTRD;
      (*fEventInfo_CacheTrackChi2)[2] += TMath::Sqrt(chi2TRD);
      (*fEventInfo_CacheTrackMatchEff)[0]+=itsOK;
      (*fEventInfo_CacheTrackChi2)[5]+= (fTrackChi2TPC>0) ? TMath::Sqrt(fTrackChi2TPC) : 2; // TPC chi2 in case prolongation to TRD
      if (nclTRD>80){
        (*fEventInfo_CacheTrackCounters)[8]++;
        (*fEventInfo_CacheTrackNcl)[8] += nclTRD;
        (*fEventInfo_CacheTrackChi2)[8]+=TMath::Min(TMath::Sqrt(chi2TRD),10.);
      }
    }
    if (tofOK) {  // TOF
      (*fEventInfo_CacheTrackCounters)[3]++;
      (*fEventInfo_CacheTrackdEdxRatio)[3] += 1;   // dummy for the moment
      (*fEventInfo_CacheTrackChi2)[3]+= 1;   //
    }
  } // end of track LOOP
  //
  // ======================================================================
  //  calculate event averages
  // ======================================================================
  //
  for (Int_t i=0; i<9; i++) if ((*fEventInfo_CacheTrackCounters)[i]>0) (*fEventInfo_CacheTrackdEdxRatio)[i]/=(*fEventInfo_CacheTrackCounters)[i];
  for (Int_t i=0; i<4; i++) if ((*fEventInfo_CacheTrackCounters)[i]>0) (*fEventInfo_CacheTrackChi2)[i]/=(*fEventInfo_CacheTrackCounters)[i];

  for (Int_t i=4; i<7; i++)  if ((*fEventInfo_CacheTrackCounters)[1]>0) (*fEventInfo_CacheTrackdEdxRatio)[i]/=(*fEventInfo_CacheTrackCounters)[1];
  //
  if ((*fEventInfo_CacheTrackCounters)[6]>0){
    for (Int_t i=0; i<27; i++)   (*fEventInfo_CacheTrackdEdxRatio)[i]/=(*fEventInfo_CacheTrackCounters)[6];
  }
  //
  // conditional matching efficiency and chi2
  if ((*fEventInfo_CacheTrackCounters)[0]>0){
    (*fEventInfo_CacheTrackMatchEff)[2]/=(*fEventInfo_CacheTrackCounters)[0];  // TRD if ITS
    (*fEventInfo_CacheTrackMatchEff)[3]/=(*fEventInfo_CacheTrackCounters)[0];  // TOF if ITS
    (*fEventInfo_CacheTrackChi2)[4]/=(*fEventInfo_CacheTrackCounters)[0];
  }
  if ((*fEventInfo_CacheTrackCounters)[2]>0) {
    (*fEventInfo_CacheTrackMatchEff)[0]/=(*fEventInfo_CacheTrackCounters)[2];
    (*fEventInfo_CacheTrackChi2)[5]/=(*fEventInfo_CacheTrackCounters)[2];
  } //ITS if TRD
  (*fEventInfo_CacheTrackCounters)[9]=event->GetNumberOfTracks();  // original number of ESDtracks
  //
  //
  for (Int_t iphi=0; iphi<36; iphi++){

    // count nonzero entries
    Int_t countNonZerosA=0;
    Int_t countNonZerosC=0;
    for (Int_t itrack=0;itrack<nNumberOfTracks;itrack++) {
      if (tpcDCAarrPhiA[iphi][itrack]>kDCAtpcNULL) countNonZerosA++;
      if (tpcDCAarrPhiC[iphi][itrack]>kDCAtpcNULL) countNonZerosC++;
    }
    //
    // create arrays from nonzero entries for A side
    Float_t tmpA[countNonZerosA];
    Int_t j=0;
    for (Int_t itrack=0;itrack<nNumberOfTracks;itrack++) {
      if (tpcDCAarrPhiA[iphi][itrack]>kDCAtpcNULL) { tmpA[j]=tpcDCAarrPhiA[iphi][itrack]; j++; }
    }
    //
    // create arrays from nonzero entries for C side
    Float_t tmpC[countNonZerosC];
    Int_t k=0;
    for (Int_t itrack=0;itrack<nNumberOfTracks;itrack++) {
      if (tpcDCAarrPhiC[iphi][itrack]>kDCAtpcNULL) { tmpC[k]=tpcDCAarrPhiC[iphi][itrack]; k++; }
    }
    (*fEventInfo_PhiTPCdcarA)[iphi]=TMath::Median(countNonZerosA,tmpA);
    (*fEventInfo_PhiTPCdcarC)[iphi]=TMath::Median(countNonZerosC,tmpC);

  }

}
//________________________________________________________________________
Int_t AliAnalysisTaskTIdentityPID::CacheTPCEventInformation()
{

  AliVEvent *event=InputEvent();
  const Int_t kNCRCut=80;
  const Double_t kDCACut=5;
  const Float_t knTrackletCut=1.5;
  // FILL DCA histograms
  fEventInfo_HisTPCVertexA->Reset();
  fEventInfo_HisTPCVertexC->Reset();
  fEventInfo_HisTPCVertexACut->Reset();
  fEventInfo_HisTPCVertexCCut->Reset();
  fEventInfo_HisTPCVertex->Reset();

  Int_t nTracks=event->GetNumberOfTracks();
  Int_t selected=0;
  for (Int_t iTrack=0; iTrack<nTracks; iTrack++){
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack * pTrack = fESD->GetTrack(iTrack);
    if (pTrack==nullptr) continue;
    if (pTrack->IsOn(AliESDtrack::kTPCin)==0) continue;
    if (pTrack->GetTPCClusterInfo(3,1)<kNCRCut) continue;
    Float_t dcaxy,dcaz;
    pTrack->GetImpactParameters(dcaxy,dcaz);
    if (TMath::Abs(dcaxy)>kDCACut) continue;
    pTrack->SetESDEvent(fESD);
    selected++;
    if ((pTrack->GetNumberOfTRDClusters()/20.+pTrack->GetNumberOfITSClusters())>knTrackletCut){
      fEventInfo_HisTPCVertex->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()>0) fEventInfo_HisTPCVertexACut->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()<0) fEventInfo_HisTPCVertexCCut->Fill(pTrack->GetTPCInnerParam()->GetZ());
    }else{
      if (pTrack->GetTgl()>0) fEventInfo_HisTPCVertexA->Fill(pTrack->GetTPCInnerParam()->GetZ());
      if (pTrack->GetTgl()<0) fEventInfo_HisTPCVertexC->Fill(pTrack->GetTPCInnerParam()->GetZ());
    }
  }
  if (fUseCouts) printf(" Info::marsland: ===== In the CacheTPCEventInformation:  %d\n",selected);
  return selected;
}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::DumpDownScaledTree()
{

  if (fUseCouts) std::cout << " Info::marsland: ===== In the DumpDownScaledTree ===== " << std::endl;

  AliVEvent *event=InputEvent();
  const Int_t kNclTPCcut=60;
  const Float_t kTglCut=1.5;
  const Float_t kPtCut=0.100;
  //
  // --------------------------------------------------------------
  //  Counter only for the TPC track multiplicity
  // --------------------------------------------------------------
  //
  // Go into track loop
  AliTPCdEdxInfo tpcdEdxInfo;
  TRandom r;
  Int_t nNumberOfTracks = event->GetNumberOfTracks();
  for (Int_t itrack=0;itrack<nNumberOfTracks;++itrack)
  {
    fTrackCutBits=0;  // reset the bits for the next track
    AliESDtrack *track = fESD->GetTrack(itrack);
    if (!track) continue;
    if (!fESDtrackCutsLoose->AcceptTrack(track))  continue;    // Loose cuts
    if (!(track->GetTPCsignalN()>0)) continue;
    Int_t tpcCrossedRows=0, tpcSignalN=0;
    Double_t eta=-100.;
    Double_t tgl  = track->Pz()/track->Pt();
    Double_t phi  = track->Phi()-TMath::Pi(); // ????
    Double_t pt   = track->Pt();
    Int_t sign    = track->GetSign();
    Double_t phi85 = track->GetParameterAtRadius(85,5,7);
    ULong64_t flag = track->GetStatus();
    eta = track->Eta();
    Bool_t isOnITS = track->IsOn(AliESDtrack::kITSrefit);
    Bool_t isOnTRD = track->IsOn(AliESDtrack::kTRDrefit);
    Bool_t isOnTPC = track->IsOn(AliESDtrack::kTPCrefit);
    Int_t nclTPC   = track->GetTPCncls(); if (nclTPC<1) nclTPC=-1;
    Int_t nclITS   = track->GetITSNcls(); if (nclITS<1) nclITS=-1;
    Int_t nclTRD   = track->GetTRDncls(); if (nclTRD<1) nclTRD=-1;
    Double_t chi2TPC = TMath::Sqrt(TMath::Abs(track->GetTPCchi2()/nclTPC));
    Double_t chi2ITS = TMath::Sqrt(TMath::Abs(track->GetITSchi2()));
    Double_t chi2TRD = TMath::Sqrt(TMath::Abs(track->GetTRDchi2()));
    Double_t itsdEdx = track->GetITSsignal();
    Double_t trddEdx = track->GetTRDsignal();
    Double_t tpcdEdx = track->GetTPCsignal();
    Double_t ptot0   = track->GetP();
    Double_t qP      = track->Charge()/track->P();
    Float_t dcaRPhi, dcaZ;
    track->GetImpactParameters(dcaRPhi, dcaZ);
    Bool_t fBit96_base   = fESDtrackCuts_Bit96->AcceptTrack(track);
    Bool_t fBit128       = fESDtrackCuts_Bit128->AcceptTrack(track);
    Bool_t fBit768       = fESDtrackCuts_Bit768->AcceptTrack(track);
    Bool_t ifDCAcutIfNoITSPixel = ApplyDCAcutIfNoITSPixel(track);
    //
    // --------------------------------------------------------------
    //      Preparation for downscaled tree
    // --------------------------------------------------------------
    //
    Double_t sharedTPCClusters=0.;
    Double_t phiTPC=-100.;
    Double_t ptotTPC=0.;
    Double_t lengthInActiveZone=0.;
    Float_t pTPC[2],covTPC[3];          // p[0]=fdTPC; p[1]=fzTPC; cov[0]=fCddTPC; cov[1]=fCdzTPC; cov[2]=fCzzTPC;
    if (track->GetInnerParam()){
      track->GetTPCdEdxInfo(tpcdEdxInfo);
      phiTPC  = track->GetInnerParam()->GetParameterAtRadius(85,5,7);
      ptotTPC = track->GetInnerParam()->GetP();
      lengthInActiveZone = track->GetLengthInActiveZone(1,3,230, track->GetBz(),0,0);
      track->GetImpactParametersTPC(pTPC,covTPC);
      tpcCrossedRows = track->GetTPCCrossedRows();
      tpcSignalN = track->GetTPCsignalN();
      Double_t closestPar[3];
      GetExpecteds(track,closestPar);
      SetCutBitsAndSomeTrackVariables(track,0);
      if(nclTPC) sharedTPCClusters = static_cast<Double_t>(track->GetTPCnclsS())/static_cast<Double_t>(nclTPC);
    }
    UChar_t itsclmap = track->GetITSClusterMap();
    Float_t pv[2],cov[3]; //  fTrackDCAxy = pv[0] and fTrackDCAz  = pv[1]
    track->GetImpactParameters(pv,cov); // p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
    Bool_t dca11h     = TMath::Abs(pv[0])<0.0105+0.0350/TMath::Power(pt,1.1);    // 10h tuned loose cut
    Bool_t dca10h     = TMath::Abs(pv[0])<0.0182+0.0350/TMath::Power(pt,1.01);    // 10h tuned loose cut
    Bool_t dcaBaseCut = TMath::Abs(pv[0])<0.0208+0.04/TMath::Power(pt,1.01);  // 10h tuned loose cut
    // Get TOF beta
    Double_t tofSignalTunedOnData = track->GetTOFsignalTunedOnData();
    Double_t length = track->GetIntegratedLength();
    Double_t tofSignal = track->GetTOFsignal();
    Double_t beta = -.05;
    if((length > 0) && (tofSignal > 0)) beta = length / 2.99792458e-2 / tofSignal;


    // High dEdx its conditon
    TDatabasePDG *pdg = TDatabasePDG::Instance();
    Double_t mproton = pdg->GetParticle(2212)->Mass();  // Double_t mproton = 9.3827199e-01;
    Bool_t itsHighDeDx = kFALSE;
    if (ptot0>0.1 && ptot0<100 && itsdEdx>0){
      itsHighDeDx = ( (TMath::Log(itsdEdx/AliExternalTrackParam::BetheBlochSolid(ptot0/mproton))>11.5) && (nclITS>4) );
    }
    //
    // --------------------------------------------------------------
    //   Fill downscaled tree
    // --------------------------------------------------------------
    //
    if (nclTPC<kNclTPCcut) continue;
    if (TMath::Abs(tgl)>kTglCut) continue;
    if (track->Pt()<kPtCut) continue;
    if ( (fRandom.Rndm()*(qP*qP) < 0.005) || itsHighDeDx )
    {
      if(!fTreeSRedirector) return;
      (*fTreeSRedirector)<<"tracksdscaled"<<
      "p=" << ptot0 << // vertex momentum momentum
      "pt=" << pt << // pT
      "ptot=" << ptotTPC << // TPC momentum
      "eta=" << eta << // eta
      "dEdx=" << tpcdEdx << // dEdx of the track
      "beta=" << beta << // TOF beta  --> beta = GetIntegratedLength / 2.99792458e-2 / tofSignal;
      "phi=" << phi << // ph
      "phiTPC=" << phiTPC <<
      "phi85=" << phi85 <<
      "sign=" << sign << // charge
      "cent=" << fCentrality << // centrality
      "qP=" << qP << // charge/momentu,
      "tgl=" << tgl << // tangent
      "centEst.=" << fEventInfo_CentralityEstimates << // track counter
      "dEdxInfo.=" << &tpcdEdxInfo << // TPC dEdx info
      //
      "gid=" << fEventGID << // global event ID
      "event=" << fEventCountInFile <<
      "cutBit=" << fTrackCutBits << // Systematic Cuts
      "run=" << fRunNo << // run number
      "intrate=" << fIntRate << // interaction rate
      "timestamp=" << fTimeStamp << // timestamp
      "tpcMult=" << fTPCMult << // TPC multiplicity
      "primMult=" << fNContributors << // #prim tracks
      "eventmult=" << fEventMult << // event multiplicity
      "vZ=" << fVz << // vertex Z
      "tpcvz=" << fTPCvZ << // TPC event vertex
      "spdvz=" << fSPDvZ <<
      //
      "defCut=" << fDefaultCuts << // default cuts tuned by hand
      "bit96=" << fBit96_base << // tight cuts of 2011 tuned data
      "bit128=" << fBit128 << // TPC only tracks cuts
      "bit768=" << fBit768 << // Hybrid track cuts
      "pixCut=" << ifDCAcutIfNoITSPixel << // cut: apply a DCAcut If No ITS Pixel
      "dcabase=" << dcaBaseCut << // dcaxy base cut
      "dca10h=" << dca10h << // dcaxy cut tuned to 10h
      "dca11h=" << dca11h << // dcaxy cut tuned to 11h
      //
      "expel=" << fDEdxEl <<
      "exppi=" << fDEdxPi <<
      "expka=" << fDEdxKa <<
      "exppr=" << fDEdxPr <<
      "expde=" << fDEdxDe <<
      //
      "eltpcpid=" << fNSigmasElTPC << // nsigma TPC for electrons
      "pitpcpid=" << fNSigmasPiTPC << // nsigma TPC for pions
      "katpcpid=" << fNSigmasKaTPC << // nsigma TPC for kaons
      "prtpcpid=" << fNSigmasPrTPC << // nsigma TPC for protons
      "detpcpid=" << fNSigmasPrTOF << // nsigma TPC for deuterons
      //
      "eltofpid=" << fNSigmasElTOF << // nsigma TOF for electrons
      "pitofpid=" << fNSigmasPiTOF << // nsigma TOF for pions
      "katofpid=" << fNSigmasKaTOF << // nsigma TOF for kaons
      "prtofpid=" << fNSigmasPrTOF << // nsigma TOF for protons
      "detofpid=" << fNSigmasDeTOF << // nsigma TOF for deuterons
      //
      "itsdEdx=" << itsdEdx <<
      "trddEdx=" << trddEdx <<
      "isOnITS=" << isOnITS <<
      "isOnTRD=" << isOnTRD <<
      "isOnTPC=" << isOnTPC <<
      "ncltpc=" << nclTPC << // #ITS clusters
      "nclits=" << nclITS << // #ITS clusters
      "ncltrd=" << nclTRD << // #TRD clusters
      "lengthInActiveZone=" << lengthInActiveZone << // fTrackLengthInActiveZone in TPC
      "sharedTPCClusters=" << sharedTPCClusters <<
      "tpcSignalN=" << tpcSignalN << // number of cl used in dEdx
      "cRows=" << tpcCrossedRows << // crossed rows
      "missCl=" << fMissingCl << // fraction of missing clusters
      "flag=" << flag <<
      "fd=" << pv[0] <<
      "fz=" << pv[1] <<
      "fCdd=" << cov[0] <<
      "fCdz=" << cov[1] <<
      "fCzz=" << cov[2] <<
      "fdTPC=" << pTPC[0] <<
      "fzTPC=" << pTPC[1] <<
      "fCddTPC=" << covTPC[0] <<
      "fCdzTPC=" << covTPC[1] <<
      "fCzzTPC=" << covTPC[2] <<
      "chi2tpc=" << chi2TPC << // TPC chi2
      "chi2its=" << chi2ITS << // ITS chi2
      "chi2trd=" << chi2TRD << // TRD chi2
      "itsclmap=" << itsclmap <<  //  vertex Z
      "itspixel01=" << fIsITSpixel01 <<
      "primRes=" << fPrimRestriction <<
      "tofSignal=" << tofSignal << // tof signal --> GetTOFsignal
      "tofSignalTOD=" << tofSignalTunedOnData << // tof signal tuneondata --> GetTOFsignalTunedOnData
      "flag=" << flag <<
      "\n";
    }

  } // end of track LOOP

}
//________________________________________________________________________
Bool_t AliAnalysisTaskTIdentityPID::CheckPsiPair(const AliESDv0* v0)
{
  // Angle between daughter momentum plane and plane
  // taken from AliESDv0KineCuts

  if(!fESD) return kFALSE;

  Float_t magField = fESD->GetMagneticField();

  // check if indices have been correctly applied
  Int_t pIndexTemp = -1;
  Int_t nIndexTemp = -1;

  pIndexTemp = v0->GetPindex();
  nIndexTemp = v0->GetNindex();

  AliESDtrack* d[2];
  d[0] = dynamic_cast<AliESDtrack*>(fESD->GetTrack(pIndexTemp));
  d[1] = dynamic_cast<AliESDtrack*>(fESD->GetTrack(nIndexTemp));

  Int_t sign[2];
  sign[0] = (int)d[0]->GetSign();
  sign[1] = (int)d[1]->GetSign();

  Int_t pIndex = 0, nIndex = 0;
  if(-1 == sign[0] && 1 == sign[1]){
    pIndex = v0->GetPindex();
    nIndex = v0->GetNindex();
  }
  else{
    pIndex = v0->GetNindex();
    nIndex = v0->GetPindex();
  }

  AliESDtrack* daughter[2];

  daughter[0] = dynamic_cast<AliESDtrack *>(fESD->GetTrack(pIndex));
  daughter[1] = dynamic_cast<AliESDtrack *>(fESD->GetTrack(nIndex));

  Double_t x, y, z;
  v0->GetXYZ(x,y,z);//Reconstructed coordinates of V0; to be replaced by Markus Rammler's method in case of conversions!

  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};


  v0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  v0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;


  Double_t deltat = 1.;
  deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) -  TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis

  Double_t radiussum = TMath::Sqrt(x*x + y*y) + 50;//radius to which tracks shall be propagated

  Double_t momPosProp[3];
  Double_t momNegProp[3];

  AliExternalTrackParam pt(*daughter[0]), nt(*daughter[1]);

  Double_t psiPair = 4.;

  if(nt.PropagateTo(radiussum,magField) == 0) psiPair = -5.;
  if(pt.PropagateTo(radiussum,magField) == 0) psiPair = -5.;
  pt.GetPxPyPz(momPosProp);//Get momentum vectors of tracks after propagation
  nt.GetPxPyPz(momNegProp);

  Double_t pEle = TMath::Sqrt(momNegProp[0]*momNegProp[0]+momNegProp[1]*momNegProp[1]+momNegProp[2]*momNegProp[2]);//absolute momentum value of negative daughter
  Double_t pPos = TMath::Sqrt(momPosProp[0]*momPosProp[0]+momPosProp[1]*momPosProp[1]+momPosProp[2]*momPosProp[2]);//absolute momentum value of positive daughter
  Double_t scalarproduct = momPosProp[0]*momNegProp[0]+momPosProp[1]*momNegProp[1]+momPosProp[2]*momNegProp[2];//scalar product of propagated positive and negative daughters' momenta

  Double_t chipair = TMath::ACos(scalarproduct/(pEle*pPos));//Angle between propagated daughter tracks

  return abs(deltat / chipair) <= 1;
}
//______________________________________________________________________________
Double_t AliAnalysisTaskTIdentityPID::GetTrackEfficiency(const Int_t& part, const Double_t& ptot, const Double_t& eta, const Int_t& setting, const Int_t& sign, const Int_t& pidCutType) {
  //
  // Get the efficiency for a given particle species, momentum, eta, centrality, syst setting, sign, and TOF setting
  Double_t ret = -1.;
  Int_t signIndex = (sign == 1) ? 0 : 1;

  // get centrality index
  Int_t centIndex = -1;
  for (Int_t iCent = 0; iCent < fEffMatrixCentBins.size(); iCent++) {
    if (fCentrality >= fEffMatrixCentBins[iCent] && fCentrality < fEffMatrixCentBins[iCent+1]) {
      centIndex = iCent;
      break;
    }
  }
  if (fCollisionType == 1) centIndex = 0;
  //
  // make projections
  if (!fEffMatrixProjections[centIndex][setting][signIndex][pidCutType]) {
    THnF* effMatrixGen = (sign == 1) ? fEffMatrixGenPos : fEffMatrixGenNeg;
    THnF* effMatrixRec = (sign == 1) ? fEffMatrixRecPos : fEffMatrixRecNeg;

    for (THnF* effMatrix: {effMatrixGen, effMatrixRec}) {
      effMatrix->GetAxis(1)->SetRangeUser(0, 1); // origin
      effMatrix->GetAxis(3)->SetRangeUser(part, part + 1);
      effMatrix->GetAxis(0)->SetRangeUser(pidCutType, pidCutType + 1);
      if (effMatrix == effMatrixRec) {
        effMatrix->GetAxis(2)->SetRangeUser(setting, setting+1);
      } else {
        effMatrix->GetAxis(2)->SetRangeUser(0, 1);
      }
      effMatrix->GetAxis(4)->SetRangeUser(fEffMatrixCentBins[centIndex], fEffMatrixCentBins[centIndex + 1]);
    }

    TH2F* tmpEffGen = (TH2F*)effMatrixGen->Projection(5, 6)->Clone();
    TH2F* tmpEffRec = (TH2F*)effMatrixRec->Projection(5, 6)->Clone();
    tmpEffRec->Divide(tmpEffGen);
    fEffMatrixProjections[centIndex][setting][signIndex][pidCutType] = tmpEffRec;
  }

  ret = fEffMatrixProjections[centIndex][setting][signIndex][pidCutType]->GetBinContent(fEffMatrixProjections[centIndex][setting][signIndex][pidCutType]->FindBin(eta, ptot));
  return ret;
};
//________________________________________________________________________
Float_t AliAnalysisTaskTIdentityPID::RelativePhi(Float_t mphi, Float_t vphi) {
  if (mphi < -TMath::Pi()){
    mphi += TMath::TwoPi();
  }
  else if (mphi > TMath::Pi()){
    mphi -= TMath::TwoPi();
  }
  if (vphi < -TMath::Pi()){
    vphi += TMath::TwoPi();
  }
  else if (vphi > TMath::Pi()){
    vphi -= TMath::TwoPi();
  }
  Float_t dphi = mphi - vphi;
  if (dphi < -TMath::Pi()){
    dphi += TMath::TwoPi();
  }
  else if (dphi > TMath::Pi()){
    dphi -= TMath::TwoPi();
  }
  return dphi;  // returns them in the range [-pi,pi]
}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::FindJetsFJ()
{
  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FindJetsFJ ===== " << std::endl;
  //
  // Create jetwrapper with the same settings used in FindJetsEMC
  float fTrackMinPt = 0.15;
  float fMaxRap = 0.9;
  float fGhostArea = 0.005;
  float bgJetAbsEtaCut = 0.7;
  float bgJetRadius = 0.2;
  //
  int nJetRadiusBins = 3;
  std::vector<float> fJetRadius{0.2,0.4,0.6};
  //
  // loop over settings and jets radius
  for (size_t iset=0; iset<fSystSettings.size(); iset++){
    Int_t setting = fSystSettings[iset];
    //
    // Check if the event selection is passed
    Int_t countTrack = 0;
    for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
      AliESDtrack* track = fESD->GetTrack(iTrack);
      if (!track) continue;
      if (fMCEvent){
        Int_t lab = TMath::Abs(track->GetLabel());
        Bool_t bPrim = fMCStack->IsPhysicalPrimary(lab);
        fIsMCPileup = IsFromPileup(lab);
        if (fIsMCPileup) continue;
        if (!bPrim) continue;
      } else {
        SetCutBitsAndSomeTrackVariables(track,0);
        if (!GetSystematicClassIndex(fTrackCutBits,setting)) continue;
      }
      countTrack++;
    }
    if (countTrack<1) return;
    //
    // radius loop
    for (int iJetRadius=0; iJetRadius<nJetRadiusBins; iJetRadius++){
      //
      // run jet finder only for set==-1,0 and 4 in case of MC
      if (fMCtrue && (setting > 0) ) continue; //
      Float_t jetAbsEtaCut = 0.9-fJetRadius[iJetRadius];   // fixed
      double particleEtaCut = 0.9;
      //
      //SOME CODE FROM NIMA
      AliFJWrapper *fFastJetWrapper;
      fFastJetWrapper = new AliFJWrapper("fFastJetWrapper","fFastJetWrapper");
      fFastJetWrapper->Clear();
      fFastJetWrapper->SetR(fJetRadius[iJetRadius]);
      fFastJetWrapper->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
      fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::E_scheme); // fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
      fFastJetWrapper->SetStrategy(fastjet::Strategy::Best);
      fFastJetWrapper->SetGhostArea(fGhostArea);
      fFastJetWrapper->SetAreaType(fastjet::AreaType::active_area);
      fFastJetWrapper->SetMaxRap(fMaxRap);
      fFastJetWrapper->SetMinJetPt(fTrackMinPt);
      std::vector<int> trackTTIndex;
      trackTTIndex.clear();
      std::vector<fastjet::PseudoJet> particlesEmbeddedSubtracted; //will be filled with your subtracted event
      std::vector<fastjet::PseudoJet> particlesEmbedded; //fill this with your event
      //
      // loop over esd tracks and add their four vector to wrapper --> identical to track container in EMC jet
      Int_t nTracksSelected = 0;
      for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
        AliESDtrack* track = fESD->GetTrack(iTrack);
        if (TMath::Abs(track->Eta()) > 0.9) continue;
        if (!track->GetInnerParam()) continue;
        if (!(track->GetTPCsignalN()>0)) continue;
        SetCutBitsAndSomeTrackVariables(track,0);
        if (!GetSystematicClassIndex(fTrackCutBits,setting)) continue;
        //
        if (fMCEvent){
          Int_t lab = TMath::Abs(track->GetLabel());
          Bool_t bPrim = fMCStack->IsPhysicalPrimary(lab);
          fIsMCPileup = IsFromPileup(lab);
          if (fIsMCPileup) continue;
          //
          // only primary particle condition for set=4 for set==-1 and set==0 jet finder runs over all selected particles
          if (setting>0 && !bPrim) continue;
        }
        //
        if (track->Pt() < fTrackMinPt || TMath::Abs(track->Eta()) >= particleEtaCut) continue;
        fFastJetWrapper->AddInputVector(track->Px(), track->Py(), track->Pz(), track->E(), iTrack);//TMath::Sqrt(track->P()*track->P()+0.13957*0.13957),iTrack);
        particlesEmbedded.push_back(fastjet::PseudoJet(track->Px(), track->Py(), track->Pz(), track->E()));// TMath::Sqrt(track->P()*track->P()+0.13957*0.13957) ) );
        nTracksSelected++;
      }
      //
      // background jet definitions
      fastjet::JetMedianBackgroundEstimator bgE;
      fastjet::Selector selectorBG = !fastjet::SelectorNHardest(2) * fastjet::SelectorAbsEtaMax(bgJetAbsEtaCut) * fastjet::SelectorPtRange(fTrackMinPt, 1000.0); //set the max eta cut on the estimator, then get rid of 2 highest pt jets
      bgE.set_selector(selectorBG);
      fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, bgJetRadius, fastjet::E_scheme, fastjet::Best); //define the kT jet finding which will do the average background estimation
      fastjet::GhostedAreaSpec ghostSpecBG(particleEtaCut, 1, fGhostArea); //this ghost area might be too small and increase processing time too much
      fastjet::AreaDefinition areaDefBG(fastjet::active_area_explicit_ghosts, ghostSpecBG);
      fastjet::ClusterSequenceArea cluster_seq_BG(particlesEmbedded, jetDefBG, areaDefBG);
      std::vector<fastjet::PseudoJet> jetsBG = sorted_by_pt(selectorBG(cluster_seq_BG.inclusive_jets())); //find the kT jets
      if (jetsBG.size() > 0) {
        bgE.set_jets(jetsBG);  // give the kT jets to the background estimator
        fRhoFJ = bgE.rho();
      }
      //
      // start of background jet loop
      if (fFillJetsBG==1){
        for (Int_t ijet=0; ijet<Int_t(jetsBG.size()); ijet++) {
          fastjet::PseudoJet jetbg = jetsBG[ijet];
          Float_t jetpt = jetbg.pt();
          Float_t jetphi = jetbg.phi();
          Float_t jeteta = jetbg.eta();
          Float_t jetArea = jetbg.area();
          Float_t jetptsub = jetpt - fRhoFJ*jetArea;
          Int_t nJets = jetsBG.size();
          (*fTreeSRedirector)<<"jetsFJBG"<<
          "ijet=" << ijet <<
          "gid=" << fEventGID << // global event ID
          "syst=" << setting << // syst setting
          "jetRadiusBG=" << bgJetRadius << // jet Radius
          "jetEtaCutBG=" << bgJetAbsEtaCut << //abs eta cut for jet
          "nJets=" << nJets <<
          "jetpt=" << jetpt <<
          "jetphi=" << jetphi <<
          "jeteta=" << jeteta <<
          "jetptsub=" << jetptsub << //bg sub jet pt (pt - rho*Area)
          "rhoFJ=" << fRhoFJ << //event rho
          "jetArea=" << jetArea << //jet area
          "cent=" << fCentrality << // centrality
          "pileupbit=" << fPileUpBit <<
          "\n";
        } // end of background jet loop
      }
      //
      // background subtraction on the constituent level TODO
      // fastjet::contrib::ConstituentSubtractor subtractorConstituent(&bgE); //add the background estimator to the correct subtractor
      // subtractorConstituent.set_common_bge_for_rho_and_rhom(true); //CHECK : should not be the case since particles have mass
      // subtractorConstituent.set_max_standardDeltaR(0.25); // set the max event wise subtraction distance
      // particlesEmbeddedSubtracted = subtractorConstituent.subtract_event(particlesEmbedded, particleEtaCut); //perform subtraction and fill the subtracted event container
      //
      // run jet finder using wrapper
      fFastJetWrapper->Run();
      std::vector<fastjet::PseudoJet> jets = fFastJetWrapper->GetInclusiveJets();
      auto nJets = jets.size();
      //
      // start of jet loop
      std::vector<Int_t> fJetConstituentLabels;
      TVectorF fShapesVar_Particles_E(nTracksSelected);
      TVectorF fShapesVar_Particles_pT(nTracksSelected);
      TVectorF fShapesVar_Particles_Phi(nTracksSelected);
      TVectorF fShapesVar_Particles_Theta(nTracksSelected);
      TVectorF fShapesVar_Particles_InJet(nTracksSelected);
      TVectorF fShapesVar_Particles_DeltaR(nTracksSelected);
      TVectorF fShapesVar_Particles_NRPhi(nTracksSelected);
      TVectorF fShapesVar_Particles_Eta(nTracksSelected);
      TVectorF fShapesVar_Particles_dEdx(nTracksSelected);
      for (Int_t ijet=0; ijet<Int_t(nJets); ijet++){
        //
        // get the jet object
        fJetConstituentLabels.clear();
        for (Int_t i=1;i<nTracksSelected;i++)
        {
          fShapesVar_Particles_E[i] = 0.;
          fShapesVar_Particles_pT[i] = 0.;
          fShapesVar_Particles_Phi[i] = 0.;
          fShapesVar_Particles_Theta[i] = 0.;
          fShapesVar_Particles_InJet[i] = 0.;
          fShapesVar_Particles_DeltaR[i] = 0.;
          fShapesVar_Particles_NRPhi[i] = 0.;
          fShapesVar_Particles_Eta[i] = 0.;
          fShapesVar_Particles_dEdx[i] = 0.;
        }
        fastjet::PseudoJet jet = jets[ijet];
        if (jet.pt() < fTrackMinPt || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= jetAbsEtaCut) continue;
        //
        // get the jet constituents
        std::vector<fastjet::PseudoJet> constituents(fFastJetWrapper->GetJetConstituents(ijet));
        Int_t nConstituents = constituents.size();
        std::vector<fastjet::PseudoJet> sorted_constituents = sorted_by_pt(constituents);
        auto leadingPt = sorted_constituents[0].perp();
        Float_t jetpx = jet.px();
        Float_t jetpy = jet.py();
        Float_t jetpz = jet.pz();
        Float_t jetE = jet.E();
        Float_t jetpt = jet.pt();
        Float_t jetphi = jet.phi();
        Float_t jeteta = jet.eta();
        Float_t jetArea = jet.area();
        Float_t jetptsub = jetpt - fRhoFJ*jetArea;
        Float_t jetMass = ((jetE * jetE) - (jetpt * jetpt) - (jetpz * jetpz) >0) ? TMath::Sqrt((jetE * jetE) - (jetpt * jetpt) - (jetpz * jetpz)) : 1; // TODO
        fJetHistptSub->Fill(jetptsub);
        if (jetpt<5) continue;
        //
        // Nsubjettiness
        // AliEmcalJetFinder::Nsubjettiness( *pJet, *pContJets,  Double_t dVtx[3], Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Int_t Option, Int_t Measure, Double_t Beta_SD, Double_t ZCut, Int_t SoftDropOn){
        // void AliAnalysisTaskJetPlanarFlow::SetTree(AliEmcalJet *jet, AliJetContainer *jetContainer, AliTrackContainer *trackContainer, Float_t jetPt, Int_t level)
        // Algorithm = 0
        // beta = 0.
        // option = 0
        // measure = 0
        // betasd = 0
        // zcut = 0.1
        // sotdro = 1
        Float_t Result_NSub1=10.0;
        Float_t Result_NSub2=-100.0;
        Float_t deltaR = -10.0;
        AliFJWrapper *fFastJetWrapperNsubjettines;
        fFastJetWrapperNsubjettines = new AliFJWrapper("fFastJetWrapperNsubjettines","fFastJetWrapperNsubjettines");
        fFastJetWrapperNsubjettines->Clear();
        fFastJetWrapperNsubjettines->SetR(fJetRadius[iJetRadius]);
        fFastJetWrapperNsubjettines->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
        fFastJetWrapperNsubjettines->SetRecombScheme(fastjet::RecombinationScheme::E_scheme); // fFastJetWrapper->SetRecombScheme(fastjet::RecombinationScheme::pt_scheme);
        fFastJetWrapperNsubjettines->SetStrategy(fastjet::Strategy::Best);
        fFastJetWrapperNsubjettines->SetGhostArea(fGhostArea);
        fFastJetWrapperNsubjettines->SetAreaType(fastjet::AreaType::active_area);
        fFastJetWrapperNsubjettines->SetMaxRap(0.9);
        fFastJetWrapperNsubjettines->SetMinJetPt(0.15);
        //
        // Add jet constituents to wrapper for the nsubjettiness calculation
        for(Int_t i = 0; i < nConstituents; i++)
        {
          fastjet::PseudoJet &constituentSubjettiness = constituents[i];
          Int_t trackIndex = constituentSubjettiness.user_index();
          AliESDtrack* trackNsubjettines = fESD->GetTrack(trackIndex);
          fFastJetWrapperNsubjettines->AddInputVector(trackNsubjettines->Px(), trackNsubjettines->Py(), trackNsubjettines->Pz(), trackNsubjettines->E(), i);
        }
        Result_NSub1 = fFastJetWrapperNsubjettines->AliFJWrapper::NSubjettiness(1,0,fJetRadius[iJetRadius], 0., 0, 0, 0., 0.1, 1);
        Result_NSub2 = fFastJetWrapperNsubjettines->AliFJWrapper::NSubjettiness(2,0,fJetRadius[iJetRadius], 0., 0, 0, 0., 0.1, 1);
        deltaR       = fFastJetWrapperNsubjettines->AliFJWrapper::NSubjettiness(2,0,fJetRadius[iJetRadius], 0., 2, 0, 0., 0.1, 1);
        Float_t tau2to1 = (Result_NSub1 != 0.) ? Result_NSub2 / Result_NSub1 : 0.; // TODO
        delete fFastJetWrapperNsubjettines;
        //
        // Apply rotations
        Float_t thetaTrack = -1.0;
        Float_t phiTrack = -1.0;
        Float_t rotationMatrix[3][3];
        Float_t rotationMatrix2D[2][2];
        Float_t jetUnitVector[3] = {Float_t(TMath::Cos(jetphi) / TMath::CosH(jeteta)), Float_t(TMath::Sin(jetphi) / TMath::CosH(jeteta)), Float_t(TMath::SinH(jeteta) / TMath::CosH(jeteta))};
        Float_t magPt = TMath::Sqrt((jetUnitVector[0] * jetUnitVector[0]) + (jetUnitVector[1] * jetUnitVector[1]));
        Float_t cosTheta = jetUnitVector[2];
        Float_t sinTheta = magPt;
        Float_t cosPhi = TMath::Cos(jetphi);
        Float_t sinPhi = TMath::Sin(jetphi);
        //
        rotationMatrix[0][0] = -1.0 * cosTheta * cosPhi;
        rotationMatrix[0][1] = -1.0 * cosTheta * sinPhi;
        rotationMatrix[0][2] = sinTheta;
        rotationMatrix[1][0] = sinPhi;
        rotationMatrix[1][1] = -1.0 * cosPhi;
        rotationMatrix[1][2] = 0.;
        rotationMatrix[2][0] = sinTheta * cosPhi;
        rotationMatrix[2][1] = sinTheta * sinPhi;
        rotationMatrix[2][2] = cosTheta;
        //
        Float_t principleMatrix[2][2];
        for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 2; j++) {
            principleMatrix[i][j] = 0.0;
          }
        }
        //
        for(Int_t i = 0; i < nConstituents; i++)
        {
          fastjet::PseudoJet &constituent = constituents[i];
          Int_t trackIndex = constituent.user_index();
          AliESDtrack* track = fESD->GetTrack(trackIndex);
          Double_t closestPar[3];
          GetExpecteds(track,closestPar);
          SetCutBitsAndSomeTrackVariables(track,0);
          if (track->Pt() < 0.15 || track->Pt() >= 1000.) continue;
          fJetConstituentLabels.push_back(trackIndex);
          //
          // Reject particle wrt efficiency matrix
          // Double_t eff = 1e-5;
          // Double_t effTOF = 1e-5;
          // if (fEffMatrixGenPos) {
          //   eff = GetTrackEfficiency(closestPar[1], fPtot, fEta, setting, fSign);
          //   effTOF = GetTrackEfficiency(closestPar[1], fPtot, fEta, setting, fSign, kTRUE);
          // }
          // if (eff != 1.0) {
          //   if (fRandom.Rndm() > eff) continue;
          // }
          Double_t nsubtrPx = track->Px();
          Double_t nsubtrPy = track->Py();
          Double_t nsubtrPz = track->Pz();
          Double_t nsubtrE  = track->E();
          double normalisation_factor = (1.0 / (nsubtrE * jetMass));
          Float_t pxRotated = (rotationMatrix[0][0] * nsubtrPx) + (rotationMatrix[0][1] * nsubtrPy) + (rotationMatrix[0][2] * nsubtrPz);
          Float_t pyRotated = (rotationMatrix[1][0] * nsubtrPx) + (rotationMatrix[1][1] * nsubtrPy) + (rotationMatrix[1][2] * nsubtrPz);
          Float_t pzRotated = (rotationMatrix[2][0] * nsubtrPx) + (rotationMatrix[2][1] * nsubtrPy) + (rotationMatrix[2][2] * nsubtrPz);
          principleMatrix[0][0] += normalisation_factor * pxRotated * pxRotated;
          principleMatrix[0][1] += normalisation_factor * pxRotated * pyRotated;
          principleMatrix[1][0] += normalisation_factor * pyRotated * pxRotated;
          principleMatrix[1][1] += normalisation_factor * pyRotated * pyRotated;
        }

        Float_t principleMatrixTrace = principleMatrix[0][0] + principleMatrix[1][1];
        Float_t PrinciplMatrixDeterminant = (principleMatrix[0][0] * principleMatrix[1][1]) - (principleMatrix[0][1] * principleMatrix[1][0]);
        Float_t eigenValue1 = 0.5 * (principleMatrixTrace + TMath::Sqrt(principleMatrixTrace * principleMatrixTrace - 4 * PrinciplMatrixDeterminant));
        Float_t eigenValue2 = 0.5 * (principleMatrixTrace - TMath::Sqrt(principleMatrixTrace * principleMatrixTrace - 4 * PrinciplMatrixDeterminant));
        Float_t planarFlowJet = (4.0 * PrinciplMatrixDeterminant) / (principleMatrixTrace * principleMatrixTrace);

        Float_t eigenVector1[2];
        Float_t eigenVector2[2];
        if (principleMatrix[1][0] == 0.0 || principleMatrix[0][1] == 0.0) {
          eigenVector1[0] = principleMatrix[0][0];
          eigenVector1[1] = principleMatrix[1][1];
          eigenVector2[0] = principleMatrix[0][0];
          eigenVector2[1] = principleMatrix[1][1];
        }
        else {
          eigenVector1[0] = eigenValue1 - principleMatrix[1][1];
          eigenVector1[1] = principleMatrix[1][0];
          eigenVector2[0] = principleMatrix[0][1];
          eigenVector2[1] = eigenValue2 - principleMatrix[0][0];
        }
        if (eigenValue1 < eigenValue2) {
          eigenVector1[0] = eigenVector2[0];
          eigenVector1[1] = eigenVector2[1];
        }
        Float_t theta = ( TMath::Abs(eigenVector1[0])>0.0001) ? TMath::ATan(eigenVector1[1] / eigenVector1[0]) : -100.; // TODO
        if (theta < 0) theta += TMath::Pi();
        rotationMatrix2D[0][0] = TMath::Cos(theta);
        rotationMatrix2D[0][1] = TMath::Sin(theta);
        rotationMatrix2D[1][0] = -TMath::Sin(theta);
        rotationMatrix2D[1][1] = TMath::Cos(theta);
        //
        // Rorate all other tracks
        if (fCollisionType==1){
          int indexCounter = 0.0;
          for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
            AliESDtrack* track = fESD->GetTrack(iTrack);
            if (TMath::Abs(track->Eta()) > 0.9) continue;
            if (!track->GetInnerParam()) continue;
            if (!(track->GetTPCsignalN()>0)) continue;
            SetCutBitsAndSomeTrackVariables(track,0);
            if (!GetSystematicClassIndex(fTrackCutBits,setting)) continue;
            //
            if (fMCEvent){
              Int_t lab = TMath::Abs(track->GetLabel());
              Bool_t bPrim = fMCStack->IsPhysicalPrimary(lab);
              fIsMCPileup = IsFromPileup(lab);
              if (fIsMCPileup) continue;
              if (setting>0 && !bPrim) continue;
            }
            //
            if (track->Pt() < fTrackMinPt || TMath::Abs(track->Eta()) >= particleEtaCut) continue;
            Int_t trackLabel = track->GetLabel();
            //
            // Check if in jet
            Float_t isInJet = 0.0;
            for (Int_t iConstituent = 0; iConstituent < nConstituents; iConstituent++) {
              // if (fJetConstituentLabels[iConstituent] == iTrack) std::cout << fJetConstituentLabels[iConstituent] << " ---- " << trackLabel << " --> " << iTrack << std::endl;
              if (fJetConstituentLabels[iConstituent] == iTrack) {  // Is this correct?
                isInJet = 1.0;
                break;
              }
            }
            //
            // Rotate tracks outside of jets
            Float_t pxRotated = (rotationMatrix[0][0] * track->Px()) + (rotationMatrix[0][1] * track->Py()) + (rotationMatrix[0][2] * track->Pz());
            Float_t pyRotated = (rotationMatrix[1][0] * track->Px()) + (rotationMatrix[1][1] * track->Py()) + (rotationMatrix[1][2] * track->Pz());
            Float_t pzRotated = (rotationMatrix[2][0] * track->Px()) + (rotationMatrix[2][1] * track->Py()) + (rotationMatrix[2][2] * track->Pz());
            Float_t pxRotatedPrincipleAxis = (rotationMatrix2D[0][0] * pxRotated) + (rotationMatrix2D[0][1] * pyRotated);
            Float_t pyRotatedPrincipleAxis = (rotationMatrix2D[1][0] * pxRotated) + (rotationMatrix2D[1][1] * pyRotated);
            if (TMath::Abs(pxRotated)>0.001 && TMath::Abs(pyRotated)>0.001 && TMath::Abs(pzRotated)>0.001){
              thetaTrack = TMath::ACos(pzRotated / TMath::Sqrt((pxRotated * pxRotated) + (pyRotated * pyRotated) + (pzRotated * pzRotated)));
              phiTrack = TMath::ATan2(pyRotatedPrincipleAxis, pxRotatedPrincipleAxis);
            }
            Float_t particleDeltaR = TMath::Sqrt(TMath::Power(RelativePhi(track->Phi(),jetphi),2)+TMath::Power((track->Eta()-jeteta),2));
            //
            fShapesVar_Particles_E[indexCounter] = track->E();
            fShapesVar_Particles_pT[indexCounter] = track->Pt();
            fShapesVar_Particles_Phi[indexCounter] = phiTrack;
            fShapesVar_Particles_Theta[indexCounter] = thetaTrack;
            fShapesVar_Particles_InJet[indexCounter] = isInJet;
            fShapesVar_Particles_DeltaR[indexCounter] = particleDeltaR;
            fShapesVar_Particles_NRPhi[indexCounter] = track->Phi();
            fShapesVar_Particles_Eta[indexCounter] = track->Eta();
            fShapesVar_Particles_dEdx[indexCounter] = track->GetTPCsignal();
            indexCounter++;
          }
        }
        //
        (*fTreeSRedirector)<<"jetsFJ"<<
        "gid=" << fEventGID << // global event ID
        "pileupbit=" << fPileUpBit << // pileup selection bit
        "syst=" << setting << // syst setting
        "cent=" << fCentrality << // centrality
        "ijet=" << ijet << // jet ID in an event --> sorted by pt
        "nparticles=" << nTracksSelected << // all particles passing the event/track selection
        "nconst=" << nConstituents << // number of constituents in a jet
        "njets=" << nJets << // number of jets in a given event
        //
        "spher=" << fSpherocity << // event spherocity
        "flat=" << fFlatenicity << // event flattenicity
        "psi2=" << fEP_2_Psi << // event shape 2nd harmonic
        "psi3=" << fEP_3_Psi << // event shape 3rd harmonic
        //
        "jetetacut=" << jetAbsEtaCut << //abs eta cut for jet
        "jetradius=" << fJetRadius[iJetRadius] << // jet Radius
        "maxpt=" << leadingPt << // leading particle pt in jet
        "jetptsub=" << jetptsub << // rho subtracted jet pT --> (pt - rho*Area)
        "rhofj=" << fRhoFJ << // event rho
        "nsub1=" << Result_NSub1 <<
        "nsub2=" << Result_NSub2 <<
        "deltar=" << deltaR <<
        "planarflowjet="<< planarFlowJet <<
        //
        "jetpx=" << jetpx <<
        "jetpy=" << jetpy <<
        "jetpz=" << jetpz <<
        "jetE=" << jetE <<
        "jetpt=" << jetpt <<
        "jetphi=" << jetphi <<
        "jeteta=" << jeteta <<
        "jetarea=" << jetArea <<
        "jetmass=" << jetMass ;
        //
        if (fCollisionType==1){
          (*fTreeSRedirector)<<"jetsFJ"<<
          "particles_E.=" << &fShapesVar_Particles_E <<
          "particles_pT.=" << &fShapesVar_Particles_pT <<
          "particles_Phi.=" << &fShapesVar_Particles_Phi <<
          "particles_Theta.=" << &fShapesVar_Particles_Theta <<
          "particles_InJet.=" << &fShapesVar_Particles_InJet <<
          "particles_DeltaR.=" << &fShapesVar_Particles_DeltaR <<
          "particles_NRPhi.=" << &fShapesVar_Particles_NRPhi <<
          "particles_Eta.=" << &fShapesVar_Particles_Eta <<
          "particles_dEdx.=" << &fShapesVar_Particles_dEdx ;
        }
        (*fTreeSRedirector)<<"jetsFJ"<<"\n";

        for(Int_t i = 0; i < nConstituents; i++)
        {
          fastjet::PseudoJet &constituent = constituents[i];
          Int_t trackIndex = constituent.user_index();
          AliESDtrack* trackConst = fESD->GetTrack(trackIndex);
          Bool_t bPrim = kFALSE;
          Int_t iPart = -10;
          if (fMCEvent){
            Int_t lab = TMath::Abs(trackConst->GetLabel());
            bPrim = fMCStack->IsPhysicalPrimary(lab);
            //
            // Identify particle wrt pdg code
            AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(lab); // TParticle *trackMC  = fMCStack->Particle(lab);
            Int_t pdg = trackMCgen->Particle()->GetPdgCode();                     // Int_t pdg = trackMC->GetPdgCode();
            //
            if (TMath::Abs(pdg) == kPDGpi) { iPart = 0; } // select pi
            if (TMath::Abs(pdg) == kPDGka) { iPart = 1; } // select ka
            if (TMath::Abs(pdg) == kPDGpr) { iPart = 2; } // select pr
            if (TMath::Abs(pdg) == kPDGde) { iPart = 3; } // select de
            if (TMath::Abs(pdg) == kPDGel) { iPart = 4; } // select el
            if (iPart == -10) continue;
          }
          //
          //Track cuts start
          fDEdxEl=-100;  fDEdxPi=-100;  fDEdxKa=-100;  fDEdxPr=-100;  fDEdxDe=-100;
          fSigmaEl=-100; fSigmaPi=-100; fSigmaKa=-100; fSigmaPr=-100; fSigmaDe=-100;
          fTrackCutBits=0;  // reset the bits for the next track
          //
          // --------------------------------------------------------------
          //      Get relevant track info and set cut bits
          // --------------------------------------------------------------
          //
          Double_t closestPar[3];
          GetExpecteds(trackConst,closestPar);
          SetCutBitsAndSomeTrackVariables(trackConst,0);
          Int_t tpcNcls = trackConst->GetTPCncls();
          UShort_t tpcFindableCls = trackConst->GetTPCNclsF();
          UShort_t tpcSharedCls   = trackConst->GetTPCnclsS();
          Double_t tofSignalTunedOnData = trackConst->GetTOFsignalTunedOnData();
          Double_t length    = trackConst->GetIntegratedLength();
          Double_t tofSignal = trackConst->GetTOFsignal();
          Double_t beta = -.05;
          if((length > 0) && (tofSignal > 0)) beta = length / 2.99792458e-2 / tofSignal;
          //
          (*fTreeSRedirector)<<"jetsFJconst"<<
          "constlabel=" << trackIndex <<
          "nsub1=" << Result_NSub1 <<
          "nsub2=" << Result_NSub2 <<
          "deltar=" << deltaR <<
          "ijet=" << ijet <<
          "syst=" << setting << //  syst setting
          "jetRadius=" << fJetRadius[iJetRadius] << // jet Radius
          "jetEtaCut=" << jetAbsEtaCut << //abs eta cut for jet
          "jetpx=" << jetpx << // global event ID
          "jetpy=" << jetpy << // global event ID
          "jetpz=" << jetpz << // global event ID
          "jetE=" << jetE << // global event ID
          "jetpt=" << jetpt << //  global event ID
          "jetphi=" << jetphi << //  global event ID
          "jeteta=" << jeteta << //  global event ID
          "jetArea=" << jetArea << //jet area
          "maxpt=" << leadingPt <<
          "nparticles=" << nTracksSelected << // all particles passing the event/track selection
          "nConst=" << nConstituents << //  global event ID
          "nJets=" << nJets << //  global event ID
          "jetptsub=" << jetptsub << //bg sub jet pt (pt - rho*Area)
          "rhoFJ=" << fRhoFJ; //event rho
          if (fMCEvent){
            (*fTreeSRedirector)<<"jetsFJconst"<<
            "prim=" << bPrim <<
            "part=" << iPart ;
          }
          (*fTreeSRedirector)<<"jetsFJconst"<<
          "gid=" << fEventGID << // global event ID
          "cutBit=" << fTrackCutBits << // Systematic Cuts
          "dEdx=" << fTPCSignal << // dEdx of the track
          "sign=" << fSign << // charge
          "ptot=" << fPtot << // TPC momentum
          "p=" << fPVertex << // momentum at vertex
          "pT=" << fPt << // transverse momentum
          "eta=" << fEta << // eta
          "phi=" << fPhi << // phi
          "cent=" << fCentrality << // centrality
          "pileupbit=" << fPileUpBit ;
          if (fFillJetsBG==2){
            (*fTreeSRedirector)<<"jetsFJconst"<<
            "tpcFindableCls=" << tpcFindableCls << // number of findable clusters
            "tpcSharedCls=" << tpcSharedCls << // number of shared clusters
            "tpcSignalN=" << fTrackTPCSignalN << // number of cl used in dEdx
            "lengthInActiveZone=" << fTrackLengthInActiveZone << // track length in active zone
            "tofSignalTOD=" << tofSignalTunedOnData <<
            "beta=" << beta <<
            "tofSignal=" << tofSignal <<
            //
            "dcaxy=" << fTrackDCAxy << // dca cut on xy plane
            "dcaz=" << fTrackDCAz << // dca cut along z
            "ncltpc=" << fNcl << // number of clusters
            "cRows=" << fTrackTPCCrossedRows << // crossed Rows in TPC
            "chi2tpc=" << fTrackChi2TPC <<  // TPC chi2
            "missCl=" << fMissingCl <<  // fraction of missing clusters
            //
            "dEdxMeanEl=" << fDEdxEl << //mean dEdx for electrons
            "dEdxSigmaEl=" << fSigmaEl << //sigma dEdx for electrons
            "dEdxMeanPi=" << fDEdxPi <<
            "dEdxSigmaPi=" << fSigmaPi <<
            "dEdxMeanKa=" << fDEdxKa <<
            "dEdxSigmaKa=" << fSigmaKa <<
            "dEdxMeanPr=" << fDEdxPr <<
            "dEdxSigmaPr=" << fSigmaPr <<
            "dEdxMeanDe=" << fDEdxDe <<
            "dEdxSigmaDe=" << fSigmaDe <<
            //
            "pitofpid=" << fNSigmasPiTOF <<
            "katofpid=" << fNSigmasKaTOF <<
            "prtofpid=" << fNSigmasPrTOF <<
            "eltofpid=" << fNSigmasElTOF <<  // nsigma TPC for electrons
            "detofpid=" << fNSigmasDeTOF <<
            "eltpcpid=" << fNSigmasElTPC <<  // nsigma TPC for electrons
            "pitpcpid=" << fNSigmasPiTPC <<
            "katpcpid=" << fNSigmasKaTPC <<
            "prtpcpid=" << fNSigmasPrTPC <<
            "detpcpid=" << fNSigmasDeTPC <<
            "closestTPCPIDtype=" << closestPar[1] << //particle type
            "closestTPCPIDmass=" << closestPar[2] ; //particle mass
          }
          (*fTreeSRedirector)<<"jetsFJconst"<<"\n";
        }
      } // end of jet loop
      delete fFastJetWrapper;
    } // endd of jet radius loop
  } // end of systematic setting loop

}
//________________________________________________________________________
void AliAnalysisTaskTIdentityPID::FindJetsFJGen()
{

  //
  if (fUseCouts) std::cout << " Info::marsland: ===== In the FindJetsFJGen ===== " << std::endl;
  //
  //
  // -----------------------------------------------------------------------------------------
  // ----------------------------   MC generated pure MC particles  --------------------------
  // -----------------------------------------------------------------------------------------
  //
  //
  // Find jets
  //
  // Create jetwrapper with the same settings used in FindJetsEMC
  int nJetRadiusBins = 3;
  float fTrackMinPt = 0.15;
  float fGhostArea = 0.005;
  float bgJetAbsEtaCut = 0.7;           // fixed
  float bgJetRadius = 0.2;              // fixed
  std::vector<float> fJetRadius{0.2,0.4,0.6};
  //
  for (int iJetRadius=0; iJetRadius<nJetRadiusBins; iJetRadius++){
    Float_t jetAbsEtaCut = 0.9-fJetRadius[iJetRadius];   // fixed
    double particleEtaCut = 0.9;
    //
    //SOME CODE FROM NIMA
    AliFJWrapper *fFastJetWrapperGen;
    fFastJetWrapperGen = new AliFJWrapper("fFastJetWrapperGen","fFastJetWrapperGen");
    fFastJetWrapperGen->Clear();
    fFastJetWrapperGen->SetR(fJetRadius[iJetRadius]);
    fFastJetWrapperGen->SetAlgorithm(fastjet::JetAlgorithm::antikt_algorithm);
    fFastJetWrapperGen->SetRecombScheme(fastjet::RecombinationScheme::E_scheme);
    fFastJetWrapperGen->SetStrategy(fastjet::Strategy::Best);
    fFastJetWrapperGen->SetGhostArea(fGhostArea);
    fFastJetWrapperGen->SetAreaType(fastjet::AreaType::active_area);
    fFastJetWrapperGen->SetMaxRap(0.9);
    fFastJetWrapperGen->SetMinJetPt(0.15);
    std::vector<int> trackTTIndex;
    trackTTIndex.clear();
    //
    std::vector<fastjet::PseudoJet> particlesEmbeddedSubtracted; //will be filled with your subtracted event
    std::vector<fastjet::PseudoJet> particlesEmbedded; //fill this with your event

    for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
    {
      // Select real trigger event and reject other pile up vertices for PbPb
      fIsMCPileup = IsFromPileup(iTrack);
      if (fIsMCPileup) continue;
      AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
      if (!trackMCgen) continue;
      if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
      //
      // Acceptance cut
      Double_t ptotMCgen = trackMCgen->Pt();
      Double_t etaMCgen  = trackMCgen->Eta();
      Bool_t etaAccMaxWindow = (etaMCgen>=-0.9  && etaMCgen<=0.9);
      Bool_t momAccMaxWindow = (ptotMCgen>=0.15 && ptotMCgen<=1000);
      if (!(etaAccMaxWindow && momAccMaxWindow) ) continue;
      fFastJetWrapperGen->AddInputVector(trackMCgen->Px(), trackMCgen->Py(), trackMCgen->Pz(), trackMCgen->E(), iTrack);
      particlesEmbedded.push_back(fastjet::PseudoJet(trackMCgen->Px(), trackMCgen->Py(), trackMCgen->Pz(), trackMCgen->E()));

    }
    //
    // background jet definitions
    fastjet::JetMedianBackgroundEstimator bgE;
    fastjet::Selector selectorBG = !fastjet::SelectorNHardest(2) * fastjet::SelectorAbsEtaMax(bgJetAbsEtaCut) * fastjet::SelectorPtRange(fTrackMinPt, 1000.0); //set the max eta cut on the estimator, then get rid of 2 highest pt jets
    bgE.set_selector(selectorBG);
    fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, bgJetRadius, fastjet::E_scheme, fastjet::Best); //define the kT jet finding which will do the average background estimation
    fastjet::GhostedAreaSpec ghostSpecBG(particleEtaCut, 1, fGhostArea); //this ghost area might be too small and increase processing time too much
    fastjet::AreaDefinition areaDefBG(fastjet::active_area_explicit_ghosts, ghostSpecBG);
    fastjet::ClusterSequenceArea cluster_seq_BG(particlesEmbedded, jetDefBG, areaDefBG);
    std::vector<fastjet::PseudoJet> jetsBG = sorted_by_pt(selectorBG(cluster_seq_BG.inclusive_jets())); //find the kT jets
    if (jetsBG.size() > 0) {
      bgE.set_jets(jetsBG);  // give the kT jets to the background estimator
      fRhoFJ = bgE.rho();
    }
    //
    // start of background jet loop
    if (fFillJetsBG==1){
      for (Int_t ijet=0; ijet<Int_t(jetsBG.size()); ijet++) {
        fastjet::PseudoJet jet = jetsBG[ijet];
        Float_t jetpt = jet.pt();
        Float_t jetphi = jet.phi();
        Float_t jeteta = jet.eta();
        Float_t jetArea = jet.area();
        Float_t jetptsub = jetpt - fRhoFJ*jetArea;
        Int_t nJets = jetsBG.size();
        (*fTreeSRedirector)<<"jetsFJBGGen"<<
        "ijet=" << ijet <<
        "gid=" << fEventGID << // global event ID
        "jetRadiusBG=" << bgJetRadius << // jet Radius
        "jetEtaCutBG=" << bgJetAbsEtaCut << //abs eta cut for jet
        "nJets=" << nJets << // global event ID
        "jetpt=" << jetpt << // global event ID
        "jetphi=" << jetphi << // global event ID
        "jeteta=" << jeteta << // global event ID
        "jetptsub=" << jetptsub << //bg sub jet pt (pt - rho*Area)
        "rhoFJ=" << fRhoFJ << //event rho
        "jetArea=" << jetArea << //jet area
        "cent=" << fCentrality  <<  //  centrality
        "pileupbit=" << fPileUpBit <<
        "\n";
      } // end of background jet loop
    }


    fFastJetWrapperGen->Run();
    std::vector<fastjet::PseudoJet> jets = fFastJetWrapperGen->GetInclusiveJets();
    auto nJets = jets.size();
    //
    Float_t Result_NSub1=10.0;
    Float_t Result_NSub2=-100.0;
    Float_t deltaR = -10.0;
    // Result_NSub1 = fFastJetWrapperGen->AliFJWrapper::NSubjettiness(1,0.,fJetRadius[iJetRadius], 0., 0, 0, 0., 0.1, 1);
    // Result_NSub2 = fFastJetWrapperGen->AliFJWrapper::NSubjettiness(2,0.,fJetRadius[iJetRadius], 0., 0, 0, 0., 0.1, 1);
    // deltaR       = fFastJetWrapperGen->AliFJWrapper::NSubjettiness(2,0.,fJetRadius[iJetRadius], 0., 2, 0, 0., 0.1, 1);

    for (Int_t ijet=0; ijet<Int_t(nJets); ijet++){
      //
      // get the jet object
      fastjet::PseudoJet jet = jets[ijet];
      if (jet.pt() < fTrackMinPt || jet.perp() > 1000.0 || TMath::Abs(jet.eta()) >= jetAbsEtaCut) continue;
      fhasAcceptedFJjet = 1;
      //
      // get the jet constituents
      std::vector<fastjet::PseudoJet> constituents(fFastJetWrapperGen->GetJetConstituents(ijet));
      Int_t nConstituents = constituents.size();
      std::vector<fastjet::PseudoJet> sorted_constituents = sorted_by_pt(constituents);
      auto leadingPt = sorted_constituents[0].perp();
      Float_t jetpx = jet.px();
      Float_t jetpy = jet.py();
      Float_t jetpz = jet.pz();
      Float_t jetE  = jet.E();
      Float_t jetpt = jet.pt();
      Float_t jetphi = jet.phi();
      Float_t jeteta = jet.eta();
      Float_t jetArea = jet.area();
      Float_t jetptsub = jetpt - fRhoFJ*jetArea;
      //
      (*fTreeSRedirector)<<"jetsFJGen"<<
      // "nsub1=" << Result_NSub1 <<
      // "nsub2=" << Result_NSub2 <<
      // "deltar=" << deltaR <<
      "ijet=" << ijet <<
      "jetEtaCut=" << jetAbsEtaCut << //abs eta cut for jet
      "gid=" << fEventGID << //  global event ID
      "jetRadius=" << fJetRadius[iJetRadius] << // jet Radius
      "rhoFJ=" << fRhoFJ << //event rho
      "jetpx=" << jetpx <<
      "jetpy=" << jetpy <<
      "jetpz=" << jetpz <<
      "jetE=" << jetE <<
      "jetpt=" << jetpt <<
      "jetphi=" << jetphi <<
      "jeteta=" << jeteta <<
      "jetArea=" << jetArea << //jet area
      "maxpt=" << leadingPt <<
      "nConst=" << nConstituents <<
      "nJets=" << nJets <<
      "jetptsub=" << jetptsub << // bg sub jet pt (pt - rho*Area)
      "cent=" << fCentrality << // centrality
      "\n";

      for(Int_t i = 0; i < nConstituents; i++)
      {
        fastjet::PseudoJet &constituent = constituents[i];
        Int_t trackIndex = constituent.user_index();
        AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(trackIndex);
        Int_t pdg  = trackMCgen->Particle()->GetPdgCode();
        //
        // select particle of interest
        Int_t iPart = -10;
        if (TMath::Abs(pdg) == kPDGpi) iPart = 0;
        if (TMath::Abs(pdg) == kPDGka) iPart = 1;
        if (TMath::Abs(pdg) == kPDGpr) iPart = 2;
        if (TMath::Abs(pdg) == kPDGde) iPart = 3;
        if (TMath::Abs(pdg) == kPDGel) iPart = 4;
        if (iPart == -10) continue;
        //
        // apply primary track and acceptance cuts
        Double_t ptotMCgen = trackMCgen->P();
        Double_t etaMCgen  = trackMCgen->Eta();
        Double_t ptMCgen   = trackMCgen->Pt();
        Float_t phiMCGen  = trackMCgen->Phi();
        //
        (*fTreeSRedirector)<<"jetsFJconstGen"<<
        // "nsub1=" << Result_NSub1 <<
        // "nsub2=" << Result_NSub2 <<
        // "deltar=" << deltaR <<
        "ijet=" << ijet <<
        "jetRadius=" << fJetRadius[iJetRadius] << // jet Radius
        "jetEtaCut=" << jetAbsEtaCut << //abs eta cut for jet
        "jetpx=" << jetpx << // global event ID
        "jetpy=" << jetpy << // global event ID
        "jetpz=" << jetpz << // global event ID
        "jetE=" << jetE << // global event ID
        "jetpt=" << jetpt << // global event ID
        "jetphi=" << jetphi << // global event ID
        "jeteta=" << jeteta << // global event ID
        "jetArea=" << jetArea << // jet area
        "maxpt=" << leadingPt <<
        "nConst=" << nConstituents << // global event ID
        "nJets=" << nJets << // global event ID
        "jetptsub=" << jetptsub << //bg sub jet pt (pt - rho*Area)
        "rhoFJ=" << fRhoFJ << //event rho
        //
        "gid=" << fEventGID << // global event ID
        "part=" << iPart << // particle type
        "sign=" << fSign << // charge
        "p=" << ptotMCgen << // momentum at vertex
        "pT=" << ptMCgen << // transverse momentum
        "eta=" << etaMCgen << // eta
        "phi=" << phiMCGen << // phi
        "cent=" << fCentrality << // centrality
        "\n";


      }
    } // end of jet loop
    delete fFastJetWrapperGen;




  } // ======= end of track loop for generated particles =======


}
//________________________________________________________________________
Int_t AliAnalysisTaskTIdentityPID::MakeEventPlane(Int_t doEP_Psi2, Int_t doEP_Psi3, Int_t setting)
{

  TVector3 Qvec_eta_pos, Qvec_eta_neg;
  vector< vector<TVector3> > vec_TV3_Qvec_eta;
  vec_TV3_Qvec_eta.resize(2);
  vec_TV3_Qvec_eta[0].clear();
  vec_TV3_Qvec_eta[1].clear();
  Int_t nTracksTPCDefaultSelection = 0;
  if (setting == -1){ // for MC gen level event plane
    if (fUseCouts) std::cout << " Info::marsland: ===== In the MakeEventPlane from MC gen info ===== " << std::endl;
    for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
    {
      AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
      if (!trackMCgen) continue;
      if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
      Float_t pTMCgen   = trackMCgen->Pt();
      Float_t pxMCgen   = trackMCgen->Px();
      Float_t pyMCgen   = trackMCgen->Py();
      Float_t etaMCgen  = trackMCgen->Eta();
      if(pTMCgen < 0.15) continue;
      if(pTMCgen > 3.0 ) continue;
      if((fabs(etaMCgen)) > fEtaMax) continue;
      nTracksTPCDefaultSelection++;
      //
      // use all tracks and accumulate them in a vector for the
      if(etaMCgen < fEtaMax && etaMCgen > 0.1)
      {
        Qvec_eta_pos.SetXYZ(pxMCgen,pyMCgen,0.0);
        vec_TV3_Qvec_eta[0].push_back(Qvec_eta_pos);
      }
      if(etaMCgen > -fEtaMax && etaMCgen < -0.1)
      {
        Qvec_eta_neg.SetXYZ(pxMCgen,pyMCgen,0.0);
        vec_TV3_Qvec_eta[1].push_back(Qvec_eta_neg);
      }
    }
  } else { // for data event plane
    if (fUseCouts) std::cout << " Info::marsland: ===== In the MakeEventPlane from ESD tracks ===== " << std::endl;
    for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
      AliESDtrack* track = fESD->GetTrack(iTrack);
      if (!fESDtrackCuts->AcceptTrack(track)) continue;
      if (!track->GetInnerParam()) continue;
      if (!(track->GetTPCsignalN()>0)) continue;
      SetCutBitsAndSomeTrackVariables(track,0);
      if(fPt < 0.15) continue;
      if(fPt > 3.0 ) continue;
      if((fabs(fEta)) > fEtaMax) continue;
      nTracksTPCDefaultSelection++;
      //
      // use all tracks and accumulate them in a vector for the
      if(fEta < fEtaMax && fEta > 0.1)
      {
        Qvec_eta_pos.SetXYZ(fPx,fPy,0.0);
        vec_TV3_Qvec_eta[0].push_back(Qvec_eta_pos);
      }
      if(fEta > -fEtaMax && fEta < -0.1)
      {
        Qvec_eta_neg.SetXYZ(fPx,fPy,0.0);
        vec_TV3_Qvec_eta[1].push_back(Qvec_eta_neg);
      }

    }
  }
  // check if the event is full
  if (vec_TV3_Qvec_eta[0].size() <1 ||  vec_TV3_Qvec_eta[1].size() < 1) return 0;
  if (fCollisionType == 1) fCentrality = nTracksTPCDefaultSelection;

  if(doEP_Psi2)
  {
    Double_t Psi_full_r  = 0.0;
    TVector3 TV3_sum_Qvec_eta[2];
    std::vector<Double_t> Psi_pos_neg = {0.0,0.0};
    Double_t Qvec_correction_qx[2] = {0.0,0.0};
    Double_t Qvec_correction_qy[2] = {0.0,0.0};
    Double_t harmonic = 2.0; // level of harmonic
    //
    // get values for Q-vector correction
    // if (h2D_input_qx_qy_for_EP_pos_eta && setting == -1){
    //   Qvec_correction_qx[0] = h2D_input_qx_qy_for_EP_pos_eta -> GetMean(1);
    //   Qvec_correction_qy[0] = h2D_input_qx_qy_for_EP_pos_eta -> GetMean(2);
    //   Qvec_correction_qx[1] = h2D_input_qx_qy_for_EP_neg_eta -> GetMean(1);
    //   Qvec_correction_qy[1] = h2D_input_qx_qy_for_EP_neg_eta -> GetMean(2);
    // } else {
    //   Qvec_correction_qx[0] = 0;
    //   Qvec_correction_qy[0] = 0;
    //   Qvec_correction_qx[1] = 0;
    //   Qvec_correction_qy[1] = 0;
    // }
    for(Int_t i_eta_pos_neg = 0; i_eta_pos_neg < 2; i_eta_pos_neg++)
    {
      TV3_sum_Qvec_eta[i_eta_pos_neg].SetXYZ(0.0,0.0,0.0);
      Double_t Psi_nom = 0.0;
      Double_t Psi_den = 0.0;
      Double_t track_phi = 0.0;
      Int_t    N_tracks = 0;
      Double_t weight  = 1.0;
      for(Int_t i_Qvec = 0; i_Qvec < (Int_t)vec_TV3_Qvec_eta[i_eta_pos_neg].size(); i_Qvec++)
      {
        TV3_sum_Qvec_eta[i_eta_pos_neg] += vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec];
        weight    = vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec].Perp(); // pt weight
        track_phi = vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec].Phi(); // track azimuthal angle
        Psi_nom += ((weight*TMath::Sin(harmonic*track_phi)) - Qvec_correction_qy[i_eta_pos_neg]);
        Psi_den += ((weight*TMath::Cos(harmonic*track_phi)) - Qvec_correction_qx[i_eta_pos_neg]);
        N_tracks++;

      }
      if (i_eta_pos_neg == 0) fEP_ntracks_pos = N_tracks;
      if (i_eta_pos_neg == 1) fEP_ntracks_neg = N_tracks;
      if (i_eta_pos_neg == 1 && fEP_ntracks_neg>0) {
        fEP_2_Qx_neg = Psi_nom/N_tracks;
        fEP_2_Qy_neg = Psi_den/N_tracks;
      }
      if (i_eta_pos_neg == 0  && fEP_ntracks_pos>0) {
        fEP_2_Qx_pos = Psi_nom/N_tracks;
        fEP_2_Qy_pos = Psi_den/N_tracks;
      }
      //
      // if (h2D_qx_qy_for_EP_pos_eta){
      //   if(i_eta_pos_neg == 0) h2D_qx_qy_for_EP_pos_eta->Fill(Psi_den/N_tracks,Psi_nom/N_tracks);
      //   if(i_eta_pos_neg == 1) h2D_qx_qy_for_EP_neg_eta->Fill(Psi_den/N_tracks,Psi_nom/N_tracks);
      // }
      Psi_pos_neg[i_eta_pos_neg] = TMath::RadToDeg()*TMath::ATan2(Psi_nom,Psi_den)/harmonic;
    }
    //
    // To be filled in the histogram
    fEP_2_Psi_pos = Psi_pos_neg[0];
    fEP_2_Psi_neg = Psi_pos_neg[1];
    if(Psi_pos_neg[0] - Psi_pos_neg[1] > 90.0) Psi_pos_neg[1]  -= 180.0;
    else
    {
      if(Psi_pos_neg[1] - Psi_pos_neg[0] > 90.0) Psi_pos_neg[0]  -= 180.0;
    }
    //
    // final observable per event --> Phi2
    Psi_full_r = (Psi_pos_neg[0] + Psi_pos_neg[1])/2.0;
    if(Psi_full_r > 90.0)  Psi_full_r -= 180.0;
    if(Psi_full_r < -90.0) Psi_full_r += 180.0;
    fEP_2_Psi = Psi_full_r;

  }
  //
  //
  if(doEP_Psi3)
  {
    TVector3 TV3_sum_Qvec_eta_Psi3[2];
    Double_t Psi_full_r_Psi3  = 0.0;
    std::vector<Double_t> Psi_pos_neg_Psi3 = {0.0,0.0};
    Double_t Qvec_correction_qx_Psi3[2] = {0.0};
    Double_t Qvec_correction_qy_Psi3[2] = {0.0};

    Double_t harmonic = 3.0;
    // Qvec_correction_qx_Psi3[0] = h2D_input_qx_qy_for_EP_pos_eta_v3 -> GetMean(1);
    // Qvec_correction_qy_Psi3[0] = h2D_input_qx_qy_for_EP_pos_eta_v3 -> GetMean(2);
    // Qvec_correction_qx_Psi3[1] = h2D_input_qx_qy_for_EP_neg_eta_v3 -> GetMean(1);
    // Qvec_correction_qy_Psi3[1] = h2D_input_qx_qy_for_EP_neg_eta_v3 -> GetMean(2);
    // Qvec_correction_qx_Psi3[0] = 0;
    // Qvec_correction_qy_Psi3[0] = 0;
    // Qvec_correction_qx_Psi3[1] = 0;
    // Qvec_correction_qy_Psi3[1] = 0;

    for(Int_t i_eta_pos_neg = 0; i_eta_pos_neg < 2; i_eta_pos_neg++)
    {
      TV3_sum_Qvec_eta_Psi3[i_eta_pos_neg].SetXYZ(0.0,0.0,0.0);
      Double_t Psi_nom = 0.0;
      Double_t Psi_den = 0.0;
      Double_t weight  = 1.0;
      Double_t track_phi = 0.0;
      Int_t    N_tracks = 0;
      for(Int_t i_Qvec = 0; i_Qvec < (Int_t)vec_TV3_Qvec_eta[i_eta_pos_neg].size(); i_Qvec++)
      {
        TV3_sum_Qvec_eta_Psi3[i_eta_pos_neg] += vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec];
        weight    = vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec].Perp(); // pt weight
        track_phi = vec_TV3_Qvec_eta[i_eta_pos_neg][i_Qvec].Phi(); // track azimuthal angle
        Psi_nom += ((weight*TMath::Sin(harmonic*track_phi))- Qvec_correction_qy_Psi3[i_eta_pos_neg]);
        Psi_den += ((weight*TMath::Cos(harmonic*track_phi))- Qvec_correction_qx_Psi3[i_eta_pos_neg]);
        N_tracks++;
      }
      Psi_pos_neg_Psi3[i_eta_pos_neg] = TMath::RadToDeg()*TMath::ATan2(Psi_nom,Psi_den)/harmonic;
      //
      if (i_eta_pos_neg == 1 && fEP_ntracks_neg>0) {
        fEP_3_Qx_neg = Psi_nom/N_tracks;
        fEP_3_Qy_neg = Psi_den/N_tracks;
      }
      if (i_eta_pos_neg == 0  && fEP_ntracks_pos>0) {
        fEP_3_Qx_pos = Psi_nom/N_tracks;
        fEP_3_Qy_pos = Psi_den/N_tracks;
      }

    }
    //
    // to fill histogram
    fEP_3_Psi_neg = Psi_pos_neg_Psi3[0];
    fEP_3_Psi_pos = Psi_pos_neg_Psi3[1];
    if(Psi_pos_neg_Psi3[0] - Psi_pos_neg_Psi3[1] > 60.0) Psi_pos_neg_Psi3[1]  -= 120.0;
    else
    {
      if(Psi_pos_neg_Psi3[1] - Psi_pos_neg_Psi3[0] > 60.0) Psi_pos_neg_Psi3[0]  -= 120.0;
    }
    //
    // final observable per event --> Phi3
    Psi_full_r_Psi3 = (Psi_pos_neg_Psi3[0] + Psi_pos_neg_Psi3[1])/2.0;
    if(Psi_full_r_Psi3 > 60.0)  Psi_full_r_Psi3 -= 120.0;
    if(Psi_full_r_Psi3 < -60.0) Psi_full_r_Psi3 += 120.0;


    fEP_3_Psi = Psi_full_r_Psi3;

  }
  return 1;

}
//______________________________________________________________________________
void AliAnalysisTaskTIdentityPID::GetFlatenicity()
{

  if (fUseCouts) std::cout << " Info::marsland: ===== In the GetFlatenicity ===== " << std::endl;
  // Get VZERO Information for multiplicity later
  AliVVZERO *lVV0 = fESD->GetVZEROData();
  //
  // Flatenicity calculation
  const Int_t nRings = 4;
  const Int_t nSectors = 8;
  Float_t minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
  Float_t maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
  Float_t maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
  Float_t minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};

  // Grid
  const Int_t nCells = nRings * 2 * nSectors;
  Float_t RhoLattice[nCells];
  Float_t multLattice[nCells];
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    RhoLattice[iCh] = 0.0;
    multLattice[iCh] = 0.0;
  }

  Int_t nringA = 0;
  Int_t nringC = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    Float_t detaV0 = -1;
    Float_t mult = lVV0->GetMultiplicity(iCh);
    if (iCh < 32) { // V0C
      if (iCh < 8) {
        nringC = 0;
      } else if (iCh >= 8 && iCh < 16) {
        nringC = 1;
      } else if (iCh >= 16 && iCh < 24) {
        nringC = 2;
      } else {
        nringC = 3;
      }
      detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
    } else { // V0A
      if (iCh < 40) {
        nringA = 0;
      } else if (iCh >= 40 && iCh < 48) {
        nringA = 1;
      } else if (iCh >= 48 && iCh < 56) {
        nringA = 2;
      } else {
        nringA = 3;
      }
      detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
    }
    RhoLattice[iCh] = mult / detaV0; // needed to consider the different eta coverage
    multLattice[iCh] = mult;
  }

  Float_t mRho = 0;
  //   Float_t multRho = 0;
  Float_t flatenicity = -1;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    mRho += RhoLattice[iCh];
    //     multRho += multLattice[iCh];
  }
  Float_t multV0Mdeta = mRho;
  //   Float_t multV0M = multRho;

  // average activity per cell
  mRho /= (1.0 * nCells);

  // get sigma
  Double_t sRho_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  Float_t sRho = TMath::Sqrt(sRho_tmp);
  if (mRho > 0) {
    fFlatenicityScaled = TMath::Sqrt(multV0Mdeta) * sRho / mRho; // scaling by absolute tot mult
    fFlatenicity = sRho / mRho;
  } else {
    fFlatenicity = -1;
    fFlatenicityScaled = -1;
  }

}
//______________________________________________________________________________
void AliAnalysisTaskTIdentityPID::GetFlatenicityMC()
{

  if (fUseCouts) std::cout << " Info::marsland: ===== In the GetFlatenicityMC ===== " << std::endl;
  // Flatenicity calculation
  const Int_t nRings = 8;
  Float_t maxEta[nRings] = {-3.2, -2.7, -2.2, -1.7, 5.1, 4.5, 3.9, 3.4};
  Float_t minEta[nRings] = {-3.7, -3.2, -2.7, -2.2, 4.5, 3.9, 3.4, 2.8};

  const Int_t nSectors = 8;
  Float_t PhiBins[nSectors + 1];
  Float_t deltaPhi = (2.0 * TMath::Pi()) / (1.0 * nSectors);
  for (int i_phi = 0; i_phi < nSectors + 1; ++i_phi) {
    PhiBins[i_phi] = 0;
    if (i_phi < nSectors) {
      PhiBins[i_phi] = i_phi * deltaPhi;
    } else {
      PhiBins[i_phi] = 2.0 * TMath::Pi();
    }
  }

  // Grid
  const Int_t nCells = nRings * nSectors;
  Float_t RhoLattice[nCells];
  Float_t multLattice[nCells];
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    RhoLattice[iCh] = 0.0;
    multLattice[iCh] = 0.0;
  }

  Int_t nMult = 0;
  for (Int_t i = 0; i < fMCEvent->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMCEvent->GetTrack(i);
    if (!particle)
    continue;
    if (!fMCEvent->IsPhysicalPrimary(i)) continue;
    if (particle->Pt() <= 0.0) continue;
    if (TMath::Abs(particle->Charge()) < 0.1) continue;
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMCEvent)) continue;

    Double_t phi = particle->Phi();
    Double_t eta = particle->Eta();

    Int_t i_segment = 0;
    for (int i_eta = 0; i_eta < nRings; ++i_eta) {

      for (int i_phi = 0; i_phi < nSectors; ++i_phi) {

        if (eta >= minEta[i_eta] && eta < maxEta[i_eta] && phi >= PhiBins[i_phi] && phi < PhiBins[i_phi + 1])
        {
          nMult++;
          RhoLattice[i_segment] += 1.0;
          multLattice[i_segment] += 1.0;
        }
        i_segment++;
      }
    }
  }

  Int_t i_seg = 0;
  for (int i_eta = 0; i_eta < nRings; ++i_eta) {
    for (int i_phi = 0; i_phi < nSectors; ++i_phi) {
      Float_t deltaEta = TMath::Abs(maxEta[i_eta] - minEta[i_eta]);
      RhoLattice[i_seg] /= deltaEta;
      i_seg++;
    }
  }

  Float_t mRho = 0;
  Float_t multRho = 0;
  Float_t flatenicity = -1;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    mRho    += RhoLattice[iCh];
    multRho += multLattice[iCh];
  }
  Float_t multiplicityV0M = mRho;
  Float_t multV0M = multRho;

  // average activity per cell
  mRho /= (1.0 * nCells);

  // get sigma
  Float_t sRho_tmp = 0;
  for (Int_t iCh = 0; iCh < nCells; iCh++) {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  Float_t sRho = TMath::Sqrt(sRho_tmp);

  if (mRho > 0) {
    fFlatenicityScaled = TMath::Sqrt(1.0 * multV0M) * sRho / mRho;
    fFlatenicity = sRho / mRho;
  } else {
    sRho = -1;
    fFlatenicity = -1;
    fFlatenicityScaled = -1;
  }

}
//______________________________________________________________________________
Double_t AliAnalysisTaskTIdentityPID::ComputeSpherocity(Int_t setting)
{
  if (fUseCouts) std::cout << " Info::marsland: ===== In the ComputeSpherocity ===== " << std::endl;
  Int_t ntracksLoop = (setting == -1) ? fMCEvent->GetNumberOfTracks() : fESD->GetNumberOfTracks();
  Float_t pFull = 0;
  Float_t Spherocity = 1000;
  Int_t minMult = 10;
  vector<Float_t> pt(ntracksLoop);
  vector<Float_t> phi(ntracksLoop);
  Int_t GoodTracks = 0;
  Float_t sumpt = 0;
  if (setting == -1){
    for (Int_t iTrack = 0; iTrack < fMCEvent->GetNumberOfTracks(); iTrack++)
    {
      AliMCParticle *trackMCgen = (AliMCParticle *)fMCEvent->GetTrack(iTrack);
      if (!trackMCgen) continue;
      if (!fMCStack->IsPhysicalPrimary(iTrack)) continue;
      Float_t pTMCgen   = trackMCgen->Pt();
      Float_t phiMCGen  = trackMCgen->Phi();
      Float_t etaMCgen  = trackMCgen->Eta();
      if((fabs(pTMCgen) < 0.15)) continue;
      if((fabs(etaMCgen)) > 0.8) continue;
      pt[GoodTracks] = pTMCgen;
      phi[GoodTracks] = phiMCGen;
      sumpt += pt[GoodTracks];
      GoodTracks++;
    }
  } else {
    for (Int_t iTrack = 0; iTrack < fESD->GetNumberOfTracks(); iTrack++) {
      AliESDtrack* track = fESD->GetTrack(iTrack);
      if (!fESDtrackCuts->AcceptTrack(track)) continue;
      if (!track->GetInnerParam()) continue;
      if (!(track->GetTPCsignalN()>0)) continue;
      Float_t pTreal   = track->Pt();
      Float_t phireal  = track->Phi();
      Float_t etareal  = track->Eta();
      if((fabs(pTreal) < 0.15)) continue;
      if((fabs(etareal)) > 0.8) continue;
      pt[GoodTracks] = pTreal;
      phi[GoodTracks] = phireal;
      sumpt += pt[GoodTracks];
      GoodTracks++;
    }
  }
  if (GoodTracks < minMult) return -10.0;
  //
  //Getting thrust
  Float_t stepSize = 0.1;
  for(Int_t i = 0; i < 360/stepSize; ++i){
    Float_t numerador = 0;
    Float_t phiparam  = 0;
    Float_t nx = 0;
    Float_t ny = 0;
    phiparam=( (TMath::Pi()) * i * 0.1 ) / 180.; // parametrization of the angle
    nx = TMath::Cos(phiparam);            // x component of an unitary vector n
    ny = TMath::Sin(phiparam);            // y component of an unitary vector n
    for(Int_t i1 = 0; i1 < GoodTracks; ++i1){
      Float_t pxA = pt[i1] * TMath::Cos( phi[i1] );
      Float_t pyA = pt[i1] * TMath::Sin( phi[i1] );
      numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
    }
    if (sumpt == 0) return -10.0;
    pFull=TMath::Power( (numerador / sumpt), 2);
    if(pFull < Spherocity) Spherocity = pFull;
  }
  if (GoodTracks >= minMult)  return ((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;
  else return -10.0;
}
