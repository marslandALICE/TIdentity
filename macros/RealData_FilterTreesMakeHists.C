/*

meld /home/marsland/Desktop/ubuntu_desktop/workdir/code/RealData_FilterTreesMakeHists.C /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_FilterTreesMakeHists.C

events->Draw("pileUp1DITS/multTPC","multTPC-pileUp1DITS<1000 && pileUp1DITS/multTPC<1.4 && pileUp1DITS/multTPC>0.005")
events->Draw("multV0/multTPC/2.55","multTPC-pileUp1DITS<1000 && pileUp1DITS/multTPC<1.4 && pileUp1DITS/multTPC>0.005","same")

clean->Draw("dEdx:ptot","abs(dEdx-50)<50&&abs(alfa)<0.5&&purity==2&&qt>0.14","")

*/



#include <TH3.h>
#include "THn.h"
#include "TCut.h"
#include "TCutG.h"
#include <TStyle.h>
#include <THnSparse.h>
#include "TMath.h"
#include "TLine.h"
#include "TVectorF.h"
#include "TH1.h"
#include "TH2.h"
#include "TF2.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TTreeStream.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TLinearFitter.h"
#include "AliNDLocalRegression.h"
#include "TF1.h"
#include "TAxis.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "AliXRDPROOFtoolkit.h"
#include "TStatToolkit.h"
#include "AliMathBase.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliTreePlayer.h"
#include "AliPID.h"
#include <fstream>
#include <iostream>
using namespace std;


void InitInitials();
void WriteHistsToFile();
void ProcessDataHists();
void ProcessCleanSamples();
void ProcessDScaledTree();
void ProcessHighPtTree();
void ProcessEventTree();
void PlotTimeSeriesPerSector(TString performanceEventList);
void PlotTimeSeriesFullPeriod(TString performanceEventList);
void CreateSplinesFromTHnSparse(Int_t icent, Int_t ieta);
Bool_t ApplyTreeSelection(Int_t syst, UInt_t cutBit);
TGraphErrors * ProfileToGraphErrors(TH2F * h2);
void SetBranchAddresses();
void ProcessExpectedHists(TString plotVS, TString cutON, Int_t systSet, TString hFile);


// ======= Modification part =======
const Int_t colors[] = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
const Int_t nParticles = 5;
const Int_t nCentBins = 9;
Bool_t   oldMap       = kFALSE;
//
const Double_t dEdxNbins = 1000;   // ??? be careful it must match with the analysis task binwidth
Double_t dEdxMin = 20;
Double_t dEdxMax = 1020;
//
const Double_t ptNbins = 150;   // mostly for real data analysis
Double_t ptotMin = 0.1;
Double_t ptotMax = 3.1;
//
const Int_t nEtaBins = 16;  // ???
Double_t etaMin = -0.8;
Double_t etaMax = 0.8;
//
const Int_t nCentDim = 10;
const Int_t nEtaDim = 17;
Float_t etaBinning[nEtaDim]   = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
Float_t centBinning[nCentDim] = { 0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
//
// Data hists
TH2F *h2Dall[nCentBins][nEtaBins];                 TString hName2Dall[nCentBins][nEtaBins];
TH2F *h2Dpos[nCentBins][nEtaBins];                 TString hName2Dpos[nCentBins][nEtaBins];
TH2F *h2Dneg[nCentBins][nEtaBins];                 TString hName2Dneg[nCentBins][nEtaBins];
//
TH2F *h2DallPrTOF[nCentBins][nEtaBins];            TString hName2DallPrTOF[nCentBins][nEtaBins];
TH2F *h2DallPrTOFPos[nCentBins][nEtaBins];         TString hName2DallPrTOFPos[nCentBins][nEtaBins];
TH2F *h2DallPrTOFNeg[nCentBins][nEtaBins];         TString hName2DallPrTOFNeg[nCentBins][nEtaBins];
//
TH2F *h2DallKaTOF[nCentBins][nEtaBins];            TString hName2DallKaTOF[nCentBins][nEtaBins];
TH2F *h2DallKaTOFPos[nCentBins][nEtaBins];         TString hName2DallKaTOFPos[nCentBins][nEtaBins];
TH2F *h2DallKaTOFNeg[nCentBins][nEtaBins];         TString hName2DallKaTOFNeg[nCentBins][nEtaBins];
//
TH2F *h2DCleanPiKineCut[nCentBins][nEtaBins];      TString hName2DCleanPiKineCut[nCentBins][nEtaBins];
TH2F *h2DCleanPrKineCut[nCentBins][nEtaBins];      TString hName2DCleanPrKineCut[nCentBins][nEtaBins];
TH2F *h2DCleanElKineCut[nCentBins][nEtaBins];      TString hName2DCleanElKineCut[nCentBins][nEtaBins];
//
TH2F *h2DCleanKaTOFTRD[nCentBins][nEtaBins];       TString hName2DCleanKaTOFTRD[nCentBins][nEtaBins];
TH2F *h2DCleanKaBayes[nCentBins][nEtaBins];        TString hName2DCleanKaBayes[nCentBins][nEtaBins];
TH2F *h2DCleanPiTight[nCentBins][nEtaBins];        TString hName2DCleanPiTight[nCentBins][nEtaBins];
TH2F *h2DCleanPiTOF[nCentBins][nEtaBins];          TString hName2DCleanPiTOF[nCentBins][nEtaBins];
//
TH2F *h2DClean[nParticles][nCentBins][nEtaBins];
TString hName2DClean[nParticles][nCentBins][nEtaBins];
//
// ----------------------------------------------------------------------------------------------------------------------
TH2F *h2DallCorr[nCentBins][nEtaBins];                 TString hName2DallCorr[nCentBins][nEtaBins];
TH2F *h2DposCorr[nCentBins][nEtaBins];                 TString hName2DposCorr[nCentBins][nEtaBins];
TH2F *h2DnegCorr[nCentBins][nEtaBins];                 TString hName2DnegCorr[nCentBins][nEtaBins];
//
TH2F *h2DallPrTOFCorr[nCentBins][nEtaBins];            TString hName2DallPrTOFCorr[nCentBins][nEtaBins];
TH2F *h2DallPrTOFPosCorr[nCentBins][nEtaBins];         TString hName2DallPrTOFPosCorr[nCentBins][nEtaBins];
TH2F *h2DallPrTOFNegCorr[nCentBins][nEtaBins];         TString hName2DallPrTOFNegCorr[nCentBins][nEtaBins];
//
TH2F *h2DallKaTOFCorr[nCentBins][nEtaBins];            TString hName2DallKaTOFCorr[nCentBins][nEtaBins];
TH2F *h2DallKaTOFPosCorr[nCentBins][nEtaBins];         TString hName2DallKaTOFPosCorr[nCentBins][nEtaBins];
TH2F *h2DallKaTOFNegCorr[nCentBins][nEtaBins];         TString hName2DallKaTOFNegCorr[nCentBins][nEtaBins];
//
TH2F *h2DCleanPiKineCutCorr[nCentBins][nEtaBins];      TString hName2DCleanPiKineCutCorr[nCentBins][nEtaBins];
TH2F *h2DCleanPrKineCutCorr[nCentBins][nEtaBins];      TString hName2DCleanPrKineCutCorr[nCentBins][nEtaBins];
TH2F *h2DCleanElKineCutCorr[nCentBins][nEtaBins];      TString hName2DCleanElKineCutCorr[nCentBins][nEtaBins];
//
TH2F *h2DCleanKaTOFTRDCorr[nCentBins][nEtaBins];       TString hName2DCleanKaTOFTRDCorr[nCentBins][nEtaBins];
TH2F *h2DCleanKaBayesCorr[nCentBins][nEtaBins];        TString hName2DCleanKaBayesCorr[nCentBins][nEtaBins];
TH2F *h2DCleanPiTightCorr[nCentBins][nEtaBins];        TString hName2DCleanPiTightCorr[nCentBins][nEtaBins];
TH2F *h2DCleanPiTOFCorr[nCentBins][nEtaBins];          TString hName2DCleanPiTOFCorr[nCentBins][nEtaBins];
//
TH2F *h2DCleanCorr[nParticles][nCentBins][nEtaBins];   TString hName2DCleanCorr[nParticles][nCentBins][nEtaBins];
// ----------------------------------------------------------------------------------------------------------------------
//
// Expecteds
TH2F *h2Expected[nParticles][nCentBins][nEtaBins];                TString hName2Expected[nParticles][nCentBins][nEtaBins];
TH2F *h2ExpectedSigma[nParticles][nCentBins][nEtaBins];           TString hName2ExpectedSigma[nParticles][nCentBins][nEtaBins];
TGraphErrors *grExpected[nParticles][nCentBins][nEtaBins];        TString grNameExpected[nParticles][nCentBins][nEtaBins];
TGraphErrors *grExpectedSigma[nParticles][nCentBins][nEtaBins];   TString grNameExpectedSigma[nParticles][nCentBins][nEtaBins];
//
//
// Data tree branches
UInt_t ftracks_cutBit=0;
ULong64_t ftracks_gid=0;
Double_t ftracks_eventtime=0;
Float_t ftracks_intrate=0;
Float_t ftracks_dEdx=0;
Int_t ftracks_sign=0;
Float_t ftracks_ptot=0;
Float_t ftracks_p=0;
Float_t ftracks_pT=0;
Float_t ftracks_eta=0;
Float_t ftracks_cent=0;
Int_t ftracks_nclTPC=0;
Float_t ftracks_dcaxy=0;
Float_t ftracks_dcaz=0;
//
// Clean Tree Branches
Int_t ffArmPodTree_purity = 0;
UInt_t ffArmPodTree_cutBit=0;
ULong64_t ffArmPodTree_gid;
Double_t ffArmPodTree_eventtime;
Float_t  ffArmPodTree_intrate;
Float_t ffArmPodTree_dEdx=0;
Int_t ffArmPodTree_sign=0;
Float_t ffArmPodTree_ptot=0;
Float_t ffArmPodTree_p=0;
Float_t ffArmPodTree_pT=0;
Float_t ffArmPodTree_eta=0;
Float_t ffArmPodTree_cent=0;
Float_t ffArmPodTree_qt;
Float_t ffArmPodTree_alfa;
Float_t ffArmPodTree_piTOFnSigma=0;
Float_t ffArmPodTree_prTOFnSigma=0;
Char_t ffArmPodTree_piFromK0=0;
Char_t ffArmPodTree_v0haspixel=0;
//
// dscaled tree
Float_t fdscaled_dcaxy = 0;
Float_t fdscaled_dcaz = 0;
UInt_t fdscaled_cutBit = 0;
Float_t fdscaled_cent = 0;
ULong64_t fdscaled_gid = 0;
Float_t fdscaled_intrate = 0;
Double_t fdscaled_eventtime = 0;
Double_t fdscaled_vz = 0;
Double_t fdscaled_sharedTPCClusters = 0;
Int_t fdscaled_nclTPC = 0;
Int_t fdscaled_tpcSignalN = 0;
Int_t fdscaled_cRows = 0;
Double_t fdscaled_lengthInActiveZone = 0;
Double_t fdscaled_phi = 0;
Double_t fdscaled_phiTPC = 0;
Double_t fdscaled_phi85 = 0;
Double_t fdscaled_itsdEdx = 0;
Double_t fdscaled_trddEdx = 0;
Float_t fdscaled_fdTPC = 0;
Float_t fdscaled_fzTPC = 0;
Float_t fdscaled_fCdd = 0;
Float_t fdscaled_fCdz = 0;
Float_t fdscaled_fCzz = 0;
Int_t fdscaled_nclITS = 0;
Int_t fdscaled_nclTRD = 0;
Double_t fdscaled_chi2TPC = 0;
Double_t fdscaled_chi2ITS = 0;
Double_t fdscaled_chi2TRD = 0;
//
// evets tree variables
Int_t fevents_run=0;
Float_t fevents_intrate=0;
Double_t fevents_bField=0;
ULong64_t fevents_gid=0;
Double_t fevents_timestamp=0;
ULong64_t fevents_triggerMask=0;
Double_t fevents_vz=0;
Double_t fevents_tpcvz=0;
Double_t fevents_spdvz=0;
Int_t fevents_tpcMult=0;
Short_t fevents_eventMult=0;
Int_t fevents_eventMultESD=0;
Int_t fevents_nTracksStored=0;
Int_t fevents_primMult=0;
Int_t fevents_tpcTrackBeforeClean=0;
Int_t fevents_itsTracklets=0;
TVectorF *fevents_centrality=0x0;
TVectorF *fevents_tZeroMult=0x0;
TVectorF *fevents_vZeroMult=0x0;
TVectorF *fevents_itsClustersPerLayer=0x0;
TVectorF *fevents_trackCounters=0x0;
TVectorF *fevents_trackdEdxRatio=0x0;
TVectorF *fevents_trackNcl=0x0;
TVectorF *fevents_trackChi2=0x0;
TVectorF *fevents_trackMatchEff=0x0;
TVectorF *fevents_trackTPCCountersZ=0x0;
TVectorF *fevents_phiCountA=0x0;
TVectorF *fevents_phiCountC=0x0;
TVectorF *fevents_phiCountAITS=0x0;
TVectorF *fevents_phiCountCITS=0x0;
TVectorF *fevents_phiCountAITSOnly=0x0;
TVectorF *fevents_phiCountCITSOnly=0x0;
TVectorF *fevents_tpcVertexInfo=0x0;
//
// highPt tree varianles
ULong64_t fhighPt_gid=0;
Int_t fhighPt_selectionPtMask=0;
Float_t fhighPt_centralityF=0;
Int_t fhighPt_runNumber=0;
Int_t fhighPt_evtTimeStamp=0;
Int_t fhighPt_evtNumberInFile=0;
Float_t fhighPt_Bz=0;
Int_t fhighPt_IRtot=0;
Int_t fhighPt_IRint2=0;
Int_t fhighPt_mult=0;
TObjString *fhighPt_triggerClass=0;
AliESDtrack *fhighPt_esdTrack=0;
TObjString *fhighPt_fileName=0;
AliESDVertex *fhighPt_vtxESD=0;
//
//
const Int_t fnCutBins=10;
Int_t fCutArr[fnCutBins];
Int_t fSystSet=0;
TString fPlotVS="";
TString fCutON="";
TString inputDataTree    = "tracks";
TString inputCleanTree   = "fArmPodTree";
TString inputEventTree   = "events";
TString inputDscaledTree = "dscaled";
TString inputHighPtTree  = "highPt";
TString inputHists     = "cleanHists";
TString inputHistsFile = "";
TFile *fdata=NULL;
TFile *fhist=NULL;
TFile *fdEdxMap=NULL;
TTree *dataTree=NULL, *armtree=NULL, *eventtree=NULL, *dscaltree=NULL, *highPttree=NULL;
THnSparse *fhnExpected=NULL;
TTreeSRedirector *treeStream=0, *histStream=0, *histStreamCorr, *debugStream=0, *splinesStream=0;
TH1D  *hEta=0x0, *hCent=0x0, *hMom=0x0;
TCutG *pionCutG=0, *antiProtonCutG=0, *protonCutG=0;
TStopwatch timer;
AliNDLocalRegression *hdEdxAShifttMNTglDist_meanGFitAll=NULL;
AliNDLocalRegression *hdEdxAShifttMNTglDist_meanGFitNoPileUp=NULL;
//
TH2F *htracks_dcaxy2D=NULL;
TH2F *htracks_dcaz2D=NULL;
TH2F *htracks_nclTPC2D=NULL;
//
TH2F *hevents_pileUpV02D=NULL;
TH2F *hevents_pileUpITS2D=NULL;
TH2F *htracks_shiftM_nPileUpPrim=NULL;
TH1F *hevents_pileUpV01D=NULL;
TH1F *hevents_pileUpITS1D=NULL;
TH1F *hevents_ITSTPCeff=NULL;
//
TH1F *hdscaled_sharedTPCClusters1D=NULL;
TH1F *hdscaled_tpcSignalN1D=NULL;
TH1F *hdscaled_lengthInActiveZone1D=NULL;
TH1F *hdscaled_cRows1D=NULL;
TH1F *hdscaled_nclITS1D=NULL;
TH1F *hdscaled_nclTRD1D=NULL;
TH1F *hdscaled_chi2TPC1D=NULL;
TH1F *hdscaled_chi2ITS1D=NULL;
TH1F *hdscaled_chi2TRD1D=NULL;
TH1F *hdscaled_vz1D=NULL;
TH1F *hdscaled_dcaxy1D=NULL;
TH1F *hdscaled_dcaz1D=NULL;
TH1F *hdscaled_nclTPC1D=NULL;
//
TH2F *hhighPt_dEdxPtot=NULL;
//
//
enum trackCutBit
{
  kNCrossedRowsTPC60=0,
  kNCrossedRowsTPC80=1,
  kNCrossedRowsTPC100=2,
  kMaxChi2PerClusterTPC3=3,
  kMaxChi2PerClusterTPC4=4,
  kMaxChi2PerClusterTPC5=5,
  kMaxDCAToVertexXYPtDepSmall=6,
  kMaxDCAToVertexXYPtDep=7,
  kMaxDCAToVertexXYPtDepLarge=8,
  kVertexZSmall=9,
  kVertexZ=10,
  kVertexZLarge=11,
  kEventVertexZSmall=12,
  kEventVertexZ=13,
  kEventVertexZLarge=14,
  kRequireITSRefit=15,
  kPixelRequirementITS=16,
  kNewITSCut=17,
  kActiveZoneSmall=18,
  kActiveZone=19,
  kActiveZoneLarge=20,
  kTPCSignalNSmall=21,
  kTPCSignalN=22,
  kTPCSignalNLarge=23,
};
TString parName[nParticles] = {"El", "Pi", "Ka", "Pr", "De"};
enum parType
{
  kEl=0,
  kPi=1,
  kKa=2,
  kPr=3,
  kDe=4
};



/*

kPileUp = 1000, 2000, 4000, 1000000;
kTimeSeriesEff = 0.80, 0.85, 0;
kVzPileup = +, -, All;



*/

Int_t    fInPileUpLow     = 0;
Int_t    fInPileUpHigh    = 0;
Double_t fInTimeSeriesEff = 0.;
Int_t    fInVzPileup      = 0;
Int_t    fPeriod          = 0;
//
TString dEdxMapStr = "";
TString timeSeriesMap = "";
TString MapDirName = "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists";
// TString MapDirName = "/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp";
// Double_t testEntries  = 2000000;
Double_t testEntries  = -1.;

void RealData_FilterTreesMakeHists(TString plotVS, TString cutON, Int_t systSet, TString dFile, Int_t period, Int_t kPileUpLow, Int_t kPileUpHigh, Double_t kTimeSeriesEff, Int_t kVzPileup)
{
  //
  // Produce all hists form MC and Real data
  // plotVS: "ptot", "pT", "p"
  //

  /*

  meld /home/marsland/Desktop/ubuntu_desktop/workdir/code/RealData_FilterTreesMakeHists.C /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_FilterTreesMakeHists.C


  cd /home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists
  aliroot -l
  TString file = "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN2/LHC15o_pass1_NoSelection_4Runs_20072019/mergedRuns/merged/000246148/Sub_0/filteredTreesAndTracks/AnalysisResults_tree2.root"
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/code/RealData_FilterTreesMakeHists.C+
  RealData_FilterTreesMakeHists("ptot","p",0,file,0  ,-100000,1000,    0.88, 0)
  RealData_FilterTreesMakeHists("ptot","p",0,file,0  ,-100000,100000,  0,    0)
  RealData_FilterTreesMakeHists("ptot","p",0,file,0  , 2000,8000,      0.88, 0)

  RealData_FilterTreesMakeHists("ptot","p",0,file,0  ,-50000,1000,    0.88, 0)   // full correction and selection
  RealData_FilterTreesMakeHists("ptot","p",0,file,0  ,-50000,50000, 0,    0)     // no correction
  ProcessExpectedHists("AnalysisResults_hist.root",)


  aliroot -l
  TString file = "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/ATO-465/data/LHC18q/LHC18q_pass1_AllRuns/merged/000295908/filteredTrees/AnalysisResults_filteredTree10.root"
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/code/RealData_FilterTreesMakeHists.C+
  RealData_FilterTreesMakeHists("ptot","p",0,file,1  ,1, 1000, 0.85, 0)   // full correction and selection


  cd /home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists
  aliroot -l
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/code/RealData_FilterTreesMakeHists.C+
  ProcessExpectedHists("ptot","p",0, "AnalysisResults_hist.root")


  */
  //
  fInPileUpLow     = kPileUpLow;
  fInPileUpHigh    = kPileUpHigh;
  fInTimeSeriesEff = kTimeSeriesEff;
  fInVzPileup      = kVzPileup;
  cout << " fInPileUpLow     = " <<  fInPileUpLow << endl;
  cout << " fInPileUpHigh    = " <<  fInPileUpHigh << endl;
  cout << " fInTimeSeriesEff = " <<  fInTimeSeriesEff << endl;
  cout << " fInVzPileup      = " <<  fInVzPileup << endl;

  //
  fSystSet         = systSet;
  fPlotVS          = plotVS;
  fCutON           = cutON;
  fPeriod          = period;
  //
  if (fPeriod==0)  dEdxMapStr = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2015/LHC15o/pass1/dEdxFit.root",MapDirName.Data());
  if (fPeriod==1)  dEdxMapStr = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18q/pass1/dEdxFit.root",MapDirName.Data());
  if (fPeriod==2)  dEdxMapStr = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18r/pass1/dEdxFit.root",MapDirName.Data());
  //
  if (oldMap)  dEdxMapStr = Form("%s/old/dEdxFit.root",MapDirName.Data());
  else dEdxMapStr = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2015/LHC15o/pass1/dEdxFit.root",MapDirName.Data());
  //
  cout << "  ------------ > dEdx Map = " << dEdxMapStr << endl;
  //
  fdata      = TFile::Open(dFile);
  armtree    = (TTree*)fdata->Get(inputCleanTree);
  dataTree   = (TTree*)fdata->Get(inputDataTree);
  eventtree  = (TTree*)fdata->Get(inputEventTree);
  dscaltree  = (TTree*)fdata->Get(inputDscaledTree);
  highPttree = (TTree*)fdata->Get(inputHighPtTree);
  //
  //
  SetBranchAddresses();
  //
  InitInitials();
  if (fPeriod==0){
    //
    if (dataTree) dataTree->BuildIndex("gid");
    if (eventtree) {
      eventtree->BuildIndex("gid");
      dataTree->AddFriend(eventtree,"event");
    }
    if (dscaltree) {
      dscaltree->BuildIndex("gid");
      dataTree->AddFriend(dscaltree,"dscaled");
    }
    if (armtree) {
      armtree->BuildIndex("gid");
      armtree->AddFriend(eventtree,"event");
    }
    std::cout << "Indexing is finished " << dFile << std::endl;
    std::cout << "get trees from root file  ---------> " << dFile << std::endl;
    ProcessDataHists();
    ProcessEventTree();
    if (testEntries<0) ProcessCleanSamples();
    if (testEntries<0) ProcessDScaledTree();
    if (testEntries<0) ProcessHighPtTree();
  }
  if (fPeriod>0) ProcessEventTree();
  WriteHistsToFile();
  delete treeStream;
  delete debugStream;
  delete histStream;
  delete histStreamCorr;


}
//____________________________________________________________________________________________________________
void ProcessDataHists()
{

  timer.Reset(); timer.Start();
  std::cout << " ========= ProcessDataHists ========= " << std::endl;
  Double_t nTreeEntriesAll = dataTree -> GetEntries();
  Double_t nTreeEntries    = (testEntries>0) ? testEntries : nTreeEntriesAll;
  cout << " Data Tree entries = " << nTreeEntriesAll << endl;
  if (nTreeEntriesAll<10) { std::cout << " === upss data tree is empty === " << std::endl; return; }
  //
  // Loop over tree entries
  Int_t eventCount=0;
  TGraph *fgrTimeSeriesEventWeighted=NULL;
  TGraph *fgrTimeSeriesNTracksWeighted=NULL;
  for(Int_t i = 0; i < nTreeEntries; ++i)
  {
    //
    // Get Track and event information
    dataTree -> GetEntry(i);
    if(i%Int_t(nTreeEntriesAll/10) == 0) {
      cout << i << "   ftracks_dEdx = " << ftracks_dEdx      << "         ftracks_eta        = " << ftracks_eta << "  ftracks_dcaxy = " << ftracks_dcaxy << endl;
      cout << "   ftracks_gid       = " << ftracks_gid       << "    -->  fevents_gid        = " <<  fevents_gid              << endl;
      cout << "   ftracks_eventtime = " << ftracks_eventtime << "    -->  fevents_timestamp  = " <<  fevents_timestamp       << endl;
      cout << "   ftracks_cent      = " << ftracks_cent      << "    -->  fevents_centrality = " <<  (*fevents_centrality)[0] << endl;
    }
    //
    // Get the time series map for a given run
    static Int_t runCache = -1;
    if (runCache!=Int_t(fevents_run)) {
      runCache=Int_t(fevents_run);
      if (fPeriod==0) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2015/LHC15o/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==1) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18q/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==2) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18r/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      cout << "run number = " << runCache << "  ------------ > timeSeries Map = " <<  timeSeriesMap << endl;
      TFile *fQAtimeSeries = TFile::Open(timeSeriesMap);
      if (!fQAtimeSeries) {
        std::cout << " yolun acik ola --> file does not exist  " << std::endl;
      }
      fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.binMedian");   // weighting wrt event
      fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.binMedian");     // weighting wrt track
    }
    //
    // Apply event selection
    static ULong64_t gidCache = -1;
    static Int_t acceptEvent = -1;
    if (gidCache!=ULong64_t(ftracks_gid)) {
      gidCache=ULong64_t(ftracks_gid);
      //
      //  TPC ITS matching eff from time series
      Double_t fITSTPCeffEvent = fgrTimeSeriesEventWeighted  ->Eval(fevents_timestamp);
      Double_t fITSTPCeffTrack = fgrTimeSeriesNTracksWeighted->Eval(fevents_timestamp);
      //
      // PileUp selection:  chain->SetAlias("multTPCITS","((multSSD+multSDD)/2.38)");   // chain->SetAlias("cut1000","tpcMult-multTPCITS<1000");
      Double_t multSSD = (*fevents_itsClustersPerLayer)[4]+(*fevents_itsClustersPerLayer)[5];
      Double_t multSDD = (*fevents_itsClustersPerLayer)[2]+(*fevents_itsClustersPerLayer)[3];
      Double_t multTPC = fevents_tpcMult;
      Double_t pileUp1DITS = (multSSD+multSDD)/2.38;
      Double_t multV0 = fevents_vZeroMult->Sum();
      Double_t multT0 = fevents_tZeroMult->Sum();
      //
      // pile up and time series cut
      if ( (multTPC-pileUp1DITS) > fInPileUpLow  && (multTPC-pileUp1DITS) < fInPileUpHigh && fITSTPCeffTrack > fInTimeSeriesEff )  acceptEvent=1;
      else acceptEvent=-1;
      //
      // Dump debug histograms
      if (multTPC>0 && acceptEvent>0){
        if (pileUp1DITS/multTPC>0.005) hevents_pileUpV01D->Fill(multV0/multTPC/2.55);
        if (pileUp1DITS/multTPC>0.005) hevents_pileUpITS1D->Fill(pileUp1DITS/multTPC);
        hevents_pileUpV02D ->Fill(multTPC,multV0);
        hevents_pileUpITS2D->Fill(multTPC,pileUp1DITS);
        hevents_ITSTPCeff  ->Fill(fITSTPCeffTrack);
      }
    }
    //
    // Event Cuts --> pile up and time series cut
    if (acceptEvent<0) continue;
    //
    // Apply dEdx correction for a given track at a given eta
    Double_t shiftM      = 0.5*((*fevents_tpcVertexInfo)[1]+(*fevents_tpcVertexInfo)[0]);
    Double_t norm        = 1-TMath::Abs(shiftM/210.);
    Double_t primMult    = fevents_primMult;
    Double_t nPileUpPrim = ((*fevents_tpcVertexInfo)[3]+(*fevents_tpcVertexInfo)[4])/norm;
    Double_t trackTgl    = TMath::Abs(TMath::SinH(ftracks_eta));
    //
    Double_t params1[]  = {shiftM,nPileUpPrim,primMult,trackTgl};
    Double_t params0[]  = {primMult,trackTgl};
    Double_t dEdxCorr0 =  hdEdxAShifttMNTglDist_meanGFitNoPileUp->Eval(params0);
    Double_t dEdxCorr1 =  hdEdxAShifttMNTglDist_meanGFitAll     ->Eval(params1);
    Double_t dEdxCorr0ND = AliNDLocalRegression::GetCorrND(1,primMult,trackTgl+0);
    Double_t dEdxCorr1ND = AliNDLocalRegression::GetCorrND(2,shiftM,nPileUpPrim,primMult,trackTgl+0);
    Double_t dEdxCorr  =  ftracks_dEdx-50*dEdxCorr1;
    //
    // Vz pileup cut
    if ( shiftM<0 && fInVzPileup== 1 ) continue;
    if ( shiftM>0 && fInVzPileup==-1 ) continue;
    htracks_shiftM_nPileUpPrim->Fill(shiftM,nPileUpPrim);
    //
    // dump tidentree
    treeStream->GetFile()->cd();
    (*treeStream)<<"tracks"<<
    // "gid="                  << ftracks_gid          <<  //  global event ID
    // "p="                    << ftracks_p            <<  //  TPC momentum
    "intrate="              << ftracks_intrate      <<  // interaction rate
    "cutBit="               << ftracks_cutBit       <<  //  Systematic Cuts
    "dEdx="                 << ftracks_dEdx         <<  //  dEdx of the track
    "dEdxCorr="             << dEdxCorr             <<  //  dEdx of the track
    "sign="                 << ftracks_sign         <<  //  charge
    "ptot="                 << ftracks_ptot         <<  //  TPC momentum
    "pT="                   << ftracks_pT           <<
    "eta="                  << ftracks_eta          <<  //  eta
    "cent="                 << ftracks_cent         <<  //  centrality
    "dEdxCorr0="            << dEdxCorr0            <<  //  centrality
    "dEdxCorr1="            << dEdxCorr1            <<  //  4D correction using Eval inteface
    "shiftM="               << shiftM               <<  // interaction rate
    "primMult="             << primMult             <<  // interaction rate
    "tgl="                  << trackTgl             <<  // interaction rate
    "nPileUpPrim="          << nPileUpPrim          <<  // interaction rate
    //
    "pileUpCutLow="         << fInPileUpLow          <<  // interaction rate
    "pileUpCutHigh="        << fInPileUpHigh          <<  // interaction rate
    "tSeriesCut="           << fInTimeSeriesEff          <<  // interaction rate
    "vzSideCut="            << fInVzPileup          <<  // interaction rate
    // "dEdxCorr1ND="            << dEdxCorr1ND            <<  // 4D correction using ND interface
    // "dEdxCorr0ND="            << dEdxCorr0ND            <<  //  centrality
    "\n";
    //
    // dump some debug histogram
    if (ftracks_dcaxy>-10 && ftracks_dcaxy<10)   htracks_dcaxy2D ->Fill(ftracks_pT,ftracks_dcaxy);
    if (ftracks_dcaz>-10  && ftracks_dcaz<10)    htracks_dcaz2D  ->Fill(ftracks_pT,ftracks_dcaz);
    if (ftracks_nclTPC>30 && ftracks_nclTPC<170) htracks_nclTPC2D->Fill(ftracks_pT,Float_t(ftracks_nclTPC));
    //
    // apply dca and nclusters cut
    if (TMath::Abs(ftracks_dcaz)>3) continue;
    if (TMath::Abs(ftracks_dcaxy)>3) continue;
    if (ftracks_nclTPC<70) continue;
    //
    // Prepare dEdx histograms
    for (Int_t icent=0; icent<nCentBins; icent++){
      for (Int_t ieta=0; ieta<nEtaBins; ieta++){

        Bool_t etaCentString = (ftracks_eta>=etaBinning[ieta] && ftracks_eta<etaBinning[ieta+1] && ftracks_cent>=centBinning[icent] && ftracks_cent<centBinning[icent+1]);
        Bool_t cleanKaCutTOFTRD = ((ftracks_cutBit >> 26) & 1);
        Bool_t cleanKaCut = ((ftracks_cutBit >> 27) & 1);
        Bool_t cleanDeCut = ((ftracks_cutBit >> 30) & 1);
        Bool_t parPos = ftracks_sign>0.;
        Bool_t parNeg = ftracks_sign<0.;
        //
        Bool_t prTOF = ((ftracks_cutBit >> 24) & 1);
        Bool_t kaTOF = ((ftracks_cutBit >> 25) & 1);
        //
        Bool_t systCut = ApplyTreeSelection(fSystSet, ftracks_cutBit);
        Bool_t vertexPcut=kFALSE;
        if(fCutON=="pT")   vertexPcut = (ftracks_pT>=ptotMin      && ftracks_pT<=ptotMax);
        if(fCutON=="ptot") vertexPcut = (ftracks_ptot>=ptotMin    && ftracks_ptot<=ptotMax);
        if(fCutON=="p")    vertexPcut = (ftracks_p>=ptotMin && ftracks_p<=ptotMax);

        if (etaCentString && vertexPcut && systCut          ) h2Dall[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && parPos) h2Dpos[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && parNeg) h2Dneg[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        //
        if (etaCentString && vertexPcut && systCut && kaTOF)           h2DallKaTOF[icent][ieta]   ->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && kaTOF && parPos) h2DallKaTOFPos[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && kaTOF && parNeg) h2DallKaTOFNeg[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        //
        if (etaCentString && vertexPcut && systCut && prTOF)           h2DallPrTOF[icent][ieta]   ->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && prTOF && parPos) h2DallPrTOFPos[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && prTOF && parNeg) h2DallPrTOFNeg[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        //
        if (etaCentString && vertexPcut && systCut && cleanKaCut)       h2DCleanKaBayes[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && cleanKaCutTOFTRD) h2DCleanKaTOFTRD[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);

        if (etaCentString && vertexPcut && systCut && cleanDeCut)       h2DClean[4][icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && cleanKaCutTOFTRD) h2DClean[2][icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        //
        //
        //
        if (etaCentString && vertexPcut && systCut          ) h2DallCorr[icent][ieta]->Fill(ftracks_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && systCut && parPos) h2DposCorr[icent][ieta]->Fill(ftracks_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && systCut && parNeg) h2DnegCorr[icent][ieta]->Fill(ftracks_ptot,dEdxCorr);
        //
        if (etaCentString && vertexPcut && systCut && kaTOF)           h2DallKaTOFCorr[icent][ieta]   ->Fill(ftracks_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && systCut && kaTOF && parPos) h2DallKaTOFPosCorr[icent][ieta]->Fill(ftracks_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && systCut && kaTOF && parNeg) h2DallKaTOFNegCorr[icent][ieta]->Fill(ftracks_ptot,dEdxCorr);
        //
        if (etaCentString && vertexPcut && systCut && prTOF)           h2DallPrTOFCorr[icent][ieta]   ->Fill(ftracks_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && systCut && prTOF && parPos) h2DallPrTOFPosCorr[icent][ieta]->Fill(ftracks_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && systCut && prTOF && parNeg) h2DallPrTOFNegCorr[icent][ieta]->Fill(ftracks_ptot,dEdxCorr);
        //
        if (etaCentString && vertexPcut && systCut && cleanKaCut)       h2DCleanKaBayesCorr[icent][ieta]->Fill(ftracks_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && systCut && cleanKaCutTOFTRD) h2DCleanKaTOFTRDCorr[icent][ieta]->Fill(ftracks_ptot,dEdxCorr);

        if (etaCentString && vertexPcut && systCut && cleanDeCut)       h2DCleanCorr[4][icent][ieta]->Fill(ftracks_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && systCut && cleanKaCutTOFTRD) h2DCleanCorr[2][icent][ieta]->Fill(ftracks_ptot,dEdxCorr);

      }
    }

  }  // tree loop
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessDataHists DONE ========= #events = " << eventCount << std::endl;

}
//____________________________________________________________________________________________________________
void ProcessEventTree()
{

  timer.Reset(); timer.Start();
  std::cout << " ========= ProcessEventTree ========= " << std::endl;
  Double_t nTreeEntriesAll = eventtree -> GetEntries();
  cout << " Event Tree entries = " << nTreeEntriesAll << endl;
  //
  // Loop over tree entries
  Int_t eventCount=0;
  TGraph *fgrTimeSeriesEventWeighted=NULL;
  TGraph *fgrTimeSeriesNTracksWeighted=NULL;
  for(Int_t i = 0; i < nTreeEntriesAll; ++i)
  {
    //
    // Get Track and event information
    eventtree -> GetEntry(i);
    if(i%Int_t(nTreeEntriesAll/10) == 0) {
      cout << i << "   fevents_run = " << fevents_run << "    fevents_primMult = " << fevents_primMult << "    fevents_intrate = " << fevents_intrate << endl;
    }
    //
    // Get the time series map for a given run
    static Int_t runCache = -1;
    if (runCache!=Int_t(fevents_run)) {
      runCache=Int_t(fevents_run);
      TString timeSeriesMap = "";
      if (fPeriod==0) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2015/LHC15o/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==1) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18q/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==2) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18r/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      cout << "timeSeries Map = " <<  timeSeriesMap << endl;
      TFile *fQAtimeSeries = TFile::Open(timeSeriesMap);
      if (!fQAtimeSeries) {
        std::cout << " yolun acik ola --> file does not exist  " << std::endl;
      }
      fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.binMedian");   // weighting wrt event
      fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.binMedian");     // weighting wrt track
    }
    //
    //
    Double_t shiftM0 = 0.5*((*fevents_tpcVertexInfo)[1]+(*fevents_tpcVertexInfo)[0]);
    Double_t norm0   = 1-TMath::Abs(shiftM0/210.);
    Double_t nPileUpPrim0 = ((*fevents_tpcVertexInfo)[3]+(*fevents_tpcVertexInfo)[4])/norm0;
    Double_t primMult0    = fevents_primMult;
    //
    //  TPC ITS matching eff from time series
    Double_t fITSTPCeffEvent = fgrTimeSeriesEventWeighted  ->Eval(fevents_timestamp);
    Double_t fITSTPCeffTrack = fgrTimeSeriesNTracksWeighted->Eval(fevents_timestamp);
    //
    // PileUp selection:  chain->SetAlias("multTPCITS","((multSSD+multSDD)/2.38)");   // chain->SetAlias("cut1000","tpcMult-multTPCITS<1000");
    Double_t multSSD = (*fevents_itsClustersPerLayer)[4]+(*fevents_itsClustersPerLayer)[5];
    Double_t multSDD = (*fevents_itsClustersPerLayer)[2]+(*fevents_itsClustersPerLayer)[3];
    Double_t multSPD = (*fevents_itsClustersPerLayer)[0]+(*fevents_itsClustersPerLayer)[1];
    Double_t multV0 = fevents_vZeroMult->Sum();
    Double_t multT0 = fevents_tZeroMult->Sum();
    Double_t multTPC = fevents_tpcMult;
    Double_t pileUp1DITS = (multSSD+multSDD)/2.38;
    //
    // dump tidentree
    treeStream->GetFile()->cd();
    (*treeStream)<<"events"<<
    "run="            << fevents_run         <<  //  global event ID
    "timeStamp="      << fevents_timestamp   <<  //  global event ID
    "bField="         << fevents_bField      <<  //  global event ID
    "gid="            << fevents_gid         <<  //  global event ID
    "vZ="             << fevents_vz          <<  //  global event ID
    "tpcvz="          << fevents_tpcvz          <<  //  global event ID
    "spdvz="          << fevents_spdvz          <<  //  global event ID
    "eventMultESD="   << fevents_eventMultESD     <<  //  global event ID
    "cent="           << (*fevents_centrality)[0] <<
    "intrate="        << fevents_intrate     <<  // interaction rate
    "multSSD="        << multSSD             <<  //  global event ID
    "multSDD="        << multSDD             <<  // interaction rate
    "multSPD="        << multSPD             <<  //  Systematic Cuts
    "multV0="         << multV0              <<  //  dEdx of the track
    "multT0="         << multT0              <<  //  charge
    "multTPC="        << multTPC             <<  //  TPC momentum
    "pileUp1DITS="    << pileUp1DITS         <<  //  TPC momentum
    "ITSTPCeffEvent=" << fITSTPCeffEvent     <<  //  TPC momentum
    "ITSTPCeffTrack=" << fITSTPCeffTrack     <<  //  TPC momentum
    "tpcTrackBeforeClean=" << fevents_tpcTrackBeforeClean     <<  //  TPC momentum
    "nTracksStored=" << fevents_nTracksStored     <<  //  TPC momentum
    //
    "shiftM="         << shiftM0             <<  //  TPC momentum
    "norm="           << norm0               <<  //  TPC momentum
    "nPileUpPrim="    << nPileUpPrim0        <<  //  TPC momentum
    "primMult="       << primMult0           <<  //  TPC momentum
    //
    "pileUpCutLow="         << fInPileUpLow          <<  // interaction rate
    "pileUpCutHigh="        << fInPileUpHigh          <<  // interaction rate
    "tSeriesCut="           << fInTimeSeriesEff          <<  // interaction rate
    "vzSideCut="            << fInVzPileup          <<  // interaction rate
    "\n";
    eventCount++;


  }  // tree loop
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessDataHists DONE ========= #events = " << eventCount << std::endl;

}
//____________________________________________________________________________________________________________
void ProcessDScaledTree()
{


  timer.Reset(); timer.Start();
  std::cout << " ========= ProcessDScaledTree ========= " << std::endl;
  Double_t nTreeEntries = dscaltree -> GetEntries();
  if (nTreeEntries<10) { std::cout << " === upss data tree is empty === " << std::endl; return; }
  //
  // Loop over tree entries
  for(Int_t i = 0; i < nTreeEntries; ++i)
  {

    dscaltree -> GetEntry(i);
    Bool_t systCut = ApplyTreeSelection(fSystSet, fdscaled_cutBit);
    if(i%Int_t(nTreeEntries/10) == 0) cout << i << " fdscaled_cRows = " << fdscaled_cRows << "  fdscaled_nclTPC = " << fdscaled_nclTPC << "  fdscaled_vz = " << fdscaled_vz << endl;
    //
    if (fdscaled_sharedTPCClusters>0    && fdscaled_sharedTPCClusters<0.5    && systCut) hdscaled_sharedTPCClusters1D ->Fill(fdscaled_sharedTPCClusters);
    if (fdscaled_tpcSignalN>0           && fdscaled_tpcSignalN<170           && systCut) hdscaled_tpcSignalN1D        ->Fill(fdscaled_tpcSignalN);
    if (fdscaled_lengthInActiveZone>0   && fdscaled_lengthInActiveZone<170   && systCut) hdscaled_lengthInActiveZone1D->Fill(fdscaled_lengthInActiveZone);
    if (fdscaled_cRows>0                && fdscaled_cRows<170                && systCut) hdscaled_cRows1D->Fill(fdscaled_cRows);
    if (fdscaled_nclITS>-2              && fdscaled_nclITS<8                 && systCut) hdscaled_nclITS1D->Fill(fdscaled_nclITS);
    if (fdscaled_nclTRD>0               && fdscaled_nclTRD<200               && systCut) hdscaled_nclTRD1D->Fill(fdscaled_nclTRD);
    if (fdscaled_chi2TPC>0              && fdscaled_chi2TPC<10               && systCut) hdscaled_chi2TPC1D->Fill(fdscaled_chi2TPC);
    if (fdscaled_chi2ITS>0              && fdscaled_chi2ITS<20               && systCut) hdscaled_chi2ITS1D->Fill(fdscaled_chi2ITS);
    if (fdscaled_chi2TRD>0              && fdscaled_chi2TRD<20               && systCut) hdscaled_chi2TRD1D->Fill(fdscaled_chi2TRD);
    if (fdscaled_vz>-20                 && fdscaled_vz<20                    && systCut) hdscaled_vz1D->Fill(fdscaled_vz);
    if (fdscaled_dcaxy>-10              && fdscaled_dcaxy<10                 && systCut) hdscaled_dcaxy1D ->Fill(fdscaled_dcaxy);
    if (fdscaled_dcaz>-10               && fdscaled_dcaz<10                  && systCut) hdscaled_dcaz1D ->Fill(fdscaled_dcaz);
    if (fdscaled_nclTPC>30              && fdscaled_nclTPC<170               && systCut) hdscaled_nclTPC1D ->Fill(fdscaled_nclTPC);

  }  // tree loop
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessDScaledTree DONE ========= " << std::endl;

}
//____________________________________________________________________________________________________________
void ProcessHighPtTree()
{


  timer.Reset(); timer.Start();
  std::cout << " ========= ProcessHighPtTree ========= " << std::endl;
  Double_t nTreeEntries = highPttree -> GetEntries();
  if (nTreeEntries<10) { std::cout << " === upss data tree is empty === " << std::endl; return; }
  //
  // Loop over tree entries
  for(Int_t i = 0; i < nTreeEntries; ++i)
  {

    highPttree -> GetEntry(i);
    hhighPt_dEdxPtot->Fill(fhighPt_esdTrack->GetInnerParam()->GetP(),fhighPt_esdTrack->GetTPCsignal());

  }  // tree loop
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessHighPtTree DONE ========= " << std::endl;

}
//____________________________________________________________________________________________________________
void ProcessExpectedHists(TString plotVS, TString cutON, Int_t systSet, TString hFile)
{
  inputHistsFile = hFile;
  //
  // Read THnSparses
  std::cout << " ========= ProcessExpectedHists ========= " << std::endl;
  //
  Bool_t fileExist = (inputHistsFile.Contains(".root") || inputHistsFile.Contains(".list"));
  if ( !fileExist ) {
    std::cout << "This is not a proper file ---------> " << inputHistsFile << std::endl;
    return;
  }

  fhist = TFile::Open(inputHistsFile);
  if (!fhist->GetListOfKeys()->Contains("cleanHists")) return;
  TList * list  = (TList*)fhist->Get("cleanHists");
  if (list){
    fhnExpected = (THnSparse*)list->FindObject(Form("hExpected_%d",systSet));
  } else {
    std::cout << "This is not a proper file ---------> " << inputHistsFile << std::endl;
    return;
  }
  //
  //
  TString splinesFileName = Form("Splines_Syst%d_PlotVS_%s_CutON_%s.root",systSet,plotVS.Data(),cutON.Data());
  splinesStream = new TTreeSRedirector(splinesFileName,"recreate");
  splinesStream->GetFile()->cd();
  //
  // Initialize hists
  for (Int_t icent=0; icent<nCentBins; icent++){
    for (Int_t ieta=0; ieta<nEtaBins; ieta++){
      for (Int_t ipart=0; ipart<nParticles; ipart++){
        h2Expected[ipart][icent][ieta]=NULL;
        h2ExpectedSigma[ipart][icent][ieta]=NULL;
        grExpected[ipart][icent][ieta]=NULL;
        grExpectedSigma[ipart][icent][ieta]=NULL;
      }
    }
  }
  //
  // Create Expected Hists
  for (Int_t icent=0; icent<nCentBins; icent++){
    for (Int_t ieta=0; ieta<1; ieta++){

      cout << icent << "    "  << ieta << endl;
      CreateSplinesFromTHnSparse(icent,ieta);

      for (Int_t ipart=0; ipart<nParticles; ipart++){
        if(h2Expected[ipart][icent][ieta])       h2Expected[ipart][icent][ieta]->Write();
        if(h2ExpectedSigma[ipart][icent][ieta])  h2ExpectedSigma[ipart][icent][ieta]->Write();
        if(grExpected[ipart][icent][ieta])       grExpected[ipart][icent][ieta]->Write();
        if(grExpectedSigma[ipart][icent][ieta])  grExpectedSigma[ipart][icent][ieta]->Write();
      }
    }
  }
  delete splinesStream;
  //
  std::cout << " ========= ProcessExpectedHists DONE ========= " << std::endl;

}
//____________________________________________________________________________________________________________
void ProcessCleanSamples()
{

  timer.Reset(); timer.Start();
  std::cout << " ========= ProcessCleanSamples ========= " << std::endl;
  Double_t nTreeEntriesAll = armtree -> GetEntries();
  Double_t nTreeEntries    = nTreeEntriesAll;
  cout << " Data Tree entries = " << nTreeEntriesAll << endl;
  if (nTreeEntriesAll<10) { std::cout << " === upss data tree is empty === " << std::endl; return; }
  //
  // Loop over tree entries
  TGraph *fgrTimeSeriesEventWeighted=NULL;
  TGraph *fgrTimeSeriesNTracksWeighted=NULL;
  for(Int_t i = 0; i < nTreeEntries; ++i)
  {

    armtree -> GetEntry(i);
    Bool_t piFromK0 = Int_t(ffArmPodTree_piFromK0);
    Bool_t v0haspixel = Int_t(ffArmPodTree_v0haspixel);
    if(i%Int_t(nTreeEntriesAll/10) == 0) cout << i << " ffArmPodTree_dEdx = " << ffArmPodTree_dEdx << "  ffArmPodTree_eta = " << ffArmPodTree_eta << "  ffArmPodTree_cent = " << ffArmPodTree_cent << "  --->  " << piFromK0 << " --  " << v0haspixel << endl;
    //
    // Get the time series map for a given run
    static Int_t runCache = -1;
    if (runCache!=Int_t(fevents_run)) {
      runCache=Int_t(fevents_run);
      if (fPeriod==0) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2015/LHC15o/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==1) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18q/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==2) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18r/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      cout << "timeSeries Map = " <<  timeSeriesMap << endl;
      TFile *fQAtimeSeries = TFile::Open(timeSeriesMap);
      if (!fQAtimeSeries) {
        std::cout << " yolun acik ola --> file does not exist  " << std::endl;
      }
      fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.binMedian");   // weighting wrt event
      fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.binMedian");     // weighting wrt track
    }
    //
    // Apply event selection
    static ULong64_t gidCache = -1;
    static Int_t acceptEvent = -1;
    if (gidCache!=ULong64_t(ffArmPodTree_gid)) {
      gidCache=ULong64_t(ffArmPodTree_gid);
      //
      //  TPC ITS matching eff from time series
      Double_t fITSTPCeffEvent = fgrTimeSeriesEventWeighted  ->Eval(fevents_timestamp);
      Double_t fITSTPCeffTrack = fgrTimeSeriesNTracksWeighted->Eval(fevents_timestamp);
      //
      // PileUp selection:  chain->SetAlias("multTPCITS","((multSSD+multSDD)/2.38)");   // chain->SetAlias("cut1000","tpcMult-multTPCITS<1000");
      Double_t multSSD = (*fevents_itsClustersPerLayer)[4]+(*fevents_itsClustersPerLayer)[5];
      Double_t multSDD = (*fevents_itsClustersPerLayer)[2]+(*fevents_itsClustersPerLayer)[3];
      Double_t multTPC = fevents_tpcMult;
      Double_t pileUp1DITS = (multSSD+multSDD)/2.38;
      //
      // pile up and time series cut
      if ( (multTPC-pileUp1DITS) > fInPileUpLow  && (multTPC-pileUp1DITS) < fInPileUpHigh && fITSTPCeffTrack > fInTimeSeriesEff )  acceptEvent=1;
      else acceptEvent=-1;
    }
    //
    // Event Cuts --> pile up and time series cut
    if (acceptEvent<0) continue;
    //
    // Apply dEdx correction for a given track at a given eta
    Double_t shiftM      = 0.5*((*fevents_tpcVertexInfo)[1]+(*fevents_tpcVertexInfo)[0]);
    Double_t norm        = 1-TMath::Abs(shiftM/210.);
    Double_t primMult    = fevents_primMult;
    Double_t nPileUpPrim = ((*fevents_tpcVertexInfo)[3]+(*fevents_tpcVertexInfo)[4])/norm;
    Double_t trackTgl    = TMath::Abs(TMath::SinH(ffArmPodTree_eta));
    //
    Double_t params1[]  = {shiftM,nPileUpPrim,primMult,trackTgl};
    Double_t params0[]  = {primMult,trackTgl};
    Double_t dEdxCorr0 =  hdEdxAShifttMNTglDist_meanGFitNoPileUp->Eval(params0);
    Double_t dEdxCorr1 =  hdEdxAShifttMNTglDist_meanGFitAll     ->Eval(params1);
    Double_t dEdxCorr0ND = AliNDLocalRegression::GetCorrND(1,primMult,trackTgl+0);
    Double_t dEdxCorr1ND = AliNDLocalRegression::GetCorrND(2,shiftM,nPileUpPrim,primMult,trackTgl+0);
    Double_t dEdxCorr  =  ffArmPodTree_dEdx-50*dEdxCorr1;
    //
    // Vz pileup cut
    if ( shiftM<0 && fInVzPileup== 1 ) continue;
    if ( shiftM>0 && fInVzPileup==-1 ) continue;
    //
    // dump tidentree
    treeStream->GetFile()->cd();
    (*treeStream)<<"clean"<<
    // "gid="                  << ffArmPodTree_gid          <<  //  global event ID
    // "p="                    << ffArmPodTree_p            <<  //  TPC momentum
    "intrate="              << ffArmPodTree_intrate      <<  // interaction rate
    "cutBit="               << ffArmPodTree_cutBit       <<  //  Systematic Cuts
    "dEdx="                 << ffArmPodTree_dEdx         <<  //  dEdx of the track
    "dEdxCorr="             << dEdxCorr             <<  //  dEdx of the track
    "sign="                 << ffArmPodTree_sign         <<  //  charge
    "ptot="                 << ffArmPodTree_ptot         <<  //  TPC momentum
    "pT="                   << ffArmPodTree_pT           <<
    "eta="                  << ffArmPodTree_eta          <<  //  eta
    "cent="                 << ffArmPodTree_cent         <<  //  centrality
    "alfa="                 << ffArmPodTree_alfa         <<  //  centrality
    "qt="                   << ffArmPodTree_qt         <<  //  centrality
    "purity="               << ffArmPodTree_purity         <<  //  centrality
    "piFromK0="             << ffArmPodTree_piFromK0         <<  //  centrality
    "dEdxCorr0="            << dEdxCorr0            <<  //  centrality
    "dEdxCorr1="            << dEdxCorr1            <<  //  4D correction using Eval inteface
    "shiftM="               << shiftM               <<  // interaction rate
    "primMult="             << primMult             <<  // interaction rate
    "tgl="                  << trackTgl             <<  // interaction rate
    "nPileUpPrim="          << nPileUpPrim          <<  // interaction rate
    //
    "pileUpCutLow="         << fInPileUpLow          <<  // interaction rate
    "pileUpCutHigh="        << fInPileUpHigh          <<  // interaction rate
    "tSeriesCut="           << fInTimeSeriesEff          <<  // interaction rate
    "vzSideCut="            << fInVzPileup          <<  // interaction rate
    // "dEdxCorr1ND="          << dEdxCorr1ND            <<  // 4D correction using ND interface
    // "dEdxCorr0ND="          << dEdxCorr0ND            <<  //  centrality
    "\n";


    for (Int_t icent=0; icent<nCentBins; icent++){
      for (Int_t ieta=0; ieta<nEtaBins; ieta++){

        Bool_t etaCentString = (ffArmPodTree_eta>=etaBinning[ieta] && ffArmPodTree_eta<etaBinning[ieta+1] && ffArmPodTree_cent>=centBinning[icent] && ffArmPodTree_cent<centBinning[icent+1]);
        Bool_t vertexPcut=kFALSE;
        if(fCutON=="pT")   vertexPcut = (ffArmPodTree_pT>=ptotMin      && ffArmPodTree_pT<=ptotMax);
        if(fCutON=="ptot") vertexPcut = (ffArmPodTree_ptot>=ptotMin    && ffArmPodTree_ptot<=ptotMax);
        if(fCutON=="p")    vertexPcut = (ffArmPodTree_p>=ptotMin && ffArmPodTree_p<=ptotMax);
        //
        // Bool_t systCut = ApplyTreeSelection(fSystSet, ftracks_cutBit);
        // Bool_t cleanCutPrTOF = TMath::Abs(ffArmPodTree_prTOFnSigma)<2. && ffArmPodTree_qt>0.07 && ffArmPodTree_qt<0.11 && (antiProtonCutG||protonCutG);
        // Bool_t cleanCutPr    = "ffArmPodTree_qt>0.07 && ffArmPodTree_qt<0.11 && (antiProtonCutG||protonCutG)";
        // Bool_t cleanCutPi    = (TMath::Abs(ffArmPodTree_alfa)<0.5 && pionCutG);
        //
        Bool_t cleanCutPiTOF = ( TMath::Abs(ffArmPodTree_piTOFnSigma)<2. );
        Bool_t cleanCutPrTOF = ( (TMath::Abs(ffArmPodTree_prTOFnSigma)<2.) && (ffArmPodTree_qt>0.07) && (ffArmPodTree_qt<0.11) );
        Bool_t cleanCutPi    = ( TMath::Abs(ffArmPodTree_alfa)<0.5 );
        Bool_t cleanCutPr    = ( (ffArmPodTree_qt>0.07) && (ffArmPodTree_qt<0.11) );
        Bool_t cleanCutEl    = ( (ffArmPodTree_qt<0.005) && (TMath::Abs(ffArmPodTree_alfa)<0.5) );
        Bool_t piK0cut       = ( Bool_t(piFromK0) );
        Bool_t piPixelcut    = ( !Bool_t(v0haspixel) );
        Bool_t cleanCutPiKineCut = ( (ffArmPodTree_qt>0.14) && (ffArmPodTree_purity==2) );
        Bool_t cleanCutElKineCut = ( (ffArmPodTree_qt<0.01) && (ffArmPodTree_purity==2) );
        Bool_t cleanCutPrKineCut = ( (ffArmPodTree_qt<0.12) && (ffArmPodTree_qt>0.01) && (ffArmPodTree_purity==2) );
        //
        //
        if (etaCentString && vertexPcut && cleanCutPi    && piK0cut && piPixelcut) h2DCleanPiTight[icent][ieta]->Fill(ffArmPodTree_ptot,ffArmPodTree_dEdx);
        if (etaCentString && vertexPcut && cleanCutPi    && cleanCutPiTOF)         h2DCleanPiTOF[icent][ieta]->Fill(ffArmPodTree_ptot,ffArmPodTree_dEdx);
        if (etaCentString && vertexPcut && cleanCutEl)                             h2DClean[0][icent][ieta]->Fill(ffArmPodTree_ptot,ffArmPodTree_dEdx);
        if (etaCentString && vertexPcut && cleanCutPi    && piK0cut)               h2DClean[1][icent][ieta]->Fill(ffArmPodTree_ptot,ffArmPodTree_dEdx);
        if (etaCentString && vertexPcut && cleanCutPrTOF)                          h2DClean[3][icent][ieta]->Fill(ffArmPodTree_ptot,ffArmPodTree_dEdx);
        //
        if (etaCentString && vertexPcut && cleanCutPiKineCut) h2DCleanPiKineCut[icent][ieta]->Fill(ffArmPodTree_ptot,ffArmPodTree_dEdx);
        if (etaCentString && vertexPcut && cleanCutElKineCut) h2DCleanElKineCut[icent][ieta]->Fill(ffArmPodTree_ptot,ffArmPodTree_dEdx);
        if (etaCentString && vertexPcut && cleanCutPrKineCut) h2DCleanPrKineCut[icent][ieta]->Fill(ffArmPodTree_ptot,ffArmPodTree_dEdx);

        //
        if (etaCentString && vertexPcut && cleanCutPi    && piK0cut && piPixelcut) h2DCleanPiTightCorr[icent][ieta]->Fill(ffArmPodTree_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && cleanCutPi    && cleanCutPiTOF)         h2DCleanPiTOFCorr[icent][ieta]->Fill(ffArmPodTree_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && cleanCutEl)                             h2DCleanCorr[0][icent][ieta]->Fill(ffArmPodTree_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && cleanCutPi    && piK0cut)               h2DCleanCorr[1][icent][ieta]->Fill(ffArmPodTree_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && cleanCutPrTOF)                          h2DCleanCorr[3][icent][ieta]->Fill(ffArmPodTree_ptot,dEdxCorr);
        //
        //
        if (etaCentString && vertexPcut && cleanCutPiKineCut) h2DCleanPiKineCutCorr[icent][ieta]->Fill(ffArmPodTree_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && cleanCutElKineCut) h2DCleanElKineCutCorr[icent][ieta]->Fill(ffArmPodTree_ptot,dEdxCorr);
        if (etaCentString && vertexPcut && cleanCutPrKineCut) h2DCleanPrKineCutCorr[icent][ieta]->Fill(ffArmPodTree_ptot,dEdxCorr);


      }
    }

  }
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessCleanSamples DONE ========= " << std::endl;

}
//____________________________________________________________________________________________________________
void CreateSplinesFromTHnSparse(Int_t icent, Int_t ieta)
{

  //
  // Make projections form THnSparse
  fhnExpected->GetAxis(2)->SetRangeUser(centBinning[icent],centBinning[icent+1]);   // centrality
  fhnExpected->GetAxis(3)->SetRangeUser(etaBinning[ieta],etaBinning[ieta+1]);     // eta
  fhnExpected->GetAxis(5)->SetRangeUser(2.,60.);           // sigma
  fhnExpected->GetAxis(6)->SetRangeUser(25.,1000.);        // mean
  //
  // get each particle PID response
  for (Int_t ipart=0;ipart<nParticles;ipart++){
    fhnExpected->GetAxis(0)->SetRangeUser(ipart,ipart+1);
    h2Expected[ipart][icent][ieta]      = (TH2F*)fhnExpected->Projection(6,4);
    h2ExpectedSigma[ipart][icent][ieta] = (TH2F*)fhnExpected->Projection(5,4);
    grExpected[ipart][icent][ieta]      = ProfileToGraphErrors(h2Expected[ipart][icent][ieta]);
    grExpectedSigma[ipart][icent][ieta] = ProfileToGraphErrors(h2ExpectedSigma[ipart][icent][ieta]);
    h2Expected[ipart][icent][ieta]     ->SetName(Form("h2Expected_%d_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"     ,ipart,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]));
    h2ExpectedSigma[ipart][icent][ieta]->SetName(Form("h2ExpectedSigma_%d_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",ipart,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]));
    grExpected[ipart][icent][ieta]     ->SetName(Form("grExpected_%d_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"     ,ipart,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]));
    grExpectedSigma[ipart][icent][ieta]->SetName(Form("grExpectedSigma_%d_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",ipart,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]));

  }

}
//____________________________________________________________________________________________________________
void InitInitials()
{

  TString outputFileNameTree     = Form("Trees_Syst%d_PlotVS_%s_CutON_%s_%d_%d_%3.2f_%d.root"    ,fSystSet,fPlotVS.Data(),fCutON.Data(),fInPileUpLow,fInPileUpHigh,fInTimeSeriesEff,fInVzPileup);
  TString outputFileNameHist     = Form("Hists_Syst%d_PlotVS_%s_CutON_%s_%d_%d_%3.2f_%d.root"    ,fSystSet,fPlotVS.Data(),fCutON.Data(),fInPileUpLow,fInPileUpHigh,fInTimeSeriesEff,fInVzPileup);
  TString outputFileNameHistCorr = Form("HistsCorr_Syst%d_PlotVS_%s_CutON_%s_%d_%d_%3.2f_%d.root",fSystSet,fPlotVS.Data(),fCutON.Data(),fInPileUpLow,fInPileUpHigh,fInTimeSeriesEff,fInVzPileup);
  TString debugFile              = Form("Debug_Syst%d_PlotVS_%s_CutON_%s_%d_%d_%3.2f_%d.root"    ,fSystSet,fPlotVS.Data(),fCutON.Data(),fInPileUpLow,fInPileUpHigh,fInTimeSeriesEff,fInVzPileup);
  treeStream      = new TTreeSRedirector(outputFileNameTree,"recreate");
  histStream      = new TTreeSRedirector(outputFileNameHist,"recreate");
  histStreamCorr  = new TTreeSRedirector(outputFileNameHistCorr,"recreate");
  debugStream = new TTreeSRedirector(debugFile,"recreate");

  //
  // Initialise histograms to be used for binning
  //
  hevents_pileUpV02D            = new TH2F("hevents_pileUpV02D"  ,"hevents_pileUpV02D"   ,1000,0.,35000., 1000 ,0., 35000. );
  hevents_pileUpITS2D           = new TH2F("hevents_pileUpITS2D" ,"hevents_pileUpITS2D"  ,1500,0.,45000., 1500 ,0., 45000. );
  htracks_shiftM_nPileUpPrim     = new TH2F("htracks_shiftM_nPileUpPrim" ,"htracks_shiftM_nPileUpPrim"  ,1000,-200,200., 1000 ,0., 10000. );
  hevents_pileUpV01D            = new TH1F("hevents_pileUpV01D"  ,"hevents_pileUpV01D"   ,1000,0.05,1.4);
  hevents_pileUpITS1D           = new TH1F("hevents_pileUpITS1D" ,"hevents_pileUpITS1D"  ,1000,0.05,1.4);
  hevents_ITSTPCeff             = new TH1F("hevents_ITSTPCeff"   ,"hevents_ITSTPCeff"    ,1000,0.,1.2);
  //
  //
  htracks_dcaxy2D              = new TH2F("htracks_dcaxy2D" ,"htracks_dcaxy2D"  ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  htracks_dcaz2D               = new TH2F("htracks_dcaz2D"  ,"htracks_dcaz2D"   ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  htracks_nclTPC2D             = new TH2F("htracks_nclTPC2D","htracks_nclTPC2D" ,ptNbins,ptotMin,ptotMax, 140 , 30., 170.);
  //
  hdscaled_sharedTPCClusters1D  = new TH1F("hdscaled_sharedTPCClusters1D" ,"hdscaled_sharedTPCClusters1D Bins" ,200   ,  0., 0.5 );
  hdscaled_tpcSignalN1D         = new TH1F("hdscaled_tpcSignalN1D"        ,"hdscaled_tpcSignalN1D Bins"        ,170   ,  0., 170. );
  hdscaled_lengthInActiveZone1D = new TH1F("hdscaled_lengthInActiveZone1D","hdscaled_lengthInActiveZone1D Bins",170   ,  0., 170. );
  hdscaled_cRows1D              = new TH1F("hdscaled_cRows1D"             ,"hdscaled_cRows1D Bins"             ,170   ,  0., 170. );
  hdscaled_nclITS1D             = new TH1F("hdscaled_nclITS1D"            ,"hdscaled_nclITS1D Bins"            ,10    , -2., 8. );
  hdscaled_nclTRD1D             = new TH1F("hdscaled_nclTRD1D"            ,"hdscaled_nclTRD1D Bins"            ,200   ,  0., 200. );
  hdscaled_chi2TPC1D            = new TH1F("hdscaled_chi2TPC1D"           ,"hdscaled_chi2TPC1D Bins"           ,400   ,  0., 10. );
  hdscaled_chi2ITS1D            = new TH1F("hdscaled_chi2ITS1D"           ,"hdscaled_chi2ITS1D Bins"           ,400   ,  0., 20. );
  hdscaled_chi2TRD1D            = new TH1F("hdscaled_chi2TRD1D"           ,"hdscaled_chi2TRD1D Bins"           ,400   ,  0., 20. );
  hdscaled_vz1D                 = new TH1F("hdscaled_vz1D"                ,"hdscaled_vz1D Bins"                ,200   ,-20., 20. );
  hdscaled_dcaxy1D            = new TH1F("hdcaxy" ,"hdcaxy"  ,400 ,-10., 10. );
  hdscaled_dcaz1D             = new TH1F("hdcaz"  ,"hdcaz"   ,400 ,-10., 10. );
  hdscaled_nclTPC1D           = new TH1F("hnclTPC","hnclTPC" ,140 , 30., 170.);
  //
  hhighPt_dEdxPtot = new TH2F("hhighPt_dEdxPtot", "hhighPt_dEdxPtot" ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
  //
  for (Int_t icent=0; icent<nCentBins; icent++){
    for (Int_t ieta=0; ieta<nEtaBins; ieta++){

      h2Dall[icent][ieta]=NULL;   h2Dpos[icent][ieta]=NULL;   h2Dneg[icent][ieta]=NULL;
      hName2Dall[icent][ieta]=Form("h2Dall_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2Dpos[icent][ieta]=Form("h2Dpos_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2Dneg[icent][ieta]=Form("h2Dneg_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      h2Dall[icent][ieta] = new TH2F(hName2Dall[icent][ieta], hName2Dall[icent][ieta] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2Dpos[icent][ieta] = new TH2F(hName2Dpos[icent][ieta], hName2Dpos[icent][ieta] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2Dneg[icent][ieta] = new TH2F(hName2Dneg[icent][ieta], hName2Dneg[icent][ieta] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DallPrTOF[icent][ieta]=NULL;      h2DallPrTOFPos[icent][ieta]=NULL;      h2DallPrTOFNeg[icent][ieta]=NULL;
      h2DallKaTOF[icent][ieta]=NULL;      h2DallKaTOFPos[icent][ieta]=NULL;      h2DallKaTOFNeg[icent][ieta]=NULL;
      hName2DallPrTOF[icent][ieta]   =Form("h2DallPrTOF_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallPrTOFPos[icent][ieta]=Form("h2DallPrTOFPos_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallPrTOFNeg[icent][ieta]=Form("h2DallPrTOFNeg_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallKaTOF[icent][ieta]   =Form("h2DallKaTOF_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallKaTOFPos[icent][ieta]=Form("h2DallKaTOFPos_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallKaTOFNeg[icent][ieta]=Form("h2DallKaTOFNeg_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      h2DallKaTOF[icent][ieta]     = new TH2F(hName2DallPrTOF[icent][ieta]    ,hName2DallPrTOF[icent][ieta]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallKaTOFPos[icent][ieta]  = new TH2F(hName2DallPrTOFPos[icent][ieta] ,hName2DallPrTOFPos[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallKaTOFNeg[icent][ieta]  = new TH2F(hName2DallPrTOFNeg[icent][ieta] ,hName2DallPrTOFNeg[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallPrTOF[icent][ieta]     = new TH2F(hName2DallKaTOF[icent][ieta]    ,hName2DallKaTOF[icent][ieta]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallPrTOFPos[icent][ieta]  = new TH2F(hName2DallKaTOFPos[icent][ieta] ,hName2DallKaTOFPos[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallPrTOFNeg[icent][ieta]  = new TH2F(hName2DallKaTOFNeg[icent][ieta] ,hName2DallKaTOFNeg[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DCleanPiKineCut[icent][ieta]=NULL;      h2DCleanPrKineCut[icent][ieta]=NULL;      h2DCleanElKineCut[icent][ieta]=NULL;
      hName2DCleanPiKineCut[icent][ieta]   =Form("h2DCleanPiKineCut_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f" ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanPrKineCut[icent][ieta]   =Form("h2DCleanPrKineCut_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f" ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanElKineCut[icent][ieta]   =Form("h2DCleanElKineCut_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f" ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      h2DCleanPiKineCut[icent][ieta]  = new TH2F(hName2DCleanPiKineCut[icent][ieta] ,hName2DCleanPiKineCut[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanPrKineCut[icent][ieta]  = new TH2F(hName2DCleanPrKineCut[icent][ieta] ,hName2DCleanPrKineCut[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanElKineCut[icent][ieta]  = new TH2F(hName2DCleanElKineCut[icent][ieta] ,hName2DCleanElKineCut[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DCleanKaTOFTRD[icent][ieta]=NULL;  h2DCleanKaBayes[icent][ieta]=NULL;  h2DCleanPiTight[icent][ieta]=NULL;  h2DCleanPiTOF[icent][ieta]=NULL;
      hName2DCleanKaTOFTRD[icent][ieta] = Form("h2DCleanKaTOFTRD_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanKaBayes[icent][ieta] = Form("h2DCleanKaBayes_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanPiTight[icent][ieta] = Form("h2DCleanPiTight_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanPiTOF[icent][ieta] = Form("h2DCleanPiTOF_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      h2DCleanKaTOFTRD[icent][ieta]     = new TH2F(hName2DCleanKaTOFTRD[icent][ieta],hName2DCleanKaTOFTRD[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanKaBayes[icent][ieta]     = new TH2F(hName2DCleanKaBayes[icent][ieta],hName2DCleanKaBayes[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanPiTight[icent][ieta]     = new TH2F(hName2DCleanPiTight[icent][ieta],hName2DCleanPiTight[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanPiTOF[icent][ieta]     = new TH2F(hName2DCleanPiTOF[icent][ieta],hName2DCleanPiTOF[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      for (Int_t ipart=0; ipart<nParticles; ipart++)
      {
        h2DClean[ipart][icent][ieta]=NULL;
        hName2DClean[ipart][icent][ieta] = Form("h2DClean%s_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f", parName[ipart].Data(), centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
        h2DClean[ipart][icent][ieta]     = new TH2F(hName2DClean[ipart][icent][ieta],hName2DClean[ipart][icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      }
      //
      h2DallCorr[icent][ieta]=NULL;   h2DposCorr[icent][ieta]=NULL;   h2DnegCorr[icent][ieta]=NULL;
      hName2DallCorr[icent][ieta]=Form("h2DallCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DposCorr[icent][ieta]=Form("h2DposCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DnegCorr[icent][ieta]=Form("h2DnegCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      h2DallCorr[icent][ieta] = new TH2F(hName2DallCorr[icent][ieta], hName2DallCorr[icent][ieta] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DposCorr[icent][ieta] = new TH2F(hName2DposCorr[icent][ieta], hName2DposCorr[icent][ieta] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DnegCorr[icent][ieta] = new TH2F(hName2DnegCorr[icent][ieta], hName2DnegCorr[icent][ieta] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DallPrTOFCorr[icent][ieta]=NULL;      h2DallPrTOFPosCorr[icent][ieta]=NULL;      h2DallPrTOFNegCorr[icent][ieta]=NULL;
      h2DallKaTOFCorr[icent][ieta]=NULL;      h2DallKaTOFPosCorr[icent][ieta]=NULL;      h2DallKaTOFNegCorr[icent][ieta]=NULL;
      hName2DallPrTOFCorr[icent][ieta]   =Form("h2DallPrTOFCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallPrTOFPosCorr[icent][ieta]=Form("h2DallPrTOFPosCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallPrTOFNegCorr[icent][ieta]=Form("h2DallPrTOFNegCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallKaTOFCorr[icent][ieta]   =Form("h2DallKaTOFCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallKaTOFPosCorr[icent][ieta]=Form("h2DallKaTOFPosCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DallKaTOFNegCorr[icent][ieta]=Form("h2DallKaTOFNegCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      h2DallKaTOFCorr[icent][ieta]     = new TH2F(hName2DallPrTOFCorr[icent][ieta]    ,hName2DallPrTOFCorr[icent][ieta]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallKaTOFPosCorr[icent][ieta]  = new TH2F(hName2DallPrTOFPosCorr[icent][ieta] ,hName2DallPrTOFPosCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallKaTOFNegCorr[icent][ieta]  = new TH2F(hName2DallPrTOFNegCorr[icent][ieta] ,hName2DallPrTOFNegCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallPrTOFCorr[icent][ieta]     = new TH2F(hName2DallKaTOFCorr[icent][ieta]    ,hName2DallKaTOFCorr[icent][ieta]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallPrTOFPosCorr[icent][ieta]  = new TH2F(hName2DallKaTOFPosCorr[icent][ieta] ,hName2DallKaTOFPosCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallPrTOFNegCorr[icent][ieta]  = new TH2F(hName2DallKaTOFNegCorr[icent][ieta] ,hName2DallKaTOFNegCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DCleanPiKineCutCorr[icent][ieta]=NULL;      h2DCleanPrKineCutCorr[icent][ieta]=NULL;      h2DCleanElKineCutCorr[icent][ieta]=NULL;
      hName2DCleanPiKineCutCorr[icent][ieta]   =Form("h2DCleanPiKineCutCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f" ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanPrKineCutCorr[icent][ieta]   =Form("h2DCleanPrKineCutCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f" ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanElKineCutCorr[icent][ieta]   =Form("h2DCleanElKineCutCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f" ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      h2DCleanPiKineCutCorr[icent][ieta]  = new TH2F(hName2DCleanPiKineCutCorr[icent][ieta] ,hName2DCleanPiKineCutCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanPrKineCutCorr[icent][ieta]  = new TH2F(hName2DCleanPrKineCutCorr[icent][ieta] ,hName2DCleanPrKineCutCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanElKineCutCorr[icent][ieta]  = new TH2F(hName2DCleanElKineCutCorr[icent][ieta] ,hName2DCleanElKineCutCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DCleanKaTOFTRDCorr[icent][ieta]=NULL;  h2DCleanKaBayesCorr[icent][ieta]=NULL;  h2DCleanPiTightCorr[icent][ieta]=NULL;  h2DCleanPiTOFCorr[icent][ieta]=NULL;
      hName2DCleanKaTOFTRDCorr[icent][ieta] = Form("h2DCleanKaTOFTRDCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanKaBayesCorr[icent][ieta] = Form("h2DCleanKaBayesCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanPiTightCorr[icent][ieta] = Form("h2DCleanPiTightCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      hName2DCleanPiTOFCorr[icent][ieta] = Form("h2DCleanPiTOFCorr_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f"   ,centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
      h2DCleanKaTOFTRDCorr[icent][ieta]     = new TH2F(hName2DCleanKaTOFTRDCorr[icent][ieta],hName2DCleanKaTOFTRDCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanKaBayesCorr[icent][ieta]     = new TH2F(hName2DCleanKaBayesCorr[icent][ieta],hName2DCleanKaBayesCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanPiTightCorr[icent][ieta]     = new TH2F(hName2DCleanPiTightCorr[icent][ieta],hName2DCleanPiTightCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanPiTOFCorr[icent][ieta]     = new TH2F(hName2DCleanPiTOFCorr[icent][ieta],hName2DCleanPiTOFCorr[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      for (Int_t ipart=0; ipart<nParticles; ipart++)
      {
        h2DCleanCorr[ipart][icent][ieta]=NULL;
        hName2DCleanCorr[ipart][icent][ieta] = Form("h2DCleanCorr%s_cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f", parName[ipart].Data(), centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
        h2DCleanCorr[ipart][icent][ieta]     = new TH2F(hName2DCleanCorr[ipart][icent][ieta],hName2DCleanCorr[ipart][icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      }


    }
  }
  //
  std::cout <<  "testEntries =        " << testEntries        << std::endl;
  std::cout <<  "ptNbins =            " << ptNbins            << std::endl;
  std::cout <<  "nEtaBins =           " << nEtaBins           << std::endl;
  // Create the histograms to be used in the binning of eta, cent and momentum
  hEta  =  new TH1D("hEta" ,"Eta Bins"       ,nEtaBins   ,etaMin, etaMax );
  hMom  =  new TH1D("hMom" ,"Mom Bins"       ,ptNbins    ,ptotMin,  ptotMax );
  hCent =  new TH1D("hCent","Centrality Bins",nCentBins  ,centBinning );
  hEta  -> FillRandom("gaus",1000);
  hCent -> FillRandom("gaus",1000);
  hMom  -> FillRandom("gaus",1000);

}
//____________________________________________________________________________________________________________
Bool_t ApplyTreeSelection(Int_t syst, UInt_t cutBit)
{

  /*
  kNCrossedRowsTPC60=0,
  kNCrossedRowsTPC80=1,
  kNCrossedRowsTPC100=2,
  kMaxChi2PerClusterTPC3=3,
  kMaxChi2PerClusterTPC4=4,
  kMaxChi2PerClusterTPC5=5,
  kMaxDCAToVertexXYPtDepSmall=6,
  kMaxDCAToVertexXYPtDep=7,
  kMaxDCAToVertexXYPtDepLarge=8,
  kVertexZSmall=9,
  kVertexZ=10,
  kVertexZLarge=11,
  kEventVertexZSmall=12,
  kEventVertexZ=13,
  kEventVertexZLarge=14,
  kRequireITSRefit=15,
  kPixelRequirementITS=16,
  kNewITSCut=17,
  kActiveZoneSmall=18,
  kActiveZone=19,
  kActiveZoneLarge=20,
  kTPCSignalNSmall=21,
  kTPCSignalN=22,
  kTPCSignalNLarge=23,
  kTrackProbPiTPC=24,
  kTrackProbKaTPC=25,
  kTrackProbPrTPC=26,
  kTrackProbDeTPC=27,
  kTrackProbPiTOF=28,
  kTrackProbKaTOF=29,
  kTrackProbPrTOF=30,
  kTrackProbDeTOF=31,
  */

  /*
  syst:
  0 -->  Reference
  1 -->  CRows60
  2 -->  CRows100
  3 -->  Chi2TPC3
  4 -->  Chi2TPC5
  5 -->  DCAXYSmall
  6 -->  DCAXYLarge
  7 -->  VZSmall
  8 -->  VZLarge
  9 -->  EventVertexZSmall
  10 --> EventVertexZLarge
  11 --> RequireITSRefit
  11 --> PixelRequirementITS
  12 --> NewITSCut
  // extra settings
  13 --> ActiveZoneSmall,
  14 --> ActiveZone,
  15 --> ActiveZoneLarge,
  16 --> TPCSignalNSmall,
  17 --> TPCSignalN,
  18 --> kTPCSignalNLarge,
  */

  Bool_t cutAll = kFALSE;
  const Int_t fnCutBins=10;
  Int_t fCutArrTmp[fnCutBins]={0};

  switch(syst) {

    case 0:   // 0 -->  Reference
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 1:  // 1 -->  CRows60
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC60,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 2:  // 2 -->  CRows100
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC100, kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 3:   // 3 -->  Chi2TPC3
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC3, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 4:   // 4 -->  Chi2TPC5
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC5, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 5:   // 5 -->  kMaxDCAToVertexXYPtDepSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDepSmall, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 6:   // 6 -->  kMaxDCAToVertexXYPtDepLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDepLarge, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 7:   // 7 -->  kVertexZSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZSmall, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 8:   // 8 -->  kVertexZLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZLarge, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 9:  // 9 -->  kEventVertexZSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZSmall, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 10:  // 10 -->  kEventVertexZLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZLarge, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 11:  // 11 -->  kRequireITSRefit
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, 1,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 12:  // 12 -->  kNewITSCut
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, 1, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 13:  // 13 -->  kPixelRequirementITS
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, 1,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 14:  // 14 -->  kTPCSignalN
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalN,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 15:  // 15 -->  kActiveZone
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,kActiveZone};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 16:  // 16 -->  kTPCSignalN + kActiveZone
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalN,kActiveZone};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 17:  // 17 -->  kTPCSignalNSmall + kActiveZoneSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalNSmall,kActiveZoneSmall};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case 18:  // 18 -->  kTPCSignalNLarge + kActiveZoneLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC4, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalNLarge,kActiveZoneLarge};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    default:
    {
      Int_t fCutArrTmp[fnCutBins] = {0};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }

  }

  //
  //
  //  Apply conditions
  cutAll = fCutArr[0]&&fCutArr[1]&&fCutArr[2]&&fCutArr[3]&&fCutArr[4]&&fCutArr[5]&&fCutArr[6]&&fCutArr[7]&&fCutArr[8]&&fCutArr[9];
  return cutAll;



}
//____________________________________________________________________________________________________________
TGraphErrors * ProfileToGraphErrors(TH2F * h2)
{

  // cout << h2->GetName() << endl;
  TProfile *hpro = (TProfile*)h2->ProfileX();
  const Int_t Nbins = hpro->GetNbinsX();
  // cout << "number of bins = " << Nbins << endl;
  Double_t x[Nbins];
  Double_t y[Nbins];
  Double_t errx[Nbins];
  Double_t erry[Nbins];

  for (Int_t i=0; i<Nbins; i++){
    y[i]=hpro->GetBinContent(i);
    x[i]=hpro->GetBinCenter(i);
    erry[i]=hpro->GetBinError(i);
    errx[i]=0;
  }

  TGraphErrors *gr = new TGraphErrors(Nbins,x,y,errx,erry);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.5);

  for (Int_t i=0; i<gr->GetN();i++){
    if (gr->GetY()[i]<1) {
      gr->RemovePoint(i);
    }
  }
  for (Int_t i=0; i<gr->GetN();i++){
    if (gr->GetY()[i]<1) {
      gr->RemovePoint(i);
    }
  }

  return gr;
}
//____________________________________________________________________________________________________________
void WriteHistsToFile()
{

  histStream->GetFile()->cd();
  for (Int_t icent=0; icent<nCentBins; icent++){
    for (Int_t ieta=0; ieta<nEtaBins; ieta++){

      if(h2Dall[icent][ieta])  h2Dall[icent][ieta]->Write();
      if(h2Dpos[icent][ieta])  h2Dpos[icent][ieta]->Write();
      if(h2Dneg[icent][ieta])  h2Dneg[icent][ieta]->Write();
      //
      if(h2DallPrTOF[icent][ieta])     h2DallPrTOF[icent][ieta]->Write();
      if(h2DallPrTOFPos[icent][ieta])  h2DallPrTOFPos[icent][ieta]->Write();
      if(h2DallPrTOFNeg[icent][ieta])  h2DallPrTOFNeg[icent][ieta]->Write();
      if(h2DallKaTOF[icent][ieta])     h2DallKaTOF[icent][ieta]->Write();
      if(h2DallKaTOFPos[icent][ieta])  h2DallKaTOFPos[icent][ieta]->Write();
      if(h2DallKaTOFNeg[icent][ieta])  h2DallKaTOFNeg[icent][ieta]->Write();
      //
      if(h2DCleanPiKineCut[icent][ieta]) h2DCleanPiKineCut[icent][ieta]->Write();
      if(h2DCleanPrKineCut[icent][ieta]) h2DCleanPrKineCut[icent][ieta]->Write();
      if(h2DCleanElKineCut[icent][ieta]) h2DCleanElKineCut[icent][ieta]->Write();
      //
      if(h2DCleanKaTOFTRD[icent][ieta]) h2DCleanKaTOFTRD[icent][ieta]->Write();
      if(h2DCleanKaBayes[icent][ieta])  h2DCleanKaBayes[icent][ieta]->Write();
      if(h2DCleanPiTight[icent][ieta])  h2DCleanPiTight[icent][ieta]->Write();
      if(h2DCleanPiTOF[icent][ieta])    h2DCleanPiTOF[icent][ieta]->Write();
      //
      for (Int_t ipart=0; ipart<nParticles; ipart++){
        if(h2Expected[ipart][icent][ieta])       h2Expected[ipart][icent][ieta]->Write();
        if(h2ExpectedSigma[ipart][icent][ieta])  h2ExpectedSigma[ipart][icent][ieta]->Write();
        if(grExpected[ipart][icent][ieta])       grExpected[ipart][icent][ieta]->Write();
        if(grExpectedSigma[ipart][icent][ieta])  grExpectedSigma[ipart][icent][ieta]->Write();
        if(h2DClean[ipart][icent][ieta])         h2DClean[ipart][icent][ieta]->Write();
      }
    }
  }


  histStreamCorr->GetFile()->cd();
  for (Int_t icent=0; icent<nCentBins; icent++){
    for (Int_t ieta=0; ieta<nEtaBins; ieta++){
      //
      if(h2DallCorr[icent][ieta])  h2DallCorr[icent][ieta]->Write();
      if(h2DposCorr[icent][ieta])  h2DposCorr[icent][ieta]->Write();
      if(h2DnegCorr[icent][ieta])  h2DnegCorr[icent][ieta]->Write();
      //
      if(h2DallPrTOFCorr[icent][ieta])     h2DallPrTOFCorr[icent][ieta]->Write();
      if(h2DallPrTOFPosCorr[icent][ieta])  h2DallPrTOFPosCorr[icent][ieta]->Write();
      if(h2DallPrTOFNegCorr[icent][ieta])  h2DallPrTOFNegCorr[icent][ieta]->Write();
      if(h2DallKaTOFCorr[icent][ieta])     h2DallKaTOFCorr[icent][ieta]->Write();
      if(h2DallKaTOFPosCorr[icent][ieta])  h2DallKaTOFPosCorr[icent][ieta]->Write();
      if(h2DallKaTOFNegCorr[icent][ieta])  h2DallKaTOFNegCorr[icent][ieta]->Write();
      //
      if(h2DCleanPiKineCutCorr[icent][ieta]) h2DCleanPiKineCutCorr[icent][ieta]->Write();
      if(h2DCleanPrKineCutCorr[icent][ieta]) h2DCleanPrKineCutCorr[icent][ieta]->Write();
      if(h2DCleanElKineCutCorr[icent][ieta]) h2DCleanElKineCutCorr[icent][ieta]->Write();
      //
      if(h2DCleanKaTOFTRDCorr[icent][ieta]) h2DCleanKaTOFTRDCorr[icent][ieta]->Write();
      if(h2DCleanKaBayesCorr[icent][ieta])  h2DCleanKaBayesCorr[icent][ieta]->Write();
      if(h2DCleanPiTightCorr[icent][ieta])  h2DCleanPiTightCorr[icent][ieta]->Write();
      if(h2DCleanPiTOFCorr[icent][ieta])    h2DCleanPiTOFCorr[icent][ieta]->Write();
      //
      for (Int_t ipart=0; ipart<nParticles; ipart++){
        if(h2DCleanCorr[ipart][icent][ieta]) h2DCleanCorr[ipart][icent][ieta]->Write();
      }
    }
  }


  //
  debugStream->GetFile()->cd();
  if(htracks_dcaxy2D)               htracks_dcaxy2D  ->Write();
  if(htracks_dcaz2D)                htracks_dcaz2D   ->Write();
  if(hevents_pileUpV02D)            hevents_pileUpV02D  ->Write();
  if(hevents_pileUpITS2D)           hevents_pileUpITS2D  ->Write();
  if(htracks_shiftM_nPileUpPrim)    htracks_shiftM_nPileUpPrim  ->Write();
  if(hevents_pileUpV01D)            hevents_pileUpV01D  ->Write();
  if(hevents_pileUpITS1D)           hevents_pileUpITS1D  ->Write();
  if(hevents_ITSTPCeff)             hevents_ITSTPCeff  ->Write();
  if(htracks_nclTPC2D)              htracks_nclTPC2D ->Write();
  if(hdscaled_sharedTPCClusters1D)  hdscaled_sharedTPCClusters1D ->Write();
  if(hdscaled_tpcSignalN1D)         hdscaled_tpcSignalN1D ->Write();
  if(hdscaled_lengthInActiveZone1D) hdscaled_lengthInActiveZone1D ->Write();
  if(hdscaled_cRows1D)              hdscaled_cRows1D ->Write();
  if(hdscaled_nclITS1D)             hdscaled_nclITS1D ->Write();
  if(hdscaled_nclTRD1D)             hdscaled_nclTRD1D ->Write();
  if(hdscaled_chi2TPC1D)            hdscaled_chi2TPC1D ->Write();
  if(hdscaled_chi2ITS1D)            hdscaled_chi2ITS1D ->Write();
  if(hdscaled_chi2TRD1D)            hdscaled_chi2TRD1D ->Write();
  if(hdscaled_vz1D)                 hdscaled_vz1D      ->Write();
  if(hdscaled_dcaxy1D)              hdscaled_dcaxy1D  ->Write();
  if(hdscaled_dcaz1D)               hdscaled_dcaz1D   ->Write();
  if(hdscaled_nclTPC1D)             hdscaled_nclTPC1D ->Write();
  if(hhighPt_dEdxPtot)              hhighPt_dEdxPtot->Write();

}
//____________________________________________________________________________________________________________
void SetBranchAddresses()
{

  //
  // Clean Samples Tree
  if (armtree){
    armtree->SetBranchAddress("gid"        ,&ffArmPodTree_gid);
    armtree->SetBranchAddress("eventtime"  ,&ffArmPodTree_eventtime);
    armtree->SetBranchAddress("intrate"    ,&ffArmPodTree_intrate);
    armtree->SetBranchAddress("cutBit"     ,&ffArmPodTree_cutBit);
    armtree->SetBranchAddress("purity"     ,&ffArmPodTree_purity);
    armtree->SetBranchAddress("dEdx"       ,&ffArmPodTree_dEdx);
    armtree->SetBranchAddress("sign"       ,&ffArmPodTree_sign);
    armtree->SetBranchAddress("ptot"       ,&ffArmPodTree_ptot);
    armtree->SetBranchAddress("p"          ,&ffArmPodTree_p);
    armtree->SetBranchAddress("pT"         ,&ffArmPodTree_pT);
    armtree->SetBranchAddress("eta"        ,&ffArmPodTree_eta);
    armtree->SetBranchAddress("cent"       ,&ffArmPodTree_cent);
    armtree->SetBranchAddress("qt"         ,&ffArmPodTree_qt);
    armtree->SetBranchAddress("alfa"       ,&ffArmPodTree_alfa);
    armtree->SetBranchAddress("piTOFnSigma",&ffArmPodTree_piTOFnSigma);
    armtree->SetBranchAddress("prTOFnSigma",&ffArmPodTree_prTOFnSigma);
    armtree->SetBranchAddress("piFromK0"   ,&ffArmPodTree_piFromK0);
    armtree->SetBranchAddress("v0haspixel" ,&ffArmPodTree_v0haspixel);
  }
  //
  // Tracks Tree
  if (dataTree){
    dataTree->SetBranchAddress("gid"      ,&ftracks_gid);
    dataTree->SetBranchAddress("eventtime",&ftracks_eventtime);
    dataTree->SetBranchAddress("intrate"  ,&ftracks_intrate);
    dataTree->SetBranchAddress("cutBit"   ,&ftracks_cutBit);
    dataTree->SetBranchAddress("dEdx"     ,&ftracks_dEdx);
    dataTree->SetBranchAddress("sign"     ,&ftracks_sign);
    dataTree->SetBranchAddress("ptot"     ,&ftracks_ptot);
    dataTree->SetBranchAddress("p"        ,&ftracks_p);
    dataTree->SetBranchAddress("pT"       ,&ftracks_pT);
    dataTree->SetBranchAddress("eta"      ,&ftracks_eta);
    dataTree->SetBranchAddress("cent"     ,&ftracks_cent);
    dataTree->SetBranchAddress("nclTPC"   ,&ftracks_nclTPC);
    dataTree->SetBranchAddress("dcaxy"    ,&ftracks_dcaxy);
    dataTree->SetBranchAddress("dcaz"     ,&ftracks_dcaz);
  }
  //
  // dscaled tree
  if (dscaltree){
    dscaltree->SetBranchAddress("dcaxy"       ,&fdscaled_dcaxy);
    dscaltree->SetBranchAddress("dcaz"        ,&fdscaled_dcaz);
    dscaltree->SetBranchAddress("cutBit"      ,&fdscaled_cutBit);
    dscaltree->SetBranchAddress("cent"        ,&fdscaled_cent);
    dscaltree->SetBranchAddress("gid"         ,&fdscaled_gid);
    dscaltree->SetBranchAddress("intrate"     ,&fdscaled_intrate);
    dscaltree->SetBranchAddress("eventtime"   ,&fdscaled_eventtime);
    dscaltree->SetBranchAddress("vz"          ,&fdscaled_vz);
    dscaltree->SetBranchAddress("sharedTPCClusters" ,&fdscaled_sharedTPCClusters);
    dscaltree->SetBranchAddress("nclTPC"      ,&fdscaled_nclTPC);
    dscaltree->SetBranchAddress("tpcSignalN"  ,&fdscaled_tpcSignalN);
    dscaltree->SetBranchAddress("cRows"       ,&fdscaled_cRows);
    dscaltree->SetBranchAddress("lengthInActiveZone" ,&fdscaled_lengthInActiveZone);
    dscaltree->SetBranchAddress("phi"         ,&fdscaled_phi);
    dscaltree->SetBranchAddress("phiTPC"      ,&fdscaled_phiTPC);
    dscaltree->SetBranchAddress("phi85"       ,&fdscaled_phi85);
    dscaltree->SetBranchAddress("itsdEdx"     ,&fdscaled_itsdEdx);
    dscaltree->SetBranchAddress("trddEdx"     ,&fdscaled_trddEdx);
    dscaltree->SetBranchAddress("fdTPC"       ,&fdscaled_fdTPC);
    dscaltree->SetBranchAddress("fzTPC"       ,&fdscaled_fzTPC);
    dscaltree->SetBranchAddress("fCdd"        ,&fdscaled_fCdd);
    dscaltree->SetBranchAddress("fCdz"        ,&fdscaled_fCdz);
    dscaltree->SetBranchAddress("fCzz"        ,&fdscaled_fCzz);
    dscaltree->SetBranchAddress("nclITS"      ,&fdscaled_nclITS);
    dscaltree->SetBranchAddress("nclTRD"      ,&fdscaled_nclTRD);
    dscaltree->SetBranchAddress("chi2TPC"     ,&fdscaled_chi2TPC);
    dscaltree->SetBranchAddress("chi2ITS"     ,&fdscaled_chi2ITS);
    dscaltree->SetBranchAddress("chi2TRD"     ,&fdscaled_chi2TRD);
  }
  //
  // event tree
  if (eventtree){
    eventtree->SetBranchAddress("run"                ,&fevents_run);
    eventtree->SetBranchAddress("intrate"            ,&fevents_intrate);
    eventtree->SetBranchAddress("bField"             ,&fevents_bField);
    eventtree->SetBranchAddress("gid"                ,&fevents_gid);
    eventtree->SetBranchAddress("timestamp"         ,&fevents_timestamp);
    eventtree->SetBranchAddress("triggerMask"        ,&fevents_triggerMask);
    eventtree->SetBranchAddress("vz"                 ,&fevents_vz);
    eventtree->SetBranchAddress("tpcvz"              ,&fevents_tpcvz);
    eventtree->SetBranchAddress("spdvz"              ,&fevents_spdvz);
    eventtree->SetBranchAddress("tpcMult"            ,&fevents_tpcMult);
    eventtree->SetBranchAddress("eventMult"          ,&fevents_eventMult);
    eventtree->SetBranchAddress("eventMultESD"       ,&fevents_eventMultESD);
    eventtree->SetBranchAddress("nTracksStored"      ,&fevents_nTracksStored);
    eventtree->SetBranchAddress("primMult"           ,&fevents_primMult);
    eventtree->SetBranchAddress("tpcTrackBeforeClean",&fevents_tpcTrackBeforeClean);
    eventtree->SetBranchAddress("itsTracklets"       ,&fevents_itsTracklets);
    eventtree->SetBranchAddress("centrality."         ,&fevents_centrality);
    eventtree->SetBranchAddress("tZeroMult."          ,&fevents_tZeroMult);
    eventtree->SetBranchAddress("vZeroMult."          ,&fevents_vZeroMult);
    eventtree->SetBranchAddress("itsClustersPerLayer.",&fevents_itsClustersPerLayer);
    eventtree->SetBranchAddress("trackCounters."      ,&fevents_trackCounters);
    eventtree->SetBranchAddress("trackdEdxRatio."     ,&fevents_trackdEdxRatio);
    eventtree->SetBranchAddress("trackNcl."           ,&fevents_trackNcl);
    eventtree->SetBranchAddress("trackChi2."          ,&fevents_trackChi2);
    eventtree->SetBranchAddress("trackMatchEff."      ,&fevents_trackMatchEff);
    eventtree->SetBranchAddress("trackTPCCountersZ."  ,&fevents_trackTPCCountersZ);
    eventtree->SetBranchAddress("phiCountA."          ,&fevents_phiCountA);
    eventtree->SetBranchAddress("phiCountC."          ,&fevents_phiCountC);
    eventtree->SetBranchAddress("phiCountAITS."       ,&fevents_phiCountAITS);
    eventtree->SetBranchAddress("phiCountCITS."       ,&fevents_phiCountCITS);
    eventtree->SetBranchAddress("phiCountAITSOnly."   ,&fevents_phiCountAITSOnly);
    eventtree->SetBranchAddress("phiCountCITSOnly."   ,&fevents_phiCountCITSOnly);
    eventtree->SetBranchAddress("tpcVertexInfo."      ,&fevents_tpcVertexInfo);
    // eventtree->SetBranchAddress("itsVertexInfo."      ,&fevents_itsVertexInfo);

    eventtree->SetAlias("shiftM"       ,"0.5*(tpcVertexInfo.fElements[1]+tpcVertexInfo.fElements[0])");
    eventtree->SetAlias("nPileUpPrim"  ,"(tpcVertexInfo.fElements[3]+tpcVertexInfo.fElements[4])");
  }

  //
  // event tree
  if (highPttree){
    highPttree->SetBranchAddress("gid"             ,&fhighPt_gid);
    highPttree->SetBranchAddress("selectionPtMask" ,&fhighPt_selectionPtMask);
    highPttree->SetBranchAddress("centralityF"     ,&fhighPt_centralityF);
    highPttree->SetBranchAddress("runNumber"       ,&fhighPt_runNumber);
    highPttree->SetBranchAddress("evtTimeStamp"    ,&fhighPt_evtTimeStamp);
    highPttree->SetBranchAddress("evtNumberInFile" ,&fhighPt_evtNumberInFile);
    highPttree->SetBranchAddress("triggerClass"    ,&fhighPt_triggerClass);
    highPttree->SetBranchAddress("Bz"              ,&fhighPt_Bz);
    highPttree->SetBranchAddress("IRtot"           ,&fhighPt_IRtot);
    highPttree->SetBranchAddress("IRint2"          ,&fhighPt_IRint2);
    highPttree->SetBranchAddress("mult"            ,&fhighPt_mult);
    highPttree->SetBranchAddress("esdTrack."       ,&fhighPt_esdTrack);
    highPttree->SetBranchAddress("fileName."       ,&fhighPt_fileName);
    highPttree->SetBranchAddress("vtxESD."         ,&fhighPt_vtxESD);
  }
  //
  // "downscaleCounter="<<downscaleCounter<<
  // "fLowPtTrackDownscaligF="<<fLowPtTrackDownscaligF<<
  // "selectionPtMask="<<selectionPtMask<<          // high pt trigger mask
  // "selectionPIDMask="<<selectionPIDMask<<         // selection PIDmask
  //
  // Initialize Map
  fdEdxMap   = TFile::Open(dEdxMapStr);
  hdEdxAShifttMNTglDist_meanGFitAll      = (AliNDLocalRegression*)fdEdxMap->Get("hdEdxAShifttMNTglDist_meanGFitA");
  hdEdxAShifttMNTglDist_meanGFitNoPileUp = (AliNDLocalRegression*)fdEdxMap->Get("hdEdxAShifttMNTglDist_meanGFit0");
  //
  // This is needed in case you want to use it in TF1 or TTreeFormula
  AliNDLocalRegression::AddVisualCorrection(hdEdxAShifttMNTglDist_meanGFitAll,2); //register correction - 4 parameters
  AliNDLocalRegression::AddVisualCorrection(hdEdxAShifttMNTglDist_meanGFitNoPileUp,1); //register correction - 2 parameters

}
//____________________________________________________________________________________________________________
void PlotTimeSeriesPerSector(TString performanceEventList)
{


  /*


  .L /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_FilterTreesMakeHists.C+
  TString file15o0="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_15o_246272.list"
  TString file15o="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_15o.list"
  TString file18q="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_18q.list"
  TString file18r="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_18r.list"
  PlotTimeSeriesPerSector(file15o0);


  */

  ::Info("LoadSumaryTrees","Load simple time trees");
  TTree *treeTimeAll=AliTreePlayer::LoadTrees(Form("cat %s",performanceEventList.Data()),"(hisTime.*All|hisTimeChi2.*|hisTimeRa.*|hisTimeNcl.*|hisTimeEff.*|hisPrimtMult.*|hisTPCMult.*)","xxx",".*","","");
  ::Info("LoadSumaryTrees","Load time trees with phi granularity");
  TTree *treePhi=AliTreePlayer::LoadTrees(Form("cat %s",performanceEventList.Data()),"hisTimephi","AllDist",".*","","");
  //
  treePhi    ->SetAlias("sectorBin","Iteration$Bin/2.");
  treePhi    ->SetAlias("T","timestampCenter");
  treeTimeAll->SetAlias("T","timestampCenter");
  treePhi    ->SetAlias("ntrITSRatio","hisTimephiCountAITSOnly0Dist.entries/(hisTimephiCountAITS0Dist.entries+hisTimephiCountAITSOnly0Dist.entries)");

  TStatToolkit::AddMetadata(treePhi,"ntrITSRatio.AxisTitle","N_{tr}(ITS_{only})/(N_{tr}(ITS_{only})+N_{tr}(ITS+TPC))");
  treePhi->SetAlias("ntrITSRatioC","hisTimephiCountCITSOnly0Dist.entries/(hisTimephiCountCITS0Dist.entries+hisTimephiCountCITSOnly0Dist.entries)");
  TStatToolkit::AddMetadata(treePhi,"ntrITSRatioC.AxisTitle","N_{trC}(ITSC_{only})/(N_{tr}(ITSC_{only})+N_{trC}(ITS+TPC))");
  treePhi    ->BuildIndex("timestampBin");
  treeTimeAll->BuildIndex("timestampBin");
  treePhi->AddFriend(treeTimeAll,"TA");

  TCanvas *canvas = new TCanvas("canvas","canvas",1200,700);
  canvas->Divide(1,2);
  canvas->cd(1);
  treePhi->Draw("ntrITSRatio:T:sectorBin","abs(sectorBin-2)<2","colz");
  canvas->cd(2);
  treeTimeAll->Draw("hisTimeEffITSDist.binMedian:T","hisTimeEffITSDist.binMedian>0.5");


}
//____________________________________________________________________________________________________________
void PlotTimeSeriesFullPeriod(TString performanceEventList)
{

  /*


  .L /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_FilterTreesMakeHists.C+
  TString file15o0="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_15o_246272.list"
  TString file15o="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_15o.list"
  TString file18q="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_18q.list"
  TString file18r="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_18r.list"
  PlotTimeSeriesFullPeriod(file15o);


  */

  TString file15o="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_15o.list";
  TString file18q="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_18q.list";
  TString file18r="/lustre/nyx/alice/users/marsland/Maps/dEdxCorrection_PileUp/timeSeries_18r.list";

  TChain *chain15o = AliXRDPROOFtoolkit::MakeChainRandom(file15o,"hisTimeEffITSDist",0,1000);
  Double_t entries0  =  chain15o->Draw("binMedian>>hisEff(300,0.5,1)","binMedian>0.5","");

  TChain *chain18q = AliXRDPROOFtoolkit::MakeChainRandom(file18q,"hisTimeEffITSDist",0,1000);
  Double_t entries1  =  chain18q->Draw("binMedian>>hisEff(300,0.5,1)","binMedian>0.5","");

  TChain *chain18r = AliXRDPROOFtoolkit::MakeChainRandom(file18r,"hisTimeEffITSDist",0,1000);
  Double_t entries2  =  chain18r->Draw("binMedian>>hisEff(300,0.5,1)","binMedian>0.5","");

  TString dEdxMapStr = "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/dEdxMap15o/dEdxFit.root";
  TFile *dEdxMap = TFile::Open(dEdxMapStr);
  AliNDLocalRegression *hdEdxAShifttMNTglDist_meanGFitAll      = hdEdxAShifttMNTglDist_meanGFitAll      = (AliNDLocalRegression*)dEdxMap->Get("hdEdxAShifttMNTglDist_meanGFitA");
  AliNDLocalRegression *hdEdxAShifttMNTglDist_meanGFitNoPileUp = hdEdxAShifttMNTglDist_meanGFitNoPileUp = (AliNDLocalRegression*)dEdxMap->Get("hdEdxAShifttMNTglDist_meanGFit0");
  AliNDLocalRegression::AddVisualCorrection(hdEdxAShifttMNTglDist_meanGFitAll,2); //register correction - 4 parameters
  AliNDLocalRegression::AddVisualCorrection(hdEdxAShifttMNTglDist_meanGFitNoPileUp,1); //register correction - 2 parameters
  //
  TF2* f2DdEdxCorrMap     = new TF2("f2DdEdxCorrMap"    ,"AliNDLocalRegression::GetCorrND(2,x,y,0,0)"   ,-100,100,0,6000); // x: pileupvZ, y: pileUpMult, z:primary mult: t:pZ/pt
  TF1* f1DdEdxCorrMap500  = new TF1("f1DdEdxCorrMap500" ,"AliNDLocalRegression::GetCorrND(2,x,500,0,0)",-100,100);
  TF1* f1DdEdxCorrMap1000 = new TF1("f1DdEdxCorrMap1000","AliNDLocalRegression::GetCorrND(2,x,1000,0,0)",-100,100);
  TF1* f1DdEdxCorrMap2000 = new TF1("f1DdEdxCorrMap2000","AliNDLocalRegression::GetCorrND(2,x,2000,0,0)",-100,100);
  TF1* f1DdEdxCorrMap4000 = new TF1("f1DdEdxCorrMap4000","AliNDLocalRegression::GetCorrND(2,x,4000,0,0)",-100,100);
  TF1* f1DdEdxCorrMap6000 = new TF1("f1DdEdxCorrMap6000","AliNDLocalRegression::GetCorrND(2,x,6000,0,0)",-100,100);
  const Int_t color[] = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
  f1DdEdxCorrMap500  -> SetLineColor(color[0]);
  f1DdEdxCorrMap1000 -> SetLineColor(color[1]);
  f1DdEdxCorrMap2000 -> SetLineColor(color[2]);
  f1DdEdxCorrMap4000 -> SetLineColor(color[3]);
  f1DdEdxCorrMap6000 -> SetLineColor(color[4]);
  //
  f1DdEdxCorrMap500  -> SetLineWidth(2);
  f1DdEdxCorrMap1000 -> SetLineWidth(2);
  f1DdEdxCorrMap2000 -> SetLineWidth(2);
  f1DdEdxCorrMap4000 -> SetLineWidth(2);
  f1DdEdxCorrMap6000 -> SetLineWidth(2);

  TF2 *f2D = (TF2*)f2DdEdxCorrMap->Clone();
  TF1 *f0 = (TF1*)f1DdEdxCorrMap500->Clone();
  TF1 *f1 = (TF1*)f1DdEdxCorrMap1000->Clone();
  TF1 *f2 = (TF1*)f1DdEdxCorrMap2000->Clone();
  TF1 *f3 = (TF1*)f1DdEdxCorrMap4000->Clone();
  TF1 *f4 = (TF1*)f1DdEdxCorrMap6000->Clone();

  //
  // TCanvas *canvas = new TCanvas("canvas","canvas",1200,700);
  // canvas->cd();
  // f2DdEdxCorrMap->Draw("colz");
  // TCanvas *canvas1 = new TCanvas("canvas1","canvas1",1200,700);
  // canvas1->cd();
  // f1DdEdxCorrMap6000->Draw();
  // f1DdEdxCorrMap500->Draw("same");
  // f1DdEdxCorrMap1000->Draw("same");
  // f1DdEdxCorrMap2000->Draw("same");
  // f1DdEdxCorrMap4000->Draw("same");
  //

  TTreeSRedirector *maps12D = new TTreeSRedirector("maps12D.root","recreate");
  maps12D->GetFile()->cd();
  f2D->Write();
  f0->Write();
  f1->Write();
  f2->Write();
  f3->Write();
  f4->Write();
  delete maps12D;

}


/*

cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/Fits/testHists

// TH2F *h1 = (TH2F*)f->Get("h2DCleanPiTOF")
// h->GetYaxis()->SetRangeUser(35,70); h1->GetYaxis()->SetRangeUser(35,70);

TFile *f = new TFile("Hists_PbPb_syst0_eta_0.00_0.20_cent_0.00_10.00.root");
TH2F *h  = (TH2F*)f->Get("h2DCleanPiTight")
h->RebinX(4)
TF1 *func = new TF1("func","gaus",45,55);
TObjArray arr;
h->Draw("colz");
h->FitSlicesY(func, 0, -1,0, "QNR", &arr);
TH1D *hh = (TH1D*)arr[1];
hh->SetLineColor(kBlack)

TH2F *h1  = (TH2F*)f->Get("h2DCleanPiTOF")
h1->RebinX(4)
TObjArray arr1;
h1->Draw("colz");
h1->FitSlicesY(func, 0, -1,0, "QNR", &arr1);
TH1D *hh1 = (TH1D*)arr1[1];
hh->SetLineColor(kGreen)

hh->Draw()
hh1->Draw("same")

*/

/*

TFile f("AnalysisResults_hist.root")
TList * list  = (TList*)f.Get("cleanHists");
fhnExpected = (THnSparse*)list->FindObject(Form("hExpected_%d",0));
fhnExpected->GetAxis(0)->SetRangeUser(1,2);       // particle type
fhnExpected->GetAxis(1)->SetRangeUser(-1,0);      // sign
fhnExpected->GetAxis(2)->SetRangeUser(0,5);       // centrality
fhnExpected->GetAxis(3)->SetRangeUser(0,0.2);     // eta
fhnExpected->Projection(6,4)->Draw()


TFile f("AnalysisResults_tree1.root");
TVectorF cent(27);
TVectorF *acent=&cent;
events->SetBranchAddress("centrality.",&acent);
events->GetEntry(2);
cout << (*acent)[2]  << endl;

TFile f("AnalysisResults_tree1.root");
TVectorF *cent=0x0;
events->SetBranchAddress("centrality.",&cent);
events->GetEntry(2);
cout << (*cent)[2]  << endl;



*/
