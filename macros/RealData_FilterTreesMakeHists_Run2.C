/*

meld /home/marsland/Desktop/ubuntu_desktop/workdir/code/RealData_FilterTreesMakeHists.C /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_FilterTreesMakeHists_Run2.C

events->Draw("pileUp1DITS/multTPC","multTPC-pileUp1DITS<1000 && pileUp1DITS/multTPC<1.4 && pileUp1DITS/multTPC>0.005")
events->Draw("multV0/multTPC/2.55","multTPC-pileUp1DITS<1000 && pileUp1DITS/multTPC<1.4 && pileUp1DITS/multTPC>0.005","same")

clean->Draw("dEdx:ptot","abs(dEdx-50)<50&&abs(alfa)<0.5&&purity==2&&qt>0.14","")

events->SetAlias("primMultITS","itsVertexInfo.fElements[2]");
events->SetAlias("secMultITS0","itsVertexInfo.fElements[3]");
events->SetAlias("secMultITS","itsVertexInfo.fElements[3]-0.01*itsVertexInfo.fElements[2]");
events->SetAlias("cut5_50","secMultITS<0.05*primMultITS+50");
events->Draw("secMultITS0:primMultITS>>his(100,0,3000,100,0,100)","","colz",1000000);
events->GetHistogram()->Fit("pol1")

//
events->SetMarkerColor();events->SetMarkerStyle(1);
events->Draw("itsClustersPerLayer.fElements[0]:vZeroMult->Sum()","spdvz!=0&&abs(spdvz)<8","",1000000);
events->SetMarkerColor(2);events->SetMarkerStyle(1);
events->Draw("itsClustersPerLayer.fElements[0]:vZeroMult->Sum()","spdvz!=0&&abs(spdvz)<8&&!cut5_50","same",10000000);



kPileUp = 1000, 2000, 4000, 1000000;
kTimeSeriesEff = 0.80, 0.85, 0;
kVzPileup = +, -, All;

// Marian:
TFile f("events.root")
treeEvent = events
makeAliasesEvent()
EventFlat.SetAlias("multTPCITS","(multSSD+multSDD)/2.38")
EventFlat.Draw("tpcMult-multTPCITS:nPileUpSum","centV0<80&&spdvz!=0","",100000)
EventFlat.Draw("multV0:tpcMult-multTPCITS","centV0<80&&spdvz!=0","",100000)
//
EventFlat.SetAlias("multTPCITS","(multSSD+multSPD)/2.38")
EventFlat.Draw("mutlTPCITS/tpcMult","centV0<20","",1000)
EventFlat.Draw("multTPCITS/tpcMult","centV0<20","",1000)
EventFlat.Draw("multTPCITS/tpcMult","centV0<50&&spdvz!=0","",1000)
EventFlat.Draw("multTPCITS/tpcMult","centV0<50&&spdvz!=0","",10000)
EventFlat.Draw("multTPCITS/tpcMult>>his(100,0,2)","centV0<50&&spdvz!=0","",10000)
EventFlat.Draw("multTPCITS/tpcMult>>his(100,0,2)","centV0<50&&spdvz!=0","",10000)
EventFlat.Draw("multTPCITS/tpcMult>>his(100,0,2)","centV0<50&&spdvz!=0","",100000)
EventFlat.Draw("multTPCITS-tpcMult:centV0","centV0<80&&spdvz!=0","",100000)
EventFlat.Draw("multTPCITS-tpcMult:centV0","multV0<80&&spdvz!=0","",100000)
EventFlat.Draw("multTPCITS-tpcMult:centV0","centV0<80&&spdvz!=0","",100000)
EventFlat.Draw("mutlV0:multTPCITS-tpcMult","centV0<80&&spdvz!=0","",100000)
EventFlat.Draw("multV0:multTPCITS-tpcMult","centV0<80&&spdvz!=0","",100000)
EventFlat.Draw("multV0:tpcMult-multTPCITS","centV0<80&&spdvz!=0","",100000)
EventFlat.Draw("tpcMult-multTPCITS:nPileUpSum","centV0<80&&spdvz!=0","",100000)
EventFlat.Draw("nPileUpSum:tpcMult-multTPCITS","centV0<80&&spdvz!=0","",100000)



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
void ProcessJetConstHists();
void ProcessCleanSamples();
void ProcessDScaledTree();
void ProcessHighPtTree();
void ProcessEventTree();
void ProcessJetEventTree();
void PlotTimeSeriesPerSector(TString performanceEventList);
void CreateSplinesFromTHnSparse(Int_t icent, Int_t ieta);
Bool_t ApplyTreeSelection(Int_t syst, UInt_t cutBit);
TGraphErrors * ProfileToGraphErrors(TH2F * h2);
void SetBranchAddresses();
void ProcessExpectedHists(TString plotVS, TString cutON, Int_t systSet, TString hFile);
Double_t TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc);
Double_t TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc);


// ======= Modification part =======
// inbunch pileup settings:
Double_t inbunchCutSlope  = 0.095;
Double_t inbunchCutOffSet = 40;
Double_t inbunchCutRMS    = 0.05;
Int_t runList[] = {246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928,
  246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847,
  246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804,
  246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750,
  246676, 246675, 246540, 246495, 246493, 246488, 246487, 246434, 246431,
  246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217,
  246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148,
  246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042,
  246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949,
  245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};
//
const Int_t colors[] = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
const Int_t nParticles = 5;
const Int_t nCentBins = 13;
Bool_t   oldMap       = kFALSE;
//
const Double_t dEdxNbins = 1000;   // ??? be careful it must match with the analysis task binwidth
Double_t dEdxMin = 20;
Double_t dEdxMax = 1020;
//
// const Double_t ptNbins = 150;   // mostly for real data analysis
// Double_t ptotMin = 0.1;
// Double_t ptotMax = 3.1;
const Double_t ptNbins = 250;   // mostly for real data analysis
Double_t ptotMin = 0.2;
Double_t ptotMax = 5.2;
//
// const Int_t nEtaBins = 16;  // ???
// Double_t etaMin = -0.8;
// Double_t etaMax = 0.8;
const Int_t nEtaBins = 8;  // ???
const Int_t njetRad = 3;
Double_t etaMin = -0.8;
Double_t etaMax = 0.8;
Double_t jet_etaMin = -0.9;
Double_t jet_etaMax = 0.9;
//
const Int_t nCentDim = 14;
// const Int_t nEtaDim = 17;
// Float_t etaBinning[nEtaDim]   = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
const Int_t nEtaDim = 19;
const Int_t njetRadDim = 3;
Float_t jetRadSizes[njetRadDim]   = {0.2,0.4,0.6};
Float_t etaBinning[nEtaDim]   = {-0.9,-0.8,-0.7,-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
Float_t jet_etaBinning[nEtaDim]   = {-0.9,-0.8,-0.7,-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
Float_t centBinning[nCentDim] = { 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 85., 90., 95., 100.};
//
TTreeSRedirector *tidenTreeStream[nCentBins];
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
TH2F *h2DallBkgTOF[nCentBins][nEtaBins];            TString hName2DallBkgTOF[nCentBins][nEtaBins];
TH2F *h2DallBkgTOFPos[nCentBins][nEtaBins];         TString hName2DallBkgTOFPos[nCentBins][nEtaBins];
TH2F *h2DallBkgTOFNeg[nCentBins][nEtaBins];         TString hName2DallBkgTOFNeg[nCentBins][nEtaBins];
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
// ----------------------------------------------------------------------------------------------------------------------
//
// Jet const hists
TH2F *jet_h2Dall[nCentBins][nEtaBins][njetRad];                 TString jet_hName2Dall[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2Dpos[nCentBins][nEtaBins][njetRad];                 TString jet_hName2Dpos[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2Dneg[nCentBins][nEtaBins][njetRad];                 TString jet_hName2Dneg[nCentBins][nEtaBins][njetRad];
//
TH2F *jet_h2DallPiTOF[nCentBins][nEtaBins][njetRad];            TString jet_hName2DallPiTOF[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2DallPiTOFPos[nCentBins][nEtaBins][njetRad];         TString jet_hName2DallPiTOFPos[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2DallPiTOFNeg[nCentBins][nEtaBins][njetRad];         TString jet_hName2DallPiTOFNeg[nCentBins][nEtaBins][njetRad];
//
TH2F *jet_h2DallPrTOF[nCentBins][nEtaBins][njetRad];            TString jet_hName2DallPrTOF[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2DallPrTOFPos[nCentBins][nEtaBins][njetRad];         TString jet_hName2DallPrTOFPos[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2DallPrTOFNeg[nCentBins][nEtaBins][njetRad];         TString jet_hName2DallPrTOFNeg[nCentBins][nEtaBins][njetRad];
//
TH2F *jet_h2DallBkgTOF[nCentBins][nEtaBins][njetRad];            TString jet_hName2DallBkgTOF[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2DallBkgTOFPos[nCentBins][nEtaBins][njetRad];         TString jet_hName2DallBkgTOFPos[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2DallBkgTOFNeg[nCentBins][nEtaBins][njetRad];         TString jet_hName2DallBkgTOFNeg[nCentBins][nEtaBins][njetRad];
//
TH2F *jet_h2DallKaTOF[nCentBins][nEtaBins][njetRad];            TString jet_hName2DallKaTOF[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2DallKaTOFPos[nCentBins][nEtaBins][njetRad];         TString jet_hName2DallKaTOFPos[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2DallKaTOFNeg[nCentBins][nEtaBins][njetRad];         TString jet_hName2DallKaTOFNeg[nCentBins][nEtaBins][njetRad];
//
TH2F *jet_h2DCleanKaTOFTRD[nCentBins][nEtaBins][njetRad];       TString jet_hName2DCleanKaTOFTRD[nCentBins][nEtaBins][njetRad];
TH2F *jet_h2DCleanKaBayes[nCentBins][nEtaBins][njetRad];        TString jet_hName2DCleanKaBayes[nCentBins][nEtaBins][njetRad];
//
TH2F *jet_h2DClean[nParticles][nCentBins][nEtaBins][njetRad];
TString jet_hName2DClean[nParticles][nCentBins][nEtaBins][njetRad];
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
UInt_t    ftracks_cutBit=0;
ULong64_t ftracks_gid=0;
Double_t  ftracks_eventtime=0;
Float_t   ftracks_intrate=0;
Float_t   ftracks_dEdx=0;
Int_t     ftracks_sign=0;
Float_t   ftracks_ptot=0;
Float_t   ftracks_p=0;
Float_t   ftracks_pT=0;
Float_t   ftracks_eta=0;
Float_t   ftracks_phi=0;
Float_t   ftracks_cent=0;
Int_t     ftracks_ncltpc=0;
Float_t   ftracks_dcaxy=0;
Float_t   ftracks_dcaz=0;
Float_t   ftracks_cRows=0;
Float_t   ftracks_chi2tpc=0;
//
// Jet constituent tree branches
UInt_t     fjetConst_cutBit=0;
ULong64_t  fjetConst_gid=0;
Float_t    fjetConst_radius=0;
Float_t    fjetConst_area=0;
Double_t   fjetConst_maxpt=0;
Float_t    fjetConst_ptsub=0;
//Double_t fjetConst_eventtime=0;
//Float_t  fjetConst_intrate=0;
Float_t    fjetConst_dEdx=0;
Int_t      fjetConst_sign=0;
Float_t    fjetConst_ptot=0;
Float_t    fjetConst_p=0;
Float_t    fjetConst_pT=0;
Float_t    fjetConst_eta=0;
Float_t    fjetConst_phi=0;
Float_t    fjetConst_cent=0;
//Int_t    fjetConst_ncltpc=0;
//Float_t  fjetConst_dcaxy=0;
//Float_t  fjetConst_dcaz=0;
//Float_t  fjetConst_cRows=0;
//Float_t  fjetConst_chi2tpc=0;
//
// Clean Tree Branches
UInt_t    ffArmPodTree_purity = 0;
ULong64_t ffArmPodTree_gid=0;
Double_t  ffArmPodTree_eventtime=0;
Float_t   ffArmPodTree_intrate=0;
Float_t   ffArmPodTree_cent=0;
Float_t   ffArmPodTree_qt=0;
Float_t   ffArmPodTree_alfa=0;
Char_t    ffArmPodTree_v0haspixel=0;
Char_t    ffArmPodTree_piFromK0=0;
//
UInt_t    ffArmPodTree_cutBit0=0;
Float_t   ffArmPodTree_dEdx0=0;
Int_t     ffArmPodTree_sign0=0;
Float_t   ffArmPodTree_ptot0=0;
Float_t   ffArmPodTree_p0=0;
Float_t   ffArmPodTree_pT0=0;
Float_t   ffArmPodTree_eta0=0;
Float_t   ffArmPodTree_phi0=0;
Float_t   ffArmPodTree_nSigmasPiTOF0=0;
Float_t   ffArmPodTree_nSigmasPrTOF0=0;
//
UInt_t  ffArmPodTree_cutBit1=0;
Float_t ffArmPodTree_dEdx1=0;
Int_t   ffArmPodTree_sign1=0;
Float_t ffArmPodTree_ptot1=0;
Float_t ffArmPodTree_p1=0;
Float_t ffArmPodTree_pT1=0;
Float_t ffArmPodTree_eta1=0;
Float_t ffArmPodTree_phi1=0;
Float_t ffArmPodTree_nSigmasPiTOF1=0;
Float_t ffArmPodTree_nSigmasPrTOF1=0;
//
// dscaled tree
Float_t   fdscaled_dcaxy = 0;
Float_t   fdscaled_dcaz = 0;
UInt_t    fdscaled_cutBit = 0;
Float_t   fdscaled_cent = 0;
ULong64_t fdscaled_gid = 0;
Float_t   fdscaled_intrate = 0;
Double_t  fdscaled_eventtime = 0;
Double_t  fdscaled_vz = 0;
Double_t  fdscaled_sharedTPCClusters = 0;
Int_t     fdscaled_ncltpc = 0;
Int_t     fdscaled_tpcSignalN = 0;
Int_t     fdscaled_cRows = 0;
Double_t  fdscaled_lengthInActiveZone = 0;
Double_t  fdscaled_phi = 0;
Double_t  fdscaled_phiTPC = 0;
Double_t  fdscaled_phi85 = 0;
Double_t  fdscaled_itsdEdx = 0;
Double_t  fdscaled_trddEdx = 0;
Float_t   fdscaled_fdTPC = 0;
Float_t   fdscaled_fzTPC = 0;
Float_t   fdscaled_fCdd = 0;
Float_t   fdscaled_fCdz = 0;
Float_t   fdscaled_fCzz = 0;
Int_t     fdscaled_nclits = 0;
Int_t     fdscaled_ncltrd = 0;
Double_t  fdscaled_chi2tpc = 0;
Double_t  fdscaled_chi2its = 0;
Double_t  fdscaled_chi2trd = 0;
//
// evets tree variables
Int_t     fevents_run=0;
Float_t   fevents_intrate=0;
Double_t  fevents_bField=0;
ULong64_t fevents_gid=0;
Double_t  fevents_timestamp=0;
ULong64_t fevents_triggerMask=0;
Double_t  fevents_vz=0;
Double_t  fevents_tpcvz=0;
Double_t  fevents_spdvz=0;
Int_t     fevents_tpcMult=0;
Short_t   fevents_eventMult=0;
Int_t     fevents_eventMultESD=0;
Int_t     fevents_nTracksStored=0;
Int_t     fevents_primMult=0;
Int_t     fevents_tpcTrackBeforeClean=0;
Int_t     fevents_itsTracklets=0;
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
TVectorF *fevents_itsVertexInfo=0x0;
//
// jet events tree variables
Int_t      fjetEvents_run=0;
Float_t    fjetEvents_cent=0;
Float_t    fjetEvents_rhoFJ=0;
Char_t     fjetEvents_accjet=0;
Char_t     fjetEvents_realjet=0;
Float_t    fjetEvents_intrate=0;
Float_t    fjetEvents_bField=0;
ULong64_t  fjetEvents_gid=0;
//Double_t fjetEvents_timestamp=0;
ULong64_t  fjetEvents_triggerMask=0;
Double_t   fjetEvents_vz=0;
Float_t    fjetEvents_tpcvz=0;
Float_t    fjetEvents_spdvz=0;
Int_t      fjetEvents_tpcMult=0;
Int_t      fjetEvents_eventMult=0;
Int_t      fjetEvents_eventMultESD=0;
//nt_t     fjetEvents_nTracksStored=0;
Float_t    fjetEvents_primMult=0;
Int_t      fjetEvents_tpcTrackBeforeClean=0;
Int_t      fjetEvents_itsTracklets=0;
/*
TVectorF *fjetEvents_centrality=0x0;
TVectorF *fjetEvents_tZeroMult=0x0;
TVectorF *fjetEvents_vZeroMult=0x0;
TVectorF *fjetEvents_itsClustersPerLayer=0x0;
TVectorF *fjetEvents_trackCounters=0x0;
TVectorF *fjetEvents_trackdEdxRatio=0x0;
TVectorF *fjetEvents_trackNcl=0x0;
TVectorF *fjetEvents_trackChi2=0x0;
TVectorF *fjetEvents_trackMatchEff=0x0;
TVectorF *fjetEvents_trackTPCCountersZ=0x0;
TVectorF *fjetEvents_phiCountA=0x0;
TVectorF *fjetEvents_phiCountC=0x0;
TVectorF *fjetEvents_phiCountAITS=0x0;
TVectorF *fjetEvents_phiCountCITS=0x0;
TVectorF *fjetEvents_phiCountAITSOnly=0x0;
TVectorF *fjetEvents_phiCountCITSOnly=0x0;
TVectorF *fjetEvents_tpcVertexInfo=0x0;
TVectorF *fjetEvents_itsVertexInfo=0x0;
*/
//
// highPt tree varianles
ULong64_t     fhighPt_gid=0;
Int_t         fhighPt_selectionPtMask=0;
Float_t       fhighPt_centralityF=0;
Int_t         fhighPt_runNumber=0;
Int_t         fhighPt_evtTimeStamp=0;
Int_t         fhighPt_evtNumberInFile=0;
Float_t       fhighPt_Bz=0;
Int_t         fhighPt_IRtot=0;
Int_t         fhighPt_IRint2=0;
Int_t         fhighPt_mult=0;
TObjString   *fhighPt_triggerClass=0;
AliESDtrack  *fhighPt_esdTrack=0;
TObjString   *fhighPt_fileName=0;
AliESDVertex *fhighPt_vtxESD=0;
//
//
const Int_t fnCutBins=10;
Int_t fCutArr[fnCutBins];
Int_t fSystSet=0;
TString fPlotVS="";
TString fCutON="";
TString inputDataTree     = "tracks";
TString inputJetConstTree = "jetsFJconst";
TString inputCleanTree    = "fArmPodTree";
TString inputEventTree    = "events";
TString inputJetEventTree = "jetsEventInfo";
TString inputDscaledTree = "dscaled";
TString inputHighPtTree  = "highPt";
TString inputHists     = "cleanHists";
TString inputHistsFile = "";
TFile *fdata=NULL;
TFile *fhist=NULL;
TFile *fdEdxMap=NULL;
TTree *dataTree=NULL, *jetConstTree=NULL, *armtree=NULL, *eventtree=NULL, *jetEventtree=NULL, *dscaltree=NULL, *highPttree=NULL;
THnSparse *fhnExpected=NULL;
TTreeSRedirector *treeStream=0, *histStream=0, *jetHistStream=0, *debugStream=0, *splinesStream=0;
TH1D  *hEta=0x0, *hCent=0x0, *hMom=0x0;
TCutG *pionCutG=0, *antiProtonCutG=0, *protonCutG=0;
TStopwatch timer;
AliNDLocalRegression *hdEdxAShifttMNTglDist_meanGFitAll=NULL;
AliNDLocalRegression *hdEdxAShifttMNTglDist_meanGFitNoPileUp=NULL;
//
TH2F *htracks_dcaxy2D=NULL;
TH2F *htracks_dcaz2D=NULL;
TH2F *htracks_ncltpc2D=NULL;
TH2F *htracks_dcaxy2D_After=NULL;
TH2F *htracks_dcaz2D_After=NULL;
/*
TH2F *hjetConst_dcaxy2D=NULL;
TH2F *hjetConst_dcaz2D=NULL;
TH2F *hjetConst_ncltpc2D=NULL;
TH2F *hjetConst_dcaxy2D_After=NULL;
TH2F *hjetConst_dcaz2D_After=NULL;
*/
//
TH2F *hevents_pileUpV02D=NULL;
TH2F *hevents_pileUpV02D_inbunchCut=NULL;
TH2F *hevents_V0M_CL0=NULL;
TH2F *hevents_V0M_CL0_After=NULL;
TH2F *hevents_itsLayer0V02D_inbunchCut=NULL;
TH2F *hevents_itsLayer0V02D=NULL;
TH2F *hevents_pileUpITS2D_inbunchCut=NULL;
TH2F *hevents_pileUpITS2D=NULL;
TH2F *hevents_secMultITS0_primMultITS=NULL;
TH1F *hevents_pileUpV01D=NULL;
TH1F *hevents_pileUpITS1D=NULL;
TH1F *hevents_ITSTPCeff=NULL;
// Jet events
TH2F *jet_hevents_pileUpV02D=NULL;
TH2F *jet_hevents_pileUpV02D_inbunchCut=NULL;
TH2F *jet_hevents_V0M_CL0=NULL;
TH2F *jet_hevents_V0M_CL0_After=NULL;
TH2F *jet_hevents_itsLayer0V02D_inbunchCut=NULL;
TH2F *jet_hevents_itsLayer0V02D=NULL;
TH2F *jet_hevents_pileUpITS2D_inbunchCut=NULL;
TH2F *jet_hevents_pileUpITS2D=NULL;
TH2F *jet_hevents_secMultITS0_primMultITS=NULL;
TH1F *jet_hevents_pileUpV01D=NULL;
TH1F *jet_hevents_pileUpITS1D=NULL;
TH1F *jet_hevents_ITSTPCeff=NULL;
//
TH1F *hdscaled_sharedTPCClusters1D=NULL;
TH1F *hdscaled_tpcSignalN1D=NULL;
TH1F *hdscaled_lengthInActiveZone1D=NULL;
TH1F *hdscaled_cRows1D=NULL;
TH1F *hdscaled_nclits1D=NULL;
TH1F *hdscaled_ncltrd1D=NULL;
TH1F *hdscaled_chi2tpc1D=NULL;
TH1F *hdscaled_chi2its1D=NULL;
TH1F *hdscaled_chi2trd1D=NULL;
TH1F *hdscaled_vz1D=NULL;
TH1F *htracks_dcaxy1D=NULL;
TH1F *htracks_dcaz1D=NULL;
TH1F *htracks_ncltpc1D=NULL;
TH1F *htracks_eta1D=NULL;
TH1F *htracks_cRows1D=NULL;
TH1F *htracks_cent1D=NULL;
TH1F *htracks_phi1D=NULL;
TH1F *htracks_chi2tpc1D=NULL;
TH1F *htracks_vz1D=NULL;
/*
TH1F *hjetConst_dcaxy1D=NULL;
TH1F *hjetConst_dcaz1D=NULL;
TH1F *hjetConst_ncltpc1D=NULL;
*/
TH1F *hjetConst_eta1D=NULL;
//TH1F *hjetConst_cRows1D=NULL;
TH1F *hjetConst_cent1D=NULL;
TH1F *hjetConst_phi1D=NULL;
//TH1F *hjetConst_chi2tpc1D=NULL;
TH1F *hjetConst_vz1D=NULL;
//
TH2F *hhighPt_dEdxPtot=NULL;
//
//
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
  kTrackProbKaTOF=18,
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
  kCutNSigmaTOFLoose=16,
  kCutNSigmaTOFLoose2=17
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


Bool_t fQAtimeSeriesExist = kTRUE;
Int_t    fInPileupTPCCut  = 0;
Double_t fInTimeSeriesEff = 0.;
Int_t    fPeriod          = 0;
//
TString dEdxMapStr = "";
TString timeSeriesMap = "";
TString MapDirName = "/home/marsland/Desktop/ubuntu_desktop/workdir/Data/IlyaEos";
Double_t testEntries  = 100000;
// Double_t testEntries  = -1.;
// Bool_t testIntegratedHist=kTRUE;
Bool_t testIntegratedHist=kFALSE;
Bool_t fillTIdenTree=kTRUE;



void RealData_FilterTreesMakeHists_Run2(TString plotVS, TString cutON, Int_t systSet, TString dFile, Int_t period, Int_t kPileupTPCCut, Double_t kTimeSeriesEff)
{
  //
  // plotVS: "ptot", "pT", "p"
  //

  /*

  meld /home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/code/RealData_FilterTreesMakeHists_Run2.C /home/marsland/Desktop/ubuntu_desktop/github/TIdentity3D/TIdentity/macros/RealData_FilterTreesMakeHists_Run2.C

  cd /home/marsland/Desktop/ubuntu_desktop/workdir/fourthMomentAnalysis/test_data
  aliroot -l
  TString file = "/home/marsland/Desktop/ubuntu_desktop/workdir/Data/IlyaEos/AnalysisResults_data.root"
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/code/RealData_FilterTreesMakeHists_Run2.C+
  RealData_FilterTreesMakeHists_Run2("ptot","p",0,file ,1, 1, 0.85)   // full correction and selection

  */
  //
  fInPileupTPCCut  = kPileupTPCCut;
  fInTimeSeriesEff = kTimeSeriesEff;
  cout << " fInPileupTPCCut     = " <<  fInPileupTPCCut << endl;
  cout << " fInTimeSeriesEff = " <<  fInTimeSeriesEff << endl;
  cout << " setting          = " <<  systSet << endl;
  //
  fSystSet         = systSet;
  fPlotVS          = plotVS;
  fCutON           = cutON;
  fPeriod          = period;
  cout << "  ------------ > dEdx Map = " << dEdxMapStr << endl;
  //
  fdata      = TFile::Open(dFile);
  armtree    = (TTree*)fdata->Get(inputCleanTree);
  dataTree   = (TTree*)fdata->Get(inputDataTree);
  jetConstTree   = (TTree*)fdata->Get(inputJetConstTree);
  eventtree  = (TTree*)fdata->Get(inputEventTree);
  jetEventtree  = (TTree*)fdata->Get(inputJetEventTree);
  dscaltree  = (TTree*)fdata->Get(inputDscaledTree);
  highPttree = (TTree*)fdata->Get(inputHighPtTree);
  //
  //
  SetBranchAddresses();
  //
  InitInitials();
  //
  if (dataTree) dataTree->BuildIndex("gid");
  if (jetConstTree) {
    jetConstTree->BuildIndex("gid");
    dataTree->AddFriend(jetConstTree,"jetConsts");
    jetConstTree->AddFriend(dataTree,"jetConsts");
  }
  if (eventtree) {
    eventtree->BuildIndex("gid");
    dataTree->AddFriend(eventtree,"event");
    jetConstTree->AddFriend(eventtree,"event");
  }
  if (jetEventtree) {
    jetEventtree->BuildIndex("gid");
    dataTree->AddFriend(jetEventtree,"jetEvent");
    jetConstTree->AddFriend(jetEventtree,"jetEvent");
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
  ProcessJetConstHists();
  ProcessEventTree();
  ProcessJetEventTree();
  ProcessCleanSamples();
  ProcessDScaledTree();
  ProcessHighPtTree();
  WriteHistsToFile();
  if (treeStream)  delete treeStream;
  if (debugStream) delete debugStream;
  if (histStream)  delete histStream;
  if (jetHistStream)  delete jetHistStream;
  //
  for (Int_t icent=0; icent<nCentBins; icent++) delete tidenTreeStream[icent];


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
    // --------------------------------------------------------------------------------
    // Get the time series map for a given run and check if the run is in good run list
    // --------------------------------------------------------------------------------
    //
    static Int_t runCache = -1;
    static Int_t acceptRun = -1;
    if (runCache!=Int_t(fevents_run)) {
      runCache=Int_t(fevents_run);
      if (fPeriod==0) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2015/LHC15o/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==1) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18q/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==2) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18r/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      cout << "run number = " << runCache << "  ------------ > timeSeries Map = " <<  timeSeriesMap << endl;
      //
      // select run number
      // if (fPeriod==0){
      //   for (Int_t irun = 0; irun<Int_t(sizeof(runList)/sizeof(runList[0])); irun++){
      //     if (runList[i]==runCache) {
      //       acceptRun=1;
      //       cout << " --- run is in the good run list from RCT --- GO AHEAD " << endl;
      //       break;
      //     }
      //   }
      // }
      //
      // get the time series
      TFile *fQAtimeSeries = TFile::Open(timeSeriesMap);
      if (!fQAtimeSeries) {
        std::cout << " yolun acik ola --> file does not exist  " << std::endl;
        fQAtimeSeriesExist = kFALSE;
      } else {
        if (fPeriod==0){
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.binMedian");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.binMedian");     // weighting wrt track
        } else {
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.mean");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.mean");     // weighting wrt track
        }
      }
    }
    //
    // Run selection
    // if (fPeriod==0 && acceptRun<0) { cout << " !!! run is not in the good run list from RCT !!! " << endl; break;}
    //
    // --------------------------------------------------------------------------------
    // Apply event selection --> tpcPileUpCut && timeSeriesCut && inbunchPileUpCut && spdVzCut
    // --------------------------------------------------------------------------------
    //
    static ULong64_t gidCache = -1;
    static Int_t acceptEvent = -1;
    if (gidCache!=ULong64_t(ftracks_gid)) {
      eventCount++;
      gidCache=ULong64_t(ftracks_gid);
      //
      //  TPC ITS matching eff from time series
      Double_t fITSTPCeffEvent = (fQAtimeSeriesExist && fgrTimeSeriesEventWeighted)   ? fgrTimeSeriesEventWeighted  ->Eval(fevents_timestamp) : 1;
      Double_t fITSTPCeffTrack = (fQAtimeSeriesExist && fgrTimeSeriesNTracksWeighted) ? fgrTimeSeriesNTracksWeighted->Eval(fevents_timestamp) : 1;
      //
      // PileUp selection:  chain->SetAlias("multTPCITS","((multSSD+multSDD)/2.38)");   // chain->SetAlias("cut1000","tpcMult-multTPCITS<1000");
      Double_t multSSD = (*fevents_itsClustersPerLayer)[4]+(*fevents_itsClustersPerLayer)[5];
      Double_t multSDD = (*fevents_itsClustersPerLayer)[2]+(*fevents_itsClustersPerLayer)[3];
      Double_t multTPC = TMath::Max(fevents_tpcMult,fevents_tpcTrackBeforeClean);
      Double_t pileUp1DITS = (multSSD+multSDD)/2.38;
      Double_t multV0 = fevents_vZeroMult->Sum();
      Double_t multT0 = fevents_tZeroMult->Sum();
      Double_t primMult    = fevents_primMult;
      Double_t trackTgl    = TMath::Abs(TMath::SinH(ftracks_eta));
      //
      // in bunch pileup cut
      Double_t primMultITS = (*fevents_itsVertexInfo)[2];
      Double_t secMultITS0 = (*fevents_itsVertexInfo)[3];
      Double_t secMultITS  = secMultITS0 - inbunchCutSlope*primMultITS;
      //
      // Define the cuts
      Bool_t tpcPileUpCut     = ( ftracks_cutBit & (1 << fInPileupTPCCut) );
      Bool_t timeSeriesCut    = (fITSTPCeffTrack > fInTimeSeriesEff);
      Bool_t inbunchPileUpCut = (secMultITS < (inbunchCutRMS*primMultITS+inbunchCutOffSet));
      Bool_t spdVzCut         = (fevents_spdvz!=0 && TMath::Abs(fevents_spdvz)<10);
      Bool_t inbunchPileUpTailCut = (pileUp1DITS<14500);
      //
      // trigger selections
      Bool_t bMB          = ((fevents_triggerMask&0x1)>0);
      Bool_t bCentral     = ((fevents_triggerMask&0x80000)>0);
      Bool_t bSemiCentral = ((fevents_triggerMask&0x100000)>0);
      Bool_t bMB10        = ((fevents_triggerMask&0x10)>0);
      Bool_t bAlltiggers = (bMB || bCentral || bSemiCentral || bMB10 );
      //
      // inbunch pile hists
      hevents_secMultITS0_primMultITS->Fill(primMultITS,secMultITS0);
      if (!inbunchPileUpCut && spdVzCut) hevents_pileUpV02D_inbunchCut   ->Fill(multTPC,multV0);
      if (!inbunchPileUpCut && spdVzCut) hevents_pileUpITS2D_inbunchCut  ->Fill(multTPC,pileUp1DITS);
      if (!inbunchPileUpCut && spdVzCut) hevents_itsLayer0V02D_inbunchCut->Fill((*fevents_itsClustersPerLayer)[0],multV0);
      hevents_pileUpITS2D_inbunchCut->SetMarkerStyle(20);
      hevents_pileUpITS2D_inbunchCut->SetMarkerSize(0.5);
      hevents_pileUpITS2D_inbunchCut->SetMarkerColor(kRed+1);
      hevents_itsLayer0V02D_inbunchCut->SetMarkerStyle(20);
      hevents_itsLayer0V02D_inbunchCut->SetMarkerSize(0.5);
      hevents_itsLayer0V02D_inbunchCut->SetMarkerColor(kRed+1);
      //
      // centrality vs centrality
      hevents_V0M_CL0->Fill((*fevents_centrality)[0],(*fevents_centrality)[1]);
      //
      // CL0 vs V0M cuts
      Bool_t dioganalCut = ( (*fevents_centrality)[1] < (82./75.)*(*fevents_centrality)[0]+8. && (*fevents_centrality)[1] > (70./82.)* ((*fevents_centrality)[0]-8.) );
      //
      // pile up and time series cut
      // if ( tpcPileUpCut && timeSeriesCut && inbunchPileUpCut && spdVzCut && bAlltiggers)  acceptEvent=1;
      if ( tpcPileUpCut && timeSeriesCut && inbunchPileUpCut && spdVzCut && bAlltiggers && inbunchPileUpTailCut && dioganalCut)  acceptEvent=1;
      else acceptEvent=-1;
      //
      // Dump debug histograms
      if (multTPC>0 && acceptEvent>0){
        if (pileUp1DITS/multTPC>0.005) hevents_pileUpV01D->Fill(multV0/multTPC/2.55);
        if (pileUp1DITS/multTPC>0.005) hevents_pileUpITS1D->Fill(pileUp1DITS/multTPC);
        hevents_pileUpV02D   ->Fill(multTPC,multV0);
        hevents_pileUpITS2D  ->Fill(multTPC,pileUp1DITS);
        hevents_itsLayer0V02D->Fill((*fevents_itsClustersPerLayer)[0],multV0);
        hevents_ITSTPCeff    ->Fill(fITSTPCeffTrack);
        hevents_V0M_CL0_After->Fill((*fevents_centrality)[0],(*fevents_centrality)[1]);
        //
        htracks_cent1D->Fill(ftracks_cent);
        htracks_vz1D->Fill(fevents_vz);

      }
    }
    //
    // Event Cuts --> pile up and time series cut
    if (acceptEvent<0) continue;
    Double_t primMult = fevents_primMult;
    Double_t trackTgl = TMath::Abs(TMath::SinH(ftracks_eta));
    //
    // dump some debug histogram
    if (ftracks_dcaxy>-10 && ftracks_dcaxy<10)   htracks_dcaxy2D ->Fill(ftracks_pT,ftracks_dcaxy);
    if (ftracks_dcaz>-10  && ftracks_dcaz<10)    htracks_dcaz2D  ->Fill(ftracks_pT,ftracks_dcaz);
    if (ftracks_ncltpc>30 && ftracks_ncltpc<170) htracks_ncltpc2D->Fill(ftracks_pT,Float_t(ftracks_ncltpc));
    //
    // Retrieve cut setting
    Bool_t systCut = ApplyTreeSelection(fSystSet, ftracks_cutBit);
    Bool_t etaAcc = (ftracks_eta >= -0.8 && ftracks_eta <= 0.8);
    Bool_t momAcc = (ftracks_p >= 0.2    && ftracks_p <= 3.2);
    //
    // dump tidentree
    if (systCut && etaAcc && momAcc && fillTIdenTree){
      //
      //
      // Fill final histograms
      if ( fSystSet==21 && TMath::Abs(ftracks_dcaxy) > 3.2 ) continue;
      if ( fSystSet==21 && TMath::Abs(ftracks_dcaz)  > 2.4 ) continue;
      htracks_dcaxy1D  ->Fill(ftracks_dcaxy);
      htracks_dcaz1D   ->Fill(ftracks_dcaz);
      htracks_ncltpc1D ->Fill(ftracks_ncltpc);
      htracks_eta1D    ->Fill(ftracks_eta);
      htracks_cRows1D  ->Fill(ftracks_cRows);
      htracks_phi1D    ->Fill(ftracks_phi);
      htracks_chi2tpc1D->Fill(ftracks_chi2tpc);
      htracks_dcaxy2D_After ->Fill(ftracks_pT,ftracks_dcaxy);
      htracks_dcaz2D_After  ->Fill(ftracks_pT,ftracks_dcaz);
      //
      // fill trees for a given cent
      for (Int_t icent=0; icent<nCentBins; icent++){
        Bool_t centAcc = (ftracks_cent>=centBinning[icent] && ftracks_cent<centBinning[icent+1]);
        if ( centAcc ){
          //
          tidenTreeStream[icent]->GetFile()->cd();
          //
          // tracks tree
          (*tidenTreeStream[icent])<<"tracks"<<
          "gid="                  << ftracks_gid          <<  //  global event ID
          "dEdx="                 << ftracks_dEdx         <<  //  dEdx of the track
          "sign="                 << ftracks_sign         <<  //  charge
          "ptot="                 << ftracks_ptot         <<  //  TPC momentum
          "eta="                  << ftracks_eta          <<  //  eta
          "cent="                 << ftracks_cent         <<  //  centrality
          "\n";
          continue;
        }
      }
    }
    if (systCut && testIntegratedHist){
      treeStream->GetFile()->cd();
      (*treeStream)<<"tracks"<<
      "gid="                  << ftracks_gid          <<  //  global event ID
      "p="                    << ftracks_p            <<  //  vertex momentum
      "intrate="              << ftracks_intrate      <<  //  interaction rate
      "cutBit="               << ftracks_cutBit       <<  //  Systematic Cuts
      "dEdx="                 << ftracks_dEdx         <<  //  dEdx of the track
      "sign="                 << ftracks_sign         <<  //  charge
      "ptot="                 << ftracks_ptot         <<  //  TPC momentum
      "pT="                   << ftracks_pT           <<  //  transverse momentum
      "eta="                  << ftracks_eta          <<  //  eta
      "phi="                  << ftracks_phi          <<  //  eta
      "cent="                 << ftracks_cent         <<  //  centrality
      //
      "fevents_tpcVertexInfo.=" << fevents_tpcVertexInfo <<
      "fevents_itsClustersPerLayer.=" << fevents_itsClustersPerLayer <<
      "primMult="             << primMult             <<  // interaction rate
      "tgl="                  << trackTgl             <<  // interaction rate
      "cRows="                << ftracks_cRows        <<  // interaction rate
      "chi2tpc="              << ftracks_chi2tpc      <<  // interaction rate
      "dcaz="                 << ftracks_dcaz         <<  // interaction rate
      "dcaxy="                << ftracks_dcaxy        <<  // interaction rate
      "\n";
    }
    //
    // apply dca and nclusters cut
    if (TMath::Abs(ftracks_dcaz)>3) continue;
    if (TMath::Abs(ftracks_dcaxy)>3) continue;
    //
    // Prepare dEdx histograms
    for (Int_t icent=0; icent<nCentBins; icent++){
      for (Int_t ieta=0; ieta<nEtaBins; ieta++){
        //
        if (testIntegratedHist && ieta>0) continue;
        //
        Bool_t etaCentString = (ftracks_eta>=etaBinning[ieta] && ftracks_eta<etaBinning[ieta+1] && ftracks_cent>=centBinning[icent] && ftracks_cent<centBinning[icent+1]);
        if (testIntegratedHist && ieta==0) {
          etaCentString = (ftracks_cent>=centBinning[icent] && ftracks_cent<centBinning[icent+1]);
        }
        //
        Bool_t cleanKaCutTOFTRD = ((ftracks_cutBit >> kCleanKaTOFTRD) & 1);
        Bool_t cleanKaCutBayes = ((ftracks_cutBit >> kTrackProbKaTOF) & 1);
        Bool_t cleanDeCutTOF = ((ftracks_cutBit >> kCleanDeTOF) & 1);
        Bool_t parPos = ftracks_sign>0.;
        Bool_t parNeg = ftracks_sign<0.;
        //
        Bool_t prTOF = ((ftracks_cutBit >> kCleanPrTOF) & 1);
        Bool_t kaTOF = ((ftracks_cutBit >> kCleanKaTOF) & 1);
        //
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
        if (etaCentString && vertexPcut && systCut && !prTOF)           h2DallBkgTOF[icent][ieta]   ->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && !prTOF && parPos) h2DallBkgTOFPos[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && !prTOF && parNeg) h2DallBkgTOFNeg[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        //
        if (etaCentString && vertexPcut && systCut && cleanKaCutBayes)  h2DCleanKaBayes[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && cleanKaCutTOFTRD) h2DCleanKaTOFTRD[icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);

        if (etaCentString && vertexPcut && systCut && cleanDeCutTOF)    h2DClean[4][icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
        if (etaCentString && vertexPcut && systCut && cleanKaCutTOFTRD) h2DClean[2][icent][ieta]->Fill(ftracks_ptot,ftracks_dEdx);
      }
    }

  }  // tree loop
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessEventTree DONE ========= #events = " << eventCount << std::endl;

}
//____________________________________________________________________________________________________________
void ProcessJetConstHists()
{

  timer.Reset(); timer.Start();
  std::cout << " ========= ProcessJetConstHists ========= " << std::endl;
  Double_t nTreeEntriesAll = jetConstTree -> GetEntries();
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
    jetConstTree -> GetEntry(i);
    if(i%Int_t(nTreeEntriesAll/10) == 0) {
      cout << i << "   fjetConst_dEdx = " << fjetConst_dEdx      << "         fjetConst_eta        = " << fjetConst_eta << endl; //"  fjetConst_dcaxy = " << fjetConst_dcaxy << endl;
      cout << "   fjetConst_gid       = " << fjetConst_gid       << "    -->  fevents_gid        = " <<  fevents_gid              << endl;
      cout << "   fjetConst_gid       = " << fjetConst_gid       << "    -->  fjetEvents_gid        = " <<  fjetEvents_gid              << endl;
      cout << "    -->  fevents_timestamp  = " <<  fevents_timestamp       << endl; //"   fjetConst_eventtime = " << fjetConst_eventtime
      cout << "   fjetConst_cent      = " << fjetConst_cent      << "    -->  fevents_centrality = " <<  (*fevents_centrality)[0] << endl;
    }

    //
    // --------------------------------------------------------------------------------
    // Get the time series map for a given run and check if the run is in good run list
    // --------------------------------------------------------------------------------
    //
    static Int_t runCache = -1;
    static Int_t acceptRun = -1;
    if (runCache!=Int_t(fevents_run)) {
      runCache=Int_t(fevents_run);
      if (fPeriod==0) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2015/LHC15o/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==1) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18q/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==2) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18r/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      cout << "run number = " << runCache << "  ------------ > timeSeries Map = " <<  timeSeriesMap << endl;
      //
      // select run number
      // if (fPeriod==0){
      //   for (Int_t irun = 0; irun<Int_t(sizeof(runList)/sizeof(runList[0])); irun++){
      //     if (runList[i]==runCache) {
      //       acceptRun=1;
      //       cout << " --- run is in the good run list from RCT --- GO AHEAD " << endl;
      //       break;
      //     }
      //   }
      // }
      //
      // get the time series
      TFile *fQAtimeSeries = TFile::Open(timeSeriesMap);
      if (!fQAtimeSeries) {
        std::cout << " yolun acik ola --> file does not exist  " << std::endl;
        fQAtimeSeriesExist = kFALSE;
      } else {
        if (fPeriod==0){
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.binMedian");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.binMedian");     // weighting wrt track
        } else {
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.mean");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.mean");     // weighting wrt track
        }
      }
    }
    //
    // Run selection
    // if (fPeriod==0 && acceptRun<0) { cout << " !!! run is not in the good run list from RCT !!! " << endl; break;}
    //
    // --------------------------------------------------------------------------------
    // Apply event selection --> tpcPileUpCut && timeSeriesCut && inbunchPileUpCut && spdVzCut
    // --------------------------------------------------------------------------------
    //
    static ULong64_t gidCache = -1;
    static Int_t acceptEvent = -1;
    if (gidCache!=ULong64_t(fjetConst_gid)) {
      eventCount++;
      gidCache=ULong64_t(fjetConst_gid);
      //
      //  TPC ITS matching eff from time series
      Double_t fITSTPCeffEvent = (fQAtimeSeriesExist && fgrTimeSeriesEventWeighted)   ? fgrTimeSeriesEventWeighted  ->Eval(fevents_timestamp) : 1;
      Double_t fITSTPCeffTrack = (fQAtimeSeriesExist && fgrTimeSeriesNTracksWeighted) ? fgrTimeSeriesNTracksWeighted->Eval(fevents_timestamp) : 1;
      //
      // PileUp selection:  chain->SetAlias("multTPCITS","((multSSD+multSDD)/2.38)");   // chain->SetAlias("cut1000","tpcMult-multTPCITS<1000");
      Double_t multSSD = (*fevents_itsClustersPerLayer)[4]+(*fevents_itsClustersPerLayer)[5];
      Double_t multSDD = (*fevents_itsClustersPerLayer)[2]+(*fevents_itsClustersPerLayer)[3];
      Double_t multTPC = TMath::Max(fevents_tpcMult,fevents_tpcTrackBeforeClean);
      Double_t pileUp1DITS = (multSSD+multSDD)/2.38;
      Double_t multV0 = fevents_vZeroMult->Sum();
      Double_t multT0 = fevents_tZeroMult->Sum();
      Double_t primMult    = fevents_primMult;
      Double_t trackTgl    = TMath::Abs(TMath::SinH(fjetConst_eta));
      //
      // in bunch pileup cut
      Double_t primMultITS = (*fevents_itsVertexInfo)[2];
      Double_t secMultITS0 = (*fevents_itsVertexInfo)[3];
      Double_t secMultITS  = secMultITS0 - inbunchCutSlope*primMultITS;
      //
      // Define the cuts
      Bool_t tpcPileUpCut     = ( fjetConst_cutBit & (1 << fInPileupTPCCut) );
      Bool_t timeSeriesCut    = (fITSTPCeffTrack > fInTimeSeriesEff);
      Bool_t inbunchPileUpCut = (secMultITS < (inbunchCutRMS*primMultITS+inbunchCutOffSet));
      Bool_t spdVzCut         = (fevents_spdvz!=0 && TMath::Abs(fevents_spdvz)<10);
      Bool_t inbunchPileUpTailCut = (pileUp1DITS<14500);
      //
      // trigger selections
      Bool_t bMB          = ((fevents_triggerMask&0x1)>0);
      Bool_t bCentral     = ((fevents_triggerMask&0x80000)>0);
      Bool_t bSemiCentral = ((fevents_triggerMask&0x100000)>0);
      Bool_t bMB10        = ((fevents_triggerMask&0x10)>0);
      Bool_t bAlltiggers = (bMB || bCentral || bSemiCentral || bMB10 );
      //
      // inbunch pile hists
      jet_hevents_secMultITS0_primMultITS->Fill(primMultITS,secMultITS0);
      if (!inbunchPileUpCut && spdVzCut) jet_hevents_pileUpV02D_inbunchCut   ->Fill(multTPC,multV0);
      if (!inbunchPileUpCut && spdVzCut) jet_hevents_pileUpITS2D_inbunchCut  ->Fill(multTPC,pileUp1DITS);
      if (!inbunchPileUpCut && spdVzCut) jet_hevents_itsLayer0V02D_inbunchCut->Fill((*fevents_itsClustersPerLayer)[0],multV0);
      jet_hevents_pileUpITS2D_inbunchCut->SetMarkerStyle(20);
      jet_hevents_pileUpITS2D_inbunchCut->SetMarkerSize(0.5);
      jet_hevents_pileUpITS2D_inbunchCut->SetMarkerColor(kRed+1);
      jet_hevents_itsLayer0V02D_inbunchCut->SetMarkerStyle(20);
      jet_hevents_itsLayer0V02D_inbunchCut->SetMarkerSize(0.5);
      jet_hevents_itsLayer0V02D_inbunchCut->SetMarkerColor(kRed+1);
      //
      // centrality vs centrality
      jet_hevents_V0M_CL0->Fill((*fevents_centrality)[0],(*fevents_centrality)[1]);
      //
      // CL0 vs V0M cuts
      Bool_t dioganalCut = ( (*fevents_centrality)[1] < (82./75.)*(*fevents_centrality)[0]+8. && (*fevents_centrality)[1] > (70./82.)* ((*fevents_centrality)[0]-8.) );
      //
      // pile up and time series cut
      // if ( tpcPileUpCut && timeSeriesCut && inbunchPileUpCut && spdVzCut && bAlltiggers)  acceptEvent=1;
      if ( tpcPileUpCut && timeSeriesCut && inbunchPileUpCut && spdVzCut && inbunchPileUpTailCut && dioganalCut && bAlltiggers)  acceptEvent=1;
      else acceptEvent=-1;
      //
      // Dump debug histograms
      if (multTPC>0 && acceptEvent>0){
        if (pileUp1DITS/multTPC>0.005) jet_hevents_pileUpV01D->Fill(multV0/multTPC/2.55);
        if (pileUp1DITS/multTPC>0.005) jet_hevents_pileUpITS1D->Fill(pileUp1DITS/multTPC);
        jet_hevents_pileUpV02D   ->Fill(multTPC,multV0);
        jet_hevents_pileUpITS2D  ->Fill(multTPC,pileUp1DITS);
        jet_hevents_itsLayer0V02D->Fill((*fevents_itsClustersPerLayer)[0],multV0);
        jet_hevents_ITSTPCeff    ->Fill(fITSTPCeffTrack);
        jet_hevents_V0M_CL0_After->Fill((*fevents_centrality)[0],(*fevents_centrality)[1]);
        //
        hjetConst_cent1D->Fill(fjetConst_cent);
        hjetConst_vz1D->Fill(fevents_vz);

      }
    }
    //
    // Event Cuts --> pile up and time series cut
    if (acceptEvent<0) continue;
    Double_t primMult = fevents_primMult;
    Double_t trackTgl = TMath::Abs(TMath::SinH(fjetConst_eta));
    //
    // dump some debug histogram
    /*
    if (fjetConst_dcaxy>-10 && fjetConst_dcaxy<10)   hjetConst_dcaxy2D ->Fill(fjetConst_pT,fjetConst_dcaxy);
    if (fjetConst_dcaz>-10  && fjetConst_dcaz<10)    hjetConst_dcaz2D  ->Fill(fjetConst_pT,fjetConst_dcaz);
    if (fjetConst_ncltpc>30 && fjetConst_ncltpc<170) hjetConst_ncltpc2D->Fill(fjetConst_pT,Float_t(fjetConst_ncltpc));
    */
    //

    // Retrieve cut setting
    Bool_t systCut = ApplyTreeSelection(fSystSet, fjetConst_cutBit);
    Bool_t etaAcc = (fjetConst_eta >= -0.8 && fjetConst_eta <= 0.8);
    Bool_t momAcc = (fjetConst_p >= 0.2    && fjetConst_p <= 3.2);
    //
    // dump tidentree
    //Check jet cuts, add to if
    Bool_t jet_ptsub_acc = ( (fjetConst_ptsub >= 40.0 && fjetConst_radius >= 0.15 && fjetConst_radius <= 0.25) ||
    (fjetConst_ptsub >= 60.0 && fjetConst_radius >= 0.35 && fjetConst_radius <= 0.45) ||
    (fjetConst_ptsub >= 80.0 && fjetConst_radius >= 0.55 && fjetConst_radius <= 0.65) );
    Bool_t jet_area_acc = (fjetConst_area >= 0.6*3.14159265359*fjetConst_radius*fjetConst_radius);
    if (systCut && etaAcc && momAcc && jet_ptsub_acc && jet_area_acc && fillTIdenTree){
      //
      //
      // Fill final histograms
      /*
      if ( fSystSet==21 && TMath::Abs(fjetConst_dcaxy) > 3.2 ) continue;
      if ( fSystSet==21 && TMath::Abs(fjetConst_dcaz)  > 2.4 ) continue;
      hjetConst_dcaxy1D  ->Fill(fjetConst_dcaxy);
      hjetConst_dcaz1D   ->Fill(fjetConst_dcaz);
      hjetConst_ncltpc1D ->Fill(fjetConst_ncltpc);
      */
      hjetConst_eta1D    ->Fill(fjetConst_eta);
      //hjetConst_cRows1D  ->Fill(fjetConst_cRows);
      hjetConst_phi1D    ->Fill(fjetConst_phi);
      //hjetConst_chi2tpc1D->Fill(fjetConst_chi2tpc);
      //hjetConst_dcaxy2D_After ->Fill(fjetConst_pT,fjetConst_dcaxy);
      //hjetConst_dcaz2D_After  ->Fill(fjetConst_pT,fjetConst_dcaz);

      //
      // fill trees for a given cent
      for (Int_t icent=0; icent<nCentBins; icent++){
        Bool_t centAcc = (fjetConst_cent>=centBinning[icent] && fjetConst_cent<centBinning[icent+1]);
        if ( centAcc ){
          //
          tidenTreeStream[icent]->GetFile()->cd();
          //
          // jetConst tree
          (*tidenTreeStream[icent])<<"jetConst"<<
          "gid="                  << fjetConst_gid          <<  //  global event ID
          "dEdx="                 << fjetConst_dEdx         <<  //  dEdx of the track
          "sign="                 << fjetConst_sign         <<  //  charge
          "ptot="                 << fjetConst_ptot         <<  //  TPC momentum
          "p="                    << fjetConst_p            <<  //  vertex momentum
          "pT="                   << fjetConst_pT           <<  //  vertex pT
          "eta="                  << fjetConst_eta          <<  //  eta
          "cent="                 << fjetConst_cent         <<  //  centrality
          "jetRadius="            << fjetConst_radius       <<  //  jet radius
          "jetArea="              << fjetConst_area         <<  //  jetArea
          "maxpt="                << fjetConst_maxpt        <<  //  jet max const pt
          "jetptsub="             << fjetConst_ptsub        <<  //  jetpt after subtraction
          "\n";
          continue;
        }
      }
    }
    if (systCut && jet_ptsub_acc && jet_area_acc && testIntegratedHist){
      treeStream->GetFile()->cd();
      (*treeStream)<<"jetConst"<<
      "gid="                  << fjetConst_gid          <<  //  global event ID
      "p="                    << fjetConst_p            <<  //  vertex momentum
      //"intrate="              << fjetConst_intrate      <<  //  interaction rate
      "cutBit="               << fjetConst_cutBit       <<  //  Systematic Cuts
      "dEdx="                 << fjetConst_dEdx         <<  //  dEdx of the track
      "sign="                 << fjetConst_sign         <<  //  charge
      "ptot="                 << fjetConst_ptot         <<  //  TPC momentum
      "p="                    << fjetConst_p            <<  //  vertex momentum
      "pT="                   << fjetConst_pT           <<  //  transverse momentum
      "eta="                  << fjetConst_eta          <<  //  eta
      "phi="                  << fjetConst_phi          <<  //  eta
      "cent="                 << fjetConst_cent         <<  //  centrality
      "jetRadius="            << fjetConst_radius       <<  //  jet radius
      "jetArea="              << fjetConst_area         <<  //  jetArea
      "maxpt="                << fjetConst_maxpt        <<  //  jet max const pt
      "jetptsub="             << fjetConst_ptsub        <<  //  jetpt after subtraction
      //
      "fevents_tpcVertexInfo.=" << fevents_tpcVertexInfo <<
      "fevents_itsClustersPerLayer.=" << fevents_itsClustersPerLayer <<
      "primMult="             << primMult             <<  // interaction rate
      "tgl="                  << trackTgl             <<  // interaction rate
      /*
      "cRows="                << fjetConst_cRows        <<  // interaction rate
      "chi2tpc="              << fjetConst_chi2tpc      <<  // interaction rate
      "dcaz="                 << fjetConst_dcaz         <<  // interaction rate
      "dcaxy="                << fjetConst_dcaxy        <<  // interaction rate
      */
      "\n";
    }
    //
    // apply dca and nclusters cut
    //if (TMath::Abs(fjetConst_dcaz)>3) continue;
    //if (TMath::Abs(fjetConst_dcaxy)>3) continue;
    //
    // Prepare dEdx histograms
    for (Int_t icent=0; icent<nCentBins; icent++){
      for (Int_t ieta=0; ieta<nEtaBins; ieta++){
        for (Int_t ijetRad=0; ijetRad<njetRad; ijetRad++){
        //
        if (testIntegratedHist && ieta>0) continue;
        if (!jet_ptsub_acc || !jet_area_acc) continue;
        //
        Bool_t jet_etaCentString = (fjetConst_eta>=jet_etaBinning[ieta] && fjetConst_eta<jet_etaBinning[ieta+1] && fjetConst_cent>=centBinning[icent] && fjetConst_cent<centBinning[icent+1] && fjetConst_radius>(jetRadSizes[ijetRad]-0.05) && fjetConst_radius<(jetRadSizes[ijetRad]+0.05));
        if (testIntegratedHist && ieta==0) {
          jet_etaCentString = (fjetConst_cent>=centBinning[icent] && fjetConst_cent<centBinning[icent+1] && fjetConst_radius>(jetRadSizes[ijetRad]-0.05) && fjetConst_radius<(jetRadSizes[ijetRad]+0.05));
        }
        //
        Bool_t cleanKaCutTOFTRD = ((fjetConst_cutBit >> kCleanKaTOFTRD) & 1);
        Bool_t cleanKaCutBayes = ((fjetConst_cutBit >> kTrackProbKaTOF) & 1);
        Bool_t cleanDeCutTOF = ((fjetConst_cutBit >> kCleanDeTOF) & 1);
        Bool_t parPos = fjetConst_sign>0.;
        Bool_t parNeg = fjetConst_sign<0.;
        //
        Bool_t prTOF = ((fjetConst_cutBit >> kCleanPrTOF) & 1);
        Bool_t kaTOF = ((fjetConst_cutBit >> kCleanKaTOF) & 1);
        //
        Bool_t vertexPcut=kFALSE;
        if(fCutON=="pT")   vertexPcut = (fjetConst_pT>=ptotMin      && fjetConst_pT<=ptotMax);
        if(fCutON=="ptot") vertexPcut = (fjetConst_ptot>=ptotMin    && fjetConst_ptot<=ptotMax);
        if(fCutON=="p")    vertexPcut = (fjetConst_p>=ptotMin && fjetConst_p<=ptotMax);

        if (jet_etaCentString && vertexPcut && systCut          ) jet_h2Dall[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && parPos) jet_h2Dpos[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && parNeg) jet_h2Dneg[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        //
        if (jet_etaCentString && vertexPcut && systCut && kaTOF)           jet_h2DallKaTOF[icent][ieta][ijetRad]   ->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && kaTOF && parPos) jet_h2DallKaTOFPos[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && kaTOF && parNeg) jet_h2DallKaTOFNeg[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        //
        if (jet_etaCentString && vertexPcut && systCut && prTOF)           jet_h2DallPrTOF[icent][ieta][ijetRad]   ->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && prTOF && parPos) jet_h2DallPrTOFPos[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && prTOF && parNeg) jet_h2DallPrTOFNeg[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        //
        if (jet_etaCentString && vertexPcut && systCut && prTOF)           jet_h2DallPiTOF[icent][ieta][ijetRad]   ->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && prTOF && parPos) jet_h2DallPiTOFPos[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && prTOF && parNeg) jet_h2DallPiTOFNeg[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        //
        if (jet_etaCentString && vertexPcut && systCut && !prTOF)           jet_h2DallBkgTOF[icent][ieta][ijetRad]   ->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && !prTOF && parPos) jet_h2DallBkgTOFPos[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && !prTOF && parNeg) jet_h2DallBkgTOFNeg[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        //
        if (jet_etaCentString && vertexPcut && systCut && cleanKaCutBayes)  jet_h2DCleanKaBayes[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && cleanKaCutTOFTRD) jet_h2DCleanKaTOFTRD[icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);

        if (jet_etaCentString && vertexPcut && systCut && cleanDeCutTOF)    jet_h2DClean[4][icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);
        if (jet_etaCentString && vertexPcut && systCut && cleanKaCutTOFTRD) jet_h2DClean[2][icent][ieta][ijetRad]->Fill(fjetConst_ptot,fjetConst_dEdx);

      }
      }
    }

  }  // tree loop
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessJetConstHists DONE ========= #events = " << eventCount << std::endl;

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
        fQAtimeSeriesExist = kFALSE;
      } else {
        if (fPeriod==0){
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.binMedian");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.binMedian");     // weighting wrt track
        } else {
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.mean");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.mean");     // weighting wrt track
        }
      }
    }
    //
    //  TPC ITS matching eff from time series
    Double_t fITSTPCeffEvent = (fQAtimeSeriesExist && fgrTimeSeriesEventWeighted)   ? fgrTimeSeriesEventWeighted  ->Eval(fevents_timestamp) : 1;
    Double_t fITSTPCeffTrack = (fQAtimeSeriesExist && fgrTimeSeriesNTracksWeighted) ? fgrTimeSeriesNTracksWeighted->Eval(fevents_timestamp) : 1;
    //
    // PileUp selection:  chain->SetAlias("multTPCITS","((multSSD+multSDD)/2.38)");   // chain->SetAlias("cut1000","tpcMult-multTPCITS<1000");
    Double_t multSSD = (*fevents_itsClustersPerLayer)[4]+(*fevents_itsClustersPerLayer)[5];
    Double_t multSDD = (*fevents_itsClustersPerLayer)[2]+(*fevents_itsClustersPerLayer)[3];
    Double_t multSPD = (*fevents_itsClustersPerLayer)[0]+(*fevents_itsClustersPerLayer)[1];
    Double_t multV0 = fevents_vZeroMult->Sum();
    Double_t multT0 = fevents_tZeroMult->Sum();
    Double_t multTPC = TMath::Max(fevents_tpcMult,fevents_tpcTrackBeforeClean);
    Double_t pileUp1DITS = (multSSD+multSDD)/2.38;
    //
    // dump tidentree
    treeStream->GetFile()->cd();
    (*treeStream)<<"events"<<
    "run="            << fevents_run         <<  //  global event ID
    "timeStamp="      << fevents_timestamp   <<  //  global event ID
    "bField="         << fevents_bField      <<  //  global event ID
    "gid="            << fevents_gid         <<  //  global event ID
    "vz="             << fevents_vz          <<  //  global event ID
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
    "pileupTPCCut="         << fInPileupTPCCut          <<  // interaction rate
    "tSeriesCut="           << fInTimeSeriesEff          <<  // interaction rate
    "\n";
    eventCount++;


  }  // tree loop
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessDScaledTree DONE ========= #events = " << eventCount << std::endl;

}
//____________________________________________________________________________________________________________
void ProcessJetEventTree()
{

  timer.Reset(); timer.Start();
  std::cout << " ========= ProcessJetEventTree ========= " << std::endl;
  Double_t nTreeEntriesAll = jetEventtree -> GetEntries();
  cout << " Jet Event Tree entries = " << nTreeEntriesAll << endl;
  //
  // Loop over tree entries
  Int_t eventCount=0;
  TGraph *fgrTimeSeriesEventWeighted=NULL;
  TGraph *fgrTimeSeriesNTracksWeighted=NULL;
  for(Int_t i = 0; i < nTreeEntriesAll; ++i)
  {
    //
    // Get Track and event information
    jetEventtree -> GetEntry(i);
    if(i%Int_t(nTreeEntriesAll/10) == 0) {
      cout << i << "   fjetEvents_run = " << fjetEvents_run << "    fjetEvents_primMult = " << fjetEvents_primMult << "    fjetEvents_intrate = " << fjetEvents_intrate << endl;
    }
    //
    // Get the time series map for a given run
    static Int_t runCache = -1;
    if (runCache!=Int_t(fjetEvents_run)) {
      runCache=Int_t(fjetEvents_run);
      TString timeSeriesMap = "";
      if (fPeriod==0) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2015/LHC15o/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==1) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18q/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      if (fPeriod==2) timeSeriesMap = Form("%s/lustre/nyx/alice/users/miranov/NOTESData/alice-tpc-notes/JIRA/PWGPP-538/alice/data/2018/LHC18r/pass1/000%d/QAtimeSeries.root",MapDirName.Data(),runCache);
      cout << "timeSeries Map = " <<  timeSeriesMap << endl;
      TFile *fQAtimeSeries = TFile::Open(timeSeriesMap);
      if (!fQAtimeSeries) {
        std::cout << " yolun acik ola --> file does not exist  " << std::endl;
        fQAtimeSeriesExist = kFALSE;
      } else {
        if (fPeriod==0){
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.binMedian");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.binMedian");     // weighting wrt track
        } else {
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.mean");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.mean");     // weighting wrt track
        }
      }
    }
    //
    //  TPC ITS matching eff from time series
    //Double_t fITSTPCeffEvent = (fQAtimeSeriesExist && fgrTimeSeriesEventWeighted)   ? fgrTimeSeriesEventWeighted  ->Eval(fjetEvents_timestamp) : 1;
    //Double_t fITSTPCeffTrack = (fQAtimeSeriesExist && fgrTimeSeriesNTracksWeighted) ? fgrTimeSeriesNTracksWeighted->Eval(fjetEvents_timestamp) : 1;
    //
    // PileUp selection:  chain->SetAlias("multTPCITS","((multSSD+multSDD)/2.38)");   // chain->SetAlias("cut1000","tpcMult-multTPCITS<1000");
    /*
    Double_t multSSD = (*fevents_itsClustersPerLayer)[4]+(*fevents_itsClustersPerLayer)[5];
    Double_t multSDD = (*fevents_itsClustersPerLayer)[2]+(*fevents_itsClustersPerLayer)[3];
    Double_t multSPD = (*fevents_itsClustersPerLayer)[0]+(*fevents_itsClustersPerLayer)[1];
    Double_t multV0 = fjetEvent_vZeroMult->Sum();
    Double_t multT0 = fjetEvent_tZeroMult->Sum();
    */
    Double_t multTPC = TMath::Max(fjetEvents_tpcMult,fjetEvents_tpcTrackBeforeClean);
    //Double_t pileUp1DITS = (multSSD+multSDD)/2.38;
    //
    // dump tidentree
    treeStream->GetFile()->cd();
    (*treeStream)<<"jetEvents"<<
    "run="            << fjetEvents_run         <<  //  global event ID
    "cent="            << fjetEvents_cent         <<  //  cent
    "rhoFJ="            << fjetEvents_rhoFJ         <<  //  rhoFJ
    "hasAcceptedFJjet="  << fjetEvents_accjet         <<  //  hasAcceptedFJjet
    "hasRealFJjet="       << fjetEvents_realjet         <<  //  hasRealFJjet - only events with jet pt sub > 40GeV. No other jet cuts
    //"timeStamp="      << fjetEvents_timestamp   <<  //  global event ID
    "bField="         << fjetEvents_bField      <<  //  global jetEvent ID
    "gid="            << fjetEvents_gid         <<  //  global jetEvent ID
    "vz="             << fjetEvents_vz          <<  //  global jetEvent ID
    "tpcvz="          << fjetEvents_tpcvz          <<  //  global jetEvent ID
    "spdvz="          << fjetEvents_spdvz          <<  //  global jetEvent ID
    "eventMultESD="   << fjetEvents_eventMultESD     <<  //  global event ID
    //"cent="           << (*fjetEvents_centrality)[0] <<
    "intrate="        << fjetEvents_intrate     <<  // interaction rate
    /*
    "multSSD="        << multSSD             <<  //  global event ID
    "multSDD="        << multSDD             <<  // interaction rate
    "multSPD="        << multSPD             <<  //  Systematic Cuts
    "multV0="         << multV0              <<  //  dEdx of the track
    "multT0="         << multT0              <<  //  charge
    */
    "multTPC="        << multTPC             <<  //  TPC momentum
    //"pileUp1DITS="    << pileUp1DITS         <<  //  TPC momentum
    //"ITSTPCeffEvent=" << fITSTPCeffEvent     <<  //  TPC momentum
    //"ITSTPCeffTrack=" << fITSTPCeffTrack     <<  //  TPC momentum
    "tpcTrackBeforeClean=" << fjetEvents_tpcTrackBeforeClean     <<  //  TPC momentum
    //"nTracksStored=" << fjetEvents_nTracksStored     <<  //  TPC momentum
    //
    "pileupTPCCut="         << fInPileupTPCCut          <<  // interaction rate
    "tSeriesCut="           << fInTimeSeriesEff          <<  // interaction rate
    "\n";
    eventCount++;


  }  // tree loop
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessJetEventTree DONE ========= #events = " << eventCount << std::endl;

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
    if(i%Int_t(nTreeEntries/10) == 0) cout << i << " fdscaled_cRows = " << fdscaled_cRows << "  fdscaled_ncltpc = " << fdscaled_ncltpc << "  fdscaled_vz = " << fdscaled_vz << endl;
    //
    if (fdscaled_sharedTPCClusters>0    && fdscaled_sharedTPCClusters<0.5    && systCut) hdscaled_sharedTPCClusters1D ->Fill(fdscaled_sharedTPCClusters);
    if (fdscaled_tpcSignalN>0           && fdscaled_tpcSignalN<170           && systCut) hdscaled_tpcSignalN1D        ->Fill(fdscaled_tpcSignalN);
    if (fdscaled_lengthInActiveZone>0   && fdscaled_lengthInActiveZone<170   && systCut) hdscaled_lengthInActiveZone1D->Fill(fdscaled_lengthInActiveZone);
    if (fdscaled_cRows>0                && fdscaled_cRows<170                && systCut) hdscaled_cRows1D->Fill(fdscaled_cRows);
    if (fdscaled_nclits>-2              && fdscaled_nclits<8                 && systCut) hdscaled_nclits1D->Fill(fdscaled_nclits);
    if (fdscaled_ncltrd>0               && fdscaled_ncltrd<200               && systCut) hdscaled_ncltrd1D->Fill(fdscaled_ncltrd);
    if (fdscaled_chi2tpc>0              && fdscaled_chi2tpc<10               && systCut) hdscaled_chi2tpc1D->Fill(fdscaled_chi2tpc);
    if (fdscaled_chi2its>0              && fdscaled_chi2its<20               && systCut) hdscaled_chi2its1D->Fill(fdscaled_chi2its);
    if (fdscaled_chi2trd>0              && fdscaled_chi2trd<20               && systCut) hdscaled_chi2trd1D->Fill(fdscaled_chi2trd);
    if (fdscaled_vz>-20                 && fdscaled_vz<20                    && systCut) hdscaled_vz1D->Fill(fdscaled_vz);
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

  /*

  cd /home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists
  aliroot -l
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_FilterTreesMakeHists_Run2.C+
  // .L /home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/code/RealData_FilterTreesMakeHists_Run2.C+
  TString input = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN2/LHC15o_pass1_NoSelection_06082019/mergedRuns/mergedHists/AnalysisResults_hist.root"
  ProcessExpectedHists("ptot","p",0, input)

  */

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
    for (Int_t ieta=0; ieta<nEtaBins; ieta++){

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
  Double_t nTreeEntries    = (testEntries>0) ? testEntries : nTreeEntriesAll;
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
    if(i%Int_t(nTreeEntriesAll/10) == 0) {
      cout << i << " ffArmPodTree_dEdx0 = " << ffArmPodTree_dEdx0 << "  ffArmPodTree_eta0 = " << ffArmPodTree_eta0 << "  ffArmPodTree_cent = " << ffArmPodTree_cent << "  --->  " << piFromK0 << " --  " << v0haspixel << endl;
      cout << i << " ffArmPodTree_dEdx1 = " << ffArmPodTree_dEdx1 << "  ffArmPodTree_eta1 = " << ffArmPodTree_eta1 << "  ffArmPodTree_cent = " << ffArmPodTree_cent << "  --->  " << piFromK0 << " --  " << v0haspixel << endl;
    }
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
        fQAtimeSeriesExist = kFALSE;
      } else {
        if (fPeriod==0){
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.binMedian");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.binMedian");     // weighting wrt track
        } else {
          fgrTimeSeriesEventWeighted   = (TGraph*)fQAtimeSeries -> Get("hisTimeEvEffITSDist.mean");   // weighting wrt event
          fgrTimeSeriesNTracksWeighted = (TGraph*)fQAtimeSeries -> Get("hisTimeEffITSDist.mean");     // weighting wrt track
        }
      }
    }
    //
    // Apply event selection
    static ULong64_t gidCache = -1;
    static Int_t acceptEvent = -1;
    if (gidCache!=ULong64_t(ffArmPodTree_gid)) {
      gidCache=ULong64_t(ffArmPodTree_gid);
      //
      //  TPC ITS matching eff from time series
      Double_t fITSTPCeffEvent = (fQAtimeSeriesExist && fgrTimeSeriesEventWeighted)   ? fgrTimeSeriesEventWeighted  ->Eval(fevents_timestamp) : 1;
      Double_t fITSTPCeffTrack = (fQAtimeSeriesExist && fgrTimeSeriesNTracksWeighted) ? fgrTimeSeriesNTracksWeighted->Eval(fevents_timestamp) : 1;
      //
      // PileUp selection:  chain->SetAlias("multTPCITS","((multSSD+multSDD)/2.38)");   // chain->SetAlias("cut1000","tpcMult-multTPCITS<1000");
      Double_t multSSD = (*fevents_itsClustersPerLayer)[4]+(*fevents_itsClustersPerLayer)[5];
      Double_t multSDD = (*fevents_itsClustersPerLayer)[2]+(*fevents_itsClustersPerLayer)[3];
      Double_t multTPC = TMath::Max(fevents_tpcMult,fevents_tpcTrackBeforeClean);
      Double_t pileUp1DITS = (multSSD+multSDD)/2.38;
      //
      // in bunch pileup cut
      Double_t primMultITS = (*fevents_itsVertexInfo)[2];
      Double_t secMultITS0 = (*fevents_itsVertexInfo)[3];
      Double_t secMultITS  = secMultITS0 - 0.01*primMultITS;
      //
      // trigger selections
      Bool_t bMB          = ((fevents_triggerMask&0x1)>0);
      Bool_t bCentral     = ((fevents_triggerMask&0x80000)>0);
      Bool_t bSemiCentral = ((fevents_triggerMask&0x100000)>0);
      Bool_t bMB10        = ((fevents_triggerMask&0x10)>0);
      Bool_t bAlltiggers = (bMB || bCentral || bSemiCentral || bMB10 );
      //
      // Define the cuts
      Bool_t tpcPileUpCut     = ( ftracks_cutBit & (1 << fInPileupTPCCut) );
      Bool_t timeSeriesCut    = (fITSTPCeffTrack > fInTimeSeriesEff);
      Bool_t inbunchPileUpCut = (secMultITS < (0.05*primMultITS+50));
      Bool_t spdVzCut         = (fevents_spdvz!=0 && TMath::Abs(fevents_spdvz)<7);
      //
      // pile up and time series cut
      if ( tpcPileUpCut && timeSeriesCut && inbunchPileUpCut && spdVzCut && bAlltiggers)  acceptEvent=1;
      else acceptEvent=-1;

    }
    //
    // Event Cuts --> pile up and time series cut
    if (acceptEvent<0) continue;
    Double_t primMult = fevents_primMult;
    Double_t trackTglPos  = TMath::Abs(TMath::SinH(ffArmPodTree_eta0));
    Double_t trackTglNeg  = TMath::Abs(TMath::SinH(ffArmPodTree_eta1));
    //
    // dump tidentree
    if (testIntegratedHist){
      treeStream->GetFile()->cd();
      (*treeStream)<<"clean"<<
      // "gid="                  << ffArmPodTree_gid          <<  //  global event ID
      // "p="                    << ffArmPodTree_p0            <<  //  TPC momentum
      "intrate="              << ffArmPodTree_intrate      <<  // interaction rate
      "cutBit="               << ffArmPodTree_cutBit0       <<  //  Systematic Cuts
      "cent="                 << ffArmPodTree_cent         <<  //  centrality
      "purity="               << ffArmPodTree_purity         <<  //  centrality
      "piFromK0="             << ffArmPodTree_piFromK0         <<  //  centrality
      "alfa="                 << ffArmPodTree_alfa         <<  //  centrality
      "qt="                   << ffArmPodTree_qt         <<  //  centralit
      //
      "tgl0="                 << trackTglPos             <<  // interaction rate
      "dEdx0="                << ffArmPodTree_dEdx0         <<  //  dEdx of the track
      "sign0="                << ffArmPodTree_sign0         <<  //  charge
      "ptot0="                << ffArmPodTree_ptot0         <<  //  TPC momentum
      "pT0="                  << ffArmPodTree_pT0           <<
      "p0="                   << ffArmPodTree_p0           <<
      "eta0="                 << ffArmPodTree_eta0          <<  //  eta
      "phi0="                 << ffArmPodTree_phi0          <<  //  eta
      //
      "tgl1="                 << trackTglNeg             <<  // interaction rate
      "dEdx1="                << ffArmPodTree_dEdx1         <<  //  dEdx of the track
      "sign1="                << ffArmPodTree_sign1         <<  //  charge
      "ptot1="                << ffArmPodTree_ptot1         <<  //  TPC momentum
      "pT1="                  << ffArmPodTree_pT1           <<
      "p1="                   << ffArmPodTree_p1           <<
      "eta1="                 << ffArmPodTree_eta1          <<  //  eta
      "phi1="                 << ffArmPodTree_phi1          <<  //  eta
      //
      "fevents_tpcVertexInfo.=" << fevents_tpcVertexInfo <<
      "fevents_itsClustersPerLayer.=" << fevents_itsClustersPerLayer <<
      "primMult="             << primMult             <<  // interaction rate
      "\n";
    }


    for (Int_t icent=0; icent<nCentBins; icent++){
      for (Int_t ieta=0; ieta<nEtaBins; ieta++){
        //
        if (testIntegratedHist && ieta>0) continue;
        //
        Bool_t etaCentString = (ffArmPodTree_eta0>=etaBinning[ieta] && ffArmPodTree_eta0<etaBinning[ieta+1] && ffArmPodTree_cent>=centBinning[icent] && ffArmPodTree_cent<centBinning[icent+1]);
        if (testIntegratedHist && ieta==0) {
          etaCentString = (ffArmPodTree_cent>=centBinning[icent] && ffArmPodTree_cent<centBinning[icent+1]);
        }
        //
        Bool_t vertexPcut=kFALSE;
        if(fCutON=="pT")   vertexPcut = (ffArmPodTree_pT0>=ptotMin   && ffArmPodTree_pT0<=ptotMax);
        if(fCutON=="ptot") vertexPcut = (ffArmPodTree_ptot0>=ptotMin && ffArmPodTree_ptot0<=ptotMax);
        if(fCutON=="p")    vertexPcut = (ffArmPodTree_p0>=ptotMin    && ffArmPodTree_p0<=ptotMax);
        //
        Bool_t cleanCutPiTOF0 = ( TMath::Abs(ffArmPodTree_nSigmasPiTOF0)<2. );
        Bool_t cleanCutPiTOF1 = ( TMath::Abs(ffArmPodTree_nSigmasPiTOF1)<2. );
        Bool_t cleanCutPrTOF0 = ( TMath::Abs(ffArmPodTree_nSigmasPrTOF0)<2. && (ffArmPodTree_qt>0.05) && (ffArmPodTree_qt<0.12) );
        Bool_t cleanCutPrTOF1 = ( TMath::Abs(ffArmPodTree_nSigmasPrTOF1)<2. && (ffArmPodTree_qt>0.05) && (ffArmPodTree_qt<0.12) );
        //
        Bool_t cleanCutPi    = ( TMath::Abs(ffArmPodTree_alfa)<0.5 );
        Bool_t cleanCutEl    = ( (ffArmPodTree_qt<0.005) && (TMath::Abs(ffArmPodTree_alfa)<0.5) );
        Bool_t piK0cut       = (  Bool_t(piFromK0) );
        Bool_t piPixelcut    = ( !Bool_t(v0haspixel) );
        Bool_t cleanCutPiKineCut = ( (ffArmPodTree_qt>0.14) && ((ffArmPodTree_purity >> 3) & 1) );
        Bool_t cleanCutElKineCut = ( (ffArmPodTree_qt<0.01) && ((ffArmPodTree_purity >> 5) & 1) );
        Bool_t cleanCutPrKineCut = ( (ffArmPodTree_qt<0.12) && (ffArmPodTree_qt>0.01) && ((ffArmPodTree_purity >> 4) & 1) && TMath::Abs(ffArmPodTree_nSigmasPrTOF0)<3. && TMath::Abs(ffArmPodTree_nSigmasPrTOF1)<3. );

        //
        // clean sample histograms
        if (etaCentString && vertexPcut && cleanCutEl){
          h2DClean[0][icent][ieta]->Fill(ffArmPodTree_ptot0,ffArmPodTree_dEdx0);
          h2DClean[0][icent][ieta]->Fill(ffArmPodTree_ptot1,ffArmPodTree_dEdx1);
        }
        if (etaCentString && vertexPcut && cleanCutPi && piK0cut && piPixelcut){
          h2DCleanPiTight[icent][ieta]->Fill(ffArmPodTree_ptot0,ffArmPodTree_dEdx0);
          h2DCleanPiTight[icent][ieta]->Fill(ffArmPodTree_ptot1,ffArmPodTree_dEdx1);
        }
        if (etaCentString && vertexPcut && cleanCutPi && piK0cut && cleanCutPiTOF0) h2DClean[1][icent][ieta]  ->Fill(ffArmPodTree_ptot1,ffArmPodTree_dEdx1);
        if (etaCentString && vertexPcut && cleanCutPi && piK0cut && cleanCutPiTOF1) h2DClean[1][icent][ieta]  ->Fill(ffArmPodTree_ptot0,ffArmPodTree_dEdx0);
        if (etaCentString && vertexPcut && cleanCutPi && cleanCutPiTOF0)            h2DCleanPiTOF[icent][ieta]->Fill(ffArmPodTree_ptot0,ffArmPodTree_dEdx0);
        if (etaCentString && vertexPcut && cleanCutPi && cleanCutPiTOF1)            h2DCleanPiTOF[icent][ieta]->Fill(ffArmPodTree_ptot1,ffArmPodTree_dEdx1);
        if (etaCentString && vertexPcut && cleanCutPrTOF0)                          h2DClean[3][icent][ieta]  ->Fill(ffArmPodTree_ptot0,ffArmPodTree_dEdx0);
        if (etaCentString && vertexPcut && cleanCutPrTOF1)                          h2DClean[3][icent][ieta]  ->Fill(ffArmPodTree_ptot1,ffArmPodTree_dEdx1);
        //
        // Kine cuts
        if (etaCentString && vertexPcut && cleanCutPiKineCut) { h2DCleanPiKineCut[icent][ieta]->Fill(ffArmPodTree_ptot0,ffArmPodTree_dEdx0); h2DCleanPiKineCut[icent][ieta]->Fill(ffArmPodTree_ptot1,ffArmPodTree_dEdx1);}
        if (etaCentString && vertexPcut && cleanCutElKineCut) { h2DCleanElKineCut[icent][ieta]->Fill(ffArmPodTree_ptot0,ffArmPodTree_dEdx0); h2DCleanElKineCut[icent][ieta]->Fill(ffArmPodTree_ptot1,ffArmPodTree_dEdx1);}
        if (etaCentString && vertexPcut && cleanCutPrKineCut) { h2DCleanPrKineCut[icent][ieta]->Fill(ffArmPodTree_ptot0,ffArmPodTree_dEdx0); h2DCleanPrKineCut[icent][ieta]->Fill(ffArmPodTree_ptot1,ffArmPodTree_dEdx1);}

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
  TString centEtaStr = Form("cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);
  //
  // get each particle PID response
  for (Int_t ipart=0;ipart<nParticles;ipart++){
    fhnExpected->GetAxis(0)->SetRangeUser(ipart,ipart+1);
    h2Expected[ipart][icent][ieta]      = (TH2F*)fhnExpected->Projection(6,4);
    h2ExpectedSigma[ipart][icent][ieta] = (TH2F*)fhnExpected->Projection(5,4);
    grExpected[ipart][icent][ieta]      = ProfileToGraphErrors(h2Expected[ipart][icent][ieta]);
    grExpectedSigma[ipart][icent][ieta] = ProfileToGraphErrors(h2ExpectedSigma[ipart][icent][ieta]);
    h2Expected[ipart][icent][ieta]     ->SetName(Form("h2Expected_%d_%s"     ,ipart,centEtaStr.Data()));
    h2ExpectedSigma[ipart][icent][ieta]->SetName(Form("h2ExpectedSigma_%d_%s",ipart,centEtaStr.Data()));
    grExpected[ipart][icent][ieta]     ->SetName(Form("grExpected_%d_%s"     ,ipart,centEtaStr.Data()));
    grExpectedSigma[ipart][icent][ieta]->SetName(Form("grExpectedSigma_%d_%s",ipart,centEtaStr.Data()));

  }

}
//____________________________________________________________________________________________________________
void InitInitials()
{

  TString outputFileNameTree     = Form("Trees_Syst%d_PlotVS_%s_CutON_%s_%d_%3.2f.root"    ,fSystSet,fPlotVS.Data(),fCutON.Data(),fInPileupTPCCut,fInTimeSeriesEff);
  TString outputFileNameHist     = Form("Hists_Syst%d_PlotVS_%s_CutON_%s_%d_%3.2f.root"    ,fSystSet,fPlotVS.Data(),fCutON.Data(),fInPileupTPCCut,fInTimeSeriesEff);
  TString jetOutputFileNameHist     = Form("Jet_Hists_Syst%d_PlotVS_%s_CutON_%s_%d_%3.2f.root"    ,fSystSet,fPlotVS.Data(),fCutON.Data(),fInPileupTPCCut,fInTimeSeriesEff);
  TString debugFile              = Form("Debug_Syst%d_PlotVS_%s_CutON_%s_%d_%3.2f.root"    ,fSystSet,fPlotVS.Data(),fCutON.Data(),fInPileupTPCCut,fInTimeSeriesEff);
  treeStream      = new TTreeSRedirector(outputFileNameTree,"recreate");
  histStream      = new TTreeSRedirector(outputFileNameHist,"recreate");
  jetHistStream   = new TTreeSRedirector(jetOutputFileNameHist,"recreate");
  debugStream     = new TTreeSRedirector(debugFile,"recreate");
  //
  //
  if (fillTIdenTree){
    for (Int_t icent=0; icent<nCentBins; icent++){
      TString idenStrName = Form("TIdenTree_Syst%d_cent_%3.2f_%3.2f_PlotVS_%s_CutON_%s_%d_%3.2f.root",fSystSet,centBinning[icent],centBinning[icent+1],fPlotVS.Data(),fCutON.Data(),fInPileupTPCCut,fInTimeSeriesEff);
      tidenTreeStream[icent] = new TTreeSRedirector(idenStrName,"recreate");
    }
  }
  //
  // Initialise histograms to be used for binning
  //
  hevents_pileUpV02D               = new TH2F("hevents_pileUpV02D"                ,"hevents_pileUpV02D"              ,1000,0.,45000., 1000 ,0., 45000. );
  hevents_pileUpV02D_inbunchCut    = new TH2F("hevents_pileUpV02D_inbunchCut"     ,"hevents_pileUpV02D_inbunchCut"   ,1000,0.,45000., 1000 ,0., 45000. );
  hevents_pileUpITS2D              = new TH2F("hevents_pileUpITS2D"               ,"hevents_pileUpITS2D"                ,1500,0.,45000., 1500 ,0., 45000. );
  hevents_pileUpITS2D_inbunchCut   = new TH2F("hevents_pileUpITS2D_inbunchCut"    ,"hevents_pileUpITS2D_inbunchCut"     ,1000,0.,45000., 1000 ,0., 45000. );
  hevents_itsLayer0V02D            = new TH2F("hevents_itsLayer0V02D"             ,"hevents_itsLayer0V02D"              ,1000,0.,45000., 1000 ,0., 45000. );
  hevents_itsLayer0V02D_inbunchCut = new TH2F("hevents_itsLayer0V02D_inbunchCut"  ,"hevents_itsLayer0V02D_inbunchCut"   ,1000,0.,45000., 1000 ,0., 45000. );
  //
  hevents_secMultITS0_primMultITS  = new TH2F("hevents_secMultITS0_primMultITS"   ,"hevents_secMultITS0_primMultITS"    ,250,0.,5000., 200 ,0., 200. );
  hevents_V0M_CL0                  = new TH2F("hevents_V0M_CL0"                   ,"hevents_V0M_CL0"                    ,100,0.,100., 100 ,0., 100. );
  hevents_V0M_CL0_After            = new TH2F("hevents_V0M_CL0_After"             ,"hevents_V0M_CL0_After"              ,100,0.,100., 100 ,0., 100. );
  hevents_pileUpV01D               = new TH1F("hevents_pileUpV01D"                ,"hevents_pileUpV01D"                 ,1000,0.05,1.4);
  hevents_pileUpITS1D              = new TH1F("hevents_pileUpITS1D"               ,"hevents_pileUpITS1D"                ,1000,0.05,1.4);
  hevents_ITSTPCeff                = new TH1F("hevents_ITSTPCeff"                 ,"hevents_ITSTPCeff"                  ,1000,0.,1.2);
  //
  // Jet event histos
  jet_hevents_pileUpV02D               = new TH2F("jet_hevents_pileUpV02D"                ,"jet_hevents_pileUpV02D"              ,1000,0.,45000., 1000 ,0., 45000. );
  jet_hevents_pileUpV02D_inbunchCut    = new TH2F("jet_hevents_pileUpV02D_inbunchCut"     ,"jet_hevents_pileUpV02D_inbunchCut"   ,1000,0.,45000., 1000 ,0., 45000. );
  jet_hevents_pileUpITS2D              = new TH2F("jet_hevents_pileUpITS2D"               ,"jet_hevents_pileUpITS2D"                ,1500,0.,45000., 1500 ,0., 45000. );
  jet_hevents_pileUpITS2D_inbunchCut   = new TH2F("jet_hevents_pileUpITS2D_inbunchCut"    ,"jet_hevents_pileUpITS2D_inbunchCut"     ,1000,0.,45000., 1000 ,0., 45000. );
  jet_hevents_itsLayer0V02D            = new TH2F("jet_hevents_itsLayer0V02D"             ,"jet_hevents_itsLayer0V02D"              ,1000,0.,45000., 1000 ,0., 45000. );
  jet_hevents_itsLayer0V02D_inbunchCut = new TH2F("jet_hevents_itsLayer0V02D_inbunchCut"  ,"jet_hevents_itsLayer0V02D_inbunchCut"   ,1000,0.,45000., 1000 ,0., 45000. );
  //
  jet_hevents_secMultITS0_primMultITS  = new TH2F("jet_hevents_secMultITS0_primMultITS"   ,"jet_hevents_secMultITS0_primMultITS"    ,250,0.,5000., 200 ,0., 200. );
  jet_hevents_V0M_CL0                  = new TH2F("jet_hevents_V0M_CL0"                   ,"jet_hevents_V0M_CL0"                    ,100,0.,100., 100 ,0., 100. );
  jet_hevents_V0M_CL0_After            = new TH2F("jet_hevents_V0M_CL0_After"             ,"jet_hevents_V0M_CL0_After"              ,100,0.,100., 100 ,0., 100. );
  jet_hevents_pileUpV01D               = new TH1F("jet_hevents_pileUpV01D"                ,"jet_hevents_pileUpV01D"                 ,1000,0.05,1.4);
  jet_hevents_pileUpITS1D              = new TH1F("jet_hevents_pileUpITS1D"               ,"jet_hevents_pileUpITS1D"                ,1000,0.05,1.4);
  jet_hevents_ITSTPCeff                = new TH1F("jet_hevents_ITSTPCeff"                 ,"jet_hevents_ITSTPCeff"                  ,1000,0.,1.2);
  //
  htracks_dcaxy2D             = new TH2F("htracks_dcaxy2D"                     ,"htracks_dcaxy2D"  ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  htracks_dcaz2D              = new TH2F("htracks_dcaz2D"                      ,"htracks_dcaz2D"   ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  htracks_ncltpc2D            = new TH2F("htracks_ncltpc2D"                    ,"htracks_ncltpc2D" ,ptNbins,ptotMin,ptotMax, 140 , 30., 170.);
  htracks_dcaxy2D_After       = new TH2F("htracks_dcaxy2D_After"               ,"htracks_dcaxy2D_After"  ,ptNbins,ptotMin,ptotMax, 1600 ,-4., 4. );
  htracks_dcaz2D_After        = new TH2F("htracks_dcaz2D_After"                ,"htracks_dcaz2D_After"   ,ptNbins,ptotMin,ptotMax, 1600 ,-4., 4. );
  htracks_dcaxy1D             = new TH1F("htracks_dcaxy1D"                     ,"htracks_dcaxy1D"   ,1600 ,-4., 4. );
  htracks_dcaz1D              = new TH1F("htracks_dcaz1D"                      ,"htracks_dcaz1D"    ,1600 ,-4., 4. );
  htracks_ncltpc1D            = new TH1F("htracks_ncltpc1D"                    ,"htracks_ncltpc1D " ,140 , 30., 170.);
  htracks_eta1D               = new TH1F("htracks_eta1D"                       ,"htracks_eta1D"     ,100 , -1., 1.);
  htracks_cRows1D             = new TH1F("htracks_cRows1D"                     ,"htracks_cRows1D"   ,170 , 0., 170.);
  htracks_cent1D              = new TH1F("htracks_cent1D"                      ,"htracks_cent1D"    ,200 , 0., 100.);
  htracks_phi1D               = new TH1F("htracks_phi1D"                       ,"htracks_phi1D"     ,200 , 4., 4.);
  htracks_chi2tpc1D           = new TH1F("htracks_chi2tpc1D"                   ,"htracks_chi2tpc1D" ,200 , -1., 10.);
  htracks_vz1D                = new TH1F("htracks_vz1D"                        ,"htracks_vz1D"      ,200 ,-20., 20. );
  //
  /*
  hjetConst_dcaxy2D             = new TH2F("hjetConst_dcaxy2D"                     ,"hjetConst_dcaxy2D"  ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  hjetConst_dcaz2D              = new TH2F("hjetConst_dcaz2D"                      ,"hjetConst_dcaz2D"   ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  hjetConst_ncltpc2D            = new TH2F("hjetConst_ncltpc2D"                    ,"hjetConst_ncltpc2D" ,ptNbins,ptotMin,ptotMax, 140 , 30., 170.);
  hjetConst_dcaxy2D_After       = new TH2F("hjetConst_dcaxy2D_After"               ,"hjetConst_dcaxy2D_After"  ,ptNbins,ptotMin,ptotMax, 1600 ,-4., 4. );
  hjetConst_dcaz2D_After        = new TH2F("hjetConst_dcaz2D_After"                ,"hjetConst_dcaz2D_After"   ,ptNbins,ptotMin,ptotMax, 1600 ,-4., 4. );
  hjetConst_dcaxy1D             = new TH1F("hjetConst_dcaxy1D"                     ,"hjetConst_dcaxy1D"   ,1600 ,-4., 4. );
  hjetConst_dcaz1D              = new TH1F("hjetConst_dcaz1D"                      ,"hjetConst_dcaz1D"    ,1600 ,-4., 4. );
  hjetConst_ncltpc1D            = new TH1F("hjetConst_ncltpc1D"                    ,"hjetConst_ncltpc1D " ,140 , 30., 170.);
  */
  hjetConst_eta1D               = new TH1F("hjetConst_eta1D"                       ,"hjetConst_eta1D"     ,100 , -1., 1.);
  //hjetConst_cRows1D             = new TH1F("hjetConst_cRows1D"                     ,"hjetConst_cRows1D"   ,170 , 0., 170.);
  hjetConst_cent1D              = new TH1F("hjetConst_cent1D"                      ,"hjetConst_cent1D"    ,200 , 0., 100.);
  hjetConst_phi1D               = new TH1F("hjetConst_phi1D"                       ,"hjetConst_phi1D"     ,200 , 4., 4.);
  //hjetConst_chi2tpc1D           = new TH1F("hjetConst_chi2tpc1D"                   ,"hjetConst_chi2tpc1D" ,200 , -1., 10.);
  hjetConst_vz1D                = new TH1F("hjetConst_vz1D"                        ,"hjetConst_vz1D"      ,200 ,-20., 20. );
  //
  hdscaled_sharedTPCClusters1D  = new TH1F("hdscaled_sharedTPCClusters1D" ,"hdscaled_sharedTPCClusters1D Bins" ,200   ,  0., 0.5 );
  hdscaled_tpcSignalN1D         = new TH1F("hdscaled_tpcSignalN1D"        ,"hdscaled_tpcSignalN1D Bins"        ,170   ,  0., 170. );
  hdscaled_lengthInActiveZone1D = new TH1F("hdscaled_lengthInActiveZone1D","hdscaled_lengthInActiveZone1D Bins",170   ,  0., 170. );
  hdscaled_cRows1D              = new TH1F("hdscaled_cRows1D"             ,"hdscaled_cRows1D Bins"             ,170   ,  0., 170. );
  hdscaled_nclits1D             = new TH1F("hdscaled_nclits1D"            ,"hdscaled_nclits1D Bins"            ,10    , -2., 8. );
  hdscaled_ncltrd1D             = new TH1F("hdscaled_ncltrd1D"            ,"hdscaled_ncltrd1D Bins"            ,200   ,  0., 200. );
  hdscaled_chi2tpc1D            = new TH1F("hdscaled_chi2tpc1D"           ,"hdscaled_chi2tpc1D Bins"           ,400   ,  0., 10. );
  hdscaled_chi2its1D            = new TH1F("hdscaled_chi2its1D"           ,"hdscaled_chi2its1D Bins"           ,400   ,  0., 20. );
  hdscaled_chi2trd1D            = new TH1F("hdscaled_chi2trd1D"           ,"hdscaled_chi2trd1D Bins"           ,400   ,  0., 20. );
  hdscaled_vz1D                 = new TH1F("hdscaled_vz1D"                ,"hdscaled_vz1D Bins"                ,200   ,-20., 20. );

  //
  hhighPt_dEdxPtot = new TH2F("hhighPt_dEdxPtot", "hhighPt_dEdxPtot" ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
  //
  for (Int_t icent=0; icent<nCentBins; icent++){
    for (Int_t ieta=0; ieta<nEtaBins; ieta++){

      TString centEtaStr = Form("cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centBinning[icent],centBinning[icent+1],etaBinning[ieta],etaBinning[ieta+1]);

      h2Dall[icent][ieta]=NULL;   h2Dpos[icent][ieta]=NULL;   h2Dneg[icent][ieta]=NULL;
      hName2Dall[icent][ieta]=Form("h2Dall_%s",centEtaStr.Data());
      hName2Dpos[icent][ieta]=Form("h2Dpos_%s",centEtaStr.Data());
      hName2Dneg[icent][ieta]=Form("h2Dneg_%s",centEtaStr.Data());
      h2Dall[icent][ieta] = new TH2F(hName2Dall[icent][ieta], hName2Dall[icent][ieta] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2Dpos[icent][ieta] = new TH2F(hName2Dpos[icent][ieta], hName2Dpos[icent][ieta] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2Dneg[icent][ieta] = new TH2F(hName2Dneg[icent][ieta], hName2Dneg[icent][ieta] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DallPrTOF[icent][ieta]=NULL;      h2DallPrTOFPos[icent][ieta]=NULL;      h2DallPrTOFNeg[icent][ieta]=NULL;
      hName2DallPrTOF[icent][ieta]   =Form("h2DallPrTOF_%s"   ,centEtaStr.Data());
      hName2DallPrTOFPos[icent][ieta]=Form("h2DallPrTOFPos_%s",centEtaStr.Data());
      hName2DallPrTOFNeg[icent][ieta]=Form("h2DallPrTOFNeg_%s",centEtaStr.Data());
      h2DallPrTOF[icent][ieta]     = new TH2F(hName2DallPrTOF[icent][ieta]    ,hName2DallPrTOF[icent][ieta]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallPrTOFPos[icent][ieta]  = new TH2F(hName2DallPrTOFPos[icent][ieta] ,hName2DallPrTOFPos[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallPrTOFNeg[icent][ieta]  = new TH2F(hName2DallPrTOFNeg[icent][ieta] ,hName2DallPrTOFNeg[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DallBkgTOF[icent][ieta]=NULL;      h2DallBkgTOFPos[icent][ieta]=NULL;      h2DallBkgTOFNeg[icent][ieta]=NULL;
      hName2DallBkgTOF[icent][ieta]   =Form("h2DallBkgTOF_%s"   ,centEtaStr.Data());
      hName2DallBkgTOFPos[icent][ieta]=Form("h2DallBkgTOFPos_%s",centEtaStr.Data());
      hName2DallBkgTOFNeg[icent][ieta]=Form("h2DallBkgTOFNeg_%s",centEtaStr.Data());
      h2DallBkgTOF[icent][ieta]     = new TH2F(hName2DallBkgTOF[icent][ieta]    ,hName2DallBkgTOF[icent][ieta]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallBkgTOFPos[icent][ieta]  = new TH2F(hName2DallBkgTOFPos[icent][ieta] ,hName2DallBkgTOFPos[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallBkgTOFNeg[icent][ieta]  = new TH2F(hName2DallBkgTOFNeg[icent][ieta] ,hName2DallBkgTOFNeg[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DallKaTOF[icent][ieta]=NULL;      h2DallKaTOFPos[icent][ieta]=NULL;      h2DallKaTOFNeg[icent][ieta]=NULL;
      hName2DallKaTOF[icent][ieta]   =Form("h2DallKaTOF_%s"   ,centEtaStr.Data());
      hName2DallKaTOFPos[icent][ieta]=Form("h2DallKaTOFPos_%s",centEtaStr.Data());
      hName2DallKaTOFNeg[icent][ieta]=Form("h2DallKaTOFNeg_%s",centEtaStr.Data());
      h2DallKaTOF[icent][ieta]     = new TH2F(hName2DallKaTOF[icent][ieta]    ,hName2DallKaTOF[icent][ieta]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallKaTOFPos[icent][ieta]  = new TH2F(hName2DallKaTOFPos[icent][ieta] ,hName2DallKaTOFPos[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallKaTOFNeg[icent][ieta]  = new TH2F(hName2DallKaTOFNeg[icent][ieta] ,hName2DallKaTOFNeg[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DCleanPiKineCut[icent][ieta]=NULL;      h2DCleanPrKineCut[icent][ieta]=NULL;      h2DCleanElKineCut[icent][ieta]=NULL;
      hName2DCleanPiKineCut[icent][ieta]   =Form("h2DCleanPi_KineCut_%s" ,centEtaStr.Data());
      hName2DCleanPrKineCut[icent][ieta]   =Form("h2DCleanPr_KineCut_%s" ,centEtaStr.Data());
      hName2DCleanElKineCut[icent][ieta]   =Form("h2DCleanEl_KineCut_%s" ,centEtaStr.Data());
      h2DCleanPiKineCut[icent][ieta]  = new TH2F(hName2DCleanPiKineCut[icent][ieta] ,hName2DCleanPiKineCut[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanPrKineCut[icent][ieta]  = new TH2F(hName2DCleanPrKineCut[icent][ieta] ,hName2DCleanPrKineCut[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanElKineCut[icent][ieta]  = new TH2F(hName2DCleanElKineCut[icent][ieta] ,hName2DCleanElKineCut[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      h2DCleanKaTOFTRD[icent][ieta]=NULL;  h2DCleanKaBayes[icent][ieta]=NULL;  h2DCleanPiTight[icent][ieta]=NULL;  h2DCleanPiTOF[icent][ieta]=NULL;
      hName2DCleanKaTOFTRD[icent][ieta] = Form("h2DCleanKa_TOFTRD_%s"   ,centEtaStr.Data());
      hName2DCleanKaBayes[icent][ieta] = Form("h2DCleanKa_Bayes_%s"   ,centEtaStr.Data());
      hName2DCleanPiTight[icent][ieta] = Form("h2DCleanPi_Tight_%s"   ,centEtaStr.Data());
      hName2DCleanPiTOF[icent][ieta] = Form("h2DCleanPi_TOF_%s"   ,centEtaStr.Data());
      h2DCleanKaTOFTRD[icent][ieta]     = new TH2F(hName2DCleanKaTOFTRD[icent][ieta],hName2DCleanKaTOFTRD[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanKaBayes[icent][ieta]     = new TH2F(hName2DCleanKaBayes[icent][ieta],hName2DCleanKaBayes[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanPiTight[icent][ieta]     = new TH2F(hName2DCleanPiTight[icent][ieta],hName2DCleanPiTight[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DCleanPiTOF[icent][ieta]     = new TH2F(hName2DCleanPiTOF[icent][ieta],hName2DCleanPiTOF[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //

      //Jet consts
      for (Int_t ijetRad=0; ijetRad<njetRad; ijetRad++){
            TString jet_centEtaStr = Form("cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f_jetRad_%3.2f",centBinning[icent],centBinning[icent+1],jet_etaBinning[ieta],jet_etaBinning[ieta+1],jetRadSizes[ijetRad]);

            jet_h2Dall[icent][ieta][ijetRad]=NULL;   jet_h2Dpos[icent][ieta][ijetRad]=NULL;   jet_h2Dneg[icent][ieta][ijetRad]=NULL;
            jet_hName2Dall[icent][ieta][ijetRad]=Form("jet_h2Dall_%s",jet_centEtaStr.Data());
            jet_hName2Dpos[icent][ieta][ijetRad]=Form("jet_h2Dpos_%s",jet_centEtaStr.Data());
            jet_hName2Dneg[icent][ieta][ijetRad]=Form("jet_h2Dneg_%s",jet_centEtaStr.Data());
            jet_h2Dall[icent][ieta][ijetRad] = new TH2F(jet_hName2Dall[icent][ieta][ijetRad], jet_hName2Dall[icent][ieta][ijetRad] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2Dpos[icent][ieta][ijetRad] = new TH2F(jet_hName2Dpos[icent][ieta][ijetRad], jet_hName2Dpos[icent][ieta][ijetRad] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2Dneg[icent][ieta][ijetRad] = new TH2F(jet_hName2Dneg[icent][ieta][ijetRad], jet_hName2Dneg[icent][ieta][ijetRad] ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            //
            jet_h2DallPrTOF[icent][ieta][ijetRad]=NULL;      jet_h2DallPrTOFPos[icent][ieta][ijetRad]=NULL;      jet_h2DallPrTOFNeg[icent][ieta][ijetRad]=NULL;
            jet_hName2DallPrTOF[icent][ieta][ijetRad]   =Form("jet_h2DallPrTOF_%s"   ,jet_centEtaStr.Data());
            jet_hName2DallPrTOFPos[icent][ieta][ijetRad]=Form("jet_h2DallPrTOFPos_%s",jet_centEtaStr.Data());
            jet_hName2DallPrTOFNeg[icent][ieta][ijetRad]=Form("jet_h2DallPrTOFNeg_%s",jet_centEtaStr.Data());
            jet_h2DallPrTOF[icent][ieta][ijetRad]     = new TH2F(jet_hName2DallPrTOF[icent][ieta][ijetRad]    ,jet_hName2DallPrTOF[icent][ieta][ijetRad]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2DallPrTOFPos[icent][ieta][ijetRad]  = new TH2F(jet_hName2DallPrTOFPos[icent][ieta][ijetRad] ,jet_hName2DallPrTOFPos[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2DallPrTOFNeg[icent][ieta][ijetRad]  = new TH2F(jet_hName2DallPrTOFNeg[icent][ieta][ijetRad] ,jet_hName2DallPrTOFNeg[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            //
            jet_h2DallBkgTOF[icent][ieta][ijetRad]=NULL;      jet_h2DallBkgTOFPos[icent][ieta][ijetRad]=NULL;      jet_h2DallBkgTOFNeg[icent][ieta][ijetRad]=NULL;
            jet_hName2DallBkgTOF[icent][ieta][ijetRad]   =Form("jet_h2DallBkgTOF_%s"   ,jet_centEtaStr.Data());
            jet_hName2DallBkgTOFPos[icent][ieta][ijetRad]=Form("jet_h2DallBkgTOFPos_%s",jet_centEtaStr.Data());
            jet_hName2DallBkgTOFNeg[icent][ieta][ijetRad]=Form("jet_h2DallBkgTOFNeg_%s",jet_centEtaStr.Data());
            jet_h2DallBkgTOF[icent][ieta][ijetRad]     = new TH2F(jet_hName2DallBkgTOF[icent][ieta][ijetRad]    ,jet_hName2DallBkgTOF[icent][ieta][ijetRad]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2DallBkgTOFPos[icent][ieta][ijetRad]  = new TH2F(jet_hName2DallBkgTOFPos[icent][ieta][ijetRad] ,jet_hName2DallBkgTOFPos[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2DallBkgTOFNeg[icent][ieta][ijetRad]  = new TH2F(jet_hName2DallBkgTOFNeg[icent][ieta][ijetRad] ,jet_hName2DallBkgTOFNeg[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            //
            jet_h2DallKaTOF[icent][ieta][ijetRad]=NULL;      jet_h2DallKaTOFPos[icent][ieta][ijetRad]=NULL;      jet_h2DallKaTOFNeg[icent][ieta][ijetRad]=NULL;
            jet_hName2DallKaTOF[icent][ieta][ijetRad]   =Form("jet_h2DallKaTOF_%s"   ,jet_centEtaStr.Data());
            jet_hName2DallKaTOFPos[icent][ieta][ijetRad]=Form("jet_h2DallKaTOFPos_%s",jet_centEtaStr.Data());
            jet_hName2DallKaTOFNeg[icent][ieta][ijetRad]=Form("jet_h2DallKaTOFNeg_%s",jet_centEtaStr.Data());
            jet_h2DallKaTOF[icent][ieta][ijetRad]     = new TH2F(jet_hName2DallKaTOF[icent][ieta][ijetRad]    ,jet_hName2DallKaTOF[icent][ieta][ijetRad]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2DallKaTOFPos[icent][ieta][ijetRad]  = new TH2F(jet_hName2DallKaTOFPos[icent][ieta][ijetRad] ,jet_hName2DallKaTOFPos[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2DallKaTOFNeg[icent][ieta][ijetRad]  = new TH2F(jet_hName2DallKaTOFNeg[icent][ieta][ijetRad] ,jet_hName2DallKaTOFNeg[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            //
            jet_h2DallPiTOF[icent][ieta][ijetRad]=NULL;      jet_h2DallPiTOFPos[icent][ieta][ijetRad]=NULL;      jet_h2DallPiTOFNeg[icent][ieta][ijetRad]=NULL;
            jet_hName2DallPiTOF[icent][ieta][ijetRad]   =Form("jet_h2DallPiTOF_%s"   ,jet_centEtaStr.Data());
            jet_hName2DallPiTOFPos[icent][ieta][ijetRad]=Form("jet_h2DallPiTOFPos_%s",jet_centEtaStr.Data());
            jet_hName2DallPiTOFNeg[icent][ieta][ijetRad]=Form("jet_h2DallPiTOFNeg_%s",jet_centEtaStr.Data());
            jet_h2DallPiTOF[icent][ieta][ijetRad]     = new TH2F(jet_hName2DallPiTOF[icent][ieta][ijetRad]    ,jet_hName2DallPiTOF[icent][ieta][ijetRad]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2DallPiTOFPos[icent][ieta][ijetRad]  = new TH2F(jet_hName2DallPiTOFPos[icent][ieta][ijetRad] ,jet_hName2DallPiTOFPos[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2DallPiTOFNeg[icent][ieta][ijetRad]  = new TH2F(jet_hName2DallPiTOFNeg[icent][ieta][ijetRad] ,jet_hName2DallPiTOFNeg[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            //
            jet_h2DCleanKaTOFTRD[icent][ieta][ijetRad]=NULL;  jet_h2DCleanKaBayes[icent][ieta][ijetRad]=NULL;
            jet_hName2DCleanKaTOFTRD[icent][ieta][ijetRad] = Form("jet_h2DCleanKa_TOFTRD_%s"   ,jet_centEtaStr.Data());
            jet_hName2DCleanKaBayes[icent][ieta][ijetRad] = Form("jet_h2DCleanKa_Bayes_%s"   ,jet_centEtaStr.Data());
            jet_h2DCleanKaTOFTRD[icent][ieta][ijetRad]     = new TH2F(jet_hName2DCleanKaTOFTRD[icent][ieta][ijetRad],jet_hName2DCleanKaTOFTRD[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            jet_h2DCleanKaBayes[icent][ieta][ijetRad]     = new TH2F(jet_hName2DCleanKaBayes[icent][ieta][ijetRad],jet_hName2DCleanKaBayes[icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            //
            for (Int_t ipart=0; ipart<nParticles; ipart++)
            {
              jet_h2DClean[ipart][icent][ieta][ijetRad]=NULL;
              jet_hName2DClean[ipart][icent][ieta][ijetRad] = Form("jet_h2DClean%s_%s", parName[ipart].Data(), jet_centEtaStr.Data());
              jet_h2DClean[ipart][icent][ieta][ijetRad]     = new TH2F(jet_hName2DClean[ipart][icent][ieta][ijetRad],jet_hName2DClean[ipart][icent][ieta][ijetRad]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
            }
      }
      for (Int_t ipart=0; ipart<nParticles; ipart++)
      {
        h2DClean[ipart][icent][ieta]=NULL;
        hName2DClean[ipart][icent][ieta] = Form("h2DClean%s_%s", parName[ipart].Data(), centEtaStr.Data());
        h2DClean[ipart][icent][ieta]     = new TH2F(hName2DClean[ipart][icent][ieta],hName2DClean[ipart][icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
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
  syst:
  0 -->  Reference
  1 -->  kNCrossedRowsTPC70
  2 -->  kNCrossedRowsTPC90
  3 -->  kActiveZone
  4 -->  kMaxChi2PerClusterTPCSmall
  5 -->  kMaxChi2PerClusterTPCLarge
  6 -->  kMaxDCAToVertexXYPtDepLarge
  7 -->  kVertexZSmall
  8 -->  kEventVertexZLarge
  9 -->  kSharedClsLoose
  10 --> kFindableClsLoose
  11 --> kFindableClsLoosest
  12 --> kPileupLoose
  13 --> kBFieldPos
  13 --> kBFieldNeg
  13 --> kTPCSignalNSmall
  13 --> kTPCSignalNLarge
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
    case kCutActiveZone:  // 3 -->  kActiveZone
    {
      fCutArr = {kActiveZone,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCSmall:   // 4 -->  kMaxChi2PerClusterTPCSmall
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCSmall, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCLarge:   // 5 -->  kMaxChi2PerClusterTPCLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCLarge, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutMaxDCAToVertexXYPtDepLarge:   // 6 -->  kMaxDCAToVertexXYPtDepLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDepLarge, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutVertexZSmall:   // 7 -->  kVertexZSmall
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZSmall, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutEventVertexZLarge:  // 8 -->  kEventVertexZLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZLarge, kPileup, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutSharedCls:   // 9 -->  kSharedClsLoose
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedClsLoose, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutFindableClsTight:   // 10 -->  kFindableClsTight
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableClsTight,kTPCSignalN};
    }
    break;
    //
    case kCutFindableClsLoose:   // 11 -->  kFindableClsLoose
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableClsLoose,kTPCSignalN};
    }
    break;
    //
    case kCutPileupLoose:   // 12 -->  kPileupLoose
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileupLoose, kSharedCls, kFindableCls,kTPCSignalN};
    }
    break;
    //
    case kCutBFieldPos:   // 13 -->  kBFieldPos
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN,kBFieldPos};
    }
    break;
    //
    case kCutBFieldNeg:   // 14 --> kBFieldNeg
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalN,kBFieldNeg};
    }
    break;
    //
    case kCutTPCSignalNSmall:   // 15 --> kTPCSignalNSmall
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalNSmall};
    }
    break;
    //
    case kCutTPCSignalNLarge:   // 16 --> kTPCSignalNLarge
    {
      fCutArr = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPileup, kSharedCls, kFindableCls,kTPCSignalNLarge};
    }
    break;
    //
    default:
    {
      fCutArr = {0};
    }

  }
  //
  //
  //  Apply conditions
  Bool_t cutAll = kTRUE;
  for (Int_t cut : fCutArr) {
    if ( (cutBit & (1 << cut)) == kFALSE ) {
      cutAll = kFALSE;
      break;
    }
  }
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
      //
      if (testIntegratedHist && ieta>0) continue;
      //
      if(h2Dall[icent][ieta])  h2Dall[icent][ieta]->Write();
      if(h2Dpos[icent][ieta])  h2Dpos[icent][ieta]->Write();
      if(h2Dneg[icent][ieta])  h2Dneg[icent][ieta]->Write();
      //
      if(h2DallPrTOF[icent][ieta])     h2DallPrTOF[icent][ieta]->Write();
      if(h2DallPrTOFPos[icent][ieta])  h2DallPrTOFPos[icent][ieta]->Write();
      if(h2DallPrTOFNeg[icent][ieta])  h2DallPrTOFNeg[icent][ieta]->Write();
      //
      if(h2DallBkgTOF[icent][ieta])     h2DallBkgTOF[icent][ieta]->Write();
      if(h2DallBkgTOFPos[icent][ieta])  h2DallBkgTOFPos[icent][ieta]->Write();
      if(h2DallBkgTOFNeg[icent][ieta])  h2DallBkgTOFNeg[icent][ieta]->Write();
      //
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
  //
  jetHistStream->GetFile()->cd();
  for (Int_t icent=0; icent<nCentBins; icent++){
    for (Int_t ieta=0; ieta<nEtaBins; ieta++){
      for (Int_t ijetRad=0; ijetRad<njetRad; ijetRad++){
      //
      if (testIntegratedHist && ieta>0) continue;
      //
      if(jet_h2Dall[icent][ieta][ijetRad])  jet_h2Dall[icent][ieta][ijetRad]->Write();
      if(jet_h2Dpos[icent][ieta][ijetRad])  jet_h2Dpos[icent][ieta][ijetRad]->Write();
      if(jet_h2Dneg[icent][ieta][ijetRad])  jet_h2Dneg[icent][ieta][ijetRad]->Write();
      //
      if(jet_h2DallPrTOF[icent][ieta][ijetRad])     jet_h2DallPrTOF[icent][ieta][ijetRad]->Write();
      if(jet_h2DallPrTOFPos[icent][ieta][ijetRad])  jet_h2DallPrTOFPos[icent][ieta][ijetRad]->Write();
      if(jet_h2DallPrTOFNeg[icent][ieta][ijetRad])  jet_h2DallPrTOFNeg[icent][ieta][ijetRad]->Write();
      //
      if(jet_h2DallBkgTOF[icent][ieta][ijetRad])     jet_h2DallBkgTOF[icent][ieta][ijetRad]->Write();
      if(jet_h2DallBkgTOFPos[icent][ieta][ijetRad])  jet_h2DallBkgTOFPos[icent][ieta][ijetRad]->Write();
      if(jet_h2DallBkgTOFNeg[icent][ieta][ijetRad])  jet_h2DallBkgTOFNeg[icent][ieta][ijetRad]->Write();
      //
      if(jet_h2DallKaTOF[icent][ieta][ijetRad])     jet_h2DallKaTOF[icent][ieta][ijetRad]->Write();
      if(jet_h2DallKaTOFPos[icent][ieta][ijetRad])  jet_h2DallKaTOFPos[icent][ieta][ijetRad]->Write();
      if(jet_h2DallKaTOFNeg[icent][ieta][ijetRad])  jet_h2DallKaTOFNeg[icent][ieta][ijetRad]->Write();
      //
      if(jet_h2DallPiTOF[icent][ieta][ijetRad])     jet_h2DallPiTOF[icent][ieta][ijetRad]->Write();
      if(jet_h2DallPiTOFPos[icent][ieta][ijetRad])  jet_h2DallPiTOFPos[icent][ieta][ijetRad]->Write();
      if(jet_h2DallPiTOFNeg[icent][ieta][ijetRad])  jet_h2DallPiTOFNeg[icent][ieta][ijetRad]->Write();
      //
      if(jet_h2DCleanKaTOFTRD[icent][ieta][ijetRad]) jet_h2DCleanKaTOFTRD[icent][ieta][ijetRad]->Write();
      if(jet_h2DCleanKaBayes[icent][ieta][ijetRad])  jet_h2DCleanKaBayes[icent][ieta][ijetRad]->Write();
      //
      /* //EMPTY Becuase we never run ProcessExpectedHists in inclusive case or here.
      for (Int_t ipart=0; ipart<nParticles; ipart++){
        if(jet_h2Expected[ipart][icent][ieta][ijetRad])       jet_h2Expected[ipart][icent][ieta][ijetRad]->Write();
        if(jet_h2ExpectedSigma[ipart][icent][ieta][ijetRad])  jet_h2ExpectedSigma[ipart][icent][ieta][ijetRad]->Write();
        if(jet_grExpected[ipart][icent][ieta][ijetRad])       jet_grExpected[ipart][icent][ieta][ijetRad]->Write();
        if(jet_grExpectedSigma[ipart][icent][ieta][ijetRad])  jet_grExpectedSigma[ipart][icent][ieta][ijetRad]->Write();
        if(jet_h2DClean[ipart][icent][ieta][ijetRad])         jet_h2DClean[ipart][icent][ieta][ijetRad]->Write(); //EMPTY BECAUSE WE DON'T HAVE CLEAN SAMPLES
      }
      */
    }
    }
  }
  //
  debugStream->GetFile()->cd();
  if(hevents_secMultITS0_primMultITS)  hevents_secMultITS0_primMultITS  ->Write();
  //
  if(hevents_itsLayer0V02D)         hevents_itsLayer0V02D  ->Write();
  if(hevents_itsLayer0V02D_inbunchCut) hevents_itsLayer0V02D_inbunchCut  ->Write();
  if(hevents_pileUpITS2D)           hevents_pileUpITS2D  ->Write();
  if(hevents_pileUpITS2D_inbunchCut)   hevents_pileUpITS2D_inbunchCut  ->Write();
  if(hevents_pileUpV02D)            hevents_pileUpV02D  ->Write();
  if(hevents_pileUpV02D_inbunchCut) hevents_pileUpV02D_inbunchCut  ->Write();
  //
  if(hevents_V0M_CL0)               hevents_V0M_CL0  ->Write();
  if(hevents_V0M_CL0_After)         hevents_V0M_CL0_After  ->Write();
  if(hevents_pileUpV01D)            hevents_pileUpV01D  ->Write();
  if(hevents_pileUpITS1D)           hevents_pileUpITS1D  ->Write();
  if(hevents_ITSTPCeff)             hevents_ITSTPCeff  ->Write();
  //Jet events
  if(jet_hevents_secMultITS0_primMultITS)  jet_hevents_secMultITS0_primMultITS  ->Write();
  //
  if(jet_hevents_itsLayer0V02D)         jet_hevents_itsLayer0V02D  ->Write();
  if(jet_hevents_itsLayer0V02D_inbunchCut) jet_hevents_itsLayer0V02D_inbunchCut  ->Write();
  if(jet_hevents_pileUpITS2D)           jet_hevents_pileUpITS2D  ->Write();
  if(jet_hevents_pileUpITS2D_inbunchCut)   jet_hevents_pileUpITS2D_inbunchCut  ->Write();
  if(jet_hevents_pileUpV02D)            jet_hevents_pileUpV02D  ->Write();
  if(jet_hevents_pileUpV02D_inbunchCut) jet_hevents_pileUpV02D_inbunchCut  ->Write();
  //
  if(jet_hevents_V0M_CL0)               jet_hevents_V0M_CL0  ->Write();
  if(jet_hevents_V0M_CL0_After)         jet_hevents_V0M_CL0_After  ->Write();
  if(jet_hevents_pileUpV01D)            jet_hevents_pileUpV01D  ->Write();
  if(jet_hevents_pileUpITS1D)           jet_hevents_pileUpITS1D  ->Write();
  if(jet_hevents_ITSTPCeff)             jet_hevents_ITSTPCeff  ->Write();

  if(hdscaled_sharedTPCClusters1D)  hdscaled_sharedTPCClusters1D ->Write();
  if(hdscaled_tpcSignalN1D)         hdscaled_tpcSignalN1D ->Write();
  if(hdscaled_lengthInActiveZone1D) hdscaled_lengthInActiveZone1D ->Write();
  if(hdscaled_cRows1D)              hdscaled_cRows1D ->Write();
  if(hdscaled_nclits1D)             hdscaled_nclits1D ->Write();
  if(hdscaled_ncltrd1D)             hdscaled_ncltrd1D ->Write();
  if(hdscaled_chi2tpc1D)            hdscaled_chi2tpc1D ->Write();
  if(hdscaled_chi2its1D)            hdscaled_chi2its1D ->Write();
  if(hdscaled_chi2trd1D)            hdscaled_chi2trd1D ->Write();
  if(hdscaled_vz1D)                 hdscaled_vz1D      ->Write();
  if(hhighPt_dEdxPtot)              hhighPt_dEdxPtot->Write();
  if(htracks_dcaxy2D)               htracks_dcaxy2D  ->Write();
  if(htracks_dcaz2D)                htracks_dcaz2D   ->Write();
  if(htracks_ncltpc2D)              htracks_ncltpc2D ->Write();
  if(htracks_dcaxy2D_After)         htracks_dcaxy2D_After->Write();
  if(htracks_dcaz2D_After)          htracks_dcaz2D_After->Write();
  //
  if(htracks_dcaxy1D)               htracks_dcaxy1D  ->Write();
  if(htracks_dcaz1D)                htracks_dcaz1D   ->Write();
  if(htracks_ncltpc1D)              htracks_ncltpc1D ->Write();
  if(htracks_eta1D)                 htracks_eta1D ->Write();
  if(htracks_cRows1D)               htracks_cRows1D ->Write();
  if(htracks_cent1D)                htracks_cent1D ->Write();
  if(htracks_phi1D)                 htracks_phi1D ->Write();
  if(htracks_chi2tpc1D)             htracks_chi2tpc1D ->Write();
  if(htracks_vz1D)                  htracks_vz1D ->Write();
  // jets
  /*
  if(hjetConst_dcaxy1D)               hjetConst_dcaxy1D  ->Write();
  if(hjetConst_dcaz1D)                hjetConst_dcaz1D   ->Write();
  if(hjetConst_ncltpc1D)              hjetConst_ncltpc1D ->Write();
  */
  if(hjetConst_eta1D)                 hjetConst_eta1D ->Write();
  //if(hjetConst_cRows1D)               hjetConst_cRows1D ->Write();
  if(hjetConst_cent1D)                hjetConst_cent1D ->Write();
  if(hjetConst_phi1D)                 hjetConst_phi1D ->Write();
  //if(hjetConst_chi2tpc1D)             hjetConst_chi2tpc1D ->Write();
  if(hjetConst_vz1D)                  hjetConst_vz1D ->Write();

}
//____________________________________________________________________________________________________________
void SetBranchAddresses()
{

  //
  // Clean Samples Tree
  if (armtree){
    armtree->SetBranchAddress("cutBit1"       ,&ffArmPodTree_cutBit1);
    armtree->SetBranchAddress("dEdx1"         ,&ffArmPodTree_dEdx1);
    armtree->SetBranchAddress("sign1"         ,&ffArmPodTree_sign1);
    armtree->SetBranchAddress("ptot1"         ,&ffArmPodTree_ptot1);
    armtree->SetBranchAddress("p1"            ,&ffArmPodTree_p1);
    armtree->SetBranchAddress("pT1"           ,&ffArmPodTree_pT1);
    armtree->SetBranchAddress("eta1"          ,&ffArmPodTree_eta1);
    armtree->SetBranchAddress("phi1"          ,&ffArmPodTree_phi1);
    armtree->SetBranchAddress("nSigmasPiTOF1" ,&ffArmPodTree_nSigmasPiTOF1);
    armtree->SetBranchAddress("nSigmasPrTOF1" ,&ffArmPodTree_nSigmasPrTOF1);
    //
    armtree->SetBranchAddress("cutBit0"       ,&ffArmPodTree_cutBit0);
    armtree->SetBranchAddress("dEdx0"         ,&ffArmPodTree_dEdx0);
    armtree->SetBranchAddress("sign0"         ,&ffArmPodTree_sign0);
    armtree->SetBranchAddress("ptot0"         ,&ffArmPodTree_ptot0);
    armtree->SetBranchAddress("p0"            ,&ffArmPodTree_p0);
    armtree->SetBranchAddress("pT0"           ,&ffArmPodTree_pT0);
    armtree->SetBranchAddress("eta0"          ,&ffArmPodTree_eta0);
    armtree->SetBranchAddress("phi0"          ,&ffArmPodTree_phi0);
    armtree->SetBranchAddress("nSigmasPiTOF0" ,&ffArmPodTree_nSigmasPiTOF0);
    armtree->SetBranchAddress("nSigmasPrTOF0" ,&ffArmPodTree_nSigmasPrTOF0);
    //
    armtree->SetBranchAddress("gid"        ,&ffArmPodTree_gid);
    armtree->SetBranchAddress("eventtime"  ,&ffArmPodTree_eventtime);
    armtree->SetBranchAddress("intrate"    ,&ffArmPodTree_intrate);
    armtree->SetBranchAddress("purity"     ,&ffArmPodTree_purity);
    armtree->SetBranchAddress("cent"       ,&ffArmPodTree_cent);
    armtree->SetBranchAddress("qt"         ,&ffArmPodTree_qt);
    armtree->SetBranchAddress("alfa"       ,&ffArmPodTree_alfa);
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
    dataTree->SetBranchAddress("phi"      ,&ftracks_phi);
    dataTree->SetBranchAddress("cent"     ,&ftracks_cent);
    dataTree->SetBranchAddress("dcaxy"    ,&ftracks_dcaxy);
    dataTree->SetBranchAddress("dcaz"     ,&ftracks_dcaz);
    dataTree->SetBranchAddress("ncltpc"   ,&ftracks_ncltpc);
    dataTree->SetBranchAddress("cRows"    ,&ftracks_cRows);
    dataTree->SetBranchAddress("chi2tpc"  ,&ftracks_chi2tpc);
  }
  // Jet constituents Tree
  if (jetConstTree){
    jetConstTree->SetBranchAddress("jetRadius"      ,&fjetConst_radius);
    jetConstTree->SetBranchAddress("jetArea"      ,&fjetConst_area);
    jetConstTree->SetBranchAddress("maxpt"      ,&fjetConst_maxpt);
    jetConstTree->SetBranchAddress("jetptsub"      ,&fjetConst_ptsub);
    jetConstTree->SetBranchAddress("gid"      ,&fjetConst_gid);
    //jetConstTree->SetBranchAddress("eventtime",&fjetConst_eventtime);
    //jetConstTree->SetBranchAddress("intrate"  ,&fjetConst_intrate);
    jetConstTree->SetBranchAddress("cutBit"   ,&fjetConst_cutBit);
    jetConstTree->SetBranchAddress("dEdx"     ,&fjetConst_dEdx);
    jetConstTree->SetBranchAddress("sign"     ,&fjetConst_sign);
    jetConstTree->SetBranchAddress("ptot"     ,&fjetConst_ptot);
    jetConstTree->SetBranchAddress("p"        ,&fjetConst_p);
    jetConstTree->SetBranchAddress("pT"       ,&fjetConst_pT);
    jetConstTree->SetBranchAddress("eta"      ,&fjetConst_eta);
    jetConstTree->SetBranchAddress("phi"      ,&fjetConst_phi);
    jetConstTree->SetBranchAddress("cent"     ,&fjetConst_cent);
    //jetConstTree->SetBranchAddress("dcaxy"    ,&fjetConst_dcaxy);
    //jetConstTree->SetBranchAddress("dcaz"     ,&fjetConst_dcaz);
    //jetConstTree->SetBranchAddress("ncltpc"   ,&fjetConst_ncltpc);
    //jetConstTree->SetBranchAddress("cRows"    ,&fjetConst_cRows);
    //jetConstTree->SetBranchAddress("chi2tpc"  ,&fjetConst_chi2tpc);
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
    dscaltree->SetBranchAddress("nclTPC"      ,&fdscaled_ncltpc);
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
    dscaltree->SetBranchAddress("nclits"      ,&fdscaled_nclits);
    dscaltree->SetBranchAddress("ncltrd"      ,&fdscaled_ncltrd);
    dscaltree->SetBranchAddress("chi2tpc"     ,&fdscaled_chi2tpc);
    dscaltree->SetBranchAddress("chi2its"     ,&fdscaled_chi2its);
    dscaltree->SetBranchAddress("chi2trd"     ,&fdscaled_chi2trd);
  }
  //
  // event tree
  if (eventtree){
    eventtree->SetBranchAddress("run"                 ,&fevents_run);
    eventtree->SetBranchAddress("intrate"             ,&fevents_intrate);
    eventtree->SetBranchAddress("bField"              ,&fevents_bField);
    eventtree->SetBranchAddress("gid"                 ,&fevents_gid);
    eventtree->SetBranchAddress("timestamp"           ,&fevents_timestamp);
    eventtree->SetBranchAddress("triggerMask"         ,&fevents_triggerMask);
    eventtree->SetBranchAddress("vz"                  ,&fevents_vz);
    eventtree->SetBranchAddress("tpcvz"               ,&fevents_tpcvz);
    eventtree->SetBranchAddress("spdvz"               ,&fevents_spdvz);
    eventtree->SetBranchAddress("tpcMult"             ,&fevents_tpcMult);
    eventtree->SetBranchAddress("eventMult"           ,&fevents_eventMult);
    eventtree->SetBranchAddress("eventMultESD"        ,&fevents_eventMultESD);
    eventtree->SetBranchAddress("nTracksStored"       ,&fevents_nTracksStored);
    eventtree->SetBranchAddress("primMult"            ,&fevents_primMult);
    eventtree->SetBranchAddress("tpcTrackBeforeClean" ,&fevents_tpcTrackBeforeClean);
    eventtree->SetBranchAddress("itsTracklets"        ,&fevents_itsTracklets);
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
    eventtree->SetBranchAddress("itsVertexInfo."      ,&fevents_itsVertexInfo);

    eventtree->SetAlias("nPileUpPrim"  ,"(tpcVertexInfo.fElements[3]+tpcVertexInfo.fElements[4])");
  }
  // jet event tree
  if (jetEventtree){
    jetEventtree->SetBranchAddress("run"                 ,&fjetEvents_run);
    jetEventtree->SetBranchAddress("cent"                 ,&fjetEvents_cent);
    jetEventtree->SetBranchAddress("rhoFJ"                 ,&fjetEvents_rhoFJ);
    jetEventtree->SetBranchAddress("hasAcceptedFJjet"       ,&fjetEvents_accjet);
    jetEventtree->SetBranchAddress("hasRealFJjet"            ,&fjetEvents_realjet);
    jetEventtree->SetBranchAddress("intrate"             ,&fjetEvents_intrate);
    jetEventtree->SetBranchAddress("bField"              ,&fjetEvents_bField);
    jetEventtree->SetBranchAddress("gid"                 ,&fjetEvents_gid);
    //jetEventtree->SetBranchAddress("timestamp"           ,&fjetEvents_timestamp);
    jetEventtree->SetBranchAddress("triggerMask"         ,&fjetEvents_triggerMask);
    jetEventtree->SetBranchAddress("vz"                  ,&fjetEvents_vz);
    jetEventtree->SetBranchAddress("tpcvz"               ,&fjetEvents_tpcvz);
    jetEventtree->SetBranchAddress("spdvz"               ,&fjetEvents_spdvz);
    jetEventtree->SetBranchAddress("tpcMult"             ,&fjetEvents_tpcMult);
    jetEventtree->SetBranchAddress("eventMult"           ,&fjetEvents_eventMult);
    jetEventtree->SetBranchAddress("eventMultESD"        ,&fjetEvents_eventMultESD);
    //jetEventtree->SetBranchAddress("nTracksStored"       ,&fjetEvents_nTracksStored);
    jetEventtree->SetBranchAddress("primMult"            ,&fjetEvents_primMult);
    jetEventtree->SetBranchAddress("tpcTrackBeforeClean" ,&fjetEvents_tpcTrackBeforeClean);
    jetEventtree->SetBranchAddress("itsTracklets"        ,&fjetEvents_itsTracklets);
    /*
    jetEventtree->SetBranchAddress("centrality."         ,&fjetEvents_centrality);
    jetEventtree->SetBranchAddress("tZeroMult."          ,&fjetEvents_tZeroMult);
    jetEventtree->SetBranchAddress("vZeroMult."          ,&fjetEvents_vZeroMult);
    jetEventtree->SetBranchAddress("itsClustersPerLayer.",&fjetEvents_itsClustersPerLayer);
    jetEventtree->SetBranchAddress("trackCounters."      ,&fjetEvents_trackCounters);
    jetEventtree->SetBranchAddress("trackdEdxRatio."     ,&fjetEvents_trackdEdxRatio);
    jetEventtree->SetBranchAddress("trackNcl."           ,&fjetEvents_trackNcl);
    jetEventtree->SetBranchAddress("trackChi2."          ,&fjetEvents_trackChi2);
    jetEventtree->SetBranchAddress("trackMatchEff."      ,&fjetEvents_trackMatchEff);
    jetEventtree->SetBranchAddress("trackTPCCountersZ."  ,&fjetEvents_trackTPCCountersZ);
    jetEventtree->SetBranchAddress("phiCountA."          ,&fjetEvents_phiCountA);
    jetEventtree->SetBranchAddress("phiCountC."          ,&fjetEvents_phiCountC);
    jetEventtree->SetBranchAddress("phiCountAITS."       ,&fjetEvents_phiCountAITS);
    jetEventtree->SetBranchAddress("phiCountCITS."       ,&fjetEvents_phiCountCITS);
    jetEventtree->SetBranchAddress("phiCountAITSOnly."   ,&fjetEvents_phiCountAITSOnly);
    jetEventtree->SetBranchAddress("phiCountCITSOnly."   ,&fjetEvents_phiCountCITSOnly);
    jetEventtree->SetBranchAddress("tpcVertexInfo."      ,&fjetEvents_tpcVertexInfo);
    jetEventtree->SetBranchAddress("itsVertexInfo."      ,&fjetEvents_itsVertexInfo);

    jetEventtree->SetAlias("nPileUpPrim"  ,"(tpcVertexInfo.fElements[3]+tpcVertexInfo.fElements[4])");
    */
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

}
//____________________________________________________________________________________________________________
void PlotTimeSeriesPerSector(TString performanceEventList)
{


  /*


  .L /home/marsland/Desktop/ubuntu_desktop/workdir/RUN_ON_GRID/Ebye/code/RealData_FilterTreesMakeHists_Run2.C+
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
Double_t TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
    return TMath::Min((0.972865 + 0.0117093 * pTmc), 1.);
}
//____________________________________________________________________________________________________________
Double_t TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
    return (1 - 0.129758 *TMath::Exp(-pTmc*0.679612));
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
