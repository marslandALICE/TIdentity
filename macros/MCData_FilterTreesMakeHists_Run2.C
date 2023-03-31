/*

meld /home/marsland/Desktop/ubuntu_desktop/workdir/code/MCData_FilterTreesMakeHists.C /u/marsland/PHD/macros/marsland_EbyeRatios/MCData_FilterTreesMakeHists.C

TFile f("AnalysisResults_trees1.root")
mcFull->Show(0)
mcFull->Draw("isample","syst==8&&etaUp>0.7&&orig==0&&pDown<0.5")
mcFull->Draw("netPiMomRec.fElements[0]/netPiMomGen.fElements[0]","syst==8&&etaUp>0.7")
mcFull->Draw("vZ:syst","etaUp>0.7&&orig==0&&pDown<0.5","*")

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
void PrepareEffMatrix(TString matrixFile, Int_t detectorCom, Int_t systSetting, Float_t pDownAcc, Float_t pUpAcc);
Bool_t ApplyTreeSelection(Int_t syst, UInt_t cutBit);
TGraphErrors * ProfileToGraphErrors(TH2F * h2);
void SetBranchAddresses();
void effMatrix1D(Int_t effType,Int_t pid, Int_t cent);
void effMatrix2D(Int_t pid, Int_t cent);
void MeanEffInAcceptance(Int_t pid);
void ModifyEffMatrix(THnF *hnF);

// ======= Modification part =======
const Int_t colors[] = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
const Int_t nParticles = 5;
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
// Acceptance arrays
const Int_t nCentAcc = 9;
const Int_t nEtaAcc = 8;
Float_t etaAccDown[nEtaAcc] = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1};
Float_t etaAccUp[nEtaAcc]   = { 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
Float_t centAccDown[nCentAcc] = {0.,5.0,10.,20.,30.,40.,50.,60.,70.};
Float_t centAccUp[nCentAcc]   = {5.,10.,20.,30.,40.,50.,60.,70.,80.};
//
// Efficiency  matrix
const Int_t nSigns = 2;
const Int_t nCentBins = 9;
const Int_t nParticlesEffMatrix = 3;
THnF *posRec, *posGen, *negRec, *negGen;
TH2F *h2EtaMom[nSigns][nParticlesEffMatrix][nCentBins];
TH2F *h2CentEta[nSigns][nParticlesEffMatrix];
//
TH1F *h1CentEta[nSigns][nParticlesEffMatrix][nEtaAcc];
TH1F *h1Mom[nSigns][nParticlesEffMatrix][nCentBins];
TH1F *h1Eta[nSigns][nParticlesEffMatrix][nCentBins];
TH1F *h1Cent[nSigns][nParticlesEffMatrix];
TH1F *h1EffAccScan[nSigns][nParticlesEffMatrix][nCentBins];
//
TTreeSRedirector *tidenTreeStream[nCentBins];
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
TH2F *h2DCleanPos[nParticles][nCentBins][nEtaBins];
TH2F *h2DCleanNeg[nParticles][nCentBins][nEtaBins];
TH2F *h2DClean[nParticles][nCentBins][nEtaBins];
//
TString hName2DCleanPos[nParticles][nCentBins][nEtaBins];
TString hName2DCleanNeg[nParticles][nCentBins][nEtaBins];
TString hName2DClean[nParticles][nCentBins][nEtaBins];

//
// ----------------------------------------------------------------------------------------------------------------------
//
// Data tree branches
UInt_t ffTreeMC_cutBit=0;
ULong64_t ffTreeMC_gid=0;
Float_t ffTreeMC_dEdx=0;
Int_t ffTreeMC_sign=0;
Int_t ffTreeMC_isample=0;
Int_t ffTreeMC_origin=0;
Int_t ffTreeMC_part=0;
Float_t ffTreeMC_ptot=0;
Float_t ffTreeMC_p=0;
Float_t ffTreeMC_pT=0;
Float_t ffTreeMC_eta=0;
Float_t ffTreeMC_cent=0;
Float_t ffTreeMC_centimp=0;
Double_t ffTreeMC_vZ=0;
Int_t ffTreeMC_ncltpc=0;
Float_t ffTreeMC_dcaxy=0;
Float_t ffTreeMC_dcaz=0;
Float_t  ffTreeMC_cRows=0;
Float_t ffTreeMC_chi2tpc=0;
Float_t ffTreeMC_phi=0;
Float_t ffTreeMC_nsigmatofka=0;
Float_t ffTreeMC_nsigmatofpr=0;
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
TVectorF *fevents_itsVertexInfo=0x0;
//
//
const Int_t fnCutBins=10;
Int_t fCutArr[fnCutBins];
Int_t fSystSet=0;
Int_t fTrackOrigin=0;
TString fPlotVS="";
TString fCutON="";
TString inputMCTree    = "fTreeMC";
TString inputHists     = "cleanHists";
TString inputHistsFile = "";
TFile *fMCdata=NULL;
TFile *fMChist=NULL;
TTree *mctree=NULL;
TTree *eventtree=NULL;
TTreeSRedirector *histStream=0, *debugStream=0, *effMatrixStream=0;;
TH1D  *hEta=0x0, *hCent=0x0, *hMom=0x0;
TStopwatch timer;
AliNDLocalRegression *hdEdxAShifttMNTglDist_meanGFitAll=NULL;
AliNDLocalRegression *hdEdxAShifttMNTglDist_meanGFitNoPileUp=NULL;
//
TH2F *hfTreeMC_dcaxy2D=NULL;
TH2F *hfTreeMC_dcaz2D=NULL;
TH2F *hfTreeMC_ncltpc2D=NULL;
TH2F *hfTreeMC_dcaxy2D_After=NULL;
TH2F *hfTreeMC_dcaz2D_After=NULL;
TH2F *hfTreeMC_ncltpc2D_After=NULL;
//
TH1F *htracks_dcaxy1D=NULL;
TH1F *htracks_dcaz1D=NULL;
TH1F *htracks_ncltpc1D=NULL;
TH1F *htracks_eta1D=NULL;
TH1F *htracks_cRows1D=NULL;
TH1F *htracks_cent1D=NULL;
TH1F *htracks_phi1D=NULL;
TH1F *htracks_chi2tpc1D=NULL;
TH1F *htracks_vz1D=NULL;
//
enum trackCutBit {
  kNCrossedRowsTPC60=0,
  kNCrossedRowsTPC80=1,
  kNCrossedRowsTPC100=2,
  kMaxChi2PerClusterTPCSmall=3,
  kMaxChi2PerClusterTPC=4,
  kMaxChi2PerClusterTPCLarge=5,
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
  kCleanPrTOF=24,
  kCleanKaTOF=25,
  kCleanKaTOFTRD=26,
  kTrackProbKaTOF=27,
  kTrackProbPrTOF=28,
  kCleanDeTOF=29,
  kEventVertexZALICE=30,
  kEventVertexZALICETight=31,
};

enum cutSettings {
  kCutReference=0,
  kCutCrossedRowsTPC60=1,
  kCutCrossedRowsTPC100=2,
  kCutMaxChi2PerClusterTPCSmall=3,
  kCutMaxChi2PerClusterTPCLarge=4,
  kCutMaxDCAToVertexXYPtDepSmall=5,
  kCutMaxDCAToVertexXYPtDepLarge=6,
  kCutVertexZSmall=7,
  kCutVertexZLarge=8,
  kCutEventVertexZSmall=9,
  kCutEventVertexZLarge=10,
  kCutRequireITSRefit=11,
  kCutPixelRequirementITS=12,
  kCutNewITSCut=13,
  kCutTPCSignalN=14,
  kCutActiveZone=15,
  kCutTPCSignalNActiveZone=16,
  kCutTPCSignalNSmallActiveZoneSmall=17,
  kCutTPCSignalNLargeActiveZoneLarge=18,
  kCutEventVertexZALICE=19,
  kCutEventVertexZALICETight=20,
};


TString parName[] = {"El", "Pi", "Ka", "Pr", "De"};
TString effTypeName[] = {"mom", "eta"};
enum parType
{
  kEl=0,
  kPi=1,
  kKa=2,
  kPr=3,
  kDe=4
};

enum binType
{
  kDet=0,
  kSyst=1,
  kPart=2,
  kCent=3,
  kMom=4,
  kEta=5
};


//
// Double_t testEntries  = 2000000;
Double_t testEntries  = -1.;
Bool_t testIntegratedHist=kFALSE;
Bool_t fillTIdenTree = kTRUE;

void MCData_FilterTreesMakeHists( TString plotVS, TString cutON, TString mcFile, TString histFile, Int_t systSet, Int_t trackOrigin)
{
  //
  // Produce all hists form MC and Real data
  // plotVS: "ptot", "pT", "p"
  //

  /*

  meld /home/marsland/Desktop/ubuntu_desktop/workdir/code/MCData_FilterTreesMakeHists.C /u/marsland/PHD/macros/marsland_EbyeRatios/MCData_FilterTreesMakeHists.C


  cd /home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists
  aliroot -l
  TString histFile = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN2/LHC16g1/LHC16g1_pass1_NoSelection_06082019/mergedRuns/mergedHists/AnalysisResults_hists.root";
  TString mcFile = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN2/LHC16g1/LHC16g1_pass1_NoSelection_06082019/mergedRuns/mergedData/runs/000246087/Sub_0/AnalysisResults_fTreeMC.root";
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/fourthMomentAnalysis/code/MCData_FilterTreesMakeHists.C+
  MCData_FilterTreesMakeHists("ptot","p",mcFile,histFile,0,0)


/media/marsland/Samsung_T5/lustre/nyx/alice/users/marsland/workdir/FlowPaper/data/PbPb2018/LHC20e3a_pass3/mergedData/AnalysisResults1.root


  aliroot -l
  TString histFile = "/media/marsland/Samsung_T5/lustre/nyx/alice/users/marsland/workdir/FlowPaper/data/PbPb2018/LHC20e3a_pass3/mergedData/AnalysisResults_hists.root";
  TString mcFile = "/media/marsland/Samsung_T5/lustre/nyx/alice/users/marsland/workdir/FlowPaper/data/PbPb2018/LHC20e3a_pass3/mergedData/AnalysisResults_fTreeMC1.root";
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/fourthMomentAnalysis/code/MCData_FilterTreesMakeHists.C+
  MCData_FilterTreesMakeHists("ptot","p",mcFile,histFile,0,0)


  */
  //
  fSystSet         = systSet;
  fPlotVS          = plotVS;
  fCutON           = cutON;
  fTrackOrigin     = trackOrigin; // 1 --> only primaries; 0 --> all particles which survive cuts
  //
  fMCdata = TFile::Open(mcFile);
  if (histFile != "") fMChist = TFile::Open(histFile);
  mctree  = (TTree*)fMCdata->Get(inputMCTree);
  //
  SetBranchAddresses();
  InitInitials();
  ProcessDataHists();
  WriteHistsToFile();
  //
  if (debugStream)     delete debugStream;
  if (histStream)      delete histStream;
  if (effMatrixStream) delete effMatrixStream;
  //
  for (Int_t icent=0; icent<nCentBins; icent++) delete tidenTreeStream[icent];

}
//____________________________________________________________________________________________________________
void ProcessDataHists()
{

  timer.Reset(); timer.Start();
  std::cout << " ========= ProcessDataHists ========= " << std::endl;
  Double_t nTreeEntriesAll = mctree -> GetEntries();
  Double_t nTreeEntries    = (testEntries>0) ? testEntries : nTreeEntriesAll;
  cout << " Data Tree entries = " << nTreeEntriesAll << endl;
  if (nTreeEntriesAll<10) { std::cout << " === upss data tree is empty === " << std::endl; return; }
  //
  // Loop over tree entries
  Int_t eventCount=0;
  for(Int_t i = 0; i < nTreeEntries; ++i)
  {
    //
    // Get Track and event information
    mctree -> GetEntry(i);
    if(i%Int_t(nTreeEntriesAll/10) == 0) {
      cout << i << "   ffTreeMC_dEdx = " << ffTreeMC_dEdx      << "         ffTreeMC_eta        = " << ffTreeMC_eta << "  ffTreeMC_dcaxy = " << ffTreeMC_dcaxy << endl;
    }
    //
    // Apply event selection
    static ULong64_t gidCache = -1;
    static Int_t acceptEvent = -1;
    if (gidCache!=ULong64_t(ffTreeMC_gid)) {
      gidCache=ULong64_t(ffTreeMC_gid);
      //
      // pile up and time series cut
      if ( 1 ) { acceptEvent=1; eventCount++;}
      else acceptEvent=-1;
    }
    //
    // Event Cuts --> pile up and time series cut
    if (acceptEvent<0) continue;
    //
    // dump some debug histogram
    if (ffTreeMC_dcaxy>-10 && ffTreeMC_dcaxy<10)   hfTreeMC_dcaxy2D ->Fill(ffTreeMC_pT,ffTreeMC_dcaxy);
    if (ffTreeMC_dcaz>-10  && ffTreeMC_dcaz<10)    hfTreeMC_dcaz2D  ->Fill(ffTreeMC_pT,ffTreeMC_dcaz);
    if (ffTreeMC_ncltpc>30 && ffTreeMC_ncltpc<170) hfTreeMC_ncltpc2D->Fill(ffTreeMC_pT,Float_t(ffTreeMC_ncltpc));
    //
    // Retrieve cut setting
    Bool_t systCut = ApplyTreeSelection(fSystSet, ffTreeMC_cutBit);
    Bool_t etaAcc = (ffTreeMC_eta >= -0.8 && ffTreeMC_eta <= 0.8);
    Bool_t momAcc = (ffTreeMC_p >= 0.2    && ffTreeMC_p <= 3.2);
    //
    // dump tidentree
    if (systCut && etaAcc && momAcc){
      //
      //
      // Fill final histograms
      if ( fSystSet==21 && TMath::Abs(ffTreeMC_dcaxy) > 3.2 ) continue;
      if ( fSystSet==21 && TMath::Abs(ffTreeMC_dcaz)  > 2.4 ) continue;
      htracks_dcaxy1D  ->Fill(ffTreeMC_dcaxy);
      htracks_dcaz1D   ->Fill(ffTreeMC_dcaz);
      htracks_ncltpc1D ->Fill(ffTreeMC_ncltpc);
      htracks_eta1D    ->Fill(ffTreeMC_eta);
      htracks_cRows1D  ->Fill(ffTreeMC_cRows);
      htracks_phi1D    ->Fill(ffTreeMC_phi);
      htracks_cent1D    ->Fill(ffTreeMC_cent);
      htracks_chi2tpc1D->Fill(ffTreeMC_chi2tpc);
      htracks_vz1D     ->Fill(ffTreeMC_vZ);
      //
      hfTreeMC_dcaxy2D_After ->Fill(ffTreeMC_pT,ffTreeMC_dcaxy);
      hfTreeMC_dcaz2D_After  ->Fill(ffTreeMC_pT,ffTreeMC_dcaz);
      hfTreeMC_ncltpc2D_After->Fill(ffTreeMC_pT,Float_t(ffTreeMC_ncltpc));
      //
      for (Int_t icent=0; icent<nCentBins; icent++){
        Bool_t centAcc = (ffTreeMC_cent>=centBinning[icent] && ffTreeMC_cent<centBinning[icent+1]);
        if ( centAcc ){

          tidenTreeStream[icent]->GetFile()->cd();
          (*tidenTreeStream[icent])<<"tracks"<<
          "gid="                  << ffTreeMC_gid          <<  //  global event ID
          "dEdx="                 << ffTreeMC_dEdx         <<  //  dEdx of the track
          "sign="                 << ffTreeMC_sign         <<  //  charge
          "ptot="                 << ffTreeMC_ptot         <<  //  TPC momentum
          "eta="                  << ffTreeMC_eta          <<  //  eta
          "cent="                 << ffTreeMC_cent         <<  //  centrality
          //
          "cutBit="               << ffTreeMC_cutBit       <<  //  Systematic Cuts
          "part="                 << ffTreeMC_part         <<  //  Systematic Cuts
          "p="                    << ffTreeMC_p            <<  //  TPC momentum
          "pT="                   << ffTreeMC_pT           <<
          "cRows="                << ffTreeMC_cRows        <<  // interaction rate
          "chi2tpc="              << ffTreeMC_chi2tpc      <<  // interaction rate
          "dcaz="                 << ffTreeMC_dcaz         <<  // interaction rate
          "dcaxy="                << ffTreeMC_dcaxy        <<  // interaction rate
          "\n";
        }
      }
    }
    //
    // apply dca and nclusters cut
    if (TMath::Abs(ffTreeMC_dcaz)>3) continue;
    if (TMath::Abs(ffTreeMC_dcaxy)>3) continue;
    //
    // Prepare dEdx histograms
    for (Int_t icent=0; icent<nCentBins; icent++){
      for (Int_t ieta=0; ieta<nEtaBins; ieta++){
        //
        if (testIntegratedHist && ieta>0) continue;
        //
        // trackOrigin = 0 --> IsPhysicalPrimary,  trackOrigin = 1 --> IsSecondaryFromMaterial, trackOrigin = 2 --> IsSecondaryFromWeakDecay
        //
        Bool_t etaCentString = (ffTreeMC_eta>=etaBinning[ieta] && ffTreeMC_eta<etaBinning[ieta+1] && ffTreeMC_cent>=centBinning[icent] && ffTreeMC_cent<centBinning[icent+1]);
        if (testIntegratedHist && ieta==0) {
          etaCentString = (ffTreeMC_cent>=centBinning[icent] && ffTreeMC_cent<centBinning[icent+1]);
        }
        //
        Bool_t parPos = ffTreeMC_sign>0.;
        Bool_t parNeg = ffTreeMC_sign<0.;
        Bool_t originCut = (fTrackOrigin) ? (ffTreeMC_origin==0) : (ffTreeMC_origin>-1);
        //
        Bool_t prTOF = ((ffTreeMC_cutBit >> kCleanPrTOF) & 1);
        Bool_t kaTOF = ((ffTreeMC_cutBit >> kCleanKaTOF) & 1);
        //
        Bool_t vertexPcut=kFALSE;
        if(fCutON=="pT")   vertexPcut = (ffTreeMC_pT>=ptotMin      && ffTreeMC_pT<=ptotMax);
        if(fCutON=="ptot") vertexPcut = (ffTreeMC_ptot>=ptotMin    && ffTreeMC_ptot<=ptotMax);
        if(fCutON=="p")    vertexPcut = (ffTreeMC_p>=ptotMin && ffTreeMC_p<=ptotMax);

        if (originCut && etaCentString && vertexPcut && systCut          ) h2Dall[icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
        if (originCut && etaCentString && vertexPcut && systCut && parPos) h2Dpos[icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
        if (originCut && etaCentString && vertexPcut && systCut && parNeg) h2Dneg[icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
        //
        if (originCut && etaCentString && vertexPcut && systCut && kaTOF)           h2DallKaTOF[icent][ieta]   ->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
        if (originCut && etaCentString && vertexPcut && systCut && kaTOF && parPos) h2DallKaTOFPos[icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
        if (originCut && etaCentString && vertexPcut && systCut && kaTOF && parNeg) h2DallKaTOFNeg[icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
        //
        if (originCut && etaCentString && vertexPcut && systCut && prTOF)           h2DallPrTOF[icent][ieta]   ->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
        if (originCut && etaCentString && vertexPcut && systCut && prTOF && parPos) h2DallPrTOFPos[icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
        if (originCut && etaCentString && vertexPcut && systCut && prTOF && parNeg) h2DallPrTOFNeg[icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
        //
        for (Int_t ipart=0; ipart<nParticles; ipart++){
          //
          // part=0 --> electron, part=1 --> pion, part=2 --> kaon, part=3 --> proton, part=4 --> deuteron
          if (originCut && etaCentString && vertexPcut && systCut && (ffTreeMC_part==ipart)) h2DClean[ipart][icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
          if (originCut && etaCentString && vertexPcut && systCut && (ffTreeMC_part==ipart) && parPos) h2DCleanPos[ipart][icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);
          if (originCut && etaCentString && vertexPcut && systCut && (ffTreeMC_part==ipart) && parNeg) h2DCleanNeg[ipart][icent][ieta]->Fill(ffTreeMC_ptot,ffTreeMC_dEdx);

        }

      }
    }

  }  // tree loop
  timer.Stop(); timer.Print();
  std::cout << " ========= ProcessDataHists DONE ========= #events = " << eventCount << std::endl;

}
//____________________________________________________________________________________________________________
void InitInitials()
{

  TString outputFileNameHist     = Form("Hists_Origin_%d_Syst%d.root" ,fTrackOrigin, fSystSet);
  TString debugFile              = Form("Debug_Origin_%d_Syst%d.root" ,fTrackOrigin, fSystSet);
  histStream      = new TTreeSRedirector(outputFileNameHist,"recreate");
  debugStream     = new TTreeSRedirector(debugFile,"recreate");
  //
  // Initialise histograms to be used for binning
  //
  if (fillTIdenTree){
    for (Int_t icent=0; icent<nCentBins; icent++){
      TString idenStrName = Form("TIdenTree_cent_%3.2f_%3.2f.root",centBinning[icent],centBinning[icent+1]);
      tidenTreeStream[icent] = new TTreeSRedirector(idenStrName,"recreate");
    }
  }
  //
  hfTreeMC_dcaxy2D              = new TH2F("hfTreeMC_dcaxy2D" ,"hfTreeMC_dcaxy2D"  ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  hfTreeMC_dcaz2D               = new TH2F("hfTreeMC_dcaz2D"  ,"hfTreeMC_dcaz2D"   ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  hfTreeMC_ncltpc2D             = new TH2F("hfTreeMC_ncltpc2D","hfTreeMC_ncltpc2D" ,ptNbins,ptotMin,ptotMax, 140 , 30., 170.);
  //
  hfTreeMC_dcaxy2D_After        = new TH2F("hfTreeMC_dcaxy2D_After" ,"hfTreeMC_dcaxy2D_After"  ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  hfTreeMC_dcaz2D_After         = new TH2F("hfTreeMC_dcaz2D_After"  ,"hfTreeMC_dcaz2D_After"   ,ptNbins,ptotMin,ptotMax, 400 ,-10., 10. );
  hfTreeMC_ncltpc2D_After       = new TH2F("hfTreeMC_ncltpc2D_After","hfTreeMC_ncltpc2D_After" ,ptNbins,ptotMin,ptotMax, 140 , 30., 170.);


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
      h2DallKaTOF[icent][ieta]=NULL;      h2DallKaTOFPos[icent][ieta]=NULL;      h2DallKaTOFNeg[icent][ieta]=NULL;
      hName2DallKaTOF[icent][ieta]   =Form("h2DallKaTOF_%s"   ,centEtaStr.Data());
      hName2DallKaTOFPos[icent][ieta]=Form("h2DallKaTOFPos_%s",centEtaStr.Data());
      hName2DallKaTOFNeg[icent][ieta]=Form("h2DallKaTOFNeg_%s",centEtaStr.Data());
      h2DallKaTOF[icent][ieta]     = new TH2F(hName2DallKaTOF[icent][ieta]    ,hName2DallKaTOF[icent][ieta]     ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallKaTOFPos[icent][ieta]  = new TH2F(hName2DallKaTOFPos[icent][ieta] ,hName2DallKaTOFPos[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      h2DallKaTOFNeg[icent][ieta]  = new TH2F(hName2DallKaTOFNeg[icent][ieta] ,hName2DallKaTOFNeg[icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
      //
      for (Int_t ipart=0; ipart<nParticles; ipart++)
      {
        h2DCleanPos[ipart][icent][ieta]=NULL;
        hName2DCleanPos[ipart][icent][ieta] = Form("h2DClean%sPos_%s", parName[ipart].Data(), centEtaStr.Data());
        h2DCleanPos[ipart][icent][ieta]     = new TH2F(hName2DCleanPos[ipart][icent][ieta],hName2DCleanPos[ipart][icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
        //
        h2DCleanNeg[ipart][icent][ieta]=NULL;
        hName2DCleanNeg[ipart][icent][ieta] = Form("h2DClean%sNeg_%s", parName[ipart].Data(), centEtaStr.Data());
        h2DCleanNeg[ipart][icent][ieta]     = new TH2F(hName2DCleanNeg[ipart][icent][ieta],hName2DCleanNeg[ipart][icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
        //
        h2DClean[ipart][icent][ieta]=NULL;
        hName2DClean[ipart][icent][ieta]    = Form("h2DClean%s_%s", parName[ipart].Data(), centEtaStr.Data());
        h2DClean[ipart][icent][ieta]        = new TH2F(hName2DClean[ipart][icent][ieta],hName2DClean[ipart][icent][ieta]  ,ptNbins,ptotMin,ptotMax,dEdxNbins,dEdxMin,dEdxMax);
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
  1 -->  CRows60
  2 -->  CRows100
  3 -->  Chi2TPCDmall
  4 -->  Chi2TPCLarge
  5 -->  DCAXYSmall
  6 -->  DCAXYLarge
  7 -->  VZSmall
  8 -->  VZLarge
  9 -->  EventVertexZSmall
  10 --> EventVertexZLarge
  11 --> RequireITSRefit
  12 --> NewITSCut
  13 --> PixelRequirementITS
  // extra settings
  14 --> kTPCSignalN,
  15 --> kActiveZone,
  16 --> kTPCSignalN+kActiveZone,
  17 --> kTPCSignalNSmall+kActiveZoneSmall,
  18 --> kTPCSignalNLarge,kActiveZoneLarge,
  //
  19 --> kTPCSignalNLarge,kActiveZoneLarge,
  20 --> kTPCSignalNLarge,kActiveZoneLarge,
  */


  Bool_t cutAll = kFALSE;
  const Int_t fnCutBins=10;
  Int_t fCutArrTmp[fnCutBins]={0};

  switch(syst) {

    case kCutReference:   // 0 -->  Reference
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutCrossedRowsTPC60:  // 1 -->  CRows60
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC60,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutCrossedRowsTPC100:  // 2 -->  CRows100
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC100, kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCSmall:   // 3 -->  Chi2TPCSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCSmall, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutMaxChi2PerClusterTPCLarge:   // 4 -->  Chi2TPCLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPCLarge, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutMaxDCAToVertexXYPtDepSmall:   // 5 -->  kMaxDCAToVertexXYPtDepSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDepSmall, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutMaxDCAToVertexXYPtDepLarge:   // 6 -->  kMaxDCAToVertexXYPtDepLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDepLarge, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutVertexZSmall:   // 7 -->  kVertexZSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZSmall, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutVertexZLarge:   // 8 -->  kVertexZLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZLarge, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutEventVertexZSmall:  // 9 -->  kEventVertexZSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZSmall, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutEventVertexZLarge:  // 10 -->  kEventVertexZLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZLarge, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutRequireITSRefit:  // 11 -->  no kRequireITSRefit
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, 1, 1, 1, 1 ,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutPixelRequirementITS:  // 12 --> no  kNewITSCut
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, 1, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutNewITSCut:  // 13 -->  kPixelRequirementITS
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, 1, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    // ----------------------------------------------------------------------------------------------------
    //                                              Dangerous cuts
    // ----------------------------------------------------------------------------------------------------
    //
    case kCutTPCSignalN:  // 14 -->  kTPCSignalN
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalN,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutActiveZone:  // 15 -->  kActiveZone
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,kActiveZone};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutTPCSignalNActiveZone:  // 16 -->  kTPCSignalN + kActiveZone
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalN,kActiveZone};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutTPCSignalNSmallActiveZoneSmall:  // 17 -->  tightest cut kTPCSignalNSmall + kActiveZoneSmall
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalNSmall,kActiveZoneSmall};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    case kCutTPCSignalNLargeActiveZoneLarge:  // 18 -->  kTPCSignalNLarge + kActiveZoneLarge
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZ, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,kTPCSignalNLarge,kActiveZoneLarge};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = ((cutBit >> fCutArrTmp[i]) & 1);
    }
    break;
    //
    // ----------------------------------------------------------------------------------------------------
    //                                              Extra vZ cuts
    // ----------------------------------------------------------------------------------------------------
    //
    //
    case kCutEventVertexZALICE:  // 19 -->  kEventVertexZALICE
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZALICE, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case kCutEventVertexZALICETight:  // 20 -->  kEventVertexZALICETight
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, kMaxDCAToVertexXYPtDep, kVertexZ, kEventVertexZALICETight, kPixelRequirementITS, kNewITSCut, kRequireITSRefit,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
    }
    break;
    //
    case 21:  // 21 -->  dcaxy<2.4 dcaz<3.2
    {
      Int_t fCutArrTmp[fnCutBins] = {kNCrossedRowsTPC80,  kMaxChi2PerClusterTPC, 1, 1, kEventVertexZ, kPixelRequirementITS, kNewITSCut, 1,1,1};
      for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
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
      if(h2DallKaTOF[icent][ieta])     h2DallKaTOF[icent][ieta]->Write();
      if(h2DallKaTOFPos[icent][ieta])  h2DallKaTOFPos[icent][ieta]->Write();
      if(h2DallKaTOFNeg[icent][ieta])  h2DallKaTOFNeg[icent][ieta]->Write();
      //
      for (Int_t ipart=0; ipart<nParticles; ipart++){
        if(h2DCleanPos[ipart][icent][ieta]) h2DCleanPos[ipart][icent][ieta]->Write();
        if(h2DCleanNeg[ipart][icent][ieta]) h2DCleanNeg[ipart][icent][ieta]->Write();
        if(h2DClean[ipart][icent][ieta])    h2DClean[ipart][icent][ieta]->Write();
      }
    }
  }
  //
  debugStream->GetFile()->cd();
  if(hfTreeMC_dcaxy2D)   hfTreeMC_dcaxy2D  ->Write();
  if(hfTreeMC_dcaz2D)    hfTreeMC_dcaz2D   ->Write();
  if(hfTreeMC_ncltpc2D)  hfTreeMC_ncltpc2D ->Write();
  if(hfTreeMC_dcaxy2D_After)   hfTreeMC_dcaxy2D_After  ->Write();
  if(hfTreeMC_dcaz2D_After)    hfTreeMC_dcaz2D_After   ->Write();
  if(hfTreeMC_ncltpc2D_After)  hfTreeMC_ncltpc2D_After ->Write();
  //
  if(htracks_dcaxy1D)    htracks_dcaxy1D  ->Write();
  if(htracks_dcaz1D)     htracks_dcaz1D   ->Write();
  if(htracks_ncltpc1D)   htracks_ncltpc1D ->Write();
  if(htracks_eta1D)      htracks_eta1D ->Write();
  if(htracks_cRows1D)    htracks_cRows1D ->Write();
  if(htracks_cent1D)     htracks_cent1D ->Write();
  if(htracks_phi1D)      htracks_phi1D ->Write();
  if(htracks_chi2tpc1D)  htracks_chi2tpc1D ->Write();
  if(htracks_vz1D)       htracks_vz1D ->Write();

}
//____________________________________________________________________________________________________________
void SetBranchAddresses()
{
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
    eventtree->SetBranchAddress("itsVertexInfo."      ,&fevents_itsVertexInfo);
  }
  //
  // mc tree
  if (mctree){

    mctree->SetBranchAddress("cutBit"          ,&ffTreeMC_cutBit);
    mctree->SetBranchAddress("gid"             ,&ffTreeMC_gid);
    mctree->SetBranchAddress("dEdx"       ,&ffTreeMC_dEdx);
    mctree->SetBranchAddress("sign"       ,&ffTreeMC_sign);
    mctree->SetBranchAddress("isample"       ,&ffTreeMC_isample);
    mctree->SetBranchAddress("origin"       ,&ffTreeMC_origin);
    mctree->SetBranchAddress("part"       ,&ffTreeMC_part);
    mctree->SetBranchAddress("ptot"       ,&ffTreeMC_ptot);
    mctree->SetBranchAddress("p"       ,&ffTreeMC_p);
    mctree->SetBranchAddress("pT"       ,&ffTreeMC_pT);
    mctree->SetBranchAddress("eta"       ,&ffTreeMC_eta);
    mctree->SetBranchAddress("cent"       ,&ffTreeMC_cent);
    mctree->SetBranchAddress("centimp"       ,&ffTreeMC_centimp);
    mctree->SetBranchAddress("vZ"       ,&ffTreeMC_vZ);
    mctree->SetBranchAddress("ncltpc"       ,&ffTreeMC_ncltpc);
    mctree->SetBranchAddress("dcaxy"       ,&ffTreeMC_dcaxy);
    mctree->SetBranchAddress("dcaz"       ,&ffTreeMC_dcaz);
    mctree->SetBranchAddress("cRows"       ,&ffTreeMC_cRows);
    mctree->SetBranchAddress("chi2tpc"       ,&ffTreeMC_chi2tpc);
    mctree->SetBranchAddress("phi"       ,&ffTreeMC_phi);
    mctree->SetBranchAddress("nsigmatofka"       ,&ffTreeMC_nsigmatofka);
    mctree->SetBranchAddress("nsigmatofpr"       ,&ffTreeMC_nsigmatofpr);

  }

}
//____________________________________________________________________________________________________________
void PrepareEffMatrix(TString matrixFile, Int_t detectorCom, Int_t systSetting, Float_t pDownAcc, Float_t pUpAcc)
{

  /*

  /media/marsland/Samsung_T5/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN2/LHC16g1/LHC16g1_pass1_EffCheckForPaper_12112020/mergedRuns/AnalysisResults_hists.root

  cd /home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/MChists
  aliroot -l
  // .L /home/marsland/Desktop/ubuntu_desktop/workdir/code/MCData_FilterTreesMakeHists.C+
  // .L /home/marsland/Desktop/ubuntu_desktop/workdir/fourthMomentAnalysis/code/MCData_FilterTreesMakeHists.C+
  // TString histFile = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN2/LHC16g1/LHC16g1_pass1_NoSelection_06082019/mergedRuns/mergedHists/AnalysisResults_hists.root";


  aliroot -l
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/gsi_marsland_EbyeRatios/MCData_FilterTreesMakeHists.C+
  TString histFile = "/home/marsland/Desktop/QM19_tmp/PbPb/MC/RUN2/LHC16g1/LHC16g1_pass1_NoSelection_06082019/mergedRuns/mergedHists/AnalysisResults_hists.root";

  PrepareEffMatrix(histFile,0,    0   ,0.6, 2.);  // TPC eff.
  PrepareEffMatrix(histFile,0,    0   ,0.6, 1.5); // TPC eff.

  PrepareEffMatrix(histFile,1,    10   ,0.6, 2.);  // TPC+TOF eff.
  PrepareEffMatrix(histFile,1,    10   ,0.6, 1.5); // TPC+TOF eff.

  aliroot -l
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/code/gsi_marsland_EbyeRatios/MCData_FilterTreesMakeHists.C+
  // TString matrixFile = "/media/marsland/Samsung_T5/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN2/LHC16g1/LHC16g1_pass1_EffCheckForPaper_12112020/mergedRuns/AnalysisResults_hists.root";
  TString matrixFile = "/media/marsland/Samsung_T5/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN2/LHC16g1/LHC16g1_pass1_NoSelection_06082019/mergedRuns/mergedHists/AnalysisResults_hists.root";
  PrepareEffMatrix(matrixFile,0,    0    ,0.6, 1.5); // TPC eff.
  PrepareEffMatrix(matrixFile,0,    11   ,0.6, 1.5); // TPC eff.
  PrepareEffMatrix(matrixFile,0,    15   ,0.6, 1.5); // TPC eff.
  PrepareEffMatrix(matrixFile,0,    21   ,0.6, 1.5); // TPC eff.
  PrepareEffMatrix(matrixFile,0,    0    ,0.6, 2.0); // TPC eff.
  PrepareEffMatrix(matrixFile,0,    11   ,0.6, 2.0); // TPC eff.
  PrepareEffMatrix(matrixFile,0,    15   ,0.6, 2.0); // TPC eff.
  PrepareEffMatrix(matrixFile,0,    21   ,0.6, 2.0); // TPC eff.

  */
  fSystSet = systSetting;
  if (systSetting>20) fSystSet=0;

  gStyle->SetPadRightMargin(0.15);
  gROOT->ProcessLine(".x ~/bin/color.C");
  TGaxis::SetMaxDigits(5);
  Double_t etaDown = -0.8;
  Double_t etaUp   =  0.8;

  TFile *fEffMatrix = new TFile(matrixFile);
  TList *list = (TList*)fEffMatrix -> Get("cleanHists");
  posRec = (THnF*)list -> FindObject("fHistPosEffMatrixScanRec");
  posGen = (THnF*)list -> FindObject("fHistPosEffMatrixScanGen");
  negRec = (THnF*)list -> FindObject("fHistNegEffMatrixScanRec");
  negGen = (THnF*)list -> FindObject("fHistNegEffMatrixScanGen");
  //
  ModifyEffMatrix(posRec);
  ModifyEffMatrix(posGen);
  ModifyEffMatrix(negRec);
  ModifyEffMatrix(negGen);
  //
  // axis = 0 -->   detector
  // axis = 1 -->   setting
  // axis = 2 -->   particle type
  // axis = 3 -->   Centrality (%)
  // axis = 4 -->   #it{p}_{T} (GeV/#it{c})
  // axis = 5 -->   #eta
  //
  // 0:detector, 1:syst, 2: particle, 3:cent, 4:momentum, 5:eta
  posRec -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  posGen -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  negRec -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  negGen -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  //
  posRec -> GetAxis(kMom)-> SetRangeUser(pDownAcc,pUpAcc);
  posGen -> GetAxis(kMom)-> SetRangeUser(pDownAcc,pUpAcc);
  negRec -> GetAxis(kMom)-> SetRangeUser(pDownAcc,pUpAcc);
  negGen -> GetAxis(kMom)-> SetRangeUser(pDownAcc,pUpAcc);
  //
  posRec -> GetAxis(kDet)-> SetRangeUser(detectorCom, detectorCom+1);
  negRec -> GetAxis(kDet)-> SetRangeUser(detectorCom, detectorCom+1);
  posGen -> GetAxis(kDet)-> SetRangeUser(detectorCom, detectorCom+1);
  negGen -> GetAxis(kDet)-> SetRangeUser(detectorCom, detectorCom+1);
  //
  posRec -> GetAxis(kSyst)-> SetRangeUser(fSystSet, fSystSet+1);
  negRec -> GetAxis(kSyst)-> SetRangeUser(fSystSet, fSystSet+1);
  posGen -> GetAxis(kSyst)-> SetRangeUser(fSystSet, fSystSet+1);
  negGen -> GetAxis(kSyst)-> SetRangeUser(fSystSet, fSystSet+1);
  //
  effMatrixStream = new TTreeSRedirector(Form("EffMatrix_Syst_%d_Det_%d_Acc_%3.2f_%3.2f.root",systSetting,detectorCom,pDownAcc,pUpAcc),"recreate");

  MeanEffInAcceptance(0);
  MeanEffInAcceptance(1);
  MeanEffInAcceptance(2);

  for (Int_t icent=0; icent<nCentAcc; icent++){
    effMatrix1D(4,0,icent);
    effMatrix1D(4,1,icent);
    effMatrix1D(4,2,icent);
    //
    effMatrix1D(5,0,icent);
    effMatrix1D(5,1,icent);
    effMatrix1D(5,2,icent);
    //
    effMatrix2D(0,icent);
    effMatrix2D(1,icent);
    effMatrix2D(2,icent);
  }

}

void ModifyEffMatrix(THnF *hnF)
{

  Int_t nDetBins = hnF->GetAxis(0)->GetNbins();
  Int_t nSetBins = hnF->GetAxis(1)->GetNbins();
  Int_t nParBins = hnF->GetAxis(2)->GetNbins();
  Int_t nCenBins = hnF->GetAxis(3)->GetNbins();
  Int_t nMomBins = hnF->GetAxis(4)->GetNbins();
  Int_t nEtaBins = hnF->GetAxis(5)->GetNbins();
  cout << " axis = 0 -->   " << hnF->GetAxis(0)->GetTitle() << "  Nbins = " << hnF->GetAxis(0)->GetNbins() << endl;
  cout << " axis = 1 -->   " << hnF->GetAxis(1)->GetTitle() << "  Nbins = " << hnF->GetAxis(1)->GetNbins() << endl;
  cout << " axis = 2 -->   " << hnF->GetAxis(2)->GetTitle() << "  Nbins = " << hnF->GetAxis(2)->GetNbins() << endl;
  cout << " axis = 3 -->   " << hnF->GetAxis(3)->GetTitle() << "  Nbins = " << hnF->GetAxis(3)->GetNbins() << endl;
  cout << " axis = 4 -->   " << hnF->GetAxis(4)->GetTitle() << "  Nbins = " << hnF->GetAxis(4)->GetNbins() << endl;
  cout << " axis = 5 -->   " << hnF->GetAxis(5)->GetTitle() << "  Nbins = " << hnF->GetAxis(5)->GetNbins() << endl;
  //
  for (Int_t idet = 1; idet <=nDetBins; idet++) {
    for (Int_t iset = 1; iset <=nSetBins; iset++) {
      for (Int_t ipar = 1; ipar <=nParBins; ipar++) {
        for (Int_t icen = 1; icen <=nCenBins; icen++) {
          for (Int_t imom = 1; imom <=nMomBins; imom++) {

            Int_t a8[6] ={idet,iset,ipar,icen,imom,8};
            Int_t a9[6] ={idet,iset,ipar,icen,imom,9};
            Int_t a7[6] ={idet,iset,ipar,icen,imom,7};
            Int_t a10[6]={idet,iset,ipar,icen,imom,10};
            hnF->SetBinContent(a8,hnF->GetBinContent(a7));
            hnF->SetBinContent(a9,hnF->GetBinContent(a10));
            // hnF->SetBinContent(a8,0.01);
            // hnF->SetBinContent(a9,0.01);

          }
        }
      }
    }
  }

}

void effMatrix1D(Int_t effType,Int_t pid, Int_t cent)
{

  //
  // pion    pid  ==  0;
  // kaon    pid  ==  1;
  // proton  pid  ==  2;
  //
  // mom dependence effType=4
  // eta dependence effType=5
  //
  Double_t etaDown = -0.8;
  Double_t etaUp   =  0.8;
  posRec -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  posGen -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  negRec -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  negGen -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);

  posRec -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);
  posGen -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);
  negRec -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);
  negGen -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);

  posGen -> GetAxis(kCent)-> SetRangeUser(centAccDown[cent],centAccUp[cent]);
  posRec -> GetAxis(kCent)-> SetRangeUser(centAccDown[cent],centAccUp[cent]);
  negGen -> GetAxis(kCent)-> SetRangeUser(centAccDown[cent],centAccUp[cent]);
  negRec -> GetAxis(kCent)-> SetRangeUser(centAccDown[cent],centAccUp[cent]);

  TH1F *histoRecPos = (TH1F*)posRec -> Projection(effType);  histoRecPos -> Sumw2();
  TH1F *histoGenPos = (TH1F*)posGen -> Projection(effType);  histoGenPos -> Sumw2();
  TH1F *histoRecNeg = (TH1F*)negRec -> Projection(effType);  histoRecNeg -> Sumw2();
  TH1F *histoGenNeg = (TH1F*)negGen -> Projection(effType);  histoGenNeg -> Sumw2();

  histoRecPos -> Divide(histoGenPos);
  histoRecNeg -> Divide(histoGenNeg);

  histoRecPos -> SetMarkerStyle(20);
  histoRecPos -> SetMarkerColor(colors[cent]);
  histoRecPos -> SetLineColor(colors[cent]);
  histoRecPos -> SetLineWidth(3);
  histoRecPos -> SetMarkerSize(1.8);

  histoRecNeg -> SetMarkerStyle(25);
  histoRecNeg -> SetMarkerColor(colors[cent]);
  histoRecNeg -> SetLineColor(colors[cent]);
  histoRecNeg -> SetLineWidth(3);
  histoRecNeg -> SetMarkerSize(1.8);

  TString xAxisNAme;
  if (effType==4) xAxisNAme = "#it{p} (GeV/#it{c})";
  if (effType==5) xAxisNAme = "#eta";

  histoRecPos -> GetYaxis() -> SetRangeUser(0,1);
  histoRecPos -> GetXaxis() -> SetTitle(xAxisNAme);
  if(pid == 0) histoRecPos -> GetYaxis() -> SetTitle("pion detection efficiency");
  if(pid == 1) histoRecPos -> GetYaxis() -> SetTitle("kaon detection efficiency");
  if(pid == 2) histoRecPos -> GetYaxis() -> SetTitle("proton detection efficiency");

  histoRecNeg -> GetYaxis() -> SetRangeUser(0,1);
  histoRecNeg -> GetXaxis() -> SetTitle(xAxisNAme);
  if(pid == 0) histoRecNeg -> GetYaxis() -> SetTitle("pion detection efficiency");
  if(pid == 1) histoRecNeg -> GetYaxis() -> SetTitle("kaon detection efficiency");
  if(pid == 2) histoRecNeg -> GetYaxis() -> SetTitle("proton detection efficiency");

  if (effType==4){
    h1Mom[0][pid][cent]=(TH1F*)histoRecPos->Clone();
    h1Mom[1][pid][cent]=(TH1F*)histoRecNeg->Clone();

    h1Mom[0][pid][cent]->SetName(Form("h1DPos_%s_%d_Cent_%d",parName[pid+1].Data(),effType,cent));
    h1Mom[1][pid][cent]->SetName(Form("h1DNeg_%s_%d_Cent_%d",parName[pid+1].Data(),effType,cent));

    effMatrixStream->GetFile()->cd();
    h1Mom[0][pid][cent]->Write();
    h1Mom[1][pid][cent]->Write();
  }

  if (effType==5){
    h1Eta[0][pid][cent]=(TH1F*)histoRecPos->Clone();
    h1Eta[1][pid][cent]=(TH1F*)histoRecNeg->Clone();

    h1Eta[0][pid][cent]->SetName(Form("h1DPos_%s_%d_Cent_%d",parName[pid+1].Data(),effType,cent));
    h1Eta[1][pid][cent]->SetName(Form("h1DNeg_%s_%d_Cent_%d",parName[pid+1].Data(),effType,cent));
    //
    effMatrixStream->GetFile()->cd();
    h1Eta[0][pid][cent]->Write();
    h1Eta[1][pid][cent]->Write();
  }

}

void effMatrix2D(Int_t pid, Int_t cent)
{

  //
  //pion    pid  ==  0;
  //kaon    pid  ==  1;
  //proton  pid  ==  2;
  //
  Double_t etaDown = -0.8;
  Double_t etaUp   =  0.8;
  posRec -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  posGen -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  negRec -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);
  negGen -> GetAxis(kEta)-> SetRangeUser(etaDown,etaUp);

  posRec -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);
  posGen -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);
  negRec -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);
  negGen -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);

  posGen -> GetAxis(kCent)-> SetRangeUser(centAccDown[cent],centAccUp[cent]);
  posRec -> GetAxis(kCent)-> SetRangeUser(centAccDown[cent],centAccUp[cent]);
  negGen -> GetAxis(kCent)-> SetRangeUser(centAccDown[cent],centAccUp[cent]);
  negRec -> GetAxis(kCent)-> SetRangeUser(centAccDown[cent],centAccUp[cent]);

  TH2F *histoRecPos = (TH2F*)posRec -> Projection(4,5);  histoRecPos -> Sumw2();
  TH2F *histoGenPos = (TH2F*)posGen -> Projection(4,5);  histoGenPos -> Sumw2();
  TH2F *histoRecNeg = (TH2F*)negRec -> Projection(4,5);  histoRecNeg -> Sumw2();
  TH2F *histoGenNeg = (TH2F*)negGen -> Projection(4,5);  histoGenNeg -> Sumw2();

  histoRecPos -> Divide(histoGenPos);
  histoRecNeg -> Divide(histoGenNeg);
  histoRecPos -> GetXaxis() -> SetTitle("#eta");
  histoRecPos -> GetYaxis() -> SetTitle("#it{p} (GeV/#it{c})");
  histoRecNeg -> GetXaxis() -> SetTitle("#eta");
  histoRecNeg -> GetYaxis() -> SetTitle("#it{p} (GeV/#it{c})");

  h2EtaMom[0][pid][cent]=(TH2F*)histoRecPos->Clone();
  h2EtaMom[1][pid][cent]=(TH2F*)histoRecNeg->Clone();

  h2EtaMom[0][pid][cent]->SetName(Form("h2D_EtaMom_Pos_%s_Cent_%d",parName[pid+1].Data(),cent));
  h2EtaMom[1][pid][cent]->SetName(Form("h2D_EtaMom_Neg_%s_Cent_%d",parName[pid+1].Data(),cent));

  effMatrixStream->GetFile()->cd();
  h2EtaMom[0][pid][cent]->Write();
  h2EtaMom[1][pid][cent]->Write();

}

void MeanEffInAcceptance(Int_t pid)
{

  //
  //pion    pid  ==  0;
  //kaon    pid  ==  1;
  //proton  pid  ==  2;
  //

  posGen -> GetAxis(kCent)-> SetRangeUser(centAccDown[0],centAccUp[8]);
  posRec -> GetAxis(kCent)-> SetRangeUser(centAccDown[0],centAccUp[8]);
  negGen -> GetAxis(kCent)-> SetRangeUser(centAccDown[0],centAccUp[8]);
  negRec -> GetAxis(kCent)-> SetRangeUser(centAccDown[0],centAccUp[8]);

  posRec -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);
  posGen -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);
  negRec -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);
  negGen -> GetAxis(kPart)-> SetRangeUser(pid,pid+1);

  TString centEtaStr_Pos = Form("h2D_CentEta_Pos_%s",parName[pid+1].Data());
  TString centEtaStr_Neg = Form("h2D_CentEta_Neg_%s",parName[pid+1].Data());

  h2CentEta[0][pid] = new TH2F(centEtaStr_Pos, centEtaStr_Pos , nEtaAcc,0.,nEtaAcc, nCentAcc ,0., nCentAcc );
  h2CentEta[1][pid] = new TH2F(centEtaStr_Neg ,centEtaStr_Neg , nEtaAcc,0.,nEtaAcc, nCentAcc ,0., nCentAcc );

  for (Int_t ieta=0; ieta<nEtaAcc; ieta++){

    // 0:detector, 1:syst, 2: particle, 3:cent, 4:momentum, 5:eta
    posRec -> GetAxis(kEta)-> SetRangeUser(etaAccDown[ieta],etaAccUp[ieta]);
    posGen -> GetAxis(kEta)-> SetRangeUser(etaAccDown[ieta],etaAccUp[ieta]);
    negRec -> GetAxis(kEta)-> SetRangeUser(etaAccDown[ieta],etaAccUp[ieta]);
    negGen -> GetAxis(kEta)-> SetRangeUser(etaAccDown[ieta],etaAccUp[ieta]);
    //
    TH1F *histoRecPos = (TH1F*)posRec -> Projection(3);  histoRecPos -> Sumw2();
    TH1F *histoGenPos = (TH1F*)posGen -> Projection(3);  histoGenPos -> Sumw2();
    TH1F *histoRecNeg = (TH1F*)negRec -> Projection(3);  histoRecNeg -> Sumw2();
    TH1F *histoGenNeg = (TH1F*)negGen -> Projection(3);  histoGenNeg -> Sumw2();

    histoRecPos -> Divide(histoGenPos);
    histoRecNeg -> Divide(histoGenNeg);

    h1CentEta[0][pid][ieta]=(TH1F*)histoRecPos->Clone();
    h1CentEta[1][pid][ieta]=(TH1F*)histoRecNeg->Clone();

    h1CentEta[0][pid][ieta]->SetName(Form("h1D_CentEta_Pos_%s_EtaAcc_%d",parName[pid+1].Data(),ieta));
    h1CentEta[1][pid][ieta]->SetName(Form("h1D_CentEta_Neg_%s_EtaAcc_%d",parName[pid+1].Data(),ieta));

    for (Int_t icent=0; icent<nCentAcc; icent++){

      cout << etaAccDown[ieta]   << " - " << etaAccUp[ieta] << "   ";
      cout << centAccDown[icent] << " - " << centAccUp[icent];
      cout << "   Eff = " << h1CentEta[0][pid][ieta]->GetBinContent(icent+1) << "  ---  " << h1CentEta[1][pid][ieta]->GetBinContent(icent+1) << endl;
      h2CentEta[0][pid] -> Fill(ieta,icent,h1CentEta[0][pid][ieta]->GetBinContent(icent+1));
      h2CentEta[1][pid] -> Fill(ieta,icent,h1CentEta[1][pid][ieta]->GetBinContent(icent+1));
      if (icent==nCentAcc-1) {
        h2CentEta[0][pid] -> Fill(ieta,icent,h1CentEta[0][pid][ieta]->GetBinContent(nCentAcc-1));
        h2CentEta[1][pid] -> Fill(ieta,icent,h1CentEta[1][pid][ieta]->GetBinContent(nCentAcc-1));
      }

    }

    effMatrixStream->GetFile()->cd();
    h1CentEta[0][pid][ieta]->Write();
    h1CentEta[1][pid][ieta]->Write();

  }

  effMatrixStream->GetFile()->cd();
  h2CentEta[0][pid]->Write();
  h2CentEta[1][pid]->Write();


}
