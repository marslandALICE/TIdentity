#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TLine.h"
#include "TVectorF.h"
#include "TClonesArray.h"
#include "TTreeStream.h"
#include "TDatabasePDG.h"
#include "TCut.h"
#include "TProfile.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using std::cout;
using std::setw;

void InitHists();
void PlotdEdx(TString datafile, TString lineShapesFile);
void PrepareTemplates();
void ReadTemplates();
Bool_t AcceptTrack(Int_t effSetting, Float_t mom, Float_t eta, Int_t sign, Int_t part, Float_t cent, Float_t event);
void SetDataTreeAliases(TTree *tree);
TGraphErrors *ConvertProfToGraph(TProfile *hh);
void ApplyEffCorrectionAllMoments(TString inputDataFile, Int_t nEvents, Float_t bgFraction, Int_t effSetting,  Int_t centSetting);

/*

cd /home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/workdir
aliroot -l
.L /home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/TOYMC_HighMoments.C+
TOYMC_HighMoments(0.1, 0, 0)

*/
Int_t nSubSample = 25;
Int_t switchOffParticle=-100;   // default = -100 for 4 particles
Int_t fUsedSign = 0;
const Int_t nParticles=4;
const Int_t nNetParticles=3;
//
Float_t elMean  = 2, prMean =10;
Float_t elBMean = 2, prBMean=10;
Int_t nTracksPerEventArr[nParticles]={0};
Double_t elParams[]={10.,1.5};
Double_t prParams[]={10.,1.5};
const Int_t nMaxTracksPerEvent = 10000;
Float_t trackdEdx[nMaxTracksPerEvent]={0.};
Float_t trackMom[nMaxTracksPerEvent]={0.};
Float_t trackEta[nMaxTracksPerEvent]={0.};
Float_t trackPhi[nMaxTracksPerEvent]={0.};
Int_t trackID[nMaxTracksPerEvent]={0};
Int_t trackSign[nMaxTracksPerEvent]={0};

TH1D *hParticles[nParticles];
TH1D *hFirstMoms[nParticles];
TF1 *fParticles[nParticles];
TF1 *fNetParticles[nNetParticles];
UInt_t cutBit=0;
const Int_t colors[] = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray LShapesArr("TF1",5000);
TClonesArray LShapesArr_Old("TF1",5000);
enum momentType{kEl=0,kPr=1,kBEl=2,kBPr=3,};
//
Int_t centDown[9] = {0,5, 10,20,30,40,50,60,70};
Int_t centUp[9]   = {5,10,20,30,40,50,60,70,80};
//
TH1D *hEta, *hMom0, *hMom1, *hPhi;
TH1D *hCent[9];
TH1D *hEffTPC_Eta[2];
TH1D *hEffTPC_Mom[2];
TH1D *hEffTOF_Eta[2];
TH1D *hEffTOF_Mom[2];
TH2D *hEffTPC_EtaMom[2];
TH2D *hEffTOF_EtaMom[2];
TH1D *hMeanEff[2];
//
TRandom randomGen;
Float_t fEffAPr = 0.;
Float_t fEffPr  = 0.;
Float_t fPrEffCorFactor = 0.91;
Float_t fAPrEffCorFactor = 0.84;
THnF *fHistPosEffMatrixScanRec = NULL;
THnF *fHistNegEffMatrixScanRec = NULL;
THnF *fHistPosEffMatrixScanGen = NULL;
THnF *fHistNegEffMatrixScanGen = NULL;
Float_t fMomDown = 0.6;
Float_t fMomUp   = 1.5;
Float_t fEtaDown = -0.8;
Float_t fEtaUp   =  0.8;
//
Float_t fBGFraction = 0.1;
ULong64_t fNevents=50000;
Int_t fEffSetting = 0;
Int_t fCentSetting = 0;
Float_t fCent = 0;
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
};

//
//
//
// ================================================================================================
//
//
void TOYMC_HighMoments(Int_t nEvents, Float_t bgFraction, Int_t effSetting,  Int_t centSetting)
{
  /*

  // bgFraction  --> fraction of background wrt protons
  // effSetting  --> selection of eff.
  // centSetting --> centSetting; 0: cent taken from templates volume fluct, 1: fixed multiplicity for a given cent bin, no volume fluct

  cd /home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/workdir
  aliroot -l
  .L /home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/TOYMC_HighMoments.C+
  TOYMC_HighMoments(50000, 0.1, 0, 0)
  PlotdEdx("DataTree_Net_50000_BG_0.10_eff_0_VF_0.root", "LineShapes_Net_50000_BG_0.10_eff_0_VF_0.root")

  */

  fBGFraction  = bgFraction;
  fNevents     = nEvents;
  fEffSetting  = effSetting;
  fCentSetting = centSetting;
  ReadTemplates();
  InitHists();
  if (fEffSetting==0){
    fAPrEffCorFactor = hMeanEff[0]->GetMean();
    fPrEffCorFactor  = hMeanEff[1]->GetMean();
    cout << "meanEffs = " <<fAPrEffCorFactor << " -- " << fPrEffCorFactor << endl;
  }
  //
  TDatabasePDG pdg;
  const Double_t masses[3]={pdg.GetParticle("pi+")->Mass(), pdg.GetParticle("K+")->Mass(), pdg.GetParticle("proton")->Mass()};
  //
  TTreeSRedirector *outputData = new TTreeSRedirector(Form("DataTree_Net_%d_BG_%3.2f_eff_%d_VF_%d.root",nEvents,bgFraction,effSetting,centSetting), "recreate");
  TTreeSRedirector *outputFits = new TTreeSRedirector(Form("LineShapes_Net_%d_BG_%3.2f_eff_%d_VF_%d.root",nEvents,bgFraction,effSetting,centSetting), "recreate");
  LShapesArr.SetOwner(kTRUE);
  LShapesArr_Old.SetOwner(kTRUE);
  for (Int_t i=0;i<nParticles;i++) {
    hParticles[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),3000,-50.,50.);
    hFirstMoms[i] = new TH1D(Form("firstMom_%d",i),Form("hist_%d",i),3000,-50.,50.);
    fParticles[i] = new TF1(Form("particle_%d",i),"gaus",-50,50);
  }

  //
  // TRandom randomGen;
  Float_t subsampleID = -1;
  for (ULong64_t ievent=0;ievent<fNevents;ievent++){  // event loop
    subsampleID = (ievent+1)%nSubSample;
    // rest particle counters
    for (Int_t j=0;j<nParticles;j++){
      nTracksPerEventArr[j]=0;
    }
    // reset track counter
    for (Int_t i=0;i<nMaxTracksPerEvent;i++) { trackdEdx[i]=0.; trackSign[i]=0; trackMom[i]=0.; trackEta[i]=0.; trackPhi[i]=0.; trackID[i]=-100; }
    //
    // Find centrality
    Int_t centIndex = randomGen.Uniform(0,9);
    fCent = (centDown[centIndex]+centUp[centIndex])/2.;
    //
    // Find particle multiplicity
    if (fCentSetting==0) {
      prMean  = hCent[centIndex]->GetRandom();
      prBMean = hCent[centIndex]->GetRandom();
    }
    if (fCentSetting==1) {
      prMean  = hCent[centIndex]->GetMean();
      prBMean = hCent[centIndex]->GetMean();
    }
    elMean  = prMean*fBGFraction;
    elBMean = prBMean*fBGFraction;
    //
    // Find particle multiplicity per event
    nTracksPerEventArr[kEl]  = randomGen.Poisson(elMean);
    nTracksPerEventArr[kPr]  = randomGen.Poisson(prMean);
    nTracksPerEventArr[kBEl] = randomGen.Poisson(elBMean);
    nTracksPerEventArr[kBPr] = randomGen.Poisson(prBMean);
    //
    // Track loop
    Int_t trCount=0;
    for (Int_t ipart=0; ipart<nParticles; ipart++){

      for (Int_t i=0; i<nTracksPerEventArr[ipart];i++){

        if (ipart == kEl) {trackdEdx[trCount]  = randomGen.Gaus(elParams[0],elParams[1]);     hParticles[ipart]->Fill(trackdEdx[trCount]); trackID[trCount] = 2;}
        if (ipart == kPr) {trackdEdx[trCount]  = randomGen.Gaus(prParams[0],prParams[1]);     hParticles[ipart]->Fill(trackdEdx[trCount]); trackID[trCount] = 0;}
        if (ipart == kBEl) {trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1])*-1.; hParticles[ipart]->Fill(trackdEdx[trCount]); trackID[trCount] = 2;}
        if (ipart == kBPr) {trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1])*-1.; hParticles[ipart]->Fill(trackdEdx[trCount]); trackID[trCount] = 1;}
        trackSign[trCount]= (ipart<nParticles/2.) ? 1 : -1 ;
        //
        // trackPhi[trCount] = randomGen.Uniform(0,2*TMath::Pi());
        // trackMom[trCount] = randomGen.Exp(0.5);
        // trackEta[trCount] = randomGen.Uniform(-0.8,0.8);
        trackPhi[trCount] = hPhi->GetRandom();
        trackMom[trCount] = hMom1->GetRandom();
        trackEta[trCount] = hEta->GetRandom();
        trCount++;
      }

    }
    //
    // Fill histograms first moms
    for (Int_t i=0;i<nParticles;i++){
      if ( !(switchOffParticle>-1 && i>=switchOffParticle) ) {
        hFirstMoms[i]->Fill(nTracksPerEventArr[i]);
      }
    }

    if(ievent%10000==0) {
      cout << ievent << "  " << trCount;
      cout <<  "  " << nTracksPerEventArr[kEl] <<  "  " << nTracksPerEventArr[kBEl] ;
      cout <<  "  " << nTracksPerEventArr[kPr] <<  "  " << nTracksPerEventArr[kBPr] ;
      cout << endl;
    }
    //
    // -------------------------------------------------------------------------
    //                       Moment calculation part
    // -------------------------------------------------------------------------
    //
    // Initialize moments
    const Int_t nMoments   = 14;
    Float_t genPos=0, genNeg=0,recPos=0, recNeg=0;
    TVectorF fMomNetPrGen(nMoments);
    TVectorF fMomNetPrRec(nMoments);
    for(Int_t i=0;i<nMoments; i++){ fMomNetPrGen[i]=0.; fMomNetPrRec[i]=0.; }
    //
    // Fill track loop
    for (Int_t itr=0;itr<=trCount;itr++){

      if (trackdEdx[itr]==0) continue;
      if ( trackID[itr]==0) genPos++;
      if ( trackID[itr]==1) genNeg++;
      Double_t xxxGenSystScan[3]={Float_t(trackID[itr]),trackMom[itr],trackEta[itr]};
      if (trackSign[itr]>0) fHistPosEffMatrixScanGen->Fill(xxxGenSystScan);
      if (trackSign[itr]<0) fHistNegEffMatrixScanGen->Fill(xxxGenSystScan);
      //
      // Track cuts
      if (AcceptTrack(effSetting,trackMom[itr],trackEta[itr],trackSign[itr],trackID[itr],fCent,ievent))
      {
        if ( trackID[itr]==0) recPos++;
        if ( trackID[itr]==1) recNeg++;
        Double_t xxxRecSystScan[3]={Float_t(trackID[itr]),trackMom[itr],trackEta[itr]};
        if (trackSign[itr]>0) fHistPosEffMatrixScanRec->Fill(xxxRecSystScan);
        if (trackSign[itr]<0) fHistNegEffMatrixScanRec->Fill(xxxRecSystScan);
      }

      //
      outputData->GetFile()->cd();
      *outputData << "tracks" <<
      "effSet="         << fEffSetting <<    // efficiency loss setting
      "centSet="        << fCentSetting <<   // cent setting, with or without volume fluct
      "gid="            << ievent <<         // unique id of an event
      "effAPr="         << fEffAPr <<        // mean eff of antiprotons for a given setting
      "effPr="          << fEffPr <<         // mean eff of protons for a given setting
      "subsample="      << subsampleID <<    // subsample index
      "dEdx="           << trackdEdx[itr] << // simulated dEdx
      "mom="            << trackMom[itr] <<
      "eta="            << trackEta[itr] <<
      "phi="            << trackPhi[itr] <<
      "cutBit="         << cutBit <<
      "sign="           << trackSign[itr] <<
      "part="           << trackID[itr] <<
      "cent="           << fCent <<
      "\n";

    }
    // Net Protons
    fMomNetPrGen[kA]    = genPos;                       fMomNetPrRec[kA]    = recPos;
    fMomNetPrGen[kB]    = genNeg;                       fMomNetPrRec[kB]    = recNeg;
    fMomNetPrGen[kAA]   = genPos*genPos;                fMomNetPrRec[kAA]   = recPos*recPos;
    fMomNetPrGen[kBB]   = genNeg*genNeg;                fMomNetPrRec[kBB]   = recNeg*recNeg;
    fMomNetPrGen[kAB]   = genPos*genNeg;                fMomNetPrRec[kAB]   = recPos*recNeg;
    fMomNetPrGen[kAAA]  = genPos*genPos*genPos;         fMomNetPrRec[kAAA]  = recPos*recPos*recPos;
    fMomNetPrGen[kBBB]  = genNeg*genNeg*genNeg;         fMomNetPrRec[kBBB]  = recNeg*recNeg*recNeg;
    fMomNetPrGen[kAAB]  = genPos*genPos*genNeg;         fMomNetPrRec[kAAB]  = recPos*recPos*recNeg;
    fMomNetPrGen[kBBA]  = genNeg*genNeg*genPos;         fMomNetPrRec[kBBA]  = recNeg*recNeg*recPos;
    fMomNetPrGen[kABBB] = genPos*genNeg*genNeg*genNeg;  fMomNetPrRec[kABBB] = recPos*recNeg*recNeg*recNeg;
    fMomNetPrGen[kAABB] = genPos*genPos*genNeg*genNeg;  fMomNetPrRec[kAABB] = recPos*recPos*recNeg*recNeg;
    fMomNetPrGen[kAAAB] = genPos*genPos*genPos*genNeg;  fMomNetPrRec[kAAAB] = recPos*recPos*recPos*recNeg;
    fMomNetPrGen[kAAAA] = genPos*genPos*genPos*genPos;  fMomNetPrRec[kAAAA] = recPos*recPos*recPos*recPos;
    fMomNetPrGen[kBBBB] = genNeg*genNeg*genNeg*genNeg;  fMomNetPrRec[kBBBB] = recNeg*recNeg*recNeg*recNeg;
    //
    // Dump moments to tree
    outputData->GetFile()->cd();
    *outputData << "events"  <<
    "bgFrac="   << fBGFraction <<     // efficiency loss setting
    "effSet="   << fEffSetting <<    // efficiency loss setting
    "centSet="  << fCentSetting <<   // cent setting, with or without volume fluct
    "gid="      << ievent <<         // unique id of an event
    "subsample="<< subsampleID <<    // subsample index
    "cent="     << fCent <<           // centrality
    "prMean="   << prMean <<         // number of protons for a given event
    "elMean="   << elMean <<         // number of background for a given event
    "netPrMomGen.="   << &fMomNetPrGen <<   // generated level: moments of protons
    "netPrMomRec.="   << &fMomNetPrRec <<   // reconstructed level: moments of protons
    "\n";

    //
  } // event loop ends

  //
  // dump line shape
  Int_t objcounter=0;
  for (Int_t i=0;i<nParticles/2.;i++) {
    hParticles[i] -> SetLineColor(colors[i+1]);
    hFirstMoms[i] -> SetLineColor(colors[i+1]);
    fParticles[i] -> SetLineColor(colors[i+1]);
    hParticles[i+nParticles/2] -> SetLineColor(colors[i+1]);
    hFirstMoms[i+nParticles/2] -> SetLineColor(colors[i+1]);
    fParticles[i+nParticles/2] -> SetLineColor(colors[i+1]);
  }
  //
  // dump line shape
  for (Int_t i=0;i<nParticles;i++) {
    fParticles[i]->SetLineWidth(2);
    fParticles[i]->SetNpx(1000);
    fParticles[i]->SetParameters(hParticles[i]->GetMean(),hParticles[i]->GetRMS());
    if (hParticles[i]) hParticles[i] -> Fit(fParticles[i],"QN");
    LShapesArr_Old[i] = (TF1*)fParticles[i];
  }
  //
  // Create lineshapes for the net-particle case
  fNetParticles[0] = (TF1*)fParticles[1]->Clone(); fNetParticles[0]->SetName("fNetParticles_0");
  fNetParticles[1] = (TF1*)fParticles[3]->Clone(); fNetParticles[1]->SetName("fNetParticles_1");
  fNetParticles[2] = new TF1("fNetParticles_2","particle_0+particle_2",-50,50);
  for (Int_t i=0;i<nNetParticles;i++) {
    fNetParticles[i] -> SetLineColor(colors[i+1]);
    fNetParticles[i] -> SetNpx(1000);
    LShapesArr[i] = (TF1*)fNetParticles[i];
  }

  outputFits->GetFile()->cd();
  LShapesArr.Write("LShapesArr",TObject::kSingleKey);          LShapesArr.Clear("C");
  LShapesArr_Old.Write("LShapesArr_Old",TObject::kSingleKey);  LShapesArr_Old.Clear("C");
  if (fHistPosEffMatrixScanRec) fHistPosEffMatrixScanRec->Write();
  if (fHistNegEffMatrixScanRec) fHistNegEffMatrixScanRec->Write();
  if (fHistPosEffMatrixScanGen) fHistPosEffMatrixScanGen->Write();
  if (fHistNegEffMatrixScanGen) fHistNegEffMatrixScanGen->Write();
  for (Int_t i=0;i<nParticles;i++) { hParticles[i] -> Write(); hFirstMoms[i] -> Write();}
  delete outputFits;
  delete outputData;

}

void PlotdEdx(TString datafile, TString lineShapesFile)
{

  gStyle->SetOptStat(0);
  TFile *f = new TFile(datafile);
  TTree *tree = (TTree*)f->Get("tracks");
  TFile *g = new TFile(lineShapesFile);
  TClonesArray *cloneArrFunc = (TClonesArray*)g->Get("LShapesArr");
  TF1 *fShape[nNetParticles];

  for (Int_t ipart = 0; ipart<nNetParticles; ipart++) {
    TString objName = Form("fNetParticles_%d",ipart);
    fShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
    //
    if (ipart<2){
      Double_t mean  = fShape[ipart]->GetParameter(1);
      Double_t sigma = fShape[ipart]->GetParameter(2);
      fShape[ipart]->SetRange(mean-5*sigma,mean+5*sigma);
    }
  }
  //
  //
  tree->Draw("dEdx>>h(3000,-50,50)");
  TH1D *h0 = (TH1D*)tree->GetHistogram()->Clone();
  h0 -> GetXaxis()->SetTitle("d#it{E}/d#it{x} Signal (arb. units)");
  h0 -> GetYaxis()->SetTitle("Counts (arb. units)");
  h0 -> GetYaxis()->SetTitleOffset(1.1);
  h0 -> SetLineColor(kBlack);
  h0 -> SetMarkerStyle(24);
  h0 -> SetStats(kFALSE);
  h0 -> Rebin(6); h0->Scale(1./6.);
  h0 -> GetXaxis()->SetRangeUser(-40,40);
  h0 -> GetYaxis()->SetRangeUser(0,20000);
  h0 -> GetYaxis()->SetNdivisions(505);
  //
  //
  TLegend *leg = new TLegend(0.8, 0.6, 0.9, 0.9);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(4000);
  leg->AddEntry(h0,        "data" ,"LP");
  leg->AddEntry(fShape[0], "p"       ,"L");
  leg->AddEntry(fShape[1], "#bar{p}" ,"L");
  leg->AddEntry(fShape[2], "BG"      ,"L");
  //
  //
  TLine *line = new TLine(0., 0.,0., 20000.);
  line -> SetLineStyle(2);
  line -> SetLineWidth(2);
  //
  // Plot all in one canvas
  TCanvas *canD = new TCanvas("canD", "canD", 1700, 700);
  h0->Draw("E");
  for (Int_t ipart = 0; ipart<nNetParticles; ipart++) fShape[ipart]->Draw("same");
  leg->Draw("same");
  line->Draw();

}

void PrepareTemplates()
{

  /*

  cd /home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/workdir/input
  aliroot -l
  .L /home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/TOYMC_HighMoments.C+
  PrepareTemplates()

  */

  TTreeSRedirector *outputTemplates = new TTreeSRedirector("Templates.root", "recreate");
  //
  // Get mom, eta and phi dist
  TFile *fDist = TFile::Open("/home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/workdir/input/AnalysisResults_mcGen22.root");
  TTree *tree = (TTree*)fDist->Get("mcGen");
  tree->Draw("eta>>hEta(40,-0.8,0.8)","origin==0&&part==3");
  TH1F *hEta = (TH1F*)tree->GetHistogram()->Clone("hEta");
  tree->Draw("p>>hMom0(80,0.4,2)","origin==0&&part==3");
  TH1F *hMom0 = (TH1F*)tree->GetHistogram()->Clone("hMom0");
  tree->Draw("p>>hMom1(45,0.6,1.5)","origin==0&&part==3");
  TH1F *hMom1 = (TH1F*)tree->GetHistogram()->Clone("hMom1");
  tree->Draw("phi>>hPhi(60,0,6.283185)","origin==0&&part==3");
  TH1F *hPhi = (TH1F*)tree->GetHistogram()->Clone("hPhi");
  cout << "Got mom, eta and phi dist" << endl;

  //
  // Get centrality dist
  TFile *fmcFull = TFile::Open("/home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/workdir/input/AnalysisResults_mcFull34.root");
  TTree *treeFull = (TTree*)fmcFull->Get("mcFull");
  TH1D *hCent[9];
  treeFull->Draw("posGen.fElements[2]>>hCent_0(100,0,100)","orig==0&&abs(etaUp-0.8)<0.001&&abs(pUp-2)<0.001&&cent<5");            hCent[0] = (TH1D*)treeFull->GetHistogram()->Clone("hCent_0");
  treeFull->Draw("posGen.fElements[2]>>hCent_1(100,0,100)","orig==0&&abs(etaUp-0.8)<0.001&&abs(pUp-2)<0.001&&cent<10&&cent>=5");  hCent[1] = (TH1D*)treeFull->GetHistogram()->Clone("hCent_1");
  treeFull->Draw("posGen.fElements[2]>>hCent_2(100,0,100)","orig==0&&abs(etaUp-0.8)<0.001&&abs(pUp-2)<0.001&&cent<20&&cent>=10"); hCent[2] = (TH1D*)treeFull->GetHistogram()->Clone("hCent_2");
  treeFull->Draw("posGen.fElements[2]>>hCent_3(100,0,100)","orig==0&&abs(etaUp-0.8)<0.001&&abs(pUp-2)<0.001&&cent<30&&cent>=20"); hCent[3] = (TH1D*)treeFull->GetHistogram()->Clone("hCent_3");
  treeFull->Draw("posGen.fElements[2]>>hCent_4(100,0,100)","orig==0&&abs(etaUp-0.8)<0.001&&abs(pUp-2)<0.001&&cent<40&&cent>=30"); hCent[4] = (TH1D*)treeFull->GetHistogram()->Clone("hCent_4");
  treeFull->Draw("posGen.fElements[2]>>hCent_5(100,0,100)","orig==0&&abs(etaUp-0.8)<0.001&&abs(pUp-2)<0.001&&cent<50&&cent>=40"); hCent[5] = (TH1D*)treeFull->GetHistogram()->Clone("hCent_5");
  treeFull->Draw("posGen.fElements[2]>>hCent_6(100,0,100)","orig==0&&abs(etaUp-0.8)<0.001&&abs(pUp-2)<0.001&&cent<60&&cent>=50"); hCent[6] = (TH1D*)treeFull->GetHistogram()->Clone("hCent_6");
  treeFull->Draw("posGen.fElements[2]>>hCent_7(100,0,100)","orig==0&&abs(etaUp-0.8)<0.001&&abs(pUp-2)<0.001&&cent<70&&cent>=60"); hCent[7] = (TH1D*)treeFull->GetHistogram()->Clone("hCent_7");
  treeFull->Draw("posGen.fElements[2]>>hCent_8(100,0,100)","orig==0&&abs(etaUp-0.8)<0.001&&abs(pUp-2)<0.001&&cent<80&&cent>=70"); hCent[8] = (TH1D*)treeFull->GetHistogram()->Clone("hCent_8");
  cout << "Got centrality dist" << endl;
  //
  // Get eff matrix dist
  TFile *fEffTPC = TFile::Open("/home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/workdir/input/EffMatrices/EffMatrix_Syst_0_Det_0_Acc_0.60_1.50.root");
  TFile *fEffTOF = TFile::Open("/home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/workdir/input/EffMatrices/EffMatrix_Syst_0_Det_1_Acc_0.60_1.50.root");
  TH1D *hEffTPC_Eta[2];
  TH1D *hEffTPC_Mom[2];
  TH1D *hEffTOF_Eta[2];
  TH1D *hEffTOF_Mom[2];
  hEffTPC_Eta[0] = (TH1D*)fEffTPC->Get("h1DNeg_Pr_5_Cent_0");
  hEffTPC_Eta[1] = (TH1D*)fEffTPC->Get("h1DPos_Pr_5_Cent_0");
  hEffTPC_Mom[0] = (TH1D*)fEffTPC->Get("h1DNeg_Pr_4_Cent_0");
  hEffTPC_Mom[1] = (TH1D*)fEffTPC->Get("h1DPos_Pr_4_Cent_0");
  //
  hEffTOF_Eta[0] = (TH1D*)fEffTOF->Get("h1DNeg_Pr_5_Cent_0");
  hEffTOF_Eta[1] = (TH1D*)fEffTOF->Get("h1DPos_Pr_5_Cent_0");
  hEffTOF_Mom[0] = (TH1D*)fEffTOF->Get("h1DNeg_Pr_4_Cent_0");
  hEffTOF_Mom[1] = (TH1D*)fEffTOF->Get("h1DPos_Pr_4_Cent_0");
  //
  //
  TH2D *hEffEtaMomTPC[2];
  TH2D *hEffEtaMomTOF[2];
  hEffEtaMomTPC[0] = (TH2D*)fEffTPC->Get("h2D_EtaMom_Neg_Pr_Cent_0");
  hEffEtaMomTPC[1] = (TH2D*)fEffTPC->Get("h2D_EtaMom_Pos_Pr_Cent_0");
  hEffEtaMomTOF[0] = (TH2D*)fEffTOF->Get("h2D_EtaMom_Neg_Pr_Cent_0");
  hEffEtaMomTOF[1] = (TH2D*)fEffTOF->Get("h2D_EtaMom_Pos_Pr_Cent_0");
  //
  hMeanEff[0] = new TH1D("meanEffApr","meanEffApr",1000, 0., 1.);
  hMeanEff[1] = new TH1D("meanEffpr" ,"meanEffpr",1000, 0., 1.);
  for (Int_t ibinx=1; ibinx<=hEffEtaMomTPC[0]->GetNbinsX();ibinx++ ){
    for (Int_t ibiny=1; ibiny<=hEffEtaMomTPC[0]->GetNbinsY();ibiny++ ){
      hMeanEff[0]->Fill(hEffEtaMomTPC[0] -> GetBinContent(ibinx,ibiny));
      hMeanEff[1]->Fill(hEffEtaMomTPC[1] -> GetBinContent(ibinx,ibiny));
    }
  }
  cout << "Got eff matrix dist  -->  meanEffs = " << hMeanEff[0]->GetMean() << " -- " << hMeanEff[1]->GetMean() << endl;

  //
  // Centrality template
  const Int_t fnCentBins = 9;
  Float_t xCentBins[] = {0, 5,  10,  20, 30, 40, 50, 60, 70, 80};
  TH1D *fhCent  = new TH1D("hCentTemp","Centrality Bins",fnCentBins ,xCentBins );


  outputTemplates->GetFile()->cd();
  hEta ->SetMinimum(0);
  hMom0->SetMinimum(0);
  hMom1->SetMinimum(0);
  hPhi ->SetMinimum(0);
  //
  hEta ->Write("hEta");
  hMom0->Write("hMom0");
  hMom1->Write("hMom1");
  hPhi ->Write("hPhi");
  //
  for (Int_t i=0; i<9; i++) hCent[i]->Write(Form("hCent_%d",i));
  //
  hEffTPC_Eta[0]->Write("hEffTPC_Eta_APr");
  hEffTPC_Eta[1]->Write("hEffTPC_Eta_Pr");
  hEffTPC_Mom[0]->Write("hEffTPC_Mom_APr");
  hEffTPC_Mom[1]->Write("hEffTPC_Mom_Pr");
  //
  hEffTOF_Eta[0]->Write("hEffTOF_Eta_APr");
  hEffTOF_Eta[1]->Write("hEffTOF_Eta_Pr");
  hEffTOF_Mom[0]->Write("hEffTOF_Mom_APr");
  hEffTOF_Mom[1]->Write("hEffTOF_Mom_Pr");
  //
  hEffEtaMomTPC[0]->Write("hEffEtaMomTPC_APr");
  hEffEtaMomTPC[1]->Write("hEffEtaMomTPC_Pr");
  hEffEtaMomTOF[0]->Write("hEffEtaMomTOF_APr");
  hEffEtaMomTOF[1]->Write("hEffEtaMomTOF_Pr");
  //
  hMeanEff[0]->Write();
  hMeanEff[1]->Write();
  delete outputTemplates;

}

void ReadTemplates()
{

  TFile *fTemplates = TFile::Open("/home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/workdir/input/Templates.root");

  hEta  = (TH1D*)fTemplates->Get("hEta");
  hMom0 = (TH1D*)fTemplates->Get("hMom0");
  hMom1 = (TH1D*)fTemplates->Get("hMom1");
  hPhi  = (TH1D*)fTemplates->Get("hPhi");
  //
  for (Int_t i=0; i<9; i++) hCent[i] = (TH1D*)fTemplates->Get(Form("hCent_%d",i));
  //
  hEffTPC_Eta[0] = (TH1D*)fTemplates->Get("hEffTPC_Eta_APr");
  hEffTPC_Eta[1] = (TH1D*)fTemplates->Get("hEffTPC_Eta_Pr");
  hEffTPC_Mom[0] = (TH1D*)fTemplates->Get("hEffTPC_Mom_APr");
  hEffTPC_Mom[1] = (TH1D*)fTemplates->Get("hEffTPC_Mom_Pr");
  //
  hEffTOF_Eta[0] = (TH1D*)fTemplates->Get("hEffTOF_Eta_APr");
  hEffTOF_Eta[1] = (TH1D*)fTemplates->Get("hEffTOF_Eta_Pr");
  hEffTOF_Mom[0] = (TH1D*)fTemplates->Get("hEffTOF_Mom_APr");
  hEffTOF_Mom[1] = (TH1D*)fTemplates->Get("hEffTOF_Mom_Pr");
  //
  hEffTPC_EtaMom[0] = (TH2D*)fTemplates->Get("hEffEtaMomTPC_APr");
  hEffTPC_EtaMom[1] = (TH2D*)fTemplates->Get("hEffEtaMomTPC_Pr");
  hEffTOF_EtaMom[0] = (TH2D*)fTemplates->Get("hEffEtaMomTOF_APr");
  hEffTOF_EtaMom[1] = (TH2D*)fTemplates->Get("hEffEtaMomTOF_Pr");
  //
  hMeanEff[0] = (TH1D*)fTemplates->Get("meanEffApr");
  hMeanEff[1] = (TH1D*)fTemplates->Get("meanEffpr");

}

Bool_t AcceptTrack(Int_t effSetting, Float_t mom, Float_t eta, Int_t sign, Int_t part, Float_t cent, Float_t event)
{

  /*

  effSetting==4 --> realistic efficiency from templates
  //
  effSetting==0 --> uniform flat efficiency --> equivalent to realistic eff. [0.90, 0.82]
  effSetting==1 --> uniform flat efficiency --> equivalent to realistic eff. [0.70, 0.62]
  effSetting==2 --> uniform flat efficiency --> equivalent to realistic eff. [0.50, 0.42]
  effSetting==3 --> uniform flat efficiency --> equivalent to realistic eff. [0.30, 0.22]
  //


  */
  //
  //
  if (effSetting==0) {
    // fEffAPr = hEffTPC_Mom[0] -> GetBinContent(hEffTPC_Mom[0]->FindBin(mom));
    // fEffPr  = hEffTPC_Mom[1] -> GetBinContent(hEffTPC_Mom[1]->FindBin(mom));
    fEffAPr = hEffTPC_EtaMom[0] -> GetBinContent(hEffTPC_EtaMom[0]->GetXaxis()->FindBin(eta), hEffTPC_EtaMom[0]->GetYaxis()->FindBin(mom));
    fEffPr  = hEffTPC_EtaMom[1] -> GetBinContent(hEffTPC_EtaMom[1]->GetXaxis()->FindBin(eta), hEffTPC_EtaMom[1]->GetYaxis()->FindBin(mom));
    fPrEffCorFactor  = 0.91;    // ??? TODO
    fAPrEffCorFactor = 0.84;
    if (part==0 && randomGen.Uniform(0,1)<=fEffPr) return kTRUE;
    else if (part==1 && randomGen.Uniform(0,1)<fEffAPr) return kTRUE;
    else { fEffAPr=1.; fEffPr=1.; return kFALSE;}
  }
  //
  //  -------------------------------------------------------------------
  //
  if (effSetting==1) {
    fEffPr  = 0.90;  fPrEffCorFactor  = fEffPr;
    fEffAPr = 0.82;  fAPrEffCorFactor = fEffAPr;
    if (part==0 && randomGen.Uniform(0,1)<=fEffPr) return kTRUE;
    else if (part==1 && randomGen.Uniform(0,1)<fEffAPr) return kTRUE;
    else { fEffAPr=1.; fEffPr=1.; return kFALSE;}
  }
  //
  //
  if (effSetting==2) {
    fEffPr  = 0.70;  fPrEffCorFactor  = fEffPr;
    fEffAPr = 0.62;  fAPrEffCorFactor = fEffAPr;
    if (part==0 && randomGen.Uniform(0,1)<=fEffPr) return kTRUE;
    else if (part==1 && randomGen.Uniform(0,1)<fEffAPr) return kTRUE;
    else { fEffAPr=1.; fEffPr=1.; return kFALSE;}
  }
  //
  //
  if (effSetting==3) {
    fEffPr  = 0.50;  fPrEffCorFactor  = fEffPr;
    fEffAPr = 0.42;  fAPrEffCorFactor = fEffAPr;
    if (part==0 && randomGen.Uniform(0,1)<=fEffPr) return kTRUE;
    else if (part==1 && randomGen.Uniform(0,1)<fEffAPr) return kTRUE;
    else { fEffAPr=1.; fEffPr=1.; return kFALSE;}
  }
  //
  //
  if (effSetting==4) {
    fEffPr  = 0.30;  fPrEffCorFactor  = fEffPr;
    fEffAPr = 0.22;  fAPrEffCorFactor = fEffAPr;
    if (part==0 && randomGen.Uniform(0,1)<=fEffPr) return kTRUE;
    else if (part==1 && randomGen.Uniform(0,1)<fEffAPr) return kTRUE;
    else { fEffAPr=1.; fEffPr=1.; return kFALSE;}
  }
  //
  //  -------------------------------------------------------------------
  //
  if (effSetting==5) {
    fEffPr  = 0.90;  fPrEffCorFactor  = fEffPr;
    fEffAPr = 0.82;  fAPrEffCorFactor = fEffAPr;
    // introduce eta gap
    Double_t gapsizePercent = ((fEtaUp-fEtaDown)*0.05)/2.;
    if ( TMath::Abs(eta-0.3)<gapsizePercent && randomGen.Uniform(0,1)>=0.05 )
    {fEffAPr=0.; fEffPr=0.; return kFALSE;}
    //
    if (part==0 && randomGen.Uniform(0,1)<=fEffPr) return kTRUE;
    else if (part==1 && randomGen.Uniform(0,1)<fEffAPr) return kTRUE;
    else { fEffAPr=1.; fEffPr=1.; return kFALSE;}
  }
  //
  //
  if (effSetting==6) {
    fEffPr  = 0.90;  fPrEffCorFactor  = fEffPr;
    fEffAPr = 0.82;  fAPrEffCorFactor = fEffAPr;
    // introduce eta gap
    Double_t gapsizePercent = ((fEtaUp-fEtaDown)*0.10)/2.;
    if ( TMath::Abs(eta-0.3)<gapsizePercent && randomGen.Uniform(0,1)>=0.05 )
    {fEffAPr=0.05; fEffPr=0.05; return kFALSE;}
    //
    if (part==0 && randomGen.Uniform(0,1)<=fEffPr) return kTRUE;
    else if (part==1 && randomGen.Uniform(0,1)<fEffAPr) return kTRUE;
    else { fEffAPr=1.; fEffPr=1.; return kFALSE;}
  }
  //
  //
  if (effSetting==7) {
    fEffPr  = 0.90;  fPrEffCorFactor  = fEffPr;
    fEffAPr = 0.82;  fAPrEffCorFactor = fEffAPr;
    // introduce eta gap
    Double_t gapsizePercent = ((fEtaUp-fEtaDown)*0.15)/2.;
    if ( TMath::Abs(eta-0.3)<gapsizePercent && randomGen.Uniform(0,1)>=0.05 )
    {fEffAPr=0.05; fEffPr=0.05; return kFALSE;}
    //
    if (part==0 && randomGen.Uniform(0,1)<=fEffPr) return kTRUE;
    else if (part==1 && randomGen.Uniform(0,1)<fEffAPr) return kTRUE;
    else { fEffAPr=1.; fEffPr=1.; return kFALSE;}
  }
  //
  //
  if (effSetting==8) {
    fEffPr  = 0.90;  fPrEffCorFactor  = fEffPr;
    fEffAPr = 0.82;  fAPrEffCorFactor = fEffAPr;
    // introduce eta gap
    Double_t gapsizePercent = ((fEtaUp-fEtaDown)*0.20)/2.;
    if ( TMath::Abs(eta-0.3)<gapsizePercent && randomGen.Uniform(0,1)>=0.05 )
    {fEffAPr=0.05; fEffPr=0.05; return kFALSE;}
    //
    if (part==0 && randomGen.Uniform(0,1)<=fEffPr) return kTRUE;
    else if (part==1 && randomGen.Uniform(0,1)<fEffAPr) return kTRUE;
    else { fEffAPr=1.; fEffPr=1.; return kFALSE;}
  }

}

void SetDataTreeAliases(TTree *tree)
{

  TString mcRec = "Rec";
  TString mcGen = "Gen";
  //
  tree->SetAlias("a"  ,Form("netPrMom%s.fElements[0]",mcRec.Data()));   // <a> --> mean proton
  tree->SetAlias("b"  ,Form("netPrMom%s.fElements[1]",mcRec.Data()));   // <b> --> mean antiproton
  tree->SetAlias("a2" ,Form("netPrMom%s.fElements[2]",mcRec.Data()));  // <aa>
  tree->SetAlias("b2" ,Form("netPrMom%s.fElements[3]",mcRec.Data()));  // <bb>
  tree->SetAlias("ab" ,Form("netPrMom%s.fElements[4]",mcRec.Data()));  // <ab>
  tree->SetAlias("a3" ,Form("netPrMom%s.fElements[5]",mcRec.Data()));  // <aaa>
  tree->SetAlias("b3" ,Form("netPrMom%s.fElements[6]",mcRec.Data()));  // <bbb>
  tree->SetAlias("a2b",Form("netPrMom%s.fElements[7]",mcRec.Data()));  // <aab>
  tree->SetAlias("b2a",Form("netPrMom%s.fElements[8]",mcRec.Data()));  // <bba>
  //
  tree->SetAlias("c"  ,Form("netPrMom%s.fElements[0]",mcGen.Data()));   // <a> --> mean proton
  tree->SetAlias("d"  ,Form("netPrMom%s.fElements[1]",mcGen.Data()));   // <b> --> mean antiproton
  tree->SetAlias("c2" ,Form("netPrMom%s.fElements[2]",mcGen.Data()));  // <aa>
  tree->SetAlias("d2" ,Form("netPrMom%s.fElements[3]",mcGen.Data()));  // <bb>
  tree->SetAlias("cd" ,Form("netPrMom%s.fElements[4]",mcGen.Data()));  // <ab>
  tree->SetAlias("c3" ,Form("netPrMom%s.fElements[5]",mcGen.Data()));  // <aaa>
  tree->SetAlias("d3" ,Form("netPrMom%s.fElements[6]",mcGen.Data()));  // <bbb>
  tree->SetAlias("c2d",Form("netPrMom%s.fElements[7]",mcGen.Data()));  // <aab>
  tree->SetAlias("d2c",Form("netPrMom%s.fElements[8]",mcGen.Data()));  // <bba>

}

void ApplyEffCorrectionAllMoments(TString inputDataFile, Int_t nEvents, Float_t bgFraction, Int_t effSetting,  Int_t centSetting )
{

  /*

  cutIndex  = [0,10] --> index to get the corresponding eff matrix

  cd /home/marsland/Desktop/QM19_tmp/results
  aliroot -l
  .L /home/marsland/Desktop/ubuntu_desktop/TOYMC_HighMoments/toyMC/TOYMC_HighMoments.C+
  ApplyEffCorrectionAllMoments("DataTree_Net_500000_BG_0.10_eff_2_VF_1.root",500000, 0.1, 2, 1)
  ApplyEffCorrectionAllMoments("DataTree_Net_500000_BG_0.10_eff_2_VF_0.root",500000, 0.1, 2, 0)

  TFile *f = TFile::Open("Moments.root");
  TTree *tree = (TTree*)f->Get("moms");
  tree->SetMarkerStyle(20)
  */
  fBGFraction  = bgFraction;
  fNevents     = nEvents;
  fEffSetting  = effSetting;
  fCentSetting = centSetting;
  TTreeSRedirector *streamCumsCent = new TTreeSRedirector(Form("Moments_%d_BG_%3.2f_eff_%d_VF_%d.root",nEvents,bgFraction,effSetting,centSetting),"recreate");
  TFile *fMC_tiden0  = new TFile(inputDataFile);
  TTree *tree = (TTree*)fMC_tiden0->Get("events");
  SetDataTreeAliases(tree);
  tree->SetMarkerStyle(20);
  tree->SetMarkerColor(colors[0]);
  tree->SetLineColor(colors[0]);
  //
  for (Int_t isamp=0; isamp<nSubSample; isamp++){
    //
    Int_t nBins    = 9;
    TString plotVS = "cent";
    TCut accCut = "";
    TCut corrCut = "";    // only for Tiden tree
    TCut sampCut = Form("subsample==%d",isamp);
    // nBins = tree->Draw("a:cent",corrCut&&accCut&&sampCut,"prof");
    cout << " nBins " << nBins << "   isample = " << isamp << endl;
    //
    // Read tree information into graphs
    tree->Draw(Form("a:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");   TProfile *prof_a = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_a   = ConvertProfToGraph(prof_a); gr_a ->SetName("gr_a");
    tree->Draw(Form("b:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");   TProfile *prof_b = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_b   = ConvertProfToGraph(prof_b); gr_b ->SetName("gr_b");
    tree->Draw(Form("a2:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_a2 = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_a2   = ConvertProfToGraph(prof_a2); gr_a2 ->SetName("gr_a2");
    tree->Draw(Form("b2:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_b2 = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_b2   = ConvertProfToGraph(prof_b2); gr_b2 ->SetName("gr_b2");
    tree->Draw(Form("ab:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_ab = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_ab   = ConvertProfToGraph(prof_ab); gr_ab ->SetName("gr_ab");
    tree->Draw(Form("a3:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_a3 = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_a3   = ConvertProfToGraph(prof_a3); gr_a3 ->SetName("gr_a3");
    tree->Draw(Form("b3:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_b3 = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_b3   = ConvertProfToGraph(prof_b3); gr_b3 ->SetName("gr_b3");
    tree->Draw(Form("a2b:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof"); TProfile *prof_a2b = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_a2b   = ConvertProfToGraph(prof_a2b); gr_a2b ->SetName("gr_a2b");
    tree->Draw(Form("b2a:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof"); TProfile *prof_b2a = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_b2a   = ConvertProfToGraph(prof_b2a); gr_b2a ->SetName("gr_b2a");
    //
    //
    // Read tree information into graphs
    tree->Draw(Form("c:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");   TProfile *prof_c = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_c   = ConvertProfToGraph(prof_c); gr_c ->SetName("gr_c");
    tree->Draw(Form("d:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");   TProfile *prof_d = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_d   = ConvertProfToGraph(prof_d); gr_d ->SetName("gr_d");
    tree->Draw(Form("c2:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_c2 = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_c2   = ConvertProfToGraph(prof_c2); gr_c2 ->SetName("gr_c2");
    tree->Draw(Form("d2:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_d2 = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_d2   = ConvertProfToGraph(prof_d2); gr_d2 ->SetName("gr_d2");
    tree->Draw(Form("cd:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_cd = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_cd   = ConvertProfToGraph(prof_cd); gr_cd ->SetName("gr_cd");
    tree->Draw(Form("c3:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_c3 = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_c3   = ConvertProfToGraph(prof_c3); gr_c3 ->SetName("gr_c3");
    tree->Draw(Form("d3:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof");  TProfile *prof_d3 = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_d3   = ConvertProfToGraph(prof_d3); gr_d3 ->SetName("gr_d3");
    tree->Draw(Form("c2d:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof"); TProfile *prof_c2d = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_c2d   = ConvertProfToGraph(prof_c2d); gr_c2d ->SetName("gr_c2d");
    tree->Draw(Form("d2c:%s",plotVS.Data()),corrCut&&accCut&&sampCut,"prof"); TProfile *prof_d2c = (TProfile*)tree->GetHistogram()->Clone(); TGraphErrors *gr_d2c   = ConvertProfToGraph(prof_d2c); gr_d2c ->SetName("gr_d2c");
    // cout << " moments are read " << endl;
    //
    // loop over cent or eta bins
    for (Int_t ibin=0; ibin<nBins; ibin++)
    {
      //
      Double_t e1 = fPrEffCorFactor;
      Double_t e2 = fAPrEffCorFactor;
      //
      // Reconstructed moments
      Double_t a   = gr_a->GetY()[ibin];
      Double_t b   = gr_b->GetY()[ibin];
      Double_t a2  = gr_a2->GetY()[ibin];
      Double_t b2  = gr_b2->GetY()[ibin];
      Double_t ab  = gr_ab->GetY()[ibin];
      Double_t a3  = gr_a3->GetY()[ibin];
      Double_t b3  = gr_b3->GetY()[ibin];
      Double_t a2b = gr_a2b->GetY()[ibin];
      Double_t b2a = gr_b2a->GetY()[ibin];
      Double_t m1 = a-b;
      Double_t m2 = a2+b2-2*ab;
      Double_t m3 = a3-3*a2b+3*b2a-b3;
      Double_t C2 = m2-m1*m1;
      Double_t C3 = m3 - 3*m2*m1 + 2*m1*m1*m1;
      Double_t S  = a+b;
      Double_t C2a = a2-a*a;
      Double_t C2b = b2-b*b;
      Double_t C3a = a3 - 3*a2*a + 2*a*a*a;
      Double_t C3b = b3 - 3*b2*b + 2*b*b*b;
      //
      // Corrected ones
      //
      Double_t A = a/e1;
      Double_t B = b/e2;
      Double_t AB = ab/(e1*e2);
      Double_t A2 = a2/(e1*e1) - a/(e1*e1) + a/e1;
      Double_t B2 = b2/(e2*e2) - b/(e2*e2) + b/e2;
      Double_t A3 = a3/(e1*e1*e1)-(3*a2)/(e1*e1*e1)+(2*a)/(e1*e1*e1)+3*A2-2*A;
      Double_t B3 = b3/(e2*e2*e2)-(3*b2)/(e2*e2*e2)+(2*b)/(e2*e2*e2)+3*B2-2*B;
      Double_t A2B = a2b/(e1*e1*e2)-ab/(e1*e1*e2)+ab/(e1*e2);
      Double_t B2A = b2a/(e2*e2*e1)-ab/(e2*e2*e1)+ab/(e1*e2);
      Double_t M1 = A-B;
      Double_t M2 = A2+B2-2*AB;
      Double_t M3 = A3-3*A2B+3*B2A-B3;
      Double_t CC2 = M2-M1*M1;
      Double_t CC3 = M3 - 3*M2*M1 + 2*M1*M1*M1;
      Double_t SS = A+B;
      Double_t C2A = A2-A*A;
      Double_t C2B = B2-B*B;
      Double_t C3A = A3 - 3*A2*A + 2*A*A*A;
      Double_t C3B = B3 - 3*B2*B + 2*B*B*B;
      //
      // Reconstructed moments
      Double_t c   = gr_c->GetY()[ibin];
      Double_t d   = gr_d->GetY()[ibin];
      Double_t c2  = gr_c2->GetY()[ibin];
      Double_t d2  = gr_d2->GetY()[ibin];
      Double_t cd  = gr_cd->GetY()[ibin];
      Double_t c3  = gr_c3->GetY()[ibin];
      Double_t d3  = gr_d3->GetY()[ibin];
      Double_t c2d = gr_c2d->GetY()[ibin];
      Double_t d2c = gr_d2c->GetY()[ibin];
      Double_t n1 = c-d;
      Double_t n2 = c2+d2-2*cd;
      Double_t n3 = c3-3*c2d+3*d2c-d3;
      Double_t K2 = n2-n1*n1;
      Double_t K3 = n3 - 3*n2*n1 + 2*n1*n1*n1;
      Double_t R  = c+d;
      Double_t K2c = c2-c*c;
      Double_t K2d = d2-d*d;
      Double_t K3c = c3 - 3*c2*c + 2*c*c*c;
      Double_t K3d = d3 - 3*d2*d + 2*d*d*d;
      //
      // tracks tree
      streamCumsCent->GetFile()->cd();
      (*streamCumsCent)<<"moms"<<
      "bgFrac="            << fBGFraction <<    // efficiency loss setting
      "effSet="            << fEffSetting <<    // efficiency loss setting
      "centSet="           << fCentSetting <<   // cent setting, with or without volume fluct
      "bin="               << ibin                   <<
      "cent="              << gr_a->GetX()[ibin]     <<
      "subsample="         << isamp                  <<
      "e1="                << e1                     <<
      "e2="                << e2                     <<
      //
      "a="                 << a         <<
      "b="                 << b         <<
      "a2="                << a2        <<
      "b2="                << b2        <<
      "ab="                << ab        <<
      "a3="                << a3        <<
      "b3="                << b3        <<
      "a2b="               << a2b       <<
      "b2a="               << b2a       <<
      "m1="                << m1         <<
      "m2="                << m2         <<
      "m3="                << m3        <<
      "C2="                << C2         <<
      "C3="                << C3         <<
      "S="                 << S        <<
      "C2a="               << C2a         <<
      "C2b="               << C2b         <<
      "C3a="               << C3a        <<
      "C3b="               << C3b        <<
      //
      //  Generated ones
      "c="                 << c         <<
      "d="                 << d         <<
      "c2="                << c2        <<
      "d2="                << d2        <<
      "cd="                << cd        <<
      "c3="                << c3        <<
      "d3="                << d3        <<
      "c2d="               << c2d       <<
      "d2c="               << d2c       <<
      "n1="                << n1         <<
      "n2="                << n2         <<
      "n3="                << n3        <<
      "K2="                << K2         <<
      "K3="                << K3         <<
      "R="                 << R        <<
      "K2c="               << K2c         <<
      "K2d="               << K2d         <<
      "K3c="               << K3c        <<
      "K3d="               << K3d        <<
      //
      // Corrected ones
      //
      "A="                 << A         <<
      "B="                 << B         <<
      "AB="                << AB        <<
      "A2="                << A2         <<
      "B2="                << B2         <<
      "A3="                << A3        <<
      "B3="                << B3         <<
      "A2B="               << A2B         <<
      "B2A="               << B2A        <<
      "M1="                << M1        <<
      "M2="                << M2         <<
      "M3="                << M3         <<
      "CC2="               << CC2        <<
      "CC3="               << CC3         <<
      "SS="                << SS        <<
      "C2A="               << C2A         <<
      "C2B="               << C2B         <<
      "C3A="               << C3A        <<
      "C3B="               << C3B        <<
      "\n";

    }

    delete gr_a;
    delete gr_b;
    delete gr_a2;
    delete gr_b2;
    delete gr_ab;
    delete gr_a3;
    delete gr_b3;
    delete gr_a2b;
    delete gr_b2a;
    delete gr_c;
    delete gr_d;
    delete gr_c2;
    delete gr_d2;
    delete gr_cd;
    delete gr_c3;
    delete gr_d3;
    delete gr_c2d;
    delete gr_d2c;


  }

  delete streamCumsCent;

}

TGraphErrors *ConvertProfToGraph(TProfile *hh)
{

  /*
  TFile f("/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN1/10h_MCclosure_cutON_p_02122018/TIdentity/TIdenResults_0_0_plotVS_ptot_cutON_ptot_0.1_3.1/AllMoments_TIdenResults_0_0_plotVS_ptot_cutON_ptot_0.1_3.1.root")
  momTree->Draw("fMoments1st.fElements[0]:centBin","","profgoff");
  TProfile *h = (TProfile*)momTree->GetHistogram()->Clone();
  for (Int_t i=0; i<100;i++) if (h->GetBinContent(i)>0) cout << h->GetBinCenter(i) << "   " << h->GetBinContent(i) << "  " << h->GetBinError(i) << endl;
  */

  Int_t nBinsProf = hh->GetNbinsX();
  TString pName = hh->GetName();
  Int_t nCentBins=0;
  for (Int_t i=0; i<nBinsProf; i++){
    if (hh->GetBinError(i) && hh->GetBinContent(i) ) nCentBins++;
  }
  Double_t xx[nCentBins] = {0.};
  Double_t yy[nCentBins] = {0.};
  Double_t exx[nCentBins] = {0.};
  Double_t eyy[nCentBins] = {0.};

  Int_t binCount=0;
  for (Int_t i=0; i<nBinsProf; i++){

    if (hh->GetBinError(i) && hh->GetBinContent(i) ) {
      // cout << hh->GetBinCenter(i) << "   " << hh->GetBinContent(i) << "  " << hh->GetBinError(i) << endl;
      xx[binCount]  = hh->GetBinCenter(i);
      yy[binCount]  = hh->GetBinContent(i);
      eyy[binCount] = hh->GetBinError(i);
      exx[binCount] = 0.;
      binCount++;
    }

  }

  TGraphErrors * gr = new TGraphErrors(nCentBins,xx,yy,exx,eyy);
  gr->SetName(Form("gr_%s",pName.Data()));
  return gr;

}

void InitHists()
{

  const Int_t ndimScan=3;
  Int_t nbinsScan[ndimScan]   = {2, 200      ,200       };
  Double_t xminScan[ndimScan] = {0, fMomDown ,fEtaDown };
  Double_t xmaxScan[ndimScan] = {2, fMomUp   ,fEtaUp   };
  fHistPosEffMatrixScanRec  =new THnF("fHistPosEffMatrixScanRec","fHistPosEffMatrixScanRec",ndimScan, nbinsScan,xminScan,xmaxScan);
  fHistNegEffMatrixScanRec  =new THnF("fHistNegEffMatrixScanRec","fHistNegEffMatrixScanRec",ndimScan, nbinsScan,xminScan,xmaxScan);
  fHistPosEffMatrixScanGen  =new THnF("fHistPosEffMatrixScanGen","fHistPosEffMatrixScanGen",ndimScan, nbinsScan,xminScan,xmaxScan);
  fHistNegEffMatrixScanGen  =new THnF("fHistNegEffMatrixScanGen","fHistNegEffMatrixScanGen",ndimScan, nbinsScan,xminScan,xmaxScan);
  TString axisNameEffScan[ndimScan]  = {"particle"      ,"momentum"                ,"eta" };
  TString axisTitleEffScan[ndimScan] = {"particle type" ,"#it{p}_{T} (GeV/#it{c})" ,"#eta"};
  for (Int_t iEff=0; iEff<ndimScan;iEff++){
    fHistPosEffMatrixScanRec->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistPosEffMatrixScanRec->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
    fHistNegEffMatrixScanRec->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistNegEffMatrixScanRec->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
    fHistPosEffMatrixScanGen->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistPosEffMatrixScanGen->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
    fHistNegEffMatrixScanGen->GetAxis(iEff)->SetName(axisNameEffScan[iEff]);  fHistNegEffMatrixScanGen->GetAxis(iEff)->SetTitle(axisTitleEffScan[iEff]);
  }

}
