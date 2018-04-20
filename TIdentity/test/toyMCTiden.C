#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TVectorF.h"
#include "TClonesArray.h"
#include "TTreeStream.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using std::cout;
using std::setw;


void PrintMoments();
void PlotdEdx(Int_t sign);
void PlotWs();

/*

cd /u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/
aliroot -l
.L toyMCTiden.C+
toyMCTiden()
PlotdEdx(0)

PlotWs()

*/
Int_t switchOffParticle=-100;   // default = -100 for 4 particles
Int_t fUsedSign = -1;
const Int_t nParticles=5;
const Int_t nSignBins =3;
const Int_t nMoments = 19;
Int_t elMean = 6, piMean=10, kaMean=6, prMean=8, deMean=4;
Int_t nTracksPerEventArr[nSignBins][nParticles]={0};
Double_t elParams[]={8.,1.5};
Double_t piParams[]={4.,1.};
Double_t kaParams[]={11.,1.5};
Double_t prParams[]={17.,1.6};
Double_t deParams[]={25.,1.8};

Float_t trackdEdx[1000]={0.};
Int_t trackSign[1000]={0};

TH1D *hParticles[nSignBins][nParticles];
TH1D *hFirstMoms[nSignBins][nParticles];
TF1 *fParticles[nSignBins][nParticles];
UInt_t cutBit=0;
Int_t sign=0;
Int_t signBin=-100;
ULong64_t nEvents=500000;
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray funcLineShapesCArr("TF1",50000);
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kDe=4,
  kElEl=5,kPiPi=6,kKaKa=7,kPrPr=8,kDeDe=9,
  kElPi=10,kElKa=11,kElPr=12,kElDe=13,
  kPiKa=14,kPiPr=15,kPiDe=16,
  kKaPr=17,kKaDe=18,};
TString momNames[nMoments] = {"kEl","kPi","kKa","kPr","kDe",
  "kElEl","kPiPi","kKaKa","kPrPr","kDeDe",
  "kElPi","kElKa","kElPr","kElDe",
  "kPiKa","kPiPr","kPiDe",
  "kKaPr","kKaDe"};
  enum momentTypeUnlikeSign {
    kPiPosPiNeg=0,
    kPiPosKaNeg=1,
    kPiPosPrNeg=2,
    kKaPosPiNeg=3,
    kKaPosKaNeg=4,
    kKaPosPrNeg=5,
    kPrPosPiNeg=6,
    kPrPosKaNeg=7,
    kPrPosPrNeg=8,
    kLaPosLaNeg=9,
    kChPosChNeg=10,
    kBaPosBaNeg=11,
  };
//
//
// ================================================================================================
//
//
void toyMCTiden(){

  TTreeSRedirector *outputData = new TTreeSRedirector("dataTree.root", "recreate");
  TTreeSRedirector *outputFits = new TTreeSRedirector("LineShapes.root", "recreate");
  funcLineShapesCArr.SetOwner(kTRUE);
  for (Int_t j=0;j<nSignBins;j++) {
    for (Int_t i=0;i<nParticles;i++) {
      hParticles[j][i] = new TH1D(Form("hist_%d_bin_%d",i,j),Form("hist_%d_bin_%d",i,j),1500,0.,50.);
      hFirstMoms[j][i] = new TH1D(Form("firstMom_%d_bin_%d",i,j),Form("firstMom_%d_bin_%d",i,j),1500,0.,50.);
      fParticles[j][i] = new TF1(Form("particle_%d_bin_%d",i,j),"gaus",0,50);
    }
  }
  if (fUsedSign==-1) signBin = 0;
  if (fUsedSign== 1) signBin = 1;
  if (fUsedSign== 0) signBin = 2;

  Double_t pidVar = 0;
  TRandom randomGen;
  TVectorF recMoments(nMoments);
  for(Int_t i=0;i<nMoments; i++) recMoments[i]=0.;

  for (ULong64_t ievent=0;ievent<nEvents;ievent++){  // event loop

    // rest particle counters
    for (Int_t i=0;i<nSignBins;i++){
      for (Int_t j=0;j<nParticles;j++){
        nTracksPerEventArr[i][j]=0;
      }
    }
    // reset track counter
    for (Int_t i=0;i<1000;i++) { trackdEdx[i]=0.; trackSign[i]=0; }
    //
    Float_t cent=randomGen.Uniform(0,10);
    Int_t trCount=0;
    nTracksPerEventArr[2][kEl] = randomGen.Poisson(elMean);
    nTracksPerEventArr[2][kPi] = randomGen.Poisson(piMean);
    nTracksPerEventArr[2][kKa] = randomGen.Poisson(kaMean);
    nTracksPerEventArr[2][kPr] = randomGen.Poisson(prMean);
    nTracksPerEventArr[2][kDe] = randomGen.Poisson(deMean);
    //
    // Switch of some particles
    for (Int_t i=0;i<nParticles;i++) {
      if ( switchOffParticle>-1 && i>=switchOffParticle ) nTracksPerEventArr[2][i]=0;
    }
    //
    // Generate electron dEdx
    for (Int_t i=0; i<nTracksPerEventArr[2][kEl];i++){
      trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]);
      hParticles[2][kEl]->Fill(trackdEdx[trCount]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kEl]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kEl]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kEl]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kEl]++; trackSign[trCount]=1; }
      trCount++;
    }
    //
    // Generate pion dEdx
    for (Int_t i=0; i<nTracksPerEventArr[2][kPi];i++){
      trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]);
      hParticles[2][kPi]->Fill(trackdEdx[trCount]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kPi]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kPi]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kPi]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kPi]++; trackSign[trCount]=1; }
      trCount++;
    }
    //
    // Generate kaon dEdx
    for (Int_t i=0; i<nTracksPerEventArr[2][kKa];i++){
      trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]);
      hParticles[2][kKa]->Fill(trackdEdx[trCount]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kKa]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kKa]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kKa]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kKa]++; trackSign[trCount]=1; }
      trCount++;
    }
    //
    // Generate proton dEdx
    for (Int_t i=0; i<nTracksPerEventArr[2][kPr];i++){
      trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]);
      hParticles[2][kPr]->Fill(trackdEdx[trCount]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kPr]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kPr]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kPr]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kPr]++; trackSign[trCount]=1; }
      trCount++;
    }
    // Generate proton dEdx
    for (Int_t i=0; i<nTracksPerEventArr[2][kDe];i++){
      trackdEdx[trCount] = randomGen.Gaus(deParams[0],deParams[1]);
      hParticles[2][kDe]->Fill(trackdEdx[trCount]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kDe]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kDe]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kDe]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kDe]++; trackSign[trCount]=1; }
      trCount++;
    }
    //
    // Fill histograms first moms
    for (Int_t i=0;i<nParticles;i++){
      if ( !(switchOffParticle>-1 && i>=switchOffParticle) ) {
        hFirstMoms[0][i]->Fill(nTracksPerEventArr[0][i]);
        hFirstMoms[1][i]->Fill(nTracksPerEventArr[1][i]);
        hFirstMoms[2][i]->Fill(nTracksPerEventArr[2][i]);
      }
    }

    if(ievent%100000==0) {
      cout << ievent << "  " << trCount << " " <<  nTracksPerEventArr[signBin][kEl];
      cout <<  "  " << nTracksPerEventArr[signBin][kPi] ;
      cout <<  "  " << nTracksPerEventArr[signBin][kKa] ;
      cout <<  "  " << nTracksPerEventArr[signBin][kPr] ;
      cout <<  "  " << nTracksPerEventArr[signBin][kDe] << endl;
    }
    //
    // Fill data tree
    for (Int_t itr=0;itr<1000;itr++){
      outputData->GetFile()->cd();
      if(trackdEdx[itr]>0){
        *outputData << "tracks"
        << "gid="   << ievent <<
        "dEdx="     << trackdEdx[itr] <<
        "cutBit="   << cutBit <<
        "sign="     << trackSign[itr] <<
        "cent="     << cent <<
        "\n";
      }
    }
    //
    // Calculate Moments
    recMoments[kEl]=Float_t(nTracksPerEventArr[signBin][kEl]);
    recMoments[kPi]=Float_t(nTracksPerEventArr[signBin][kPi]);
    recMoments[kKa]=Float_t(nTracksPerEventArr[signBin][kKa]);
    recMoments[kPr]=Float_t(nTracksPerEventArr[signBin][kPr]);
    recMoments[kDe]=Float_t(nTracksPerEventArr[signBin][kDe]);

    recMoments[kElEl]=recMoments[kEl]*recMoments[kEl];
    recMoments[kPiPi]=recMoments[kPi]*recMoments[kPi];
    recMoments[kKaKa]=recMoments[kKa]*recMoments[kKa];
    recMoments[kPrPr]=recMoments[kPr]*recMoments[kPr];
    recMoments[kDeDe]=recMoments[kDe]*recMoments[kDe];

    recMoments[kElPi]=recMoments[kEl]*recMoments[kPi];
    recMoments[kElKa]=recMoments[kEl]*recMoments[kKa];
    recMoments[kElPr]=recMoments[kEl]*recMoments[kPr];
    recMoments[kElDe]=recMoments[kEl]*recMoments[kDe];

    recMoments[kPiKa]=recMoments[kPi]*recMoments[kKa];
    recMoments[kPiPr]=recMoments[kPi]*recMoments[kPr];
    recMoments[kPiDe]=recMoments[kPi]*recMoments[kDe];

    recMoments[kKaPr]=recMoments[kKa]*recMoments[kPr];
    recMoments[kKaDe]=recMoments[kKa]*recMoments[kDe];
    //
    // Dump moments to tree
    outputData->GetFile()->cd();
    *outputData << "events"
    << "gid=" << ievent <<
    "moment.=" << &recMoments <<             // second moments for particle+antiparticle
    "\n";

    //
  } // event loop ends
  //
  // dump line shape
  Int_t objcounter=0;
  for (Int_t j=0;j<nSignBins;j++) {
    for (Int_t i=0;i<nParticles;i++) {
      hParticles[j][i] -> SetLineColor(colors[i+1]);
      hFirstMoms[j][i] -> SetLineColor(colors[i+1]);
      fParticles[j][i] -> SetLineColor(colors[i+1]);
      fParticles[j][i] -> SetLineWidth(2);
      fParticles[j][i]->SetNpx(1000);
      fParticles[j][i]->SetParameters(hParticles[j][i]->GetMean(),hParticles[j][i]->GetRMS());
      if (hParticles[j][i]) hParticles[j][i] -> Fit(fParticles[j][i],"QN");
      funcLineShapesCArr[objcounter] = (TF1*)fParticles[j][i];
      objcounter++;
    }
  }
  outputFits->GetFile()->cd();
  funcLineShapesCArr . Write("funcLineShapesCArr",TObject::kSingleKey);
  funcLineShapesCArr.Clear("C");
  for (Int_t j=0;j<nSignBins;j++){
    for (Int_t i=0;i<nParticles;i++){
      hParticles[j][i] -> Write();
    }
  }
  for (Int_t j=0;j<nSignBins;j++){
    for (Int_t i=0;i<nParticles;i++){
      hFirstMoms[j][i] -> Write();
    }
  }
  // for (Int_t i=0;i<4;i++) fParticles[i] -> Write();
  delete outputFits;
  delete outputData;

  PrintMoments();

}

void PrintMoments(){

  cout << " ================================================================================== "<<endl;
  cout << " ================================================================================== "<<endl;
  TFile *f = new TFile("dataTree.root");
  TTree *tree = (TTree*)f->Get("events");
  TH1D *h = new TH1D("hGen","hGen",nMoments,0.,nMoments);
  TH1D *htemp[nMoments];
  for (Int_t i=1;i<nMoments+1;i++) {
    tree->Draw(Form("moment.fElements[%d]",i-1),"","goff");
    htemp[i-1] = (TH1D*)tree->GetHistogram()->Clone();
    htemp[i-1]->SetName(Form("htmp_%d",i-1));
    cout << momNames[i-1] << "  " << htemp[i-1]->GetMean() << endl;
    h->SetBinContent(i,htemp[i-1]->GetMean());
  }
  cout << " ================================================================================== "<<endl;
  cout << " ================================================================================== "<<endl;
  for (Int_t i=1;i<nMoments+1;i++) h->GetXaxis()->SetBinLabel(i,momNames[i-1]);
  TFile *outFile = new TFile("toyMC_Moments_Gen.root","recreate");
  h->Write();
  outFile -> Close();
  delete outFile; //yeni eklave etdim.

}

void PlotdEdx(Int_t sign){


  TFile *f = new TFile("dataTree.root");
  TTree *tree = (TTree*)f->Get("tracks");
  TFile *g = new TFile("LineShapes.root");
  TClonesArray *cloneArrFunc = (TClonesArray*)g->Get("funcLineShapesCArr");
  TF1 *fShape[nParticles];

  Int_t binIndex;
  if (sign==-1) binIndex = 0;
  if (sign== 1) binIndex = 1;
  if (sign== 0) binIndex = 2;

  for (Int_t ipart = 0; ipart<nParticles; ipart++) {
    TString objName = Form("particle_%d_bin_%d",ipart,binIndex);
    fShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
  }
  TString cutSign = (TMath::Abs(sign)==1) ? Form("sign==%d",sign) : "";
  tree->Draw("dEdx>>h(1500,0,50)",cutSign);
  for (Int_t ipart = 0; ipart<nParticles; ipart++){
    fShape[ipart]->Draw("same");
  }

}

void PlotWs(){


  TFile *f = new TFile("dataTree_10M.root");
  TTree *tree = (TTree*)f->Get("tracks");
  TFile *g = new TFile("LineShapes_10M.root");
  TFile *h = new TFile("TIdenDebug_10M.root");
  TClonesArray *cloneArrFunc = (TClonesArray*)g->Get("funcLineShapesCArr");
  TF1 *fShape[4];
  TH1D *hW[4];
  TH1D *hN[4];

  for (Int_t ipart = 0; ipart<4; ipart++) {
    TString objName = Form("particle_%d",ipart);
    fShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
    hN[ipart]     = (TH1D*)g->Get(Form("mom1_%d",ipart));
    hW[ipart]     = (TH1D*)h->Get(Form("hW_%d",ipart));

    hW[ipart] -> SetLineColor(colors[ipart+1]);
    hN[ipart] -> SetLineColor(colors[ipart+1]);
  }
  //
  // tree->Draw("dEdx>>h(900,0,30)");
  // for (Int_t ipart = 0; ipart<4; ipart++){
  //   fShape[ipart]->Draw("same");
  // }

  TCanvas *can = new TCanvas("can", "can", 1200, 600);
  can->Divide(4,2);
  can->cd(1); hW[0] -> Draw();
  can->cd(5); hN[0] -> Draw();

  can->cd(2); hW[1] -> Draw();
  can->cd(6); hN[1] -> Draw();

  can->cd(3); hW[2] -> Draw();
  can->cd(7); hN[2] -> Draw();

  can->cd(4); hW[3] -> Draw();
  can->cd(8); hN[3] -> Draw();




}
