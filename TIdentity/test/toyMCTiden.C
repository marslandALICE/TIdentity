#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
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


/*

cd /u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/
aliroot -l
.L toyMCTiden.C+
toyMCTiden()

TFile f("dataTree.root")
TTree *tree = (TTree*)f.Get("tracks")
TFile g("LineShapes.root")
TClonesArray *cloneArrFunc   = (TClonesArray*)g.Get("funcLineShapesCArr");
TF1 *fLineShape[4]
for (Int_t ipart = 0; ipart<4; ipart++) {
TString objName = Form("particle_%d",ipart);
fLineShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
}
tree->Draw("dEdx>>h(900,0,30)")
fLineShape[0]->Draw("same")
fLineShape[1]->Draw("same")
fLineShape[2]->Draw("same")
fLineShape[3]->Draw("same")

*/
Int_t switchOffParticle=3;   // default = -100 for 4 particles
Int_t nParticles=4;
Int_t elMean = 6, piMean=10, kaMean=6, prMean=8;
Int_t nTracksPerEventArr[4]={0};
Double_t elParams[]={8.,1.5};
Double_t piParams[]={4.,1.};
Double_t kaParams[]={12.,1.5};
Double_t prParams[]={20.,1.6};
Float_t trackdEdx[1000]={0.};
TH1D *hParticles[4];
TH1D *hFirstMoms[4];
TF1 *fParticles[4];
UInt_t cutBit=0;
Int_t sign=0;
ULong64_t nEvents=500000;
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray funcLineShapesCArr("TF1",50000);
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kElEl=4,kPiPi=5,kKaKa=6,kPrPr=7,kElPi=8,kElKa=9,kElPr=10,kPiKa=11,kPiPr=12,kKaPr=13,};
TString momNames[14] = {"El1","Pi1","Ka1","Pr1","El2","Pi2","Ka2","Pr2","ElPi","ElKa","ElPr","PiKa","PiPr","KaPr"};

//
//
// ================================================================================================
//
//
void toyMCTiden(){

  TTreeSRedirector *outputData = new TTreeSRedirector("dataTree.root", "recreate");
  TTreeSRedirector *outputFits = new TTreeSRedirector("LineShapes.root", "recreate");
  funcLineShapesCArr.SetOwner(kTRUE);
  for (Int_t i=0;i<nParticles;i++) {
    hParticles[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),900,0.,30.);
    hFirstMoms[i] = new TH1D(Form("firstMom_%d",i),Form("firstMom_%d",i),900,0.,30.);
    fParticles[i] = new TF1(Form("particle_%d",i),"gaus",0,30);
  }

  Double_t pidVar = 0;
  TRandom randomGen;
  const Int_t nMoments = 14;
  TVectorF recMoments(nMoments);
  for(Int_t i=0;i<nMoments; i++) recMoments[i]=0.;

  for (ULong64_t ievent=0;ievent<nEvents;ievent++){

    Float_t cent=randomGen.Uniform(0,10);
    for (Int_t i=0;i<1000;i++) trackdEdx[i]=0.;
    Int_t trCount=0;
    nTracksPerEventArr[kEl] = randomGen.Poisson(elMean);  hFirstMoms[0]->Fill(nTracksPerEventArr[kEl]);
    nTracksPerEventArr[kPi] = randomGen.Poisson(piMean);  hFirstMoms[1]->Fill(nTracksPerEventArr[kPi]);
    nTracksPerEventArr[kKa] = randomGen.Poisson(kaMean);  hFirstMoms[2]->Fill(nTracksPerEventArr[kKa]);
    nTracksPerEventArr[kPr] = randomGen.Poisson(prMean);  hFirstMoms[3]->Fill(nTracksPerEventArr[kPr]);
    if (switchOffParticle>-1) nTracksPerEventArr[switchOffParticle]=0;
    //
    // Generate electron dEdx
    for (Int_t i=0; i<nTracksPerEventArr[kEl];i++){
      trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]);
      hParticles[0]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    //
    // Generate pion dEdx
    for (Int_t i=0; i<nTracksPerEventArr[kPi];i++){
      trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]);
      hParticles[1]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    //
    // Generate kaon dEdx
    for (Int_t i=0; i<nTracksPerEventArr[kKa];i++){
      trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]);
      hParticles[2]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    //
    // Generate proton dEdx
    for (Int_t i=0; i<nTracksPerEventArr[kPr];i++){
      trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]);
      hParticles[3]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    if(ievent%100000==0) cout << ievent << "  " << trCount << " " <<  nTracksPerEventArr[kEl] << "  " << nTracksPerEventArr[kPi] << "  " << nTracksPerEventArr[kKa] << "  " << nTracksPerEventArr[kPr] << endl;
    //
    // Fill data tree
    for (Int_t itr=0;itr<1000;itr++){
      outputData->GetFile()->cd();
      if(trackdEdx[itr]>0){
        *outputData << "tracks"
        << "gid=" << ievent <<
        "dEdx=" << trackdEdx[itr] <<
        "cutBit=" << cutBit <<
        "sign=" << sign <<
        "cent=" << cent <<
        "\n";
      }
    }
    //
    // Calculate Moments
    recMoments[kEl]=Float_t(nTracksPerEventArr[kEl]);
    recMoments[kPi]=Float_t(nTracksPerEventArr[kPi]);
    recMoments[kKa]=Float_t(nTracksPerEventArr[kKa]);
    recMoments[kPr]=Float_t(nTracksPerEventArr[kPr]);
    recMoments[kElEl]=recMoments[kEl]*recMoments[kEl];
    recMoments[kPiPi]=recMoments[kPi]*recMoments[kPi];
    recMoments[kKaKa]=recMoments[kKa]*recMoments[kKa];
    recMoments[kPrPr]=recMoments[kPr]*recMoments[kPr];
    recMoments[kElPi]=recMoments[kEl]*recMoments[kPi];
    recMoments[kElKa]=recMoments[kEl]*recMoments[kKa];
    recMoments[kElPr]=recMoments[kEl]*recMoments[kPr];
    recMoments[kPiKa]=recMoments[kPi]*recMoments[kKa];
    recMoments[kPiPr]=recMoments[kPi]*recMoments[kPr];
    recMoments[kKaPr]=recMoments[kKa]*recMoments[kPr];
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
  // dump line shapes
  for (Int_t i=0;i<nParticles;i++) {
    hParticles[i] -> SetLineColor(colors[i+1]);
    hFirstMoms[i] -> SetLineColor(colors[i+1]);
    fParticles[i] -> SetLineColor(colors[i+1]);
    fParticles[i] -> SetLineWidth(2);
    fParticles[i]->SetNpx(1000);
    fParticles[i]->SetParameters(hParticles[i]->GetMean(),hParticles[i]->GetRMS());
    funcLineShapesCArr[i] = (TF1*)fParticles[i];
    if (hParticles[i]) hParticles[i] -> Fit(fParticles[i],"QN");
  }
  outputFits->GetFile()->cd();
  funcLineShapesCArr . Write("funcLineShapesCArr",TObject::kSingleKey);
  funcLineShapesCArr.Clear("C");
  for (Int_t i=0;i<nParticles;i++) hParticles[i] -> Write();
  for (Int_t i=0;i<nParticles;i++) hFirstMoms[i] -> Write();
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
  TH1D *h = new TH1D("hGen","hGen",14,0.,14.);
  TH1D *htemp[14];
  for (Int_t i=1;i<15;i++) {
    tree->Draw(Form("moment.fElements[%d]",i-1),"","goff");
    htemp[i-1] = (TH1D*)tree->GetHistogram()->Clone();
    htemp[i-1]->SetName(Form("htmp_%d",i-1));
    cout << momNames[i-1] << "  " << htemp[i-1]->GetMean() << endl;
    h->SetBinContent(i,htemp[i-1]->GetMean());
  }
  cout << " ================================================================================== "<<endl;
  cout << " ================================================================================== "<<endl;
  for (Int_t i=1;i<15;i++) h->GetXaxis()->SetBinLabel(i,momNames[i-1]);
  TFile *outFile = new TFile("toyMC_Moments_Gen.root","recreate");
  h->Write();
  outFile -> Close();
  delete outFile; //yeni eklave etdim.

}
