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
for (Int_t ipart = 0; ipart<5; ipart++) {
TString objName = Form("particle_%d",ipart);
fLineShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
}
tree->Draw("dEdx>>h(1500,0,50)")
fLineShape[0]->Draw("same")
fLineShape[1]->Draw("same")
fLineShape[2]->Draw("same")
fLineShape[3]->Draw("same")
fLineShape[4]->Draw("same")

*/
Int_t switchOffParticle=-100;   // default = -100 for 4 particles
const Int_t nParticles=5;
const Int_t nMoments = 19;
Int_t elMean = 6, piMean=10, kaMean=6, prMean=8, deMean=4;
Int_t nTracksPerEventArr[4]={0};
Double_t elParams[]={8.,1.5};
Double_t piParams[]={4.,1.};
Double_t kaParams[]={11.,1.5};
Double_t prParams[]={17.,1.6};
Double_t deParams[]={25.,1.8};

Float_t trackdEdx[1000]={0.};
TH1D *hParticles[nParticles];
TH1D *hFirstMoms[nParticles];
TF1 *fParticles[nParticles];
UInt_t cutBit=0;
Int_t sign=0;
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
    hParticles[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),1500,0.,50.);
    hFirstMoms[i] = new TH1D(Form("firstMom_%d",i),Form("firstMom_%d",i),1500,0.,50.);
    fParticles[i] = new TF1(Form("particle_%d",i),"gaus",0,50);
  }

  Double_t pidVar = 0;
  TRandom randomGen;
  TVectorF recMoments(nMoments);
  for(Int_t i=0;i<nMoments; i++) recMoments[i]=0.;

  for (ULong64_t ievent=0;ievent<nEvents;ievent++){

    Float_t cent=randomGen.Uniform(0,10);
    for (Int_t i=0;i<1000;i++) trackdEdx[i]=0.;
    Int_t trCount=0;
    nTracksPerEventArr[kEl] = randomGen.Poisson(elMean);
    nTracksPerEventArr[kPi] = randomGen.Poisson(piMean);
    nTracksPerEventArr[kKa] = randomGen.Poisson(kaMean);
    nTracksPerEventArr[kPr] = randomGen.Poisson(prMean);
    nTracksPerEventArr[kDe] = randomGen.Poisson(deMean);

    for (Int_t i=0;i<nParticles;i++){
      if ( !(switchOffParticle>-1 && i>=switchOffParticle) ) hFirstMoms[i]->Fill(nTracksPerEventArr[i]);
    }
    //
    // Switch of some particles
    for (Int_t i=0;i<nParticles;i++) {
      if ( switchOffParticle>-1 && i>=switchOffParticle ) nTracksPerEventArr[i]=0;
    }
    //
    // Generate electron dEdx
    for (Int_t i=0; i<nTracksPerEventArr[kEl];i++){
      trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
      hParticles[0]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    //
    // Generate pion dEdx
    for (Int_t i=0; i<nTracksPerEventArr[kPi];i++){
      trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
      hParticles[1]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    //
    // Generate kaon dEdx
    for (Int_t i=0; i<nTracksPerEventArr[kKa];i++){
      trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]);
      hParticles[2]->Fill(trackdEdx[trCount]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
      trCount++;
    }
    //
    // Generate proton dEdx
    for (Int_t i=0; i<nTracksPerEventArr[kPr];i++){
      trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]);
      hParticles[3]->Fill(trackdEdx[trCount]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
      trCount++;
    }
    // Generate proton dEdx
    for (Int_t i=0; i<nTracksPerEventArr[kDe];i++){
      trackdEdx[trCount] = randomGen.Gaus(deParams[0],deParams[1]);
      hParticles[4]->Fill(trackdEdx[trCount]);
      sign = (i%2==0) ? sign = 1 : sign = -1;
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
    recMoments[kDe]=Float_t(nTracksPerEventArr[kDe]);

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
