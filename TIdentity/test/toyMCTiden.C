#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TClonesArray.h"
#include "TTreeStream.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::cout;

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
tree->Draw("dEdx>>h(600,0,30)")
fLineShape[0]->Draw("same")
fLineShape[1]->Draw("same")
fLineShape[2]->Draw("same")
fLineShape[3]->Draw("same")

*/

Int_t elMean = 6, piMean=10, kaMean=6, prMean=8;
Double_t elParams[]={8.,1.5};
Double_t piParams[]={4.,1.};
Double_t kaParams[]={12.,1.5};
Double_t prParams[]={20.,1.6};
Float_t trackdEdx[1000]={0.};
TH1D *hParticles[4];
TF1 *fParticles[4];
UInt_t cutBit=0;
Int_t sign=0;
ULong64_t nEvents=1000000;
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray funcLineShapesCArr("TF1",50000);
//
//
// ================================================================================================
//
//
void toyMCTiden(){

  TTreeSRedirector *outputData = new TTreeSRedirector("dataTree.root", "recreate");
  TTreeSRedirector *outputFits = new TTreeSRedirector("LineShapes.root", "recreate");
  funcLineShapesCArr.SetOwner(kTRUE);
  for (Int_t i=0;i<4;i++) {
    hParticles[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),600,0.,30.);
    fParticles[i] = new TF1(Form("particle_%d",i),"gaus",0,30);
  }

  Double_t pidVar = 0;
  TRandom randomGen;
  for (ULong64_t ievent=0;ievent<nEvents;ievent++){

    for (Int_t i=0;i<1000;i++) trackdEdx[i]=0.;
    Int_t trCount=0;
    Int_t nEl = randomGen.Poisson(elMean);
    Int_t nPi = randomGen.Poisson(piMean);
    Int_t nKa = randomGen.Poisson(kaMean);
    Int_t nPr = randomGen.Poisson(prMean);
    //
    // Generate electron dEdx
    for (Int_t i=0; i<nEl;i++){
      trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]);
      hParticles[0]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    //
    // Generate pion dEdx
    for (Int_t i=0; i<nPi;i++){
      trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]);
      hParticles[1]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    //
    // Generate kaon dEdx
    for (Int_t i=0; i<nKa;i++){
      trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]);
      hParticles[2]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    //
    // Generate proton dEdx
    for (Int_t i=0; i<nPr;i++){
      trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]);
      hParticles[3]->Fill(trackdEdx[trCount]);
      trCount++;
    }
    if(ievent%100000==0) cout << ievent << "  " << trCount << " " <<  nEl << "  " << nPi << "  " << nKa << "  " << nPr << endl;

    //
    // Fill data tree
    for (Int_t itr=0;itr<nEl+nPi+nKa+nPr;itr++){
      outputData->GetFile()->cd();
      if(trackdEdx[itr]>0){
        *outputData << "tracks"
        << "gid=" << ievent <<
        "dEdx=" << trackdEdx[itr] <<
        "cutBit=" << cutBit <<
        "sign=" << sign <<
        "\n";
      }
    }

  } // event loop ends

  //
  // dump line shapes
  for (Int_t i=0;i<4;i++) {
    hParticles[i] -> SetLineColor(colors[i+1]);
    fParticles[i] -> SetLineColor(colors[i+1]);
    fParticles[i] -> SetLineWidth(2);
    fParticles[i]->SetNpx(1000);
    funcLineShapesCArr[i] = (TF1*)fParticles[i];
    hParticles[i] -> Fit(fParticles[i],"QN");
  }
  outputFits->GetFile()->cd();
  funcLineShapesCArr . Write("funcLineShapesCArr",TObject::kSingleKey);
  funcLineShapesCArr.Clear("C");
  for (Int_t i=0;i<4;i++) hParticles[i] -> Write();
  // for (Int_t i=0;i<4;i++) fParticles[i] -> Write();
  delete outputFits;
  delete outputData;

}
