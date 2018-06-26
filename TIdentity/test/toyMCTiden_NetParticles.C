#include "TFile.h"
#include "TH1.h"
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
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using std::cout;
using std::setw;


void PlotdEdx();

/*

cd /u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/
aliroot -l
.L toyMCTiden_NetParticles.C+
toyMCTiden_NetParticles()
PlotdEdx()

*/
Int_t nSubSample = 25;
Int_t switchOffParticle=-100;   // default = -100 for 4 particles
Int_t fUsedSign = 0;
ULong64_t nEvents=100000;
const Int_t nParticles=10;
const Int_t nMoments = 19;
Int_t elMean  = 8, piMean =14, kaMean =8, prMean =10, deMean =6;
Int_t elBMean = 6, piBMean=10, kaBMean=6, prBMean=8, deBMean=5;

Int_t nTracksPerEventArr[nParticles]={0};
Double_t elParams[]={8.,1.5};
Double_t piParams[]={4.,1.};
Double_t kaParams[]={11.,1.5};
Double_t prParams[]={17.,1.6};
Double_t deParams[]={25.,1.8};
const Int_t nMaxTracksPerEvent = 10000;
Float_t trackdEdx[nMaxTracksPerEvent]={0.};
Int_t trackSign[nMaxTracksPerEvent]={0};

TH1D *hParticles[nParticles];
TH1D *hFirstMoms[nParticles];
TF1 *fParticles[nParticles];
UInt_t cutBit=0;
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray funcLineShapesCArr("TF1",50000);
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kDe=4,
  kBEl=5,kBPi=6,kBKa=7,kBPr=8,kBDe=9,};
//
//
// ================================================================================================
//
//
void toyMCTiden_NetParticles()
{

  TTreeSRedirector *outputData = new TTreeSRedirector("DataTree_Net.root", "recreate");
  TTreeSRedirector *outputFits = new TTreeSRedirector("LineShapes_Net.root", "recreate");
  funcLineShapesCArr.SetOwner(kTRUE);
  for (Int_t i=0;i<nParticles;i++) {
    hParticles[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),3000,-50.,50.);
    hFirstMoms[i] = new TH1D(Form("firstMom_%d",i),Form("hist_%d",i),3000,-50.,50.);
    fParticles[i] = new TF1(Form("particle_%d",i),"gaus",-50,50);
  }

  Double_t pidVar = 0;
  TRandom randomGen;
  TVectorF recMoments1st(nParticles);
  TVectorF recMoments2nd(nParticles*nParticles);
  for(Int_t i=0;i<nParticles; i++) recMoments1st[i]=0.;
  for(Int_t i=0;i<nParticles*nParticles; i++) recMoments2nd[i]=0.;

  Float_t subsampleID = -1;
  for (ULong64_t ievent=0;ievent<nEvents;ievent++){  // event loop
    subsampleID = (ievent+1)%nSubSample;
    // rest particle counters
    for (Int_t j=0;j<nParticles;j++){
      nTracksPerEventArr[j]=0;
    }
    // reset track counter
    for (Int_t i=0;i<nMaxTracksPerEvent;i++) { trackdEdx[i]=0.; trackSign[i]=0; }
    //
    Float_t cent=randomGen.Uniform(0,10);
    Int_t trCount=0;
    nTracksPerEventArr[kEl] = randomGen.Poisson(elMean);
    nTracksPerEventArr[kPi] = randomGen.Poisson(piMean);
    nTracksPerEventArr[kKa] = randomGen.Poisson(kaMean);
    nTracksPerEventArr[kPr] = randomGen.Poisson(prMean);
    nTracksPerEventArr[kDe] = randomGen.Poisson(deMean);

    nTracksPerEventArr[kBEl] = randomGen.Poisson(elBMean);
    nTracksPerEventArr[kBPi] = randomGen.Poisson(piBMean);
    nTracksPerEventArr[kBKa] = randomGen.Poisson(kaBMean);
    nTracksPerEventArr[kBPr] = randomGen.Poisson(prBMean);
    nTracksPerEventArr[kBDe] = randomGen.Poisson(deBMean);
    //
    // put extra correlation
    nTracksPerEventArr[kEl]=nTracksPerEventArr[kKa];

    if(ievent%100000==0) {
      cout << ievent << "  " << trCount;
      cout <<  "  " << nTracksPerEventArr[kEl] <<  "  " << nTracksPerEventArr[kBEl] ;
      cout <<  "  " << nTracksPerEventArr[kPi] <<  "  " << nTracksPerEventArr[kBPi] ;
      cout <<  "  " << nTracksPerEventArr[kKa] <<  "  " << nTracksPerEventArr[kBKa] ;
      cout <<  "  " << nTracksPerEventArr[kPr] <<  "  " << nTracksPerEventArr[kBPr] ;
      cout <<  "  " << nTracksPerEventArr[kDe] <<  "  " << nTracksPerEventArr[kBDe] ;
      cout << endl;
    }

    //
    // Switch of some particles
    for (Int_t i=0;i<nParticles;i++) {
      if ( switchOffParticle>-1 && i>=switchOffParticle ) nTracksPerEventArr[i]=0;
    }
    //
    // Generate electron dEdx
    for (Int_t ipart=0; ipart<nParticles; ipart++){

      for (Int_t i=0; i<nTracksPerEventArr[ipart];i++){

        if (ipart == kEl) {trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kPi) {trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kKa) {trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kPr) {trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kDe) {trackdEdx[trCount] = randomGen.Gaus(deParams[0],deParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}

        if (ipart == kBEl) {trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1])*-1.; hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kBPi) {trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1])*-1.; hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kBKa) {trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1])*-1.; hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kBPr) {trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1])*-1.; hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kBDe) {trackdEdx[trCount] = randomGen.Gaus(deParams[0],deParams[1])*-1.; hParticles[ipart]->Fill(trackdEdx[trCount]);}

        // hParticles[ipart]->Fill(trackdEdx[trCount]);
        trackSign[trCount]= (ipart<nParticles/2.) ? 1 : -1 ;
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

    if(ievent%100000==0) {
      cout << ievent << "  " << trCount;
      cout <<  "  " << nTracksPerEventArr[kEl] <<  "  " << nTracksPerEventArr[kBEl] ;
      cout <<  "  " << nTracksPerEventArr[kPi] <<  "  " << nTracksPerEventArr[kBPi] ;
      cout <<  "  " << nTracksPerEventArr[kKa] <<  "  " << nTracksPerEventArr[kBKa] ;
      cout <<  "  " << nTracksPerEventArr[kPr] <<  "  " << nTracksPerEventArr[kBPr] ;
      cout <<  "  " << nTracksPerEventArr[kDe] <<  "  " << nTracksPerEventArr[kBDe] ;
      cout << endl;
    }
    //
    // Fill data tree
    for (Int_t itr=0;itr<10000;itr++){
      outputData->GetFile()->cd();
      if(trackdEdx[itr]!=0){
        *outputData << "tracks"
        << "gid="         << ievent <<
        "subsampleindex=" << subsampleID <<
        "dEdx="           << trackdEdx[itr] <<
        "cutBit="         << cutBit <<
        "sign="           << trackSign[itr] <<
        "cent="           << cent <<
        "\n";
      }
    }
    //
    // Calculate 1st Moments
    for (Int_t i=0; i<nParticles;i++) recMoments1st[i]=Float_t(nTracksPerEventArr[i]);
    //
    // Calculate 2nd Moments
    Int_t k=0;
    for (Int_t i=0; i<nParticles;i++) {
      for (Int_t j=0; j<nParticles;j++) {
        recMoments2nd[k]=recMoments1st[i]*recMoments1st[j];
        k++;
      }
    }
    //
    // Dump moments to tree
    outputData->GetFile()->cd();
    *outputData << "events"
    << "gid=" << ievent <<
    "moment1.=" << &recMoments1st <<             // second moments for particle+antiparticle
    "moment2.=" << &recMoments2nd <<             // second moments for particle+antiparticle
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

  for (Int_t i=0;i<nParticles;i++) {
    fParticles[i]->SetLineWidth(2);
    fParticles[i]->SetNpx(1000);
    fParticles[i]->SetParameters(hParticles[i]->GetMean(),hParticles[i]->GetRMS());
    if (hParticles[i]) hParticles[i] -> Fit(fParticles[i],"QN");
    funcLineShapesCArr[objcounter] = (TF1*)fParticles[i];
    objcounter++;
  }

  outputFits->GetFile()->cd();
  funcLineShapesCArr.Write("funcLineShapesCArr",TObject::kSingleKey);
  funcLineShapesCArr.Clear("C");

  for (Int_t i=0;i<nParticles;i++) { hParticles[i] -> Write(); hFirstMoms[i] -> Write();}

  delete outputFits;
  delete outputData;

}

void PlotdEdx()
{


  TFile *f = new TFile("DataTree_Net.root");
  TTree *tree = (TTree*)f->Get("tracks");
  TFile *g = new TFile("LineShapes_Net.root");
  TClonesArray *cloneArrFunc = (TClonesArray*)g->Get("funcLineShapesCArr");
  TF1 *fShape[nParticles];

  for (Int_t ipart = 0; ipart<nParticles; ipart++) {
    TString objName = Form("particle_%d",ipart);
    fShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
  }

  tree->Draw("dEdx>>h(3000,-50,50)");
  for (Int_t ipart = 0; ipart<nParticles; ipart++){
    fShape[ipart]->Draw("same");
  }

}
