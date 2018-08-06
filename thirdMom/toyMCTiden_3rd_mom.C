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

cd /home/marsland/Desktop/TIdentityCMake_Dev/test/thirdMom
aliroot -l
.L toyMCTiden_3rd_mom.C+
toyMCTiden_3rd_mom()
PlotdEdx()


n=10.; cout << n*n*n+3*n*n+n << "  " << n*n+n << endl;
n=14.; cout << n*n*n+3*n*n+n << "  " << n*n+n << endl;
n=8.;  cout << n*n*n+3*n*n+n << "  " << n*n+n << endl;

1310  110  -->  1308.51  109.914
3346  210  -->  3344.59  209.864
712  72    -->  711.713  71.986

*/

Int_t fUsedSign = 0;
ULong64_t nEvents=300000;
const Int_t nParticles=4;
Int_t elMean  = 0, piMean =14, kaMean =8, prMean =10;

Int_t nTracksPerEventArr[nParticles]={0};
Double_t elParams[]={8.,1.5};
Double_t piParams[]={4.,1.};
Double_t kaParams[]={11.,1.5};
Double_t prParams[]={17.,1.6};

// Double_t elParams[]={4.,1.5};
// Double_t piParams[]={15.,1.};
// Double_t kaParams[]={19.,1.5};
// Double_t prParams[]={26.,1.6};

const Int_t nMaxTracksPerEvent = 10000;
Float_t trackdEdx[nMaxTracksPerEvent]={0.};
Int_t trackSign[nMaxTracksPerEvent]={0};

TH1D *hParticles[nParticles];
TH1D *hFirstMoms[nParticles];
TF1 *fParticles[nParticles];
UInt_t cutBit=0;
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray funcLineShapesCArr("TF1",50000);
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,};
//
//
// ================================================================================================
//
//
void toyMCTiden_3rd_mom()
{

  TTreeSRedirector *outputData = new TTreeSRedirector("DataTree_Net.root", "recreate");
  TTreeSRedirector *outputFits = new TTreeSRedirector("LineShapes_Net.root", "recreate");
  funcLineShapesCArr.SetOwner(kTRUE);
  for (Int_t i=0;i<nParticles;i++) {
    hParticles[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),1500,0.,50.);
    hFirstMoms[i] = new TH1D(Form("firstMom_%d",i),Form("hist_%d",i),1500,0.,50.);
    fParticles[i] = new TF1(Form("particle_%d",i),"gaus",0,50);
  }

  Double_t pidVar = 0;
  TRandom randomGen;
  //
  for (ULong64_t ievent=0;ievent<nEvents;ievent++){  // event loop

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
    //
    // put extra correlation
    // nTracksPerEventArr[kEl]=nTracksPerEventArr[kKa];

    if(ievent%10000==0) {
      cout << " event = " << ievent << "  ----  " << trCount;
      cout <<  "  " << nTracksPerEventArr[kEl] ;
      cout <<  "  " << nTracksPerEventArr[kPi] ;
      cout <<  "  " << nTracksPerEventArr[kKa] ;
      cout <<  "  " << nTracksPerEventArr[kPr] ;
      cout << endl;
    }

    //
    // Generate electron dEdx
    for (Int_t ipart=0; ipart<nParticles; ipart++){

      for (Int_t i=0; i<nTracksPerEventArr[ipart];i++){

        if (ipart == kEl) {trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kPi) {trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kKa) {trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kPr) {trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        // hParticles[ipart]->Fill(trackdEdx[trCount]);
        trackSign[trCount]= 1;
        trCount++;
      }

    }
    //
    // Fill histograms first moms
    for (Int_t i=0;i<nParticles;i++) hFirstMoms[i]->Fill(nTracksPerEventArr[i]);
    //
    // Fill data tree
    for (Int_t itr=0;itr<10000;itr++){
      outputData->GetFile()->cd();
      if(trackdEdx[itr]!=0){
        *outputData << "tracks"
        << "gid="         << ievent <<
        "dEdx="           << trackdEdx[itr] <<
        "cutBit="         << cutBit <<
        "sign="           << trackSign[itr] <<
        "cent="           << cent <<
        "\n";
      }
    }

  } // event loop ends
  //
  // dump line shape
  Int_t objcounter=0;
  for (Int_t i=0;i<nParticles;i++) {
    hParticles[i] -> SetLineColor(colors[i+1]);
    hFirstMoms[i] -> SetLineColor(colors[i+1]);
    fParticles[i] -> SetLineColor(colors[i+1]);
    hParticles[i] -> SetDirectory(0);
    hFirstMoms[i] -> SetDirectory(0);
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

  tree->Draw("dEdx>>h(1500,0,50)");
  for (Int_t ipart = 0; ipart<nParticles; ipart++){
    fShape[ipart]->Draw("same");
  }

}
