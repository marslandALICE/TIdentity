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

pr=10.;
pi=14.;
k=8.+7;
{
cout << " pr3 " << pr*pr*pr+3*pr*pr+pr << "  " << pr*pr+pr << endl;
cout << " pi3 " << pi*pi*pi+3*pi*pi+pi << "  " << pi*pi+pi << endl;
cout << " K3  " << k*k*k+3*k*k+k << "  " << k*k+k << endl;
cout << " Pr2Pi  " << (pr*pr + pr)*pi << endl;
cout << " Pr2K   " << (pr*pr + pr)*k << endl;
cout << " Pi2K   " << (pi*pi + pi)*k << endl;
cout << " Pi2Pr  " << (pi*pi + pi)*pr << endl;
cout << " K2Pr   " << (k*k + k)*pr << endl;
cout << " K2Pi   " << (k*k + k)*pi << endl;
cout << " PiPrK  " << pi*pr*k << endl;
}

pr3 1310  110
pi3 3346  210
K3  712  72
Pr2Pi  1540
Pr2K   880
Pi2K   1680
Pi2Pr  2100
K2Pr   720
K2Pi   1008
PiPrK  1120


 pr3 1310  110
 pi3 3346  210
 K3  4065  240
 Pr2Pi  1540
 Pr2K   1650
 Pi2K   3150
 Pi2Pr  2100
 K2Pr   2400
 K2Pi   3360
 PiPrK  2100

*/

Int_t fUsedSign = 0;
ULong64_t nEvents=500000;
const Int_t nParticles=4;
Int_t piMean =14, kaMean =8, prMean =10, elMean  = 0;

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
TF1 *fSumElKa;
UInt_t cutBit=0;
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray funcLineShapesCArr("TF1",50000);
enum momentType{kPr=0,kPi=1,kKa=2,kEl=3,};
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
    nTracksPerEventArr[kPr] = randomGen.Poisson(prMean);
    nTracksPerEventArr[kPi] = randomGen.Poisson(piMean);
    nTracksPerEventArr[kKa] = randomGen.Poisson(kaMean);
    nTracksPerEventArr[kEl] = randomGen.Poisson(elMean);

    //
    // put extra correlation
    // nTracksPerEventArr[kEl]=nTracksPerEventArr[kKa];

    if(ievent%10000==0) {
      cout << " event = " << ievent << "  ----  " << trCount;
      cout <<  "  " << nTracksPerEventArr[kPr] ;
      cout <<  "  " << nTracksPerEventArr[kPi] ;
      cout <<  "  " << nTracksPerEventArr[kKa] ;
      cout <<  "  " << nTracksPerEventArr[kEl] ;
      cout << endl;
    }

    //
    // Generate electron dEdx
    for (Int_t ipart=0; ipart<nParticles; ipart++){

      for (Int_t i=0; i<nTracksPerEventArr[ipart];i++){

        if (ipart == kPr) {trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kPi) {trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kKa) {trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kEl) {trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]); hParticles[ipart]->Fill(trackdEdx[trCount]);}
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

  }

  Int_t objcounter=0;
  for (Int_t i=0;i<nParticles;i++) {

    if (i==2 && elMean>0) {
      fSumElKa = new TF1("particle_2","particle_2+particle_3",0,50);
      funcLineShapesCArr[objcounter] = (TF1*)fSumElKa;
    } else {
      funcLineShapesCArr[objcounter] = (TF1*)fParticles[i];
    }
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
