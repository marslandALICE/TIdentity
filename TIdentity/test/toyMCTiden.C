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


void PrintMoments();
void PlotdEdx(Int_t sign);
void PlotWs();
void PlotRecGenRatio();
void PlotSubSample();

/*

cd /u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/
aliroot -l
.L toyMCTiden.C+
toyMCTiden()
PlotdEdx(0)

PlotWs()

*/
Int_t nSubSample = 25;
Int_t switchOffParticle=-100;   // default = -100 for 4 particles
Int_t fUsedSign = 0;
ULong64_t nEvents=10000000;
const Int_t nParticles=5;
const Int_t nSignBins =3;
const Int_t nMoments = 19;
Int_t elMean = 8, piMean=14, kaMean=8, prMean=10, deMean=6;
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
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray funcLineShapesCArr("TF1",50000);
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kDe=4,
  kElEl=5,kPiPi=6,kKaKa=7,kPrPr=8,kDeDe=9,
  kElPi=10,kElKa=11,kElPr=12,kElDe=13,
  kPiKa=14,kPiPr=15,kPiDe=16,
  kKaPr=17,kKaDe=18,};
TString fMomNames[nMoments] = {"#LT e #GT","#LT #pi #GT","#LT K #GT","#LT p #GT","#LT d #GT",
  "#LT e^{2} #GT","#LT #pi^{2} #GT","#LT K^{2} #GT","#LT p^{2} #GT","#LT d^{2} #GT",
  "#LT e,#pi #GT","#LT e,K #GT","#LT e,p #GT","#LT e,d #GT",
  "#LT #pi,K #GT","#LT #pi,p #GT","#LT #pi,d #GT",
  "#LT K,p #GT","#LT K,d #GT"};
//
//
// ================================================================================================
//
//
void toyMCTiden()
{

  TTreeSRedirector *outputData = new TTreeSRedirector("DataTree.root", "recreate");
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

  Float_t subsampleID = -1;
  for (ULong64_t ievent=0;ievent<nEvents;ievent++){  // event loop
    subsampleID = (ievent+1)%nSubSample;
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
    // put extra correlation
    nTracksPerEventArr[2][kEl]=nTracksPerEventArr[2][kKa];
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
      sign = (i%3==0) ? sign = 1 : sign = -1;
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
      sign = (i%3==0) ? sign = 1 : sign = -1;
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

void PrintMoments()
{

  cout << " ================================================================================== "<<endl;
  cout << " ================================================================================== "<<endl;
  TFile *f = new TFile("DataTree.root");
  TTree *tree = (TTree*)f->Get("events");
  TH1D *h = new TH1D("hGen","hGen",nMoments,0.,nMoments);
  TH1D *htemp[nMoments];
  for (Int_t i=1;i<nMoments+1;i++) {
    tree->Draw(Form("moment.fElements[%d]",i-1),"","goff");
    htemp[i-1] = (TH1D*)tree->GetHistogram()->Clone();
    htemp[i-1]->SetName(Form("htmp_%d",i-1));
    cout << fMomNames[i-1] << "  " << htemp[i-1]->GetMean() << endl;
    h->SetBinContent(i,htemp[i-1]->GetMean());
  }
  cout << " ================================================================================== "<<endl;
  cout << " ================================================================================== "<<endl;
  for (Int_t i=1;i<nMoments+1;i++) h->GetXaxis()->SetBinLabel(i,fMomNames[i-1]);
  TFile *outFile = new TFile("ToyMC_Moments.root","recreate");
  h->Write();
  outFile -> Close();
  delete outFile; //yeni eklave etdim.

}

void PlotdEdx(Int_t sign)
{


  TFile *f = new TFile("DataTree.root");
  TTree *tree = (TTree*)f->Get("tracks");
  TFile *g = new TFile("LineShapes.root");
  TClonesArray *cloneArrFunc = (TClonesArray*)g->Get("funcLineShapesCArr");
  TF1 *fShape[nParticles];

  Int_t binIndex=0;
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

  /*

  cd /u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/
  aliroot -l
  .L toyMCTiden.C+
  PlotWs()

  */

  gStyle->SetTextFont(62);
  TGaxis::SetMaxDigits(3);
  const Int_t npars = 5;
  TFile *f = new TFile("DataTree.root");
  TTree *tree = (TTree*)f->Get("tracks");
  TFile *g = new TFile("LineShapes.root");
  TFile *h = new TFile("TIdenDebug.root");
  TClonesArray *cloneArrFunc = (TClonesArray*)g->Get("funcLineShapesCArr");
  TF1 *fShape[npars];
  TH1D *hW[npars];
  TH1D *hN[npars];
  TH1D *hO[npars];

  TString WNames[npars] = {"W_{e}","W_{#pi}","W_{K}","W_{p}","W_{d}"};
  TString momNames[npars] = {"N_{e}","N_{#pi}","N_{K}","N_{p}","N_{d}"};
  TString OmegaNames[npars] = {"#omega_{e}","#omega_{#pi}","#omega_{K}","#omega_{p}","#omega_{e}"};

  for (Int_t ipart = 0; ipart<npars; ipart++) {

    TString objName = Form("particle_%d_bin_0",ipart);
    fShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
    hN[ipart]     = (TH1D*)g->Get(Form("firstMom_%d_bin_2",ipart));
    hW[ipart]     = (TH1D*)h->Get(Form("hW_%d",ipart));
    hO[ipart]     = (TH1D*)h->Get(Form("hOmega_%d",ipart));

    hN[ipart]->GetXaxis()->SetRangeUser(0,30);
    hW[ipart]->GetXaxis()->SetRangeUser(0,30);

    hN[ipart]->GetYaxis()->SetTitle("Counts (arb. units)");
    hW[ipart]->GetYaxis()->SetTitle("Counts (arb. units)");
    hO[ipart]->GetYaxis()->SetTitle("Counts (arb. units)");

    hN[ipart]->GetXaxis()->SetTitle(momNames[ipart]);
    hW[ipart]->GetXaxis()->SetTitle(WNames[ipart]);
    hO[ipart]->GetXaxis()->SetTitle(OmegaNames[ipart]);


    hN[ipart]->GetYaxis()->SetTitleOffset(1.2);
    hN[ipart]->GetXaxis()->SetTitleOffset(1.);
    hN[ipart]->GetYaxis()->SetTitleSize(0.06);
    hN[ipart]->GetXaxis()->SetTitleSize(0.07);
    hN[ipart]->GetXaxis()->SetLabelSize(0.06);
    hN[ipart]->GetYaxis()->SetLabelSize(0.06);
    hN[ipart]->GetYaxis()->SetLabelOffset(0.01);

    hW[ipart]->GetYaxis()->SetTitleOffset(1.2);
    hW[ipart]->GetXaxis()->SetTitleOffset(1.);
    hW[ipart]->GetYaxis()->SetTitleSize(0.06);
    hW[ipart]->GetXaxis()->SetTitleSize(0.07);
    hW[ipart]->GetXaxis()->SetLabelSize(0.06);
    hW[ipart]->GetYaxis()->SetLabelSize(0.06);
    hW[ipart]->GetYaxis()->SetLabelOffset(0.01);

    hO[ipart]->GetYaxis()->SetTitleOffset(1.15);
    hO[ipart]->GetXaxis()->SetTitleOffset(1.);
    hO[ipart]->GetYaxis()->SetTitleSize(0.06);
    hO[ipart]->GetYaxis()->SetLabelSize(0.06);
    hO[ipart]->GetXaxis()->SetTitleSize(0.07);
    hO[ipart]->GetXaxis()->SetLabelSize(0.06);
    hO[ipart]->GetYaxis()->SetLabelOffset(0.01);
    hO[ipart]->SetStats(kFALSE);

    hO[ipart]->GetXaxis()->SetNdivisions(505);
    hW[ipart]->GetXaxis()->SetNdivisions(505);
    hN[ipart]->GetXaxis()->SetNdivisions(505);

    hW[ipart] -> SetLineColor(colors[ipart+1]);
    hN[ipart] -> SetLineColor(colors[ipart+1]);
    hO[ipart] -> SetLineColor(colors[ipart+1]);

    // hN[ipart]->SetLineWidth(2);
    // hW[ipart]->SetLineWidth(2);
    // hO[ipart]->SetLineWidth(2);

    hW[ipart]->Scale(1./hW[ipart]->GetMaximum());
    hN[ipart]->Scale(1./hN[ipart]->GetMaximum());

    hN[ipart]->GetYaxis()->SetRangeUser(0.01,1.1);
    hW[ipart]->GetYaxis()->SetRangeUser(0.01,1.1);

  }

  TCanvas *can = new TCanvas("can", "can", 900, 1500);
  can->Divide(3,npars);
  gStyle->SetOptStat(1100);
  gStyle->SetStatY(0.93);
  gStyle->SetStatX(0.93);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.32);
  can->cd(1)->SetLogy();  gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hO[0] -> Draw();
  can->cd(2);  gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hW[0] -> Draw();
  can->cd(3);  gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hN[0] -> Draw();

  can->cd(4)->SetLogy();;  gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hO[1] -> Draw();
  can->cd(5);  gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hW[1] -> Draw();
  can->cd(6);  gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hN[1] -> Draw();

  can->cd(7)->SetLogy();;  gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.01); hO[2] -> Draw();
  can->cd(8);  gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hW[2] -> Draw();
  can->cd(9);  gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hN[2] -> Draw();

  can->cd(10)->SetLogy();;  gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.01); hO[3] -> Draw();
  can->cd(11); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hW[3] -> Draw();
  can->cd(12); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hN[3] -> Draw();

  can->cd(13)->SetLogy();; gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.01); hO[4] -> Draw();
  can->cd(14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hW[4] -> Draw();
  can->cd(15); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.01); hN[4] -> Draw();

  //
  //

  TH1D *h0pars[5],*h1pars[5];
  TFile *fGenvsRec = new TFile("TIdentity_Moments.root");
  TH1D *h0 = (TH1D*)fGenvsRec->Get("hDedxDebug_0");
  h0 -> GetXaxis()->SetTitle("d#it{E}/d#it{x} Signal (arb. units)");
  h0 -> GetYaxis()->SetTitle("Counts");
  h0 -> GetYaxis()->SetTitleOffset(1.3);
  h0 -> SetLineColor(kBlack);
  h0 -> SetMarkerStyle(24);
  h0 -> SetStats(kFALSE);
  h0 -> Rebin(6); h0->Scale(1./6.);
  h0 -> GetXaxis()->SetRangeUser(0,40);
  TH1D *h1 = (TH1D*)fGenvsRec->Get("hDedxDebug_1");
  h1 -> GetXaxis()->SetTitle("d#it{E}/d#it{x} Signal (arb. units)");
  h1 -> GetYaxis()->SetTitle("Counts");
  h1 -> GetYaxis()->SetTitleOffset(1.3);
  h1 -> SetLineColor(kBlack);
  h1 -> SetStats(kFALSE);
  h1 -> SetMarkerStyle(24);
  h1 -> Rebin(6); h1->Scale(1./6.);
  h1 -> GetXaxis()->SetRangeUser(0,40);
  for (Int_t i=0; i<npars;i++){
    h0pars[i] = (TH1D*)fGenvsRec->Get(Form("LineShape_sign_0_part_%d",i));
    h1pars[i] = (TH1D*)fGenvsRec->Get(Form("LineShape_sign_1_part_%d",i));
    h1pars[i]->SetLineWidth(1);
    h0pars[i]->SetLineWidth(1);
  }

  TLegend *leg = new TLegend(0.8, 0.6, 0.95, 0.9);
  leg->SetTextFont(62);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(4000);
  leg->AddEntry(h0,        "data" ,"LP");
  leg->AddEntry(h0pars[0], "e"    ,"L");
  leg->AddEntry(h0pars[1], "#pi"  ,"L");
  leg->AddEntry(h0pars[2], "K"    ,"L");
  leg->AddEntry(h0pars[3], "p"    ,"L");
  leg->AddEntry(h0pars[4], "d"    ,"L");

  TCanvas *canD = new TCanvas("canD", "canD", 1700, 700);
  canD->Divide(2,1);
  canD->cd(1);  gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.09); gPad->SetLeftMargin(0.09);
  h0->Draw("E"); for (Int_t i=0; i<npars;i++) h0pars[i]->Draw("same");
  leg->Draw("same");

  canD->cd(2);  gPad->SetRightMargin(0.02); gPad->SetTopMargin(0.04); gPad->SetBottomMargin(0.09); gPad->SetLeftMargin(0.09);
  h1->Draw("E"); for (Int_t i=0; i<npars;i++) h1pars[i]->Draw("same");
  leg->Draw("same");
  //
  //
  TH1D *hR = (TH1D*)fGenvsRec->Get("hRatio");
  for (Int_t i=1;i<nMoments+1;i++) hR->GetXaxis()->SetBinLabel(i,fMomNames[i-1]);
  hR->SetStats(kFALSE);
  TCanvas *canR = new TCanvas("canR", "canR", 1700, 700);
  canR->cd()->SetGridx();   gPad->SetRightMargin(0.01); gPad->SetTopMargin(0.02); gPad->SetBottomMargin(0.09); gPad->SetLeftMargin(0.12);
  hR->GetYaxis()->SetNdivisions(505);
  hR->GetYaxis()->SetTitle("Generated/Reconstructed");
  hR->GetYaxis()->SetTitleOffset(1.1);
  hR->SetMarkerColor(kBlack);
  hR->SetMarkerStyle(20);
  hR->Draw();

  can ->SaveAs("can_wDists.pdf");
  canD->SaveAs("can_dEdx.pdf");
  canR->SaveAs("can_ratio.pdf");


}

void PlotRecGenRatio(){

  gStyle->SetOptStat(0);
  TFile *outFile = new TFile("Ratio_GenVsRec.root","recreate");
  TFile *fGen = new TFile("ToyMC_Moments.root");
  TFile *fRec = new TFile("TIdentity_Moments.root");

  TTree *momTree = (TTree*)fRec->Get("momTree");
  TH1D *hGen   = (TH1D*)fGen->Get("hGen"); hGen->SetLineColor(kBlack);
  TH1D *hRec   = (TH1D*)hGen->Clone();   hRec->SetName("hRec"); hRec->SetLineColor(kRed+1);
  TH1D *hRatio = (TH1D*)hGen->Clone();   hRatio->SetName("hRatio"); hRatio->GetYaxis()->SetTitle("Generated/Reconstructed");

  outFile->cd();
  for (Int_t i=0;i<nMoments;i++) {

    momTree->Draw(Form("moment.fElements[%d]:subsampleindex",i),"","goff");
    TGraphErrors *gr = new TGraphErrors(momTree->GetSelectedRows(),momTree->GetV2(),momTree->GetV1());
    gr -> SetMarkerStyle(20);
    gr -> SetName(Form("mom_%d",i));
    gr->Write();
    // Calculate error for a given moment
    Double_t mean = gr->GetMean(2);
    Int_t N = gr->GetN();
    Double_t errsum = 0.;
    for (Int_t ibin=0;ibin<N;ibin++) errsum = (gr->GetY()[ibin]-mean)*(gr->GetY()[ibin]-mean);
    Double_t err = TMath::Sqrt(errsum/(N*(N-1)));
    hRec->SetBinContent(i+1,mean);
    hRec->SetBinError(i+1,err);
    cout << i << "  " << mean << "  " << err << endl;


    Double_t xx    = hGen->GetBinContent(i+1);
    Double_t errxx = 0.;

    Double_t yy    = mean;
    Double_t erryy = err;

    Double_t ratio = xx/yy;
    Double_t errRatio = TMath::Sqrt((errxx*errxx)/(yy*yy)+(xx*xx*erryy*erryy)/(yy*yy*yy*yy));

    hRatio->SetBinContent(i+1,ratio);
    hRatio->SetBinError(i+1,errRatio);
    hRatio->SetMarkerStyle(20);



  }
  hRatio->Write();
  hRatio->Draw();
  delete outFile;
  // Read generated histogram
  // for (Int_t i=1;i<nMoments+1;i++) {
  //   Double_t gen = hGen->GetBinContent(i);
  //   Double_t rec = hRec->GetBinContent(i);
  //   Double_t ratio = hGen->GetBinContent(i)/hRec->GetBinContent(i);
  //   cout << std::setw(10) << fMomNames[i-1] << " -->  gen:   ";
  //   cout << std::setw(10) << gen           << " -->  ";
  //   cout << std::setw(10) << rec           << "    ///  ";
  //   cout << std::setw(10) << ratio << endl;
  //   hRatio->SetBinContent(i,ratio);
  // }
}

void PlotSubSample()
{

  TFile *fTiden = new TFile("TIdentity_Moments.root");
  TFile *fToymc = new TFile("ToyMC_Moments.root");

  TTree *tree = (TTree*)fTiden->Get("momTree");
  TH1D *hist = (TH1D*)fToymc->Get("hGen");

  tree->Draw("moment.fElements[2]:subsampleindex+1","");
  TGraphErrors *gr = new TGraphErrors(tree->GetSelectedRows(),tree->GetV2(),tree->GetV1());
  Double_t mean = gr->GetMean(2);
  gr->SetMarkerStyle(20);
  gr->GetYaxis()->SetNdivisions(505);
  gr->GetYaxis()->SetTitle("#LT K #GT");
  gr->GetXaxis()->SetTitle("subsample index");
  gr->Draw("ap");
  //
  Double_t rec = hist->GetBinContent(3);
  TLine* lineAna = new TLine(0., mean,27.5, mean);
  lineAna -> SetLineColor(kBlack); lineAna -> SetLineWidth(2); lineAna -> SetLineStyle(1);
  lineAna->Draw("same");
  TLine* lineRec = new TLine(0., rec,27.5, rec);
  lineRec -> SetLineColor(kGreen+2); lineRec -> SetLineWidth(4); lineRec -> SetLineStyle(2);
  lineRec->Draw("same");

}
