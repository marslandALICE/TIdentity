#include <TFile.h>
#include <TH3.h>
#include "TSystem.h"
#include "TROOT.h"
#include "THn.h"
#include "TCut.h"
#include "TGaxis.h"
#include "TMinuit.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <THnSparse.h>
#include "TMath.h"
#include "TLine.h"
#include "TStopwatch.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTreeStream.h"
#include "TGraphSmooth.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TError.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TLinearFitter.h"
#include "TVirtualFitter.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TLegend.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::cout;
using std::setw;


// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************* Modification region Start ********************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************


// For free Kurtosis and skewness

// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************* Modification region Ends *********************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************


Int_t fSign=0;
Int_t fNSlices = 0;
Double_t parErrors[30]={};
// Define amp mean and sigma parlimits


// Defien max position in each expected particle mean bin
Double_t maxBin,maxBin0,maxBin1,maxBin2,maxBin3,maxBin4;

// Containers
TClonesArray fitResArr("TH1D", 1000);
TClonesArray fitResiduals("TH1D", 1000);
TClonesArray histLineShapesCArr("TH1D",50000);
TClonesArray funcLineShapesCArr("TF1",50000);

TObjArray    cleanResArr(1000);
TObjArray    cleanResArrFreeKS(1000);


// Debuggers and files
const Int_t   fnEtaBins       = 16;///16;    MC: 8, Data:16
const Float_t fEtaRangeDown   = -0.8;///16;
const Float_t fEtaRangeUp     = 0.8;///16;
const Int_t   fnMomBins       = 150;
const Float_t fMomRangeDown   = 0.2;
const Float_t fMomRangeUp     = 3.2;
const Int_t   fnCentBins      = 9;
Float_t     xCentBins[] = {0, 5,  10,  20, 30, 40, 50, 60, 70, 80};
Double_t skewMin       = 0.;
Double_t skewMax       = 1.5;
Double_t kurtosisMin   = 1.7;
Double_t kurtosisMax   = 2.2;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

enum ParticleType {
  kElectron = 0,
  kPion     = 1,
  kKaon     = 2,
  kProton   = 3,
  kDeuteron = 4,
} pType;

//TString fitFunctionGenGaus = "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
TString fitFunctionGenGaus = "[0]*exp(-(abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/1.414213))";

Float_t fEtaDown;
Float_t fCentDown;
Int_t corrType=-1;
Int_t setting = 0;

Int_t particleBin,centBin,etaBin,momBin,signBin;

TString filesCorr[]={
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_2000_8000_timeSeries_0.88_vZPileUp_0/HistsCorr.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_-100000_100000_timeSeries_0_vZPileUp_0/HistsCorr.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_-100000_1000_timeSeries_0.88_vZPileUp_0/HistsCorr.root", //
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_2000_100000_timeSeries_0.88_vZPileUp_-1/HistsCorr.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_2000_100000_timeSeries_0.88_vZPileUp_1/HistsCorr.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_8000_100000_timeSeries_0.88_vZPileUp_-1/HistsCorr.root",  //
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_8000_100000_timeSeries_0.88_vZPileUp_1/HistsCorr.root"    //
};

TString files[]={
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_2000_8000_timeSeries_0.88_vZPileUp_0/Hists.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_-100000_100000_timeSeries_0_vZPileUp_0/Hists.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_-100000_1000_timeSeries_0.88_vZPileUp_0/Hists.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_2000_100000_timeSeries_0.88_vZPileUp_-1/Hists.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_2000_100000_timeSeries_0.88_vZPileUp_1/Hists.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_8000_100000_timeSeries_0.88_vZPileUp_-1/Hists.root",
  "/home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits/mergedHists/Settings_pileUp_8000_100000_timeSeries_0.88_vZPileUp_1/Hists.root"
};

const Int_t nPars = 10;
const Int_t nCentDim = 10;
const Int_t nEtaDim = 17;
Float_t etaBinning[nEtaDim]   = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
Float_t centBinning[nCentDim] = { 0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
const Int_t colors[] = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TTreeSRedirector *histsStream=0;

// Helper functions for severel computations
void FitParticleSample(TH1D *hClean, Double_t *parSample, Double_t *parError);
void Fit_APWPerFile(TString corrStrr, TString fileData, const Int_t nSlice, Double_t ptMin, Double_t ptMax, Double_t etaMin, Double_t etaMax,Double_t centMin, Double_t centMax);

// Set functions for the parameters

void Fit_APW(TString corrStrIn, const Int_t nSlice, Double_t ptMin, Double_t ptMax)
{

  //
  // Apply gaussian fits to the ptot slices --> Main Function
  /*

  meld /u/marsland/PHD/macros/marsland_EbyeRatios/Fit_APW.C /home/marsland/Desktop/ubuntu_desktop/workdir/code/Fit_APW.C

  cd /home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits
  aliroot -l
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/code/Fit_APW.C+
  Fit_APW("Corr", 100,0.2,2.2); 2> err.log


  */
  // Double_t etaMin=-0.3, etaMax=-0.2, centMin=0, centMax=5;
  histsStream = new TTreeSRedirector(Form("Debug_%s.root",corrStrIn.Data()),"recreate");

  for (Int_t icent=0; icent<nCentDim-1; icent++){
    // for (Int_t ieta=0; ieta<1; ieta++){
      for (Int_t ieta=0; ieta<nEtaDim-1; ieta++){
      cout << "eta = " << etaBinning[ieta] <<" - " << etaBinning[ieta+1] << "     cent = " << centBinning[icent] << " - " << centBinning[icent+1] << endl;
      //
      for (Int_t i=0; i<7; i++){
        setting = i;
        if (corrStrIn=="Corr") { corrType = 1; Fit_APWPerFile("Corr", filesCorr[i], nSlice, ptMin, ptMax, etaBinning[ieta], etaBinning[ieta+1], centBinning[icent], centBinning[icent+1]);}
        if (corrStrIn=="")     { corrType = 0; Fit_APWPerFile("",     files[i],     nSlice, ptMin, ptMax, etaBinning[ieta], etaBinning[ieta+1], centBinning[icent], centBinning[icent+1]);}
      }
    }
  }

  delete histsStream;

}


// -------------------------------------------------------------------------------------------------------
void Fit_APWPerFile(TString corrStrr, TString fileData, const Int_t nSlice, Double_t ptMin, Double_t ptMax, Double_t etaMin, Double_t etaMax,Double_t centMin, Double_t centMax)
{

  //
  // Apply gaussian fits to the ptot slices --> Main Function
  /*

  cd /home/marsland/Desktop/ubuntu_desktop/workdir/TEST/filterTreesMakeHists/fits
  aliroot -l
  .L /home/marsland/Desktop/ubuntu_desktop/workdir/code/Fit_APW.C+
  Fit_APWPerFile("Corr","HistsCorr.root",100,0.2,2.2,  -0.3,-0.2,  0,5); 2> err.log


  */
  // Double_t etaMin=-0.3, etaMax=-0.2, centMin=0, centMax=5;
  // TString fileData = "HistsCorr.root";
  // TString corrStrr = "Corr";
  //

  // dump all info to tree
  TString centEtaStr =  Form("cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centMin,centMax,etaMin,etaMax);
  TTreeSRedirector *debugFile = new TTreeSRedirector(Form("ParamTrees_%s_%s_setting_%d.root",corrStrr.Data(), centEtaStr.Data(),setting),"recreate");

  TFile *inputFile = TFile::Open(fileData);
  TH2D *h2Pi = (TH2D *)inputFile->Get(Form("h2DCleanPiKineCut%s_%s",corrStrr.Data(),centEtaStr.Data()));
  TH2D *h2El = (TH2D *)inputFile->Get(Form("h2DCleanElKineCut%s_%s",corrStrr.Data(),centEtaStr.Data()));
  if (setting==0 && centMin<5 && etaMin>-0.05 &&  etaMin<0.05 ) {
    histsStream -> GetFile()->cd();
    h2Pi->Write();
    h2El->Write();
  }
  //
  //
  Int_t ptDownBin = 0;
  Int_t ptUpBin   = 0;
  Double_t ptStep = 0.;
  Double_t ptDown = 0.;
  Double_t ptUp   = 0.;
  Double_t pt     = 0.;
  for (Int_t islice = 0; islice<nSlice; islice++){

    ptStep = ((ptMax-ptMin)/(Double_t)nSlice);
    ptDown = ptMin+islice*ptStep;
    ptUp   = ptDown+ptStep;
    pt     = (ptDown+ptUp)/2.;

    // find bins to be used in the projection
    ptDownBin = (h2Pi->GetXaxis()->FindBin(ptDown+0.001));                                // TO FIX
    ptUpBin   = (h2Pi->GetXaxis()->FindBin(ptUp+0.001));
    momBin = ptUpBin;
    //
    TH1D * h1DPi = h2Pi->ProjectionY(Form("h1DPi_%s_%s_slice_%d",corrStrr.Data(),centEtaStr.Data(),islice),  ptDownBin,ptDownBin);
    TH1D * h1DEl = h2El->ProjectionY(Form("h1DEl_%s_%s_slice_%d",corrStrr.Data(),centEtaStr.Data(),islice),  ptDownBin,ptDownBin);
    Double_t binWidth = h1DPi->GetXaxis()->GetBinWidth(50);
    h1DPi->GetXaxis()->SetRangeUser(0,150);
    h1DEl->GetXaxis()->SetRangeUser(0,150);
    h1DPi->Scale(1./binWidth);
    h1DEl->Scale(1./binWidth);

    Double_t * elParClean = new Double_t[nPars];
    Double_t * piParClean = new Double_t[nPars];
    Double_t * elParError = new Double_t[nPars];
    Double_t * piParError = new Double_t[nPars];

    for (Int_t ipar=0; ipar<nPars; ipar++) {
      elParClean[ipar]  = 0;
      piParClean[ipar]  = 0;
      elParError[ipar]  = 0;
      piParError[ipar]  = 0;
    }
    //
    FitParticleSample(h1DPi, piParClean, elParError);
    FitParticleSample(h1DEl, elParClean, piParError);
    if (setting==0 && centMin<5 && etaMin>-0.05 &&  etaMin<0.05 ) {
      histsStream -> GetFile()->cd();
      h1DPi->Write();
      h1DEl->Write();
    }
    //
    // dump params
    debugFile -> GetFile()->cd();
    *debugFile << "params" <<
    "setting="      << setting      <<
    "corr="         << corrType     <<
    "pt="           << pt           <<
    "etaMin="       << etaMin       <<
    "etaMax="       << etaMax       <<
    "centMin="      << centMin      <<
    "centMax="      << centMax      <<
    "elamp="        << elParClean[0] <<
    "elmean="       << elParClean[1] <<
    "elsigma="      << elParClean[2] <<
    "elskew="       << elParClean[3] <<
    "elkurt="       << elParClean[4] <<
    "elint="        << elParClean[5] <<
    "elchi2="       << elParClean[6] <<
    "elmax="        << elParClean[7] <<
    "elmaxbincenter=" << elParClean[9] <<
    //
    "piamp="        << piParClean[0] <<
    "pimean="       << piParClean[1] <<
    "pisigma="      << piParClean[2] <<
    "piskew="       << piParClean[3] <<
    "pikurt="       << piParClean[4] <<
    "piint="        << piParClean[5] <<
    "pichi2="       << piParClean[6] <<
    "pimax="        << piParClean[7] <<
    "pimaxhist="    << piParClean[8] <<
    "pimaxbincenter=" << piParClean[9] <<
    //
    "elampErr="        << elParError[0] <<
    "elmeanErr="       << elParError[1] <<
    "elsigmaErr="      << elParError[2] <<
    "elskewErr="       << elParError[3] <<
    "elkurtErr="       << elParError[4] <<
    //
    "piampErr="        << piParError[0] <<
    "pimeanErr="       << piParError[1] <<
    "pisigmaErr="      << piParError[2] <<
    "piskewErr="       << piParError[3] <<
    "pikurtErr="       << piParError[4] <<
    "\n";

    delete [] piParClean;
    delete [] elParClean;
    delete [] elParError;
    delete [] piParError;

  }

  delete debugFile;
  inputFile->Close();

}
// -------------------------------------------------------------------------------------------------------
void FitParticleSample(TH1D *hClean, Double_t *parSample, Double_t *parError)
{

  //
  // Fit Clean sample slice
  //

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter::SetMaxIterations(100000);
  maxBin = hClean->GetBinContent(hClean->GetMaximumBin());
  Double_t mean  = hClean->GetMean();
  Double_t maxBin = hClean->GetBinContent(hClean->GetMaximumBin());
  Double_t maxBinCenter = hClean->GetBinCenter(hClean->GetMaximumBin());
  Double_t rms  = hClean->GetRMS();
  Double_t sampleFitWindow = 1.5;
  Double_t integral = hClean->Integral(hClean->FindBin(maxBinCenter-sampleFitWindow*rms),hClean->FindBin(maxBinCenter+sampleFitWindow*rms));

  Bool_t histCheck = ( (mean>20 && mean<1020) && (rms>0 && rms<200) && (maxBin>=1 && maxBin<1e+30) && (integral>=5 && integral<1e+30));
  if (!histCheck) {   // there is too few entries, fit would fail
    for (Int_t ipar=0; ipar<nPars; ipar++) {
      parSample[ipar] = 1;
      parError[ipar]  = 0;
    }
    return;
  }

  // fit function
  TF1 *asyGaus = new TF1("asyGaus",fitFunctionGenGaus,maxBinCenter-sampleFitWindow*rms,maxBinCenter+sampleFitWindow*1.2*rms);
  asyGaus->SetParNames("Amplitude","Mean","Sigma","Kurtosis","Skewness");
  asyGaus ->SetLineColor(2);
  asyGaus ->SetLineWidth(2);

  // Set Parameters
  asyGaus->SetParLimits(0,maxBin/20.,maxBin);
  asyGaus->SetParLimits(1,maxBinCenter-rms,maxBinCenter+rms);
  asyGaus->SetParLimits(2,rms/2. ,3.*rms);
  asyGaus->SetParLimits(3,kurtosisMin,kurtosisMax);
  asyGaus->SetParLimits(4,skewMin,skewMax);
  // asyGaus->FixParameter(3, 0.);
  // asyGaus->FixParameter(4, 2.);

  // retrieve fit parameters
  hClean->Fit(asyGaus,"QNMR");
  hClean->GetListOfFunctions()->Add(asyGaus);
  parSample[0] = (asyGaus->GetParameter(0)>1.) ? (asyGaus->GetParameter(0)) : 1.;      // avoid floating point exception for muons
  parSample[1] = asyGaus->GetParameter(1);
  parSample[2] = asyGaus->GetParameter(2);
  parSample[3] = asyGaus->GetParameter(3);
  parSample[4] = asyGaus->GetParameter(4);
  parSample[5] = asyGaus->Integral(parSample[1]-parSample[2]*4.,parSample[1]+parSample[2]*4.);
  if (asyGaus->GetNDF()>1e-4) parSample[6]=asyGaus->GetChisquare()/asyGaus->GetNDF();
  parSample[7] = asyGaus->GetMaximumX();
  parSample[8] = maxBin;
  parSample[9] = maxBinCenter;
  //
  parError[0] = asyGaus->GetParError(0);
  parError[1] = asyGaus->GetParError(1);
  parError[2] = asyGaus->GetParError(2);
  parError[3] = asyGaus->GetParError(3);
  parError[4] = asyGaus->GetParError(4);

  //
  //
  // Plot lines for electrons
  TLine* line = new TLine(maxBinCenter, 0., maxBinCenter, maxBin);
  line->SetLineColor(4);
  line->SetLineWidth(2);
  if (maxBinCenter>0 && maxBinCenter<1020) hClean->GetListOfFunctions()->Add(line);
  //
  // Add the legend
  TLegend *leg = new TLegend(0.65, 0.55, 0.85, 0.85);
  leg->SetTextFont(62); leg->SetTextSize(0.03); leg->SetFillColor(0); leg->SetNColumns(2);
  leg->AddEntry((TObject*)0,"#chi^{2}","");
  if (asyGaus->GetNDF()>0) leg->AddEntry((TObject*)0,Form(" = %4.2f",asyGaus->GetChisquare()/asyGaus->GetNDF()),"");
  else leg->AddEntry((TObject*)0,"0","");
  leg->AddEntry((TObject*)0,"Abundance","");
  leg->AddEntry((TObject*)0,Form(" = %4.2f",asyGaus->GetParameter(0)),"");
  leg->AddEntry((TObject*)0,"Mean"     ,"");
  leg->AddEntry((TObject*)0,Form(" = %4.2f",asyGaus->GetParameter(1)),"");
  leg->AddEntry((TObject*)0,"Sigma"    ,"");
  leg->AddEntry((TObject*)0,Form(" = %4.2f",asyGaus->GetParameter(2)),"");
  leg->AddEntry((TObject*)0,"Kurtosis" ,"");
  leg->AddEntry((TObject*)0,Form(" = %4.2f",asyGaus->GetParameter(3)),"");
  leg->AddEntry((TObject*)0,"Skewness" ,"");
  leg->AddEntry((TObject*)0,Form(" = %4.2f",asyGaus->GetParameter(4)),"");
  hClean->GetListOfFunctions()->Add(leg);
  hClean->GetXaxis()->SetTitle("TPC d#it{E}/d#it{x} Signal (a.u.)");
  hClean->GetYaxis()->SetTitle("entries");


}
