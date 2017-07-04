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
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTreeStream.h"
//#include "/hera/alice/marsland/software/aliroot/trunk_private/STEER/STEERBase/TTreeStream.h"
#include "TGraphSmooth.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TError.h"
#include "TGraphErrors.h"
#include "TClonesArray.h"
#include "TLinearFitter.h"
#include "TVirtualFitter.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TF1.h"
#include "TLegend.h"
#include "AliXRDPROOFtoolkit.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::cout;
using std::setw;


// Initial Settings
Bool_t Systematic             = kFALSE;     // for the systematics study jump directly to the 6th iter
Bool_t lowPfit                = kFALSE;
Bool_t dump                   = kFALSE;
Bool_t MC                     = kFALSE;      // to analyse LHC11a10a_bis or LHC10h
Bool_t normalisedCleanSamples = kFALSE;
Bool_t automaticKS            = kTRUE;
Bool_t pp                     = kFALSE;


// For Fixed kurtosis and skewness
Bool_t fixedK   = kTRUE;       // incase kTRUE modify if in  " SetFixedKSparameterFor*()" functions 
Bool_t fixedS   = kTRUE;

Double_t skewnessFixFit[6]   = {0., 0., 0., 0., 0., 0.};
Double_t skewnessFixClean[6] = {0., 0., 0., 0., 0., 0.};
Double_t skewnessFixMC[6]    = {0., 0., 0., 0., 0., 0.};
Double_t skewnessFix[6]      = {0., 0., 0., 0., 0., 0.};

Double_t kurtosisFixFit[6]   = {2., 2., 2., 2., 2., 2.};
Double_t kurtosisFixClean[6] = {2., 2., 2., 2., 2., 2.};
Double_t kurtosisFixMC[6]    = {2., 2., 2., 2., 2., 2.};
Double_t kurtosisFix[6]      = {2., 2., 2., 2., 2., 2.};


// For free Kurtosis and skewness
Double_t skewMin       = 0.;
Double_t skewMax       = 2.;
Double_t kurtosisMin   = 1.5;
Double_t kurtosisMax   = 2.5;

Int_t assump           = 0;        // choose a value 0,1,2
Int_t particleSign     = 0;        // choose a value 0,1,2 --> 0; all particles, 1; positive particles, 2; negative particles
Int_t freedomForPion   = 6;        // up to freedomForPion iteration fix all particle Kurtosis and skewness then let pion KSfree
TString smoothFunction = "pol2";   // outlier removal fit function to be used in TLinearFitter
Double_t windowSetting = 4.;       // in window arrangment for mean "distToClosestParticle/windowSetting"
Double_t smoothSeg     = 5.;       // outlier removal range setting for the TLinearFitter (nSlice/smoothSeg)
Double_t smoothLevel   = 9.;       // Increases the iterations in smoothing (10 is too much) for TGraphSmooth
Double_t highPbinWidth = 0.1;      // high momentum bin width which is used after 2 GeV/c on
Double_t analysisRange = 5.;      // Maximum p value used in the analysis
Double_t highPskewness = 0.;
Double_t highPkurtosis = 1.92;

// Define amp mean and sigma parlimits
Double_t elMeanMin ,piMeanMin ,kaMeanMin ,prMeanMin ,deMeanMin ,muMeanMin;
Double_t elMeanMax ,piMeanMax ,kaMeanMax ,prMeanMax ,deMeanMax ,muMeanMax; 
Double_t elAmpMin  ,piAmpMin  ,kaAmpMin  ,prAmpMin  ,deAmpMin  ,muAmpMin;
Double_t elAmpMax  ,piAmpMax  ,kaAmpMax  ,prAmpMax  ,deAmpMax  ,muAmpMax;
Double_t elSigmaMin,piSigmaMin,kaSigmaMin,prSigmaMin,deSigmaMin,muSigmaMin;
Double_t elSigmaMax,piSigmaMax,kaSigmaMax,prSigmaMax,deSigmaMax,muSigmaMax;

// Defien max position in each expected particle mean bin
Double_t maxBin,maxBin0,maxBin1,maxBin2,maxBin3;
Double_t muAmpMCscaled;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Containers
TClonesArray fitResArr("TH1D", 1000);
TClonesArray fitResiduals("TH1D", 1000);
TObjArray    cleanResArr(1000);
TObjArray    MCResArr(1000);

// Debuggers and files
TTreeSRedirector *debugFile     = 0;
TTreeSRedirector *smoothResults = 0;
TTreeSRedirector *histsFile     = 0;
TTreeSRedirector *fitFile       = 0;
TTreeSRedirector *sampleFile    = 0;
TGraphSmooth *gs;


//   Some objects  
TObjArray    * Fit;          
TObjArray    * Clean;       
TObjArray    * MCs;           
TObjArray    * ScaledClean;  
TObjArray    * ScaledMC;     
TObjArray    * ScaledEx;      
TObjArray    * OutlierSmooth; 
TObjArray    * Windows;  
TObjArray    * HighPFit;
TGraphErrors * grpiToMuRatio; 

// helper graphs 
TGraphErrors *grelFitAmpSmooth=0x0, *grelFitMeanSmooth=0x0, *grelFitSigmaSmooth=0x0;
TGraphErrors *grpiFitAmpSmooth=0x0, *grpiFitMeanSmooth=0x0, *grpiFitSigmaSmooth=0x0;
TGraphErrors *grkaFitAmpSmooth=0x0, *grkaFitMeanSmooth=0x0, *grkaFitSigmaSmooth=0x0;
TGraphErrors *grprFitAmpSmooth=0x0, *grprFitMeanSmooth=0x0, *grprFitSigmaSmooth=0x0;

  // MC Scaled graphs
TGraphErrors *grpiMCAmpScaledSmooth=0x0, *grpiMCMeanScaledSmooth=0x0, *grpiMCSigmaScaledSmooth=0x0;
TGraphErrors *grkaMCAmpScaledSmooth=0x0, *grkaMCMeanScaledSmooth=0x0, *grkaMCSigmaScaledSmooth=0x0;
TGraphErrors *grprMCAmpScaledSmooth=0x0, *grprMCMeanScaledSmooth=0x0, *grprMCSigmaScaledSmooth=0x0;

  // Extra graphs
TGraphErrors *grelMeanExScaled=0x0,               *grelSigmaExScaled=0x0,                *grelFitAmpOutlierSmooth=0x0;
TGraphErrors *grpiCleanMeanSmooth=0x0,            *grpiCleanSigmaScaledSmooth=0x0,       *grpiCleanMeanScaledSmooth=0x0, *grpiSigmaExScaled=0x0, *grpiMeanExScaled=0x0;
TGraphErrors *grkaCleanMean=0x0, *grkaMeanExScaled=0x0, *grkaSigmaExScaled=0x0, *grkaCleanSigma=0x0, *grkaCleanAmpScaledSmooth=0x0, *grkaFitAmpOutlierSmooth=0x0, *grkaFitSigmaOutlierSmooth=0x0, *grkaCleanSigmaScaled=0x0;
TGraphErrors *grprCleanMean=0x0, *grprMeanExScaled=0x0, *grprSigmaExScaled=0x0, *grprCleanSigma=0x0, *grprCleanAmpScaledSmooth=0x0, *grprCleanSigmaScaledSmooth=0x0, *grprFitAmpOutlierSmooth=0x0;
  
  // High momentum part smooth graphs
TGraphErrors *grpiFitAmpSmoothHP=0x0,  *grkaFitAmpSmoothHP=0x0,  *grprFitAmpSmoothHP=0x0;
TGraphErrors *grpiFitMeanSmoothHP=0x0, *grkaFitMeanSmoothHP=0x0, *grprFitMeanSmoothHP=0x0;
 
  // Muons are special
TGraphErrors *grmuMeanExScaled=0x0, *grmuSigmaExScaled=0x0;

  // For h2DMC analysis
TGraphErrors *grelMCSigmaScaled=0x0, *grpiMCSigmaScaled=0x0, *grkaMCSigmaScaled=0x0, *grprMCSigmaScaled=0x0, *grmuMCSigmaScaled=0x0; 
TGraphErrors *grelMCMeanScaled=0x0,  *grpiMCMeanScaled=0x0,  *grkaMCMeanScaled=0x0,  *grprMCMeanScaled=0x0,  *grmuMCMeanScaled=0x0;  
TGraphErrors *grelMCAmpScaled=0x0,   *grpiMCAmpScaled=0x0,   *grkaMCAmpScaled=0x0,   *grprMCAmpScaled=0x0,   *grmuMCAmpScaled=0x0;              
       
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

enum ParticleType {
  kElectron = 0,
  kPion     = 1,
  kKaon     = 2,
  kProton   = 3,
  kDeuteron = 4,
  kMuon     = 5
} pType;

TString fitFunctionGenGaus = "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";


Int_t dEdxMin        = 20;
Int_t dEdxMax        = 1020;
Double_t ptRangeDown = 0.2;
Double_t ptRangeUp   = 10.2;
Int_t dEdxNbins      = 4000;
Int_t ptNbins        = 2000;

Float_t etaBin;
Int_t centBin;

// Main Functions for iterative fitting procedure 
void          SmoothAmplitudes(TString fileFit, TString fileSample, TString fileOut, const Int_t nSlice, Int_t iIter);
void          GetExandSampleParameters(Int_t islice, TFile *sFile, Double_t *arrMean, Double_t *arrSigma, Double_t *cleanParMean, Double_t *MCPars);
void          GetHelperGraphs(Int_t iIter, TFile *ressFile);
void          GetCleanExParams(TH1D *hClean, Double_t *arrSigma, Double_t *arrCleanSigma, Double_t *arrCleanMean, ParticleType pSpecy);

void          EstimateKS(Int_t iks, TH2D *h2DKS, TH2D **hCleanSamples, TH2D **hMCSamples, TH2D **hExpectedMean, TH2D **hExpectedSigma);
void          FitAllSamples(const Int_t nSlice, const Double_t ptMin, const Double_t ptMax,TH2D *h2D, TH2D **hCleanSamples, TH2D **hMCSamples, TH2D **hExpectedMean, TH2D **hExpectedSigma, TString sampleFileName);
void          IterativeFitting(Int_t iIter, const Int_t nSlice, const Double_t ptMin, const Double_t ptMax, TH2D *h2D, TString fileIn, TString readSmooth, TString readSamp);
void          ProduceGraphs(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grCleanTmp, TGraph **grCleanSmoothTmp);
void          ProduceHighPGraphs(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grTmp, TGraph **grSmoothTmp, TGraphErrors **grAllFit);
void          ApplyOutlierSmooth(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grTmp, TGraph **grSmoothTmp);

// Helper functions for severel computations
void          Fit1DSliceForKSEstimate(Int_t iks, TH1D *h1DKS, Double_t *pars, Double_t *arrMean, Double_t *arrSigma);
void          FitParticleSampleForKSEstimate(Int_t iks, TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy, Int_t sampleType);
TMatrixD     *FitParticleSample(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy, Int_t sampleType);
TF1          *FitParamInRange(TGraphErrors* gr, TString funcType, Double_t min, Double_t max);
TH1D         *FuncToHist(TF1 *f, TH1D * h, Int_t islice, Double_t *arrMean, Double_t *arrSigma);
TH1D         *GraphShade(TGraphErrors *grFitParam, Double_t percent);
TGraphErrors *RemoveOutliers(TGraphErrors *grsmooth, Int_t fitWindow);
TGraphErrors *ApplyScalingWrtExpected(const Int_t nSlice, TString exp, TString fit, const Double_t safeP, TObjArray * arr, TTree *t);
TGraphErrors *ApplyScalingWrtSample(const Int_t nSlice, TString sample, TString fit, const Double_t safeP, TObjArray * arr, TObjArray * arrSample, TTree *tSample, TGraph **grScaledSmooth,Int_t k);
TGraphErrors *ProduceWindowGraphs(Int_t nSlice,Int_t index, TString str, TTree *tree);


Double_t      ComputeGaussIntegral(Double_t parAmp,Double_t parMean,Double_t parSigma,Double_t parKurtosis, Double_t parSkew);
TH1D         *GetClean1DSlice(TH2D *h2Clean, TString parName, Int_t ptDownBin, Int_t ptUpBin);
Int_t         GetMaxIndex(TGraphErrors *grIndex);
Int_t         GetPtIndex(TGraphErrors *grIndex,Double_t p);
Double_t      GetClosestParticleMean(Int_t islice, Double_t sig, Double_t ref, Double_t gr1, Double_t gr2, Double_t gr3);
TGraphErrors *GraphSmooth(TGraphErrors * gr);
TGraphErrors *HistToGraphErrors(TH1D * h);
TCanvas      *GetFitResCanvas(TObjArray *arr);
Int_t         CalculateNBinsP(Int_t nSlice, Int_t ptMin, Int_t ptMax);


// Set functions for the parameters
void          SetTotalFitParameters(Int_t iIter, Int_t islice, Int_t nSlice, TF1* total, Double_t *arrMean, Double_t *arrSigma, Double_t *arrMeanWindow, Double_t *cleanParMean, Double_t *MCPars, TH1D* h1D, Double_t pt);
void          CheckIfMeanIsZero(TF1 *total, Double_t *arrMean);
void          SetParNamesOfTotalFit(TF1 *total);
void          SetParamsForKSEstimate(Int_t iks, TF1 *total, TH1D *h1DKS, Double_t *arrMean, Double_t *arrSigma);
void          SetFixedKSparameterForTotalFit(TF1 *total);
void          SetFixedKSparameterForTotalFitForPions(TF1 *total);
void          SetFixedKSparameterForCleanSampleFit(TF1 *asyGaus, ParticleType pSpecy);
void          SetFixedKSparameterForMCSampleFit(TF1 *asyGaus, ParticleType pSpecy);
void          SetFixedKSparameterForIndividualFits(TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4, TF1 *g5, TF1 *g6);
void          ResetParametersOfEachFunction(TF1 *g1,TF1 *g2,TF1 *g3,TF1 *g4,TF1 *g5,TF1 *g6,TF1 *total,TH1D *h1D,const Double_t pt);
void          SetParamsMC(Int_t islice, Double_t pt);
void          SetParams1stIteration(Double_t pt, TH1D *h1D, Double_t *arrMean, Double_t *arrSigma,Double_t *arrMeanWindow, Double_t *cleanParMean);
void          SetParams1stIterationMC(Double_t pt, TH1D *h1D, Double_t *arrMean, Double_t *arrSigma, Double_t * MCPars);
void          SetHighPFitParams(TH1D *h1D, Double_t *arrMean, Double_t *arrSigma,Double_t *arrMeanWindow);

// Fit Function
// Double_t fitFunctionGenGaus(Double_t *x, Double_t *par)
// {
//   //
//   // Generalised gauss function --> "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
//   // Skew-normal distribution   --> "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]/TMath::Sqrt(2)))";
//   // par[0] --> Amplitude
//   // par[1] --> Mean
//   // par[2] --> Sigma
//   // par[3] --> Kurtosis
//   // par[4] --> Skewness
//   //
//   
//   Double_t fun = par[0]*exp(-TMath::Power((TMath::Abs(x[0]-par[1])/par[2]),par[3]))
//       *(1+TMath::Erf(par[4]*(x[0]-par[1])/par[2]/TMath::Sqrt(2)));
//   return fun;
// }
// -------------------------------------------------------------------------------------------------------
void PIDIterativeFittingMCandReal(TString fileData, const Int_t nSlice, Double_t ptMin, Double_t ptMax, Int_t maxIter=5){
    
  //
  // Apply gaussian fits to the ptot slices --> Main Function
  /*
  
  cd /hera/alice/marsland/pFluct/files/analysis/tests
  aliroot -l -b
  .L ~/PHD/macros/marsland_EbyeRatios/PIDIterativeFittingMCandReal.C+
  PIDIterativeFittingMCandReal("/hera/alice/marsland/pFluct/files/analysis/IdentityFiles/AllHists_PbPb_eta_0_0.2_cent_60_70.root",10,0.5,0.6,0);  2> err.log
  
  .L ~/PHD/macros/marsland_EbyeRatios/PIDIterativeFittingMCandReal.C+
  PIDIterativeFittingMCandReal("/hera/alice/marsland/pFluct/files/analysis/rootFilesPbPb/fullStat/ForAnar/PbPb_MC/PbPb_MergedHists/AllHists_PbPb_eta_0_0.2_cent_0_5.root",10,0.5,0.6,0);  2> err.log
  */
  //
    
//   gStyle->SetOptStat(0);   // remove statistics from the plot and check if the range is given correctly
  
  // Get eta and cent info
  TObjArray *objArr1  = fileData.Tokenize("/");
//   TString objFileName = (TString)((objArr1->At((Int_t)objArr1->GetLast()))->GetName());
  TString objFileName = ((objArr1->At((Int_t)objArr1->GetLast()))->GetName());

  TObjArray *objArr2  = objFileName.Tokenize("_");
  etaBin      = atof((objArr2->At(3))->GetName());
  centBin     = atoi((objArr2->At(6))->GetName());
  
  // Copy files for backup
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/PIDIterativeFittingMCandReal.C .");
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllHistsOnFarm.C .");
    
  // Create an debug file
  debugFile = new TTreeSRedirector(Form("debugFile_%d_%3.1f_%3.1f_%2.1f_%d.root",nSlice,ptMin,ptMax,etaBin,centBin),"recreate");
  
  // Cache histograms in memory 
  TH2D ** hCleanSamples   = new TH2D *[6];
  TH2D ** hMCSamples      = new TH2D *[6];
  TH2D ** hExpectedMean   = new TH2D *[6];
  TH2D ** hExpectedSigma  = new TH2D *[6];
  
  for (Int_t icache=0; icache<6; icache++)
  {
    hMCSamples[icache]    = NULL;
    hCleanSamples[icache] = NULL;
    hExpectedMean[icache] = NULL;
    hExpectedSigma[icache]= NULL;
  }
  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  cout << " ================================== "    << endl;
  cout << " fixedK         = " << fixedK            << endl;
  cout << " ================================== "    << endl;
  cout << " fixedS         = " << fixedS            << endl;
  cout << " ================================== "    << endl;
  cout << " MC             = " << MC                << endl;
  cout << " ================================== "    << endl;
  cout << " eta            = " << etaBin            << endl;
  cout << " ================================== "    << endl;
  cout << " cent           = " << centBin           << endl;
  cout << " ================================== "    << endl;
  cout << " file Path      = " << fileData          << endl;
  cout << " ================================== "    << endl;
  cout << " file name      = " << objFileName       << endl;
  cout << " ================================== "    << endl;
  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // open file and retrive all histograms to be processed
  TFile *f = TFile::Open(fileData);
  
  // Take main 2D hists
  cout << " take Files " << endl;
  TH2D  * h2D;
  
  switch (particleSign) 
  { 
    case 0: 
      h2D = (MC) ? (TH2D *)f->Get("h2DMCall") : (TH2D *)f->Get("h2Dall");
      break; 
    case 1: 
      h2D = (MC) ? (TH2D *)f->Get("h2DMCpos") : (TH2D *)f->Get("h2Dpos");
      break; 
    case 2: 
      h2D = (MC) ? (TH2D *)f->Get("h2DMCneg") : (TH2D *)f->Get("h2Dneg"); 
      break;  
  } 
  
  hCleanSamples[0]  = (TH2D *)f->Get("h2CleanElectron");
  hCleanSamples[1]  = (TH2D *)f->Get("h2CleanPion");
  hCleanSamples[2]  = (TH2D *)f->Get("h2CleanKaon");
  hCleanSamples[3]  = (TH2D *)f->Get("h2CleanProton");
  hCleanSamples[4]  = (TH2D *)f->Get("h2CleanPionXTOF");
  hMCSamples[0]     = (TH2D *)f->Get("h2MCelectron");
  hMCSamples[1]     = (TH2D *)f->Get("h2MCpion");
  hMCSamples[2]     = (TH2D *)f->Get("h2MCkaon");
  hMCSamples[3]     = (TH2D *)f->Get("h2MCproton");
  hMCSamples[5]     = (TH2D *)f->Get("h2MCmuon");
  hExpectedMean[0]  = (TH2D *)f->Get("h2ExpectedEl");
  hExpectedMean[1]  = (TH2D *)f->Get("h2ExpectedPi");
  hExpectedMean[2]  = (TH2D *)f->Get("h2ExpectedKa");
  hExpectedMean[3]  = (TH2D *)f->Get("h2ExpectedPr");
  hExpectedMean[4]  = (TH2D *)f->Get("h2ExpectedDe");
  hExpectedMean[5]  = (TH2D *)f->Get("h2ExpectedMu");
  hExpectedSigma[0] = (TH2D *)f->Get("h2ExpectedSigmaEl");
  hExpectedSigma[1] = (TH2D *)f->Get("h2ExpectedSigmaPi");
  hExpectedSigma[2] = (TH2D *)f->Get("h2ExpectedSigmaKa");
  hExpectedSigma[3] = (TH2D *)f->Get("h2ExpectedSigmaPr");
  hExpectedSigma[4] = (TH2D *)f->Get("h2ExpectedSigmaDe");
  hExpectedSigma[5] = (TH2D *)f->Get("h2ExpectedSigmaMu");
 
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cout << " =========== Write histograms into the debug file ============ " << endl;
  debugFile -> GetFile()->cd();  
  h2D               -> Write("h2DReal");
  hCleanSamples[0]  -> Write("h2Electron");
  hCleanSamples[1]  -> Write("h2Pion");
  hCleanSamples[2]  -> Write("h2Kaon");
  hCleanSamples[3]  -> Write("h2Proton");
  hCleanSamples[4]  -> Write("h2PionTOF");
  hMCSamples[0]     -> Write("h2MCelectron");
  hMCSamples[1]     -> Write("h2MCpion");
  hMCSamples[2]     -> Write("h2MCkaon");
  hMCSamples[3]     -> Write("h2MCproton");
  hMCSamples[5]     -> Write("h2MCmuon");
  hExpectedMean[0]  -> Write("h2ExpectedEl");
  hExpectedMean[1]  -> Write("h2ExpectedPi");
  hExpectedMean[2]  -> Write("h2ExpectedKa");
  hExpectedMean[3]  -> Write("h2ExpectedPr");
  hExpectedMean[4]  -> Write("h2ExpectedDe");
  hExpectedMean[5]  -> Write("h2ExpectedMu");

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  // for proper error calculations
  cout << " =========== Sumw2 ============ " << endl;
  h2D               -> Sumw2();
  hMCSamples[0]     -> Sumw2();
  hMCSamples[1]     -> Sumw2();
  hMCSamples[2]     -> Sumw2();
  hMCSamples[3]     -> Sumw2();
  hMCSamples[5]     -> Sumw2();
  hCleanSamples[0]  -> Sumw2();
  hCleanSamples[1]  -> Sumw2();
  hCleanSamples[2]  -> Sumw2();
  hCleanSamples[3]  -> Sumw2();
  hCleanSamples[4]  -> Sumw2();
   
  //  Estimate kurtosis and Skewness for kaon, proton and pion by fitting 100 slices in between [0.4,1.4]GeV/c
  TH2D *h2DKS = (TH2D*)h2D->Clone();
  if (!Systematic){
    for (Int_t iks = 0; iks<2; iks++){      
      cout << " =========== KS estimate  ============  Iteration = " << iks << endl;
      EstimateKS(iks,h2DKS,hCleanSamples,hMCSamples,hExpectedMean,hExpectedSigma);
    }
  }
  delete h2DKS;
  delete debugFile;
  
  // First fit the samples 
  cout << " =========== first make sample file and clean them up =========== " << endl;
  TString fileSample  = Form("Samples_%d_%3.1f_%3.1f_%2.1f_%d.root",nSlice,ptMin,ptMax,etaBin,0);      // care only eta dependence
  FitAllSamples(nSlice,ptMin,ptMax,h2D,hCleanSamples,hMCSamples,hExpectedMean,hExpectedSigma,fileSample);
 
  // delete unnecessary histograms in memory 
  for (Int_t icache=0; icache<6; icache++)
  {
    delete hMCSamples[icache];
    delete hCleanSamples[icache];
    delete hExpectedMean[icache];
    delete hExpectedSigma[icache];
  }
  delete [] hCleanSamples;
  delete [] hMCSamples;
  delete [] hExpectedMean;
  delete [] hExpectedSigma;
  
  // iteratively apply fits
  for (Int_t iIter = 0; iIter<10; iIter++){
    if(iIter>maxIter) break;
    if( MC && (iIter == 2 || iIter == 3 || iIter == 4 || iIter == 5)) continue;
    TString fileFit    = Form("Trees_Iter%d_%d_%3.1f_%3.1f_%2.1f_%d.root",iIter,nSlice,ptMin,ptMax,etaBin,centBin);
    TString fileOut    = Form("Results_Iter%d_%d_%3.1f_%3.1f_%2.1f_%d.root",iIter,nSlice,ptMin,ptMax,etaBin,centBin);
    Int_t previousIter = (MC && iIter==6) ? 1 : iIter-1;
    TString readSmooth = Form("Results_Iter%d_%d_%3.1f_%3.1f_%2.1f_%d.root",previousIter,nSlice,ptMin,ptMax,etaBin,centBin);
    
    cout << " =========== Iterative fit procedure is being started  ============  Iteration = " << iIter << endl;
    IterativeFitting(iIter,nSlice,ptMin,ptMax,h2D,fileFit,readSmooth,fileSample);
    
    cout << " =========== Iteration = " << iIter << "   has finished go to iteration = " << iIter+1 << endl;
    if (nSlice>=40) SmoothAmplitudes(fileFit,fileSample,fileOut,nSlice,iIter);
  }
  
  if (maxIter>1) gSystem->Exec(Form("hadd AllTrees_%d_%3.1f_%3.1f_%2.1f_%d.root Trees_*.root",nSlice,ptMin,ptMax,etaBin,centBin));
  
  delete objArr1;
  delete objArr2;
  f -> Close();
    
}
// -------------------------------------------------------------------------------------------------------
void SmoothAmplitudes(TString fileFit, TString fileSample, TString fileOut, const Int_t nSlice, Int_t iIter){

  //
  // Analyse first results and smooth the amplitude graphs
  /*
  cd /hera/alice/marsland/pFluct/files/analysis/tests/trash
  aliroot -l 
  .L ~/PHD/macros/marsland_EbyeRatios/PIDIterativeFittingMCandReal.C+  
  SmoothAmplitudes("Trees_Iter4_300_0.2_3.2_-0.2_0.root","Samples_300_0.2_3.2_-0.2_0.root","SmoothTest.root",300,0)
  */
  //
    
  smoothResults = new TTreeSRedirector(fileOut,"recreate");
    
//   Open "Trees_Iteration0_200_0.2_2.2_0.0_-2.root" file and read trees
  TFile *ftrees = TFile::Open(fileFit);
  TFile *strees = TFile::Open(fileSample);
  
  TTree *treeFit   = (TTree*)ftrees->Get("FitResults");   
  TTree *treeClean = (TTree*)strees->Get("CleanSamples");
  TTree *treeMC    = (TTree*)strees->Get("MCSamples");   
  
//   GetList of branches for each tree
  TObjArray *branchArrFit   =  (TObjArray*)treeFit  ->GetListOfBranches();
  TObjArray *branchArrClean =  (TObjArray*)treeClean->GetListOfBranches();
  TObjArray *branchArrMC    =  (TObjArray*)treeMC   ->GetListOfBranches();
  
  Int_t nArrEntries = branchArrFit->GetEntriesFast();

//   Output TObjArrays
  TObjArray *GrArrClean        = new TObjArray(2*nArrEntries);  GrArrClean         -> SetOwner(kTRUE); 
  TObjArray *GrArrMC           = new TObjArray(2*nArrEntries);  GrArrMC            -> SetOwner(kTRUE);
  TObjArray *GrArrFit          = new TObjArray(2*nArrEntries);  GrArrFit           -> SetOwner(kTRUE);
  TObjArray *GrArrScaledMC     = new TObjArray(2*nArrEntries);  GrArrScaledMC      -> SetOwner(kTRUE);
  TObjArray *GrArrScaledClean  = new TObjArray(2*nArrEntries);  GrArrScaledClean   -> SetOwner(kTRUE);
  TObjArray *GrArrScaledEx     = new TObjArray(2*nArrEntries);  GrArrScaledEx      -> SetOwner(kTRUE);
  TObjArray *GrArrHighPFit     = new TObjArray(2*nArrEntries);  GrArrHighPFit      -> SetOwner(kTRUE);
  
  TObjArray *GrArrOutlierSmooth = new TObjArray(2*nArrEntries); GrArrOutlierSmooth -> SetOwner(kTRUE);
  TObjArray *GrArrWindows       = new TObjArray(2*nArrEntries); GrArrWindows       -> SetOwner(kTRUE);

//   Graph Arrays
  TGraphErrors **grClean         = new TGraphErrors *[nArrEntries];
  TGraphErrors **grMC            = new TGraphErrors *[nArrEntries];
  TGraphErrors **grFit           = new TGraphErrors *[nArrEntries];
  TGraphErrors **grHigPFit       = new TGraphErrors *[nArrEntries];
  TGraphErrors **grAllFit        = new TGraphErrors *[nArrEntries];
  TGraphErrors **grOutlierCut    = new TGraphErrors *[nArrEntries];
  TGraph       **grFitSmooth     = new TGraph *[nArrEntries];
  TGraph       **grHigPFitSmooth = new TGraph *[nArrEntries];
  TGraph       **grCleanSmooth   = new TGraph *[nArrEntries];
  TGraph       **grMCSmooth      = new TGraph *[nArrEntries];
  TGraph       **grOutlierSmooth = new TGraph *[nArrEntries];

  TGraphErrors **grScaledMC          = new TGraphErrors *[nArrEntries];
  TGraph       **grScaledMCSmooth    = new TGraph *[nArrEntries];
  TGraphErrors **grScaledClean       = new TGraphErrors *[nArrEntries];
  TGraph       **grScaledCleanSmooth = new TGraph *[nArrEntries];

  TGraphErrors **grScaledEx          = new TGraphErrors *[nArrEntries];
  
  TGraphErrors * grWindowsAroundMean[nArrEntries];

//   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ProduceHighPGraphs(nSlice, treeFit, GrArrHighPFit, branchArrFit, grHigPFit, grHigPFitSmooth, grAllFit);
  
//   Produce smooth and normal graphs
  ProduceGraphs(nSlice, treeFit  , GrArrFit  , branchArrFit  , grFit  , grFitSmooth);  
  ProduceGraphs(nSlice, treeClean, GrArrClean, branchArrClean, grClean, grCleanSmooth); 
  ProduceGraphs(nSlice, treeMC   , GrArrMC   , branchArrMC   , grMC   , grMCSmooth);    
  
//   Apply outlier removal to fits
  if (!MC) ApplyOutlierSmooth(nSlice, treeFit, GrArrOutlierSmooth, branchArrFit, grOutlierCut, grOutlierSmooth);   
  
//   Apply Scalings for expected mean and sigma
  grScaledEx[0] = ApplyScalingWrtExpected(nSlice, "elMean","elFitMean",0.32,GrArrFit,treeFit);
  grScaledEx[1] = ApplyScalingWrtExpected(nSlice, "piMean","piFitMean",0.6,GrArrFit,treeFit);
  grScaledEx[2] = ApplyScalingWrtExpected(nSlice, "kaMean","kaFitMean",0.4,GrArrFit,treeFit);
  grScaledEx[3] = ApplyScalingWrtExpected(nSlice, "prMean","prFitMean",0.6,GrArrFit,treeFit);
  grScaledEx[4] = ApplyScalingWrtExpected(nSlice, "muMean","muFitMean",0.6,GrArrFit,treeFit);
  grScaledEx[5] = ApplyScalingWrtExpected(nSlice, "elSigma","elFitSigma",0.4,GrArrFit,treeFit);
  grScaledEx[6] = ApplyScalingWrtExpected(nSlice, "piSigma","piFitSigma",0.6,GrArrFit,treeFit);
  grScaledEx[7] = ApplyScalingWrtExpected(nSlice, "kaSigma","kaFitSigma",0.4,GrArrFit,treeFit);
  grScaledEx[8] = ApplyScalingWrtExpected(nSlice, "prSigma","prFitSigma",0.6,GrArrFit,treeFit);
  grScaledEx[9] = ApplyScalingWrtExpected(nSlice, "muSigma","muFitSigma",0.6,GrArrFit,treeFit);
  for (Int_t i=0; i<10; i++) GrArrScaledEx->AddAt(grScaledEx[i],i);

//   Apply Scalings for MC data mean and sigma
  grScaledMC[0]  = ApplyScalingWrtSample(nSlice,"elMCAmp","elFitAmp",0.32,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,0);
  grScaledMC[1]  = ApplyScalingWrtSample(nSlice,"piMCAmp","piFitAmp",0.6,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,1);
  grScaledMC[2]  = ApplyScalingWrtSample(nSlice,"kaMCAmp","kaFitAmp",0.7,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,2);
  grScaledMC[3]  = ApplyScalingWrtSample(nSlice,"prMCAmp","prFitAmp",1.25,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,3);
  grScaledMC[4]  = ApplyScalingWrtSample(nSlice,"elMCMean","elFitMean",0.32,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,4);
  grScaledMC[5]  = ApplyScalingWrtSample(nSlice,"piMCMean","piFitMean",0.55,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,5);
  grScaledMC[6]  = ApplyScalingWrtSample(nSlice,"kaMCMean","kaFitMean",0.32,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,6);
  grScaledMC[7]  = ApplyScalingWrtSample(nSlice,"prMCMean","prFitMean",0.6,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,7);
  grScaledMC[8]  = ApplyScalingWrtSample(nSlice,"elMCSigma","elFitSigma",0.32,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,8);
  grScaledMC[9]  = ApplyScalingWrtSample(nSlice,"piMCSigma","piFitSigma",0.55,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,9);
  grScaledMC[10] = ApplyScalingWrtSample(nSlice,"kaMCSigma","kaFitSigma",0.32,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,10);
  grScaledMC[11] = ApplyScalingWrtSample(nSlice,"prMCSigma","prFitSigma",0.6,GrArrFit,GrArrMC,treeMC,grScaledMCSmooth,11);
  
  for (Int_t i=0; i<12; i++) {
    GrArrScaledMC->AddAt(grScaledMC[i],2*i);
    GrArrScaledMC->AddAt(grScaledMCSmooth[i],2*i+1);
  }
  
//   Apply Scalings for MC data mean and sigma KAON and PROTON AMPLITUDE !!!!!!
  grScaledClean[0]  = ApplyScalingWrtSample(nSlice,"elCleanAmp","elFitAmp"  ,0.4,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,0);
  grScaledClean[1]  = ApplyScalingWrtSample(nSlice,"piCleanAmp","piFitAmp"  ,0.6,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,1);
  grScaledClean[2]  = ApplyScalingWrtSample(nSlice,"kaCleanAmp","kaFitAmp"  ,0.65,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,2);
  grScaledClean[3]  = ApplyScalingWrtSample(nSlice,"prCleanAmp","prFitAmp"  ,1.25,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,3);
  
  grScaledClean[4]  = ApplyScalingWrtSample(nSlice,"elCleanMean","elFitMean",0.4,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,4);
  grScaledClean[5]  = ApplyScalingWrtSample(nSlice,"piCleanMean","piFitMean",0.6,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,5);
  grScaledClean[6]  = ApplyScalingWrtSample(nSlice,"kaCleanMean","kaFitMean",0.65,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,6);
  grScaledClean[7]  = ApplyScalingWrtSample(nSlice,"prCleanMean","prFitMean",0.8,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,7);
  
  grScaledClean[8]  = ApplyScalingWrtSample(nSlice,"elCleanSigma","elFitSigma",0.4,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,8);
  grScaledClean[9]  = ApplyScalingWrtSample(nSlice,"piCleanSigma","piFitSigma",0.6,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,9);
  grScaledClean[10] = ApplyScalingWrtSample(nSlice,"kaCleanSigma","kaFitSigma",0.65,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,10);
  grScaledClean[11] = ApplyScalingWrtSample(nSlice,"prCleanSigma","prFitSigma",0.75,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,11);
  
  for (Int_t i=0; i<12; i++) {
    GrArrScaledClean->AddAt(grScaledClean[i],2*i);
    GrArrScaledClean->AddAt(grScaledCleanSmooth[i],2*i+1);
  }
  
//   Make windows around mean of particles 
  grWindowsAroundMean[0]  = ProduceWindowGraphs(nSlice, 1  ,"(elMean-elWindow):p",treeFit);
  grWindowsAroundMean[1]  = ProduceWindowGraphs(nSlice, 2  ,"(elMean+elWindow):p",treeFit);
  grWindowsAroundMean[2]  = ProduceWindowGraphs(nSlice, 3  ,"elMean:p"           ,treeFit);
  grWindowsAroundMean[3]  = ProduceWindowGraphs(nSlice, 4  ,"(piMean-piWindow):p",treeFit);
  grWindowsAroundMean[4]  = ProduceWindowGraphs(nSlice, 5  ,"(piMean+piWindow):p",treeFit);
  grWindowsAroundMean[5]  = ProduceWindowGraphs(nSlice, 6  ,"piMean:p"           ,treeFit);
  grWindowsAroundMean[6]  = ProduceWindowGraphs(nSlice, 7  ,"(kaMean-kaWindow):p",treeFit);
  grWindowsAroundMean[7]  = ProduceWindowGraphs(nSlice, 8  ,"(kaMean+kaWindow):p",treeFit);
  grWindowsAroundMean[8]  = ProduceWindowGraphs(nSlice, 9  ,"kaMean:p"           ,treeFit);
  grWindowsAroundMean[9]  = ProduceWindowGraphs(nSlice, 10 ,"(prMean-prWindow):p",treeFit);
  grWindowsAroundMean[10] = ProduceWindowGraphs(nSlice, 11 ,"(prMean+prWindow):p",treeFit);
  grWindowsAroundMean[11] = ProduceWindowGraphs(nSlice, 12 ,"prMean:p"           ,treeFit);
  for (Int_t i=0; i<12; i++) GrArrWindows->AddAt(grWindowsAroundMean[i],i);
  

// +++++++++++++++++++++++++++++   Some special graphs to be saved  +++++++++++++++++++++++++++++
//   piToMuRatio from MC
  treeMC->Draw("piMCAmp/muMCAmp:p","","goff");
  TGraphErrors * piToMuRatio = new TGraphErrors(nSlice,treeMC->GetV2(),treeMC->GetV1());
  piToMuRatio->SetName("piToMuRatio");
  
//   Chi2 of the total fit
  treeFit->Draw("totalChi2/normNDF:p","","goff");
  TGraphErrors *totChi2 = new TGraphErrors(nSlice,treeFit->GetV2(),treeFit->GetV1());
  totChi2->SetMarkerStyle(7);
  
  treeFit->Draw("maxRes:p","","goff");
  TGraphErrors *maxRes = new TGraphErrors(nSlice,treeFit->GetV2(),treeFit->GetV1());
  maxRes->SetMarkerStyle(7);
  
  treeFit->Draw("maxPer:p","","goff");
  TGraphErrors *maxPer = new TGraphErrors(nSlice,treeFit->GetV2(),treeFit->GetV1());
  maxPer->SetMarkerStyle(7);
  
//   Canvas of All Fit Results
  TCanvas *canResults = GetFitResCanvas(GrArrFit);
  canResults->SaveAs(Form("summary_iter%d_eta%2.1f_cent%d.png",iIter,etaBin,centBin));
  
//   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  smoothResults->GetFile()->cd();
  // Helper graphs
  canResults         -> Write("AllFitParams");
  piToMuRatio        -> Write("piToMuRatio");
  totChi2            -> Write("totChi2");
  maxPer             -> Write("maxDeviation");
  maxRes             -> Write("maxResidual");

  // Arrays
  GrArrFit           -> Write("Fit"         ,TObject::kSingleKey);
  GrArrClean         -> Write("Clean"       ,TObject::kSingleKey);
  GrArrMC            -> Write("MC"          ,TObject::kSingleKey);
  GrArrScaledEx      -> Write("ScaledEx"    ,TObject::kSingleKey);
  GrArrScaledMC      -> Write("ScaledMC"    ,TObject::kSingleKey);
  GrArrScaledClean   -> Write("ScaledClean" ,TObject::kSingleKey);
  GrArrWindows       -> Write("Windows"     ,TObject::kSingleKey);
  GrArrOutlierSmooth -> Write("OulierSmooth",TObject::kSingleKey);
  GrArrHighPFit      -> Write("HighPFit"    ,TObject::kSingleKey);
  
//   delete arrays
  GrArrFit           -> Delete();
  GrArrClean         -> Delete();
  GrArrMC            -> Delete();
  GrArrScaledEx      -> Delete();
  GrArrScaledMC      -> Delete();
  GrArrScaledClean   -> Delete();
  GrArrOutlierSmooth -> Delete();
  GrArrHighPFit      -> Delete();
  delete GrArrFit;
  delete GrArrClean;
  delete GrArrMC;
  delete GrArrScaledEx;
  delete GrArrScaledMC;
  delete GrArrScaledClean;
  delete GrArrOutlierSmooth;
  delete GrArrHighPFit;
  
//   delete helper graphs
  delete canResults;
  delete piToMuRatio;
  delete totChi2;
  delete maxRes;
  delete maxPer;
  
//   delete outputfile
  delete smoothResults;
  ftrees->Close();
  strees->Close();
  
}
// -------------------------------------------------------------------------------------------------------
void GetExandSampleParameters(Int_t islice, TFile *sFile, Double_t *arrMean, Double_t *arrSigma, Double_t *cleanParMean, Double_t *MCPars){

  //
  // Read Sample file
  //
  
  TTree *treeSample = (TTree*)sFile->Get("CleanSamples"); 
  TTree *treeMC     = (TTree*)sFile->Get("MCSamples"); 
 
  Double_t elMean     =0., elSigma = 0., elCleanMean = 0., elMCMean = 0., elMCSigma = 0.;
  Double_t piMean     =0., piSigma = 0., piCleanMean = 0., piMCMean = 0., piMCSigma = 0.;
  Double_t kaMean     =0., kaSigma = 0., kaCleanMean = 0., kaMCMean = 0., kaMCSigma = 0.;
  Double_t prMean     =0., prSigma = 0., prCleanMean = 0., prMCMean = 0., prMCSigma = 0.;
  Double_t deMean     =0., deSigma = 0.,                   deMCMean = 0., deMCSigma = 0.;
  Double_t muMean     =0., muSigma = 0.;
        
  // Get MC mean and sigma from MC tree
  treeMC->SetBranchAddress("elMCMean" ,&elMCMean);
  treeMC->SetBranchAddress("piMCMean" ,&piMCMean);
  treeMC->SetBranchAddress("kaMCMean" ,&kaMCMean);
  treeMC->SetBranchAddress("prMCMean" ,&prMCMean);
//   treeMC->SetBranchAddress("deMCMean" ,&deMCMean);
  
  treeMC->SetBranchAddress("elMCSigma" ,&elMCSigma);
  treeMC->SetBranchAddress("piMCSigma" ,&piMCSigma);
  treeMC->SetBranchAddress("kaMCSigma" ,&kaMCSigma);
  treeMC->SetBranchAddress("prMCSigma" ,&prMCSigma);
//   treeMC->SetBranchAddress("deMCSigma" ,&deMCSigma);
  
  treeMC->GetEntry(islice);
  MCPars[0] = elMCMean;
  MCPars[1] = piMCMean;
  MCPars[2] = kaMCMean;
  MCPars[3] = prMCMean;
  MCPars[4] = deMCMean;
  
  MCPars[5] = elMCSigma;
  MCPars[6] = piMCSigma;
  MCPars[7] = kaMCSigma;
  MCPars[8] = prMCSigma;
  MCPars[9] = deMCSigma;

  
  
  // get Clean Mean and expected mean-sigma from Clean Tree
  treeSample->SetBranchAddress("elMean" ,&elMean);
  treeSample->SetBranchAddress("piMean" ,&piMean);
  treeSample->SetBranchAddress("kaMean" ,&kaMean);
  treeSample->SetBranchAddress("prMean" ,&prMean);
  treeSample->SetBranchAddress("deMean" ,&deMean);
  treeSample->SetBranchAddress("muMean" ,&muMean);
  
  treeSample->SetBranchAddress("elSigma" ,&elSigma);
  treeSample->SetBranchAddress("piSigma" ,&piSigma);
  treeSample->SetBranchAddress("kaSigma" ,&kaSigma);
  treeSample->SetBranchAddress("prSigma" ,&prSigma);
  treeSample->SetBranchAddress("deSigma" ,&deSigma);
  treeSample->SetBranchAddress("muSigma" ,&muSigma);
  
  treeSample->SetBranchAddress("elCleanMean" ,&elCleanMean);
  treeSample->SetBranchAddress("piCleanMean" ,&piCleanMean);
  treeSample->SetBranchAddress("kaCleanMean" ,&kaCleanMean);
  treeSample->SetBranchAddress("prCleanMean" ,&prCleanMean);
  
  treeSample->GetEntry(islice);
  arrMean[0] = elMean;  arrSigma[0] = elSigma; 
  arrMean[1] = piMean;  arrSigma[1] = piSigma; 
  arrMean[2] = kaMean;  arrSigma[2] = kaSigma; 
  arrMean[3] = prMean;  arrSigma[3] = prSigma; 
  arrMean[4] = deMean;  arrSigma[4] = deSigma; 
  arrMean[5] = muMean;  arrSigma[5] = muSigma; 
  
  cleanParMean[0] = elCleanMean;
  cleanParMean[1] = piCleanMean;
  cleanParMean[2] = kaCleanMean;
  cleanParMean[3] = prCleanMean;
  
}
// -------------------------------------------------------------------------------------------------------
void EstimateKS(Int_t iks, TH2D *h2DKS, TH2D **hCleanSamples, TH2D **hMCSamples, TH2D **hExpectedMean, TH2D **hExpectedSigma){
  
  //
  // apply 100 fits to FIT, MC and CLEAN samples with 10 MeV steps in between [0.4, 1.4]
  // 
  
  TMinuit g;
  g.SetMaxIterations(10000);
  gMinuit->SetMaxIterations(10000);
  Int_t nbins = 100;
  Double_t pMin = 0.4;
  Double_t pMax = 1.4;
  
  TH1D hCleanKurtosis[6];
  TH1D hMCKurtosis[6];
  TH1D hFitKurtosis[6];
  TH1D hCleanSkewness[6];
  TH1D hMCSkewness[6];
  TH1D hFitSkewness[6];
  
  for (Int_t i = 0; i<6; i++){
    hCleanKurtosis[i] = TH1D(Form("hCleanKurtosis_%d",i),Form("hCleanKurtosis_%d",i),100,1.5,3.5);
    hMCKurtosis[i]    = TH1D(Form("hMCKurtosis_%d"   ,i),Form("hMCKurtosis_%d"   ,i),100,1.5,3.5);
    hFitKurtosis[i]   = TH1D(Form("hFitKurtosis_%d"  ,i),Form("hFitKurtosis_%d"  ,i),100,1.5,3.5);
    hCleanSkewness[i] = TH1D(Form("hCleanSkewness_%d",i),Form("hCleanSkewness_%d",i),100,0.,2.);
    hMCSkewness[i]    = TH1D(Form("hMCSkewness_%d"   ,i),Form("hMCSkewness_%d"   ,i),100,0.,2.);
    hFitSkewness[i]   = TH1D(Form("hFitSkewness_%d"  ,i),Form("hFitSkewness_%d"  ,i),100,0.,2.);
  }
  

  Int_t ptDownBin = 0;
  Int_t ptUpBin   = 0; 
  Double_t ptStep = 0.;
  Double_t ptDown = 0.;
  Double_t ptUp   = 0.;
  Double_t pt     = 0.;      
  for (Int_t islice = 0; islice<nbins; islice++){
        
    ptStep = ((pMax-pMin)/(Double_t)nbins);
    ptDown = pMin+islice*ptStep;
    ptUp   = ptDown+ptStep;
    pt     = (ptDown+ptUp)/2.;
  
    // find bins to be used in the projection
    ptDownBin = (h2DKS->GetXaxis()->FindBin(pMin)-1)+2*islice+1;                               
    ptUpBin   = (h2DKS->GetXaxis()->FindBin(pMin)-1)+2*islice+2;                               
    
    // initialise 1D histograms tobe fitted
    TH1D * h1CleanSamples[6];
    TH1D * h1MCSamples[6];
    TH1D * h1ExpectedMean[6];
    TH1D * h1ExpectedSigma[6];
         
     // Real particle parameters
    Double_t * arrSigma      = new Double_t[10];
    Double_t * arrMean       = new Double_t[10];
    Double_t * arrCleanSigma = new Double_t[10];
    Double_t * arrCleanMean  = new Double_t[10];
    
    // Clean Particle parameters
    Double_t * elParClean    = new Double_t[10];
    Double_t * piParClean    = new Double_t[10];
    Double_t * kaParClean    = new Double_t[10];
    Double_t * prParClean    = new Double_t[10];
    // MC Particle parameters
    Double_t * elParMC       = new Double_t[10];
    Double_t * piParMC       = new Double_t[10];
    Double_t * kaParMC       = new Double_t[10];
    Double_t * prParMC       = new Double_t[10];
    Double_t * muParMC       = new Double_t[10];
    
    // Total fit parameters
    Double_t * pars          = new Double_t[31]; for (Int_t j=0; j<31; j++) pars[j] = 0;
        
    for (Int_t ipar=0; ipar<10; ipar++)
    {
      arrMean[ipar]       = 0;
      arrSigma[ipar]      = 0;
      arrCleanMean[ipar]  = 0;
      arrCleanSigma[ipar] = 0;
      piParClean[ipar]    = 0;
      kaParClean[ipar]    = 0;
      prParClean[ipar]    = 0;
      elParClean[ipar]    = 0;
      piParMC[ipar]       = 0;
      kaParMC[ipar]       = 0;
      prParMC[ipar]       = 0;
      elParMC[ipar]       = 0;
      muParMC[ipar]       = 0;
    }
    
    // Get 1d slice
    TH1D * h1DKS = h2DKS->ProjectionY(Form("h1D_KS_%d",islice),ptDownBin,ptUpBin);
    Double_t binWidth = h1DKS->GetXaxis()->GetBinWidth(50);
    h1DKS->Scale(1./binWidth);
    
    // Get Clean Samples 1d slice
    h1CleanSamples[0] = GetClean1DSlice(hCleanSamples[0],"h1Electron_slice_%f" ,ptDownBin,ptUpBin);
    h1CleanSamples[1] = GetClean1DSlice(hCleanSamples[1],"h1Pion_slice_%f"     ,ptDownBin,ptUpBin);
    h1CleanSamples[2] = GetClean1DSlice(hCleanSamples[2],"h1Kaon_slice_%f"     ,ptDownBin,ptUpBin);
    h1CleanSamples[3] = GetClean1DSlice(hCleanSamples[3],"h1Proton_slice_%f"   ,ptDownBin,ptUpBin);
    h1CleanSamples[4] = GetClean1DSlice(hCleanSamples[4],"h1PionXTOF_slice_%f" ,ptDownBin,ptUpBin);
    
    // Get Clean MC Samples 1d slice
    h1MCSamples[0]    = GetClean1DSlice(hMCSamples[0]   ,"h1ElectronMC_slice_%f",ptDownBin,ptUpBin);
    h1MCSamples[1]    = GetClean1DSlice(hMCSamples[1]   ,"h1PionMC_slice_%f"    ,ptDownBin,ptUpBin);
    h1MCSamples[2]    = GetClean1DSlice(hMCSamples[2]   ,"h1KaonMC_slice_%f"    ,ptDownBin,ptUpBin);
    h1MCSamples[3]    = GetClean1DSlice(hMCSamples[3]   ,"h1ProtonMC_slice_%f"  ,ptDownBin,ptUpBin);
    h1MCSamples[5]    = GetClean1DSlice(hMCSamples[5]   ,"h1MuonMC_slice_%f"    ,ptDownBin,ptUpBin);
    
    // Get Expected 1d Slice
    h1ExpectedMean[0]   = GetClean1DSlice(hExpectedMean[0],"h1ElectronExMean_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedMean[1]   = GetClean1DSlice(hExpectedMean[1],"h1PionExMean_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedMean[2]   = GetClean1DSlice(hExpectedMean[2],"h1KaonExMean_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedMean[3]   = GetClean1DSlice(hExpectedMean[3],"h1ProtonExMean_slice_%f"  ,ptDownBin,ptUpBin);
    h1ExpectedMean[4]   = GetClean1DSlice(hExpectedMean[4],"h1DeuteronExMean_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedMean[5]   = GetClean1DSlice(hExpectedMean[5],"h1MuonExMean_slice_%f"    ,ptDownBin,ptUpBin);
    
    h1ExpectedSigma[0]   = GetClean1DSlice(hExpectedSigma[0],"h1ElectronExSigma_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedSigma[1]   = GetClean1DSlice(hExpectedSigma[1],"h1PionExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedSigma[2]   = GetClean1DSlice(hExpectedSigma[2],"h1KaonExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedSigma[3]   = GetClean1DSlice(hExpectedSigma[3],"h1ProtonExSigma_slice_%f"  ,ptDownBin,ptUpBin);
    h1ExpectedSigma[4]   = GetClean1DSlice(hExpectedSigma[4],"h1DeuteronExSigma_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedSigma[5]   = GetClean1DSlice(hExpectedSigma[5],"h1MuonExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    
    // Get Expected Mean and Sigma 
    arrMean[0] = h1ExpectedMean[0]->GetMean();
    arrMean[1] = h1ExpectedMean[1]->GetMean();
    arrMean[2] = h1ExpectedMean[2]->GetMean();
    arrMean[3] = h1ExpectedMean[3]->GetMean();
    arrMean[4] = h1ExpectedMean[4]->GetMean();
    arrMean[5] = h1ExpectedMean[5]->GetMean();
    arrSigma[0] = h1ExpectedSigma[0]->GetMean();
    arrSigma[1] = h1ExpectedSigma[1]->GetMean();
    arrSigma[2] = h1ExpectedSigma[2]->GetMean();
    arrSigma[3] = h1ExpectedSigma[3]->GetMean();
    arrSigma[4] = h1ExpectedSigma[4]->GetMean();
    arrSigma[5] = h1ExpectedSigma[5]->GetMean();
    if (arrMean[4]>0 && arrSigma[4]<3 && iks<nbins) arrSigma[4]=20.;
    if (arrMean[3]>0 && arrSigma[3]<3 && iks<nbins) arrSigma[3]=15.;
    if (arrMean[4]>0 && arrSigma[4]<3 && iks<nbins) arrSigma[4]=3.8;
    if (arrMean[3]>0 && arrSigma[3]<3 && iks<nbins) arrSigma[3]=3.;
    
    //  Make the fits        
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[0] , elParClean, arrMean, arrSigma, kElectron, 0);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[1] , piParClean, arrMean, arrSigma, kPion    , 0);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[2] , kaParClean, arrMean, arrSigma, kKaon    , 0);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[3] , prParClean, arrMean, arrSigma, kProton  , 0);
    
    FitParticleSampleForKSEstimate(iks, h1MCSamples[0]    , elParMC   , arrMean, arrSigma, kElectron, 1);
    FitParticleSampleForKSEstimate(iks, h1MCSamples[1]    , piParMC   , arrMean, arrSigma, kPion    , 1);
    FitParticleSampleForKSEstimate(iks, h1MCSamples[2]    , kaParMC   , arrMean, arrSigma, kKaon    , 1);
    FitParticleSampleForKSEstimate(iks, h1MCSamples[3]    , prParMC   , arrMean, arrSigma, kProton  , 1);
    FitParticleSampleForKSEstimate(iks, h1MCSamples[5]    , muParMC   , arrMean, arrSigma, kMuon    , 1);   
     
    Fit1DSliceForKSEstimate(iks, h1DKS, pars, arrMean, arrSigma);
    
    cout.setf(ios::fixed);
    cout << " ******************************************************************************************* " << endl;
    cout << " iks = "            << setprecision(2) << iks            << setw(12);
    cout << " ptDownBin = "      << setprecision(2) << ptDownBin      << setw(12);
    cout << " ptUpBin = "        << setprecision(2) << ptUpBin        << setw(12);
    cout << " islice = "         << setprecision(2) << islice         << setw(12); 
    cout << " p = "              << setprecision(3) << pt             << setw(12);
    cout << " chi2 = "           << setprecision(2) << pars[30]       << setw(12)   << endl;
    cout << " =================================================================== " << endl;
    cout << " *** piKClean = "   << setprecision(2) << piParClean[3]  << setw(12);
    cout << " kaKClean = "       << setprecision(2) << kaParClean[3]  << setw(12);
    cout << " prKClean = "       << setprecision(2) << prParClean[3]  << setw(12);
    
    cout << " *** piKMC = "      << setprecision(2) << piParMC[3]     << setw(12);
    cout << " kaKMC = "          << setprecision(2) << kaParMC[3]     << setw(12);
    cout << " prKMC = "          << setprecision(2) << prParMC[3]     << setw(12);
    
    cout << " *** piKFit = "     << setprecision(2) << pars[8]        << setw(12);
    cout << " kaKFit = "         << setprecision(2) << pars[13]       << setw(12);
    cout << " prKFit = "         << setprecision(2) << pars[18]       << setw(12)   << endl;
    
    cout << " *** piSClean = "   << setprecision(2) << piParClean[4]  << setw(12);
    cout << " kaSClean = "       << setprecision(2) << kaParClean[4]  << setw(12);
    cout << " prSClean = "       << setprecision(2) << prParClean[4]  << setw(12);
    
    cout << " *** piSMC = "      << setprecision(2) << piParMC[4]     << setw(12);
    cout << " kaSMC = "          << setprecision(2) << kaParMC[4]     << setw(12);
    cout << " prSMC = "          << setprecision(2) << prParMC[4]     << setw(12);
    
    cout << " *** piSFit = "     << setprecision(2) << pars[9]        << setw(12);
    cout << " kaSFit = "         << setprecision(2) << pars[14]       << setw(12);
    cout << " prSFit = "         << setprecision(2) << pars[19]       << setw(12) << endl;
    cout << " ******************************************************************************************* " << endl;
    
    if (iks==0){
      
    // Safe region Kurtosis 
      if (pt>=0.6 && pt<=0.8)   hCleanKurtosis[kPion]  .Fill(piParClean[3]);
      if (pt>=0.8 && pt<=1.0)   hCleanKurtosis[kKaon]  .Fill(kaParClean[3]);
      if (pt>=1.1 && pt<=1.4)   hCleanKurtosis[kProton].Fill(prParClean[3]);
      
      if (pt>=0.4 && pt<=0.6)   hMCKurtosis[kPion]     .Fill(piParMC[3]);
      if (pt>=0.6 && pt<=0.8)   hMCKurtosis[kKaon]     .Fill(kaParMC[3]);
      if (pt>=0.8 && pt<=1. )   hMCKurtosis[kProton]   .Fill(prParMC[3]);
      
      if (pt>=0.4  && pt<=0.6)  hFitKurtosis[kPion]    .Fill(pars[8]);
      if (pt>=0.41 && pt<=0.46) hFitKurtosis[kKaon]    .Fill(pars[13]);
      if (pt>=0.7  && pt<=0.8)  hFitKurtosis[kProton]  .Fill(pars[18]);
      
    } else {
    
    // Safe region Skewness
      if (pt>=0.6 && pt<=0.8)   hCleanSkewness[kPion]  .Fill(piParClean[4]);
      if (pt>=0.8 && pt<=1.0)   hCleanSkewness[kKaon]  .Fill(kaParClean[4]);
      if (pt>=1.1 && pt<=1.4)   hCleanSkewness[kProton].Fill(prParClean[4]);
      
      if (pt>=0.4 && pt<=0.6)   hMCSkewness[kPion]     .Fill(piParMC[4]);
      if (pt>=0.6 && pt<=0.8)   hMCSkewness[kKaon]     .Fill(kaParMC[4]);
      if (pt>=1.0 && pt<=1.4)   hMCSkewness[kProton]   .Fill(prParMC[4]);
      
      if (pt>=0.4  && pt<=0.6)  hFitSkewness[kPion]    .Fill(pars[9]);
      if (pt>=0.41 && pt<=0.46) hFitSkewness[kKaon]    .Fill(pars[14]);
      if (pt>=0.7  && pt<=0.8)  hFitSkewness[kProton]  .Fill(pars[19]);

    }
        
    // dump all info to tree
    debugFile -> GetFile()->cd();
    *debugFile << "ks" <<
        
        "iks="          << iks                <<
        "p="            << pt                 <<
        "eta="          << etaBin             <<
        "cent="         << centBin            <<
        
        "elMean="       << arrMean[0]         <<
        "piMean="       << arrMean[1]         <<
        "kaMean="       << arrMean[2]         <<
        "prMean="       << arrMean[3]         <<
        "deMean="       << arrMean[4]         <<
        "muMean="       << arrMean[5]         <<
        
        "elSigma="      << arrSigma[0]        <<
        "piSigma="      << arrSigma[1]        <<
        "kaSigma="      << arrSigma[2]        <<
        "prSigma="      << arrSigma[3]        <<
        "deSigma="      << arrSigma[4]        <<
        "muSigma="      << arrSigma[5]        ;
     
    *debugFile << "ks" <<
        
        "elCleanAmp="      << elParClean[0]    << 
        "elCleanMean="     << elParClean[1]    << 
        "elCleanSigma="    << elParClean[2]    << 
        "elCleanKurtosis=" << elParClean[3]    << 
        "elCleanSkew="     << elParClean[4]    << 
        "elCleanInt="      << elParClean[5]    << 
        "elCleanChi2="     << elParClean[6]    << 
        
        "piCleanAmp="      << piParClean[0]    << 
        "piCleanMean="     << piParClean[1]    << 
        "piCleanSigma="    << piParClean[2]    << 
        "piCleanKurtosis=" << piParClean[3]    << 
        "piCleanSkew="     << piParClean[4]    << 
        "piCleanInt="      << piParClean[5]    << 
        "piCleanChi2="     << piParClean[6]    << 
            
        "kaCleanAmp="      << kaParClean[0]    << 
        "kaCleanMean="     << kaParClean[1]    << 
        "kaCleanSigma="    << kaParClean[2]    << 
        "kaCleanKurtosis=" << kaParClean[3]    << 
        "kaCleanSkew="     << kaParClean[4]    << 
        "kaCleanInt="      << kaParClean[5]    << 
        "kaCleanChi2="     << kaParClean[6]    << 

        "prCleanAmp="      << prParClean[0]    << 
        "prCleanMean="     << prParClean[1]    << 
        "prCleanSigma="    << prParClean[2]    << 
        "prCleanKurtosis=" << prParClean[3]    << 
        "prCleanSkew="     << prParClean[4]    << 
        "prCleanInt="      << prParClean[5]    << 
        "prCleanChi2="     << prParClean[6]    ;
     
    *debugFile << "ks" <<
        
        "elMCAmp="         << elParMC[0]       << 
        "elMCMean="        << elParMC[1]       << 
        "elMCSigma="       << elParMC[2]       << 
        "elMCKurtosis="    << elParMC[3]       << 
        "elMCSkew="        << elParMC[4]       << 
        "elMCInt="         << elParMC[5]       << 
        "elMCChi2="        << elParMC[6]       << 
        
        "piMCAmp="         << piParMC[0]       << 
        "piMCMean="        << piParMC[1]       << 
        "piMCSigma="       << piParMC[2]       << 
        "piMCKurtosis="    << piParMC[3]       << 
        "piMCSkew="        << piParMC[4]       << 
        "piMCInt="         << piParMC[5]       << 
        "piMCChi2="        << piParMC[6]       << 
            
        "kaMCAmp="         << kaParMC[0]       << 
        "kaMCMean="        << kaParMC[1]       << 
        "kaMCSigma="       << kaParMC[2]       << 
        "kaMCKurtosis="    << kaParMC[3]       << 
        "kaMCSkew="        << kaParMC[4]       << 
        "kaMCInt="         << kaParMC[5]       << 
        "kaMCChi2="        << kaParMC[6]       << 

        "prMCAmp="         << prParMC[0]       << 
        "prMCMean="        << prParMC[1]       << 
        "prMCSigma="       << prParMC[2]       << 
        "prMCKurtosis="    << prParMC[3]       << 
        "prMCSkew="        << prParMC[4]       << 
        "prMCInt="         << prParMC[5]       << 
        "prMCChi2="        << prParMC[6]       ;
     
    *debugFile << "ks" <<
        
        "elFitAmp="        << pars[0]       << 
        "elFitMean="       << pars[1]       << 
        "elFitSigma="      << pars[2]       << 
        "elFitKurtosis="   << pars[3]       << 
        "elFitSkew="       << pars[4]       << 
        
        "piFitAmp="        << pars[5]       << 
        "piFitMean="       << pars[6]       << 
        "piFitSigma="      << pars[7]       << 
        "piFitKurtosis="   << pars[8]       << 
        "piFitSkew="       << pars[9]       << 
        
        "kaFitAmp="        << pars[10]       << 
        "kaFitMean="       << pars[11]       << 
        "kaFitSigma="      << pars[12]       << 
        "kaFitKurtosis="   << pars[13]       << 
        "kaFitSkew="       << pars[14]       << 
        
        "prFitAmp="        << pars[15]       << 
        "prFitMean="       << pars[16]       << 
        "prFitSigma="      << pars[17]       << 
        "prFitKurtosis="   << pars[18]       << 
        "prFitSkew="       << pars[19]       << 
        
        "deFitAmp="        << pars[20]       << 
        "deFitMean="       << pars[21]       << 
        "deFitSigma="      << pars[22]       << 
        "deFitKurtosis="   << pars[23]       << 
        "deFitSkew="       << pars[24]       << 
        
        "muFitAmp="        << pars[25]       << 
        "muFitMean="       << pars[26]       << 
        "muFitSigma="      << pars[27]       << 
        "muFitKurtosis="   << pars[28]       << 
        "muFitSkew="       << pars[29]       << 
        
        "totChi2="         << pars[30]       << 

        "\n";
    
    
    delete [] arrSigma;
    delete [] arrMean;
    delete [] arrCleanSigma;
    delete [] arrCleanMean;
    delete [] piParClean;
    delete [] elParClean;
    delete [] kaParClean;
    delete [] prParClean;    
    delete [] piParMC;
    delete [] elParMC;
    delete [] kaParMC;
    delete [] prParMC;
    delete [] pars;
        
  }
  
  // Set the kurtosis and skewness to be used in further iterations
  if (iks==0){
   
    kurtosisFixClean[kPion]   = hCleanKurtosis[kPion]  .GetMean();
    kurtosisFixClean[kKaon]   = hCleanKurtosis[kKaon]  .GetMean();
    kurtosisFixClean[kProton] = hCleanKurtosis[kProton].GetMean();
    
    kurtosisFixMC[kPion]      = hMCKurtosis[kPion]  .GetMean();
    kurtosisFixMC[kKaon]      = hMCKurtosis[kKaon]  .GetMean();
    kurtosisFixMC[kProton]    = hMCKurtosis[kProton].GetMean();
    
    kurtosisFixFit[kPion]     = hFitKurtosis[kPion]  .GetMean();
    kurtosisFixFit[kKaon]     = hFitKurtosis[kKaon]  .GetMean();
    kurtosisFixFit[kProton]   = hFitKurtosis[kProton].GetMean();
  } 
   
  // dump final fixed KS into debug file
  debugFile -> GetFile()->cd();
  *debugFile << "fixedKS" <<
      
      "iks="             << iks                         <<
      "eta="             << etaBin                      <<
      "cent="            << centBin                     <<
  
      "FitSkewPi="       << skewnessFixFit[kPion]       << 
      "FitSkewKa="       << skewnessFixFit[kKaon]       << 
      "FitSkewPr="       << skewnessFixFit[kProton]     << 
        
      "MCSkewPi="        << skewnessFixMC[kPion]        << 
      "MCSkewKa="        << skewnessFixMC[kKaon]        << 
      "MCSkewPr="        << skewnessFixMC[kProton]      << 
      
      "CleanSkewPi="     << skewnessFixClean[kPion]     << 
      "CleanSkewKa="     << skewnessFixClean[kKaon]     << 
      "CleanSkewPr="     << skewnessFixClean[kProton]   << 
      
      "FitKurtPi="       << kurtosisFixFit[kPion]       << 
      "FitKurtKa="       << kurtosisFixFit[kKaon]       << 
      "FitKurtPr="       << kurtosisFixFit[kProton]     << 
      
      "MCKurtPi="        << kurtosisFixMC[kPion]       << 
      "MCKurtKa="        << kurtosisFixMC[kKaon]       << 
      "MCKurtPr="        << kurtosisFixMC[kProton]     << 
      
      "CleanKurtPi="     << kurtosisFixClean[kPion]       << 
      "CleanKurtKa="     << kurtosisFixClean[kKaon]       << 
      "CleanKurtPr="     << kurtosisFixClean[kProton]     << 
        
      "\n";
  
  debugFile -> GetFile()->cd();
  for (Int_t i = 1; i<4; i++){
    
    if (iks==0) hFitKurtosis[i]  .Write(Form("Fit_Kurtosis_%d"  ,i));
    if (iks==1) hFitSkewness[i]  .Write(Form("Fit_Skewness_%d"  ,i));

    if (iks==0) hCleanKurtosis[i].Write(Form("Clean_Kurtosis_%d",i));
    if (iks==1) hCleanSkewness[i].Write(Form("Clean_Skewness_%d",i));

    if (iks==0) hMCKurtosis[i]   .Write(Form("MC_Kurtosis_%d"   ,i));
    if (iks==1) hMCSkewness[i]   .Write(Form("MC_Skewness_%d"   ,i));

  }
  
  if (iks==1) {
    
    if (automaticKS) {
      skewnessFixFit[kPion]     = hFitSkewness[kPion]  .GetMean();
      skewnessFixFit[kKaon]     = hFitSkewness[kKaon]  .GetMean();
      skewnessFixFit[kProton]   = hFitSkewness[kProton].GetMean();
        
      skewnessFixClean[kPion]   = hCleanSkewness[kPion]  .GetMean();
      skewnessFixClean[kKaon]   = hCleanSkewness[kKaon]  .GetMean();
      skewnessFixClean[kProton] = hCleanSkewness[kProton].GetMean();
    
      skewnessFixMC[kPion]      = hMCSkewness[kPion]  .GetMean();
      skewnessFixMC[kKaon]      = hMCSkewness[kKaon]  .GetMean();
      skewnessFixMC[kProton]    = hMCSkewness[kProton].GetMean();
    } else {
      skewnessFixFit[kPion]     = TMath::Min(hFitSkewness[kPion]  .GetMean(),1.);
      skewnessFixFit[kKaon]     = TMath::Min(hFitSkewness[kKaon]  .GetMean(),0.7);
      skewnessFixFit[kProton]   = TMath::Min(hFitSkewness[kProton].GetMean(),0.6);
              
      skewnessFixClean[kPion]   = TMath::Min(hCleanSkewness[kPion]  .GetMean(),1.);
      skewnessFixClean[kKaon]   = TMath::Min(hCleanSkewness[kKaon]  .GetMean(),0.7);
      skewnessFixClean[kProton] = TMath::Min(hCleanSkewness[kProton].GetMean(),0.6);
      
      skewnessFixMC[kPion]      = TMath::Min(hMCSkewness[kPion]  .GetMean(),1.);
      skewnessFixMC[kKaon]      = TMath::Min(hMCSkewness[kKaon]  .GetMean(),0.7);
      skewnessFixMC[kProton]    = TMath::Min(hMCSkewness[kProton].GetMean(),0.6);
    }
  }
  
  // Set the final kurtosis and skewness params
  if (iks==1){
    if (MC) {
      skewnessFix[kPion]     = skewnessFixMC[kPion];
      skewnessFix[kKaon]     = skewnessFixMC[kKaon];
      skewnessFix[kProton]   = skewnessFixMC[kProton];
      kurtosisFix[kPion]     = kurtosisFixMC[kPion];
      kurtosisFix[kKaon]     = kurtosisFixMC[kKaon];
      kurtosisFix[kProton]   = kurtosisFixMC[kProton];
    } else {
      skewnessFix[kPion]     = skewnessFixFit[kPion];
      skewnessFix[kKaon]     = skewnessFixFit[kKaon];
      skewnessFix[kProton]   = skewnessFixFit[kProton];
      kurtosisFix[kPion]     = kurtosisFixFit[kPion];
      kurtosisFix[kKaon]     = kurtosisFixFit[kKaon];
      kurtosisFix[kProton]   = kurtosisFixFit[kProton];
    }
  }
 
}
// -------------------------------------------------------------------------------------------------------
void FitAllSamples(const Int_t nSlice, const Double_t ptMin, const Double_t ptMax, TH2D *h2D, TH2D **hCleanSamples, TH2D **hMCSamples, TH2D **hExpectedMean, TH2D **hExpectedSigma, TString sampleFileName){
  
  //
  // Fit only samples
  // 
  
  TMinuit g;
  g.SetMaxIterations(10000);
  gMinuit->SetMaxIterations(10000);
  sampleFile   = new TTreeSRedirector(sampleFileName,"recreate");
  cleanResArr  . SetOwner(kTRUE);
  MCResArr     . SetOwner(kTRUE);
  
  
  Int_t ptDownBin = 0;
  Int_t ptUpBin   = 0; 
  Double_t ptStep = 0.;
  Double_t ptDown = 0.;
  Double_t ptUp   = 0.;
  Double_t pt     = 0.;
  Int_t ihslice = 0;
  for (Int_t islice = 0; islice<nSlice+1000; islice++){
    
    if (islice%10 == 0) cout << " slice = " << islice << endl;
    
    if (islice<nSlice){
      ptStep = ((ptMax-ptMin)/(Double_t)nSlice);
      ptDown = ptMin+islice*ptStep;
      ptUp   = ptDown+ptStep;
      pt     = (ptDown+ptUp)/2.;
    } else {
      ptStep = highPbinWidth;
      ptDown = ptMax+ihslice*ptStep;
      ptUp   = ptDown+ptStep;
      pt     = (ptDown+ptUp)/2.;
      ihslice++;
    }
    
    // for small tests cut the loop short
    if (lowPfit   && islice==nSlice) break;
    if (nSlice<20 && islice==nSlice) break;
    
    // Fit procedure should stop at "analysisRange" which is maximum p value to be chosen
    if (ptUp>analysisRange) break;
    
    // find bins to be used in the projection
    if (islice<nSlice){
      ptDownBin = (h2D->GetXaxis()->FindBin(ptMin)-1)+2*islice+1;                                // TO FIX
      ptUpBin   = (h2D->GetXaxis()->FindBin(ptMin)-1)+2*islice+2;                                // TO FIX
    } else {
      ptDownBin = h2D->GetXaxis()->FindBin(ptDown);
      ptUpBin   = h2D->GetXaxis()->FindBin(ptUp)-1;
    }
        
    // Storage for the histograms
    TObjArray * cleanSampArr = new TObjArray(6);  cleanSampArr -> SetOwner(kTRUE);
    TObjArray * MCSampArr    = new TObjArray(6);  MCSampArr    -> SetOwner(kTRUE); 
    cleanSampArr->SetName(Form("cleanSamples_%d_%4.3f",islice,pt));
    MCSampArr   ->SetName(Form("MCSamples_%d_%4.3f",islice,pt));
    
    TH1D * h1CleanSamples[6];
    TH1D * h1MCSamples[6];
    TH1D * h1ExpectedMean[6];
    TH1D * h1ExpectedSigma[6];
         
     // Real particle parameters
    Double_t * arrSigma      = new Double_t[10];
    Double_t * arrMean       = new Double_t[10];
    Double_t * arrCleanSigma = new Double_t[10];
    Double_t * arrCleanMean  = new Double_t[10];
    
    // Clean Particle parameters
    Double_t * elParClean    = new Double_t[10];
    Double_t * piParClean    = new Double_t[10];
    Double_t * kaParClean    = new Double_t[10];
    Double_t * prParClean    = new Double_t[10];
    Double_t * piParCleanTOF = new Double_t[10];
    // MC Particle parameters
    Double_t * elParMC       = new Double_t[10];
    Double_t * piParMC       = new Double_t[10];
    Double_t * kaParMC       = new Double_t[10];
    Double_t * prParMC       = new Double_t[10];
    Double_t * muParMC       = new Double_t[10];
        
    for (Int_t ipar=0; ipar<10; ipar++)
    {
      arrMean[ipar]       = 0;
      arrSigma[ipar]      = 0;
      arrCleanMean[ipar]  = 0;
      arrCleanSigma[ipar] = 0;
      piParClean[ipar]    = 0;
      kaParClean[ipar]    = 0;
      prParClean[ipar]    = 0;
      elParClean[ipar]    = 0;
      piParCleanTOF[ipar] = 0;
      piParMC[ipar]       = 0;
      kaParMC[ipar]       = 0;
      prParMC[ipar]       = 0;
      elParMC[ipar]       = 0;
      muParMC[ipar]       = 0;
    }
    
    // Get Clean Samples 1d slice
    h1CleanSamples[0] = GetClean1DSlice(hCleanSamples[0],"h1Electron_slice_%f" ,ptDownBin,ptUpBin);
    h1CleanSamples[1] = GetClean1DSlice(hCleanSamples[1],"h1Pion_slice_%f"     ,ptDownBin,ptUpBin);
    h1CleanSamples[2] = GetClean1DSlice(hCleanSamples[2],"h1Kaon_slice_%f"     ,ptDownBin,ptUpBin);
    h1CleanSamples[3] = GetClean1DSlice(hCleanSamples[3],"h1Proton_slice_%f"   ,ptDownBin,ptUpBin);
    h1CleanSamples[4] = GetClean1DSlice(hCleanSamples[4],"h1PionXTOF_slice_%f" ,ptDownBin,ptUpBin);
    
    // Get Clean MC Samples 1d slice
    h1MCSamples[0]    = GetClean1DSlice(hMCSamples[0]   ,"h1ElectronMC_slice_%f",ptDownBin,ptUpBin);
    h1MCSamples[1]    = GetClean1DSlice(hMCSamples[1]   ,"h1PionMC_slice_%f"    ,ptDownBin,ptUpBin);
    h1MCSamples[2]    = GetClean1DSlice(hMCSamples[2]   ,"h1KaonMC_slice_%f"    ,ptDownBin,ptUpBin);
    h1MCSamples[3]    = GetClean1DSlice(hMCSamples[3]   ,"h1ProtonMC_slice_%f"  ,ptDownBin,ptUpBin);
    h1MCSamples[5]    = GetClean1DSlice(hMCSamples[5]   ,"h1MuonMC_slice_%f"    ,ptDownBin,ptUpBin);
    
    // Get Expected 1d Slice
    h1ExpectedMean[0]   = GetClean1DSlice(hExpectedMean[0],"h1ElectronExMean_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedMean[1]   = GetClean1DSlice(hExpectedMean[1],"h1PionExMean_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedMean[2]   = GetClean1DSlice(hExpectedMean[2],"h1KaonExMean_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedMean[3]   = GetClean1DSlice(hExpectedMean[3],"h1ProtonExMean_slice_%f"  ,ptDownBin,ptUpBin);
    h1ExpectedMean[4]   = GetClean1DSlice(hExpectedMean[4],"h1DeuteronExMean_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedMean[5]   = GetClean1DSlice(hExpectedMean[5],"h1MuonExMean_slice_%f"    ,ptDownBin,ptUpBin);
    
    h1ExpectedSigma[0]   = GetClean1DSlice(hExpectedSigma[0],"h1ElectronExSigma_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedSigma[1]   = GetClean1DSlice(hExpectedSigma[1],"h1PionExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedSigma[2]   = GetClean1DSlice(hExpectedSigma[2],"h1KaonExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedSigma[3]   = GetClean1DSlice(hExpectedSigma[3],"h1ProtonExSigma_slice_%f"  ,ptDownBin,ptUpBin);
    h1ExpectedSigma[4]   = GetClean1DSlice(hExpectedSigma[4],"h1DeuteronExSigma_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedSigma[5]   = GetClean1DSlice(hExpectedSigma[5],"h1MuonExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    
    // Get Expected Mean and Sigma 
    arrMean[0] = h1ExpectedMean[0]->GetMean();
    arrMean[1] = h1ExpectedMean[1]->GetMean();
    arrMean[2] = h1ExpectedMean[2]->GetMean();
    arrMean[3] = h1ExpectedMean[3]->GetMean();
    arrMean[4] = h1ExpectedMean[4]->GetMean();
    arrMean[5] = h1ExpectedMean[5]->GetMean();
    arrSigma[0] = h1ExpectedSigma[0]->GetMean();
    arrSigma[1] = h1ExpectedSigma[1]->GetMean();
    arrSigma[2] = h1ExpectedSigma[2]->GetMean();
    arrSigma[3] = h1ExpectedSigma[3]->GetMean();
    arrSigma[4] = h1ExpectedSigma[4]->GetMean();
    arrSigma[5] = h1ExpectedSigma[5]->GetMean();
    if (arrMean[4]>0 && arrSigma[4]<3 && islice<nSlice/2.) arrSigma[4]=20.;
    if (arrMean[3]>0 && arrSigma[3]<3 && islice<nSlice/2.) arrSigma[3]=15.;
    if (arrMean[4]>0 && arrSigma[4]<3 && islice>nSlice/2.) arrSigma[4]=3.8;
    if (arrMean[3]>0 && arrSigma[3]<3 && islice>nSlice/2.) arrSigma[3]=3.;
    
    //  Make the fits        
    TMatrixD *elCleanMatrix = FitParticleSample(h1CleanSamples[0] , elParClean, arrMean, arrSigma, kElectron, 0);
    TMatrixD *piCleanMatrix = FitParticleSample(h1CleanSamples[1] , piParClean, arrMean, arrSigma, kPion    , 0);
    TMatrixD *kaCleanMatrix = FitParticleSample(h1CleanSamples[2] , kaParClean, arrMean, arrSigma, kKaon    , 0);
    TMatrixD *prCleanMatrix = FitParticleSample(h1CleanSamples[3] , prParClean, arrMean, arrSigma, kProton  , 0);
    TMatrixD *piCleanMatrixTOF = FitParticleSample(h1CleanSamples[4] , piParCleanTOF, arrMean, arrSigma, kPion  , 0);
    
    TMatrixD *elMCMatrix    = FitParticleSample(h1MCSamples[0]    , elParMC   , arrMean, arrSigma, kElectron, 1);
    TMatrixD *piMCMatrix    = FitParticleSample(h1MCSamples[1]    , piParMC   , arrMean, arrSigma, kPion    , 1);
    TMatrixD *kaMCMatrix    = FitParticleSample(h1MCSamples[2]    , kaParMC   , arrMean, arrSigma, kKaon    , 1);
    TMatrixD *prMCMatrix    = FitParticleSample(h1MCSamples[3]    , prParMC   , arrMean, arrSigma, kProton  , 1);
    TMatrixD *muMCMatrix    = FitParticleSample(h1MCSamples[5]    , muParMC   , arrMean, arrSigma, kMuon    , 1);   
     
    GetCleanExParams(h1CleanSamples[0],arrSigma,arrCleanSigma,arrCleanMean,kElectron);
    GetCleanExParams(h1CleanSamples[1],arrSigma,arrCleanSigma,arrCleanMean,kPion);
    GetCleanExParams(h1CleanSamples[2],arrSigma,arrCleanSigma,arrCleanMean,kKaon);
    GetCleanExParams(h1CleanSamples[3],arrSigma,arrCleanSigma,arrCleanMean,kProton);

    // Dump histograms into arrays
    h1CleanSamples[0] -> SetName(Form("CleanElectron_%d_%4.3f" ,islice,pt));
    h1CleanSamples[1] -> SetName(Form("CleanPion_%d_%4.3f"     ,islice,pt));
    h1CleanSamples[2] -> SetName(Form("CleanKaon_%d_%4.3f"     ,islice,pt));
    h1CleanSamples[3] -> SetName(Form("CleanProton_%d_%4.3f"   ,islice,pt));
    h1CleanSamples[4] -> SetName(Form("CleanPionTOF_%d_%4.3f"   ,islice,pt));
    h1MCSamples[0]    -> SetName(Form("MC_Electron_%d_%4.3f" ,islice,pt));
    h1MCSamples[1]    -> SetName(Form("MC_Pion_%d_%4.3f"     ,islice,pt));
    h1MCSamples[2]    -> SetName(Form("MC_Kaon_%d_%4.3f"     ,islice,pt));
    h1MCSamples[3]    -> SetName(Form("MC_Proton_%d_%4.3f"   ,islice,pt));
    h1MCSamples[5]    -> SetName(Form("MC_Muon_%d_%4.3f"     ,islice,pt));
            
    cleanSampArr -> AddAt(h1CleanSamples[0],0);
    cleanSampArr -> AddAt(h1CleanSamples[1],1);
    cleanSampArr -> AddAt(h1CleanSamples[2],2);
    cleanSampArr -> AddAt(h1CleanSamples[3],3);
    cleanSampArr -> AddAt(h1CleanSamples[4],4);
    MCSampArr    -> AddAt(h1MCSamples[0],0);
    MCSampArr    -> AddAt(h1MCSamples[1],1);
    MCSampArr    -> AddAt(h1MCSamples[2],2);
    MCSampArr    -> AddAt(h1MCSamples[3],3);
    MCSampArr    -> AddAt(h1MCSamples[5],5);
    
    cleanResArr  . AddAt(cleanSampArr,islice);
    MCResArr     . AddAt(MCSampArr   ,islice);
        
    // Clean Samples dump tree
    sampleFile->GetFile()->cd();
    *sampleFile << "CleanSamples" <<
        
        "elMean="       << arrMean[0]         <<
        "piMean="       << arrMean[1]         <<
        "kaMean="       << arrMean[2]         <<
        "prMean="       << arrMean[3]         <<
        "deMean="       << arrMean[4]         <<
        "muMean="       << arrMean[5]         <<
        
        "elSigma="      << arrSigma[0]        <<
        "piSigma="      << arrSigma[1]        <<
        "kaSigma="      << arrSigma[2]        <<
        "prSigma="      << arrSigma[3]        <<
        "deSigma="      << arrSigma[4]        <<
        "muSigma="      << arrSigma[5]        ;
    
    *sampleFile << "CleanSamples" <<
        
        "elExCMean="    << arrCleanMean[0]    <<
        "piExCMean="    << arrCleanMean[1]    <<
        "kaExCMean="    << arrCleanMean[2]    <<
        "prExCMean="    << arrCleanMean[3]    <<
        
        "elExCSigma="   << arrCleanSigma[0]   <<
        "piExCSigma="   << arrCleanSigma[1]   <<
        "kaExCSigma="   << arrCleanSigma[2]   <<
        "prExCSigma="   << arrCleanSigma[3]   ;
    
    *sampleFile << "CleanSamples" <<
        
        "elSigma="      << arrSigma[0]        <<
        "piSigma="      << arrSigma[1]        <<
        "kaSigma="      << arrSigma[2]        <<
        "prSigma="      << arrSigma[3]        <<
        "deSigma="      << arrSigma[4]        <<
        "muSigma="      << arrSigma[5]        <<
        
        "slice="        << islice             <<
        "p="            << pt                 <<
        "eta="          << etaBin             <<
        "cent="         << centBin            << 
        
        "elCleanAmp="      << elParClean[0]    << 
        "elCleanMean="     << elParClean[1]    << 
        "elCleanSigma="    << elParClean[2]    << 
        "elCleanKurtosis=" << elParClean[3]    << 
        "elCleanSkew="     << elParClean[4]    << 
        "elCleanInt="      << elParClean[5]    << 
        "elCleanChi2="     << elParClean[6]    << 
        "elCleanMatrix.="  << elCleanMatrix    ;
    
    *sampleFile << "CleanSamples" <<
        
        "piCleanAmp="      << piParClean[0]    << 
        "piCleanMean="     << piParClean[1]    << 
        "piCleanSigma="    << piParClean[2]    << 
        "piCleanKurtosis=" << piParClean[3]    << 
        "piCleanSkew="     << piParClean[4]    << 
        "piCleanInt="      << piParClean[5]    << 
        "piCleanChi2="     << piParClean[6]    << 
        "piCleanMatrix.="  << piCleanMatrix    <<
    
        "piCleanAmpTOF="      << piParCleanTOF[0]    << 
        "piCleanMeanTOF="     << piParCleanTOF[1]    << 
        "piCleanSigmaTOF="    << piParCleanTOF[2]    << 
        "piCleanKurtosisTOF=" << piParCleanTOF[3]    << 
        "piCleanSkewTOF="     << piParCleanTOF[4]    << 
        "piCleanIntTOF="      << piParCleanTOF[5]    << 
        "piCleanChi2TOF="     << piParCleanTOF[6]    << 
        "piCleanMatrixTOF.="  << piCleanMatrixTOF    ;
    
    *sampleFile << "CleanSamples" <<
        
        "kaCleanAmp="      << kaParClean[0]    << 
        "kaCleanMean="     << kaParClean[1]    << 
        "kaCleanSigma="    << kaParClean[2]    << 
        "kaCleanKurtosis=" << kaParClean[3]    << 
        "kaCleanSkew="     << kaParClean[4]    << 
        "kaCleanInt="      << kaParClean[5]    << 
        "kaCleanChi2="     << kaParClean[6]    << 
        "kaCleanMatrix.="  << kaCleanMatrix    <<

        "prCleanAmp="      << prParClean[0]    << 
        "prCleanMean="     << prParClean[1]    << 
        "prCleanSigma="    << prParClean[2]    << 
        "prCleanKurtosis=" << prParClean[3]    << 
        "prCleanSkew="     << prParClean[4]    << 
        "prCleanInt="      << prParClean[5]    << 
        "prCleanChi2="     << prParClean[6]    << 
        "prCleanMatrix.="  << prCleanMatrix    <<
        
        "\n";
    
    // MC samples dump tree
    *sampleFile << "MCSamples" <<
        
        "elMean="       << arrMean[0]         <<
        "piMean="       << arrMean[1]         <<
        "kaMean="       << arrMean[2]         <<
        "prMean="       << arrMean[3]         <<
        "deMean="       << arrMean[4]         <<
        "muMean="       << arrMean[5]         <<
        
        "elSigma="      << arrSigma[0]        <<
        "piSigma="      << arrSigma[1]        <<
        "kaSigma="      << arrSigma[2]        <<
        "prSigma="      << arrSigma[3]        <<
        "deSigma="      << arrSigma[4]        <<
        "muSigma="      << arrSigma[5]     ;
    
    *sampleFile << "MCSamples" <<
        
        "slice="       << islice         <<
        "p="           << pt             <<
        "eta="         << etaBin         <<
        "cent="        << centBin        << 
        
        "elMCAmp="      << elParMC[0]    << 
        "elMCMean="     << elParMC[1]    << 
        "elMCSigma="    << elParMC[2]    << 
        "elMCKurtosis=" << elParMC[3]    << 
        "elMCSkew="     << elParMC[4]    << 
        "elMCInt="      << elParMC[5]    << 
        "elMCChi2="     << elParMC[6]    << 
        "elMCMatrix.="  << elMCMatrix    <<
        
        "piMCAmp="      << piParMC[0]    << 
        "piMCMean="     << piParMC[1]    << 
        "piMCSigma="    << piParMC[2]    << 
        "piMCKurtosis=" << piParMC[3]    << 
        "piMCSkew="     << piParMC[4]    << 
        "piMCInt="      << piParMC[5]    << 
        "piMCChi2="     << piParMC[6]    << 
        "piMCMatrix.="  << piMCMatrix    ;
    
    *sampleFile << "MCSamples" <<
        
        "kaMCAmp="      << kaParMC[0]    << 
        "kaMCMean="     << kaParMC[1]    << 
        "kaMCSigma="    << kaParMC[2]    << 
        "kaMCKurtosis=" << kaParMC[3]    << 
        "kaMCSkew="     << kaParMC[4]    << 
        "kaMCInt="      << kaParMC[5]    << 
        "kaMCChi2="     << kaParMC[6]    << 
        "kaMCMatrix.="  << kaMCMatrix    << 
    
        "prMCAmp="      << prParMC[0]    << 
        "prMCMean="     << prParMC[1]    << 
        "prMCSigma="    << prParMC[2]    << 
        "prMCKurtosis=" << prParMC[3]    << 
        "prMCSkew="     << prParMC[4]    << 
        "prMCInt="      << prParMC[5]    << 
        "prMCChi2="     << prParMC[6]    << 
        "prMCMatrix.="  << prMCMatrix    ;
    
    *sampleFile << "MCSamples" <<
        
        "muMCAmp="      << muParMC[0]    << 
        "muMCMean="     << muParMC[1]    << 
        "muMCSigma="    << muParMC[2]    << 
        "muMCKurtosis=" << muParMC[3]    << 
        "muMCSkew="     << muParMC[4]    << 
        "muMCInt="      << muParMC[5]    << 
        "muMCChi2="     << muParMC[6]    << 
        "muMCMatrix.="  << muMCMatrix    <<
        
        "\n";
    
    delete [] arrSigma;
    delete [] arrMean;
    delete [] arrCleanSigma;
    delete [] arrCleanMean;
    delete [] piParClean;
    delete [] elParClean;
    delete [] kaParClean;
    delete [] prParClean;
    delete [] piParCleanTOF;
    delete [] piParMC;
    delete [] elParMC;
    delete [] kaParMC;
    delete [] prParMC;
   
    delete elCleanMatrix; 
    delete piCleanMatrix; 
    delete kaCleanMatrix; 
    delete prCleanMatrix; 
    delete piCleanMatrixTOF; 

    delete elMCMatrix; 
    delete piMCMatrix; 
    delete kaMCMatrix; 
    delete prMCMatrix; 
    delete muMCMatrix;
        
  }
    
  sampleFile -> GetFile()->cd();
  cleanResArr . Write("cleanResArr" ,TObject::kSingleKey);
  MCResArr    . Write("MCResArr"    ,TObject::kSingleKey);
      
  // delete arrays
  cleanResArr  . Delete();
  MCResArr     . Delete();
    
  delete sampleFile;
  
}
// -------------------------------------------------------------------------------------------------------
void IterativeFitting(Int_t iIter, const Int_t nSlice, const Double_t ptMin, const Double_t ptMax, TH2D *h2D, TString fileFit, TString readSmooth, TString readSamp){
    
  //
  // Apply fits without any amplitude restriction
  //
  
  // Open files for reading samples and ttrees from previous results
  TFile * ressFile = 0x0;
  TFile * sampFile = 0x0;
  if (iIter>0) ressFile = TFile::Open(readSmooth);
  sampFile = TFile::Open(readSamp);
  GetHelperGraphs(iIter,ressFile);
  
  // calculate bin number analysed
  Int_t allBins = CalculateNBinsP(nSlice,ptMin,ptMax);
  
  TMinuit g;
  g.SetMaxIterations(10000);
  gMinuit->SetMaxIterations(10000);
  TString histFileName  = Form("Hists_Iteration%d_%d_%3.1f_%3.1f_%2.1f_%d.root",iIter,nSlice,ptMin,ptMax,etaBin,centBin);
  histsFile = new TTreeSRedirector(histFileName,"recreate");
  fitFile   = new TTreeSRedirector(fileFit     ,"recreate");
  
  
  // OutPut Arrays  
  fitResArr    . SetOwner(kTRUE);
  fitResiduals . SetOwner(kTRUE);
    
    // loop over slices
  Double_t elIntegralFit[allBins], elFitAmp[allBins], elFitMean[allBins], elFitSigma[allBins];
  Double_t piIntegralFit[allBins], piFitAmp[allBins], piFitMean[allBins], piFitSigma[allBins];
  Double_t kaIntegralFit[allBins], kaFitAmp[allBins], kaFitMean[allBins], kaFitSigma[allBins];
  Double_t prIntegralFit[allBins], prFitAmp[allBins], prFitMean[allBins], prFitSigma[allBins];
  Double_t deIntegralFit[allBins], deFitAmp[allBins], deFitMean[allBins], deFitSigma[allBins];
  Double_t muIntegralFit[allBins], muFitAmp[allBins], muFitMean[allBins], muFitSigma[allBins];
  
  Double_t elFitKurtosis[allBins], elFitSkew[allBins], elChi2Fit[allBins];
  Double_t piFitKurtosis[allBins], piFitSkew[allBins], piChi2Fit[allBins];
  Double_t kaFitKurtosis[allBins], kaFitSkew[allBins], kaChi2Fit[allBins];
  Double_t prFitKurtosis[allBins], prFitSkew[allBins], prChi2Fit[allBins];
  Double_t deFitKurtosis[allBins], deFitSkew[allBins], deChi2Fit[allBins];
  Double_t muFitKurtosis[allBins], muFitSkew[allBins], muChi2Fit[allBins];  
    
  Int_t ptDownBin = 0;
  Int_t ptUpBin   = 0; 
  Double_t ptStep = 0.;
  Double_t ptDown = 0.;
  Double_t ptUp   = 0.;
  Double_t pt     = 0.;
  Int_t ihslice = 0;
  for (Int_t islice = 0; islice<nSlice+1000; islice++){
    
    
    if (islice<nSlice){
      ptStep = ((ptMax-ptMin)/(Double_t)nSlice);
      ptDown = ptMin+islice*ptStep;
      ptUp   = ptDown+ptStep;
      pt     = (ptDown+ptUp)/2.;
    } else {
      ptStep = highPbinWidth;
      ptDown = ptMax+ihslice*ptStep;
      ptUp   = ptDown+ptStep;
      pt     = (ptDown+ptUp)/2.;
      ihslice++;
    }
    
    // For MC analysis only upto nSlice (e.g. 300)
//     if (islice>=nSlice && MC) break; 
    
    // for small tests cut the loop short
    if (lowPfit   && islice==nSlice) break;
    if (nSlice<=20 && islice==nSlice) break;
    
    // Fit procedure should stop at "analysisRange" which is maximum p value to be chosen
    if (ptUp>analysisRange) break;
    
    // find bins to be used in the projection
    if (islice<nSlice){
      ptDownBin = (h2D->GetXaxis()->FindBin(ptMin)-1)+2*islice+1;                                // TO FIX
      ptUpBin   = (h2D->GetXaxis()->FindBin(ptMin)-1)+2*islice+2;                                // TO FIX
    } else {
      ptDownBin = h2D->GetXaxis()->FindBin(ptDown);
      ptUpBin   = h2D->GetXaxis()->FindBin(ptUp)-1;
    }    
   
    // main 1d slice used in the fitting 
    TH1D * h1D = h2D->ProjectionY(Form("h1D_slice_%d",islice),ptDownBin,ptUpBin);
    Double_t binWidth = h1D->GetXaxis()->GetBinWidth(50);
    h1D->Scale(1./binWidth);
    TH1D * hResidual = (TH1D*)h1D->Clone();

    Double_t * arrSigma      = new Double_t[10];
    Double_t * arrMean       = new Double_t[10];
    Double_t * arrMeanWindow = new Double_t[10];
    Double_t * cleanParMean  = new Double_t[10];
    Double_t * MCPars        = new Double_t[10];
        
    for (Int_t ipar=0; ipar<10; ipar++)
    {
      arrMean[ipar]       = 0;
      arrSigma[ipar]      = 0;
      arrMeanWindow[ipar] = 0;
      cleanParMean[ipar]  = 0;
      MCPars[ipar]        = 0;
    }
    
    GetExandSampleParameters(islice,sampFile,arrMean,arrSigma,cleanParMean,MCPars);
    
    // Fit part
    maxBin = h1D->GetBinContent(h1D->GetMaximumBin());   
    Double_t fitWinMean  = 1.2;
    
    elMeanMin = TMath::Max(30.,(arrMean[0]-2.*arrSigma[0]));                                   // electron
    piMeanMin = TMath::Max(30.,(arrMean[1]-fitWinMean*arrSigma[1]));                           // pion
    kaMeanMin = TMath::Max(30.,(arrMean[2]-fitWinMean*arrSigma[2]));                           // kaon
    prMeanMin = TMath::Max(30.,(arrMean[3]-fitWinMean*arrSigma[3]));                           // proton
    deMeanMin = TMath::Max(30.,(arrMean[4]-fitWinMean*arrSigma[4]));                           // deuteron
    muMeanMin = TMath::Max(30.,(arrMean[5]-fitWinMean*arrSigma[5]));                           // deuteron
        
    elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+2.*arrSigma[0]));
    piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+fitWinMean*arrSigma[1]));
    kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+fitWinMean*arrSigma[2]));
    prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+fitWinMean*arrSigma[3]));
    deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[4]+fitWinMean*arrSigma[4]));
    muMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[5]+fitWinMean*arrSigma[5]));
    
    if (pt<0.3) {
      prMeanMin = TMath::Max(30.              ,(arrMean[3]-arrSigma[3]*3.));                         
      prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+arrSigma[3]*3.));

      deMeanMin = TMath::Max(30.              ,(arrMean[4]-arrSigma[4]*3.));                          
      deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[4]+arrSigma[4]*3.));
    }
    
    // Parameters:
    // 0 --> Amplitude
    // 1 --> Mean
    // 2 --> Sigma
    // 3 --> kurtosis
    // 4 --> skewness
    
    TF1 *g1    = new TF1("g1",fitFunctionGenGaus,elMeanMin,elMeanMax);        g1->SetNpx(1000);
    TF1 *g2    = new TF1("g2",fitFunctionGenGaus,piMeanMin,piMeanMax);        g2->SetNpx(1000);
    TF1 *g3    = new TF1("g3",fitFunctionGenGaus,kaMeanMin,kaMeanMax);        g3->SetNpx(1000);
    TF1 *g4    = new TF1("g4",fitFunctionGenGaus,prMeanMin,prMeanMax);        g4->SetNpx(1000);
    TF1 *g5    = new TF1("g5",fitFunctionGenGaus,deMeanMin,deMeanMax);        g5->SetNpx(1000);
    TF1 *g6    = new TF1("g6",fitFunctionGenGaus,muMeanMin,muMeanMax);        g6->SetNpx(1000);
        
    g1->SetLineColor(4);          // Blue
    g2->SetLineColor(1);          // Black
    g3->SetLineColor(3);          // Green
    g4->SetLineColor(kOrange);    // Orange
    g5->SetLineColor(6);          // Magenta
    g6->SetLineColor(2);          // Red
        
    // restrict the parameters of gaus functions of individual particles
    g1->SetParLimits(0,0.,maxBin);
    g2->SetParLimits(0,0.,maxBin);
    g3->SetParLimits(0,0.,maxBin);
    g4->SetParLimits(0,0.,maxBin);
    g5->SetParLimits(0,0.,maxBin);
    g6->SetParLimits(0,0.,maxBin);
        
    g1->SetParLimits(1,elMeanMin,elMeanMax);
    g2->SetParLimits(1,piMeanMin,piMeanMax);
    g3->SetParLimits(1,kaMeanMin,kaMeanMax);
    g4->SetParLimits(1,prMeanMin,prMeanMax);
    g5->SetParLimits(1,deMeanMin,deMeanMax);
    g6->SetParLimits(1,muMeanMin,muMeanMax);   
    
    g1->SetParLimits(2,arrSigma[0]/5.,arrSigma[0]*2.);
    g2->SetParLimits(2,arrSigma[1]/5.,arrSigma[1]*2.);
    g3->SetParLimits(2,arrSigma[2]/5.,arrSigma[2]*2.);
    g4->SetParLimits(2,arrSigma[3]/5.,arrSigma[3]*2.);
    g5->SetParLimits(2,arrSigma[4]/5.,arrSigma[4]*3.);
    g6->SetParLimits(2,arrSigma[5]/5.,arrSigma[5]*3.);
    
    g1->SetParLimits(3,kurtosisMin,kurtosisMax);
    g2->SetParLimits(3,kurtosisMin,kurtosisMax);
    g3->SetParLimits(3,kurtosisMin,kurtosisMax);
    g4->SetParLimits(3,kurtosisMin,kurtosisMax);
    g5->SetParLimits(3,kurtosisMin,kurtosisMax);
    g6->SetParLimits(3,kurtosisMin,kurtosisMax);

    g1->SetParLimits(4,skewMin,skewMax);
    g2->SetParLimits(4,skewMin,skewMax);
    g3->SetParLimits(4,skewMin,skewMax);
    g4->SetParLimits(4,skewMin,skewMax);
    g5->SetParLimits(4,skewMin,skewMax);
    g6->SetParLimits(4,skewMin,skewMax);
    
    if (fixedK || fixedS) SetFixedKSparameterForIndividualFits(g1,g2,g3,g4,g5,g6);
            
        // Fit h1D for individual particles
    if (arrMean[0]>=10) h1D->Fit(g1,"QNR");
    if (arrMean[1]>=10) h1D->Fit(g2,"QNR+");
    if (arrMean[2]>=10) h1D->Fit(g3,"QNR+");
    if (arrMean[3]>=10) h1D->Fit(g4,"QNR+");
    if (arrMean[4]>=10) h1D->Fit(g5,"QNR+");
    if (arrMean[5]>=10) h1D->Fit(g6,"QNR+");
    
    // Get Fit parameters from the individual fits and set for the total fit
    Double_t par[30] = {0};
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[5]);
    g3->GetParameters(&par[10]);
    g4->GetParameters(&par[15]);
    g5->GetParameters(&par[20]);
    g5->GetParameters(&par[25]);
        
    TF1 *total = new TF1("total","g1+g2+g3+g4+g5+g6",dEdxMin,dEdxMax);
    total->SetNpx(1000);
    total->SetLineWidth(3);
    total->SetLineColor(2);
    total->SetParameters(par);
            
    // Some setters for the total fit 
    SetTotalFitParameters(iIter,islice,nSlice,total,arrMean,arrSigma,arrMeanWindow,cleanParMean,MCPars,h1D,pt); 
    SetParNamesOfTotalFit(total);
    CheckIfMeanIsZero(total,arrMean);
    
        // Apply total fit
//     gStyle->SetOptFit(0);
//     gStyle->SetOptFit(1111);
    h1D->Fit(total,"QMR+");
    TH1D *hFitKS = (TH1D*)h1D->Clone();
    TH1D *hTotalFit = FuncToHist(total,hFitKS,islice,arrMean,arrSigma);
    hTotalFit->SetName(Form("FitHist_%d",islice));
    hFitKS   ->SetName(Form("HistKS_%d",islice));
    hTotalFit->SetMarkerColor(kGreen);
    hFitKS   ->SetMarkerColor(kBlue);
    hTotalFit->SetLineColor(kRed);
    hFitKS   ->SetLineColor(kBlack);
//     Double_t binWidthTotal = hTotalFit->GetXaxis()->GetBinWidth(50);
//     hTotalFit->Scale(1./binWidthTotal);
    
    // Fill the residuals tree
    TVectorD *resVector  = new TVectorD(hResidual->GetNbinsX());
    TVectorD *perVector  = new TVectorD(hResidual->GetNbinsX());
    TVectorD *dEdxVector = new TVectorD(hResidual->GetNbinsX());
    Double_t res;
    Double_t per;
    for (Int_t ires = 1;ires<hResidual->GetNbinsX();ires++){
      if(hResidual->GetBinContent(ires)<1e-4) continue;
      per=(hResidual->GetBinContent(ires)-total->Eval(hResidual->GetBinCenter(ires)))/hResidual->GetBinContent(ires);
      res=(hResidual->GetBinContent(ires)-total->Eval(hResidual->GetBinCenter(ires)))/hResidual->GetBinError(ires);
      hResidual->SetBinContent(ires,res);
      perVector -> GetMatrixArray()[ires] = per;
      resVector -> GetMatrixArray()[ires] = res; 
      dEdxVector-> GetMatrixArray()[ires] = hResidual->GetBinCenter(ires);
    }
    Float_t maxPer = TMath::Max( TMath::Abs(perVector->Max()),TMath::Abs(perVector->Min()) ) ;
    Float_t maxRes = TMath::Max( TMath::Abs(resVector->Max()),TMath::Abs(resVector->Min()) );
    hResidual->SetOption("hist");
    
    // Set the final parameters of individual particle fit
    ResetParametersOfEachFunction(g1,g2,g3,g4,g5,g6,total,h1D,pt);
    
    // Dump histograms into arrays
    h1D       -> SetName(Form("FinalFit_%d_%4.3f"  ,islice,pt));
    hResidual -> SetName(Form("Residual_%d_%4.3f"  ,islice,pt));
    
    fitResArr[islice]    = (TH1D*)h1D;
    fitResiduals[islice] = (TH1D*)hResidual;
        
    // Dump some info into TTreeStream
        
    elFitAmp[islice]   = g1->GetParameter(0);
    piFitAmp[islice]   = g2->GetParameter(0);
    kaFitAmp[islice]   = g3->GetParameter(0);
    prFitAmp[islice]   = g4->GetParameter(0);
    deFitAmp[islice]   = g5->GetParameter(0);
    muFitAmp[islice]   = g6->GetParameter(0);
        
    elFitMean[islice]  = g1->GetParameter(1);
    piFitMean[islice]  = g2->GetParameter(1);
    kaFitMean[islice]  = g3->GetParameter(1);
    prFitMean[islice]  = g4->GetParameter(1);
    deFitMean[islice]  = g5->GetParameter(1);
    muFitMean[islice]  = g6->GetParameter(1);   
    
    elFitSigma[islice] = g1->GetParameter(2);
    piFitSigma[islice] = g2->GetParameter(2);
    kaFitSigma[islice] = g3->GetParameter(2);
    prFitSigma[islice] = g4->GetParameter(2);
    deFitSigma[islice] = g5->GetParameter(2);
    muFitSigma[islice] = g6->GetParameter(2);
    
    elFitKurtosis[islice]  = g1->GetParameter(3);
    piFitKurtosis[islice]  = g2->GetParameter(3);
    kaFitKurtosis[islice]  = g3->GetParameter(3);
    prFitKurtosis[islice]  = g4->GetParameter(3);
    deFitKurtosis[islice]  = g5->GetParameter(3);
    muFitKurtosis[islice]  = g6->GetParameter(3);
    
    elFitSkew[islice]  = g1->GetParameter(4);
    piFitSkew[islice]  = g2->GetParameter(4);
    kaFitSkew[islice]  = g3->GetParameter(4);
    prFitSkew[islice]  = g4->GetParameter(4);
    deFitSkew[islice]  = g5->GetParameter(4);
    muFitSkew[islice]  = g6->GetParameter(4);
        
    elChi2Fit[islice]  = g1->GetChisquare();
    piChi2Fit[islice]  = g2->GetChisquare();
    kaChi2Fit[islice]  = g3->GetChisquare();
    prChi2Fit[islice]  = g4->GetChisquare();
    deChi2Fit[islice]  = g5->GetChisquare();
    muChi2Fit[islice]  = g6->GetChisquare();
   
    elIntegralFit[islice] = ComputeGaussIntegral(elFitAmp[islice],elFitMean[islice],elFitSigma[islice],g1->GetParameter(3),g1->GetParameter(4));
    piIntegralFit[islice] = ComputeGaussIntegral(piFitAmp[islice],piFitMean[islice],piFitSigma[islice],g2->GetParameter(3),g2->GetParameter(4));
    kaIntegralFit[islice] = ComputeGaussIntegral(kaFitAmp[islice],kaFitMean[islice],kaFitSigma[islice],g3->GetParameter(3),g3->GetParameter(4));
    prIntegralFit[islice] = ComputeGaussIntegral(prFitAmp[islice],prFitMean[islice],prFitSigma[islice],g4->GetParameter(3),g4->GetParameter(4));
    deIntegralFit[islice] = ComputeGaussIntegral(deFitAmp[islice],deFitMean[islice],deFitSigma[islice],g5->GetParameter(3),g5->GetParameter(4));
    muIntegralFit[islice] = ComputeGaussIntegral(muFitAmp[islice],muFitMean[islice],muFitSigma[islice],g6->GetParameter(3),g6->GetParameter(4));
    
    Double_t totalChi2    = total->GetChisquare();
    Double_t normNDF      = total->GetNDF();
    Double_t totChi2      = (normNDF<1) ? 1 : totalChi2/normNDF;
    Double_t pValue       = total->GetProb();
    Double_t ksValue      = hFitKS->KolmogorovTest(hTotalFit);
//     Double_t chi2Value    = hFitKS->Chi2Test(hTotalFit);
    
    
    
    cout << " iter = "          << iIter                     ;
    cout << " ptDownBin  = "    << ptDownBin                 ;
    cout << " ptUpBin = "       << ptUpBin                   ;
    cout << " islice = "        << islice                    ; 
    cout << " p = "             << pt                 << endl;
    
    cout << " Pion Skewness = "   << skewnessFix[1]         ;
    cout << " Pion Kurtosis = "   << kurtosisFix[1]         ;
    cout << " Kaon Skewness = "   << skewnessFix[2]         ;
    cout << " Kaon Kurtosis = "   << kurtosisFix[2]         ;
    cout << " Proton Skewness = " << skewnessFix[3]         ;
    cout << " Proton Kurtosis = " << kurtosisFix[3]  << endl;
    
    cout << " ---------------------------------------------------------------------------------------------------> kstest = " << ksValue;
    cout << " chi2 = "          << totChi2            << endl;
    
        
    // Fit results dump ttrees
    fitFile->GetFile()->cd();
    *fitFile << "IdMethodInput" <<
  
        "iter="      << iIter                 <<
        "slice="     << islice                <<
        "p="         << pt                    <<
        "eta="       << etaBin                <<
        "cent="      << centBin               << 
        
        "totalChi2=" << totalChi2             <<
        "totChi2="   << totChi2               <<
        "ksTest="    << ksValue               <<
//         "chi2Test="  << chi2Value             <<
        "pValue="    << pValue                <<
        "normNDF="   << normNDF               ;
    
    *fitFile << "IdMethodInput" <<
        
        "elAmp="     << elFitAmp[islice]      <<
        "piAmp="     << piFitAmp[islice]      <<
        "kaAmp="     << kaFitAmp[islice]      <<
        "prAmp="     << prFitAmp[islice]      <<
        "deAmp="     << deFitAmp[islice]      <<
        "muAmp="     << muFitAmp[islice]      <<
        
        "elMean="    << elFitMean[islice]     <<
        "piMean="    << piFitMean[islice]     <<
        "kaMean="    << kaFitMean[islice]     <<
        "prMean="    << prFitMean[islice]     <<
        "deMean="    << deFitMean[islice]     <<
        "muMean="    << muFitMean[islice]     <<
        
        "elSigma="   << elFitSigma[islice]    <<
        "piSigma="   << piFitSigma[islice]    <<
        "kaSigma="   << kaFitSigma[islice]    <<
        "prSigma="   << prFitSigma[islice]    <<
        "deSigma="   << deFitSigma[islice]    <<
        "muSigma="   << muFitSigma[islice]     ;
        
    *fitFile << "IdMethodInput" <<
        
        "elKurtosis="    << elFitKurtosis[islice]     <<
        "piKurtosis="    << piFitKurtosis[islice]     <<
        "kaKurtosis="    << kaFitKurtosis[islice]     <<
        "prKurtosis="    << prFitKurtosis[islice]     <<
        "deKurtosis="    << deFitKurtosis[islice]     <<
        "muKurtosis="    << muFitKurtosis[islice]     <<
        
        "elSkew="    << elFitSkew[islice]     <<
        "piSkew="    << piFitSkew[islice]     <<
        "kaSkew="    << kaFitSkew[islice]     <<
        "prSkew="    << prFitSkew[islice]     <<
        "deSkew="    << deFitSkew[islice]     <<
        "muSkew="    << muFitSkew[islice]     <<
        
        "elInt="     << elIntegralFit[islice] << 
        "piInt="     << piIntegralFit[islice] << 
        "kaInt="     << kaIntegralFit[islice] << 
        "prInt="     << prIntegralFit[islice] << 
        "deInt="     << deIntegralFit[islice] << 
        "muInt="     << muIntegralFit[islice] << 
        
        "\n";
    
    // Fit results dump ttrees
    *fitFile << "FitResults" <<
        
        "elMean="       << arrMean[0]         <<
        "piMean="       << arrMean[1]         <<
        "kaMean="       << arrMean[2]         <<
        "prMean="       << arrMean[3]         <<
        "deMean="       << arrMean[4]         <<
        "muMean="       << arrMean[5]         <<
        
        "elSigma="      << arrSigma[0]        <<
        "piSigma="      << arrSigma[1]        <<
        "kaSigma="      << arrSigma[2]        <<
        "prSigma="      << arrSigma[3]        <<
        "deSigma="      << arrSigma[4]        <<
        "muSigma="      << arrSigma[5]        <<
        
        "elWindow="     << arrMeanWindow[0]   <<
        "piWindow="     << arrMeanWindow[1]   <<
        "kaWindow="     << arrMeanWindow[2]   <<
        "prWindow="     << arrMeanWindow[3]   ;
    
    *fitFile << "FitResults" <<
        
        "iter="         << iIter              <<
        "slice="        << islice             <<
        "p="            << pt                 <<
        "eta="          << etaBin             <<
        "cent="         << centBin            << 
        "totalChi2="    << totalChi2          <<
        "totChi2="      << totChi2            <<
        "ksTest="       << ksValue            <<
//         "chi2Test="     << chi2Value          <<
        "pValue="       << pValue             <<
        "normNDF="      << normNDF            <<
        "maxRes="       << maxRes             <<
        "maxPer="       << maxPer             <<
        "resVector.="   << resVector          <<
        "perVector.="   << perVector          << 
        "dEdx.="        << dEdxVector         <<
        
        "elFitInt="     << elIntegralFit[islice] << 
        "piFitInt="     << piIntegralFit[islice] << 
        "kaFitInt="     << kaIntegralFit[islice] << 
        "prFitInt="     << prIntegralFit[islice] << 
        "deFitInt="     << deIntegralFit[islice] << 
        "muFitInt="     << muIntegralFit[islice] <<
        
        "elFitMean="    << elFitMean[islice]  <<
        "piFitMean="    << piFitMean[islice]  <<
        "kaFitMean="    << kaFitMean[islice]  <<
        "prFitMean="    << prFitMean[islice]  <<
        "deFitMean="    << deFitMean[islice]  <<
        "muFitMean="    << muFitMean[islice]  ;
    
    *fitFile << "FitResults" <<
        
        "elFitSigma="   << elFitSigma[islice] <<
        "piFitSigma="   << piFitSigma[islice] <<
        "kaFitSigma="   << kaFitSigma[islice] <<
        "prFitSigma="   << prFitSigma[islice] <<
        "muFitSigma="   << muFitSigma[islice] <<
        
        "elFitAmp="     << elFitAmp[islice]   <<
        "piFitAmp="     << piFitAmp[islice]   <<
        "kaFitAmp="     << kaFitAmp[islice]   <<
        "prFitAmp="     << prFitAmp[islice]   <<
        "deFitAmp="     << deFitAmp[islice]   <<
        "muFitAmp="     << muFitAmp[islice]   <<
        
        "elKurtosis="    << elFitKurtosis[islice]     <<
        "piKurtosis="    << piFitKurtosis[islice]     <<
        "kaKurtosis="    << kaFitKurtosis[islice]     <<
        "prKurtosis="    << prFitKurtosis[islice]     <<
        "deKurtosis="    << deFitKurtosis[islice]     <<
        "muKurtosis="    << muFitKurtosis[islice]  ;
    
    *fitFile << "FitResults" <<
        
        "elFitSkew="    << elFitSkew[islice]  <<
        "piFitSkew="    << piFitSkew[islice]  <<
        "kaFitSkew="    << kaFitSkew[islice]  <<
        "prFitSkew="    << prFitSkew[islice]  <<
        "deFitSkew="    << deFitSkew[islice]  <<
        "muFitSkew="    << muFitSkew[islice]  <<
        
        "elChi2Fit="    << elChi2Fit[islice]  <<
        "piChi2Fit="    << piChi2Fit[islice]  <<
        "kaChi2Fit="    << kaChi2Fit[islice]  <<
        "prChi2Fit="    << prChi2Fit[islice]  <<
        "deChi2Fit="    << deChi2Fit[islice]  <<
        "muChi2Fit="    << muChi2Fit[islice]  <<
      
        "\n";
    
    // write hTotalFit for debugging
    histsFile -> GetFile()->cd();
    if (islice>10 && islice<20) {
      hTotalFit->Write(Form("FitHist_%d",islice));
      hFitKS->Write(Form("HistKS_%d",islice));
    }
    
    delete [] arrSigma;
    delete [] arrMean;
    delete [] arrMeanWindow;
    delete [] cleanParMean;
    delete [] MCPars;
    delete resVector;
    delete perVector;
    delete dEdxVector;
    delete hTotalFit;
    delete hFitKS;
    
    
  }
    
  histsFile -> GetFile()->cd();
//   if ((iIter==0 && nSlice<20) || (nSlice>100 && iIter>=3)) fitResArr . Write("fitResArr",TObject::kSingleKey);
  fitResArr . Write("fitResArr",TObject::kSingleKey);
  if (iIter >=3) fitResiduals . Write("fitResiduals",TObject::kSingleKey);
      
  // delete arrays
  fitResArr    . Clear("C");
  fitResiduals . Clear("C");
    
  delete histsFile;
  delete fitFile;
  
  if (ressFile!=0x0) ressFile->Close();
  if (sampFile!=0x0) sampFile->Close();
    
}
// -------------------------------------------------------------------------------------------------------
TGraphErrors *RemoveOutliers(TGraphErrors *grsmooth, Int_t fitWindow){
    
    
  Int_t startSmoothing = 0;
  Int_t entries = grsmooth->GetN()-startSmoothing;
    
  Double_t *xxx    = grsmooth->GetX(); 
  Double_t *yyy    = grsmooth->GetY();
//   Double_t *yyy    = new Double_t[entries];
  Double_t *yyyErr = new Double_t[entries];
  // avoid from negatif labels, they include some ga
  cout << " Smooth Fuction is being processed " << endl;
    
  for (Int_t ip=startSmoothing; ip<entries; ip++){
    
    TLinearFitter fpol2(3,smoothFunction);
    yyy[ip]=grsmooth->GetY()[ip];
    yyyErr[ip]=0;
    
    for (Int_t idelta=-fitWindow; idelta<fitWindow; idelta++){
      Int_t index = ip+idelta;
      if (index<startSmoothing || index>=entries) continue;
      Double_t err=1+TMath::Sqrt(TMath::Abs(idelta)/fitWindow);
      err=1+(idelta>1)*fitWindow;
      Double_t xx[4]={idelta,idelta,idelta,idelta};
      fpol2.AddPoint(xx,grsmooth->GetY()[index],err);
    }
    fpol2.Eval();
    fpol2.EvalRobust(0.65);
    yyy[ip]=fpol2.GetParameter(0);
    xxx[ip]=grsmooth->GetX()[ip];
  }
   
  TGraphErrors * grout = new TGraphErrors(entries,xxx,yyy,0,yyyErr);
  delete [] yyy;
  delete [] yyyErr;
  return grout;
    
}
// -------------------------------------------------------------------------------------------------------
void CheckIfMeanIsZero(TF1 *total, Double_t *arrMean){

  //
  // If a particle is not seen in the pt slice fix its fit parameters to 0
  //

  if(arrMean[0]<=10){
    total->FixParameter(0,0);
    total->FixParameter(1,0);
    total->FixParameter(2,0);
    total->FixParameter(3,0);
    total->FixParameter(4,0);
  }
  if(arrMean[1]<=10){
    total->FixParameter(5,0);
    total->FixParameter(6,0);
    total->FixParameter(7,0);
    total->FixParameter(8,0);
    total->FixParameter(9,0);
  }
  if(arrMean[2]<=10){
    total->FixParameter(10,0);
    total->FixParameter(11,0);
    total->FixParameter(12,0);
    total->FixParameter(13,0);
    total->FixParameter(14,0);
  }
  if(arrMean[3]<=10){
    total->FixParameter(15,0);
    total->FixParameter(16,0);
    total->FixParameter(17,0);
    total->FixParameter(18,0);
    total->FixParameter(19,0);
  }
  if(arrMean[4]<=10){
    total->FixParameter(20,0);
    total->FixParameter(21,0);
    total->FixParameter(22,0);
    total->FixParameter(23,0);
    total->FixParameter(24,0);
  }
  if(arrMean[5]<=10){
    total->FixParameter(25,0);
    total->FixParameter(26,0);
    total->FixParameter(27,0);
    total->FixParameter(28,0);
    total->FixParameter(29,0);
  }

}
// -------------------------------------------------------------------------------------------------------
void SetParNamesOfTotalFit(TF1 *total){
  
  //
  // Set Parameter Names of total function
  //
  
  total->SetParName(0,"elAmplitude");
  total->SetParName(1,"elMean");
  total->SetParName(2,"elSigma");
  total->SetParName(3,"elKurtosis");
  total->SetParName(4,"elSkewness");
  
  total->SetParName(5,"piAmplitude");
  total->SetParName(6,"piMean");
  total->SetParName(7,"piSigma");
  total->SetParName(8,"piKurtosis");
  total->SetParName(9,"piSkewness");
  
  total->SetParName(10,"kaAmplitude");
  total->SetParName(11,"kaMean");
  total->SetParName(12,"kaSigma");
  total->SetParName(13,"kaKurtosis");
  total->SetParName(14,"kaSkewness");
  
  total->SetParName(15,"prAmplitude");
  total->SetParName(16,"prMean");
  total->SetParName(17,"prSigma");
  total->SetParName(18,"prKurtosis");
  total->SetParName(19,"prSkewness");
  
  total->SetParName(20,"deAmplitude");
  total->SetParName(21,"deMean");
  total->SetParName(22,"deSigma");
  total->SetParName(23,"deKurtosis");
  total->SetParName(24,"deSkewness");
  
  total->SetParName(25,"muAmplitude");
  total->SetParName(26,"muMean");
  total->SetParName(27,"muSigma");
  total->SetParName(28,"muKurtosis");
  total->SetParName(29,"muSkewness");
}
// --------------------------------------------------------------------------------------------
void SetFixedKSparameterForTotalFit(TF1 *total){

  // Set kurtosis 
  if (fixedK){
    total->FixParameter(3, kurtosisFix[kElectron]); 
    total->FixParameter(8, kurtosisFix[kPion]    ); 
    total->FixParameter(13,kurtosisFix[kKaon]    ); 
    total->FixParameter(18,kurtosisFix[kProton]  ); 
    total->FixParameter(23,kurtosisFix[kDeuteron]); 
    total->FixParameter(28,kurtosisFix[kMuon]    ); 
  }
  
  // Set skewness
  if (fixedS) {
    total->FixParameter(4, skewnessFix[kElectron]); 
    total->FixParameter(9, skewnessFix[kPion]    ); 
    total->FixParameter(14,skewnessFix[kKaon]    ); 
    total->FixParameter(19,skewnessFix[kProton]  ); 
    total->FixParameter(24,skewnessFix[kDeuteron]); 
    total->FixParameter(29,skewnessFix[kMuon]    ); 
  }

  
}
// --------------------------------------------------------------------------------------------
void SetFixedKSparameterForTotalFitForPions(TF1 *total){
  
  // Set kurtosis 
  total->FixParameter(3, kurtosisFix[kElectron]); 
  
  total->SetParLimits(8,kurtosisMin,kurtosisMax);
  
  total->FixParameter(13,kurtosisFix[kKaon]    ); 
  total->FixParameter(18,kurtosisFix[kProton]  ); 
  total->FixParameter(23,kurtosisFix[kDeuteron]); 
  total->FixParameter(28,kurtosisFix[kMuon]    ); 
  
  // Set skewness
  total->FixParameter(4, skewnessFix[kElectron]); 
  
  total->SetParLimits(9,skewMin,skewMax);
  
  total->FixParameter(14,skewnessFix[kKaon]    ); 
  total->FixParameter(19,skewnessFix[kProton]  ); 
  total->FixParameter(24,skewnessFix[kDeuteron]); 
  total->FixParameter(29,skewnessFix[kMuon]    ); 

}
// --------------------------------------------------------------------------------------------
void SetFixedKSparameterForIndividualFits(TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4, TF1 *g5, TF1 *g6){

  // Set kurtosis
  if (fixedK) {
    g1->FixParameter(3,kurtosisFix[kElectron]); 
    g2->FixParameter(3,kurtosisFix[kPion]    ); 
    g3->FixParameter(3,kurtosisFix[kKaon]    ); 
    g4->FixParameter(3,kurtosisFix[kProton]  ); 
    g5->FixParameter(3,kurtosisFix[kDeuteron]); 
    g6->FixParameter(3,kurtosisFix[kMuon]    ); 
  }

  
  // Set Skewness 
  if (fixedS) {
    g1->FixParameter(4,skewnessFix[kElectron]); 
    g2->FixParameter(4,skewnessFix[kPion]    ); 
    g3->FixParameter(4,skewnessFix[kKaon]    ); 
    g4->FixParameter(4,skewnessFix[kProton]  ); 
    g5->FixParameter(4,skewnessFix[kDeuteron]); 
    g6->FixParameter(4,skewnessFix[kMuon]    ); 
  }

}
// -------------------------------------------------------------------------------------------------------
void SetFixedKSparameterForCleanSampleFit(TF1 *asyGaus, ParticleType pSpecy){

 
      // Fix kurtoris and skewness for a given particle specy
  switch (pSpecy) 
  { 
    case kElectron: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kElectron]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kElectron]);
    }
    break; 
    case kPion: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kPion]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kPion]);
    }
    break; 
    case kKaon: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kKaon]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kKaon]);
    }
    break;  
    case kProton: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kProton]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kProton]);
    }
    break; 
    case kDeuteron: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kDeuteron]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kDeuteron]);
    }
    break; 
    case kMuon: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kMuon]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kMuon]);
    }
    break; 
  } 
  
}
// -------------------------------------------------------------------------------------------------------
void SetFixedKSparameterForMCSampleFit(TF1 *asyGaus, ParticleType pSpecy){

 
      // Fix kurtoris and skewness for a given particle specy
  switch (pSpecy) 
  { 
    case kElectron: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kElectron]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kElectron]);
    }
    break; 
    case kPion: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kPion]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kPion]);
    }
    break; 
    case kKaon: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kKaon]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kKaon]);
    }
    break;  
    case kProton: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kProton]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kProton]);
    }
    break; 
    case kDeuteron: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kDeuteron]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kDeuteron]);
    }
    break; 
    case kMuon: 
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kMuon]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kMuon]);
    }
    break; 
  } 
  
}
// -------------------------------------------------------------------------------------------------------
void GetHelperGraphs(Int_t iIter, TFile *ressFile){
  
  //
  // Read helper graphs from the smooth tree
  //

  if (iIter>0) {
//     get the file and objarrays of previous fit results
    Fit           = (TObjArray*)ressFile->Get("Fit");
    Clean         = (TObjArray*)ressFile->Get("Clean");
    MCs           = (TObjArray*)ressFile->Get("MC");
    ScaledClean   = (TObjArray*)ressFile->Get("ScaledClean");
    ScaledMC      = (TObjArray*)ressFile->Get("ScaledMC");
    ScaledEx      = (TObjArray*)ressFile->Get("ScaledEx");
    OutlierSmooth = (TObjArray*)ressFile->Get("OulierSmooth");
    Windows       = (TObjArray*)ressFile->Get("Windows");
    HighPFit      = (TObjArray*)ressFile->Get("HighPFit");
    grpiToMuRatio   = (TGraphErrors*)ressFile->Get("piToMuRatio");
    
//     graphs for electron
    grelFitAmpSmooth           = (TGraphErrors*)Fit          ->FindObject("elFitAmpSmooth");
    grelFitMeanSmooth          = (TGraphErrors*)Fit          ->FindObject("elFitMeanSmooth");
    grelFitSigmaSmooth         = (TGraphErrors*)Fit          ->FindObject("elFitSigmaSmooth");
    grelMeanExScaled           = (TGraphErrors*)ScaledEx     ->FindObject("elMeanExScaled");
    grelSigmaExScaled          = (TGraphErrors*)ScaledEx     ->FindObject("elSigmaExScaled");
    grelFitAmpOutlierSmooth    = (TGraphErrors*)OutlierSmooth->FindObject("elFitAmpOutlierSmooth");
    
//     graphs for pion 
    grpiFitAmpSmooth           = (TGraphErrors*)Fit          ->FindObject("piFitAmpSmooth");
    grpiFitMeanSmooth          = (TGraphErrors*)Fit          ->FindObject("piFitMeanSmooth");
    grpiFitSigmaSmooth         = (TGraphErrors*)Fit          ->FindObject("piFitSigmaSmooth");
    grpiCleanMeanSmooth        = (TGraphErrors*)Clean        ->FindObject("piCleanMeanSmooth");
    grpiCleanSigmaScaledSmooth = (TGraphErrors*)ScaledClean  ->FindObject("piCleanSigmaScaledSmooth");
    grpiCleanMeanScaledSmooth  = (TGraphErrors*)ScaledClean  ->FindObject("piCleanMeanScaledSmooth");
    grpiSigmaExScaled          = (TGraphErrors*)ScaledEx     ->FindObject("piSigmaExScaled");
    grpiMeanExScaled           = (TGraphErrors*)ScaledEx     ->FindObject("piMeanExScaled");


//     graphs for kaon 
    grkaFitAmpSmooth           = (TGraphErrors*)Fit          ->FindObject("kaFitAmpSmooth");
    grkaFitMeanSmooth          = (TGraphErrors*)Fit          ->FindObject("kaFitMeanSmooth");
    grkaFitSigmaSmooth         = (TGraphErrors*)Fit          ->FindObject("kaFitSigmaSmooth");
    grkaCleanMean              = (TGraphErrors*)Clean        ->FindObject("kaCleanMean");
    grkaCleanSigma             = (TGraphErrors*)Clean        ->FindObject("kaCleanSigma");
    grkaMeanExScaled           = (TGraphErrors*)ScaledEx     ->FindObject("kaMeanExScaled");
    grkaSigmaExScaled          = (TGraphErrors*)ScaledEx     ->FindObject("kaSigmaExScaled");
    grkaCleanAmpScaledSmooth   = (TGraphErrors*)ScaledClean  ->FindObject("kaCleanAmpScaledSmooth");
    grkaCleanSigmaScaled       = (TGraphErrors*)ScaledClean  ->FindObject("kaCleanSigmaScaled");
    grkaFitAmpOutlierSmooth    = (TGraphErrors*)OutlierSmooth->FindObject("kaFitAmpOutlierSmooth");
    grkaFitSigmaOutlierSmooth  = (TGraphErrors*)OutlierSmooth->FindObject("kaFitSigmaOutlierSmooth");

    
//     graphs for proton 
    grprFitAmpSmooth           = (TGraphErrors*)Fit          ->FindObject("prFitAmpSmooth");
    grprFitMeanSmooth          = (TGraphErrors*)Fit          ->FindObject("prFitMeanSmooth");
    grprFitSigmaSmooth         = (TGraphErrors*)Fit          ->FindObject("prFitSigmaSmooth");
    grprCleanMean              = (TGraphErrors*)Clean        ->FindObject("prCleanMean");
    grprCleanSigma             = (TGraphErrors*)Clean        ->FindObject("prCleanSigma");
    grprSigmaExScaled          = (TGraphErrors*)ScaledEx     ->FindObject("prSigmaExScaled");
    grprMeanExScaled           = (TGraphErrors*)ScaledEx     ->FindObject("prMeanExScaled");
    grprCleanAmpScaledSmooth   = (TGraphErrors*)ScaledClean  ->FindObject("prCleanAmpScaledSmooth");
    grprCleanSigmaScaledSmooth = (TGraphErrors*)ScaledClean  ->FindObject("prCleanSigmaScaledSmooth");
    grprFitAmpOutlierSmooth    = (TGraphErrors*)OutlierSmooth->FindObject("prFitAmpOutlierSmooth");

//     graphs for muon 
    grmuMeanExScaled           = (TGraphErrors*)ScaledEx     ->FindObject("muMeanExScaled");
    grmuSigmaExScaled          = (TGraphErrors*)ScaledEx     ->FindObject("muSigmaExScaled");

    //     MC scaled graphs
    grpiMCAmpScaledSmooth      = (TGraphErrors*)ScaledMC     ->FindObject("piMCAmpScaledSmooth");
    grpiMCMeanScaledSmooth     = (TGraphErrors*)ScaledMC     ->FindObject("piMCMeanScaledSmooth");
    grpiMCSigmaScaledSmooth    = (TGraphErrors*)ScaledMC     ->FindObject("piMCSigmaScaledSmooth");
    
    grkaMCAmpScaledSmooth      = (TGraphErrors*)ScaledMC     ->FindObject("kaMCAmpScaledSmooth");
    grkaMCMeanScaledSmooth     = (TGraphErrors*)ScaledMC     ->FindObject("kaMCMeanScaledSmooth");
    grkaMCSigmaScaledSmooth    = (TGraphErrors*)ScaledMC     ->FindObject("kaMCSigmaScaledSmooth");
    
    grprMCAmpScaledSmooth      = (TGraphErrors*)ScaledMC     ->FindObject("prMCAmpScaledSmooth");
    grprMCMeanScaledSmooth     = (TGraphErrors*)ScaledMC     ->FindObject("prMCMeanScaledSmooth");
    grprMCSigmaScaledSmooth    = (TGraphErrors*)ScaledMC     ->FindObject("prMCSigmaScaledSmooth");
    
//     High momentum graphs
    grpiFitAmpSmoothHP         = (TGraphErrors*)HighPFit     ->FindObject("piFitAmpSmooth");
    grkaFitAmpSmoothHP         = (TGraphErrors*)HighPFit     ->FindObject("kaFitAmpSmooth");
    grprFitAmpSmoothHP         = (TGraphErrors*)HighPFit     ->FindObject("prFitAmpSmooth");
    
    grpiFitMeanSmoothHP        = (TGraphErrors*)HighPFit     ->FindObject("piFitMeanSmooth");
    grkaFitMeanSmoothHP        = (TGraphErrors*)HighPFit     ->FindObject("kaFitMeanSmooth");
    grprFitMeanSmoothHP        = (TGraphErrors*)HighPFit     ->FindObject("prFitMeanSmooth");
    
//     For h2DMC analysis
    grelMCSigmaScaled          = (TGraphErrors*)ScaledMC->FindObject("elMCSigmaScaled");
    grelMCMeanScaled           = (TGraphErrors*)ScaledMC->FindObject("elMCMeanScaled");
    grelMCAmpScaled            = (TGraphErrors*)ScaledMC->FindObject("elMCAmpScaled");

    grpiMCSigmaScaled          = (TGraphErrors*)ScaledMC->FindObject("piMCSigmaScaled");
    grpiMCMeanScaled           = (TGraphErrors*)ScaledMC->FindObject("piMCMeanScaled");
    grpiMCAmpScaled            = (TGraphErrors*)ScaledMC->FindObject("piMCAmpScaled");
  
    grkaMCSigmaScaled          = (TGraphErrors*)ScaledMC->FindObject("kaMCSigmaScaled");
    grkaMCMeanScaled           = (TGraphErrors*)ScaledMC->FindObject("kaMCMeanScaled");
  
    grprMCSigmaScaled          = (TGraphErrors*)ScaledMC->FindObject("prMCSigmaScaled");
    grprMCMeanScaled           = (TGraphErrors*)ScaledMC->FindObject("prMCMeanScaled");
  
    grmuSigmaExScaled          = (TGraphErrors*)ScaledEx->FindObject("muSigmaExScaled");
    grmuMeanExScaled           = (TGraphErrors*)ScaledEx->FindObject("muMeanExScaled");

  }

}
// -------------------------------------------------------------------------------------------------------
void SetTotalFitParameters(Int_t iIter, Int_t islice, Int_t nSlice, TF1* total, Double_t *arrMean, Double_t *arrSigma, Double_t *arrMeanWindow, Double_t *cleanParMean, Double_t * MCPars, TH1D* h1D, Double_t pt){

  //
  // Set total Fit Parameters for each iteration
  //
    
  if (islice<nSlice) {
    
    // Set the mean position wrt closest particle
    arrMeanWindow[0] = GetClosestParticleMean(islice,arrSigma[0],arrMean[0],arrMean[1],arrMean[2],arrMean[3]);
    arrMeanWindow[1] = GetClosestParticleMean(islice,arrSigma[1],arrMean[1],arrMean[0],arrMean[2],arrMean[3]);
    arrMeanWindow[2] = GetClosestParticleMean(islice,arrSigma[2],arrMean[2],arrMean[1],arrMean[0],arrMean[3]);
    arrMeanWindow[3] = GetClosestParticleMean(islice,arrSigma[3],arrMean[3],arrMean[1],arrMean[2],arrMean[0]);

    //   set 1 iter params
    if (iIter>=0) SetParams1stIteration(pt,h1D,arrMean,arrSigma,arrMeanWindow,cleanParMean);
         
    if (iIter>=0 && MC) SetParams1stIterationMC(pt,h1D,arrMean,arrSigma,MCPars);

    
    if (iIter>=1) { 
      
      muAmpMCscaled = (grpiToMuRatio->GetY()[islice]!=0) ? grpiFitAmpSmooth->GetY()[islice]/(grpiToMuRatio->GetY()[islice]) : 1.;
      
      Double_t elAssump, piAssump, kaAssump, prAssump;
      if (assump==0){
        elAssump = grelMeanExScaled  -> GetY()[islice];
        piAssump = ((pt>=0.6 && pt<=1.8)) ? grpiCleanMeanScaledSmooth -> GetY()[islice]: grpiFitMeanSmooth -> GetY()[islice];
        kaAssump = grkaMeanExScaled  -> GetY()[islice];
        prAssump = grprMeanExScaled  -> GetY()[islice];
      } else if (assump==1){
        elAssump = grelMeanExScaled  -> GetY()[islice];
        piAssump = ((pt>=0.6 && pt<=1.8)) ? grpiCleanMeanScaledSmooth -> GetY()[islice]: grpiFitMeanSmooth -> GetY()[islice];
        kaAssump = grkaMeanExScaled  -> GetY()[islice];
        prAssump = grprCleanMean     -> GetY()[islice];
      } else if (assump==2) {
        elAssump = grelMeanExScaled  -> GetY()[islice];
        piAssump = ((pt>=0.6 && pt<=1.8)) ? grpiCleanMeanScaledSmooth -> GetY()[islice]: grpiFitMeanSmooth -> GetY()[islice];
        kaAssump = (grkaMeanExScaled -> GetY()[islice] + grkaCleanMean-> GetY()[islice])/2.;
        prAssump = (grprMeanExScaled -> GetY()[islice] + grprCleanMean-> GetY()[islice])/2.;
      }
      
    // Recalculate the windows with the clean samples 
      arrMeanWindow[0] = GetClosestParticleMean(islice, arrSigma[0], elAssump, piAssump, kaAssump, prAssump);
      arrMeanWindow[1] = GetClosestParticleMean(islice, arrSigma[1], piAssump, elAssump, kaAssump, prAssump);
      arrMeanWindow[2] = GetClosestParticleMean(islice, arrSigma[2], kaAssump, piAssump, elAssump, prAssump);
      arrMeanWindow[3] = GetClosestParticleMean(islice, arrSigma[3], prAssump, piAssump, kaAssump, elAssump);
    
      if(MC) { 
        SetParamsMC(islice,pt);
      } else { 
    
    // ***************  Electron ***************                 good  
        elMeanMin  = grelMeanExScaled  ->GetY()[islice] - arrMeanWindow[0];
        elMeanMax  = grelMeanExScaled  ->GetY()[islice] + arrMeanWindow[0];
        elSigmaMin = grelSigmaExScaled ->GetY()[islice] - grelSigmaExScaled->GetY()[islice]*0.1;
        elSigmaMax = grelSigmaExScaled ->GetY()[islice] + grelSigmaExScaled->GetY()[islice]*0.1;
    
    // ***************  Pion  *************** 
        if (pt>=0.6){   
          piMeanMin  = grpiFitMeanSmooth  ->GetY()[islice] - arrMeanWindow[1];
          piMeanMax  = grpiFitMeanSmooth  ->GetY()[islice] + arrMeanWindow[1];
          piSigmaMin = grpiFitSigmaSmooth ->GetY()[islice] - grpiFitSigmaSmooth->GetY()[islice]*0.05;
          piSigmaMax = grpiFitSigmaSmooth ->GetY()[islice] + grpiFitSigmaSmooth->GetY()[islice]*0.05;
          
          if (pt>=0.6 && pt<=1.8){   // in this interval clean samples are ok  
            piMeanMin  = grpiCleanMeanScaledSmooth  ->GetY()[islice] - arrMeanWindow[1];
            piMeanMax  = grpiCleanMeanScaledSmooth  ->GetY()[islice] + arrMeanWindow[1];
          }
 
        }
    // ***************  kaon  *************** 
        if (pt>=0.7){
          kaMeanMin  = grkaMeanExScaled  ->GetY()[islice] - arrMeanWindow[2];
          kaMeanMax  = grkaMeanExScaled  ->GetY()[islice] + arrMeanWindow[2];
          kaSigmaMin = grkaFitSigmaSmooth ->GetY()[islice] - grkaFitSigmaSmooth->GetY()[islice]*0.1;
          kaSigmaMax = grkaFitSigmaSmooth ->GetY()[islice] + grkaFitSigmaSmooth->GetY()[islice]*0.1;
        }
        
    // ***************  Proton  *************** 
        if (pt>=0.8){ 
          prMeanMin  = grprMeanExScaled  ->GetY()[islice]  - arrMeanWindow[3];
          prMeanMax  = grprMeanExScaled  ->GetY()[islice]  + arrMeanWindow[3];
          prSigmaMin = grprCleanSigmaScaledSmooth ->GetY()[islice] - grprCleanSigmaScaledSmooth ->GetY()[islice]*0.1;
          prSigmaMax = grprCleanSigmaScaledSmooth ->GetY()[islice] + grprCleanSigmaScaledSmooth ->GetY()[islice]*0.1;
        }
        
        if (pp && (pt>=0.8 && pt<=1.3)){
          prSigmaMin = grprCleanSigma ->GetY()[islice] - grprCleanSigma ->GetY()[islice]*0.1;
          prSigmaMax = grprCleanSigma ->GetY()[islice] + grprCleanSigma ->GetY()[islice]*0.1;
        }
        
    // ***************  Muon --> Fixed  ***************                  
        muMeanMin  = grmuMeanExScaled ->GetY()[islice] - grmuMeanExScaled ->GetY()[islice]*0.01;
        muMeanMax  = grmuMeanExScaled ->GetY()[islice] + grmuMeanExScaled ->GetY()[islice]*0.01;
        muSigmaMin = grmuSigmaExScaled->GetY()[islice] - grmuSigmaExScaled->GetY()[islice]*0.01;
        muSigmaMax = grmuSigmaExScaled->GetY()[islice] + grmuSigmaExScaled->GetY()[islice]*0.01;
        muAmpMin   = muAmpMCscaled                     - muAmpMCscaled*0.005;
        muAmpMax   = muAmpMCscaled                     + muAmpMCscaled*0.005; 
        
        //============================= SOME SAFE PATCHES ===============================================
      
      // electron amplitude tends to rise at el-K and el-pr overlap --> restrict elAmpMax
        if ((pt>=0.45 && pt<=0.65) || (pt>=0.8 && pt<=1.3) ){       
          elAmpMax  = grelFitAmpSmooth->GetY()[islice];     
        }
        
      }
       
    }

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

    if (iIter>=2 && !MC) {  
            
      //  el-K and el-pr overlap
      if ((pt>=0.45 && pt<=0.65) || (pt>=0.8 && pt<=1.3) ){  
             
        elAmpMax  = grelFitAmpSmooth->GetY()[islice];  
        
        kaMeanMin  = grkaMeanExScaled ->GetY()[islice] - grkaMeanExScaled ->GetY()[islice]*0.005;
        kaMeanMax  = grkaMeanExScaled ->GetY()[islice] + grkaMeanExScaled ->GetY()[islice]*0.005;
        prMeanMin  = grprMeanExScaled ->GetY()[islice] - grprMeanExScaled ->GetY()[islice]*0.005;
        prMeanMax  = grprMeanExScaled ->GetY()[islice] + grprMeanExScaled ->GetY()[islice]*0.005;
      }
      
      // pi-K and pi-pr overlap
      if ((pt>=0.8 && pt<=1.4) || (pt>=1.4 && pt<=2.)){  
                       
        piMeanMin  = grpiCleanMeanScaledSmooth  ->GetY()[islice] - grpiCleanMeanScaledSmooth ->GetY()[islice]*0.01;
        piMeanMax  = grpiCleanMeanScaledSmooth  ->GetY()[islice] + grpiCleanMeanScaledSmooth ->GetY()[islice]*0.01;
        kaMeanMin  = grkaMeanExScaled ->GetY()[islice] - grkaMeanExScaled ->GetY()[islice]*0.01;
        kaMeanMax  = grkaMeanExScaled ->GetY()[islice] + grkaMeanExScaled ->GetY()[islice]*0.01;
        prMeanMin  = grprMeanExScaled ->GetY()[islice] - grprMeanExScaled ->GetY()[islice]*0.01;
        prMeanMax  = grprMeanExScaled ->GetY()[islice] + grprMeanExScaled ->GetY()[islice]*0.01;
        
      }
      
       //     Kaon Amp
//       if (pt>=0.8 && pt<=1.4) kaAmpMax   = grkaFitAmpSmooth -> GetY()[islice];
      
      //     Proton Amp
//       if ((pt>=1.4 && pt<=2.)) prAmpMax   = grprFitAmpSmooth -> GetY()[islice];
       
    }
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    
    if (iIter>=3 && !MC) {  
      
     // ***************  electron merges with kaon and proton ***************   
      if ((pt>=0.45 && pt<=0.65) || (pt>=0.8 && pt<=1.3) ){       
        elAmpMin  = grelFitAmpSmooth->GetY()[islice] - grelFitAmpSmooth->GetY()[islice]*0.01;
        elAmpMax  = grelFitAmpSmooth->GetY()[islice] + grelFitAmpSmooth->GetY()[islice]*0.01;
      }
       
       //     Kaon
      if (pt>=0.85) {
        kaAmpMin   = grkaFitAmpSmooth   ->GetY()[islice] - grkaFitAmpSmooth   ->GetY()[islice]*0.1;
        kaAmpMax   = grkaFitAmpSmooth   ->GetY()[islice] + grkaFitAmpSmooth   ->GetY()[islice]*0.1;
        kaMeanMin  = grkaFitMeanSmooth  ->GetY()[islice] - grkaFitMeanSmooth  ->GetY()[islice]*0.01;
        kaMeanMax  = grkaFitMeanSmooth  ->GetY()[islice] + grkaFitMeanSmooth  ->GetY()[islice]*0.01;
      }
     
    //     Proton
      if (pt>1.) {
        prAmpMin   = grprFitAmpSmooth   ->GetY()[islice] - grprFitAmpSmooth   ->GetY()[islice]*0.1;
        prAmpMax   = grprFitAmpSmooth   ->GetY()[islice] + grprFitAmpSmooth   ->GetY()[islice]*0.1;
        prMeanMin  = grprFitMeanSmooth  ->GetY()[islice] - grprFitMeanSmooth  ->GetY()[islice]*0.01;
        prMeanMax  = grprFitMeanSmooth  ->GetY()[islice] + grprFitMeanSmooth  ->GetY()[islice]*0.01;
      } 
                               
    }
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

    if (iIter>=4 && !MC) {  
      
      //  pion fix 
      if (pt>=0.85) {
        piAmpMin   = grpiFitAmpSmooth   ->GetY()[islice] - grpiFitAmpSmooth   ->GetY()[islice]*0.01;
        piAmpMax   = grpiFitAmpSmooth   ->GetY()[islice] + grpiFitAmpSmooth   ->GetY()[islice]*0.01;
        piMeanMin  = grpiFitMeanSmooth  ->GetY()[islice] - grpiFitMeanSmooth  ->GetY()[islice]*0.01;
        piMeanMax  = grpiFitMeanSmooth  ->GetY()[islice] + grpiFitMeanSmooth  ->GetY()[islice]*0.01;
        piSigmaMin = grpiFitSigmaSmooth ->GetY()[islice] - grpiFitSigmaSmooth ->GetY()[islice]*0.01;
        piSigmaMax = grpiFitSigmaSmooth ->GetY()[islice] + grpiFitSigmaSmooth ->GetY()[islice]*0.01;
      }
       
    }   
      
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    
    if (iIter>=5 && !MC) {   
    
//     Electron 
      if (pt>=0.4){
        elAmpMin   = grelFitAmpSmooth   ->GetY()[islice] - grelFitAmpSmooth   ->GetY()[islice]*0.005;
        elAmpMax   = grelFitAmpSmooth   ->GetY()[islice] + grelFitAmpSmooth   ->GetY()[islice]*0.005;
        elSigmaMin = grelFitSigmaSmooth ->GetY()[islice] - grelFitSigmaSmooth ->GetY()[islice]*0.05;
        elSigmaMax = grelFitSigmaSmooth ->GetY()[islice] + grelFitSigmaSmooth ->GetY()[islice]*0.05;
      }
    
//     Pion
      if (pt>0.75) {
        piAmpMin   = grpiFitAmpSmooth   ->GetY()[islice] - grpiFitAmpSmooth   ->GetY()[islice]*0.005;
        piAmpMax   = grpiFitAmpSmooth   ->GetY()[islice] + grpiFitAmpSmooth   ->GetY()[islice]*0.005;
        piMeanMin  = grpiFitMeanSmooth  ->GetY()[islice] - grpiFitMeanSmooth  ->GetY()[islice]*0.005;
        piMeanMax  = grpiFitMeanSmooth  ->GetY()[islice] + grpiFitMeanSmooth  ->GetY()[islice]*0.005;
      }
    
//     Kaon
      if (pt>=0.75) {
        kaAmpMin   = grkaFitAmpSmooth   ->GetY()[islice] - grkaFitAmpSmooth   ->GetY()[islice]*0.05;
        kaAmpMax   = grkaFitAmpSmooth   ->GetY()[islice] + grkaFitAmpSmooth   ->GetY()[islice]*0.05;
        kaMeanMin  = grkaFitMeanSmooth  ->GetY()[islice] - grkaFitMeanSmooth  ->GetY()[islice]*0.005;
        kaMeanMax  = grkaFitMeanSmooth  ->GetY()[islice] + grkaFitMeanSmooth  ->GetY()[islice]*0.005;
        kaSigmaMin = grkaFitSigmaSmooth ->GetY()[islice] - grkaFitSigmaSmooth ->GetY()[islice]*0.05;
        kaSigmaMax = grkaFitSigmaSmooth ->GetY()[islice] + grkaFitSigmaSmooth ->GetY()[islice]*0.05;
      }
     
//     Proton
      if (pt>1.3) {
        prAmpMin   = grprFitAmpSmooth   ->GetY()[islice] - grprFitAmpSmooth   ->GetY()[islice]*0.05;
        prAmpMax   = grprFitAmpSmooth   ->GetY()[islice] + grprFitAmpSmooth   ->GetY()[islice]*0.05;
        prMeanMin  = grprFitMeanSmooth  ->GetY()[islice] - grprFitMeanSmooth  ->GetY()[islice]*0.005;
        prMeanMax  = grprFitMeanSmooth  ->GetY()[islice] + grprFitMeanSmooth  ->GetY()[islice]*0.005;
        prSigmaMin = grprFitSigmaSmooth ->GetY()[islice] - grprFitSigmaSmooth ->GetY()[islice]*0.05;
        prSigmaMax = grprFitSigmaSmooth ->GetY()[islice] + grprFitSigmaSmooth ->GetY()[islice]*0.05;
      }
     
    }
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

    if (iIter>=6 && !MC) {   
    
//     Electron 
      if (pt>=0.4){
        elAmpMin   = grelFitAmpSmooth   ->GetY()[islice] - grelFitAmpSmooth   ->GetY()[islice]*0.003;
        elAmpMax   = grelFitAmpSmooth   ->GetY()[islice] + grelFitAmpSmooth   ->GetY()[islice]*0.003;
        elSigmaMin = grelFitSigmaSmooth ->GetY()[islice] - grelFitSigmaSmooth ->GetY()[islice]*0.05;
        elSigmaMax = grelFitSigmaSmooth ->GetY()[islice] + grelFitSigmaSmooth ->GetY()[islice]*0.05;
      }
    
//     Pion
      if (pt>0.75) {
        piAmpMin   = grpiFitAmpSmooth   ->GetY()[islice] - grpiFitAmpSmooth   ->GetY()[islice]*0.003;
        piAmpMax   = grpiFitAmpSmooth   ->GetY()[islice] + grpiFitAmpSmooth   ->GetY()[islice]*0.003;
        piMeanMin  = grpiFitMeanSmooth  ->GetY()[islice] - grpiFitMeanSmooth  ->GetY()[islice]*0.005;
        piMeanMax  = grpiFitMeanSmooth  ->GetY()[islice] + grpiFitMeanSmooth  ->GetY()[islice]*0.005;
      }
    
//     Kaon
      if (pt>=0.75) {
        kaAmpMin   = grkaFitAmpSmooth   ->GetY()[islice] - grkaFitAmpSmooth   ->GetY()[islice]*0.01;
        kaAmpMax   = grkaFitAmpSmooth   ->GetY()[islice] + grkaFitAmpSmooth   ->GetY()[islice]*0.01;
        kaMeanMin  = grkaFitMeanSmooth  ->GetY()[islice] - grkaFitMeanSmooth  ->GetY()[islice]*0.005;
        kaMeanMax  = grkaFitMeanSmooth  ->GetY()[islice] + grkaFitMeanSmooth  ->GetY()[islice]*0.005;
        kaSigmaMin = grkaFitSigmaSmooth ->GetY()[islice] - grkaFitSigmaSmooth ->GetY()[islice]*0.05;
        kaSigmaMax = grkaFitSigmaSmooth ->GetY()[islice] + grkaFitSigmaSmooth ->GetY()[islice]*0.05;
      }
     
//     Proton
      if (pt>1.3) {
        prAmpMin   = grprFitAmpSmooth   ->GetY()[islice] - grprFitAmpSmooth   ->GetY()[islice]*0.01;
        prAmpMax   = grprFitAmpSmooth   ->GetY()[islice] + grprFitAmpSmooth   ->GetY()[islice]*0.01;
        prMeanMin  = grprFitMeanSmooth  ->GetY()[islice] - grprFitMeanSmooth  ->GetY()[islice]*0.005;
        prMeanMax  = grprFitMeanSmooth  ->GetY()[islice] + grprFitMeanSmooth  ->GetY()[islice]*0.005;
        prSigmaMin = grprFitSigmaSmooth ->GetY()[islice] - grprFitSigmaSmooth ->GetY()[islice]*0.05;
        prSigmaMax = grprFitSigmaSmooth ->GetY()[islice] + grprFitSigmaSmooth ->GetY()[islice]*0.05;
      }
     
    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

    if (iIter>=6 && MC) {   
    
//     Electron 
      if (pt>=0.35){
        elAmpMin   = grelFitAmpSmooth   ->GetY()[islice] - grelFitAmpSmooth   ->GetY()[islice]*0.05;
        elAmpMax   = grelFitAmpSmooth   ->GetY()[islice] + grelFitAmpSmooth   ->GetY()[islice]*0.05;
      }
     
//     Pion
      if (pt>=1.) {
        piAmpMin   = grpiFitAmpSmooth   ->GetY()[islice] - grpiFitAmpSmooth   ->GetY()[islice]*0.01;
        piAmpMax   = grpiFitAmpSmooth   ->GetY()[islice] + grpiFitAmpSmooth   ->GetY()[islice]*0.01;
        piSigmaMin = grpiFitSigmaSmooth ->GetY()[islice] - grpiFitSigmaSmooth ->GetY()[islice]*0.05;
        piSigmaMax = grpiFitSigmaSmooth ->GetY()[islice] + grpiFitSigmaSmooth ->GetY()[islice]*0.05;
      }
      
//     Kaon
      if (pt>=1.) {
        kaAmpMin   = grkaFitAmpSmooth   ->GetY()[islice] - grkaFitAmpSmooth   ->GetY()[islice]*0.01;
        kaAmpMax   = grkaFitAmpSmooth   ->GetY()[islice] + grkaFitAmpSmooth   ->GetY()[islice]*0.01;
        kaSigmaMin = grkaFitSigmaSmooth ->GetY()[islice] - grkaFitSigmaSmooth ->GetY()[islice]*0.05;
        kaSigmaMax = grkaFitSigmaSmooth ->GetY()[islice] + grkaFitSigmaSmooth ->GetY()[islice]*0.05;
      }
     
//     Proton
      if (pt>1.) {
        prAmpMin   = grprFitAmpSmooth   ->GetY()[islice] - grprFitAmpSmooth   ->GetY()[islice]*0.05;
        prAmpMax   = grprFitAmpSmooth   ->GetY()[islice] + grprFitAmpSmooth   ->GetY()[islice]*0.05;
        prSigmaMin = grprFitSigmaSmooth ->GetY()[islice] - grprFitSigmaSmooth ->GetY()[islice]*0.05;
        prSigmaMax = grprFitSigmaSmooth ->GetY()[islice] + grprFitSigmaSmooth ->GetY()[islice]*0.05;
      }
     
    }

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    
    
// !!!!!!!!!!!!!!!! iterations for Systematic uncertainty estimation !!!!!!!!!!!!!!!!!!!!!!!
    if (iIter>=7 && !MC) {   
      
      
      
      
      
      
      
      
      
      
      
      
      
    }
    
   
  }  // close settings for 10MeV bin width 
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

  //   Almost no restriction in low momentum
  if (pt<0.3 && !MC){
    
    Double_t tmpmaxBin  = h1D->GetBinContent(h1D->GetMaximumBin());
    
    elAmpMin = 1e-4;
    piAmpMin = tmpmaxBin/1.5;
    kaAmpMin = 1e-4;
    prAmpMin = 1e-4;
    deAmpMin = 1e-4;
    
    elAmpMax = tmpmaxBin/2.;
    piAmpMax = tmpmaxBin*1.1;
    kaAmpMax = tmpmaxBin;
    prAmpMax = tmpmaxBin;
    deAmpMax = tmpmaxBin;
    
    elMeanMin = TMath::Max(30.,(arrMean[0]-4.*arrSigma[0]));            
    piMeanMin = TMath::Max(30.,(arrMean[1]-4.*arrSigma[1]));            
    kaMeanMin = TMath::Max(30.,(arrMean[2]-4.*arrSigma[2]));                              
        
    elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+4.*arrSigma[0]));
    piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+4.*arrSigma[1]));
    kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+4.*arrSigma[2]));
    
    // pr and de is special for low momentum
    prMeanMin = TMath::Max(30.,(arrMean[3]-4.*arrSigma[3]));                
    deMeanMin = TMath::Max(30.,(arrMean[4]-4.*arrSigma[4]));  
    prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+4.*arrSigma[3]));
    deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[4]+4.*arrSigma[4]));
    
    elSigmaMin = arrSigma[0]/5.;
    piSigmaMin = arrSigma[1]/5.;
    kaSigmaMin = arrSigma[2]/5.;
    prSigmaMin = arrSigma[3]/5.;
    deSigmaMin = arrSigma[4]/5.;
    
    elSigmaMax = arrSigma[0]*3.;
    piSigmaMax = arrSigma[1]*3.;
    kaSigmaMax = arrSigma[2]*3.;
    prSigmaMax = arrSigma[3]*3.;
    deSigmaMax = arrSigma[4]*3.; 
    
    // Muons are special
    if (iIter>=1){
      muSigmaMin = grmuSigmaExScaled->GetY()[islice]-grmuSigmaExScaled->GetY()[islice]*0.05;
      muSigmaMax = grmuSigmaExScaled->GetY()[islice]+grmuSigmaExScaled->GetY()[islice]*0.05;
      muMeanMin  = grmuMeanExScaled ->GetY()[islice]-grmuMeanExScaled ->GetY()[islice]*0.001;
      muMeanMax  = grmuMeanExScaled ->GetY()[islice]+grmuMeanExScaled ->GetY()[islice]*0.001;
      muAmpMin   = muAmpMCscaled-muAmpMCscaled*0.05;
      muAmpMax   = muAmpMCscaled+muAmpMCscaled*0.05;
    }
    
  }
  
  // ???????
  muAmpMin   = 1e-4;
  muAmpMax   = 1e-3;
  
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
//   Finally set the params
  
  // set the Kurtosis and Skewness as free, fix them if needed (incase modify fixedK and fixedS)
  total->SetParLimits(4,skewMin,skewMax);
  total->SetParLimits(9,skewMin,skewMax);
  total->SetParLimits(14,skewMin,skewMax);
  total->SetParLimits(19,skewMin,skewMax);
  total->SetParLimits(24,skewMin,skewMax);
  total->SetParLimits(29,skewMin,skewMax);
    
  total->SetParLimits(3,kurtosisMin,kurtosisMax);
  total->SetParLimits(8,kurtosisMin,kurtosisMax);
  total->SetParLimits(13,kurtosisMin,kurtosisMax);
  total->SetParLimits(18,kurtosisMin,kurtosisMax);
  total->SetParLimits(23,kurtosisMin,kurtosisMax);
  total->SetParLimits(28,kurtosisMin,kurtosisMax);
  
  // if needed fix Kurtosis and Skewness
  if ( (fixedK || fixedS) && iIter<=freedomForPion ) SetFixedKSparameterForTotalFit(total);
  
  // smaller than 0.6 GeV/c pions should be fitted with free KS
  if (iIter>=freedomForPion && pt<0.6) {   
      
    SetFixedKSparameterForTotalFitForPions(total);
        
    piAmpMin   = h1D->GetBinContent(h1D->GetMaximumBin())*0.7;
    piAmpMax   = h1D->GetBinContent(h1D->GetMaximumBin())*1.2;
        
    piMeanMin = TMath::Min((Double_t)dEdxMax,(arrMean[1] - 4.*arrSigma[1]));
    piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1] + 4.*arrSigma[1]));
    
    piSigmaMin = arrSigma[1]/5.;
    piSigmaMax = arrSigma[1]*3.;
           
  }
  
  //if the momentum larger than a thresholold fix KS and apply different strategy
  if (islice>=nSlice && !MC) {
    
    // use pid response for mean and sigma
    if (iIter>=1 && MC) {
      arrMean[0] = grelMeanExScaled->GetY()[islice]-4. ; arrSigma[0] = grelSigmaExScaled->GetY()[islice]+1.;
      arrMean[1] = grpiMeanExScaled->GetY()[islice]+1.5; arrSigma[0] = grpiSigmaExScaled->GetY()[islice]-1.;
      arrMean[2] = grkaMeanExScaled->GetY()[islice]+5. ; arrSigma[0] = grkaSigmaExScaled->GetY()[islice]+1.;
      arrMean[3] = grprMeanExScaled->GetY()[islice]+5. ; arrSigma[0] = grprSigmaExScaled->GetY()[islice]+0.8;
      
      arrMeanWindow[0] = GetClosestParticleMean(islice,arrSigma[0],arrMean[0],arrMean[1],arrMean[2],arrMean[3]);
      arrMeanWindow[1] = GetClosestParticleMean(islice,arrSigma[1],arrMean[1],arrMean[0],arrMean[2],arrMean[3]);
      arrMeanWindow[2] = GetClosestParticleMean(islice,arrSigma[2],arrMean[2],arrMean[1],arrMean[0],arrMean[3]);
      arrMeanWindow[3] = GetClosestParticleMean(islice,arrSigma[3],arrMean[3],arrMean[1],arrMean[2],arrMean[0]);  
    } else {
      arrMeanWindow[0] = GetClosestParticleMean(islice,arrSigma[0],arrMean[0],arrMean[1],arrMean[2],arrMean[3]);
      arrMeanWindow[1] = GetClosestParticleMean(islice,arrSigma[1],arrMean[1],arrMean[0],arrMean[2],arrMean[3]);
      arrMeanWindow[2] = GetClosestParticleMean(islice,arrSigma[2],arrMean[2],arrMean[1],arrMean[0],arrMean[3]);
      arrMeanWindow[3] = GetClosestParticleMean(islice,arrSigma[3],arrMean[3],arrMean[1],arrMean[2],arrMean[0]); 
    }
    
    SetHighPFitParams(h1D,arrMean,arrSigma,arrMeanWindow);
    
    // use symmetric gauss for the fit
    if (fixedK || fixedS){
      // Set skewness to 0
      total->FixParameter(4,highPskewness); 
      total->FixParameter(9,highPskewness); 
      total->FixParameter(14,highPskewness); 
      total->FixParameter(19,highPskewness); 
      total->FixParameter(24,highPskewness); 
      total->FixParameter(29,highPskewness); 
      // set kurtosis to ~1.95
      total->FixParameter(3,highPkurtosis); 
      total->FixParameter(8,highPkurtosis); 
      total->FixParameter(13,highPkurtosis); 
      total->FixParameter(18,highPkurtosis); 
      total->FixParameter(23,highPkurtosis); 
      total->FixParameter(28,highPkurtosis);
    }
    
    // Apply smooting to amplitudes only
    if (iIter>=1){
      
      piAmpMin   = TMath::Max(1e-4,grpiFitAmpSmoothHP->GetY()[islice-nSlice] - (grpiFitAmpSmoothHP->GetY()[islice-nSlice])*0.05);
      piAmpMax   = TMath::Max(1e-3,grpiFitAmpSmoothHP->GetY()[islice-nSlice] + (grpiFitAmpSmoothHP->GetY()[islice-nSlice])*0.05);
      piMeanMin   = TMath::Max(1e-4,grpiFitMeanSmoothHP->GetY()[islice-nSlice] - (grpiFitMeanSmoothHP->GetY()[islice-nSlice])*0.01);
      piMeanMax   = TMath::Max(1e-3,grpiFitMeanSmoothHP->GetY()[islice-nSlice] + (grpiFitMeanSmoothHP->GetY()[islice-nSlice])*0.01);
      
      kaAmpMin   = TMath::Max(1e-4,grkaFitAmpSmoothHP->GetY()[islice-nSlice] - (grkaFitAmpSmoothHP->GetY()[islice-nSlice])*0.02);
      kaAmpMax   = TMath::Max(1e-3,grkaFitAmpSmoothHP->GetY()[islice-nSlice] + (grkaFitAmpSmoothHP->GetY()[islice-nSlice])*0.02);
      kaMeanMin   = TMath::Max(1e-4,grkaFitMeanSmoothHP->GetY()[islice-nSlice] - (grkaFitMeanSmoothHP->GetY()[islice-nSlice])*0.01);
      kaMeanMax   = TMath::Max(1e-3,grkaFitMeanSmoothHP->GetY()[islice-nSlice] + (grkaFitMeanSmoothHP->GetY()[islice-nSlice])*0.01);
      
      prAmpMin   = TMath::Max(1e-4,grprFitAmpSmoothHP->GetY()[islice-nSlice] - (grprFitAmpSmoothHP->GetY()[islice-nSlice])*0.02);
      prAmpMax   = TMath::Max(1e-3,grprFitAmpSmoothHP->GetY()[islice-nSlice] + (grprFitAmpSmoothHP->GetY()[islice-nSlice])*0.02);
      prMeanMin   = TMath::Max(1e-4,grprFitMeanSmoothHP->GetY()[islice-nSlice] - (grprFitMeanSmoothHP->GetY()[islice-nSlice])*0.01);
      prMeanMax   = TMath::Max(1e-3,grprFitMeanSmoothHP->GetY()[islice-nSlice] + (grprFitMeanSmoothHP->GetY()[islice-nSlice])*0.01);
    
    }
    
  }  
  
  // after some point do not fit muon and deuteron
  if (pt>1.5) { deAmpMin = 1e-4; deAmpMax = 1e-3; }
  if (pt>0.6) { muAmpMin = 1e-4; muAmpMax = 1e-3; }
  
  // Set Final total fit params
  total->SetParLimits(0,elAmpMin,elAmpMax);
  total->SetParLimits(1,elMeanMin,elMeanMax);
  total->SetParLimits(2,elSigmaMin,elSigmaMax);
        
  total->SetParLimits(5,piAmpMin,piAmpMax);
  total->SetParLimits(6,piMeanMin,piMeanMax);
  total->SetParLimits(7,piSigmaMin,piSigmaMax);
    
  total->SetParLimits(10,kaAmpMin,kaAmpMax);
  total->SetParLimits(11,kaMeanMin,kaMeanMax);
  total->SetParLimits(12,kaSigmaMin,kaSigmaMax);
        
  total->SetParLimits(15,prAmpMin,prAmpMax);
  total->SetParLimits(16,prMeanMin,prMeanMax);
  total->SetParLimits(17,prSigmaMin,prSigmaMax);
        
  total->SetParLimits(20,deAmpMin,deAmpMax);
  total->SetParLimits(21,deMeanMin,deMeanMax);
  total->SetParLimits(22,deSigmaMin,deSigmaMax);
  
  total->SetParLimits(25,muAmpMin,muAmpMax);
  total->SetParLimits(26,muMeanMin,muMeanMax);
  total->SetParLimits(27,muSigmaMin,muSigmaMax);
 
}
// -------------------------------------------------------------------------------------------------------
void  SetHighPFitParams(TH1D *h1D, Double_t *arrMean, Double_t *arrSigma,Double_t *arrMeanWindow){
  
  //   
  //    Set mean amp and sigma par limits for the first iteration
  //   
  
//   Part to make setting for fit params 
  maxBin  = h1D->GetBinContent(h1D->GetMaximumBin());
  
    
  // Restrictions on the total fit 
  elAmpMin = 1e-4;
  elAmpMax = maxBin;
  piAmpMin = 1e-4;
  piAmpMax = maxBin;
  kaAmpMin = 1e-4;
  kaAmpMax = maxBin;
  prAmpMin = 1e-4;
  prAmpMax = maxBin;
  
  elMeanMin = TMath::Max(30.,(arrMean[0]-arrMeanWindow[0]));            
  piMeanMin = TMath::Max(30.,(arrMean[1]-arrMeanWindow[1]));            
  kaMeanMin = TMath::Max(30.,(arrMean[2]-arrMeanWindow[2]));                
  prMeanMin = TMath::Max(30.,(arrMean[3]-arrMeanWindow[3]));                
  
        
  elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+arrMeanWindow[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+arrMeanWindow[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+arrMeanWindow[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+arrMeanWindow[3]));
  
    
  elSigmaMin = arrSigma[0]*0.5;
  piSigmaMin = arrSigma[1]*0.5;      
  kaSigmaMin = arrSigma[2]*0.5;
  prSigmaMin = arrSigma[3]*0.5;
  
    
  elSigmaMax = arrSigma[0]*1.5;
  piSigmaMax = arrSigma[1]*1.5;
  kaSigmaMax = arrSigma[2]*1.5;
  prSigmaMax = arrSigma[3]*1.5;
  
  // do not fit deuteron and muon
  deMeanMin = TMath::Max(30.,(arrMean[4]-arrSigma[4])); 
  deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[4]+arrSigma[4])); 
  deSigmaMin = arrSigma[4]/0.8;
  deSigmaMax = arrSigma[4]*1.2; 
  deAmpMin = 1e-4;
  deAmpMax = 1e-3;
  // Muons are special 
  muSigmaMin = arrSigma[5]/5.;
  muSigmaMax = arrSigma[5]*2.;
  muMeanMin  = arrMean[5]-arrSigma[5];
  muMeanMax  = arrMean[5]+arrSigma[5];
  muAmpMin = 1e-4;
  muAmpMax = 1e-3;
}

void  SetParamsForKSEstimate(Int_t iks, TF1 *total, TH1D *h1DKS, Double_t *arrMean, Double_t *arrSigma){
  
  //   
  //    Set mean amp and sigma par limits for KS estimate
  //   
  
  
//   Part to make setting for fit params 
  maxBin  = h1DKS->GetBinContent(h1DKS->GetMaximumBin());
  maxBin0 = TMath::Min(h1DKS->GetBinContent(h1DKS->FindBin(arrMean[0])),maxBin);
  maxBin1 = TMath::Min(h1DKS->GetBinContent(h1DKS->FindBin(arrMean[1])),maxBin);
  maxBin2 = TMath::Min(h1DKS->GetBinContent(h1DKS->FindBin(arrMean[2])),maxBin);
  maxBin3 = TMath::Min(h1DKS->GetBinContent(h1DKS->FindBin(arrMean[3])),maxBin);
  
  Double_t fitWinMaxBin      = 0.2;                      // only used for 0th iteration   // TO FIX
    
  // Restrictions on the total fit 
  elAmpMin = 1e-4;
  elAmpMax = maxBin/10.;
  piAmpMin = 1e-4;
  piAmpMax = maxBin1 + maxBin1*fitWinMaxBin;
  kaAmpMin = 1e-4;
  kaAmpMax = maxBin1 + maxBin1*fitWinMaxBin;
  prAmpMin = 1e-4;
  prAmpMax = maxBin1 + maxBin1*fitWinMaxBin;
  deAmpMin = 0.;
  deAmpMax = maxBin0 + maxBin0*fitWinMaxBin;
  
  elMeanMin = TMath::Max(30.,(arrMean[0]-2.*arrSigma[0]));            
  piMeanMin = TMath::Max(30.,(arrMean[1]-2.*arrSigma[1]));            
  kaMeanMin = TMath::Max(30.,(arrMean[2]-2.*arrSigma[2]));                
  prMeanMin = TMath::Max(30.,(arrMean[3]-2.*arrSigma[3]));                
  deMeanMin = TMath::Max(30.,(arrMean[4]-2.*arrSigma[4]));  
        
  elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+2.*arrSigma[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+2.*arrSigma[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+2.*arrSigma[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+2.*arrSigma[3]));
  deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[4]+2.*arrSigma[4]));
        
  elSigmaMin = arrSigma[0];
  piSigmaMin = arrSigma[1];      
  kaSigmaMin = arrSigma[2];
  prSigmaMin = arrSigma[3];
  deSigmaMin = arrSigma[4];
    
  elSigmaMax = arrSigma[0]*2.;
  piSigmaMax = arrSigma[1]*2.;
  kaSigmaMax = arrSigma[2]*2.;
  prSigmaMax = arrSigma[3]*2.;
  deSigmaMax = arrSigma[4]*2.; 
  
  // Muons are special 
  muSigmaMin = arrSigma[5];
  muSigmaMax = arrSigma[5]*2.;
  muMeanMin  = arrMean[5]-arrSigma[5];
  muMeanMax  = arrMean[5]+arrSigma[5];
  muAmpMin = 1e-4;
  muAmpMax = 1e-3;
  
  // Set total fit Params
  total->SetParLimits(0,elAmpMin,elAmpMax);
  total->SetParLimits(1,elMeanMin,elMeanMax);
  total->SetParLimits(2,elSigmaMin,elSigmaMax);
  (iks==0) ? total->SetParLimits(3,kurtosisMin,kurtosisMax) : total->FixParameter(3,kurtosisFix[kElectron]);
  total->SetParLimits(4,skewMin,skewMax);

  total->SetParLimits(5,piAmpMin,piAmpMax);
  total->SetParLimits(6,piMeanMin,piMeanMax);
  total->SetParLimits(7,piSigmaMin,piSigmaMax);
  (iks==0) ? total->SetParLimits(8,kurtosisMin,kurtosisMax) : total->FixParameter(8,kurtosisFix[kPion]);
  total->SetParLimits(9,skewMin,skewMax);
    
  total->SetParLimits(10,kaAmpMin,kaAmpMax);
  total->SetParLimits(11,kaMeanMin,kaMeanMax);
  total->SetParLimits(12,kaSigmaMin,kaSigmaMax);
  (iks==0) ? total->SetParLimits(13,kurtosisMin,kurtosisMax) : total->FixParameter(13,kurtosisFix[kKaon]);
  total->SetParLimits(14,skewMin,skewMax);
        
  total->SetParLimits(15,prAmpMin,prAmpMax);
  total->SetParLimits(16,prMeanMin,prMeanMax);
  total->SetParLimits(17,prSigmaMin,prSigmaMax);
  (iks==0) ? total->SetParLimits(18,kurtosisMin,kurtosisMax) : total->FixParameter(18,kurtosisFix[kProton]);
  total->SetParLimits(19,skewMin,skewMax);
        
  total->SetParLimits(20,deAmpMin,deAmpMax);
  total->SetParLimits(21,deMeanMin,deMeanMax);
  total->SetParLimits(22,deSigmaMin,deSigmaMax);
  (iks==0) ? total->SetParLimits(23,kurtosisMin,kurtosisMax) : total->FixParameter(23,kurtosisFix[kDeuteron]);
  total->SetParLimits(24,skewMin,skewMax);
  
  total->SetParLimits(25,muAmpMin,muAmpMax);
  total->SetParLimits(26,muMeanMin,muMeanMax);
  total->SetParLimits(27,muSigmaMin,muSigmaMax);
  (iks==0) ? total->SetParLimits(28,kurtosisMin,kurtosisMax) : total->FixParameter(28,kurtosisFix[kMuon]);
  total->SetParLimits(29,skewMin,skewMax);
  
}

// -------------------------------------------------------------------------------------------------------
void  SetParams1stIteration(Double_t pt, TH1D *h1D, Double_t *arrMean, Double_t *arrSigma,Double_t *arrMeanWindow, Double_t *cleanParMean){
  
  //   
    //    Set mean amp and sigma par limits for the first iteration
  //   
  
  Double_t meantmparr[5];
  for (Int_t i=0;i<5;i++) meantmparr[i] = arrMean[i];
  
  
  meantmparr[0] = ( (pt>0.35) && TMath::Abs(arrMean[0]-cleanParMean[0])>12 ) ?  arrMean[0] : cleanParMean[0]-1.5;
  meantmparr[1] = ( (pt>0.35) && TMath::Abs(arrMean[1]-cleanParMean[1])>12 ) ?  arrMean[1] : cleanParMean[1];
  meantmparr[2] = ( (pt>0.40) && TMath::Abs(arrMean[2]-cleanParMean[2])>12 ) ?  arrMean[2] : cleanParMean[2];
  meantmparr[3] = ( (pt>0.50) && TMath::Abs(arrMean[3]-cleanParMean[3])>12 ) ?  arrMean[3] : cleanParMean[3];
  
//   Part to make setting for fit params 
  maxBin  = h1D->GetBinContent(h1D->GetMaximumBin());
  maxBin0 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[0])),maxBin);
  maxBin1 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[1])),maxBin);
  maxBin2 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[2])),maxBin);
  maxBin3 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[3])),maxBin);
  
  Double_t fitWinMean        = 1.2;                      // used for every iteration
  Double_t fitWinMaxBin      = 0.2;                      // only used for 0th iteration   // TO FIX
    
  // Restrictions on the total fit 
  elAmpMin = 1e-4;
  elAmpMax = maxBin/10.;
  piAmpMin = 1e-4;
  piAmpMax = maxBin1 + maxBin1*fitWinMaxBin;
  kaAmpMin = 1e-4;
  kaAmpMax = maxBin1 + maxBin1*fitWinMaxBin;
  prAmpMin = 1e-4;
  prAmpMax = maxBin1 + maxBin1*fitWinMaxBin;
  deAmpMin = 0.;
  deAmpMax = maxBin0 + maxBin0*fitWinMaxBin;
  
  elMeanMin = TMath::Max(30.,(meantmparr[0]-arrMeanWindow[0]));            
  piMeanMin = TMath::Max(30.,(meantmparr[1]-arrMeanWindow[1]));            
  kaMeanMin = TMath::Max(30.,(meantmparr[2]-arrMeanWindow[2]));                
  prMeanMin = TMath::Max(30.,(meantmparr[3]-arrMeanWindow[3]));                
  deMeanMin = TMath::Max(30.,(meantmparr[4]-fitWinMean*arrSigma[4]));  
        
  elMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[0]+arrMeanWindow[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[1]+arrMeanWindow[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[2]+arrMeanWindow[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[3]+arrMeanWindow[3]));
  deMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[4]+fitWinMean*arrSigma[4]));
        
  elSigmaMin = (Systematic) ? arrSigma[0]/2. : arrSigma[0]+2;
  piSigmaMin = (Systematic) ? arrSigma[1]/2. : arrSigma[1]+2;      
  kaSigmaMin = (Systematic) ? arrSigma[2]/2. : arrSigma[2]+1;
  prSigmaMin = (Systematic) ? arrSigma[3]/2. : arrSigma[3]+1;
  deSigmaMin = (Systematic) ? arrSigma[4]/2. : arrSigma[4]+2;
    
  elSigmaMax = arrSigma[0]*2.;
  piSigmaMax = arrSigma[1]*2.;
  kaSigmaMax = arrSigma[2]*2.;
  prSigmaMax = arrSigma[3]*2.;
  deSigmaMax = arrSigma[4]*2.; 
  
  // Muons are special 
  muSigmaMin = arrSigma[5];
  muSigmaMax = arrSigma[5]*2.;
  muMeanMin  = arrMean[5]-arrSigma[5];
  muMeanMax  = arrMean[5]+arrSigma[5];
  muAmpMin = 1e-4;
  muAmpMax = 1e-3;
}
// -------------------------------------------------------------------------------------------------------
void  SetParams1stIterationMC(Double_t pt, TH1D *h1D, Double_t *arrMean, Double_t *arrSigma, Double_t * MCPars){
  
  //   
  //    Set mean amp and sigma par limits for the first iteration
  //   
  
  //   Part to make setting for fit params 
  maxBin  = h1D->GetBinContent(h1D->GetMaximumBin());
  maxBin0 = TMath::Min(h1D->GetBinContent(h1D->FindBin(MCPars[0])),maxBin);
  maxBin1 = TMath::Min(h1D->GetBinContent(h1D->FindBin(MCPars[1])),maxBin);
  maxBin2 = TMath::Min(h1D->GetBinContent(h1D->FindBin(MCPars[2])),maxBin);
  maxBin3 = TMath::Min(h1D->GetBinContent(h1D->FindBin(MCPars[3])),maxBin);
  
  Double_t fitWinMean        = (pt<0.4) ? 2. : 1.2;                      // used for every iteration
    
  // Restrictions on the total fit 
  elAmpMin = 1e-4;
  elAmpMax = maxBin/10.;
  piAmpMin = 1e-4;
  piAmpMax = maxBin + maxBin*0.2;
  kaAmpMin = 1e-4;
  kaAmpMax = maxBin;
  prAmpMin = 1e-4;
  prAmpMax = maxBin;
  
  
  elMeanMin = TMath::Max(30.,(MCPars[0]-fitWinMean*MCPars[5]));            
  piMeanMin = TMath::Max(30.,(MCPars[1]-fitWinMean*MCPars[6]));            
  kaMeanMin = TMath::Max(30.,(MCPars[2]-fitWinMean*MCPars[7]));                
  prMeanMin = TMath::Max(30.,(MCPars[3]-fitWinMean*MCPars[8]));                
    
        
  elMeanMax = TMath::Min((Double_t)dEdxMax,(MCPars[0]+fitWinMean*MCPars[5]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(MCPars[1]+fitWinMean*MCPars[6]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(MCPars[2]+fitWinMean*MCPars[7]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(MCPars[3]+fitWinMean*MCPars[8]));
  
        
  elSigmaMin = MCPars[5]/2.;
  piSigmaMin = MCPars[6]/2.;      
  kaSigmaMin = MCPars[7]/2.;
  prSigmaMin = MCPars[8]/2.;
  
    
  elSigmaMax = MCPars[5]*2.;
  piSigmaMax = MCPars[6]*2.;
  kaSigmaMax = MCPars[7]*2.;
  prSigmaMax = MCPars[8]*2.;
  
  TH1D *htmppr = (TH1D*)h1D->Clone();
  TH1D *htmpka = (TH1D*)h1D->Clone();
  htmppr->GetXaxis()->SetRangeUser(250,800);
  htmpka->GetXaxis()->SetRangeUser(100,250);

  if (pt<0.3) {
    kaMeanMin = (htmpka->GetMean())-arrSigma[3];                
    kaMeanMax = (htmpka->GetMean())+arrSigma[3];
    kaSigmaMin = arrSigma[2]/3.;
    kaSigmaMax = arrSigma[2]*3.;
    
    prMeanMin = (htmppr->GetMean())-2.*arrSigma[3];                
    prMeanMax = (htmppr->GetMean())+2.*arrSigma[3];
    prSigmaMin = arrSigma[3]/3.;
    prSigmaMax = arrSigma[3]*3.;
  }
  delete htmppr;
  delete htmpka;
  
  // No deuterons
  deSigmaMin = MCPars[8]/2.;
  deSigmaMax = MCPars[8]*2.; 
  deMeanMin = TMath::Max(30.,(MCPars[3]-MCPars[8]));
  deMeanMax = TMath::Min((Double_t)dEdxMax,(MCPars[3]+MCPars[8]));
  deAmpMin = 1e-4;
  deAmpMax = 1e-3;
  
  // No muons 
  muSigmaMin = MCPars[6];
  muSigmaMax = MCPars[6]*2.;
  muMeanMin  = MCPars[1]-MCPars[6];
  muMeanMax  = MCPars[1]+MCPars[6];
  muAmpMin = 1e-4;
  muAmpMax = 1e-3;
  
}
// -------------------------------------------------------------------------------------------------------
void SetParamsMC(Int_t islice, Double_t pt){
  
  //   
//  Set Params for the iterative of the h2DMC
  //   
  
  Double_t windowRange = 0.001;
   
//   ***************  Electron Patches ***************
  if (pt>0.3) {  
    elSigmaMin = grelMCSigmaScaled->GetY()[islice]-grelMCSigmaScaled->GetY()[islice]*0.01;
    elSigmaMax = grelMCSigmaScaled->GetY()[islice]+grelMCSigmaScaled->GetY()[islice]*0.01;
    elMeanMin  = grelMCMeanScaled->GetY()[islice] -grelMCMeanScaled->GetY()[islice]*windowRange;
    elMeanMax  = grelMCMeanScaled->GetY()[islice] +grelMCMeanScaled->GetY()[islice]*windowRange;
    elAmpMin   = grelMCAmpScaled->GetY()[islice]  -(grelMCAmpScaled->GetY()[islice])*0.05;
    elAmpMax   = grelMCAmpScaled->GetY()[islice]  +(grelMCAmpScaled->GetY()[islice])*0.05;  
  }
    
  if (pt>0.7) {
    elSigmaMin = grelSigmaExScaled->GetY()[islice]-grelSigmaExScaled->GetY()[islice]*0.1;
    elSigmaMax = grelSigmaExScaled->GetY()[islice]+grelSigmaExScaled->GetY()[islice]*0.1;
    elMeanMin  = grelMeanExScaled->GetY()[islice] -grelMeanExScaled->GetY()[islice]*0.05;
    elMeanMax  = grelMeanExScaled->GetY()[islice] +grelMeanExScaled->GetY()[islice]*0.05;
  }
   
//   ***************  Pion Patches ***************  
  if ((pt>=0.6)){   // pion mean patch
    piMeanMin = grpiMCMeanScaled->GetY()[islice]-(grpiMCMeanScaled->GetY()[islice])*windowRange;
    piMeanMax = grpiMCMeanScaled->GetY()[islice]+(grpiMCMeanScaled->GetY()[islice])*windowRange;
    piSigmaMin = grpiMCSigmaScaled->GetY()[islice]-(grpiMCSigmaScaled->GetY()[islice])*0.03;
    piSigmaMax = grpiMCSigmaScaled->GetY()[islice]+(grpiMCSigmaScaled->GetY()[islice])*0.03;  
    piAmpMin = grpiMCAmpScaled->GetY()[islice]-(grpiMCAmpScaled->GetY()[islice])*0.005;
    piAmpMax = grpiMCAmpScaled->GetY()[islice]+(grpiMCAmpScaled->GetY()[islice])*0.005;  
  }
  
//   ***************  kaon  patches ***************  
  if (pt>=0.35) {
    kaSigmaMin = grkaMCSigmaScaled->GetY()[islice]-(grkaMCSigmaScaled->GetY()[islice])*windowRange;
    kaSigmaMax = grkaMCSigmaScaled->GetY()[islice]+(grkaMCSigmaScaled->GetY()[islice])*windowRange;
    kaMeanMin = grkaMCMeanScaled->GetY()[islice]-(grkaMCMeanScaled->GetY()[islice])*windowRange;
    kaMeanMax = grkaMCMeanScaled->GetY()[islice]+(grkaMCMeanScaled->GetY()[islice])*windowRange;
  } 
    
// ***************  Proton Patches ***************      
  if ((pt>=0.5)){   // proton sigma patch
    prSigmaMin = grprMCSigmaScaled->GetY()[islice]-(grprMCSigmaScaled->GetY()[islice])*0.01;
    prSigmaMax = grprMCSigmaScaled->GetY()[islice]+(grprMCSigmaScaled->GetY()[islice])*0.01;
    prMeanMin = grprMCMeanScaled->GetY()[islice]-(grprMCMeanScaled->GetY()[islice])*windowRange;
    prMeanMax = grprMCMeanScaled->GetY()[islice]+(grprMCMeanScaled->GetY()[islice])*windowRange;
  }
   
//   ***************  Muon Patches ***************  
  muSigmaMin = grmuSigmaExScaled->GetY()[islice]-(grmuSigmaExScaled->GetY()[islice])*0.01;
  muSigmaMax = grmuSigmaExScaled->GetY()[islice]+(grmuSigmaExScaled->GetY()[islice])*0.01;
  muMeanMin  = grmuMeanExScaled->GetY()[islice]-(grmuMeanExScaled->GetY()[islice])*windowRange;
  muMeanMax  = grmuMeanExScaled->GetY()[islice]+(grmuMeanExScaled->GetY()[islice])*windowRange;
  muAmpMin   = 1e-4;
  muAmpMax   = 1e-3;
      
}      
// -------------------------------------------------------------------------------------------------------
TF1  *FitParamInRange(TGraphErrors* gr, TString funcType, Double_t min, Double_t max){

  //
  // Fit proton TGraph with pol2
  //
  
  TF1 *f = new TF1("f",funcType,min,max);
  gr->Fit(f,"QNR+");
  f->SetRange(0.,ptRangeUp);
  return f;

}
// -------------------------------------------------------------------------------------------------------
TH1D *FuncToHist(TF1 *f, TH1D * h, Int_t islice, Double_t *arrMean, Double_t *arrSigma){

  // Convert function to hist

  Int_t n = h->GetNbinsX();
  TH1D *hFunc = new TH1D("hFunc","hFunc",n,dEdxMin,dEdxMax);
  hFunc->SetName(Form("hFunc_%d",islice));

  for (Int_t i = 0; i<n; i++){
    hFunc->SetBinContent(i,f->Eval(h->GetBinCenter(i)));
    hFunc->SetBinError(i,h->GetBinError(i));
  }
    
    // sigmas of el, pt, ka, pr
//   Double_t s1 = arrSigma[0];
  Double_t s2 = arrSigma[1];
  Double_t s3 = arrSigma[2];
  Double_t s4 = arrSigma[3];
    
    // means of el, pt, ka, pr
//   Double_t m1 = arrMean[0];
  Double_t m2 = arrMean[1];
  Double_t m3 = arrMean[2];
  Double_t m4 = arrMean[3];
  Double_t x; 
  Double_t sigmacut = 2.;
    
    // put a cut of 3*sigma around the fit functions
  for (Int_t i = 0; i<n; i++){
      
    x = h->GetBinCenter(i);
    if ( !(TMath::Abs(x-m2)<sigmacut*s2 || TMath::Abs(x-m3)<sigmacut*s3 || TMath::Abs(x-m4)<sigmacut*s4) )
    {
      hFunc->SetBinContent(i,0.);
      hFunc->SetBinError(i,0.);
      h    ->SetBinContent(i,0.);  
    }
  }
 
  return hFunc;
}
// -------------------------------------------------------------------------------------------------------
Double_t ComputeGaussIntegral(Double_t parAmp,Double_t parMean,Double_t parSigma, Double_t parKurtosis, Double_t parSkew){
 
  // Calculate integral of the gauss functions

  Double_t min = TMath::Max(0.,parMean-parSigma*5);
  Double_t max = TMath::Min((Double_t)dEdxMax,parMean+parSigma*5);
  
  TF1 f("f",fitFunctionGenGaus,min,max);
//   f.SetNpx(1000);
  f.SetParameter(0,parAmp);
  f.SetParameter(1,parMean);
  f.SetParameter(2,parSigma);
  f.SetParameter(3,parKurtosis);
  f.SetParameter(4,parSkew);
  
  return f.Integral(min,max);
  
}
// -------------------------------------------------------------------------------------------------------
TH1D *GetClean1DSlice(TH2D *h2Clean, TString parName, Int_t ptDownBin, Int_t ptUpBin){
  
  //
  //  get 1D momentum slice
  //
  
  TH1D *h1Clean   = (TH1D*)h2Clean->ProjectionY(Form(parName,ptDownBin),ptDownBin,ptUpBin);
  Double_t mean = h1Clean->GetMean();
  Double_t rms  = h1Clean->GetRMS();
  Double_t rangeDown = TMath::Max((Double_t)dEdxMin,mean-5.*rms);
  Double_t rangeUp   = TMath::Min((Double_t)dEdxMax,mean+5.*rms);
  h1Clean->GetXaxis()->SetRangeUser(rangeDown,rangeUp);
  
  return h1Clean;
}
// -------------------------------------------------------------------------------------------------------
void FitParticleSampleForKSEstimate(Int_t iks, TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy, Int_t sampleType){
 
  //
    // Fit Muon and Gauss in K0s results
    // sampleType = 0 --> CleanSamples, 
    // sampleType = 1 --> MCsamples
  //

  TVirtualFitter::SetMaxIterations(100000);
  
    // helper numbers for setparams
  maxBin                   = hClean->GetBinContent(hClean->GetMaximumBin());  if ((maxBin<5.)) return;   
  Double_t mean            = (sampleType == 0) ? arrMean[pSpecy] : hClean->GetMean();
  Double_t rms             = (sampleType == 0) ? arrSigma[pSpecy]: hClean->GetRMS();
  Double_t sampleFitWindow = (sampleType == 0) ? 3.5 : 4.2; 
        
    // fit function
  TF1 *asyGaus = new TF1("asyGaus",fitFunctionGenGaus,mean-sampleFitWindow*rms,mean+4.2*rms);   
  asyGaus->SetParLimits(0,maxBin/20.,maxBin);
  asyGaus->SetParLimits(1,mean-rms,mean+rms);
  asyGaus->SetParLimits(2,rms/2. ,3.*rms);
  if (iks==0){
    asyGaus->SetParLimits(3,kurtosisMin,kurtosisMax);
    asyGaus->SetParLimits(4,skewMin,skewMax);
  } else {
    (sampleType==0) ? asyGaus->FixParameter(3,kurtosisFixClean[pSpecy]): asyGaus->FixParameter(3,kurtosisFixMC[pSpecy]);
    asyGaus->SetParLimits(4,skewMin,skewMax);
  }
    
  // retrieve fit parameters
  hClean->Fit(asyGaus,"QNMR");
  parSample[0] = (asyGaus->GetParameter(0)>1.) ? (asyGaus->GetParameter(0)) : 1.;      // avoid floating point exception for muons
  parSample[1] = asyGaus->GetParameter(1);
  parSample[2] = asyGaus->GetParameter(2);
  parSample[3] = asyGaus->GetParameter(3);
  parSample[4] = asyGaus->GetParameter(4);
  parSample[5] = asyGaus->Integral(parSample[1]-parSample[2]*4.,parSample[1]+parSample[2]*4.);
  if (asyGaus->GetNDF()>1e-4) parSample[5] =  asyGaus->GetChisquare()/asyGaus->GetNDF();
    
  delete asyGaus;
}
// -------------------------------------------------------------------------------------------------------
void Fit1DSliceForKSEstimate(Int_t iks, TH1D *h1DKS, Double_t *pars, Double_t *arrMean, Double_t *arrSigma){
  
    // Fit part
  maxBin = h1DKS->GetBinContent(h1DKS->GetMaximumBin());   
  Double_t fitWinMean  = 1.2;
    
  elMeanMin = TMath::Max(30.,(arrMean[0]-2.*arrSigma[0]));                                   // electron
  piMeanMin = TMath::Max(30.,(arrMean[1]-fitWinMean*arrSigma[1]));                           // pion
  kaMeanMin = TMath::Max(30.,(arrMean[2]-fitWinMean*arrSigma[2]));                           // kaon
  prMeanMin = TMath::Max(30.,(arrMean[3]-fitWinMean*arrSigma[3]));                           // proton
  deMeanMin = TMath::Max(30.,(arrMean[4]-fitWinMean*arrSigma[4]));                           // deuteron
  muMeanMin = TMath::Max(30.,(arrMean[5]-fitWinMean*arrSigma[5]));                           // deuteron
        
  elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+2.*arrSigma[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+fitWinMean*arrSigma[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+fitWinMean*arrSigma[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+fitWinMean*arrSigma[3]));
  deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[4]+fitWinMean*arrSigma[4]));
  muMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[5]+fitWinMean*arrSigma[5]));
    
   
  TF1 *g1    = new TF1("g1",fitFunctionGenGaus,elMeanMin,elMeanMax);        g1->SetNpx(1000);
  TF1 *g2    = new TF1("g2",fitFunctionGenGaus,piMeanMin,piMeanMax);        g2->SetNpx(1000);
  TF1 *g3    = new TF1("g3",fitFunctionGenGaus,kaMeanMin,kaMeanMax);        g3->SetNpx(1000);
  TF1 *g4    = new TF1("g4",fitFunctionGenGaus,prMeanMin,prMeanMax);        g4->SetNpx(1000);
  TF1 *g5    = new TF1("g5",fitFunctionGenGaus,deMeanMin,deMeanMax);        g5->SetNpx(1000);
  TF1 *g6    = new TF1("g6",fitFunctionGenGaus,muMeanMin,muMeanMax);        g6->SetNpx(1000);
            
    // restrict the parameters of gaus functions of individual particles
  g1->SetParLimits(0,0.,maxBin);
  g2->SetParLimits(0,0.,maxBin);
  g3->SetParLimits(0,0.,maxBin);
  g4->SetParLimits(0,0.,maxBin);
  g5->SetParLimits(0,0.,maxBin);
  g6->SetParLimits(0,0.,maxBin);
        
  g1->SetParLimits(1,elMeanMin,elMeanMax);
  g2->SetParLimits(1,piMeanMin,piMeanMax);
  g3->SetParLimits(1,kaMeanMin,kaMeanMax);
  g4->SetParLimits(1,prMeanMin,prMeanMax);
  g5->SetParLimits(1,deMeanMin,deMeanMax);
  g6->SetParLimits(1,muMeanMin,muMeanMax);   
    
  g1->SetParLimits(2,arrSigma[0]/2.,arrSigma[0]*2.);
  g2->SetParLimits(2,arrSigma[1]/2.,arrSigma[1]*2.);
  g3->SetParLimits(2,arrSigma[2]/2.,arrSigma[2]*2.);
  g4->SetParLimits(2,arrSigma[3]/2.,arrSigma[3]*2.);
  g5->SetParLimits(2,arrSigma[4]/2.,arrSigma[4]*3.);
  g6->SetParLimits(2,arrSigma[5]/2.,arrSigma[5]*3.);
    
  (iks==0)? g1->SetParLimits(3,kurtosisMin,kurtosisMax) : g1->FixParameter(3,kurtosisFix[kElectron]);
  (iks==0)? g2->SetParLimits(3,kurtosisMin,kurtosisMax) : g2->FixParameter(3,kurtosisFix[kPion]    );
  (iks==0)? g3->SetParLimits(3,kurtosisMin,kurtosisMax) : g3->FixParameter(3,kurtosisFix[kKaon]    );
  (iks==0)? g4->SetParLimits(3,kurtosisMin,kurtosisMax) : g4->FixParameter(3,kurtosisFix[kProton]  );
  (iks==0)? g5->SetParLimits(3,kurtosisMin,kurtosisMax) : g5->FixParameter(3,kurtosisFix[kDeuteron]);
  (iks==0)? g6->SetParLimits(3,kurtosisMin,kurtosisMax) : g6->FixParameter(3,kurtosisFix[kMuon]    );

  g1->SetParLimits(4,skewMin,skewMax);
  g2->SetParLimits(4,skewMin,skewMax);
  g3->SetParLimits(4,skewMin,skewMax);
  g4->SetParLimits(4,skewMin,skewMax);
  g5->SetParLimits(4,skewMin,skewMax);
  g6->SetParLimits(4,skewMin,skewMax);
                
    // Fit h1D for individual particles
  if (arrMean[0]>=10) h1DKS->Fit(g1,"QNR");
  if (arrMean[1]>=10) h1DKS->Fit(g2,"QNR+");
  if (arrMean[2]>=10) h1DKS->Fit(g3,"QNR+");
  if (arrMean[3]>=10) h1DKS->Fit(g4,"QNR+");
  if (arrMean[4]>=10) h1DKS->Fit(g5,"QNR+");
  if (arrMean[5]>=10) h1DKS->Fit(g6,"QNR+");
    
    // Get Fit parameters from the individual fits and set for the total fit
  Double_t par[30] = {0};
  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[5]);
  g3->GetParameters(&par[10]);
  g4->GetParameters(&par[15]);
  g5->GetParameters(&par[20]);
  g5->GetParameters(&par[25]);
        
  TF1 *total = new TF1("total","g1+g2+g3+g4+g5+g6",dEdxMin,dEdxMax);
  total->SetNpx(1000);
  total->SetParameters(par);
            
    // Some setters for the total fit 
  SetParamsForKSEstimate(iks,total,h1DKS,arrMean,arrSigma);
  CheckIfMeanIsZero(total,arrMean);
    
    // Apply total fit
  h1DKS->Fit(total,"QMR+");
  Double_t totalChi2    = total->GetChisquare();
  Double_t normNDF      = total->GetNDF();
  Double_t totChi2      = (normNDF<1) ? 1 : totalChi2/normNDF;
    
  for (Int_t i=0;i<30;i++) pars[i] = total->GetParameter(i);
  pars[30] = totChi2;
  
  delete total;
  delete g1;
  delete g2;
  delete g3;
  delete g4;
  delete g5;
  delete g6;
} 
// -------------------------------------------------------------------------------------------------------
TMatrixD *FitParticleSample(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy, Int_t sampleType){
 
  //
    // Fit Muon and Gauss in K0s results
    // sampleType = 0 --> CleanSamples, 
    // sampleType = 1 --> MCsamples
  //

//   TVirtualFitter::SetDefaultFitter("Minuit"); 
  TVirtualFitter::SetMaxIterations(100000);
  
  maxBin = hClean->GetBinContent(hClean->GetMaximumBin());  
  if ((maxBin<5.)) {   // there is too few entries, fit would fail
    parSample[0] = 1.;
    parSample[1] = 1.;
    parSample[2] = 1.;
    parSample[3] = 1.;
    parSample[4] = 1.;
    parSample[5] = 1.;
    return new TMatrixD(5,5);   
  }
    
  Double_t mean            = (sampleType == 0) ? arrMean[pSpecy] : hClean->GetMean();
  Double_t rms             = (sampleType == 0) ? arrSigma[pSpecy]: hClean->GetRMS();
  Double_t sampleFitWindow = (sampleType == 0) ? 3.5 : 4.2; 
  Double_t integral        = hClean->Integral(hClean->FindBin(mean-sampleFitWindow*rms),hClean->FindBin(mean+4.2*rms));
    
  // prepare normalised samples
  if (normalisedCleanSamples && integral > 1e-4) {
    hClean->Scale(1./integral);
    maxBin /= integral;
  }
  
  // fit function
  TF1 *asyGaus = new TF1("asyGaus",fitFunctionGenGaus,mean-sampleFitWindow*rms,mean+4.2*rms);   
  asyGaus->SetParNames("Amplitude","Mean","Sigma","Kurtosis","Skewness");
//   asyGaus->SetNpx(1000); 
  asyGaus ->SetLineColor(2);
  asyGaus ->SetLineWidth(2);
  
  // Set Parameters
  asyGaus->SetParLimits(0,maxBin/20.,maxBin);
  asyGaus->SetParLimits(1,mean-rms,mean+rms);
  asyGaus->SetParLimits(2,rms/2. ,3.*rms);
  asyGaus->SetParLimits(3,kurtosisMin,kurtosisMax);
  asyGaus->SetParLimits(4,skewMin,skewMax);
  if ((fixedK || fixedS) && sampleType==0) SetFixedKSparameterForCleanSampleFit(asyGaus,pSpecy);
  if ((fixedK || fixedS) && sampleType==1) SetFixedKSparameterForMCSampleFit(asyGaus,pSpecy);

  TH1D *hMatrix = (TH1D*)hClean->Clone();
  TF1  *fMatrix = (TF1*)asyGaus->Clone();
  
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
  
  
    // Retrive the Covariance Matrix
  TFitResultPtr fitCheck = hMatrix->Fit(fMatrix,"SRQ");
  TVirtualFitter *fitter = TVirtualFitter::GetFitter();
  TMatrixD *matrix       = (Int_t(fitCheck)==0) ? new TMatrixD(5,5,fitter->GetCovarianceMatrix()): new TMatrixD(5,5)  ;
    // fitCheck->Print("V"); 
  
  delete hMatrix;
  delete fMatrix;
  return matrix;
}
// -------------------------------------------------------------------------------------------------------
void GetCleanExParams(TH1D *hClean, Double_t *arrSigma, Double_t *arrCleanSigma, Double_t *arrCleanMean, ParticleType pSpecy){

  // 
    // Estimate the expected mean and sigma from the clean samples in order to compare with pid response
  // 

  TVirtualFitter::SetMaxIterations(10000);
  
  maxBin = hClean->GetBinContent(hClean->GetMaximumBin());  
  if ((maxBin<10.)) return;
    
  Double_t mean            = hClean->GetBinCenter(hClean->GetMaximumBin());
  Double_t rms             = arrSigma[pSpecy];
  Double_t sampleFitWindow = 1.5; 
   
  
  
  // fit function
  TF1 *gaus = new TF1("gaus","gaus",mean-sampleFitWindow*rms,mean+sampleFitWindow*rms);   
  gaus->SetParLimits(0,maxBin/2.,maxBin*1.1);
  gaus->SetParLimits(1,mean-rms,mean+rms);
  gaus->SetParLimits(2,rms/2. ,3.*rms);
 
  // retrieve fit parameters
  hClean->Fit(gaus,"QNMR");
  arrCleanMean[pSpecy]  = gaus->GetParameter(1);
  arrCleanSigma[pSpecy] = gaus->GetParameter(2);
    
  delete gaus;
  
}
// -------------------------------------------------------------------------------------------------------
Int_t GetMaxIndex(TGraphErrors *grIndex){

  // Find index of maximum point in the graph
  
  Double_t max = 0;
  Int_t maxIndex = 0;
  for (Int_t i=0;i<grIndex->GetN();i++){
        
    if (grIndex->GetY()[i]>max) {
      max = grIndex->GetY()[i];
      maxIndex = i;
    }
  }

  return maxIndex;
  
}
// -------------------------------------------------------------------------------------------------------
Int_t GetPtIndex(TGraphErrors *grIndex,Double_t p){

  // Find the index for a given momentum
  
  Int_t index = 0;
  for (Int_t i=0;i<grIndex->GetN();i++){ 
    if (grIndex->GetX()[i]>=p) {
      index = i;
      break;
    }
  }

  return index;
}
// -------------------------------------------------------------------------------------------------------
void ResetParametersOfEachFunction(TF1 *g1,TF1 *g2,TF1 *g3,TF1 *g4,TF1 *g5,TF1 *g6,TF1 *total,TH1D *h1D,const Double_t pt){
    
  //
  // Reset the parameters of individual particle after the total fit
  //
    
  // Electron
  g1->SetParameter(0,total->GetParameter(0));
  g1->SetParameter(1,total->GetParameter(1));
  g1->SetParameter(2,total->GetParameter(2));
  g1->SetParameter(3,total->GetParameter(3));
  g1->SetParameter(4,total->GetParameter(4));
   
    // Pion
  g2->SetParameter(0,total->GetParameter(5));
  g2->SetParameter(1,total->GetParameter(6));
  g2->SetParameter(2,total->GetParameter(7));
  g2->SetParameter(3,total->GetParameter(8));
  g2->SetParameter(4,total->GetParameter(9));
        
    // Kaon
  g3->SetParameter(0,total->GetParameter(10));
  g3->SetParameter(1,total->GetParameter(11));
  g3->SetParameter(2,total->GetParameter(12));
  g3->SetParameter(3,total->GetParameter(13));
  g3->SetParameter(4,total->GetParameter(14));
        
    // Proton
  g4->SetParameter(0,total->GetParameter(15));
  g4->SetParameter(1,total->GetParameter(16));
  g4->SetParameter(2,total->GetParameter(17));
  g4->SetParameter(3,total->GetParameter(18));
  g4->SetParameter(4,total->GetParameter(19));
        
    // Deuteron
  g5->SetParameter(0,total->GetParameter(20));
  g5->SetParameter(1,total->GetParameter(21));
  g5->SetParameter(2,total->GetParameter(22));
  g5->SetParameter(3,total->GetParameter(23));
  g5->SetParameter(4,total->GetParameter(24));
    
    // Muons
  g6->SetParameter(0,total->GetParameter(25));
  g6->SetParameter(1,total->GetParameter(26));
  g6->SetParameter(2,total->GetParameter(27));
  g6->SetParameter(3,total->GetParameter(28));
  g6->SetParameter(4,total->GetParameter(29));
  
  
  // Set colors
  g1->SetLineColor(4);
  g2->SetLineColor(1);
  g3->SetLineColor(3);
  g4->SetLineColor(kOrange);
  g5->SetLineColor(6);
  g6->SetLineColor(5); 
  
  // Set ranges 
  g1->SetRange(0,h1D->GetNbinsX());
  g2->SetRange(0,h1D->GetNbinsX());
  g3->SetRange(0,h1D->GetNbinsX());
  g4->SetRange(0,h1D->GetNbinsX());
  g5->SetRange(0,h1D->GetNbinsX());
  g6->SetRange(0,h1D->GetNbinsX());
  
  // add TLegend and Tlines
  TLegend *leg = new TLegend(0.7, 0.6, 0.9, 0.9);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);  // White
  leg->SetHeader(Form("p=%4.3f",pt));
  leg->AddEntry(g1,"Electron");
  leg->AddEntry(g2,"Pion");
  leg->AddEntry(g3,"Kaon");
  leg->AddEntry(g4,"Proton");
  leg->AddEntry(g5,"Deuteron");
  leg->AddEntry(g6,"Muon"); 

  // Attach functions and TLegend to the histogram
  h1D->GetListOfFunctions()->Add(g1);
  h1D->GetListOfFunctions()->Add(g2);
  h1D->GetListOfFunctions()->Add(g3);
  h1D->GetListOfFunctions()->Add(g4);
  h1D->GetListOfFunctions()->Add(g5);
  h1D->GetListOfFunctions()->Add(g6);
  h1D-> GetListOfFunctions()->Add(leg);
  
}
// -------------------------------------------------------------------------------------------------------
Double_t  GetClosestParticleMean(Int_t islice, Double_t sig, Double_t ref, Double_t gr1, Double_t gr2, Double_t gr3){
  //
  // Find the mean position of closest particle for a given slice. This distance will be used as 
  // window range for the setparlimits
  //
  
  Double_t distArr[3];
  distArr[0]=TMath::Abs(ref-gr1);
  distArr[1]=TMath::Abs(ref-gr2);
  distArr[2]=TMath::Abs(ref-gr3);
  Double_t minElement = TMath::MinElement(3,distArr);
  minElement = minElement/windowSetting;
 
  if (minElement > sig)   minElement = sig;
  if (minElement > 2*sig) minElement = 2*sig;
  if (minElement < 0.1)   minElement = 0.1;
    
//   cout << " =============================== returned value =  " << minElement << "  sig = " << sig << endl;
  return minElement;


}
// -------------------------------------------------------------------------------------------------------
TGraphErrors * GraphSmooth(TGraphErrors * gr){
    
    //   convert TH1D to TGraphErrors
    
  const Int_t N = gr->GetN();
  Double_t binFirst = gr->GetX()[0];
  Double_t binLast  = gr->GetX()[N-1];
  
  TH1D * h = new TH1D("h","h",N,binFirst,binLast);
  for (Int_t i=0; i<N+1; i++) {
    h->SetBinContent(i,gr->GetY()[i]);
  }
  
  Double_t lowp = 0; 
  for (Int_t j=0; j<N+1; j++) {
    if (h->GetBinContent(j)>0) {
      lowp=h->GetBinCenter(j+5);
      break;
    }
  }
  
  Double_t upp = 0; 
  for (Int_t j=N; j>0; j--) {
    if (h->GetBinContent(j)>0) {
      upp=h->GetBinCenter(j-5);
      break;
    }
  }
  
  lowp = TMath::Max(0.3,lowp);
  h->GetXaxis()->SetRangeUser(lowp,upp);
  h->Smooth(100000,"R");
  h->SetLineColor(kRed);
 
  TGraphErrors *grS = HistToGraphErrors(h);
  
  return grS;
}
// -------------------------------------------------------------------------------------------------------
TH1D *GraphShade(TGraphErrors *grFitParam, Double_t percent) {
  
//   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
//   if (!c1) c1 = new TCanvas("c1","");
//  c1->Clear();
   
  Int_t n = grFitParam->GetN();
  Int_t nMax = TMath::MaxElement(n,grFitParam->GetY());
  TH1D * hShade= new TH1D("hShade","hShade",n,0,2);
  hShade->GetYaxis()->SetRangeUser(0,2*nMax);
  
  Double_t x[n], y[n],ymin[n], ymax[n];
  Int_t i;
  for (i=0;i<n;i++) {
    x[i]    = grFitParam->GetX()[i];
    ymax[i] = grFitParam->GetY()[i]+percent*grFitParam->GetY()[i];
    ymin[i] = grFitParam->GetY()[i]-percent*grFitParam->GetY()[i];
    y[i]    = grFitParam->GetY()[i];
  }
  TGraph *grmin = new TGraph(n,x,ymin);
  TGraph *grmax = new TGraph(n,x,ymax);
  TGraph *gr    = new TGraph(n,x,y);
  TGraph *grshade = new TGraph(2*n);
  for (i=0;i<n;i++) {
    grshade->SetPoint(i,x[i],ymax[i]);
    grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
  }
  grshade->SetFillStyle(3013);
  grshade->SetFillColor(16);
  grshade->Draw("f");
  grmin->Draw("l");
  grmax->Draw("l");
  gr->SetLineWidth(4);
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->Draw("CP");
  
  gr->SetLineColor(kRed);
  grshade->SetLineColor(kRed);
  grmin->SetLineColor(kRed);
  grmax->SetLineColor(kRed);
  grmin->SetLineStyle(2);
  grmax->SetLineStyle(2);
  
  hShade->GetListOfFunctions()->Add(grshade);
  hShade->GetListOfFunctions()->Add(grmin);
  hShade->GetListOfFunctions()->Add(grmax);
  hShade->GetListOfFunctions()->Add(gr);

  return hShade;
}
// -------------------------------------------------------------------------------------------------------
TGraphErrors * HistToGraphErrors(TH1D * h){
    
    //   convert TH1D to TGraphErrors
    
  const Int_t N = h->GetNbinsX();
  Double_t x[N];
  Double_t y[N];
  Double_t errx[N];
  Double_t erry[N];
    
  for (Int_t i=0; i<N; i++){
    y[i]=h->GetBinContent(i);
    x[i]=i;
    errx[i]=0;
    erry[i]=0;
  }
    
  TGraphErrors *gr = new TGraphErrors(N,x,y,errx,erry);
  return gr;
}
// -------------------------------------------------------------------------------------------------------
void ProduceHighPGraphs(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grTmp, TGraph **grSmoothTmp, TGraphErrors **grAllFit){
  
  //
  // Produce graphs for High momentum part of the fits
  //
  
  TCut hihgpcut = Form("slice>%d",nSlice-1);
  Int_t nEntries = barr->GetEntriesFast();
  cout << "  in the  ProduceHighPGraphs " << nEntries << endl;  
  Int_t objIndex = 0;
  for (Int_t i=0; i<nEntries; ++i){
    
    // Get the brach name 
    TString allfitName = barr->At(i)->GetName();
    TString branchName = barr->At(i)->GetName();
    TString drawStr    = barr->At(i)->GetName();
    drawStr.Append(":p");
    
    // look at only amplitude Mean and Sigma
    if ( !(branchName.Contains("Amp") || branchName.Contains("Mean") || branchName.Contains("Sigma") || branchName.Contains("Int") || branchName.Contains("totalChi2") || branchName.Contains("normNDF") ||  branchName.Contains("totChi2") ||  branchName.Contains("max") || branchName.Contains("ksTest"))  ) continue;
       
    // create the graph object
    t->Draw(drawStr,hihgpcut,"goff");
    grTmp[objIndex] = new TGraphErrors(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    grTmp[objIndex]->SetName(branchName);
    grTmp[objIndex]->GetXaxis()->SetTitle("p(GeV/c)");
    grTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grTmp[objIndex]->SetLineColor(kBlue);
    grTmp[objIndex]->SetMarkerStyle(7);
    grTmp[objIndex]->SetMarkerColor(kBlack);
    
    t->Draw(drawStr,"","goff");
    grAllFit[objIndex] = new TGraphErrors(t->GetSelectedRows(),t->GetV2(),t->GetV1());
    grAllFit[objIndex]->SetName(allfitName.Append("AllFit"));
    grAllFit[objIndex]->GetXaxis()->SetTitle("p(GeV/c)");
    grAllFit[objIndex]->GetYaxis()->SetTitle(allfitName);
    grAllFit[objIndex]->SetLineColor(kBlack);
    grAllFit[objIndex]->SetMarkerStyle(7);
    grAllFit[objIndex]->SetMarkerColor(kBlack);
    
    // create smooth graph 
    gs      = new TGraphSmooth("normal");
    grSmoothTmp[objIndex] = gs->SmoothSuper((TGraph*)grTmp[objIndex],"",3);
    grSmoothTmp[objIndex]->SetName(branchName.Append("Smooth"));
    grSmoothTmp[objIndex]->GetXaxis()->SetTitle("p(GeV/c)");
    grSmoothTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grSmoothTmp[objIndex]->SetLineColor(kRed);
    grSmoothTmp[objIndex]->SetMarkerStyle(7);
    grSmoothTmp[objIndex]->SetMarkerColor(kRed);
        
    // Add graphs to ObjArrays
    arr       -> AddAt(grTmp[objIndex]       ,3*objIndex);
    arr       -> AddAt(grSmoothTmp[objIndex] ,3*objIndex+1);
    arr       -> AddAt(grAllFit[objIndex]    ,3*objIndex+2);
    
    objIndex++;
    
  }
 
}
void ProduceGraphs(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grTmp, TGraph **grSmoothTmp){
  //
  // Produce graphs for MC fit and Clean
  //
  
  Int_t nEntries = barr->GetEntriesFast();
  cout << "  in the  ProduceGraphs " << nEntries << endl;  
  Int_t objIndex = 0;
  for (Int_t i=0; i<nEntries; ++i){
    
    // Get the brach name 
    TString branchName = barr->At(i)->GetName();
    TString drawStr    = barr->At(i)->GetName();
    drawStr.Append(":p");
    
    // look at only amplitude Mean and Sigma
    if ( !(branchName.Contains("Amp") || branchName.Contains("Mean") || branchName.Contains("Sigma") || branchName.Contains("Int") || branchName.Contains("totalChi2") || branchName.Contains("normNDF") || branchName.Contains("Kurtosis") || branchName.Contains("Skew") ||  branchName.Contains("totChi2") ||  branchName.Contains("max") ||  branchName.Contains("ksTest")) ) continue;
       
    Double_t bassLevel = smoothLevel;
    if (branchName.Contains("elFitAmp")) bassLevel=5.;
    if (branchName.Contains("kaFitAmp")) bassLevel=8.;
    if (branchName.Contains("kaFitSigma")) bassLevel=9.;

    // create the graph object
    t->Draw(drawStr,"","goff");
    grTmp[objIndex] = new TGraphErrors(nSlice,t->GetV2(),t->GetV1());
    grTmp[objIndex]->SetName(branchName);
    grTmp[objIndex]->GetXaxis()->SetTitle("p(GeV/c)");
    grTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grTmp[objIndex]->SetLineColor(kBlack);
    grTmp[objIndex]->SetMarkerStyle(7);
    grTmp[objIndex]->SetMarkerColor(kBlack);
    
    // create smooth graph 
    gs      = new TGraphSmooth("normal");
    grSmoothTmp[objIndex] = gs->SmoothSuper((TGraph*)grTmp[objIndex],"",bassLevel);
    grSmoothTmp[objIndex]->SetName(branchName.Append("Smooth"));
    grSmoothTmp[objIndex]->GetXaxis()->SetTitle("p(GeV/c)");
    grSmoothTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grSmoothTmp[objIndex]->SetLineColor(kRed);
    grSmoothTmp[objIndex]->SetMarkerStyle(7);
    grSmoothTmp[objIndex]->SetMarkerColor(kRed);
        
    // Add graphs to ObjArrays
    arr       -> AddAt(grTmp[objIndex]       ,2*objIndex);
    arr       -> AddAt(grSmoothTmp[objIndex] ,2*objIndex+1);
    
    objIndex++;
    
  }
  
}
// -------------------------------------------------------------------------------------------------------
void ApplyOutlierSmooth(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grTmp, TGraph **grSmoothTmp){
  //
  // Produce graphs for MC fit and Clean
  //
  
  
  Int_t nEntries = barr->GetEntriesFast();
  cout << "  in the  ApplyOutlierSmooth " << nEntries << endl;  
  Int_t objIndex = 0;
  for (Int_t i=0; i<nEntries; i++){
    
    TString branchName = barr->At(i)->GetName();
    TString drawStr    = barr->At(i)->GetName();
    drawStr.Append(":p");
    
    // look at only amplitude Mean and Sigma
    if (branchName.Contains("de")) continue;
    if (branchName.Contains("mu")) continue;
    if ( !(branchName.Contains("FitAmp") || branchName.Contains("FitMean") || branchName.Contains("FitSigma")) ) continue;
//     if ( !(branchName.Contains("FitAmp")) ) continue;
   
    if (branchName.Contains("elFitAmp")) smoothSeg=10.;
    if (branchName.Contains("kaFitSigma")) smoothSeg=4.;
      
    // create the graph object
    t->Draw(drawStr,"","goff");
    grTmp[objIndex] = new TGraphErrors(nSlice,t->GetV2(),t->GetV1());
    // Aplly outlier removal
    grTmp[objIndex] = RemoveOutliers(grTmp[objIndex],TMath::Nint((Double_t)nSlice/smoothSeg));
    grTmp[objIndex] = RemoveOutliers(grTmp[objIndex],TMath::Nint((Double_t)nSlice/smoothSeg));
    
    grTmp[objIndex]->SetName(branchName.Append("Outlier"));
    grTmp[objIndex]->GetXaxis()->SetTitle("p(GeV/c)");
    grTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grTmp[objIndex]->SetLineColor(kBlack);
    grTmp[objIndex]->SetMarkerStyle(7);
    grTmp[objIndex]->SetMarkerColor(kBlack);
    
    // create smooth graph 
    gs      = new TGraphSmooth("normal");
    grSmoothTmp[objIndex] = gs->SmoothSuper((TGraph*)grTmp[objIndex],"",9);
    grSmoothTmp[objIndex]->SetName(branchName.Append("Smooth"));
    grSmoothTmp[objIndex]->GetXaxis()->SetTitle("p(GeV/c)");
    grSmoothTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grSmoothTmp[objIndex]->SetLineColor(kRed);
    grSmoothTmp[objIndex]->SetMarkerStyle(7);
    grSmoothTmp[objIndex]->SetMarkerColor(kRed);

    cout << grSmoothTmp[objIndex]->GetName() << endl;
        
    // Add graphs to ObjArrays
    arr -> AddAt(grTmp[objIndex]         ,2*objIndex);
    arr -> AddAt(grSmoothTmp[objIndex]   ,2*objIndex+1);

    objIndex++;
    
  }
  
}
// -------------------------------------------------------------------------------------------------------
TGraphErrors * ApplyScalingWrtExpected(const Int_t nSlice, TString exp, TString fit, const Double_t safeP, TObjArray * arr, TTree *t){
  //
  // Return scaled graph wrt to Expected
  //
  
  // get the graphs to be used in scaling
  TGraphErrors *grexp = (TGraphErrors*)arr->FindObject(exp);
  TGraphErrors *grfit = (TGraphErrors*)arr->FindObject(fit);
  
  // Find the index corresponding to safeP value
  Int_t index = GetPtIndex(grexp,safeP);
  
  // calculate scaling factor  exmean/fitmean (muons are special)
  Double_t scaleFactor;
  TGraphErrors *grmufit,*grmuexp;
  if (exp.Contains("muMean")){
    grmufit     = (TGraphErrors*)arr->FindObject("piFitMean");
    grmuexp     = (TGraphErrors*)arr->FindObject("piMean");
    scaleFactor = grmuexp->GetY()[index]/grmufit->GetY()[index];
  } else if (exp.Contains("muSigma")){
    grmufit     = (TGraphErrors*)arr->FindObject("piFitSigma");
    grmuexp     = (TGraphErrors*)arr->FindObject("piSigma");
    scaleFactor = grmuexp->GetY()[index]/grmufit->GetY()[index];
  } else {
    scaleFactor = grexp->GetY()[index]/grfit->GetY()[index];
  }
  
  // get a new string for newname assignment
  TString newname  = exp;
  
  // draw scaled graph and return it
  exp.Append("/%f:p");
  TString drawstr  = Form(exp,scaleFactor);
  newname.Append("ExScaled");
  
  t->Draw(drawstr,"","goff");
  TGraphErrors *grScaled = new TGraphErrors(nSlice,t->GetV2(),t->GetV1());
  grScaled->SetName(newname);
  grScaled->SetLineColor(kGreen);
  grScaled->SetMarkerStyle(7);
  grScaled->SetMarkerColor(kGreen);
  
  return grScaled;
}
// -------------------------------------------------------------------------------------------------------
TGraphErrors * ApplyScalingWrtSample(const Int_t nSlice, TString sample, TString fit, const Double_t safeP, TObjArray * arr, TObjArray * arrSample, TTree *tSample, TGraph ** grScaledSmooth,Int_t k){
  //
  // Return scaled graph wrt to Expected
  //
  
//   get the graphs to be used in scaling
  TGraphErrors *grsample  = (TGraphErrors*)arrSample -> FindObject(sample);
  TGraphErrors *grfit     = (TGraphErrors*)arr       -> FindObject(fit);
  
//   Find the index corresponding to safeP value
  Int_t index = GetPtIndex(grsample,safeP);
  
//   calculate scaling factor  samplemean/fitmean
  Double_t scaleFactor = grfit->GetY()[index]/grsample->GetY()[index];
  
//   get a new string for newname assignment
  TString newname  = sample;
  
//   draw scaled graph and return it
  sample.Append("*%f:p");
  TString drawstr  = Form(sample,scaleFactor);
  newname.Append("Scaled");
  
//   create scaled graph
  tSample->Draw(drawstr,"","goff");
  TGraphErrors *grScaled = new TGraphErrors(nSlice,tSample->GetV2(),tSample->GetV1());
  grScaled->SetName(newname);
  grScaled->GetXaxis()->SetTitle("p(GeV/c)");
  grScaled->GetYaxis()->SetTitle(newname);
  grScaled->SetLineColor(kGreen);
  grScaled->SetMarkerStyle(7);
  grScaled->SetMarkerColor(kGreen);
  
  gs = new TGraphSmooth("normal");
  grScaledSmooth[k] = gs->SmoothSuper((TGraph*)grScaled,"",4);
  grScaledSmooth[k]->SetName(newname.Append("Smooth"));
  grScaledSmooth[k]->GetXaxis()->SetTitle("p(GeV/c)");
  grScaledSmooth[k]->GetYaxis()->SetTitle(newname);
  grScaledSmooth[k]->SetLineColor(kRed);
  grScaledSmooth[k]->SetMarkerStyle(7);
  grScaledSmooth[k]->SetMarkerColor(kRed);
  
  
  return grScaled;
}
// -------------------------------------------------------------------------------------------------------
TCanvas *GetFitResCanvas(TObjArray *arr){
  
  //   
//   Fit Results Canvas
  //   
  
  TCanvas *can = new TCanvas();
  can->Divide(5,4);
  
  // ++++++++++++++++++  Sigmas +++++++++++++++++
  can->cd(1);
  can->GetPad(1)->SetGrid();
  TGraphErrors *elsigma       = (TGraphErrors*)arr->FindObject("elFitSigma");
  TGraphErrors *elsigmaSmooth = (TGraphErrors*)arr->FindObject("elFitSigmaSmooth");
  elsigma->GetXaxis()->SetTitle("p GeV/c");
  elsigma->GetYaxis()->SetTitle("electron sigma");
  elsigma->GetYaxis()->SetTitleOffset(1.6);
  elsigmaSmooth->SetLineColor(kRed);
  elsigma->Draw("alp");
  elsigmaSmooth->Draw("lp");
      
  can->cd(2);
  can->GetPad(2)->SetGrid();
  TGraphErrors *pisigma       = (TGraphErrors*)arr->FindObject("piFitSigma");
  TGraphErrors *pisigmaSmooth = (TGraphErrors*)arr->FindObject("piFitSigmaSmooth");
  pisigma->GetXaxis()->SetTitle("p GeV/c");
  pisigma->GetYaxis()->SetTitle("pion sigma");
  pisigma->GetYaxis()->SetTitleOffset(1.6);
  pisigmaSmooth->SetLineColor(kRed);
  pisigma->Draw("alp");
  pisigmaSmooth->Draw("lp");

  can->cd(3);
  can->GetPad(3)->SetGrid();
  TGraphErrors *kasigma       = (TGraphErrors*)arr->FindObject("kaFitSigma");
  TGraphErrors *kasigmaSmooth = (TGraphErrors*)arr->FindObject("kaFitSigmaSmooth");
  kasigma->GetXaxis()->SetTitle("p GeV/c");
  kasigma->GetYaxis()->SetTitle("kaon sigma");
  kasigma->GetYaxis()->SetTitleOffset(1.6);
  kasigmaSmooth->SetLineColor(kRed);
  kasigma->Draw("alp");
  kasigmaSmooth->Draw("lp");
  
  can->cd(4);
  can->GetPad(4)->SetGrid();
  TGraphErrors *prsigma       = (TGraphErrors*)arr->FindObject("prFitSigma");
  TGraphErrors *prsigmaSmooth = (TGraphErrors*)arr->FindObject("prFitSigmaSmooth");
  prsigma->GetXaxis()->SetTitle("p GeV/c");
  prsigma->GetYaxis()->SetTitle("proton sigma");
  prsigma->GetYaxis()->SetTitleOffset(1.6);
  prsigmaSmooth->SetLineColor(kRed);
  prsigma->Draw("alp");
  prsigmaSmooth->Draw("lp");
  
  can->cd(5);
  can->GetPad(5)->SetGrid();
  TGraphErrors *chi2          = (TGraphErrors*)arr->FindObject("totChi2");
  chi2->GetXaxis()->SetTitle("p GeV/c");
  chi2->GetYaxis()->SetTitle("Chi2");
  chi2->GetYaxis()->SetTitleOffset(1.6);
  chi2->SetMarkerStyle(7);
  chi2->SetMarkerColor(kRed);
  chi2->SetLineColor(kRed);
  chi2->Draw("alp");
  
  // ++++++++++++++++++  Means +++++++++++++++++
  can->cd(6);
  can->GetPad(6)->SetGrid();
  TGraphErrors *elmean       = (TGraphErrors*)arr->FindObject("elFitMean");
  TGraphErrors *elmeanSmooth = (TGraphErrors*)arr->FindObject("elFitMeanSmooth");
  elmean->GetXaxis()->SetTitle("p GeV/c");
  elmean->GetYaxis()->SetTitle("electron mean");
  elmean->GetYaxis()->SetTitleOffset(1.6);
  elmeanSmooth->SetLineColor(kRed);
  elmean->Draw("alp");
  elmeanSmooth->Draw("lp");
  
  can->cd(7);
  can->GetPad(7)->SetGrid();
  TGraphErrors *pimean       = (TGraphErrors*)arr->FindObject("piFitMean");
  TGraphErrors *pimeanSmooth = (TGraphErrors*)arr->FindObject("piFitMeanSmooth");
  pimean->GetXaxis()->SetTitle("p GeV/c");
  pimean->GetYaxis()->SetTitle("pion mean");
  pimean->GetYaxis()->SetTitleOffset(1.6);
  pimeanSmooth->SetLineColor(kRed);
  pimean->Draw("alp");
  pimeanSmooth->Draw("lp");
  
  can->cd(8);
  can->GetPad(8)->SetGrid();
  TGraphErrors *kamean       = (TGraphErrors*)arr->FindObject("kaFitMean");
  TGraphErrors *kameanSmooth = (TGraphErrors*)arr->FindObject("kaFitMeanSmooth");
  kamean->GetXaxis()->SetTitle("p GeV/c");
  kamean->GetYaxis()->SetTitle("kaon mean");
  kamean->GetYaxis()->SetTitleOffset(1.6);
  kameanSmooth->SetLineColor(kRed);
  kamean->Draw("alp");
  kameanSmooth->Draw("lp");
  
  can->cd(9);
  can->GetPad(9)->SetGrid();
  TGraphErrors *prmean       = (TGraphErrors*)arr->FindObject("prFitMean");
  TGraphErrors *prmeanSmooth = (TGraphErrors*)arr->FindObject("prFitMeanSmooth");
  prmean->GetXaxis()->SetTitle("p GeV/c");
  prmean->GetYaxis()->SetTitle("proton mean");
  prmean->GetYaxis()->SetTitleOffset(1.6);
  prmeanSmooth->SetLineColor(kRed);
  prmean->Draw("alp");
  prmeanSmooth->Draw("lp");

  can->cd(10);
  can->GetPad(10)->SetGrid();
  TGraphErrors *maxRes          = (TGraphErrors*)arr->FindObject("maxRes");
  maxRes->GetXaxis()->SetTitle("p GeV/c");
  maxRes->GetYaxis()->SetTitle("Maximum Residual per Slice");
  maxRes->GetYaxis()->SetTitleOffset(1.6);
  maxRes->SetMarkerStyle(7);
  maxRes->SetMarkerColor(kBlue);
  maxRes->SetLineColor(kBlue);
  maxRes->Draw("alp");

  // ++++++++++++++++++  Amplitudes +++++++++++++++++
  can->cd(11);
  can->GetPad(11)->SetGrid();
  TGraphErrors *elamp       = (TGraphErrors*)arr->FindObject("elFitAmp");
  TGraphErrors *elampSmooth = (TGraphErrors*)arr->FindObject("elFitAmpSmooth");
  elamp->GetXaxis()->SetTitle("p GeV/c");
  elamp->GetYaxis()->SetTitle("electron amplitude");
  elamp->GetYaxis()->SetTitleOffset(1.6);
  elampSmooth->SetLineColor(kRed);
  elamp->Draw("alp");
  elampSmooth->Draw("lp");
  
  can->cd(12);
  can->GetPad(12)->SetGrid();
  TGraphErrors *piamp       = (TGraphErrors*)arr->FindObject("piFitAmp");
  TGraphErrors *piampSmooth = (TGraphErrors*)arr->FindObject("piFitAmpSmooth");
  piamp->GetXaxis()->SetTitle("p GeV/c");
  piamp->GetYaxis()->SetTitle("pion amplitude");
  piamp->GetYaxis()->SetTitleOffset(1.6);
  piampSmooth->SetLineColor(kRed);
  piamp->Draw("alp");
  piampSmooth->Draw("lp");
  
  can->cd(13);
  can->GetPad(13)->SetGrid();
  TGraphErrors *kaamp       = (TGraphErrors*)arr->FindObject("kaFitAmp");
  TGraphErrors *kaampSmooth = (TGraphErrors*)arr->FindObject("kaFitAmpSmooth");
  kaamp->GetXaxis()->SetTitle("p GeV/c");
  kaamp->GetYaxis()->SetTitle("kaon amplitude");
  kaamp->GetYaxis()->SetTitleOffset(1.6);
  kaampSmooth->SetLineColor(kRed);
  kaamp->Draw("alp");
  kaampSmooth->Draw("lp");
  
  can->cd(14);
  can->GetPad(14)->SetGrid();
  TGraphErrors *pramp       = (TGraphErrors*)arr->FindObject("prFitAmp");
  TGraphErrors *prampSmooth = (TGraphErrors*)arr->FindObject("prFitAmpSmooth");
  pramp->GetXaxis()->SetTitle("p GeV/c");
  pramp->GetYaxis()->SetTitle("proton amplitude");
  pramp->GetYaxis()->SetTitleOffset(1.6);
  prampSmooth->SetLineColor(kRed);
  pramp->Draw("alp");
  prampSmooth->Draw("lp");
    
  can->cd(15);
  can->GetPad(15)->SetGrid();
  TGraphErrors *maxPer          = (TGraphErrors*)arr->FindObject("maxPer");
  maxPer->GetXaxis()->SetTitle("p GeV/c");
  maxPer->GetYaxis()->SetTitle("|data-fit|/data");
  maxPer->GetYaxis()->SetTitleOffset(1.6);
  maxPer->SetMarkerStyle(7);
  maxPer->SetMarkerColor(kGreen);
  maxPer->SetLineColor(kGreen);
  maxPer->Draw("alp");
    
    // ++++++++++++++++++  Yields +++++++++++++++++
  can->cd(16);
  can->GetPad(16)->SetGrid();
  TGraphErrors *elint       = (TGraphErrors*)arr->FindObject("elFitInt");
  elint->GetXaxis()->SetTitle("p GeV/c");
  elint->GetYaxis()->SetTitle("electron yield");
  elint->GetYaxis()->SetTitleOffset(1.6);
  elint->Draw("alp");
  
  can->cd(17);
  can->GetPad(17)->SetGrid();
  TGraphErrors *piint       = (TGraphErrors*)arr->FindObject("piFitInt");
  piint->GetXaxis()->SetTitle("p GeV/c");
  piint->GetYaxis()->SetTitle("pion yield");
  piint->GetYaxis()->SetTitleOffset(1.6);
  piint->Draw("alp");
  
  can->cd(18);
  can->GetPad(18)->SetGrid();
  TGraphErrors *kaint       = (TGraphErrors*)arr->FindObject("kaFitInt");
  kaint->GetXaxis()->SetTitle("p GeV/c");
  kaint->GetYaxis()->SetTitle("kaon yield");
  kaint->GetYaxis()->SetTitleOffset(1.6);
  kaint->Draw("alp");
  
  can->cd(19);
  can->GetPad(19)->SetGrid();
  TGraphErrors *print       = (TGraphErrors*)arr->FindObject("prFitInt");
  print->GetXaxis()->SetTitle("p GeV/c");
  print->GetYaxis()->SetTitle("proton yield");
  print->GetYaxis()->SetTitleOffset(1.6);
  print->Draw("alp");
  
  can->cd(20);
  can->GetPad(20)->SetGrid();
  TGraphErrors *ksValue     = (TGraphErrors*)arr->FindObject("ksTest");
  ksValue->GetXaxis()->SetTitle("p GeV/c");
  ksValue->GetYaxis()->SetTitle("ksTest");
  ksValue->GetYaxis()->SetTitleOffset(1.6);
  ksValue->SetMarkerStyle(7);
  ksValue->SetMarkerColor(kMagenta);
  ksValue->SetLineColor(kMagenta);
  ksValue->Draw("alp");
    
  return can;
  
}
// -------------------------------------------------------------------------------------------------------
TGraphErrors *ProduceWindowGraphs(Int_t nSlice, Int_t index, TString str, TTree *tree){

  // 
  //  Prepare graphs around the expected mean of each particle
  //   
  
  tree->Draw(str,"","goff");
  TGraphErrors *grwindow = new TGraphErrors(nSlice,tree->GetV2(),tree->GetV1());
  grwindow->SetName(str);
  grwindow->GetXaxis()->SetTitle("p(GeV/c)");
  grwindow->GetYaxis()->SetTitle(str);
  grwindow->SetMarkerStyle(7);
  if (index%3!=0)  grwindow->SetLineColor(kRed);
  return grwindow;

}
// -------------------------------------------------------------------------------------------------------
Int_t CalculateNBinsP(Int_t nSlice, Int_t ptMin, Int_t ptMax){

  // Calculate all bins 

  Int_t ptDownBin = 0;
  Int_t ptUpBin   = 0; 
  Double_t ptStep = 0.;
  Double_t ptDown = 0.;
  Double_t ptUp   = 0.;
  Double_t pt     = 0.;
  
  TH2D * h = new TH2D("h","h",ptNbins,ptRangeDown,ptRangeUp,dEdxNbins,dEdxMin,dEdxMax);
  
  Int_t ihslice = 0;
  Int_t ibins = 0;
  for (Int_t islice = 0; islice<nSlice+100; islice++){
    
    if (islice<nSlice){
      ptStep = ((ptMax-ptMin)/(Double_t)nSlice);
      ptDown = ptMin+islice*ptStep;
      ptUp   = ptDown+ptStep;
      pt     = (ptDown+ptUp)/2.;
    } else {
      ptStep = highPbinWidth;
      ptDown = ptMax+ihslice*ptStep;
      ptUp   = ptDown+ptStep;
      pt     = (ptDown+ptUp)/2.;
      ihslice++;
    }

    if (ptUp>analysisRange) break;
    
    // retreive the TH1 for each slice
    if (islice<nSlice){
      ptDownBin = (h->GetXaxis()->FindBin(ptMin)-1)+2*islice+1;                                // TO FIX
      ptUpBin   = (h->GetXaxis()->FindBin(ptMin)-1)+2*islice+2;                                // TO FIX
    } else {
      ptDownBin = h->GetXaxis()->FindBin(ptDown);
      ptUpBin   = h->GetXaxis()->FindBin(ptUp)-1;
    }
    
    ibins++;
  }
  
  delete h;
  return ibins;
  
}  
