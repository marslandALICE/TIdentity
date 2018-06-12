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
// ******************************* Modification region **************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************

// Initial Settings
Bool_t MCclosure              = kFALSE;   // kTRUE + KSestimation=kTRUE --> best fits || kTRUE + KSestimation=kFALSE --> gauss fits
Bool_t KSestimation           = kTRUE;   // make KS estimation from samples --> if false use gauss fits
Bool_t AutoK_FixedS           = kTRUE;    // if true only kurtosis is being estimated, output skewness values should be absolute
                                          // otherwise should be <input>*estimatedSkewness
Int_t AutoK_FixedS_Ksetting   = 0;        // kurtosis setting --> 0; pi,K,pr auto, 1; pi auto, K,pr=piKurtosis, 2; pi auto, K,pr = 2
Bool_t fixedK                 = kTRUE;   // fix to auto found (KSestimation=kTRUE) or gauss pars (KSestimation=kFALSE) --> keep true
Bool_t fixedS                 = kTRUE;   // fix to auto found (KSestimation=kTRUE) or gauss pars (KSestimation=kFALSE) --> keep true
Bool_t automaticKS            = kFALSE;   // if true kurtosis and skewness set automatically
Bool_t takeKurtosisFromTable  = kFALSE;
Bool_t specialCareForLargeEta = kFALSE;   // apply only gauss fit for large eta --> >0.6

Int_t nBinsKSestimate         = 50;      // number of bins to be used in the 0.4-1.4 window for KS estimate ; 200 or 50
Bool_t useSafeRegion          = kTRUE;   // use either safe regions or clean samples for KS estimation
Double_t sysAmp               = 0.002;
Double_t sysMean              = 0.002;
Int_t assump                  = 1;        // choose a value 0,1: 0--> use only expecteds 1--> use clean + expected

Double_t piKurtosisScan       = 1.;    // factor to be multiplied with the automatically found kurtosis
Double_t kaKurtosisScan       = 1.;    // factor to be multiplied with the automatically found kurtosis
Double_t prKurtosisScan       = 1.;    // factor to be multiplied with the automatically found kurtosis

Double_t piSkewnessScan       = 0.7;
Double_t kaSkewnessScan       = 0.5;
Double_t prSkewnessScan       = 0.5;

Double_t rebinFactor          = 1.;     // !!!!! dikkat
Int_t nBinsInLookUpTable      =1000;

TString kurtosisTableFile     = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics_cRows_80_16EtaBin_mombin20MeV/Fits/KurtosisTable/KurtosisTable.root";

// dumping bools
Bool_t normalisedCleanSamples = kFALSE;
Bool_t dump                   = kFALSE;  // dump additional plots, not necessary when running for all fits
Bool_t exWRTMean              = kTRUE;   // find

const Int_t nParticle  = 4;
Double_t skewnessFixFit[nParticle]   = {0., 0., 0., 0.};
Double_t skewnessFixClean[nParticle] = {0., 0., 0., 0.};

Double_t kurtosisFixFit[nParticle]   = {2., 2., 2., 2.};
Double_t kurtosisFixClean[nParticle] = {2., 2., 2., 2.};

Double_t kurtosisFix[nParticle]      = {2., 2., 2., 2.};
Double_t skewnessFix[nParticle]      = {0., 0., 0., 0.};

TString smoothFunction = "pol2";   // outlier removal fit function to be used in TLinearFitter
Double_t windowSetting = 4.;       // in window arrangment for mean "distToClosestParticle/windowSetting"
Double_t smoothSeg     = 5.;       // outlier removal range setting for the TLinearFitter (nSlice/smoothSeg)
Double_t smoothLevel   = 9.;       // Increases the iterations in smoothing (10 is too much) for TGraphSmooth
Double_t analysisRange = 5.;       // Maximum p value used in the analysis

Int_t dEdxMin        = 20;
Int_t dEdxMax        = 1020;
Double_t ptRangeDown = 0.2;
Double_t ptRangeUp   = 3.2;

Float_t piFitSafeMin = 0.50,   piFitSafeMax = 0.56;
Float_t kaFitSafeMin = 0.34,   kaFitSafeMax = 0.4;
Float_t prFitSafeMin = 0.62,    prFitSafeMax = 0.7;

Float_t piCleanSafeMin = 0.7,  piCleanSafeMax = 0.8;
Float_t kaCleanSafeMin = 0.8,  kaCleanSafeMax = 1.;
Float_t prCleanSafeMin = 1.3,  prCleanSafeMax = 1.4;

// For free Kurtosis and skewness
Double_t skewMin       = 0.;
Double_t skewMax       = 1.5;
Double_t kurtosisMin   = 1.7;
Double_t kurtosisMax   = 2.2;

// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************* Modification region **************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************


Int_t fSign=0;
Int_t fNSlices = 0;
Double_t parErrors[30]={};
// Define amp mean and sigma parlimits
Double_t elMeanMin ,piMeanMin ,kaMeanMin ,prMeanMin;
Double_t elMeanMax ,piMeanMax ,kaMeanMax ,prMeanMax;
Double_t elAmpMin  ,piAmpMin  ,kaAmpMin  ,prAmpMin;
Double_t elAmpMax  ,piAmpMax  ,kaAmpMax  ,prAmpMax;
Double_t elSigmaMin,piSigmaMin,kaSigmaMin,prSigmaMin;
Double_t elSigmaMax,piSigmaMax,kaSigmaMax,prSigmaMax;

TVectorD *elExMean=0x0, *elExSigma=0x0;
TVectorD *piExMean=0x0, *piExSigma=0x0;
TVectorD *kaExMean=0x0, *kaExSigma=0x0;
TVectorD *prExMean=0x0, *prExSigma=0x0;

// Defien max position in each expected particle mean bin
Double_t maxBin,maxBin0,maxBin1,maxBin2,maxBin3;

// Containers
TClonesArray fitResArr("TH1D", 1000);
TClonesArray fitResiduals("TH1D", 1000);
TClonesArray histLineShapesCArr("TH1D",50000);
TClonesArray funcLineShapesCArr("TF1",50000);
TObjArray funcLineShapesObjArr(50000);

TObjArray    cleanResArr(1000);
TObjArray    cleanResArrFreeKS(1000);


// Debuggers and files
TTreeSRedirector *debugFile     = 0;
TTreeSRedirector *smoothResults = 0;
TTreeSRedirector *histsFile     = 0;
TTreeSRedirector *fitFile       = 0;
TTreeSRedirector *lineShapesLookUp = 0;
TTreeSRedirector *sampleFile    = 0;
TGraphSmooth *gs;

//   Some objects
TObjArray    * Fit;
TObjArray    * Clean;
TObjArray    * ScaledClean;
TObjArray    * ScaledEx;
TObjArray    * OutlierSmooth;
TObjArray    * Windows;

// helper graphs
TGraphErrors *grelFitAmpSmooth=0x0, *grelFitMeanSmooth=0x0, *grelFitSigmaSmooth=0x0;
TGraphErrors *grpiFitAmpSmooth=0x0, *grpiFitMeanSmooth=0x0, *grpiFitSigmaSmooth=0x0;
TGraphErrors *grkaFitAmpSmooth=0x0, *grkaFitMeanSmooth=0x0, *grkaFitSigmaSmooth=0x0;
TGraphErrors *grprFitAmpSmooth=0x0, *grprFitMeanSmooth=0x0, *grprFitSigmaSmooth=0x0;

TGraphErrors *grelFitAmp=0x0, *grelFitMean=0x0, *grelFitSigma=0x0;
TGraphErrors *grpiFitAmp=0x0, *grpiFitMean=0x0, *grpiFitSigma=0x0;
TGraphErrors *grkaFitAmp=0x0, *grkaFitMean=0x0, *grkaFitSigma=0x0;
TGraphErrors *grprFitAmp=0x0, *grprFitMean=0x0, *grprFitSigma=0x0;

  // Extra graphs
TGraphErrors *grelMeanExScaled=0x0,               *grelSigmaExScaled=0x0,                *grelFitAmpOutlierSmooth=0x0;
TGraphErrors *grpiCleanMeanSmooth=0x0,            *grpiCleanSigmaScaledSmooth=0x0,       *grpiCleanMeanScaledSmooth=0x0, *grpiSigmaExScaled=0x0, *grpiMeanExScaled=0x0;
TGraphErrors *grkaCleanMean=0x0, *grkaMeanExScaled=0x0, *grkaSigmaExScaled=0x0, *grkaCleanSigma=0x0, *grkaCleanAmpScaledSmooth=0x0, *grkaFitAmpOutlierSmooth=0x0, *grkaFitSigmaOutlierSmooth=0x0, *grkaCleanSigmaScaled=0x0;
TGraphErrors *grprCleanMean=0x0, *grprMeanExScaled=0x0, *grprSigmaExScaled=0x0, *grprCleanSigma=0x0, *grprCleanAmpScaledSmooth=0x0, *grprCleanSigmaScaledSmooth=0x0, *grprCleanMeanScaledSmooth=0x0, *grprFitAmpOutlierSmooth=0x0;

TGraphErrors *grelCleanSigmaScaled=0x0, *grelCleanMeanScaled=0x0, *grelCleanAmpScaled=0x0,  *grpiCleanMeanScaled=0x0;
TGraphErrors *grpiCleanSigmaScaled=0x0, *grpiCleanAmpScaled=0x0,  *grkaCleanMeanScaled=0x0, *grprCleanSigmaScaled=0x0;
TGraphErrors *grprCleanMeanScaled=0x0;
TGraphErrors *grkaCleanAmpScaled=0x0, *grprCleanAmpScaled=0x0;

  // High momentum part smooth graphs
TGraphErrors *grpiFitAmpSmoothHP=0x0,  *grkaFitAmpSmoothHP=0x0,  *grprFitAmpSmoothHP=0x0;
TGraphErrors *grpiFitMeanSmoothHP=0x0, *grkaFitMeanSmoothHP=0x0, *grprFitMeanSmoothHP=0x0;

// graphs for MC closure test
TGraphErrors *grelFreeKSAmp=0x0, *grelFreeKSMean=0x0, *grelFreeKSSigma=0x0, *grelFreeKSKurtosis=0x0, *grelFreeKSSkew=0x0;
TGraphErrors *grpiFreeKSAmp=0x0, *grpiFreeKSMean=0x0, *grpiFreeKSSigma=0x0, *grpiFreeKSKurtosis=0x0, *grpiFreeKSSkew=0x0;
TGraphErrors *grkaFreeKSAmp=0x0, *grkaFreeKSMean=0x0, *grkaFreeKSSigma=0x0, *grkaFreeKSKurtosis=0x0, *grkaFreeKSSkew=0x0;
TGraphErrors *grprFreeKSAmp=0x0, *grprFreeKSMean=0x0, *grprFreeKSSigma=0x0, *grprFreeKSKurtosis=0x0, *grprFreeKSSkew=0x0;

TGraphErrors *grelCleanAmp=0x0, *grpiCleanAmp=0x0, *grkaCleanAmp=0x0,  *grprCleanAmp=0x0;


TH2D ** hCleanSamples=NULL;
TH2D ** hExpectedMean=NULL;
TH2D ** hExpectedSigma=NULL;
TH2D  * h2D=NULL;
TFile * inputFile=NULL;

const Int_t   fnEtaBins       = 8;///16;    MC: 8, Data:16
const Float_t fEtaRangeDown   = -0.8;///16;
const Float_t fEtaRangeUp     = 0.8;///16;
const Int_t   fnMomBins       = 150;
const Float_t fMomRangeDown   = 0.2;
const Float_t fMomRangeUp     = 3.2;
const Int_t   fnCentBins      = 9;
Float_t     xCentBins[] = {0, 5,  10,  20, 30, 40, 50, 60, 70, 80};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

enum ParticleType {
  kElectron = 0,
  kPion     = 1,
  kKaon     = 2,
  kProton   = 3,
} pType;

//TString fitFunctionGenGaus = "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
TString fitFunctionGenGaus = "[0]*exp(-(abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/1.414213))";

Float_t fEtaDown;
Float_t fCentDown;

Int_t particleBin,centBin,etaBin,momBin,signBin;

// Main Functions for iterative fitting procedure
void          SmoothAmplitudes(TString fileFit, TString fileSample, TString fileOut, const Int_t nSlice, Int_t iIter);
void          GetExandSampleParameters(Int_t islice, TFile *sFile, Double_t *arrMean, Double_t *arrSigma, Double_t *cleanParams);
void          GetHelperGraphs(Int_t iIter, TFile *ressFile);
void          GetCleanExParams(TH1D *hClean, Double_t *arrSigma, Double_t *arrCleanSigma, Double_t *arrCleanMean, ParticleType pSpecy);

void          EstimateKS(Int_t iks, TH2D *h2DKS);
void          FitAllSamples(const Int_t nSlice, const Double_t ptMin, const Double_t ptMax, TString sampleFileName);
void          IterativeFitting(Int_t iIter, const Int_t nSlice, const Double_t ptMin, const Double_t ptMax, TString fileIn, TString readSmooth, TString readSamp,TString lineShapesFile);
void          ProduceGraphs(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grCleanTmp, TGraph **grCleanSmoothTmp);
void          ApplyOutlierSmooth(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grTmp, TGraph **grSmoothTmp);
void          ApplyInitialFits(Int_t islice, Double_t pt, TH1D *h1D, Double_t *arrSigma, Double_t *arrMean);
void          SetStyleFor1DSlice(TH1D *h);
void          FixParamsMC(Int_t islice, Double_t pt);
void          SetParamsMC(Int_t islice, Double_t pt);
void          ReadHistograms();
void          GetKurtosisFromFile();


// Helper functions for severel computations
void          Fit1DSliceForKSEstimate(Int_t iks, TH1D *h1DKS, Double_t *pars, Double_t *arrMean, Double_t *arrSigma);
void          FitParticleSampleForKSEstimate(Int_t iks, TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy);
TMatrixD     *FitParticleSample(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy);
void          FitParticleSampleFreeKS(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy);
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
Double_t      GetClosestParticleMean(Double_t sig, Double_t ref, Double_t gr1, Double_t gr2, Double_t gr3);
TGraphErrors *GraphSmooth(TGraphErrors * gr);
TGraphErrors *HistToGraphErrors(TH1D * h);
TCanvas      *GetFitResCanvas(TObjArray *arr);

// Set functions for the parameters
void          SetTotalFitParameters(Int_t iIter, Int_t islice, Int_t nSlice, TF1* total, Double_t *arrMean, Double_t *arrSigma, Double_t *arrMeanWindow, Double_t *cleanParams, TH1D* h1D, Double_t pt);
void          CheckIfMeanIsZero(TF1 *total, Double_t *arrMean);
void          SetParNamesOfTotalFit(TF1 *total);
void          SetParamsForKSEstimate(Int_t iks, TF1 *total, TH1D *h1DKS, Double_t *arrMean, Double_t *arrSigma);
void          SetFixedKSparameterForTotalFit(TF1 *total);
void          SetFixedKSparameterForTotalFitForPions(TF1 *total);
void          SetFixedKSparameterForCleanSampleFit(TF1 *asyGaus, ParticleType pSpecy);
void          SetFixedKSparameterForIndividualFits(TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4);
void          PrepareLineShapes(TH1D *h,TF1 *ftot,TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4, Int_t centbin, Int_t etabin, Int_t mombin, Int_t signbin, Int_t iIter);
void          ResetParametersOfEachFunction(TF1 *g1,TF1 *g2,TF1 *g3,TF1 *g4,TF1 *total,TH1D *h1D,const Double_t pt);
void          SetParams1stIteration(Double_t pt, TH1D *h1D, Double_t *arrMean, Double_t *arrSigma,Double_t *arrMeanWindow, Double_t *cleanParams);
void          SetParams1stIterationMC(Double_t pt, TH1D *h1D, Double_t * cleanParams);



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
void PIDIterativeFitting(Int_t sign, TString fileData, const Int_t nSlice, Double_t ptMin, Double_t ptMax, Int_t maxIter=8, Double_t piSkewnessScanIn=0.8, Double_t kaSkewnessScanIn=0.5, Double_t prSkewnessScanIn=0.5, Double_t piKurtosisScanIn=1., Double_t kaKurtosisScanIn=1., Double_t prKurtosisScanIn=1.)
{

  //
  // Apply gaussian fits to the ptot slices --> Main Function
  /*

  cd  /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/Fits/test
  aliroot -l -b
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/PIDIterativeFitting.C+
  PIDIterativeFitting(1,"../../Hists/Hists_PbPb_eta_-0.10_0.00_cent_30.00_40.00.root",100,0.2,2.2,1,0.7,0.5,0.5,1,1,1);  2> err.log

  */
  //
  fSign=sign;
  fNSlices=nSlice;
  // Apply settings for the scan
  piSkewnessScan = piSkewnessScanIn;
  kaSkewnessScan = kaSkewnessScanIn;
  prSkewnessScan = prSkewnessScanIn;

  piKurtosisScan = piKurtosisScanIn;
  kaKurtosisScan = kaKurtosisScanIn;
  prKurtosisScan = prKurtosisScanIn;

  TH1D *fhEta   = new TH1D("fhEta" ,"Eta Bins"       ,fnEtaBins ,fEtaRangeDown, fEtaRangeUp );
  TH1D *fhPtot  = new TH1D("fhPtot","Momentum Bins"  ,fnMomBins ,fMomRangeDown, fMomRangeUp );
  TH1D *fhCent  = new TH1D("fhCent","Centrality Bins",fnCentBins ,xCentBins );


  if (!KSestimation && !automaticKS && !MCclosure) {

    skewnessFixFit[kElectron] = 0.;
    skewnessFixFit[kPion]     = piSkewnessScanIn;
    skewnessFixFit[kKaon]     = kaSkewnessScanIn;
    skewnessFixFit[kProton]   = prSkewnessScanIn;

    skewnessFixClean[kElectron] = 0.;
    skewnessFixClean[kPion]     = piSkewnessScanIn;
    skewnessFixClean[kKaon]     = kaSkewnessScanIn;
    skewnessFixClean[kProton]   = prSkewnessScanIn;

    skewnessFix[kElectron] = 0.;
    skewnessFix[kPion]     = piSkewnessScanIn;
    skewnessFix[kKaon]     = kaSkewnessScanIn;
    skewnessFix[kProton]   = prSkewnessScanIn;

    kurtosisFixFit[kElectron] = 2.;
    kurtosisFixFit[kPion]     = piKurtosisScanIn*2.;
    kurtosisFixFit[kKaon]     = kaKurtosisScanIn*2.;
    kurtosisFixFit[kProton]   = prKurtosisScanIn*2.;

    kurtosisFixClean[kElectron] = 2.;
    kurtosisFixClean[kPion]     = piKurtosisScanIn*2.;
    kurtosisFixClean[kKaon]     = kaKurtosisScanIn*2.;
    kurtosisFixClean[kProton]   = prKurtosisScanIn*2.;

    kurtosisFix[kElectron] = 2.;
    kurtosisFix[kPion]     = piKurtosisScanIn*2.;
    kurtosisFix[kKaon]     = kaKurtosisScanIn*2.;
    kurtosisFix[kProton]   = prKurtosisScanIn*2.;

  }

  // Get eta and cent info
  TObjArray *objArr1  = fileData.Tokenize("/");
  TString objFileName = ((objArr1->At((Int_t)objArr1->GetLast()))->GetName());
  TObjArray *objArr2  = objFileName.Tokenize("_");
  fEtaDown     = atof((objArr2->At(3))->GetName());
  fCentDown    = atof((objArr2->At(6))->GetName());
  Double_t etaBinCenter = (atof((objArr2->At(3))->GetName()) + atof((objArr2->At(4))->GetName()))/2.;

  etaBin    = fhEta -> FindBin(fEtaDown - 0.0000001)  + 1;
  centBin   = fhCent-> FindBin(fCentDown - 0.0000001) + 1;
  signBin   = sign+1;


   // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "    << endl;
  cout << " ================================== "          << endl;
  cout << " KSestimation   = " << KSestimation            << endl;
  cout << " ================================== "          << endl;
  cout << " AutoK_FixedS   = " << AutoK_FixedS            << endl;
  cout << " ================================== "          << endl;
  cout << " automaticKS    = " << automaticKS             << endl;
  cout << " ================================== "          << endl;
  cout << " CareLargeEta   = " << specialCareForLargeEta  << endl;
  cout << " ================================== "          << endl;
  cout << " useSafeRegion  = " << useSafeRegion           << endl;
  cout << " ================================== "          << endl;
  cout << " nBinsKSestimate= " << nBinsKSestimate         << endl;
  cout << " ================================== "    << endl;
  cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "    << endl;
  cout << " ================================== "    << endl;
  cout << " fixedK         = " << fixedK            << endl;
  cout << " ================================== "    << endl;
  cout << " fixedS         = " << fixedS            << endl;
  cout << " ================================== "    << endl;
  cout << " eta            = " << fEtaDown           << endl;
  cout << " ================================== "    << endl;
  cout << " cent           = " << fCentDown           << endl;
  cout << " ================================== "    << endl;
  cout << " file Path      = " << fileData          << endl;
  cout << " ================================== "    << endl;
  cout << " file name      = " << objFileName       << endl;
  cout << " ================================== "    << endl;
  cout << " piSkewness     = " << piSkewnessScan    << endl;
  cout << " ================================== "    << endl;
  cout << " kaSkewness     = " << kaSkewnessScan    << endl;
  cout << " ================================== "    << endl;
  cout << " prSkewness     = " << prSkewnessScan    << endl;
  cout << " ================================== "    << endl;
  cout << " piKurtScan     = " << piKurtosisScan    << endl;
  cout << " ================================== "    << endl;
  cout << " kaKurtScan     = " << kaKurtosisScan    << endl;
  cout << " ================================== "    << endl;
  cout << " prKurtScan     = " << prKurtosisScan    << endl;
  cout << " ================================== "    << endl;
  cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "    << endl;

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // Create an debug file
  debugFile = new TTreeSRedirector(Form("debugFile_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown),"recreate");
  inputFile = TFile::Open(fileData);
  ReadHistograms();

  //  Estimate kurtosis and Skewness for kaon, proton and pion by fitting 100 slices in between [0.4,1.4]GeV/c
  cout << "Get The Clone of h2D " << h2D->GetNbinsX() << "   " << h2D->GetName() << endl;
  TH2D *h2DKS = (TH2D*)h2D->Clone();
  if (KSestimation){
    TStopwatch timer; timer.Start();
    for (Int_t iks = 0; iks<2; iks++){
      cout << " =========== KS estimate  ============  Iteration = " << iks << endl;
      EstimateKS(iks,h2DKS);
    }
    timer.Stop(); timer.Print();

    // semi auto KS estimation skewness is given from outside
    if (AutoK_FixedS) {

      // set skewwness
      skewnessFix[kElectron] = 0.;
      skewnessFix[kPion]     = piSkewnessScanIn;
      skewnessFix[kKaon]     = kaSkewnessScanIn;
      skewnessFix[kProton]   = prSkewnessScanIn;
      // set kurtosis
      if (AutoK_FixedS_Ksetting==0){
	cout << " piKurtosisScanIn = " << piKurtosisScanIn <<"    " << kurtosisFixFit[kPion] << endl;
	cout << " kaKurtosisScanIn = " << kaKurtosisScanIn <<"    " << kurtosisFixClean[kKaon]<< endl;
	cout << " prKurtosisScanIn = " << prKurtosisScanIn <<"    " << kurtosisFixFit[kProton]<< endl;

	kurtosisFix[kPion]     = kurtosisFixFit[kPion]*piKurtosisScanIn;
	kurtosisFix[kKaon]     = kurtosisFixClean[kKaon]*kaKurtosisScanIn;
	kurtosisFix[kProton]   = kurtosisFixFit[kProton]*prKurtosisScanIn;
	kurtosisFix[kElectron] = 2.;
      }
      if (AutoK_FixedS_Ksetting==1){
	kurtosisFix[kPion]     = kurtosisFixFit[kPion]*piKurtosisScanIn;
	kurtosisFix[kKaon]     = kurtosisFix[kPion];
	kurtosisFix[kProton]   = kurtosisFix[kPion];
	kurtosisFix[kElectron] = 2.;
      }
      if (AutoK_FixedS_Ksetting==2){
	kurtosisFix[kPion]     = kurtosisFixFit[kPion]*piKurtosisScanIn;
	kurtosisFix[kKaon]     = 2.;
	kurtosisFix[kProton]   = 2.;
	kurtosisFix[kElectron] = 2.;
      }
      if (MCclosure){
	kurtosisFix[kPion]     = kurtosisFixClean[kPion]*piKurtosisScanIn;
	kurtosisFix[kKaon]     = kurtosisFixClean[kKaon]*kaKurtosisScanIn;
	kurtosisFix[kProton]   = kurtosisFixClean[kProton]*prKurtosisScanIn;
	kurtosisFix[kElectron] = 2.;
      }
    }

  }

  cout << " ================================== "    << endl;
  cout << "skewnessFixElectron = " <<  skewnessFix[kElectron] << endl;
  cout << "skewnessFixPion     = " <<  skewnessFix[kPion]     << endl;
  cout << "skewnessFixKaon     = " <<  skewnessFix[kKaon]     << endl;
  cout << "skewnessFixProton   = " <<  skewnessFix[kProton]   << endl;
  cout << " ================================== "    << endl;
  cout << "kurtosisFixElectron = " <<  kurtosisFix[kElectron] << endl;
  cout << "kurtosisFixPion     = " <<  kurtosisFix[kPion]     << endl;
  cout << "kurtosisFixKaon     = " <<  kurtosisFix[kKaon]     << endl;
  cout << "kurtosisFixProton   = " <<  kurtosisFix[kProton]   << endl;
  cout << " ================================== "    << endl;

  // For systenatic study use kurtosis from the table
  if (takeKurtosisFromTable) GetKurtosisFromFile();

  // First fit the samples
  cout << " =========== first make sample file and clean them up =========== " << endl;
  TString fileSample  = Form("Samples_%d_%d_%3.2f_%3.2f_%3.2f_%d.root",fSign,nSlice,ptMin,ptMax,fEtaDown,0);      // care only eta dependence
  FitAllSamples(nSlice,ptMin,ptMax,fileSample);

  // delete unnecessary histograms in memory
  for (Int_t icache=0; icache<nParticle; icache++)
  {
    delete hCleanSamples[icache];
    delete hExpectedMean[icache];
    delete hExpectedSigma[icache];
  }
  delete [] hCleanSamples;
  delete [] hExpectedMean;
  delete [] hExpectedSigma;
  delete h2DKS;
  if (!dump) delete debugFile;

  // for large eta bin do not apply auto ks estimation
  if (specialCareForLargeEta && (TMath::Abs(etaBinCenter)>=0.6) && KSestimation )
  {
    for (Int_t i=0;i<nParticle;i++) {skewnessFix[i] = 0.; kurtosisFix[i] = 2.;}
  }

  // Check the fit params just incase kurt and skewness 0
  for (Int_t i=0;i<nParticle;i++) if ((skewnessFix[i]<=-0.1) || (kurtosisFix[i]<1.)) {skewnessFix[i] = 0.; kurtosisFix[i] = 2.;}

  // iteratively apply fits
  TStopwatch timer; timer.Start();
  for (Int_t iIter = 0; iIter<10; iIter++){
    if(iIter>maxIter) break;
    TString fileFit    = Form("Trees_Iter%d_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",iIter,fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown);
    TString lineShapesFile = Form("LineShapes_%d_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",iIter,fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown);
    TString fileOut    = Form("Results_Iter%d_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",iIter,fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown);
    Int_t previousIter = (MCclosure && iIter==6) ? 1 : iIter-1;
//     Int_t previousIter = iIter-1;
    TString readSmooth = Form("Results_Iter%d_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",previousIter,fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown);


    cout << " =========== Iterative fit procedure is being started  ============  Iteration = " << iIter << endl;
    IterativeFitting(iIter,nSlice,ptMin,ptMax,fileFit,readSmooth,fileSample,lineShapesFile);

    cout << " =========== Iteration = " << iIter << "   has finished go to iteration = " << iIter+1 << endl;
    if (nSlice>50) SmoothAmplitudes(fileFit,fileSample,fileOut,nSlice,iIter);
  }
  timer.Stop(); timer.Print();
  delete objArr1;
  delete objArr2;

}
// -------------------------------------------------------------------------------------------------------
void SmoothAmplitudes(TString fileFit, TString fileSample, TString fileOut, const Int_t nSlice, Int_t iIter)
{

  //
  // Analyse first results and smooth the amplitude graphs
  /*
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/tests/trash
  aliroot -l
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/PIDIterativeFitting.C+
  SmoothAmplitudes("Trees_Iter4_300_0.2_3.2_-0.2_0.root","Samples_300_0.2_3.2_-0.2_0.root","SmoothTest.root",300,0)
  */
  //

  smoothResults = new TTreeSRedirector(fileOut,"recreate");

//   Open "Trees_Iteration0_200_0.2_2.2_0.0_-2.root" file and read trees
  TFile *ftrees = TFile::Open(fileFit);
  TFile *strees = TFile::Open(fileSample);

  TTree *treeFit   = (TTree*)ftrees->Get("FitResults");
  TTree *treeClean = (TTree*)strees->Get("CleanSamples");

//   GetList of branches for each tree
  TObjArray *branchArrFit   =  (TObjArray*)treeFit  ->GetListOfBranches();
  TObjArray *branchArrClean =  (TObjArray*)treeClean->GetListOfBranches();

  Int_t nArrEntries      = branchArrFit->GetEntriesFast();
  Int_t nArrEntriesClean = branchArrClean->GetEntriesFast();


//   Output TObjArrays
  TObjArray *GrArrClean        = new TObjArray(2*nArrEntriesClean);  GrArrClean         -> SetOwner(kTRUE);
  TObjArray *GrArrFit          = new TObjArray(2*nArrEntries);       GrArrFit           -> SetOwner(kTRUE);
  TObjArray *GrArrScaledClean  = new TObjArray(2*nArrEntriesClean);  GrArrScaledClean   -> SetOwner(kTRUE);
  TObjArray *GrArrScaledEx     = new TObjArray(2*nArrEntries);       GrArrScaledEx      -> SetOwner(kTRUE);

  TObjArray *GrArrOutlierSmooth = new TObjArray(2*nArrEntries); GrArrOutlierSmooth -> SetOwner(kTRUE);
  TObjArray *GrArrWindows       = new TObjArray(2*nArrEntries); GrArrWindows       -> SetOwner(kTRUE);

//   Graph Arrays
  TGraphErrors **grClean         = new TGraphErrors *[nArrEntriesClean];
  TGraphErrors **grFit           = new TGraphErrors *[nArrEntries];
  TGraphErrors **grOutlierCut    = new TGraphErrors *[nArrEntries];
  TGraph       **grFitSmooth     = new TGraph *[nArrEntries];
  TGraph       **grCleanSmooth   = new TGraph *[nArrEntriesClean];
  TGraph       **grOutlierSmooth = new TGraph *[nArrEntries];

  TGraphErrors **grScaledClean       = new TGraphErrors *[nArrEntries];
  TGraph       **grScaledCleanSmooth = new TGraph *[nArrEntries];

  TGraphErrors **grScaledEx          = new TGraphErrors *[nArrEntries];

  TGraphErrors * grWindowsAroundMean[nArrEntries];

//   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//   Produce smooth and normal graphs
  ProduceGraphs(nSlice, treeFit  , GrArrFit  , branchArrFit  , grFit  , grFitSmooth);
  ProduceGraphs(nSlice, treeClean, GrArrClean, branchArrClean, grClean, grCleanSmooth);

//   Apply outlier removal to fits
  ApplyOutlierSmooth(nSlice, treeFit, GrArrOutlierSmooth, branchArrFit, grOutlierCut, grOutlierSmooth);

//   Apply Scalings for expected mean and sigma
  Double_t scalingPointKaon = (MCclosure) ? 0.32 : 0.4;
  Double_t scalingPointKaon2 = (MCclosure) ? 0.32 : 0.65;
  grScaledEx[0] = ApplyScalingWrtExpected(nSlice, "elMean","elFitMean",0.32,GrArrFit,treeFit);
  grScaledEx[1] = ApplyScalingWrtExpected(nSlice, "piMean","piFitMean",0.55,GrArrFit,treeFit);
  grScaledEx[2] = ApplyScalingWrtExpected(nSlice, "kaMean","kaFitMean",scalingPointKaon,GrArrFit,treeFit);
  grScaledEx[3] = ApplyScalingWrtExpected(nSlice, "prMean","prFitMean",0.7,GrArrFit,treeFit);
  grScaledEx[4] = ApplyScalingWrtExpected(nSlice, "elSigma","elFitSigma",0.32,GrArrFit,treeFit);
  grScaledEx[5] = ApplyScalingWrtExpected(nSlice, "piSigma","piFitSigma",0.55,GrArrFit,treeFit);
  grScaledEx[6] = ApplyScalingWrtExpected(nSlice, "kaSigma","kaFitSigma",scalingPointKaon,GrArrFit,treeFit);
  grScaledEx[7] = ApplyScalingWrtExpected(nSlice, "prSigma","prFitSigma",0.7,GrArrFit,treeFit);
  for (Int_t i=0; i<8; i++) GrArrScaledEx->AddAt(grScaledEx[i],i);

//   Apply Scalings clean data mean and sigma KAON and PROTON AMPLITUDE !!!!!!
  grScaledClean[0]  = ApplyScalingWrtSample(nSlice,"elCleanAmp","elFitAmp"  ,0.25,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,0);
  grScaledClean[1]  = ApplyScalingWrtSample(nSlice,"piCleanAmp","piFitAmp"  ,0.55,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,1);
  grScaledClean[2]  = ApplyScalingWrtSample(nSlice,"kaCleanAmp","kaFitAmp"  ,scalingPointKaon2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,2);
  grScaledClean[3]  = ApplyScalingWrtSample(nSlice,"prCleanAmp","prFitAmp"  ,1.2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,3);

  grScaledClean[4]  = ApplyScalingWrtSample(nSlice,"elCleanMean","elFitMean",0.25,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,4);
  grScaledClean[5]  = ApplyScalingWrtSample(nSlice,"piCleanMean","piFitMean",0.55,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,5);
  grScaledClean[6]  = ApplyScalingWrtSample(nSlice,"kaCleanMean","kaFitMean",scalingPointKaon2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,6);
  grScaledClean[7]  = ApplyScalingWrtSample(nSlice,"prCleanMean","prFitMean",1.2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,7);

  grScaledClean[8]  = ApplyScalingWrtSample(nSlice,"elCleanSigma","elFitSigma",0.25,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,8);
  grScaledClean[9]  = ApplyScalingWrtSample(nSlice,"piCleanSigma","piFitSigma",0.55,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,9);
  grScaledClean[10] = ApplyScalingWrtSample(nSlice,"kaCleanSigma","kaFitSigma",scalingPointKaon2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,10);
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
  canResults->SaveAs(Form("summary_iter%d_%d_eta%2.1f_cent%3.2f.png",iIter,fSign,fEtaDown,fCentDown));

//   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  smoothResults->GetFile()->cd();
  // Helper graphs
  canResults         -> Write("AllFitParams");
  totChi2            -> Write("totChi2");
  maxPer             -> Write("maxDeviation");
  maxRes             -> Write("maxResidual");

  // Arrays
  GrArrFit           -> Write("Fit"         ,TObject::kSingleKey);
  GrArrClean         -> Write("Clean"       ,TObject::kSingleKey);
  GrArrScaledEx      -> Write("ScaledEx"    ,TObject::kSingleKey);
  GrArrScaledClean   -> Write("ScaledClean" ,TObject::kSingleKey);
  GrArrWindows       -> Write("Windows"     ,TObject::kSingleKey);
  GrArrOutlierSmooth -> Write("OulierSmooth",TObject::kSingleKey);

//   delete arrays
  GrArrFit           -> Delete();
  GrArrClean         -> Delete();
  GrArrScaledEx      -> Delete();
  GrArrScaledClean   -> Delete();
  GrArrOutlierSmooth -> Delete();
  delete GrArrFit;
  delete GrArrClean;
  delete GrArrScaledEx;
  delete GrArrScaledClean;
  delete GrArrOutlierSmooth;

//   delete helper graphs
  delete canResults;
  delete totChi2;
  delete maxRes;
  delete maxPer;

//   delete outputfile
  delete smoothResults;
  ftrees->Close();
  strees->Close();

}
// -------------------------------------------------------------------------------------------------------
void  SetParams1stIterationMC(Double_t pt, TH1D *h1D, Double_t * cleanParams)
{

  //
  //    Set mean amp and sigma par limits for the first iteration
  //

  //   Part to make setting for fit params
  maxBin  = h1D->GetBinContent(h1D->GetMaximumBin());
  maxBin0 = TMath::Min(h1D->GetBinContent(h1D->FindBin(cleanParams[0])),maxBin);
  maxBin1 = TMath::Min(h1D->GetBinContent(h1D->FindBin(cleanParams[1])),maxBin);
  maxBin2 = TMath::Min(h1D->GetBinContent(h1D->FindBin(cleanParams[2])),maxBin);
  maxBin3 = TMath::Min(h1D->GetBinContent(h1D->FindBin(cleanParams[3])),maxBin);

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


  elMeanMin = TMath::Max(30.,(cleanParams[0]-fitWinMean*cleanParams[4]));
  piMeanMin = TMath::Max(30.,(cleanParams[1]-fitWinMean*cleanParams[5]));
  kaMeanMin = TMath::Max(30.,(cleanParams[2]-fitWinMean*cleanParams[6]));
  prMeanMin = TMath::Max(30.,(cleanParams[3]-fitWinMean*cleanParams[7]));


  elMeanMax = TMath::Min((Double_t)dEdxMax,(cleanParams[0]+fitWinMean*cleanParams[4]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(cleanParams[1]+fitWinMean*cleanParams[5]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(cleanParams[2]+fitWinMean*cleanParams[6]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(cleanParams[3]+fitWinMean*cleanParams[7]));


  elSigmaMin = cleanParams[4]/2.;
  piSigmaMin = cleanParams[5]/2.;
  kaSigmaMin = cleanParams[6]/2.;
  prSigmaMin = cleanParams[7]/2.;


  elSigmaMax = cleanParams[4]*2.;
  piSigmaMax = cleanParams[5]*2.;
  kaSigmaMax = cleanParams[6]*2.;
  prSigmaMax = cleanParams[7]*2.;

  TH1D *htmppr = (TH1D*)h1D->Clone();
  TH1D *htmpka = (TH1D*)h1D->Clone();
  htmppr->GetXaxis()->SetRangeUser(250,800);
  htmpka->GetXaxis()->SetRangeUser(100,250);

  delete htmppr;
  delete htmpka;

}
// -------------------------------------------------------------------------------------------------------
void SetParamsMC(Int_t islice, Double_t pt)
{

  //
//  Set Params for the iterative of the h2DMC
  //

  Double_t windowRange = 0.001;

  //   ***************  Electron Patches ***************
  if (pt>0.3) {
    elSigmaMin = grelSigmaExScaled->GetY()[islice]-grelSigmaExScaled->GetY()[islice]*0.05;
    elSigmaMax = grelSigmaExScaled->GetY()[islice]+grelSigmaExScaled->GetY()[islice]*0.05;
    elMeanMin  = grelMeanExScaled->GetY()[islice] -grelMeanExScaled->GetY()[islice]*0.05;
    elMeanMax  = grelMeanExScaled->GetY()[islice] +grelMeanExScaled->GetY()[islice]*0.05;
    elAmpMin   = 0.5;
  }

  //   ***************  Pion Patches ***************
  if ((pt>=0.6)){
    piMeanMin  = grpiCleanMeanScaled->GetY()[islice] -grpiCleanMeanScaled->GetY()[islice]*windowRange;
    piMeanMax  = grpiCleanMeanScaled->GetY()[islice] +grpiCleanMeanScaled->GetY()[islice]*windowRange;
    piSigmaMin = grpiCleanSigmaScaled->GetY()[islice]-grpiCleanSigmaScaled->GetY()[islice]*0.01;
    piSigmaMax = grpiCleanSigmaScaled->GetY()[islice]+grpiCleanSigmaScaled->GetY()[islice]*0.01;
  }

  //   ***************  kaon  patches ***************
  if (pt>=0.3) {
    kaMeanMin  = grkaCleanMeanScaled->GetY()[islice] -grkaCleanMeanScaled->GetY()[islice]*windowRange;
    kaMeanMax  = grkaCleanMeanScaled->GetY()[islice] +grkaCleanMeanScaled->GetY()[islice]*windowRange;
    kaSigmaMin = grkaCleanSigmaScaled->GetY()[islice]-grkaCleanSigmaScaled->GetY()[islice]*0.01;
    kaSigmaMax = grkaCleanSigmaScaled->GetY()[islice]+grkaCleanSigmaScaled->GetY()[islice]*0.01;
  }

  // ***************  Proton Patches ***************
  if ((pt>=0.8)){
    prMeanMin  = grprCleanMeanScaled->GetY()[islice] -grprCleanMeanScaled->GetY()[islice]*windowRange;
    prMeanMax  = grprCleanMeanScaled->GetY()[islice] +grprCleanMeanScaled->GetY()[islice]*windowRange;
    prSigmaMin = grprCleanSigmaScaled->GetY()[islice]-grprCleanSigmaScaled->GetY()[islice]*0.01;
    prSigmaMax = grprCleanSigmaScaled->GetY()[islice]+grprCleanSigmaScaled->GetY()[islice]*0.01;
  }

}
// -------------------------------------------------------------------------------------------------------
void FixParamsMC(Int_t islice, Double_t pt)
{

  //
//  Set Params for the iterative of the h2DMC
  //

  Double_t windowRange = 0.001;

  elSigmaMin = grelSigmaExScaled->GetY()[islice]-grelSigmaExScaled->GetY()[islice]*0.1;
  elSigmaMax = grelSigmaExScaled->GetY()[islice]+grelSigmaExScaled->GetY()[islice]*0.1;
  elMeanMin  = grelMeanExScaled->GetY()[islice] -grelMeanExScaled->GetY()[islice]*0.1;
  elMeanMax  = grelMeanExScaled->GetY()[islice] +grelMeanExScaled->GetY()[islice]*0.1;
  elAmpMin   = 0.5;

  piSigmaMin = grpiFreeKSSigma->GetY()[islice]-grpiFreeKSSigma->GetY()[islice]*windowRange;
  piSigmaMax = grpiFreeKSSigma->GetY()[islice]+grpiFreeKSSigma->GetY()[islice]*windowRange;
  piMeanMin  = grpiFreeKSMean->GetY()[islice] -grpiFreeKSMean->GetY()[islice]*windowRange;
  piMeanMax  = grpiFreeKSMean->GetY()[islice] +grpiFreeKSMean->GetY()[islice]*windowRange;
  piAmpMin   = grpiFreeKSAmp->GetY()[islice]  -grpiFreeKSAmp->GetY()[islice]*windowRange;
  piAmpMax   = grpiFreeKSAmp->GetY()[islice]  +grpiFreeKSAmp->GetY()[islice]*windowRange;

  kaSigmaMin = grkaFreeKSSigma->GetY()[islice]-grkaFreeKSSigma->GetY()[islice]*windowRange;
  kaSigmaMax = grkaFreeKSSigma->GetY()[islice]+grkaFreeKSSigma->GetY()[islice]*windowRange;
  kaMeanMin  = grkaFreeKSMean->GetY()[islice] -grkaFreeKSMean->GetY()[islice]*windowRange;
  kaMeanMax  = grkaFreeKSMean->GetY()[islice] +grkaFreeKSMean->GetY()[islice]*windowRange;
  kaAmpMin   = grkaFreeKSAmp->GetY()[islice]  -grkaFreeKSAmp->GetY()[islice]*windowRange;
  kaAmpMax   = grkaFreeKSAmp->GetY()[islice]  +grkaFreeKSAmp->GetY()[islice]*windowRange;

  prSigmaMin = grprFreeKSSigma->GetY()[islice]-grprFreeKSSigma->GetY()[islice]*windowRange;
  prSigmaMax = grprFreeKSSigma->GetY()[islice]+grprFreeKSSigma->GetY()[islice]*windowRange;
  prMeanMin  = grprFreeKSMean->GetY()[islice] -grprFreeKSMean->GetY()[islice]*windowRange;
  prMeanMax  = grprFreeKSMean->GetY()[islice] +grprFreeKSMean->GetY()[islice]*windowRange;
  prAmpMin   = grprFreeKSAmp->GetY()[islice]  -grprFreeKSAmp->GetY()[islice]*windowRange;
  prAmpMax   = grprFreeKSAmp->GetY()[islice]  +grprFreeKSAmp->GetY()[islice]*windowRange;

}
// -------------------------------------------------------------------------------------------------------
void GetExandSampleParameters(Int_t islice, TFile *sFile, Double_t *arrMean, Double_t *arrSigma, Double_t *cleanParams)
{

  //
  // Read Sample file
  //

  TTree *treeSample = (TTree*)sFile->Get("CleanSamples");

  Double_t elMean =0., elSigma = 0., elCleanMean = 0., elCleanSigma = 0.;
  Double_t piMean =0., piSigma = 0., piCleanMean = 0., piCleanSigma = 0.;
  Double_t kaMean =0., kaSigma = 0., kaCleanMean = 0., kaCleanSigma = 0.;
  Double_t prMean =0., prSigma = 0., prCleanMean = 0., prCleanSigma = 0.;

  Double_t elCleanSkew =0., elCleanKurtosis = 0., elCleanAmp = 0.;
  Double_t piCleanSkew =0., piCleanKurtosis = 0., piCleanAmp = 0.;
  Double_t kaCleanSkew =0., kaCleanKurtosis = 0., kaCleanAmp = 0.;
  Double_t prCleanSkew =0., prCleanKurtosis = 0., prCleanAmp = 0.;


  // get Clean Mean and expected mean-sigma from Clean Tree
  treeSample->SetBranchAddress("elMean" ,&elMean);
  treeSample->SetBranchAddress("piMean" ,&piMean);
  treeSample->SetBranchAddress("kaMean" ,&kaMean);
  treeSample->SetBranchAddress("prMean" ,&prMean);

  treeSample->SetBranchAddress("elSigma" ,&elSigma);
  treeSample->SetBranchAddress("piSigma" ,&piSigma);
  treeSample->SetBranchAddress("kaSigma" ,&kaSigma);
  treeSample->SetBranchAddress("prSigma" ,&prSigma);

  treeSample->SetBranchAddress("elCleanMean" ,&elCleanMean);
  treeSample->SetBranchAddress("piCleanMean" ,&piCleanMean);
  treeSample->SetBranchAddress("kaCleanMean" ,&kaCleanMean);
  treeSample->SetBranchAddress("prCleanMean" ,&prCleanMean);

  treeSample->SetBranchAddress("elCleanSigma" ,&elCleanSigma);
  treeSample->SetBranchAddress("piCleanSigma" ,&piCleanSigma);
  treeSample->SetBranchAddress("kaCleanSigma" ,&kaCleanSigma);
  treeSample->SetBranchAddress("prCleanSigma" ,&prCleanSigma);

  treeSample->SetBranchAddress("elCleanSkew" ,&elCleanSkew);
  treeSample->SetBranchAddress("piCleanSkew" ,&piCleanSkew);
  treeSample->SetBranchAddress("kaCleanSkew" ,&kaCleanSkew);
  treeSample->SetBranchAddress("prCleanSkew" ,&prCleanSkew);

  treeSample->SetBranchAddress("elCleanKurtosis" ,&elCleanKurtosis);
  treeSample->SetBranchAddress("piCleanKurtosis" ,&piCleanKurtosis);
  treeSample->SetBranchAddress("kaCleanKurtosis" ,&kaCleanKurtosis);
  treeSample->SetBranchAddress("prCleanKurtosis" ,&prCleanKurtosis);

  treeSample->SetBranchAddress("elCleanAmp" ,&elCleanAmp);
  treeSample->SetBranchAddress("piCleanAmp" ,&piCleanAmp);
  treeSample->SetBranchAddress("kaCleanAmp" ,&kaCleanAmp);
  treeSample->SetBranchAddress("prCleanAmp" ,&prCleanAmp);

  treeSample->GetEntry(islice);
  cleanParams[0] = elCleanMean;
  cleanParams[1] = piCleanMean;
  cleanParams[2] = kaCleanMean;
  cleanParams[3] = prCleanMean;

  cleanParams[4] = elCleanSigma;
  cleanParams[5] = piCleanSigma;
  cleanParams[6] = kaCleanSigma;
  cleanParams[7] = prCleanSigma;

  cleanParams[8] = elCleanSkew;
  cleanParams[9] = piCleanSkew;
  cleanParams[10] = kaCleanSkew;
  cleanParams[11] = prCleanSkew;

  cleanParams[12] = elCleanKurtosis;
  cleanParams[13] = piCleanKurtosis;
  cleanParams[14] = kaCleanKurtosis;
  cleanParams[15] = prCleanKurtosis;

  cleanParams[16] = elCleanAmp;
  cleanParams[17] = piCleanAmp;
  cleanParams[18] = kaCleanAmp;
  cleanParams[19] = prCleanAmp;




  arrMean[0] = elMean;  arrSigma[0] = elSigma;
  arrMean[1] = piMean;  arrSigma[1] = piSigma;
  arrMean[2] = kaMean;  arrSigma[2] = kaSigma;
  arrMean[3] = prMean;  arrSigma[3] = prSigma;

  // If there is a missing PIDres point the take the next slice's parameters
  for (Int_t i = islice; i<islice+10; i++) {
    treeSample->GetEntry(i);
    if (elMean!=0 && piMean!=0 && kaMean!=0 && prMean!=0 && elSigma!=0 && piSigma!=0 && kaSigma!=0 && prSigma!=0){
      arrMean[0] = elMean;  arrSigma[0] = elSigma;
      arrMean[1] = piMean;  arrSigma[1] = piSigma;
      arrMean[2] = kaMean;  arrSigma[2] = kaSigma;
      arrMean[3] = prMean;  arrSigma[3] = prSigma;
      break;
    }
  }

  // If there is a missing PIDres point and next parameters does not help use previous ones
  if (arrMean[0]==0||arrMean[1]==0||arrMean[2]==0||arrMean[3]==0||arrSigma[0]==0||arrSigma[1]==0||arrSigma[2]==0||arrMean[3]==0){
    for (Int_t i = islice; i>0; i--) {
      treeSample->GetEntry(i);
      if (arrMean[0]==0 || arrSigma[0]==0) { arrMean[0] = elMean;  arrSigma[0] = elSigma;}
      if (arrMean[1]==0 || arrSigma[1]==0) { arrMean[1] = elMean;  arrSigma[1] = elSigma;}
      if (arrMean[2]==0 || arrSigma[2]==0) { arrMean[2] = elMean;  arrSigma[2] = elSigma;}
      if (arrMean[3]==0 || arrSigma[3]==0) { arrMean[3] = elMean;  arrSigma[3] = elSigma;}
      if (elMean!=0 && piMean!=0 && kaMean!=0 && prMean!=0 && elSigma!=0 && piSigma!=0 && kaSigma!=0 && prSigma!=0) break;
    }
  }


}
// -------------------------------------------------------------------------------------------------------
void EstimateKS(Int_t iks, TH2D *h2DKS)
{

  //
  // apply 100 fits to FIT CLEAN samples with 10 MeV steps in between [0.4, 1.4]
  //

  TMinuit g;
  g.SetMaxIterations(10000);
  gMinuit->SetMaxIterations(10000);
  Double_t pMin = 0.4;
  Double_t pMax = 1.4;

  TH1D hCleanKurtosis[nParticle];
  TH1D hFitKurtosis[nParticle];
  TH1D hCleanSkewness[nParticle];
  TH1D hFitSkewness[nParticle];

  for (Int_t i = 0; i<nParticle; i++){
    hCleanKurtosis[i] = TH1D(Form("hCleanKurtosis_%d",i),Form("hCleanKurtosis_%d",i),100,1.5,3.5);
    hFitKurtosis[i]   = TH1D(Form("hFitKurtosis_%d"  ,i),Form("hFitKurtosis_%d"  ,i),100,1.5,3.5);
    hCleanSkewness[i] = TH1D(Form("hCleanSkewness_%d",i),Form("hCleanSkewness_%d",i),100,0.,2.);
    hFitSkewness[i]   = TH1D(Form("hFitSkewness_%d"  ,i),Form("hFitSkewness_%d"  ,i),100,0.,2.);
  }

  Int_t ptDownBin = 0;
  Int_t ptUpBin   = 0;
  Double_t ptStep = 0.;
  Double_t ptDown = 0.;
  Double_t ptUp   = 0.;
  Double_t pt     = 0.;
  for (Int_t islice = 0; islice<nBinsKSestimate; islice++){

    ptStep = ((pMax-pMin)/(Double_t)nBinsKSestimate);
    ptDown = pMin+islice*ptStep;
    ptUp   = ptDown+ptStep;
    pt     = (ptDown+ptUp)/2.;

    // find bins to be used in the projection

    ptDownBin = (h2D->GetXaxis()->FindBin(ptDown+0.0001));                                // TO FIX
    ptUpBin   = (h2D->GetXaxis()->FindBin(ptUp+0.0001));

    // initialise 1D histograms tobe fitted
    TH1D * h1CleanSamples[nParticle];
    TH1D * h1ExpectedMean[nParticle];
    TH1D * h1ExpectedSigma[nParticle];

     // Real particle parameters
    Double_t * arrSigma      = new Double_t[10];
    Double_t * arrMean       = new Double_t[10];
    Double_t * arrCleanSigma = new Double_t[10];
    Double_t * arrCleanMean  = new Double_t[10];
    Double_t * cleanParams   = new Double_t[30];

    // Clean Particle parameters
    Double_t * elParClean    = new Double_t[10];
    Double_t * piParClean    = new Double_t[10];
    Double_t * kaParClean    = new Double_t[10];
    Double_t * prParClean    = new Double_t[10];

    // Total fit parameters
    Double_t * pars          = new Double_t[31]; for (Int_t j=0; j<31; j++) pars[j] = 0;

    for (Int_t ipar=0; ipar<30; ipar++) cleanParams[ipar]  = 0;
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
    }

    // Get 1d slice
    TH1D * h1DKS = h2DKS->ProjectionY(Form("h1D_KS_%d",islice),ptDownBin,ptDownBin);
    h1DKS->Rebin(rebinFactor); h1DKS->Scale(1./rebinFactor);
    Double_t binWidth = h1DKS->GetXaxis()->GetBinWidth(50);
    h1DKS->Scale(1./binWidth);

    // Get Clean Samples 1d slice
    h1CleanSamples[0]  = GetClean1DSlice(hCleanSamples[0],"h1Electron_slice_%f" ,ptDownBin,ptUpBin);
    h1CleanSamples[1]  = GetClean1DSlice(hCleanSamples[1],"h1Pion_slice_%f"     ,ptDownBin,ptUpBin);
    h1CleanSamples[2]  = GetClean1DSlice(hCleanSamples[2],"h1Kaon_slice_%f"     ,ptDownBin,ptUpBin);
    h1CleanSamples[3]  = GetClean1DSlice(hCleanSamples[3],"h1Proton_slice_%f"   ,ptDownBin,ptUpBin);

    // Get Expected 1d Slice
    h1ExpectedMean[0]  = GetClean1DSlice(hExpectedMean[0],"h1ElectronExMean_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedMean[1]  = GetClean1DSlice(hExpectedMean[1],"h1PionExMean_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedMean[2]  = GetClean1DSlice(hExpectedMean[2],"h1KaonExMean_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedMean[3]  = GetClean1DSlice(hExpectedMean[3],"h1ProtonExMean_slice_%f"  ,ptDownBin,ptUpBin);

    h1ExpectedSigma[0] = GetClean1DSlice(hExpectedSigma[0],"h1ElectronExSigma_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedSigma[1] = GetClean1DSlice(hExpectedSigma[1],"h1PionExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedSigma[2] = GetClean1DSlice(hExpectedSigma[2],"h1KaonExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedSigma[3] = GetClean1DSlice(hExpectedSigma[3],"h1ProtonExSigma_slice_%f"  ,ptDownBin,ptUpBin);

    // Get Expected Mean and Sigma
    arrMean[0]  = (exWRTMean) ? h1ExpectedMean[0]->GetMean(): h1ExpectedMean[0]->GetMaximum();
    arrMean[1]  = (exWRTMean) ? h1ExpectedMean[1]->GetMean(): h1ExpectedMean[1]->GetMaximum();
    arrMean[2]  = (exWRTMean) ? h1ExpectedMean[2]->GetMean(): h1ExpectedMean[2]->GetMaximum();
    arrMean[3]  = (exWRTMean) ? h1ExpectedMean[3]->GetMean(): h1ExpectedMean[3]->GetMaximum();
    arrSigma[0] = (exWRTMean) ? h1ExpectedSigma[0]->GetMean(): h1ExpectedSigma[0]->GetMaximum();
    arrSigma[1] = (exWRTMean) ? h1ExpectedSigma[1]->GetMean(): h1ExpectedSigma[1]->GetMaximum();
    arrSigma[2] = (exWRTMean) ? h1ExpectedSigma[2]->GetMean(): h1ExpectedSigma[2]->GetMaximum();
    arrSigma[3] = (exWRTMean) ? h1ExpectedSigma[3]->GetMean(): h1ExpectedSigma[3]->GetMaximum();

    //  Make the fits
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[0] , elParClean, arrMean, arrSigma, kElectron);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[1] , piParClean, arrMean, arrSigma, kPion    );
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[2] , kaParClean, arrMean, arrSigma, kKaon    );
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[3] , prParClean, arrMean, arrSigma, kProton  );
    Fit1DSliceForKSEstimate(iks, h1DKS, pars, arrMean, arrSigma);

    cout.setf(ios::fixed);
//     cout << " ******************************************************************************************* " << endl;
    cout << " iks = "            << setprecision(2) << iks            << setw(12);
    cout << " ptDownBin = "      << setprecision(2) << ptDownBin      << setw(12);
    cout << " ptUpBin = "        << setprecision(2) << ptUpBin        << setw(12);
    cout << " p = "              << setprecision(3) << pt             << setw(12);
    cout << " islice = "         << setprecision(2) << islice         << setw(12);
    cout << " chi2 = "           << setprecision(2) << pars[20]       << setw(12)   << endl;
    cout << " =================================================================== " << endl;
    cout << " *** piKClean = "   << setprecision(2) << piParClean[3]  << setw(12);
    cout << " kaKClean = "       << setprecision(2) << kaParClean[3]  << setw(12);
    cout << " prKClean = "       << setprecision(2) << prParClean[3]  << setw(12);

    cout << " -------- piKFit = "<< setprecision(2) << pars[8]        << setw(12);
    cout << " kaKFit = "         << setprecision(2) << pars[13]       << setw(12);
    cout << " prKFit = "         << setprecision(2) << pars[18]       << setw(12)   << endl;

    cout << " *** piSClean = "   << setprecision(2) << piParClean[4]  << setw(12);
    cout << " kaSClean = "       << setprecision(2) << kaParClean[4]  << setw(12);
    cout << " prSClean = "       << setprecision(2) << prParClean[4]  << setw(12);

    cout << " -------- piSFit = "<< setprecision(2) << pars[9]        << setw(12);
    cout << " kaSFit = "         << setprecision(2) << pars[14]       << setw(12);
    cout << " prSFit = "         << setprecision(2) << pars[19]       << setw(12) << endl;

    cout << " =*=*=* piCleanChi2 = "<< setprecision(2) << piParClean[6]  << setw(12);
    cout << "     kaCleanChi2 = "   << setprecision(2) << kaParClean[6]  << setw(12);
    cout << "     prCleanChi2 = "   << setprecision(2) << prParClean[6]  << setw(12) << endl;
    cout << " ******************************************************************************************* " << endl;

    if (iks==0){

    // Safe region Kurtosis
      if (pt>=piCleanSafeMin && pt<=piCleanSafeMax) hCleanKurtosis[kPion]  .Fill(piParClean[3]);
      if (pt>=kaCleanSafeMin && pt<=kaCleanSafeMax) hCleanKurtosis[kKaon]  .Fill(kaParClean[3]);
      if (pt>=prCleanSafeMin && pt<=prCleanSafeMax) hCleanKurtosis[kProton].Fill(prParClean[3]);

      if (pt>=piFitSafeMin && pt<=piFitSafeMax) hFitKurtosis[kPion]  .Fill(pars[8]);
      if (pt>=kaFitSafeMin && pt<=kaFitSafeMax) hFitKurtosis[kKaon]  .Fill(pars[13]);
      if (pt>=prFitSafeMin && pt<=prFitSafeMax) hFitKurtosis[kProton].Fill(pars[18]);

    } else {

    // Safe region Skewness
      if (pt>=piCleanSafeMin && pt<=piCleanSafeMax) hCleanSkewness[kPion]  .Fill(piParClean[4]);
      if (pt>=kaCleanSafeMin && pt<=kaCleanSafeMax) hCleanSkewness[kKaon]  .Fill(kaParClean[4]);
      if (pt>=prCleanSafeMin && pt<=prCleanSafeMax) hCleanSkewness[kProton].Fill(prParClean[4]);

      if (pt>=piFitSafeMin && pt<=piFitSafeMax) hFitSkewness[kPion]  .Fill(pars[9]);
      if (pt>=kaFitSafeMin && pt<=kaFitSafeMax) hFitSkewness[kKaon]  .Fill(pars[14]);
      if (pt>=prFitSafeMin && pt<=prFitSafeMax) hFitSkewness[kProton].Fill(pars[19]);

    }

    // dump all info to tree
    debugFile -> GetFile()->cd();
    *debugFile << "ks" <<

        "iks="          << iks                <<
        "p="            << pt                 <<
        "eta="          << fEtaDown             <<
        "cent="         << fCentDown            <<

        "elMean="       << arrMean[0]         <<
        "piMean="       << arrMean[1]         <<
        "kaMean="       << arrMean[2]         <<
        "prMean="       << arrMean[3]         <<

        "elSigma="      << arrSigma[0]        <<
        "piSigma="      << arrSigma[1]        <<
        "kaSigma="      << arrSigma[2]        <<
        "prSigma="      << arrSigma[3]        ;

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

        "totChi2="         << pars[20]       <<

        "\n";

    delete [] arrSigma;
    delete [] arrMean;
    delete [] arrCleanSigma;
    delete [] arrCleanMean;
    delete [] cleanParams;

    delete [] piParClean;
    delete [] elParClean;
    delete [] kaParClean;
    delete [] prParClean;
    delete [] pars;

  }

  // Set the kurtosis and skewness to be used in further iterations
  if (iks==0){

    kurtosisFixClean[kPion]   = hCleanKurtosis[kPion]  .GetMean();
    kurtosisFixClean[kKaon]   = hCleanKurtosis[kKaon]  .GetMean();
    kurtosisFixClean[kProton] = hCleanKurtosis[kProton].GetMean();

    kurtosisFixFit[kPion]     = hFitKurtosis[kPion].GetMean();
    kurtosisFixFit[kKaon]     = hFitKurtosis[kKaon].GetMean();
    kurtosisFixFit[kProton]   = hFitKurtosis[kProton].GetMean();

  }

  // dump final fixed KS into debug file
  debugFile -> GetFile()->cd();
  *debugFile << "fixedKS" <<

      "iks="             << iks                         <<
      "eta="             << fEtaDown                      <<
      "cent="            << fCentDown                     <<

      "FitSkewPi="       << skewnessFixFit[kPion]       <<
      "FitSkewKa="       << skewnessFixFit[kKaon]       <<
      "FitSkewPr="       << skewnessFixFit[kProton]     <<

      "CleanSkewPi="     << skewnessFixClean[kPion]     <<
      "CleanSkewKa="     << skewnessFixClean[kKaon]     <<
      "CleanSkewPr="     << skewnessFixClean[kProton]   <<

      "FitKurtPi="       << kurtosisFixFit[kPion]       <<
      "FitKurtKa="       << kurtosisFixFit[kKaon]       <<
      "FitKurtPr="       << kurtosisFixFit[kProton]     <<

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
  }

  if (iks==1) {

    // decide whether to let automaticly found KS or by eye
    if (automaticKS) {
      skewnessFixFit[kPion]     = hFitSkewness[kPion]  .GetMean()*piSkewnessScan;  // multiplied with piSkewnessScan since skewness decrease with p
      skewnessFixFit[kKaon]     = hFitSkewness[kKaon]  .GetMean()*kaSkewnessScan;
      skewnessFixFit[kProton]   = hFitSkewness[kProton].GetMean()*prSkewnessScan;

      skewnessFixClean[kPion]   = hCleanSkewness[kPion]  .GetMean()*piSkewnessScan;
      skewnessFixClean[kKaon]   = hCleanSkewness[kKaon]  .GetMean()*kaSkewnessScan;
      skewnessFixClean[kProton] = hCleanSkewness[kProton].GetMean()*prSkewnessScan;
    } else {
      skewnessFixFit[kPion]     = TMath::Min(hFitSkewness[kPion].GetMean(),piSkewnessScan);
      skewnessFixFit[kKaon]     = TMath::Min(hFitSkewness[kKaon].GetMean(),kaSkewnessScan);
      skewnessFixFit[kProton]   = TMath::Min(hFitSkewness[kProton].GetMean(),prSkewnessScan);

      skewnessFixClean[kPion]   = TMath::Min(hCleanSkewness[kPion].GetMean(),piSkewnessScan);
      skewnessFixClean[kKaon]   = TMath::Min(hCleanSkewness[kKaon].GetMean(),kaSkewnessScan);
      skewnessFixClean[kProton] = TMath::Min(hCleanSkewness[kProton].GetMean(),prSkewnessScan);
    }

  }

  // Set the final kurtosis and skewness params --> take the estimation form clean samples
  if (iks==1){

    if (useSafeRegion){
      skewnessFix[kPion]     = skewnessFixFit[kPion];
      skewnessFix[kKaon]     = skewnessFixClean[kKaon];
      skewnessFix[kProton]   = skewnessFixFit[kProton];
      kurtosisFix[kPion]     = kurtosisFixFit[kPion]*piKurtosisScan;
      kurtosisFix[kKaon]     = kurtosisFixClean[kKaon]*kaKurtosisScan;
      kurtosisFix[kProton]   = kurtosisFixFit[kProton]*prKurtosisScan;
    } else {
      skewnessFix[kPion]     = skewnessFixClean[kPion];
      skewnessFix[kKaon]     = skewnessFixClean[kKaon];
      skewnessFix[kProton]   = skewnessFixClean[kProton];
      kurtosisFix[kPion]     = kurtosisFixClean[kPion];
      kurtosisFix[kKaon]     = kurtosisFixClean[kKaon];
      kurtosisFix[kProton]   = kurtosisFixClean[kProton];
    }

  }

}
// -------------------------------------------------------------------------------------------------------
void FitAllSamples(const Int_t nSlice, const Double_t ptMin, const Double_t ptMax, TString sampleFileName)
{

  //
  // Fit only samples
  //

  TMinuit g;
  g.SetMaxIterations(10000);
  gMinuit->SetMaxIterations(10000);
  sampleFile   = new TTreeSRedirector(sampleFileName,"recreate");
  cleanResArr  . SetOwner(kTRUE);
  cleanResArrFreeKS.SetOwner(kTRUE);


  Int_t ptDownBin = 0;
  Int_t ptUpBin   = 0;
  Double_t ptStep = 0.;
  Double_t ptDown = 0.;
  Double_t ptUp   = 0.;
  Double_t pt     = 0.;
  for (Int_t islice = 0; islice<nSlice; islice++){

    if (islice%10 == 0) cout << " slice = " << islice << endl;

    ptStep = ((ptMax-ptMin)/(Double_t)nSlice);
    ptDown = ptMin+islice*ptStep;
    ptUp   = ptDown+ptStep;
    pt     = (ptDown+ptUp)/2.;

    // find bins to be used in the projection
    ptDownBin = (h2D->GetXaxis()->FindBin(ptDown+0.001));                                // TO FIX
    ptUpBin   = (h2D->GetXaxis()->FindBin(ptUp+0.001));


    // Storage for the histograms
    TObjArray * cleanSampArr = new TObjArray(6);  cleanSampArr -> SetOwner(kTRUE);
    cleanSampArr->SetName(Form("cleanSamples_%d_%4.3f",islice,pt));

    TObjArray * cleanSampArrFreeKS = new TObjArray(6);  cleanSampArrFreeKS -> SetOwner(kTRUE);
    cleanSampArrFreeKS->SetName(Form("cleanSamplesFreeKS_%d_%4.3f",islice,pt));

    TH1D * h1CleanSamples[nParticle];
    TH1D * h1CleanSamplesFreeKS[nParticle];
    TH1D * h1ExpectedMean[nParticle];
    TH1D * h1ExpectedSigma[nParticle];

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

    Double_t * elParCleanFreeKS    = new Double_t[10];
    Double_t * piParCleanFreeKS    = new Double_t[10];
    Double_t * kaParCleanFreeKS    = new Double_t[10];
    Double_t * prParCleanFreeKS    = new Double_t[10];

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
      piParCleanFreeKS[ipar] = 0;
      kaParCleanFreeKS[ipar] = 0;
      prParCleanFreeKS[ipar] = 0;
      elParCleanFreeKS[ipar] = 0;
    }

    // Get Clean Samples 1d slice
    h1CleanSamples[0] = GetClean1DSlice(hCleanSamples[0],"h1Electron_slice_%f" ,ptDownBin,ptUpBin);
    h1CleanSamples[1] = GetClean1DSlice(hCleanSamples[1],"h1Pion_slice_%f"     ,ptDownBin,ptUpBin);
    h1CleanSamples[2] = GetClean1DSlice(hCleanSamples[2],"h1Kaon_slice_%f"     ,ptDownBin,ptUpBin);
    h1CleanSamples[3] = GetClean1DSlice(hCleanSamples[3],"h1Proton_slice_%f"   ,ptDownBin,ptUpBin);

    h1CleanSamplesFreeKS[0] = (TH1D*)h1CleanSamples[0]->Clone();
    h1CleanSamplesFreeKS[1] = (TH1D*)h1CleanSamples[1]->Clone();
    h1CleanSamplesFreeKS[2] = (TH1D*)h1CleanSamples[2]->Clone();
    h1CleanSamplesFreeKS[3] = (TH1D*)h1CleanSamples[3]->Clone();


    // Get Expected 1d Slice
    h1ExpectedMean[0]   = GetClean1DSlice(hExpectedMean[0],"h1ElectronExMean_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedMean[1]   = GetClean1DSlice(hExpectedMean[1],"h1PionExMean_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedMean[2]   = GetClean1DSlice(hExpectedMean[2],"h1KaonExMean_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedMean[3]   = GetClean1DSlice(hExpectedMean[3],"h1ProtonExMean_slice_%f"  ,ptDownBin,ptUpBin);

    h1ExpectedSigma[0]   = GetClean1DSlice(hExpectedSigma[0],"h1ElectronExSigma_slice_%f",ptDownBin,ptUpBin);
    h1ExpectedSigma[1]   = GetClean1DSlice(hExpectedSigma[1],"h1PionExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedSigma[2]   = GetClean1DSlice(hExpectedSigma[2],"h1KaonExSigma_slice_%f"    ,ptDownBin,ptUpBin);
    h1ExpectedSigma[3]   = GetClean1DSlice(hExpectedSigma[3],"h1ProtonExSigma_slice_%f"  ,ptDownBin,ptUpBin);

    // Get Expected Mean and Sigma
    arrMean[0]  = (exWRTMean) ? h1ExpectedMean[0]->GetMean(): h1ExpectedMean[0]->GetMaximum();
    arrMean[1]  = (exWRTMean) ? h1ExpectedMean[1]->GetMean(): h1ExpectedMean[1]->GetMaximum();
    arrMean[2]  = (exWRTMean) ? h1ExpectedMean[2]->GetMean(): h1ExpectedMean[2]->GetMaximum();
    arrMean[3]  = (exWRTMean) ? h1ExpectedMean[3]->GetMean(): h1ExpectedMean[3]->GetMaximum();
    arrSigma[0] = (exWRTMean) ? h1ExpectedSigma[0]->GetMean(): h1ExpectedSigma[0]->GetMaximum();
    arrSigma[1] = (exWRTMean) ? h1ExpectedSigma[1]->GetMean(): h1ExpectedSigma[1]->GetMaximum();
    arrSigma[2] = (exWRTMean) ? h1ExpectedSigma[2]->GetMean(): h1ExpectedSigma[2]->GetMaximum();
    arrSigma[3] = (exWRTMean) ? h1ExpectedSigma[3]->GetMean(): h1ExpectedSigma[3]->GetMaximum();

    //  Make the fits
    TMatrixD *elCleanMatrix = FitParticleSample(h1CleanSamples[0] , elParClean, arrMean, arrSigma, kElectron);
    TMatrixD *piCleanMatrix = FitParticleSample(h1CleanSamples[1] , piParClean, arrMean, arrSigma, kPion    );
    TMatrixD *kaCleanMatrix = FitParticleSample(h1CleanSamples[2] , kaParClean, arrMean, arrSigma, kKaon    );
    TMatrixD *prCleanMatrix = FitParticleSample(h1CleanSamples[3] , prParClean, arrMean, arrSigma, kProton  );

    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[0] , elParCleanFreeKS, arrMean, arrSigma, kElectron);
    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[1] , piParCleanFreeKS, arrMean, arrSigma, kPion    );
    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[2] , kaParCleanFreeKS, arrMean, arrSigma, kKaon    );
    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[3] , prParCleanFreeKS, arrMean, arrSigma, kProton  );

    GetCleanExParams(h1CleanSamples[0],arrSigma,arrCleanSigma,arrCleanMean,kElectron);
    GetCleanExParams(h1CleanSamples[1],arrSigma,arrCleanSigma,arrCleanMean,kPion);
    GetCleanExParams(h1CleanSamples[2],arrSigma,arrCleanSigma,arrCleanMean,kKaon);
    GetCleanExParams(h1CleanSamples[3],arrSigma,arrCleanSigma,arrCleanMean,kProton);

    // Dump histograms into arrays
    h1CleanSamples[0] -> SetName(Form("CleanElectron_%d_%4.3f" ,islice,pt));
    h1CleanSamples[1] -> SetName(Form("CleanPion_%d_%4.3f"     ,islice,pt));
    h1CleanSamples[2] -> SetName(Form("CleanKaon_%d_%4.3f"     ,islice,pt));
    h1CleanSamples[3] -> SetName(Form("CleanProton_%d_%4.3f"   ,islice,pt));

    h1CleanSamplesFreeKS[0] -> SetName(Form("CleanElectronFreeKS_%d_%4.3f" ,islice,pt));
    h1CleanSamplesFreeKS[1] -> SetName(Form("CleanPionFreeKS_%d_%4.3f"     ,islice,pt));
    h1CleanSamplesFreeKS[2] -> SetName(Form("CleanKaonFreeKS_%d_%4.3f"     ,islice,pt));
    h1CleanSamplesFreeKS[3] -> SetName(Form("CleanProtonFreeKS_%d_%4.3f"   ,islice,pt));


    cleanSampArr -> AddAt(h1CleanSamples[0],0);
    cleanSampArr -> AddAt(h1CleanSamples[1],1);
    cleanSampArr -> AddAt(h1CleanSamples[2],2);
    cleanSampArr -> AddAt(h1CleanSamples[3],3);

    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[0],0);
    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[1],1);
    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[2],2);
    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[3],3);

    cleanResArr      . AddAt(cleanSampArr,islice);
    cleanResArrFreeKS. AddAt(cleanSampArrFreeKS,islice);

    // Clean Samples dump tree
    sampleFile->GetFile()->cd();
    *sampleFile << "CleanSamples" <<

        "slice="        << islice             <<
        "p="            << pt                 <<
        "eta="          << fEtaDown             <<
        "cent="         << fCentDown            <<

        "elMean="       << arrMean[0]         <<
        "piMean="       << arrMean[1]         <<
        "kaMean="       << arrMean[2]         <<
        "prMean="       << arrMean[3]         <<

        "elSigma="      << arrSigma[0]        <<
        "piSigma="      << arrSigma[1]        <<
        "kaSigma="      << arrSigma[2]        <<
        "prSigma="      << arrSigma[3]        ;

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

        "elCleanAmp="      << elParClean[0]    <<
        "elCleanMean="     << elParClean[1]    <<
        "elCleanSigma="    << elParClean[2]    <<
        "elCleanKurtosis=" << elParClean[3]    <<
        "elCleanSkew="     << elParClean[4]    <<
        "elCleanInt="      << elParClean[5]    <<
        "elCleanChi2="     << elParClean[6]    <<
        "elCleanMatrix.="  << elCleanMatrix    <<

        "piCleanAmp="      << piParClean[0]    <<
        "piCleanMean="     << piParClean[1]    <<
        "piCleanSigma="    << piParClean[2]    <<
        "piCleanKurtosis=" << piParClean[3]    <<
        "piCleanSkew="     << piParClean[4]    <<
        "piCleanInt="      << piParClean[5]    <<
        "piCleanChi2="     << piParClean[6]    <<
        "piCleanMatrix.="  << piCleanMatrix    ;

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
        "prCleanMatrix.="  << prCleanMatrix    ;

	*sampleFile << "CleanSamples" <<

        "elCleanAmpFreeKS="      << elParCleanFreeKS[0]    <<
        "elCleanMeanFreeKS="     << elParCleanFreeKS[1]    <<
        "elCleanSigmaFreeKS="    << elParCleanFreeKS[2]    <<
        "elCleanKurtosisFreeKS=" << elParCleanFreeKS[3]    <<
        "elCleanSkewFreeKS="     << elParCleanFreeKS[4]    <<
        "elCleanIntFreeKS="      << elParCleanFreeKS[5]    <<
        "elCleanChi2FreeKS="     << elParCleanFreeKS[6]    <<

        "piCleanAmpFreeKS="      << piParCleanFreeKS[0]    <<
        "piCleanMeanFreeKS="     << piParCleanFreeKS[1]    <<
        "piCleanSigmaFreeKS="    << piParCleanFreeKS[2]    <<
        "piCleanKurtosisFreeKS=" << piParCleanFreeKS[3]    <<
        "piCleanSkewFreeKS="     << piParCleanFreeKS[4]    <<
        "piCleanIntFreeKS="      << piParCleanFreeKS[5]    <<
        "piCleanChi2FreeKS="     << piParCleanFreeKS[6]    ;

	*sampleFile << "CleanSamples" <<

        "kaCleanAmpFreeKS="      << kaParCleanFreeKS[0]    <<
        "kaCleanMeanFreeKS="     << kaParCleanFreeKS[1]    <<
        "kaCleanSigmaFreeKS="    << kaParCleanFreeKS[2]    <<
        "kaCleanKurtosisFreeKS=" << kaParCleanFreeKS[3]    <<
        "kaCleanSkewFreeKS="     << kaParCleanFreeKS[4]    <<
        "kaCleanIntFreeKS="      << kaParCleanFreeKS[5]    <<
        "kaCleanChi2FreeKS="     << kaParCleanFreeKS[6]    <<

        "prCleanAmpFreeKS="      << prParCleanFreeKS[0]    <<
        "prCleanMeanFreeKS="     << prParCleanFreeKS[1]    <<
        "prCleanSigmaFreeKS="    << prParCleanFreeKS[2]    <<
        "prCleanKurtosisFreeKS=" << prParCleanFreeKS[3]    <<
        "prCleanSkewFreeKS="     << prParCleanFreeKS[4]    <<
        "prCleanIntFreeKS="      << prParCleanFreeKS[5]    <<
        "prCleanChi2FreeKS="     << prParCleanFreeKS[6]    <<

        "\n";

    delete [] arrSigma;
    delete [] arrMean;
    delete [] arrCleanSigma;
    delete [] arrCleanMean;
    delete [] piParClean;
    delete [] elParClean;
    delete [] kaParClean;
    delete [] prParClean;
    delete [] piParCleanFreeKS;
    delete [] elParCleanFreeKS;
    delete [] kaParCleanFreeKS;
    delete [] prParCleanFreeKS;

    delete elCleanMatrix;
    delete piCleanMatrix;
    delete kaCleanMatrix;
    delete prCleanMatrix;

  }

  sampleFile -> GetFile()->cd();
  cleanResArr       . Write("cleanResArr" ,TObject::kSingleKey);
  cleanResArrFreeKS . Write("cleanResArrFreeKS" ,TObject::kSingleKey);

  // delete arrays
  cleanResArr  . Delete();
  cleanResArrFreeKS . Delete();
  delete sampleFile;

}
// -------------------------------------------------------------------------------------------------------
void IterativeFitting(Int_t iIter, const Int_t nSlice, const Double_t ptMin, const Double_t ptMax, TString fileFit, TString readSmooth, TString readSamp, TString lineShapesFile)
{

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
  Int_t allBins = nSlice;

  TMinuit g;
  g.SetMaxIterations(100000);
  gMinuit->SetMaxIterations(100000);
  TString histFileName  = Form("Hists_Iteration%d_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",iIter,fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown);
  histsFile = new TTreeSRedirector(histFileName,"recreate");
  fitFile   = new TTreeSRedirector(fileFit     ,"recreate");
  lineShapesLookUp = new TTreeSRedirector(lineShapesFile     ,"recreate");

  // OutPut Arrays
  fitResArr    . SetOwner(kTRUE);
  fitResiduals . SetOwner(kTRUE);
  histLineShapesCArr.SetOwner(kTRUE);
  funcLineShapesCArr.SetOwner(kTRUE);
  funcLineShapesObjArr.SetOwner(kTRUE);


    // loop over slices
  Double_t elIntegralFit[allBins], elFitAmp[allBins], elFitMean[allBins], elFitSigma[allBins];
  Double_t piIntegralFit[allBins], piFitAmp[allBins], piFitMean[allBins], piFitSigma[allBins];
  Double_t kaIntegralFit[allBins], kaFitAmp[allBins], kaFitMean[allBins], kaFitSigma[allBins];
  Double_t prIntegralFit[allBins], prFitAmp[allBins], prFitMean[allBins], prFitSigma[allBins];

  Double_t elFitKurtosis[allBins], elFitSkew[allBins], elChi2Fit[allBins];
  Double_t piFitKurtosis[allBins], piFitSkew[allBins], piChi2Fit[allBins];
  Double_t kaFitKurtosis[allBins], kaFitSkew[allBins], kaChi2Fit[allBins];
  Double_t prFitKurtosis[allBins], prFitSkew[allBins], prChi2Fit[allBins];

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
    ptDownBin = (h2D->GetXaxis()->FindBin(ptDown+0.001));                                // TO FIX
    ptUpBin   = (h2D->GetXaxis()->FindBin(ptUp+0.001));
    momBin = ptUpBin;

    cout << " centBin etaBin momBin signBin " << centBin << "  " << etaBin << " " <<  momBin << "  " <<  signBin << endl;

    // main 1d slice used in the fitting
    TH1D * h1D = h2D->ProjectionY(Form("h1D_slice_%d",islice),ptDownBin,ptDownBin);
    h1D->Rebin(rebinFactor); h1D->Scale(1./rebinFactor);
    Double_t binWidth = h1D->GetXaxis()->GetBinWidth(50);
    h1D->Scale(1./binWidth);
    TH1D * hResidual = (TH1D*)h1D->Clone();
    TH1D * hLineShape = (TH1D*)h1D->Clone();

    Double_t * arrSigma      = new Double_t[10];
    Double_t * arrMean       = new Double_t[10];
    Double_t * arrMeanWindow = new Double_t[10];
    Double_t * cleanParams   = new Double_t[30];

    for (Int_t ipar=0; ipar<10; ipar++)
    {
      arrMean[ipar]       = 0;
      arrSigma[ipar]      = 0;
      arrMeanWindow[ipar] = 0;
    }
    for (Int_t ipar=0; ipar<30; ipar++) cleanParams[ipar]  = 0;


    GetExandSampleParameters(islice,sampFile,arrMean,arrSigma,cleanParams);
    if (dump) ApplyInitialFits(islice,pt,h1D,arrSigma,arrMean);

    // Fit part
    maxBin = h1D->GetBinContent(h1D->GetMaximumBin());
    Double_t fitWinMean  = 1.2;

    elMeanMin = TMath::Max(30.,(arrMean[kElectron]-2.*arrSigma[kElectron]));                                   // electron
    piMeanMin = TMath::Max(30.,(arrMean[kPion]-fitWinMean*arrSigma[kPion]));                           // pion
    kaMeanMin = TMath::Max(30.,(arrMean[kKaon]-fitWinMean*arrSigma[kKaon]));                           // kaon
    prMeanMin = TMath::Max(30.,(arrMean[kProton]-fitWinMean*arrSigma[kProton]));                           // proton

    elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kElectron]+2.*arrSigma[kElectron]));
    piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kPion]+fitWinMean*arrSigma[kPion]));
    kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kKaon]+fitWinMean*arrSigma[kKaon]));
    prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kProton]+fitWinMean*arrSigma[kProton]));

    if (pt<0.3) {
      prMeanMin = TMath::Max(30.              ,(arrMean[kProton]-arrSigma[kProton]*3.));
      prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kProton]+arrSigma[kProton]*3.));
    }

    // Parameters: 0 --> Amplitude, 1 --> Mean, 2 --> Sigma, 3 --> kurtosis, 4 --> skewness
    //     TF1 *g1 = new TF1("g1",fitFunctionGenGaus,dEdxMin,dEdxMax);        //g1->SetNpx(nBinsInLookUpTable);
    //     TF1 *g2 = new TF1("g2",fitFunctionGenGaus,dEdxMin,dEdxMax);        //g2->SetNpx(nBinsInLookUpTable);
    //     TF1 *g3 = new TF1("g3",fitFunctionGenGaus,dEdxMin,dEdxMax);        //g3->SetNpx(nBinsInLookUpTable);
    //     TF1 *g4 = new TF1("g4",fitFunctionGenGaus,dEdxMin,dEdxMax);        //g4->SetNpx(nBinsInLookUpTable);

    TF1 *g1    = new TF1("g1",fitFunctionGenGaus,elMeanMin,elMeanMax);        //g1->SetNpx(nBinsInLookUpTable);
    TF1 *g2    = new TF1("g2",fitFunctionGenGaus,piMeanMin,piMeanMax);        //g2->SetNpx(nBinsInLookUpTable);
    TF1 *g3    = new TF1("g3",fitFunctionGenGaus,kaMeanMin,kaMeanMax);        //g3->SetNpx(nBinsInLookUpTable);
    TF1 *g4    = new TF1("g4",fitFunctionGenGaus,prMeanMin,prMeanMax);        //g4->SetNpx(nBinsInLookUpTable);


    g1->SetParameters(maxBin,arrMean[kElectron],arrSigma[kElectron],2.,0.);
    g2->SetParameters(maxBin,arrMean[kPion],arrSigma[kPion],2.,0.);
    g3->SetParameters(maxBin,arrMean[kKaon],arrSigma[kKaon],2.,0.);
    g4->SetParameters(maxBin,arrMean[kProton],arrSigma[kProton],2.,0.);

    g1->SetLineColor(4);          // Blue
    g2->SetLineColor(kOrange-1);    // Black
    g3->SetLineColor(kGreen+2);   // Green
    g4->SetLineColor(kMagenta+1);          // Orange

    // restrict the parameters of gaus functions of individual particles
    g1->SetParLimits(0,0.,maxBin);
    g2->SetParLimits(0,0.,maxBin);
    g3->SetParLimits(0,0.,maxBin);
    g4->SetParLimits(0,0.,maxBin);

    g1->SetParLimits(1,elMeanMin,elMeanMax);
    g2->SetParLimits(1,piMeanMin,piMeanMax);
    g3->SetParLimits(1,kaMeanMin,kaMeanMax);
    g4->SetParLimits(1,prMeanMin,prMeanMax);

    g1->SetParLimits(2,arrSigma[0]/5.,arrSigma[0]*2.);
    g2->SetParLimits(2,arrSigma[1]/5.,arrSigma[1]*2.);
    g3->SetParLimits(2,arrSigma[2]/5.,arrSigma[2]*2.);
    g4->SetParLimits(2,arrSigma[3]/5.,arrSigma[3]*2.);

    g1->SetParLimits(3,kurtosisMin,kurtosisMax);
    g2->SetParLimits(3,kurtosisMin,kurtosisMax);
    g3->SetParLimits(3,kurtosisMin,kurtosisMax);
    g4->SetParLimits(3,kurtosisMin,kurtosisMax);

    g1->SetParLimits(4,skewMin,skewMax);
    g2->SetParLimits(4,skewMin,skewMax);
    g3->SetParLimits(4,skewMin,skewMax);
    g4->SetParLimits(4,skewMin,skewMax);

    if (fixedK || fixedS) SetFixedKSparameterForIndividualFits(g1,g2,g3,g4);

    cout << " ============================================ " << endl;
    cout << " ElMean  = "   << arrMean[0]          ;
    cout << " ElSigma = "   << arrSigma[0]          << endl;
    cout << " PiMean  = "   << arrMean[1]          ;
    cout << " PiSigma = "   << arrSigma[1]          << endl;
    cout << " KaMean  = "   << arrMean[2]          ;
    cout << " KaSigma = "   << arrSigma[2]          << endl;
    cout << " PrMean  = "   << arrMean[3]          ;
    cout << " PrSigma = "   << arrSigma[3]          << endl;


    // Clone histogram for total fit
    TH1D *hClone = (TH1D*)h1D->Clone();
    SetStyleFor1DSlice(hClone);
    SetStyleFor1DSlice(h1D);

        // Fit h1D for individual particles
    if (arrMean[0]>=10) h1D->Fit(g1,"QNR");
    if (arrMean[1]>=10) h1D->Fit(g2,"QNR+");
    if (arrMean[2]>=10) h1D->Fit(g3,"QNR+");
    if (arrMean[3]>=10) h1D->Fit(g4,"QNR+");

    // Get Fit parameters from the individual fits and set for the total fit
    Double_t par[20] = {0};
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[5]);
    g3->GetParameters(&par[10]);
    g4->GetParameters(&par[15]);

    TF1 *total = new TF1("total","g1+g2+g3+g4",dEdxMin,dEdxMax); total->SetLineWidth(3); total->SetLineColor(2); total->SetParameters(par);
    // total->SetNpx(nBinsInLookUpTable);

    // Some setters for the total fit
    SetTotalFitParameters(iIter,islice,nSlice,total,arrMean,arrSigma,arrMeanWindow,cleanParams,h1D,pt);
    SetParNamesOfTotalFit(total);
    CheckIfMeanIsZero(total,arrMean);

    // Apply total fit
    // gStyle->SetOptFit(1111);
    h1D->Fit(total,"QMR");
    hClone->Fit(total,"QMR+");

    // Add total function to legend of hClone
    TLegend *l = new TLegend(0.75, 0.75, 0.85, 0.85);
    l->SetTextFont(62);
    l->SetTextSize(0.03);
    l->SetFillColor(0);
    l->AddEntry(hClone," Data","LPE");
    l->AddEntry(total," Total Fit","L");
    hClone->GetListOfFunctions()->Add(l);
    if (dump) {debugFile->GetFile()->cd(); hClone->Write(Form("TotalFit_%d",islice)); }

    // new hists for KS test
    // TH1D *hFitKS = new TH1D(); h1D->Copy(*hFitKS);
    TH1D *hFitKS    = (TH1D*)h1D->Clone();
    TH1D *hTotalFit = FuncToHist(total,hFitKS,islice,arrMean,arrSigma);
    hTotalFit->SetName(Form("FitHist_%d",islice)); hFitKS ->SetName(Form("HistKS_%d",islice));
    hTotalFit->SetMarkerColor(kRed);               hFitKS ->SetMarkerColor(kBlue);
    hTotalFit->SetLineColor(kRed);                 hFitKS ->SetLineColor(kBlue);
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
    ResetParametersOfEachFunction(g1,g2,g3,g4,total,h1D,pt);
    PrepareLineShapes(hLineShape,total,g1,g2,g3,g4,centBin,etaBin,momBin,signBin,islice);
    //     PrepareLineShapes(h1D,fitTot,g1,g2,g3,g4,centBin,etaBin,momBin,signBin,islice);

    // Dump histograms into arrays
    h1D       -> SetName(Form("FinalFit_%d_%4.3f"  ,islice,pt));
    hResidual -> SetName(Form("Residual_%d_%4.3f"  ,islice,pt));

    if (dump) { debugFile->GetFile()->cd(); h1D -> Write(Form("FinalFit_%d_%4.3f",islice,pt)); }

    fitResArr[islice]    = (TH1D*)h1D;
    fitResiduals[islice] = (TH1D*)hResidual;

    // Dump some info into TTreeStream

    elFitAmp[islice]   = g1->GetParameter(0);
    piFitAmp[islice]   = g2->GetParameter(0);
    kaFitAmp[islice]   = g3->GetParameter(0);
    prFitAmp[islice]   = g4->GetParameter(0);

    elFitMean[islice]  = g1->GetParameter(1);
    piFitMean[islice]  = g2->GetParameter(1);
    kaFitMean[islice]  = g3->GetParameter(1);
    prFitMean[islice]  = g4->GetParameter(1);

    elFitSigma[islice] = g1->GetParameter(2);
    piFitSigma[islice] = g2->GetParameter(2);
    kaFitSigma[islice] = g3->GetParameter(2);
    prFitSigma[islice] = g4->GetParameter(2);

    elFitKurtosis[islice]  = g1->GetParameter(3);
    piFitKurtosis[islice]  = g2->GetParameter(3);
    kaFitKurtosis[islice]  = g3->GetParameter(3);
    prFitKurtosis[islice]  = g4->GetParameter(3);

    elFitSkew[islice]  = g1->GetParameter(4);
    piFitSkew[islice]  = g2->GetParameter(4);
    kaFitSkew[islice]  = g3->GetParameter(4);
    prFitSkew[islice]  = g4->GetParameter(4);

    elChi2Fit[islice]  = g1->GetChisquare();
    piChi2Fit[islice]  = g2->GetChisquare();
    kaChi2Fit[islice]  = g3->GetChisquare();
    prChi2Fit[islice]  = g4->GetChisquare();

    elIntegralFit[islice] = ComputeGaussIntegral(elFitAmp[islice],elFitMean[islice],elFitSigma[islice],g1->GetParameter(3),g1->GetParameter(4));
    piIntegralFit[islice] = ComputeGaussIntegral(piFitAmp[islice],piFitMean[islice],piFitSigma[islice],g2->GetParameter(3),g2->GetParameter(4));
    kaIntegralFit[islice] = ComputeGaussIntegral(kaFitAmp[islice],kaFitMean[islice],kaFitSigma[islice],g3->GetParameter(3),g3->GetParameter(4));
    prIntegralFit[islice] = ComputeGaussIntegral(prFitAmp[islice],prFitMean[islice],prFitSigma[islice],g4->GetParameter(3),g4->GetParameter(4));

    Double_t totalChi2    = total->GetChisquare();
    Double_t normNDF      = total->GetNDF();
    Double_t totChi2      = (normNDF<1) ? 1 : totalChi2/normNDF;
    Double_t pValue       = total->GetProb();
    Double_t ksValue      = hFitKS->KolmogorovTest(hTotalFit);
//     Double_t chi2Value    = hFitKS->Chi2Test(hTotalFit);

    // Get the fit param errors
    for (Int_t i=0;i<total->GetNpar();i++) parErrors[i] = total->GetParError(i);


    if (MCclosure && iIter>=4) {

      elFitMean[islice]  = cleanParams[0];
      piFitMean[islice]  = cleanParams[1];
      kaFitMean[islice]  = cleanParams[2];
      prFitMean[islice]  = cleanParams[3];

      elFitSigma[islice] = cleanParams[4];
      piFitSigma[islice] = cleanParams[5];
      kaFitSigma[islice] = cleanParams[6];
      prFitSigma[islice] = cleanParams[7];

      elFitSkew[islice]  = cleanParams[8];
      piFitSkew[islice]  = cleanParams[9];
      kaFitSkew[islice]  = cleanParams[10];
      prFitSkew[islice]  = cleanParams[11];

      elFitKurtosis[islice]  = cleanParams[12];
      piFitKurtosis[islice]  = cleanParams[13];
      kaFitKurtosis[islice]  = cleanParams[14];
      prFitKurtosis[islice]  = cleanParams[15];

      elFitAmp[islice]   = cleanParams[16];
      piFitAmp[islice]   = cleanParams[17];
      kaFitAmp[islice]   = cleanParams[18];
      prFitAmp[islice]   = cleanParams[19];

    }


    cout << " PiSkew = "   << skewnessFix[1]         ;
    cout << " PiKurt = "   << kurtosisFix[1]  << endl;
    cout << " KaSkew = "   << skewnessFix[2]         ;
    cout << " KaKurt = "   << kurtosisFix[2]  << endl;
    cout << " PrSkew = "   << skewnessFix[3]         ;
    cout << " PrKurt = "   << kurtosisFix[3]  << endl;
    cout << " ============================================  iter = " << iIter                     ;
    cout << " ptDownBin  = "           << ptDownBin                 ;
    cout << " ptUpBin = "              << ptUpBin                   ;
    cout << " islice = "               << islice                    ;
    cout << " p = "                    << pt                 << endl;
    cout << " --------------------------------------------------------------------------------------> kstest = " << ksValue;
    cout << " chi2 = "          << totChi2            << endl;


    // Fit results dump ttrees
    fitFile->GetFile()->cd();
    *fitFile << "IdMethodInput" <<

        "fSign="     << fSign                  <<
        "iter="      << iIter                 <<
        "slice="     << islice                <<
        "p="         << pt                    <<
        "eta="       << fEtaDown                <<
        "cent="      << fCentDown               <<

        "totalChi2=" << totalChi2             <<
        "totChi2="   << totChi2               <<
        "ksTest="    << ksValue               <<
        "pValue="    << pValue                <<
        "normNDF="   << normNDF               <<

        "piSkewScan="<< piSkewnessScan        <<
        "kaSkewScan="<< kaSkewnessScan        <<
        "prSkewScan="<< prSkewnessScan        <<
        "piKurtScan="<< piKurtosisScan        <<
        "kaKurtScan="<< kaKurtosisScan        <<
        "prKurtScan="<< prKurtosisScan        ;

	*fitFile << "IdMethodInput" <<

        "elAmp="     << elFitAmp[islice]      <<
        "piAmp="     << piFitAmp[islice]      <<
        "kaAmp="     << kaFitAmp[islice]      <<
        "prAmp="     << prFitAmp[islice]      <<

        "elMean="    << elFitMean[islice]     <<
        "piMean="    << piFitMean[islice]     <<
        "kaMean="    << kaFitMean[islice]     <<
        "prMean="    << prFitMean[islice]     <<

        "elSigma="   << elFitSigma[islice]    <<
        "piSigma="   << piFitSigma[islice]    <<
        "kaSigma="   << kaFitSigma[islice]    <<
        "prSigma="   << prFitSigma[islice]    ;

	*fitFile << "IdMethodInput" <<

        "elKurtosis="    << elFitKurtosis[islice]     <<
        "piKurtosis="    << piFitKurtosis[islice]     <<
        "kaKurtosis="    << kaFitKurtosis[islice]     <<
        "prKurtosis="    << prFitKurtosis[islice]     <<

        "elSkew="        << elFitSkew[islice]     <<
        "piSkew="        << piFitSkew[islice]     <<
        "kaSkew="        << kaFitSkew[islice]     <<
        "prSkew="        << prFitSkew[islice]     <<

        "elInt="         << elIntegralFit[islice] <<
        "piInt="         << piIntegralFit[islice] <<
        "kaInt="         << kaIntegralFit[islice] <<
        "prInt="         << prIntegralFit[islice] <<

        "elAmpErr="      << parErrors[0]      <<
        "piAmpErr="      << parErrors[5]      <<
        "kaAmpErr="      << parErrors[10]      <<
        "prAmpErr="      << parErrors[15]      ;

	*fitFile << "IdMethodInput" <<

        "elMeanErr="     << parErrors[1]      <<
        "piMeanErr="     << parErrors[6]      <<
        "kaMeanErr="     << parErrors[11]      <<
        "prMeanErr="     << parErrors[16]      <<

        "elSigmaErr="    << parErrors[2]      <<
        "piSigmaErr="    << parErrors[7]      <<
        "kaSigmaErr="    << parErrors[12]      <<
        "prSigmaErr="    << parErrors[17]      <<

        "elKurtosisErr=" << parErrors[3]      <<
        "piKurtosisErr=" << parErrors[8]      <<
        "kaKurtosisErr=" << parErrors[13]      <<
        "prKurtosisErr=" << parErrors[18]      <<

        "elSkewErr="     << parErrors[4]      <<
        "piSkewErr="     << parErrors[9]      <<
        "kaSkewErr="     << parErrors[14]      <<
        "prSkewErr="     << parErrors[19]      <<

        "\n";

    // Fit results dump ttrees
    *fitFile << "FitResults" <<

        "elMean="       << arrMean[0]         <<
        "piMean="       << arrMean[1]         <<
        "kaMean="       << arrMean[2]         <<
        "prMean="       << arrMean[3]         <<

        "elSigma="      << arrSigma[0]        <<
        "piSigma="      << arrSigma[1]        <<
        "kaSigma="      << arrSigma[2]        <<
        "prSigma="      << arrSigma[3]        <<

        "elWindow="     << arrMeanWindow[0]   <<
        "piWindow="     << arrMeanWindow[1]   <<
        "kaWindow="     << arrMeanWindow[2]   <<
        "prWindow="     << arrMeanWindow[3]   ;

    *fitFile << "FitResults" <<

        "fSign="         << fSign               <<
        "iter="         << iIter              <<
        "slice="        << islice             <<
        "p="            << pt                 <<
        "eta="          << fEtaDown             <<
        "cent="         << fCentDown            <<
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

        "elFitMean="    << elFitMean[islice]  <<
        "piFitMean="    << piFitMean[islice]  <<
        "kaFitMean="    << kaFitMean[islice]  <<
        "prFitMean="    << prFitMean[islice]  ;

    *fitFile << "FitResults" <<

        "elFitSigma="   << elFitSigma[islice] <<
        "piFitSigma="   << piFitSigma[islice] <<
        "kaFitSigma="   << kaFitSigma[islice] <<
        "prFitSigma="   << prFitSigma[islice] <<

        "elFitAmp="     << elFitAmp[islice]   <<
        "piFitAmp="     << piFitAmp[islice]   <<
        "kaFitAmp="     << kaFitAmp[islice]   <<
        "prFitAmp="     << prFitAmp[islice]   <<

        "elKurtosis="    << elFitKurtosis[islice]     <<
        "piKurtosis="    << piFitKurtosis[islice]     <<
        "kaKurtosis="    << kaFitKurtosis[islice]     <<
        "prKurtosis="    << prFitKurtosis[islice]     ;

	*fitFile << "FitResults" <<

        "elFitSkew="    << elFitSkew[islice]  <<
        "piFitSkew="    << piFitSkew[islice]  <<
        "kaFitSkew="    << kaFitSkew[islice]  <<
        "prFitSkew="    << prFitSkew[islice]  <<

        "elChi2Fit="    << elChi2Fit[islice]  <<
        "piChi2Fit="    << piChi2Fit[islice]  <<
        "kaChi2Fit="    << kaChi2Fit[islice]  <<
        "prChi2Fit="    << prChi2Fit[islice]  <<

        "elAmpErr="      << parErrors[0]      <<
        "piAmpErr="      << parErrors[5]      <<
        "kaAmpErr="      << parErrors[10]      <<
        "prAmpErr="      << parErrors[15]      ;

	*fitFile << "FitResults" <<

        "elMeanErr="     << parErrors[1]      <<
        "piMeanErr="     << parErrors[6]      <<
        "kaMeanErr="     << parErrors[11]      <<
        "prMeanErr="     << parErrors[16]      <<

        "elSigmaErr="    << parErrors[2]      <<
        "piSigmaErr="    << parErrors[7]      <<
        "kaSigmaErr="    << parErrors[12]      <<
        "prSigmaErr="    << parErrors[17]      <<

        "elKurtosisErr=" << parErrors[3]      <<
        "piKurtosisErr=" << parErrors[8]      <<
        "kaKurtosisErr=" << parErrors[13]      <<
        "prKurtosisErr=" << parErrors[18]      <<

        "elSkewErr="     << parErrors[4]      <<
        "piSkewErr="     << parErrors[9]      <<
        "kaSkewErr="     << parErrors[14]      <<
        "prSkewErr="     << parErrors[19]      <<

        "\n";

    // write hTotalFit for debugging
    histsFile -> GetFile()->cd();
//     if (islice>10 && islice<20) {
    if(nSlice<=10){
      hTotalFit->Write(Form("FitHist_%d",islice));
      hFitKS   ->Write(Form("HistKS_%d",islice));
    }

    delete [] arrSigma;
    delete [] arrMean;
    delete [] arrMeanWindow;
    delete [] cleanParams;
    delete hLineShape;
    delete resVector;
    delete perVector;
    delete dEdxVector;
    delete hTotalFit;
    delete hFitKS;
    delete hClone;

  }


  histsFile -> GetFile()->cd();
  if (iIter ==0 || iIter >=3) fitResArr . Write("fitResArr",TObject::kSingleKey);
  if (iIter >=3)  fitResiduals . Write("fitResiduals",TObject::kSingleKey);

  lineShapesLookUp->GetFile()->cd();
  histLineShapesCArr . Write("histLineShapesCArr",TObject::kSingleKey);
  histLineShapesCArr.Clear("C");
  funcLineShapesCArr . Write("funcLineShapesCArr",TObject::kSingleKey);
  funcLineShapesCArr.Clear("C");
  funcLineShapesObjArr . Write("funcLineShapesObjArr",TObject::kSingleKey);
  funcLineShapesObjArr.Delete();



  // delete arrays
  fitResArr    . Clear("C");
  fitResiduals . Clear("C");



  delete histsFile;
  delete fitFile;
  delete lineShapesLookUp;

  if (ressFile!=0x0) ressFile->Close();
  if (sampFile!=0x0) sampFile->Close();

}
// -------------------------------------------------------------------------------------------------------
TGraphErrors *RemoveOutliers(TGraphErrors *grsmooth, Int_t fitWindow)
{


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
      Double_t xx[4]={Double_t(idelta),Double_t(idelta),Double_t(idelta),Double_t(idelta)};
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
void CheckIfMeanIsZero(TF1 *total, Double_t *arrMean)
{

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

}
// -------------------------------------------------------------------------------------------------------
void SetParNamesOfTotalFit(TF1 *total)
{

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

}
// --------------------------------------------------------------------------------------------
void SetFixedKSparameterForTotalFit(TF1 *total)
{

  // Set kurtosis
  if (fixedK){
    total->FixParameter(3, kurtosisFix[kElectron]);
    total->FixParameter(8, kurtosisFix[kPion]    );
    total->FixParameter(13,kurtosisFix[kKaon]    );
    total->FixParameter(18,kurtosisFix[kProton]  );
  }

  // Set skewness
  if (fixedS) {
    total->FixParameter(4, skewnessFix[kElectron]);
    total->FixParameter(9, skewnessFix[kPion]    );
    total->FixParameter(14,skewnessFix[kKaon]    );
    total->FixParameter(19,skewnessFix[kProton]  );
  }


}
// --------------------------------------------------------------------------------------------
void SetFixedKSparameterForTotalFitForPions(TF1 *total)
{

  // Set kurtosis
  total->FixParameter(3, kurtosisFix[kElectron]);

  total->SetParLimits(8,kurtosisMin,kurtosisMax);

  total->FixParameter(13,kurtosisFix[kKaon]    );
  total->FixParameter(18,kurtosisFix[kProton]  );

  // Set skewness
  total->FixParameter(4, skewnessFix[kElectron]);

  total->SetParLimits(9,skewMin,skewMax);

  total->FixParameter(14,skewnessFix[kKaon]    );
  total->FixParameter(19,skewnessFix[kProton]  );

}
// --------------------------------------------------------------------------------------------
void SetFixedKSparameterForIndividualFits(TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4)
{

  // Set kurtosis
  if (fixedK) {
    g1->FixParameter(3,kurtosisFix[kElectron]);
    g2->FixParameter(3,kurtosisFix[kPion]    );
    g3->FixParameter(3,kurtosisFix[kKaon]    );
    g4->FixParameter(3,kurtosisFix[kProton]  );
  }


  // Set Skewness
  if (fixedS) {
    g1->FixParameter(4,skewnessFix[kElectron]);
    g2->FixParameter(4,skewnessFix[kPion]    );
    g3->FixParameter(4,skewnessFix[kKaon]    );
    g4->FixParameter(4,skewnessFix[kProton]  );
  }

}
// -------------------------------------------------------------------------------------------------------
void SetFixedKSparameterForCleanSampleFit(TF1 *asyGaus, ParticleType pSpecy)
{


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
  }

}
// -------------------------------------------------------------------------------------------------------
void GetHelperGraphs(Int_t iIter, TFile *ressFile)
{

  //
  // Read helper graphs from the smooth tree
  //

  if (iIter>0) {
//     get the file and objarrays of previous fit results
    Fit           = (TObjArray*)ressFile->Get("Fit");
    Clean         = (TObjArray*)ressFile->Get("Clean");
    ScaledClean   = (TObjArray*)ressFile->Get("ScaledClean");
    ScaledEx      = (TObjArray*)ressFile->Get("ScaledEx");
    OutlierSmooth = (TObjArray*)ressFile->Get("OulierSmooth");
    Windows       = (TObjArray*)ressFile->Get("Windows");

    // previous fir params
    grelFitAmp    = (TGraphErrors*)Fit ->FindObject("elFitAmp");
    grpiFitAmp    = (TGraphErrors*)Fit ->FindObject("piFitAmp");
    grkaFitAmp    = (TGraphErrors*)Fit ->FindObject("kaFitAmp");
    grprFitAmp    = (TGraphErrors*)Fit ->FindObject("prFitAmp");

    grelFitMean   = (TGraphErrors*)Fit ->FindObject("elFitMean");
    grpiFitMean   = (TGraphErrors*)Fit ->FindObject("piFitMean");
    grkaFitMean   = (TGraphErrors*)Fit ->FindObject("kaFitMean");
    grprFitMean   = (TGraphErrors*)Fit ->FindObject("prFitMean");

    grelFitSigma  = (TGraphErrors*)Fit ->FindObject("elFitSigma");
    grpiFitSigma  = (TGraphErrors*)Fit ->FindObject("piFitSigma");
    grkaFitSigma  = (TGraphErrors*)Fit ->FindObject("kaFitSigma");
    grprFitSigma  = (TGraphErrors*)Fit ->FindObject("prFitSigma");

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
    grprCleanMeanScaledSmooth  = (TGraphErrors*)ScaledClean  ->FindObject("prCleanMeanScaledSmooth");
    grprFitAmpOutlierSmooth    = (TGraphErrors*)OutlierSmooth->FindObject("prFitAmpOutlierSmooth");

//     Graphsto be used MC analysis case
    grelCleanSigmaScaled = (TGraphErrors*)ScaledClean  ->FindObject("elCleanSigmaScaled");
    grelCleanMeanScaled  = (TGraphErrors*)ScaledClean  ->FindObject("elCleanMeanScaled");
    grelCleanAmpScaled   = (TGraphErrors*)ScaledClean  ->FindObject("elCleanAmpScaled");
    grpiCleanMeanScaled  = (TGraphErrors*)ScaledClean  ->FindObject("piCleanMeanScaled");
    grpiCleanSigmaScaled = (TGraphErrors*)ScaledClean  ->FindObject("piCleanSigmaScaled");
    grpiCleanAmpScaled   = (TGraphErrors*)ScaledClean  ->FindObject("piCleanAmpScaled");
    grkaCleanAmpScaled   = (TGraphErrors*)ScaledClean  ->FindObject("kaCleanAmpScaled");
    grprCleanAmpScaled   = (TGraphErrors*)ScaledClean  ->FindObject("prCleanAmpScaled");

    grkaCleanMeanScaled  = (TGraphErrors*)ScaledClean  ->FindObject("kaCleanMeanScaled");
    grprCleanSigmaScaled = (TGraphErrors*)ScaledClean  ->FindObject("prCleanSigmaScaled");
    grprCleanMeanScaled  = (TGraphErrors*)ScaledClean  ->FindObject("prCleanMeanScaled");

    //     Graphsto be used MC closure test
    grelFreeKSAmp      = (TGraphErrors*)Clean  ->FindObject("elCleanAmpFreeKS");
    grpiFreeKSAmp      = (TGraphErrors*)Clean  ->FindObject("piCleanAmpFreeKS");
    grkaFreeKSAmp      = (TGraphErrors*)Clean  ->FindObject("kaCleanAmpFreeKS");
    grprFreeKSAmp      = (TGraphErrors*)Clean  ->FindObject("prCleanAmpFreeKS");

    grelFreeKSMean     = (TGraphErrors*)Clean  ->FindObject("elCleanMeanFreeKS");
    grpiFreeKSMean     = (TGraphErrors*)Clean  ->FindObject("piCleanMeanFreeKS");
    grkaFreeKSMean     = (TGraphErrors*)Clean  ->FindObject("kaCleanMeanFreeKS");
    grprFreeKSMean     = (TGraphErrors*)Clean  ->FindObject("prCleanMeanFreeKS");

    grelFreeKSSigma    = (TGraphErrors*)Clean  ->FindObject("elCleanSigmaFreeKS");
    grpiFreeKSSigma    = (TGraphErrors*)Clean  ->FindObject("piCleanSigmaFreeKS");
    grkaFreeKSSigma    = (TGraphErrors*)Clean  ->FindObject("kaCleanSigmaFreeKS");
    grprFreeKSSigma    = (TGraphErrors*)Clean  ->FindObject("prCleanSigmaFreeKS");

    grelFreeKSKurtosis = (TGraphErrors*)Clean  ->FindObject("elCleanKurtosisFreeKS");
    grpiFreeKSKurtosis = (TGraphErrors*)Clean  ->FindObject("piCleanKurtosisFreeKS");
    grkaFreeKSKurtosis = (TGraphErrors*)Clean  ->FindObject("kaCleanKurtosisFreeKS");
    grprFreeKSKurtosis = (TGraphErrors*)Clean  ->FindObject("prCleanKurtosisFreeKS");

    grelFreeKSSkew     = (TGraphErrors*)Clean  ->FindObject("elCleanSkewFreeKS");
    grpiFreeKSSkew     = (TGraphErrors*)Clean  ->FindObject("piCleanSkewFreeKS");
    grkaFreeKSSkew     = (TGraphErrors*)Clean  ->FindObject("kaCleanSkewFreeKS");
    grprFreeKSSkew     = (TGraphErrors*)Clean  ->FindObject("prCleanSkewFreeKS");

    grelCleanAmp       = (TGraphErrors*)Clean  ->FindObject("elCleanAmp");
    grpiCleanAmp       = (TGraphErrors*)Clean  ->FindObject("piCleanAmp");
    grkaCleanAmp       = (TGraphErrors*)Clean  ->FindObject("kaCleanAmp");
    grprCleanAmp       = (TGraphErrors*)Clean  ->FindObject("prCleanAmp");

  }

}
// -------------------------------------------------------------------------------------------------------
void SetTotalFitParameters(Int_t iIter, Int_t islice, Int_t nSlice, TF1* total, Double_t *arrMean, Double_t *arrSigma, Double_t *arrMeanWindow, Double_t *cleanParams, TH1D* h1D, Double_t pt)
{

  //
  // Set total Fit Parameters for each iteration
  //

  if (islice<nSlice) {

    // Set the mean position wrt closest particle
    arrMeanWindow[0] = GetClosestParticleMean(arrSigma[0],arrMean[0],arrMean[1],arrMean[2],arrMean[3]);
    arrMeanWindow[1] = GetClosestParticleMean(arrSigma[1],arrMean[1],arrMean[0],arrMean[2],arrMean[3]);
    arrMeanWindow[2] = GetClosestParticleMean(arrSigma[2],arrMean[2],arrMean[1],arrMean[0],arrMean[3]);
    arrMeanWindow[3] = GetClosestParticleMean(arrSigma[3],arrMean[3],arrMean[1],arrMean[2],arrMean[0]);

    if (iIter>=1) { // Fix only mean positions
      Double_t elAssump=0., piAssump=0., kaAssump=0., prAssump=0.;
      if (assump==0){  // use only PID response
        elAssump = grelMeanExScaled  -> GetY()[islice];
        piAssump = ((pt>=0.6 && pt<=2.)) ? grpiCleanMeanScaledSmooth -> GetY()[islice]: grpiFitMeanSmooth -> GetY()[islice];
        kaAssump = grkaMeanExScaled  -> GetY()[islice];
        prAssump = grprMeanExScaled  -> GetY()[islice];
      } else if (assump==1){  // use clean samples
        elAssump = grelMeanExScaled  -> GetY()[islice];
        piAssump = ((pt>=0.6 && pt<=2.)) ? grpiCleanMeanScaledSmooth -> GetY()[islice]: grpiFitMeanSmooth -> GetY()[islice];
        kaAssump = grkaMeanExScaled  -> GetY()[islice];
        prAssump = grprCleanMeanScaledSmooth->GetY()[islice];
      }
    // Recalculate the windows with the clean samples
      arrMeanWindow[0] = GetClosestParticleMean(arrSigma[0], elAssump, piAssump, kaAssump, prAssump);
      arrMeanWindow[1] = GetClosestParticleMean(arrSigma[1], piAssump, elAssump, kaAssump, prAssump);
      arrMeanWindow[2] = GetClosestParticleMean(arrSigma[2], kaAssump, piAssump, elAssump, prAssump);
      arrMeanWindow[3] = GetClosestParticleMean(arrSigma[3], prAssump, piAssump, kaAssump, elAssump);
    }

     //   set initials params using PIDresponse
    if (iIter>=0) SetParams1stIteration(pt,h1D,arrMean,arrSigma,arrMeanWindow,cleanParams);

    // Fix only mean positions
    if (iIter>=1) {

      // ***************  Electron ***************                 good
       elMeanMin  = grelMeanExScaled->GetY()[islice] - arrMeanWindow[0];
       elMeanMax  = grelMeanExScaled->GetY()[islice] + arrMeanWindow[0];
//      elMeanMin  = grelMeanExScaled->GetY()[islice]*0.99;
//      elMeanMax  = grelMeanExScaled->GetY()[islice]*0.01;
      // ***************  Pion  ***************
      if (pt>=0.6 && pt<=2.){   // in this interval clean samples are ok
        piMeanMin  = grpiCleanMeanScaledSmooth->GetY()[islice]*0.995;
        piMeanMax  = grpiCleanMeanScaledSmooth->GetY()[islice]*1.005;
      }
      // ***************  kaon  ***************
      if (pt>=0.7){
        kaMeanMin  = grkaMeanExScaled->GetY()[islice] - arrMeanWindow[2];
        kaMeanMax  = grkaMeanExScaled->GetY()[islice] + arrMeanWindow[2];
      }
      // ***************  Proton  ***************
      if (pt>=0.75){
        prMeanMin  = grprCleanMeanScaledSmooth->GetY()[islice]*0.997;
        prMeanMax  = grprCleanMeanScaledSmooth->GetY()[islice]*1.003;
      }
      //============================= SOME SAFE PATCHES ===============================================
      // electron amplitude tends to rise at el-K and el-pr overlap --> restrict elAmpMax
      if ((pt>=0.45 && pt<=0.65) || (pt>=0.8 && pt<=1.3) ){
        elAmpMax  = grelFitAmpSmooth->GetY()[islice];
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=2) {

      // %%%%%%%%%%%%%%%  el-K and  %%%%%%%%%%%%%%%%
      if ((pt>=0.45 && pt<=0.65)){
        elAmpMax   = grelFitAmpSmooth->GetY()[islice];
        kaMeanMin  = grkaMeanExScaled ->GetY()[islice] - grkaMeanExScaled ->GetY()[islice]*0.005;
        kaMeanMax  = grkaMeanExScaled ->GetY()[islice] + grkaMeanExScaled ->GetY()[islice]*0.005;
      }
      // %%%%%%%%%%%%%%%  el-K and  %%%%%%%%%%%%%%%%
      //
      // el-pr overlap
      if ((pt>=0.8 && pt<=1.3)){
        elAmpMax   = grelFitAmpSmooth->GetY()[islice];
        prMeanMin  = grprMeanExScaled ->GetY()[islice] - grprMeanExScaled ->GetY()[islice]*0.001;
        prMeanMax  = grprMeanExScaled ->GetY()[islice] + grprMeanExScaled ->GetY()[islice]*0.001;
      }

      if (pt>=0.75){
        prSigmaMin = grprCleanSigmaScaledSmooth->GetY()[islice]*0.95;
        prSigmaMax = grprCleanSigmaScaledSmooth->GetY()[islice]*1.05;
      }

      if (pt>=1.6 && pt<=1.9){
        prAmpMin   = grprFitAmpOutlierSmooth->GetY()[islice]*0.9;
        prAmpMax   = grprFitAmpOutlierSmooth->GetY()[islice]*1.05;
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=3) {

      // el amp fix
      if ((pt>=0.45 && pt<=0.65) || (pt>=0.8 && pt<=1.3) ){
	elAmpMin  = grelFitAmpSmooth->GetY()[islice] - grelFitAmpSmooth->GetY()[islice]*0.005;
	elAmpMax  = grelFitAmpSmooth->GetY()[islice] + grelFitAmpSmooth->GetY()[islice]*0.005;
      }
      //  kaon sigma issue
      //if (pt>=0.5 && pt<=1.){
      //  kaSigmaMin  = grkaCleanSigma->GetY()[islice]*0.95;
      //  kaSigmaMax  = grkaCleanSigma->GetY()[islice]*0.05;
      //}

      //  pion fix
      if (pt>=0.65) {
	piAmpMin   = grpiFitAmpSmooth   ->GetY()[islice] - grpiFitAmpSmooth   ->GetY()[islice]*0.05;
	piAmpMax   = grpiFitAmpSmooth   ->GetY()[islice] + grpiFitAmpSmooth   ->GetY()[islice]*0.05;
	piMeanMin  = grpiFitMeanSmooth  ->GetY()[islice] - grpiFitMeanSmooth  ->GetY()[islice]*0.002;
	piMeanMax  = grpiFitMeanSmooth  ->GetY()[islice] + grpiFitMeanSmooth  ->GetY()[islice]*0.002;
	piSigmaMin = grpiFitSigmaSmooth ->GetY()[islice] - grpiFitSigmaSmooth ->GetY()[islice]*0.05;
	piSigmaMax = grpiFitSigmaSmooth ->GetY()[islice] + grpiFitSigmaSmooth ->GetY()[islice]*0.05;
      }

      // proton
      if (pt>=1.6 && pt<=1.9){
	prAmpMin   = grprFitAmpOutlierSmooth->GetY()[islice]*0.9;
	prAmpMax   = grprFitAmpOutlierSmooth->GetY()[islice];
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=4) {

       //     Electron
      if (pt>=0.4){
	elSigmaMin = grelFitSigmaSmooth ->GetY()[islice] - grelFitSigmaSmooth ->GetY()[islice]*0.05;
        elSigmaMax = grelFitSigmaSmooth ->GetY()[islice] + grelFitSigmaSmooth ->GetY()[islice]*0.05;
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
        prAmpMax   = grprFitAmpSmooth   ->GetY()[islice] + grprFitAmpSmooth   ->GetY()[islice]*0.01;
        prMeanMin  = grprFitMeanSmooth  ->GetY()[islice] - grprFitMeanSmooth  ->GetY()[islice]*0.01;
        prMeanMax  = grprFitMeanSmooth  ->GetY()[islice] + grprFitMeanSmooth  ->GetY()[islice]*0.01;
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=5) {

      //     Electron
      if (pt>=0.4){
        elAmpMin   = grelFitAmpSmooth   ->GetY()[islice] - grelFitAmpSmooth   ->GetY()[islice]*0.005;
	elAmpMax   = grelFitAmpSmooth   ->GetY()[islice] + grelFitAmpSmooth   ->GetY()[islice]*0.005;
	elSigmaMin = grelFitSigmaSmooth ->GetY()[islice] - grelFitSigmaSmooth ->GetY()[islice]*0.05;
        elSigmaMax = grelFitSigmaSmooth ->GetY()[islice] + grelFitSigmaSmooth ->GetY()[islice]*0.05;
      }

      //     Pion
      if (pt>0.65) {
	piAmpMin   = grpiFitAmpSmooth   ->GetY()[islice] - grpiFitAmpSmooth   ->GetY()[islice]*0.005;
	piAmpMax   = grpiFitAmpSmooth   ->GetY()[islice] + grpiFitAmpSmooth   ->GetY()[islice]*0.005;
        piMeanMin  = grpiFitMeanSmooth  ->GetY()[islice] - grpiFitMeanSmooth  ->GetY()[islice]*0.005;
        piMeanMax  = grpiFitMeanSmooth  ->GetY()[islice] + grpiFitMeanSmooth  ->GetY()[islice]*0.005;
      }

    //     Kaon
      if (pt>=0.65) {
        kaAmpMin   = grkaFitAmpSmooth   ->GetY()[islice] - grkaFitAmpSmooth   ->GetY()[islice]*0.05;
        kaAmpMax   = grkaFitAmpSmooth   ->GetY()[islice] + grkaFitAmpSmooth   ->GetY()[islice]*0.05;
       }

//     Proton
      if (pt>1.3) {
        prMeanMin  = grprFitMeanSmooth  ->GetY()[islice] - grprFitMeanSmooth  ->GetY()[islice]*0.005;
        prMeanMax  = grprFitMeanSmooth  ->GetY()[islice] + grprFitMeanSmooth  ->GetY()[islice]*0.005;
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=6) {

//     Electron
      if (pt>=0.4){
        elAmpMin   = grelFitAmpSmooth   ->GetY()[islice] - grelFitAmpSmooth   ->GetY()[islice]*0.005;
        elAmpMax   = grelFitAmpSmooth   ->GetY()[islice] + grelFitAmpSmooth   ->GetY()[islice]*0.005;
        elSigmaMin = grelFitSigmaSmooth ->GetY()[islice] - grelFitSigmaSmooth ->GetY()[islice]*0.05;
        elSigmaMax = grelFitSigmaSmooth ->GetY()[islice] + grelFitSigmaSmooth ->GetY()[islice]*0.05;
      }

//     Pion
      if (pt>0.75) {
        piAmpMin   = grpiFitAmpSmooth   ->GetY()[islice] - grpiFitAmpSmooth   ->GetY()[islice]*sysAmp;
        piAmpMax   = grpiFitAmpSmooth   ->GetY()[islice] + grpiFitAmpSmooth   ->GetY()[islice]*sysAmp;
        piMeanMin  = grpiFitMeanSmooth  ->GetY()[islice] - grpiFitMeanSmooth  ->GetY()[islice]*sysMean;
        piMeanMax  = grpiFitMeanSmooth  ->GetY()[islice] + grpiFitMeanSmooth  ->GetY()[islice]*sysMean;
      }

//     Kaon
      if (pt>=0.8) {
        kaAmpMin   = grkaFitAmpSmooth   ->GetY()[islice] - grkaFitAmpSmooth   ->GetY()[islice]*sysAmp;
        kaAmpMax   = grkaFitAmpSmooth   ->GetY()[islice] + grkaFitAmpSmooth   ->GetY()[islice]*sysAmp;
        kaMeanMin  = grkaFitMeanSmooth  ->GetY()[islice] - grkaFitMeanSmooth  ->GetY()[islice]*sysMean;
        kaMeanMax  = grkaFitMeanSmooth  ->GetY()[islice] + grkaFitMeanSmooth  ->GetY()[islice]*sysMean;
        kaSigmaMin = grkaFitSigmaSmooth ->GetY()[islice] - grkaFitSigmaSmooth ->GetY()[islice]*0.05;
        kaSigmaMax = grkaFitSigmaSmooth ->GetY()[islice] + grkaFitSigmaSmooth ->GetY()[islice]*0.05;
      }

//     Proton
      if (pt>1.) {
        prAmpMin   = grprFitAmpSmooth   ->GetY()[islice] - grprFitAmpSmooth   ->GetY()[islice]*sysAmp;
        prAmpMax   = grprFitAmpSmooth   ->GetY()[islice] + grprFitAmpSmooth   ->GetY()[islice]*sysAmp;
        prMeanMin  = grprFitMeanSmooth  ->GetY()[islice] - grprFitMeanSmooth  ->GetY()[islice]*sysMean;
        prMeanMax  = grprFitMeanSmooth  ->GetY()[islice] + grprFitMeanSmooth  ->GetY()[islice]*sysMean;
        prSigmaMin = grprFitSigmaSmooth ->GetY()[islice] - grprFitSigmaSmooth ->GetY()[islice]*0.05;
        prSigmaMax = grprFitSigmaSmooth ->GetY()[islice] + grprFitSigmaSmooth ->GetY()[islice]*0.05;
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  }  // close settings for 20MeV bin width
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //   Almost no restriction in low momentum
  if (pt<0.3){

    Double_t tmpmaxBin  = h1D->GetBinContent(h1D->GetMaximumBin());

    elAmpMin = 1e-4;
    piAmpMin = tmpmaxBin/1.5;
    kaAmpMin = 1e-4;
    prAmpMin = 1e-4;

    elAmpMax = tmpmaxBin/2.;
    piAmpMax = tmpmaxBin*1.1;
    kaAmpMax = tmpmaxBin;
    prAmpMax = tmpmaxBin;

    elMeanMin = TMath::Max(30.,(arrMean[0]-4.*arrSigma[0]));
    piMeanMin = TMath::Max(30.,(arrMean[1]-4.*arrSigma[1]));
    kaMeanMin = TMath::Max(30.,(arrMean[2]-4.*arrSigma[2]));

    elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+4.*arrSigma[0]));
    piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+4.*arrSigma[1]));
    kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+4.*arrSigma[2]));

    // pr and de is special for low momentum
    prMeanMin = TMath::Max(30.,(arrMean[3]-4.*arrSigma[3]));
    prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+4.*arrSigma[3]));

    elSigmaMin = arrSigma[0]/5.;
    piSigmaMin = arrSigma[1]/5.;
    kaSigmaMin = arrSigma[2]/5.;
    prSigmaMin = arrSigma[3]/5.;

    elSigmaMax = arrSigma[0]*3.;
    piSigmaMax = arrSigma[1]*3.;
    kaSigmaMax = arrSigma[2]*3.;
    prSigmaMax = arrSigma[3]*3.;

  }

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //   Let protons free below 0.7 gev
  if (pt<0.7 && !MCclosure){

    // ------- proton -------
    prAmpMin = 1e-4;
    prAmpMax = h1D->GetBinContent(h1D->GetMaximumBin());
    prMeanMin = TMath::Max(30.,(arrMean[3]-3.*arrSigma[3]));
    prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+3.*arrSigma[3]));
    prSigmaMin = arrSigma[3]/5.;
    prSigmaMax = arrSigma[3]*3.;

  }

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //   Finally set the params

  // set the Kurtosis and Skewness as free, fix them if needed (incase modify fixedK and fixedS)
  total->SetParLimits(4,skewMin,skewMax);    // electron
  total->SetParLimits(9,skewMin,skewMax);    // pion
  total->SetParLimits(14,skewMin,skewMax);   // kaon
  total->SetParLimits(19,skewMin,skewMax);   // proton

  total->SetParLimits(3,kurtosisMin,kurtosisMax);   // electron
  total->SetParLimits(8,kurtosisMin,kurtosisMax);   // pion
  total->SetParLimits(13,kurtosisMin,kurtosisMax);  // kaon
  total->SetParLimits(18,kurtosisMin,kurtosisMax);  // proton

  // if needed fix Kurtosis and Skewness
  if ( (fixedK || fixedS) ) SetFixedKSparameterForTotalFit(total);

   // ++++++++++++++++++++++++ more freedom for some particles +++++++++++++++++++

  // let pion free for p<0.6
  if ( iIter==7 ) {

    total->SetParLimits(9,skewMin,skewMax);     // let skewness free
    total->SetParLimits(8,kurtosisMin,kurtosisMax); // let the kurtosis free

    kaAmpMin   = grkaFitAmp   ->GetY()[islice] - grkaFitAmp   ->GetY()[islice]*0.001;
    kaAmpMax   = grkaFitAmp   ->GetY()[islice] + grkaFitAmp   ->GetY()[islice]*0.001;
    kaMeanMin  = grkaFitMean  ->GetY()[islice] - grkaFitMean  ->GetY()[islice]*0.001;
    kaMeanMax  = grkaFitMean  ->GetY()[islice] + grkaFitMean  ->GetY()[islice]*0.001;
    kaSigmaMin = grkaFitSigma ->GetY()[islice] - grkaFitSigma ->GetY()[islice]*0.001;
    kaSigmaMax = grkaFitSigma ->GetY()[islice] + grkaFitSigma ->GetY()[islice]*0.001;

    elAmpMin   = grelFitAmp   ->GetY()[islice] - grelFitAmp   ->GetY()[islice]*0.001;
    elAmpMax   = grelFitAmp   ->GetY()[islice] + grelFitAmp   ->GetY()[islice]*0.001;
    elMeanMin  = grelFitMean  ->GetY()[islice] - grelFitMean  ->GetY()[islice]*0.001;
    elMeanMax  = grelFitMean  ->GetY()[islice] + grelFitMean  ->GetY()[islice]*0.001;
    elSigmaMin = grelFitSigma ->GetY()[islice] - grelFitSigma ->GetY()[islice]*0.001;
    elSigmaMax = grelFitSigma ->GetY()[islice] + grelFitSigma ->GetY()[islice]*0.001;

    piAmpMin   = grpiFitAmp   ->GetY()[islice] - grpiFitAmp   ->GetY()[islice]*0.2;
    piAmpMax   = grpiFitAmp   ->GetY()[islice] + grpiFitAmp   ->GetY()[islice]*0.2;
    piMeanMin  = grpiFitMean  ->GetY()[islice] - grpiFitMean  ->GetY()[islice]*0.001;
    piMeanMax  = grpiFitMean  ->GetY()[islice] + grpiFitMean  ->GetY()[islice]*0.001;
    piSigmaMin = grpiFitSigma ->GetY()[islice] - grpiFitSigma ->GetY()[islice]*0.2;
    piSigmaMax = grpiFitSigma ->GetY()[islice] + grpiFitSigma ->GetY()[islice]*0.2;

    prAmpMin   = grprFitAmp   ->GetY()[islice] - grprFitAmp   ->GetY()[islice]*0.001;
    prAmpMax   = grprFitAmp   ->GetY()[islice] + grprFitAmp   ->GetY()[islice]*0.001;
    prMeanMin  = grprFitMean  ->GetY()[islice] - grprFitMean  ->GetY()[islice]*0.001;
    prMeanMax  = grprFitMean  ->GetY()[islice] + grprFitMean  ->GetY()[islice]*0.001;
    prSigmaMin = grprFitSigma ->GetY()[islice] - grprFitSigma ->GetY()[islice]*0.001;
    prSigmaMax = grprFitSigma ->GetY()[islice] + grprFitSigma ->GetY()[islice]*0.001;

  }

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // let kaon free after 7th iter
  if ( iIter==8 && pt>=0.5) {

    total->SetParLimits(14,skewnessFix[kKaon]*0.9,skewnessFix[kKaon]*1.1);
    total->SetParLimits(13,kurtosisFix[kKaon]*0.9,kurtosisFix[kKaon]*1.1);

    kaAmpMin   = grkaFitAmp   ->GetY()[islice] - grkaFitAmp   ->GetY()[islice]*0.1;
    kaAmpMax   = grkaFitAmp   ->GetY()[islice] + grkaFitAmp   ->GetY()[islice]*0.1;
    kaMeanMin  = grkaFitMean  ->GetY()[islice] - grkaFitMean  ->GetY()[islice]*0.001;
    kaMeanMax  = grkaFitMean  ->GetY()[islice] + grkaFitMean  ->GetY()[islice]*0.001;
    kaSigmaMin = grkaFitSigma ->GetY()[islice] - grkaFitSigma ->GetY()[islice]*0.1;
    kaSigmaMax = grkaFitSigma ->GetY()[islice] + grkaFitSigma ->GetY()[islice]*0.1;

    elAmpMin   = grelFitAmp   ->GetY()[islice] - grelFitAmp   ->GetY()[islice]*0.001;
    elAmpMax   = grelFitAmp   ->GetY()[islice] + grelFitAmp   ->GetY()[islice]*0.001;
    elMeanMin  = grelFitMean  ->GetY()[islice] - grelFitMean  ->GetY()[islice]*0.001;
    elMeanMax  = grelFitMean  ->GetY()[islice] + grelFitMean  ->GetY()[islice]*0.001;
    elSigmaMin = grelFitSigma ->GetY()[islice] - grelFitSigma ->GetY()[islice]*0.001;
    elSigmaMax = grelFitSigma ->GetY()[islice] + grelFitSigma ->GetY()[islice]*0.001;

    piAmpMin   = grpiFitAmp   ->GetY()[islice] - grpiFitAmp   ->GetY()[islice]*0.001;
    piAmpMax   = grpiFitAmp   ->GetY()[islice] + grpiFitAmp   ->GetY()[islice]*0.001;
    piMeanMin  = grpiFitMean  ->GetY()[islice] - grpiFitMean  ->GetY()[islice]*0.001;
    piMeanMax  = grpiFitMean  ->GetY()[islice] + grpiFitMean  ->GetY()[islice]*0.001;
    piSigmaMin = grpiFitSigma ->GetY()[islice] - grpiFitSigma ->GetY()[islice]*0.001;
    piSigmaMax = grpiFitSigma ->GetY()[islice] + grpiFitSigma ->GetY()[islice]*0.001;

    prAmpMin   = grprFitAmp   ->GetY()[islice] - grprFitAmp   ->GetY()[islice]*0.001;
    prAmpMax   = grprFitAmp   ->GetY()[islice] + grprFitAmp   ->GetY()[islice]*0.001;
    prMeanMin  = grprFitMean  ->GetY()[islice] - grprFitMean  ->GetY()[islice]*0.001;
    prMeanMax  = grprFitMean  ->GetY()[islice] + grprFitMean  ->GetY()[islice]*0.001;
    prSigmaMin = grprFitSigma ->GetY()[islice] - grprFitSigma ->GetY()[islice]*0.001;
    prSigmaMax = grprFitSigma ->GetY()[islice] + grprFitSigma ->GetY()[islice]*0.001;

  }

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (MCclosure) {
    //   set 1 iter params
    if (iIter>=0) SetParams1stIterationMC(pt,h1D,cleanParams);
    if (iIter>=1) SetParamsMC(islice,pt);
    if (iIter>=2){
      if ((pt>=0.32 && pt<=0.6) || (pt>=0.8 && pt<=1.3) ) elAmpMax  = grelFitAmpSmooth->GetY()[islice];
      if (pt>=0.3) {
	piAmpMin = grpiCleanAmp->GetY()[islice] - grpiCleanAmp->GetY()[islice]*0.1;
	piAmpMax = grpiCleanAmp->GetY()[islice] + grpiCleanAmp->GetY()[islice]*0.1;
	kaAmpMin = grkaCleanAmp->GetY()[islice] - grkaCleanAmp->GetY()[islice]*0.1;
	kaAmpMax = grkaCleanAmp->GetY()[islice] + grkaCleanAmp->GetY()[islice]*0.1;
	prAmpMin = grprCleanAmp->GetY()[islice] - grprCleanAmp->GetY()[islice]*0.1;
	prAmpMax = grprCleanAmp->GetY()[islice] + grprCleanAmp->GetY()[islice]*0.1;
      }
    }

    if (iIter>=3) {
      elAmpMin = grelFitAmp->GetY()[islice] - grelFitAmp->GetY()[islice]*0.01;
      elAmpMax = grelFitAmp->GetY()[islice] + grelFitAmp->GetY()[islice]*0.01;
      piAmpMin = grpiFitAmp->GetY()[islice] - grpiFitAmp->GetY()[islice]*0.001;
      piAmpMax = grpiFitAmp->GetY()[islice] + grpiFitAmp->GetY()[islice]*0.001;
      kaAmpMin = grkaFitAmp->GetY()[islice] - grkaFitAmp->GetY()[islice]*0.001;
      kaAmpMax = grkaFitAmp->GetY()[islice] + grkaFitAmp->GetY()[islice]*0.001;
      prAmpMin = grprFitAmp->GetY()[islice] - grprFitAmp->GetY()[islice]*0.001;
      prAmpMax = grprFitAmp->GetY()[islice] + grprFitAmp->GetY()[islice]*0.001;
    }

    // Use FreeKS from MC samples to get precise fits
    if (iIter>=4) {
      FixParamsMC(islice,pt);

      total->SetParLimits(4,skewMin,skewMax);    // electron
      total->FixParameter(9,grpiFreeKSSkew->GetY()[islice]);    // pion
      total->FixParameter(14,grkaFreeKSSkew->GetY()[islice]);   // kaon
      total->FixParameter(19,grprFreeKSSkew->GetY()[islice]);   // proton

      total->SetParLimits(3,kurtosisMin,kurtosisMax);   // electron
      total->FixParameter(8,grpiFreeKSKurtosis->GetY()[islice]);   // pion
      total->FixParameter(13,grkaFreeKSKurtosis->GetY()[islice]);  // kaon
      total->FixParameter(18,grprFreeKSKurtosis->GetY()[islice]);  // proton
    }
  }

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

}
// -------------------------------------------------------------------------------------------------------
void  SetParamsForKSEstimate(Int_t iks, TF1 *total, TH1D *h1DKS, Double_t *arrMean, Double_t *arrSigma)
{

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

  elMeanMin = TMath::Max(30.,(arrMean[0]-2.*arrSigma[0]));
  piMeanMin = TMath::Max(30.,(arrMean[1]-2.*arrSigma[1]));
  kaMeanMin = TMath::Max(30.,(arrMean[2]-2.*arrSigma[2]));
  prMeanMin = TMath::Max(30.,(arrMean[3]-2.*arrSigma[3]));

  elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+2.*arrSigma[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+2.*arrSigma[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+2.*arrSigma[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+2.*arrSigma[3]));

  elSigmaMin = arrSigma[0]/2.;
  piSigmaMin = arrSigma[1]/2.;
  kaSigmaMin = arrSigma[2]/2.;
  prSigmaMin = arrSigma[3]/2.;

  elSigmaMax = arrSigma[0]*2.;
  piSigmaMax = arrSigma[1]*2.;
  kaSigmaMax = arrSigma[2]*2.;
  prSigmaMax = arrSigma[3]*2.;

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

}
// -------------------------------------------------------------------------------------------------------
void  SetParams1stIteration(Double_t pt, TH1D *h1D, Double_t *arrMean, Double_t *arrSigma,Double_t *arrMeanWindow, Double_t *cleanParams)
{

  //
    //    Set mean amp and sigma par limits for the first iteration
  //

  Double_t meantmparr[nParticle];
  for (Int_t i=0;i<nParticle;i++) meantmparr[i] = arrMean[i];


  meantmparr[0] = ( (pt>0.35) && TMath::Abs(arrMean[0]-cleanParams[0])>12 ) ?  arrMean[0] : cleanParams[0]-1.5;
  meantmparr[1] = ( (pt>0.35) && TMath::Abs(arrMean[1]-cleanParams[1])>12 ) ?  arrMean[1] : cleanParams[1];
  meantmparr[2] = ( (pt>0.40) && TMath::Abs(arrMean[2]-cleanParams[2])>12 ) ?  arrMean[2] : cleanParams[2];
  meantmparr[3] = ( (pt>0.50) && TMath::Abs(arrMean[3]-cleanParams[3])>12 ) ?  arrMean[3] : cleanParams[3];

//   Part to make setting for fit params
  maxBin  = h1D->GetBinContent(h1D->GetMaximumBin());
  maxBin0 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[0])),maxBin);
  maxBin1 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[1])),maxBin);
  maxBin2 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[2])),maxBin);
  maxBin3 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[3])),maxBin);

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

  elMeanMin = TMath::Max(30.,(meantmparr[0]-arrMeanWindow[0]));
  piMeanMin = TMath::Max(30.,(meantmparr[1]-arrMeanWindow[1]));
  kaMeanMin = TMath::Max(30.,(meantmparr[2]-arrMeanWindow[2]));
  prMeanMin = TMath::Max(30.,(meantmparr[3]-arrMeanWindow[3]));

  elMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[0]+arrMeanWindow[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[1]+arrMeanWindow[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[2]+arrMeanWindow[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[3]+arrMeanWindow[3]));

  elSigmaMin = (KSestimation) ? arrSigma[0]/2. : arrSigma[0]+2;
  piSigmaMin = (KSestimation) ? arrSigma[1]/2. : arrSigma[1]+2;
  kaSigmaMin = (KSestimation) ? arrSigma[2]/2. : arrSigma[2]+1;
  prSigmaMin = (KSestimation) ? arrSigma[3]/2. : arrSigma[3]+1;

  elSigmaMax = arrSigma[0]*2.;
  piSigmaMax = arrSigma[1]*2.;
  kaSigmaMax = arrSigma[2]*2.;
  prSigmaMax = arrSigma[3]*2.;

}
// -------------------------------------------------------------------------------------------------------
TF1  *FitParamInRange(TGraphErrors* gr, TString funcType, Double_t min, Double_t max)
{

  //
  // Fit proton TGraph with pol2
  //

  TF1 *f = new TF1("f",funcType,min,max);
  gr->Fit(f,"QNR+");
  f->SetRange(0.2,ptRangeUp);
  return f;

}
// -------------------------------------------------------------------------------------------------------
TH1D *FuncToHist(TF1 *f, TH1D * h, Int_t islice, Double_t *arrMean, Double_t *arrSigma)
{

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
Double_t ComputeGaussIntegral(Double_t parAmp,Double_t parMean,Double_t parSigma, Double_t parKurtosis, Double_t parSkew)
{

  // Calculate integral of the gauss functions

  Double_t min = TMath::Max(0.,parMean-parSigma*5);
  Double_t max = TMath::Min((Double_t)dEdxMax,parMean+parSigma*5);

  TF1 f("f",fitFunctionGenGaus,min,max);
  // f.SetNpx(nBinsInLookUpTable);
  f.SetParameter(0,parAmp);
  f.SetParameter(1,parMean);
  f.SetParameter(2,parSigma);
  f.SetParameter(3,parKurtosis);
  f.SetParameter(4,parSkew);

  return f.Integral(min,max);

}
// -------------------------------------------------------------------------------------------------------
TH1D *GetClean1DSlice(TH2D *h2Clean, TString parName, Int_t ptDownBin, Int_t ptUpBin)
{

  //
  //  get 1D momentum slice
  //

  TH1D *h1Clean   = (TH1D*)h2Clean->ProjectionY(Form(parName,ptDownBin),ptDownBin,ptDownBin);
  h1Clean->Rebin(rebinFactor); h1Clean->Scale(1./rebinFactor);
  Double_t binWidth = h1Clean->GetXaxis()->GetBinWidth(50);
  h1Clean->Scale(1./binWidth);
  Double_t mean = h1Clean->GetMean();
  Double_t rms  = h1Clean->GetRMS();
  Double_t rangeDown = TMath::Max((Double_t)dEdxMin,mean-5.*rms);
  Double_t rangeUp   = TMath::Min((Double_t)dEdxMax,mean+5.*rms);
  h1Clean->GetXaxis()->SetRangeUser(rangeDown,rangeUp);

  return h1Clean;
}
// -------------------------------------------------------------------------------------------------------
void FitParticleSampleForKSEstimate(Int_t iks, TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy)
{

  //
  // Fit Clean sample Slice
  //

  TVirtualFitter::SetMaxIterations(100000);

    // helper numbers for setparams
  maxBin        = hClean->GetBinContent(hClean->GetMaximumBin());  if ((maxBin<20.)) return;
  Double_t mean = (MCclosure) ? hClean->GetMean() : arrMean[pSpecy];
  Double_t rms  = (MCclosure) ? hClean->GetRMS()  : arrSigma[pSpecy];
  Double_t sampleFitWindow = 3.5;

    // fit function
  TF1 *asyGaus = new TF1("asyGaus",fitFunctionGenGaus,mean-sampleFitWindow*rms,mean+4.2*rms);
  asyGaus->SetParNames("Abundance","Mean","Sigma","Kurtosis","Skewness");
  asyGaus->SetParLimits(0,maxBin/20.,maxBin);
  asyGaus->SetParLimits(1,mean-rms,mean+rms);
  asyGaus->SetParLimits(2,rms/2. ,3.*rms);

  if (iks==0){
    asyGaus->SetParLimits(3,kurtosisMin,kurtosisMax);
    asyGaus->SetParLimits(4,skewMin,skewMax);
  } else {
    asyGaus->FixParameter(3,kurtosisFixClean[pSpecy]);
    asyGaus->SetParLimits(4,skewMin,skewMax);
  }

  // retrieve fit parameters
  hClean->Fit(asyGaus,"QNMR");
  asyGaus->SetRange(mean-5*rms,mean+5*rms);
  hClean->GetXaxis()->SetTitle("TPC d#it{E}/d#it{x} Signal (a.u.)");
  hClean->GetYaxis()->SetTitle("entries");
  parSample[0] = (asyGaus->GetParameter(0)>1.) ? (asyGaus->GetParameter(0)) : 1.;      // avoid floating point exception for muons
  parSample[1] = asyGaus->GetParameter(1);
  parSample[2] = asyGaus->GetParameter(2);
  parSample[3] = asyGaus->GetParameter(3);
  parSample[4] = asyGaus->GetParameter(4);
  parSample[5] = asyGaus->Integral(parSample[1]-parSample[2]*4.,parSample[1]+parSample[2]*4.);
  if (asyGaus->GetNDF()>1e-4) parSample[6] =  asyGaus->GetChisquare()/asyGaus->GetNDF();

  delete asyGaus;
}
// -------------------------------------------------------------------------------------------------------
void Fit1DSliceForKSEstimate(Int_t iks, TH1D *h1DKS, Double_t *pars, Double_t *arrMean, Double_t *arrSigma)
{

    // Fit part
  maxBin = h1DKS->GetBinContent(h1DKS->GetMaximumBin());
  Double_t fitWinMean  = 1.2;

  elMeanMin = TMath::Max(30.,(arrMean[0]-2.*arrSigma[0]));                                   // electron
  piMeanMin = TMath::Max(30.,(arrMean[1]-fitWinMean*arrSigma[1]));                           // pion
  kaMeanMin = TMath::Max(30.,(arrMean[2]-fitWinMean*arrSigma[2]));                           // kaon
  prMeanMin = TMath::Max(30.,(arrMean[3]-fitWinMean*arrSigma[3]));                           // proton

  elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+2.*arrSigma[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+fitWinMean*arrSigma[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+fitWinMean*arrSigma[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+fitWinMean*arrSigma[3]));

  TF1 *g1   = new TF1("g1",fitFunctionGenGaus,elMeanMin,elMeanMax);  // g1->SetNpx(nBinsInLookUpTable);
  TF1 *g2   = new TF1("g2",fitFunctionGenGaus,piMeanMin,piMeanMax);  // g2->SetNpx(nBinsInLookUpTable);
  TF1 *g3   = new TF1("g3",fitFunctionGenGaus,kaMeanMin,kaMeanMax);  // g3->SetNpx(nBinsInLookUpTable);
  TF1 *g4   = new TF1("g4",fitFunctionGenGaus,prMeanMin,prMeanMax);  // g4->SetNpx(nBinsInLookUpTable);

    // restrict the parameters of gaus functions of individual particles
  g1->SetParLimits(0,0.,maxBin);
  g2->SetParLimits(0,0.,maxBin);
  g3->SetParLimits(0,0.,maxBin);
  g4->SetParLimits(0,0.,maxBin);

  g1->SetParLimits(1,elMeanMin,elMeanMax);
  g2->SetParLimits(1,piMeanMin,piMeanMax);
  g3->SetParLimits(1,kaMeanMin,kaMeanMax);
  g4->SetParLimits(1,prMeanMin,prMeanMax);

  g1->SetParLimits(2,arrSigma[0]/2.,arrSigma[0]*2.);
  g2->SetParLimits(2,arrSigma[1]/2.,arrSigma[1]*2.);
  g3->SetParLimits(2,arrSigma[2]/2.,arrSigma[2]*2.);
  g4->SetParLimits(2,arrSigma[3]/2.,arrSigma[3]*2.);

  (iks==0)? g1->SetParLimits(3,kurtosisMin,kurtosisMax) : g1->FixParameter(3,kurtosisFix[kElectron]);
  (iks==0)? g2->SetParLimits(3,kurtosisMin,kurtosisMax) : g2->FixParameter(3,kurtosisFix[kPion]    );
  (iks==0)? g3->SetParLimits(3,kurtosisMin,kurtosisMax) : g3->FixParameter(3,kurtosisFix[kKaon]    );
  (iks==0)? g4->SetParLimits(3,kurtosisMin,kurtosisMax) : g4->FixParameter(3,kurtosisFix[kProton]  );

  g1->SetParLimits(4,skewMin,skewMax);
  g2->SetParLimits(4,skewMin,skewMax);
  g3->SetParLimits(4,skewMin,skewMax);
  g4->SetParLimits(4,skewMin,skewMax);

    // Fit h1D for individual particles
  if (arrMean[0]>=10) h1DKS->Fit(g1,"QNR");
  if (arrMean[1]>=10) h1DKS->Fit(g2,"QNR+");
  if (arrMean[2]>=10) h1DKS->Fit(g3,"QNR+");
  if (arrMean[3]>=10) h1DKS->Fit(g4,"QNR+");

    // Get Fit parameters from the individual fits and set for the total fit
  Double_t par[20] = {0};
  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[5]);
  g3->GetParameters(&par[10]);
  g4->GetParameters(&par[15]);

  TF1 *total = new TF1("total","g1+g2+g3+g4",dEdxMin,dEdxMax);
  // total->SetNpx(nBinsInLookUpTable);
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
  pars[20] = totChi2;

  delete total;
  delete g1;
  delete g2;
  delete g3;
  delete g4;
}
// -------------------------------------------------------------------------------------------------------
TMatrixD *FitParticleSample(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy)
{

  //
  // Fit Clean sample slice
  //

//   TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter::SetMaxIterations(100000);
  maxBin = hClean->GetBinContent(hClean->GetMaximumBin());
  Double_t mean = (MCclosure) ? hClean->GetMean() : arrMean[pSpecy];
  Double_t rms  = (MCclosure) ? hClean->GetRMS()  : arrSigma[pSpecy];
  Double_t sampleFitWindow = 3.5;
  Double_t integral        = hClean->Integral(hClean->FindBin(mean-sampleFitWindow*rms),hClean->FindBin(mean+4.2*rms));

  Bool_t histCheck = ( (mean>20 && mean<1020) && (rms>0 && rms<200) && (maxBin>=1 && maxBin<1e+30) && (integral>=5 && integral<1e+30));
  if (!histCheck) {   // there is too few entries, fit would fail
    parSample[0] = 0.;
    parSample[1] = 1.;
    parSample[2] = 1.;
    parSample[3] = 1.;
    parSample[4] = 1.;
    parSample[5] = 1.;
    return new TMatrixD(5,5);
  }

  // prepare normalised samples
  if (normalisedCleanSamples && integral > 1e-4) {
    hClean->Scale(1./integral);
    maxBin /= integral;
  }

  // fit function
  TF1 *asyGaus = new TF1("asyGaus",fitFunctionGenGaus,mean-sampleFitWindow*rms,mean+4.2*rms);
  asyGaus->SetParNames("Amplitude","Mean","Sigma","Kurtosis","Skewness");
  // asyGaus->SetNpx(nBinsInLookUpTable);
  asyGaus ->SetLineColor(2);
  asyGaus ->SetLineWidth(2);

  // Set Parameters
  asyGaus->SetParLimits(0,maxBin/20.,maxBin);
  asyGaus->SetParLimits(1,mean-rms,mean+rms);
  asyGaus->SetParLimits(2,rms/2. ,3.*rms);
  asyGaus->SetParLimits(3,kurtosisMin,kurtosisMax);
  asyGaus->SetParLimits(4,skewMin,skewMax);
  if ((fixedK || fixedS)) SetFixedKSparameterForCleanSampleFit(asyGaus,pSpecy);

  // TH1D *hMatrix = new TH1D(); hClean->Copy(*hMatrix);
  TH1D *hMatrix = (TH1D*)hClean->Clone();
  // TF1 *fMatrix = new TF1(); asyGaus->Copy(*fMatrix);
  TF1 *fMatrix = (TF1*)asyGaus->Clone();

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
void FitParticleSampleFreeKS(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy)
{

  //
  // Fit Clean sample slice
  //

  TVirtualFitter::SetMaxIterations(100000);
  maxBin = hClean->GetBinContent(hClean->GetMaximumBin());
  Double_t mean = (MCclosure) ? hClean->GetMean() : arrMean[pSpecy];
  Double_t rms  = (MCclosure) ? hClean->GetRMS()  : arrSigma[pSpecy];
  Double_t sampleFitWindow = 3.5;
  Double_t integral        = hClean->Integral(hClean->FindBin(mean-sampleFitWindow*rms),hClean->FindBin(mean+4.2*rms));

  Bool_t histCheck = ( (mean>20 && mean<1020) && (rms>0 && rms<200) && (maxBin>=1 && maxBin<1e+30) && (integral>=5 && integral<1e+30));
  if (!histCheck) {   // there is too few entries, fit would fail
    parSample[0] = 0.;
    parSample[1] = 1.;
    parSample[2] = 1.;
    parSample[3] = 1.;
    parSample[4] = 1.;
    parSample[5] = 1.;
    return;
  }

  // fit function
  TF1 *asyGaus = new TF1("asyGaus",fitFunctionGenGaus,mean-sampleFitWindow*rms,mean+4.2*rms);
  asyGaus->SetParNames("Amplitude","Mean","Sigma","Kurtosis","Skewness");
  // asyGaus->SetNpx(nBinsInLookUpTable);
  asyGaus ->SetLineColor(kBlue);
  asyGaus ->SetLineWidth(2);

  // Set Parameters
  asyGaus->SetParLimits(0,maxBin/20.,maxBin);
  asyGaus->SetParLimits(1,mean-2*rms,mean+2*rms);
  asyGaus->SetParLimits(2,rms/3. ,3.*rms);
  asyGaus->SetParLimits(3,kurtosisMin,kurtosisMax);
  asyGaus->SetParLimits(4,skewMin,skewMax);

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
// -------------------------------------------------------------------------------------------------------
void GetCleanExParams(TH1D *hClean, Double_t *arrSigma, Double_t *arrCleanSigma, Double_t *arrCleanMean, ParticleType pSpecy)
{

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
Int_t GetMaxIndex(TGraphErrors *grIndex)
{

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
Int_t GetPtIndex(TGraphErrors *grIndex,Double_t p)
{

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
void ResetParametersOfEachFunction(TF1 *g1,TF1 *g2,TF1 *g3,TF1 *g4,TF1 *total,TH1D *h1D,const Double_t pt)
{

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

  // Set colors
  g1->SetLineColor(4);
  g2->SetLineColor(kOrange-1);
  g3->SetLineColor(kGreen+2);
  g4->SetLineColor(kMagenta+1);

  g1->SetFillColor(4);
  g2->SetFillColor(kOrange-1);
  g3->SetFillColor(kGreen+2);
  g4->SetFillColor(kMagenta+1);

  g1->SetFillStyle(3004); g1->Draw("fchist");
  g2->SetFillStyle(3005); g2->Draw("fchist");
  g3->SetFillStyle(3006); g3->Draw("fchist");
  g4->SetFillStyle(3007); g4->Draw("fchist");

  // Set ranges
  g1->SetRange(0,h1D->GetNbinsX());
  g2->SetRange(0,h1D->GetNbinsX());
  g3->SetRange(0,h1D->GetNbinsX());
  g4->SetRange(0,h1D->GetNbinsX());

  // add TLegend and Tlines
  TLegend *leg = new TLegend(0.75, 0.6, 0.85, 0.85);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);  // White
  leg->SetBorderSize(0); // No border
  leg->SetHeader(Form("#it{p} = %3.2f GeV/#it{c}",pt));

  leg->AddEntry(h1D," Data","LPE");
  leg->AddEntry(total," Total Fit","L");
  leg->AddEntry(g1   ," e"  ,"F");
  leg->AddEntry(g2   ," #pi","F");
  leg->AddEntry(g3   ," K"  ,"F");
  leg->AddEntry(g4   ,"p"   ,"F");

  // Attach functions and TLegend to the histogram
  h1D->GetListOfFunctions()->Add(g1);
  h1D->GetListOfFunctions()->Add(g2);
  h1D->GetListOfFunctions()->Add(g3);
  h1D->GetListOfFunctions()->Add(g4);
  h1D->GetListOfFunctions()->Add(leg);

}
// -------------------------------------------------------------------------------------------------------
Double_t GetClosestParticleMean(Double_t sig, Double_t ref, Double_t gr1, Double_t gr2, Double_t gr3)
{
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
TGraphErrors * GraphSmooth(TGraphErrors * gr)
{

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
TH1D *GraphShade(TGraphErrors *grFitParam, Double_t percent)
{

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
TGraphErrors * HistToGraphErrors(TH1D * h)
{

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
void ProduceGraphs(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grTmp, TGraph **grSmoothTmp)
{
  //
  // Produce graphs for fit and Clean
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
    if ( branchName.Contains("normNDF") || branchName.Contains("Matrix") || branchName.Contains("ExC") ) continue;
    if ( branchName == "slice"   || branchName == "p" ) continue;
    if ( branchName == "eta"     || branchName == "cent" ) continue;

    // if ( branchName.Contains("eta") || branchName.Contains("cent") || branchName.Contains("slice") ) continue;

    Double_t bassLevel = smoothLevel;
    if (branchName.Contains("elFitAmp")) bassLevel=5.;
    if (branchName.Contains("kaFitAmp")) bassLevel=8.;
    if (branchName.Contains("kaFitSigma")) bassLevel=9.;

    // create the graph object
    t->Draw(drawStr,"","goff");
    grTmp[objIndex] = new TGraphErrors(nSlice,t->GetV2(),t->GetV1());
    grTmp[objIndex]->SetName(branchName);
    grTmp[objIndex]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    grTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grTmp[objIndex]->SetLineColor(kBlack);
    grTmp[objIndex]->SetMarkerStyle(7);
    grTmp[objIndex]->SetMarkerColor(kBlack);

    // create smooth graph
    gs      = new TGraphSmooth("normal");
    if ( branchName.Contains("Kurtosis") || branchName.Contains("Skew") ) {
      grSmoothTmp[objIndex] = (TGraph*)grTmp[objIndex]->Clone();
    } else {
      grSmoothTmp[objIndex] = gs->SmoothSuper((TGraph*)grTmp[objIndex],"",bassLevel);
    }
    grSmoothTmp[objIndex]->SetName(branchName.Append("Smooth"));
    grSmoothTmp[objIndex]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
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
void ApplyOutlierSmooth(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grTmp, TGraph **grSmoothTmp)
{
  //
  // Produce graphs fit and Clean
  //


  Int_t nEntries = barr->GetEntriesFast();
  cout << "  in the  ApplyOutlierSmooth " << nEntries << endl;
  Int_t objIndex = 0;
  for (Int_t i=0; i<nEntries; i++){

    TString branchName = barr->At(i)->GetName();
    TString drawStr    = barr->At(i)->GetName();
    drawStr.Append(":p");

    // look at only amplitude Mean and Sigma
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
    grTmp[objIndex]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    grTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grTmp[objIndex]->SetLineColor(kBlack);
    grTmp[objIndex]->SetMarkerStyle(7);
    grTmp[objIndex]->SetMarkerColor(kBlack);

    // create smooth graph
    gs      = new TGraphSmooth("normal");
    grSmoothTmp[objIndex] = gs->SmoothSuper((TGraph*)grTmp[objIndex],"",9);
    grSmoothTmp[objIndex]->SetName(branchName.Append("Smooth"));
    grSmoothTmp[objIndex]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
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
TGraphErrors * ApplyScalingWrtExpected(const Int_t nSlice, TString exp, TString fit, const Double_t safeP, TObjArray * arr, TTree *t)
{
  //
  // Return scaled graph wrt to Expected
  //

  // get the graphs to be used in scaling
  TGraphErrors *grexp = (TGraphErrors*)arr->FindObject(exp);
  TGraphErrors *grfit = (TGraphErrors*)arr->FindObject(fit);

  // Find the index corresponding to safeP value
  Int_t index = GetPtIndex(grexp,safeP);

  // calculate scaling factor
  Double_t scaleFactor = grexp->GetY()[index]/grfit->GetY()[index];

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
TGraphErrors * ApplyScalingWrtSample(const Int_t nSlice, TString sample, TString fit, const Double_t safeP, TObjArray * arr, TObjArray * arrSample, TTree *tSample, TGraph ** grScaledSmooth,Int_t k)
{
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
    grScaled->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    grScaled->GetYaxis()->SetTitle(newname);
    grScaled->SetLineColor(kGreen);
    grScaled->SetMarkerStyle(7);
    grScaled->SetMarkerColor(kGreen);

    gs = new TGraphSmooth("normal");
    grScaledSmooth[k] = gs->SmoothSuper((TGraph*)grScaled,"",4);
    grScaledSmooth[k]->SetName(newname.Append("Smooth"));
    grScaledSmooth[k]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    grScaledSmooth[k]->GetYaxis()->SetTitle(newname);
    grScaledSmooth[k]->SetLineColor(kRed);
    grScaledSmooth[k]->SetMarkerStyle(7);
    grScaledSmooth[k]->SetMarkerColor(kRed);


    return grScaled;
}
// -------------------------------------------------------------------------------------------------------
TCanvas *GetFitResCanvas(TObjArray *arr)
{

  //
//   Fit Results Canvas
  //

  TCanvas *can = new TCanvas("can", "can", 1400, 900);
  can->Divide(5,4);

  // ++++++++++++++++++  Sigmas +++++++++++++++++
  can->cd(1);
  can->GetPad(1)->SetGrid();
  TGraphErrors *elsigma       = (TGraphErrors*)arr->FindObject("elFitSigma");
  TGraphErrors *elsigmaSmooth = (TGraphErrors*)arr->FindObject("elFitSigmaSmooth");
  elsigma->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  elsigma->GetYaxis()->SetTitle("electron sigma");
  elsigma->GetYaxis()->SetTitleOffset(1.6);
  elsigmaSmooth->SetLineColor(kRed);
  elsigma->Draw("alp");
  elsigmaSmooth->Draw("lp");

  can->cd(2);
  can->GetPad(2)->SetGrid();
  TGraphErrors *pisigma       = (TGraphErrors*)arr->FindObject("piFitSigma");
  TGraphErrors *pisigmaSmooth = (TGraphErrors*)arr->FindObject("piFitSigmaSmooth");
  pisigma->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  pisigma->GetYaxis()->SetTitle("pion sigma");
  pisigma->GetYaxis()->SetTitleOffset(1.6);
  pisigmaSmooth->SetLineColor(kRed);
  pisigma->Draw("alp");
  pisigmaSmooth->Draw("lp");

  can->cd(3);
  can->GetPad(3)->SetGrid();
  TGraphErrors *kasigma       = (TGraphErrors*)arr->FindObject("kaFitSigma");
  TGraphErrors *kasigmaSmooth = (TGraphErrors*)arr->FindObject("kaFitSigmaSmooth");
  kasigma->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  kasigma->GetYaxis()->SetTitle("kaon sigma");
  kasigma->GetYaxis()->SetTitleOffset(1.6);
  kasigmaSmooth->SetLineColor(kRed);
  kasigma->Draw("alp");
  kasigmaSmooth->Draw("lp");

  can->cd(4);
  can->GetPad(4)->SetGrid();
  TGraphErrors *prsigma       = (TGraphErrors*)arr->FindObject("prFitSigma");
  TGraphErrors *prsigmaSmooth = (TGraphErrors*)arr->FindObject("prFitSigmaSmooth");
  prsigma->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  prsigma->GetYaxis()->SetTitle("proton sigma");
  prsigma->GetYaxis()->SetTitleOffset(1.6);
  prsigmaSmooth->SetLineColor(kRed);
  prsigma->Draw("alp");
  prsigmaSmooth->Draw("lp");

  can->cd(5);
  can->GetPad(5)->SetGrid();
  TGraphErrors *chi2          = (TGraphErrors*)arr->FindObject("totChi2");
  chi2->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
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
  elmean->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  elmean->GetYaxis()->SetTitle("electron mean");
  elmean->GetYaxis()->SetTitleOffset(1.6);
  elmeanSmooth->SetLineColor(kRed);
  elmean->Draw("alp");
  elmeanSmooth->Draw("lp");

  can->cd(7);
  can->GetPad(7)->SetGrid();
  TGraphErrors *pimean       = (TGraphErrors*)arr->FindObject("piFitMean");
  TGraphErrors *pimeanSmooth = (TGraphErrors*)arr->FindObject("piFitMeanSmooth");
  pimean->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  pimean->GetYaxis()->SetTitle("pion mean");
  pimean->GetYaxis()->SetTitleOffset(1.6);
  pimeanSmooth->SetLineColor(kRed);
  pimean->Draw("alp");
  pimeanSmooth->Draw("lp");

  can->cd(8);
  can->GetPad(8)->SetGrid();
  TGraphErrors *kamean       = (TGraphErrors*)arr->FindObject("kaFitMean");
  TGraphErrors *kameanSmooth = (TGraphErrors*)arr->FindObject("kaFitMeanSmooth");
  kamean->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  kamean->GetYaxis()->SetTitle("kaon mean");
  kamean->GetYaxis()->SetTitleOffset(1.6);
  kameanSmooth->SetLineColor(kRed);
  kamean->Draw("alp");
  kameanSmooth->Draw("lp");

  can->cd(9);
  can->GetPad(9)->SetGrid();
  TGraphErrors *prmean       = (TGraphErrors*)arr->FindObject("prFitMean");
  TGraphErrors *prmeanSmooth = (TGraphErrors*)arr->FindObject("prFitMeanSmooth");
  prmean->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  prmean->GetYaxis()->SetTitle("proton mean");
  prmean->GetYaxis()->SetTitleOffset(1.6);
  prmeanSmooth->SetLineColor(kRed);
  prmean->Draw("alp");
  prmeanSmooth->Draw("lp");

  can->cd(10);
  can->GetPad(10)->SetGrid();
  TGraphErrors *maxRes          = (TGraphErrors*)arr->FindObject("maxRes");
  maxRes->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
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
  elamp->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  elamp->GetYaxis()->SetTitle("electron amplitude");
  elamp->GetYaxis()->SetTitleOffset(1.6);
  elampSmooth->SetLineColor(kRed);
  elamp->Draw("alp");
  elampSmooth->Draw("lp");

  can->cd(12);
  can->GetPad(12)->SetGrid();
  TGraphErrors *piamp       = (TGraphErrors*)arr->FindObject("piFitAmp");
  TGraphErrors *piampSmooth = (TGraphErrors*)arr->FindObject("piFitAmpSmooth");
  piamp->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  piamp->GetYaxis()->SetTitle("pion amplitude");
  piamp->GetYaxis()->SetTitleOffset(1.6);
  piampSmooth->SetLineColor(kRed);
  piamp->Draw("alp");
  piampSmooth->Draw("lp");

  can->cd(13);
  can->GetPad(13)->SetGrid();
  TGraphErrors *kaamp       = (TGraphErrors*)arr->FindObject("kaFitAmp");
  TGraphErrors *kaampSmooth = (TGraphErrors*)arr->FindObject("kaFitAmpSmooth");
  kaamp->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  kaamp->GetYaxis()->SetTitle("kaon amplitude");
  kaamp->GetYaxis()->SetTitleOffset(1.6);
  kaampSmooth->SetLineColor(kRed);
  kaamp->Draw("alp");
  kaampSmooth->Draw("lp");

  can->cd(14);
  can->GetPad(14)->SetGrid();
  TGraphErrors *pramp       = (TGraphErrors*)arr->FindObject("prFitAmp");
  TGraphErrors *prampSmooth = (TGraphErrors*)arr->FindObject("prFitAmpSmooth");
  pramp->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  pramp->GetYaxis()->SetTitle("proton amplitude");
  pramp->GetYaxis()->SetTitleOffset(1.6);
  prampSmooth->SetLineColor(kRed);
  pramp->Draw("alp");
  prampSmooth->Draw("lp");

  can->cd(15);
  can->GetPad(15)->SetGrid();
  TGraphErrors *maxPer          = (TGraphErrors*)arr->FindObject("maxPer");
  maxPer->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
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
  elint->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  elint->GetYaxis()->SetTitle("electron yield");
  elint->GetYaxis()->SetTitleOffset(1.6);
  elint->Draw("alp");

  can->cd(17);
  can->GetPad(17)->SetGrid();
  TGraphErrors *piint       = (TGraphErrors*)arr->FindObject("piFitInt");
  piint->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  piint->GetYaxis()->SetTitle("pion yield");
  piint->GetYaxis()->SetTitleOffset(1.6);
  piint->Draw("alp");

  can->cd(18);
  can->GetPad(18)->SetGrid();
  TGraphErrors *kaint       = (TGraphErrors*)arr->FindObject("kaFitInt");
  kaint->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  kaint->GetYaxis()->SetTitle("kaon yield");
  kaint->GetYaxis()->SetTitleOffset(1.6);
  kaint->Draw("alp");

  can->cd(19);
  can->GetPad(19)->SetGrid();
  TGraphErrors *print       = (TGraphErrors*)arr->FindObject("prFitInt");
  print->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  print->GetYaxis()->SetTitle("proton yield");
  print->GetYaxis()->SetTitleOffset(1.6);
  print->Draw("alp");

  can->cd(20);
  can->GetPad(20)->SetGrid();
  TGraphErrors *ksValue     = (TGraphErrors*)arr->FindObject("ksTest");
  ksValue->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  ksValue->GetYaxis()->SetTitle("ksTest");
  ksValue->GetYaxis()->SetTitleOffset(1.6);
  ksValue->SetMarkerStyle(7);
  ksValue->SetMarkerColor(kMagenta);
  ksValue->SetLineColor(kMagenta);
  ksValue->Draw("alp");

  return can;

}
// -------------------------------------------------------------------------------------------------------
TGraphErrors *ProduceWindowGraphs(Int_t nSlice, Int_t index, TString str, TTree *tree)
{

  //
  //  Prepare graphs around the expected mean of each particle
  //

  tree->Draw(str,"","goff");
  TGraphErrors *grwindow = new TGraphErrors(nSlice,tree->GetV2(),tree->GetV1());
  grwindow->SetName(str);
  grwindow->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  grwindow->GetYaxis()->SetTitle(str);
  grwindow->SetMarkerStyle(7);
  if (index%3!=0)  grwindow->SetLineColor(kRed);
  return grwindow;

}
// -------------------------------------------------------------------------------------------------------
void ApplyInitialFits(Int_t islice, Double_t pt, TH1D *h1D, Double_t *arrSigma, Double_t *arrMean)
{

  //
    // Apply initial fits and dump output to streamer
  //

  maxBin = h1D->GetBinContent(h1D->GetMaximumBin());

  maxBin0 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[0])),maxBin);
  maxBin1 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[1])),maxBin);
  maxBin2 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[2])),maxBin);
  maxBin3 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[3])),maxBin);

  TH1D * h1FitWindows  = (TH1D*)h1D->Clone(); h1FitWindows->SetName(Form("FitWindows_%d",islice));
  SetStyleFor1DSlice(h1FitWindows);
  TH1D * h1InitialFits = (TH1D*)h1D->Clone(); h1InitialFits->SetName(Form("InitialFits_%d",islice));
  SetStyleFor1DSlice(h1InitialFits);

  Double_t fitWinMean  = 1.2;
  Double_t fitWinSigma = 0.5;

    // Fix the mean position between +-1 sigma
  elMeanMin = TMath::Max(30.,(arrMean[0]-fitWinMean*arrSigma[0]));            // electron
  piMeanMin = TMath::Max(30.,(arrMean[1]-fitWinMean*arrSigma[1]));            // pion
  kaMeanMin = TMath::Max(30.,(arrMean[2]-fitWinMean*arrSigma[2]));            // kaon
  prMeanMin = TMath::Max(30.,(arrMean[3]-fitWinMean*arrSigma[3]));            // proton

  elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+fitWinMean*arrSigma[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+fitWinMean*arrSigma[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+fitWinMean*arrSigma[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+fitWinMean*arrSigma[3]));

  TF1 *g1b = new TF1("g1b","[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]/TMath::Sqrt(2)))",elMeanMin,elMeanMax);
  TF1 *g2b = new TF1("g2b","[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]/TMath::Sqrt(2)))",piMeanMin,piMeanMax);
  TF1 *g3b = new TF1("g3b","[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]/TMath::Sqrt(2)))",kaMeanMin,kaMeanMax);
  TF1 *g4b = new TF1("g4b","[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]/TMath::Sqrt(2)))",prMeanMin,prMeanMax);

  // g1b->SetNpx(nBinsInLookUpTable);
  // g2b->SetNpx(nBinsInLookUpTable);
  // g3b->SetNpx(nBinsInLookUpTable);
  // g4b->SetNpx(nBinsInLookUpTable);

  g1b->SetLineColor(4);         // Blue
  g2b->SetLineColor(kOrange-1);
  g3b->SetLineColor(kGreen+2);
  g4b->SetLineColor(kMagenta+1);         // Black

  g1b->SetParLimits(0,0.,maxBin);
  g2b->SetParLimits(0,0.,maxBin);
  g3b->SetParLimits(0,0.,maxBin);
  g4b->SetParLimits(0,0.,maxBin);

  g1b->SetParLimits(1,elMeanMin,elMeanMax);
  g2b->SetParLimits(1,piMeanMin,piMeanMax);
  g3b->SetParLimits(1,kaMeanMin,kaMeanMax);
  g4b->SetParLimits(1,prMeanMin,prMeanMax);

  g1b->SetParLimits(2,arrSigma[0]-fitWinSigma*arrSigma[0],arrSigma[0]+fitWinSigma*arrSigma[0]);
  g2b->SetParLimits(2,arrSigma[1]-fitWinSigma*arrSigma[1],arrSigma[1]+fitWinSigma*arrSigma[1]);
  g3b->SetParLimits(2,arrSigma[2]-fitWinSigma*arrSigma[2],arrSigma[2]+fitWinSigma*arrSigma[2]);
  g4b->SetParLimits(2,arrSigma[3]-fitWinSigma*arrSigma[3],arrSigma[3]+fitWinSigma*arrSigma[3]);


  g1b->FixParameter(3,0.);
  g2b->FixParameter(3,0.);
  g3b->FixParameter(3,0.);
  g4b->FixParameter(3,0.);


  if (arrMean[0]>0) h1FitWindows->Fit(g1b,"QR");
  if (arrMean[1]>0) h1FitWindows->Fit(g2b,"QR+");
  if (arrMean[2]>0) h1FitWindows->Fit(g3b,"QR+");
  if (arrMean[3]>0) h1FitWindows->Fit(g4b,"QR+");

  ///////////////////////////////////////////////////////////////////////////////////
    // add TLegend
  TLegend *leg = new TLegend(0.75, 0.6, 0.85, 0.85);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);  // White
//   leg->SetHeader(Form("p=%4.3f",pt));
  TLegend *legClone = (TLegend*)leg->Clone();

    // Plot lines for electrons
  TLine* lineg1 = new TLine(arrMean[0], 0., arrMean[0], maxBin0);
  lineg1     ->SetLineColor(4);
  lineg1     ->SetLineWidth(2);
  TLine* lineg1up = new TLine(elMeanMax, 0., elMeanMax, maxBin0);
  lineg1up   ->SetLineColor(4);
  lineg1up   ->SetLineWidth(1);
  lineg1up   ->SetLineStyle(2);
  TLine* lineg1down = new TLine(elMeanMin, 0., elMeanMin, maxBin0);
  lineg1down ->SetLineColor(4);
  lineg1down ->SetLineWidth(1);
  lineg1down ->SetLineStyle(2);
  if (arrMean[0]>0 && arrMean[0]<dEdxMax) {
    h1FitWindows->GetListOfFunctions()->Add(lineg1);
    h1FitWindows->GetListOfFunctions()->Add(lineg1up);
    h1FitWindows->GetListOfFunctions()->Add(lineg1down);
  }

    // Plot lines for Pions
  TLine* lineg2 = new TLine(arrMean[1], 0., arrMean[1], maxBin1);
  lineg2      ->SetLineColor(kOrange-1);
  lineg2      ->SetLineWidth(2);
  TLine* lineg2up = new TLine(piMeanMax, 0., piMeanMax, maxBin1);
  lineg2up    ->SetLineColor(kOrange-1);
  lineg2up    ->SetLineWidth(1);
  lineg2up    ->SetLineStyle(2);
  TLine* lineg2down = new TLine(piMeanMin, 0., piMeanMin, maxBin1);
  lineg2down  ->SetLineColor(kOrange-1);
  lineg2down  ->SetLineWidth(1);
  lineg2down  ->SetLineStyle(2);
  if (arrMean[1]>0 && arrMean[1]<dEdxMax) {
    h1FitWindows->GetListOfFunctions()->Add(lineg2);
    h1FitWindows->GetListOfFunctions()->Add(lineg2up);
    h1FitWindows->GetListOfFunctions()->Add(lineg2down);
  }

    // Plot lines for Kaon
  TLine* lineg3 = new TLine(arrMean[2], 0., arrMean[2], maxBin2);
  lineg3      ->SetLineColor(kGreen+2);
  lineg3      ->SetLineWidth(2);
  TLine* lineg3up = new TLine(kaMeanMax, 0., kaMeanMax, maxBin2);
  lineg3up    ->SetLineColor(kGreen+2);
  lineg3up    ->SetLineWidth(1);
  lineg3up    ->SetLineStyle(2);
  TLine* lineg3down = new TLine(kaMeanMin, 0., kaMeanMin, maxBin2);
  lineg3down  ->SetLineColor(kGreen+2);
  lineg3down  ->SetLineWidth(1);
  lineg3down  ->SetLineStyle(2);
  if (arrMean[2]>0 && arrMean[2]<dEdxMax) {
    h1FitWindows->GetListOfFunctions()->Add(lineg3);
    h1FitWindows->GetListOfFunctions()->Add(lineg3up);
    h1FitWindows->GetListOfFunctions()->Add(lineg3down);
  }

    // Plot lines for Proton
  TLine* lineg4 = new TLine(arrMean[3], 0., arrMean[3], maxBin3);
  lineg4     ->SetLineColor(kMagenta+1);
  lineg4     ->SetLineWidth(2);
  TLine* lineg4up = new TLine(prMeanMax, 0., prMeanMax, maxBin3);
  lineg4up   ->SetLineColor(kMagenta+1);
  lineg4up   ->SetLineWidth(1);
  lineg4up   ->SetLineStyle(2);
  TLine* lineg4down = new TLine(prMeanMin, 0., prMeanMin, maxBin3);
  lineg4down ->SetLineColor(kMagenta+1);
  lineg4down ->SetLineWidth(1);
  lineg4down ->SetLineStyle(2);
  if (arrMean[3]>0 && arrMean[3]<dEdxMax) {
    h1FitWindows->GetListOfFunctions()->Add(lineg4);
    h1FitWindows->GetListOfFunctions()->Add(lineg4up);
    h1FitWindows->GetListOfFunctions()->Add(lineg4down);
  }

  leg->AddEntry(h1FitWindows,"Data","LEP");
  leg->AddEntry(g1b,"e","L");
  leg->AddEntry(g2b,"#pi","L");
  leg->AddEntry(g3b,"K","L");
  leg->AddEntry(g4b,"p","L");

  h1FitWindows->GetListOfFunctions()->Add(g1b);
  h1FitWindows->GetListOfFunctions()->Add(g2b);
  h1FitWindows->GetListOfFunctions()->Add(g3b);
  h1FitWindows->GetListOfFunctions()->Add(g4b);
  h1FitWindows->GetListOfFunctions()->Add(leg);

  ///////////////////////////////////////////////////////////////////////////////////

  // Initial fits
  TF1 * g1bClone = (TF1*)g1b->Clone();
  TF1 * g2bClone = (TF1*)g2b->Clone();
  TF1 * g3bClone = (TF1*)g3b->Clone();
  TF1 * g4bClone = (TF1*)g4b->Clone();
  //   TF1 *g1bClone = new TF1(); g1b->Copy(*g1bClone);
  //   TF1 *g2bClone = new TF1(); g2b->Copy(*g2bClone);
  //   TF1 *g3bClone = new TF1(); g3b->Copy(*g3bClone);
  //   TF1 *g4bClone = new TF1(); g4b->Copy(*g4bClone);


  g1bClone->SetRange(dEdxMin,dEdxMax);
  g2bClone->SetRange(dEdxMin,dEdxMax);
  g3bClone->SetRange(dEdxMin,dEdxMax);
  g4bClone->SetRange(dEdxMin,dEdxMax);


  legClone->AddEntry(h1InitialFits,"Data","LEP");
  legClone->AddEntry(g1bClone,"e","L");
  legClone->AddEntry(g2bClone,"#pi","L");
  legClone->AddEntry(g3bClone,"K","L");
  legClone->AddEntry(g4bClone,"p","L");


  h1InitialFits->GetListOfFunctions()->Add(g1bClone);
  h1InitialFits->GetListOfFunctions()->Add(g2bClone);
  h1InitialFits->GetListOfFunctions()->Add(g3bClone);
  h1InitialFits->GetListOfFunctions()->Add(g4bClone);
  h1InitialFits->GetListOfFunctions()->Add(legClone);

  ///////////////////////////////////////////////////////////////////////////////////

  debugFile->GetFile()->cd();
  h1FitWindows  -> Write(Form("Windows_%d_%4.3f",islice,pt));
  h1InitialFits -> Write(Form("InitialFits_%d_%4.3f",islice,pt));

  delete g1bClone;
  delete g2bClone;
  delete g3bClone;
  delete g4bClone;
  delete h1FitWindows;
  delete h1InitialFits;


}
// -------------------------------------------------------------------------------------------------------
void SetStyleFor1DSlice(TH1D *h)
{

  //
  // 1D slice plots should have the same style
  //

  h->SetMarkerStyle(24);
  h->SetLineColor(kBlack);
  h->GetXaxis()->SetTitle("TPC d#it{E}/d#it{x} Signal (a.u.)");
  h->GetYaxis()->SetTitle("entries");
  h->SetMarkerSize(0.7);

}
// -------------------------------------------------------------------------------------------------------
void ReadHistograms()
{

  //
  //  Read expected Clean and h2D form the histsFile
  //

  // Copy files for backup
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/PIDIterativeFitting.C .");
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllHistsOnFarm.C .");
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllParamPlots.C .");
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C .");

  // Cache histograms in memory
  hCleanSamples   = new TH2D *[nParticle];
  hExpectedMean   = new TH2D *[nParticle];
  hExpectedSigma  = new TH2D *[nParticle];
  h2D = new TH2D();

  for (Int_t icache=0; icache<nParticle; icache++)
  {
    hCleanSamples[icache] = NULL;
    hExpectedMean[icache] = NULL;
    hExpectedSigma[icache]= NULL;
  }

  // Take main 2D hists
  cout << " take Files " << endl;

  if (fSign==0)       h2D = (TH2D *)inputFile->Get("h2Dall");
  else if (fSign==1 ) h2D = (TH2D *)inputFile->Get("h2Dpos");
  else if (fSign==-1) h2D = (TH2D *)inputFile->Get("h2Dneg");

  if (MCclosure){
    if (fSign==1){
      hCleanSamples[0]  = (TH2D *)inputFile->Get("h2CleanElectronPos");
      hCleanSamples[1]  = (TH2D *)inputFile->Get("h2CleanPionPos");
      hCleanSamples[2]  = (TH2D *)inputFile->Get("h2CleanKaonPos");
      hCleanSamples[3]  = (TH2D *)inputFile->Get("h2CleanProtonPos");
    }
    if (fSign==-1){
      hCleanSamples[0]  = (TH2D *)inputFile->Get("h2CleanElectronNeg");
      hCleanSamples[1]  = (TH2D *)inputFile->Get("h2CleanPionNeg");
      hCleanSamples[2]  = (TH2D *)inputFile->Get("h2CleanKaonNeg");
      hCleanSamples[3]  = (TH2D *)inputFile->Get("h2CleanProtonNeg");
    }
    if (fSign==0){
      hCleanSamples[0]  = (TH2D *)inputFile->Get("h2CleanElectron");
      hCleanSamples[1]  = (TH2D *)inputFile->Get("h2CleanPion");
      hCleanSamples[2]  = (TH2D *)inputFile->Get("h2CleanKaon");
      hCleanSamples[3]  = (TH2D *)inputFile->Get("h2CleanProton");
    }
  } else {
     hCleanSamples[0]  = (TH2D *)inputFile->Get("h2CleanElectron");
     hCleanSamples[1]  = (TH2D *)inputFile->Get("h2CleanPion");
     hCleanSamples[2]  = (TH2D *)inputFile->Get("h2CleanKaon");
     hCleanSamples[3]  = (TH2D *)inputFile->Get("h2CleanProton");
  }

  hExpectedMean[0]  = (TH2D *)inputFile->Get("h2ExpectedEl");
  hExpectedMean[1]  = (TH2D *)inputFile->Get("h2ExpectedPi");
  hExpectedMean[2]  = (TH2D *)inputFile->Get("h2ExpectedKa");
  hExpectedMean[3]  = (TH2D *)inputFile->Get("h2ExpectedPr");
  hExpectedSigma[0] = (TH2D *)inputFile->Get("h2ExpectedSigmaEl");
  hExpectedSigma[1] = (TH2D *)inputFile->Get("h2ExpectedSigmaPi");
  hExpectedSigma[2] = (TH2D *)inputFile->Get("h2ExpectedSigmaKa");
  hExpectedSigma[3] = (TH2D *)inputFile->Get("h2ExpectedSigmaPr");

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cout << " =========== Write histograms into the debug file ============ " << endl;
  debugFile -> GetFile()->cd();
  h2D               -> Write("h2Dall");
  hCleanSamples[0]  -> Write("h2Electron");
  hCleanSamples[1]  -> Write("h2Pion");
  hCleanSamples[2]  -> Write("h2Kaon");
  hCleanSamples[3]  -> Write("h2Proton");
  hExpectedMean[0]  -> Write("h2ExpectedEl");
  hExpectedMean[1]  -> Write("h2ExpectedPi");
  hExpectedMean[2]  -> Write("h2ExpectedKa");
  hExpectedMean[3]  -> Write("h2ExpectedPr");
  hExpectedSigma[0] -> Write("h2ExpectedSigmaEl");
  hExpectedSigma[1] -> Write("h2ExpectedSigmaPi");
  hExpectedSigma[2] -> Write("h2ExpectedSigmaKa");
  hExpectedSigma[3] -> Write("h2ExpectedSigmaPr");

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // for proper error calculations
  cout << " =========== Sumw2 ============ " << endl;
  h2D               -> Sumw2();
  hCleanSamples[0]  -> Sumw2();
  hCleanSamples[1]  -> Sumw2();
  hCleanSamples[2]  -> Sumw2();
  hCleanSamples[3]  -> Sumw2();

}
// -------------------------------------------------------------------------------------------------------
void GetKurtosisFromFile()
{

  //
  //  Read Kurtosis values from the reference fits
  //
  cout << " ============================================================================" << endl;
  cout << " ================== Kurtosis is being taken from the table ==================" << endl;
  cout << " ============================================================================" << endl;
  cout << " ============================ Current info is ===============================" << endl;
  cout << " sign     = " << fSign     << endl;
  cout << " etaDown  = " << fEtaDown  << endl;
  cout << " centDown = " << fCentDown << endl;
  cout << " elKurt   = " << kurtosisFix[kElectron] << endl;
  cout << " piKurt   = " << kurtosisFix[kPion]     << endl;
  cout << " kaKurt   = " << kurtosisFix[kKaon]     << endl;
  cout << " prKurt   = " << kurtosisFix[kProton]   << endl;
  cout << " ============================================================================" << endl;
  cout << " ============================== New info is =================================" << endl;

  // Read kurtosis tabel tree
  TFile fKurtosisTable(kurtosisTableFile);
  TTree * treeKurtosisTable = (TTree*)fKurtosisTable.Get("IdMethodInput");

  Int_t   signTable        = 0;
  Float_t etaTable         = 0.;
  Float_t centTable        = 0.;
  Double_t elKurtosisTable = 0.;
  Double_t piKurtosisTable = 0.;
  Double_t kaKurtosisTable = 0.;
  Double_t prKurtosisTable = 0.;

  // kinematic variables
  treeKurtosisTable->SetBranchAddress("fSign"      ,&signTable);
  treeKurtosisTable->SetBranchAddress("eta"        ,&etaTable);
  treeKurtosisTable->SetBranchAddress("cent"       ,&centTable);
  treeKurtosisTable->SetBranchAddress("elKurtosis" ,&elKurtosisTable);
  treeKurtosisTable->SetBranchAddress("piKurtosis" ,&piKurtosisTable);
  treeKurtosisTable->SetBranchAddress("kaKurtosis" ,&kaKurtosisTable);
  treeKurtosisTable->SetBranchAddress("prKurtosis" ,&prKurtosisTable);
  Int_t N = treeKurtosisTable->GetEntries();
  for (Int_t i=0; i<N; ++i) {
    treeKurtosisTable->GetEntry(i);
    if (fSign == signTable && TMath::Abs(fEtaDown-etaTable)<0.01 && TMath::Abs(fCentDown-centTable)<0.01) {
      kurtosisFix[kElectron] = elKurtosisTable;
      kurtosisFix[kPion]     = piKurtosisTable;
      kurtosisFix[kKaon]     = kaKurtosisTable;
      kurtosisFix[kProton]   = prKurtosisTable;
      cout << " sign     = " << signTable     << endl;
      cout << " etaDown  = " << etaTable  << endl;
      cout << " centDown = " << centTable << endl;
      cout << " elKurt   = " << kurtosisFix[kElectron] << endl;
      cout << " piKurt   = " << kurtosisFix[kPion]     << endl;
      cout << " kaKurt   = " << kurtosisFix[kKaon]     << endl;
      cout << " prKurt   = " << kurtosisFix[kProton]   << endl;
      break;
    }
  }

}
// -------------------------------------------------------------------------------------------------------
void PrepareLineShapes(TH1D *h,TF1 *ftot,TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4, Int_t centbin, Int_t etabin, Int_t mombin, Int_t signbin, Int_t iIter){

    //
    // Dump histogram and corresponding functions into lookup table
    //
    //
    // clone of the original histogram
    TH1D *htot = (TH1D*)h->Clone();
    htot->SetMarkerStyle(24);
    htot->SetLineColor(kBlack);
    htot->GetXaxis()->SetTitle("TPC d#it{E}/d#it{x} Signal (a.u.)");
    htot->GetYaxis()->SetTitle("entries");
    htot->SetMarkerSize(0.7);
    htot->SetName(Form("h_Total_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin));
    //
    //
    // TF1 *fitTot     = new TF1(); ftot->Copy(*fitTot);     fitTot    ->SetNpx(nBinsInLookUpTable);
    // TF1 *fitTotCArr = new TF1(); ftot->Copy(*fitTotCArr); fitTotCArr->SetNpx(nBinsInLookUpTable);
    TF1 *fitTot     = (TF1*)ftot->Clone();  fitTot    ->SetNpx(nBinsInLookUpTable);
    TF1 *fitTotCArr = (TF1*)ftot->Clone();  fitTotCArr->SetNpx(nBinsInLookUpTable);
    fitTot    ->SetName(Form("f_Total_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin));
    fitTotCArr->SetName(Form("fC_Total_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin));
    //
    //
    g1->SetRange(20.,1020.);  g1->SetNpx(nBinsInLookUpTable);
    g2->SetRange(20.,1020.);  g2->SetNpx(nBinsInLookUpTable);
    g3->SetRange(20.,1020.);  g3->SetNpx(nBinsInLookUpTable);
    g4->SetRange(20.,1020.);  g4->SetNpx(nBinsInLookUpTable);
    //
    //
    TString elName = Form("particle_0_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin);
    TString piName = Form("particle_1_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin);
    TString kaName = Form("particle_2_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin);
    TString prName = Form("particle_3_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin);
    //
    //
    g1->SetName(Form("Oparticle_0_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin));
    g2->SetName(Form("Oparticle_1_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin));
    g3->SetName(Form("Oparticle_2_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin));
    g4->SetName(Form("Oparticle_3_bin_%d_bin_%d_bin_%d_bin_%d",centbin,etabin,mombin,signbin));
    //
    //
    funcLineShapesObjArr.AddAt(fitTotCArr,5*iIter+0);
    funcLineShapesObjArr.AddAt(g1,    5*iIter+1);
    funcLineShapesObjArr.AddAt(g2,    5*iIter+2);
    funcLineShapesObjArr.AddAt(g3,    5*iIter+3);
    funcLineShapesObjArr.AddAt(g4,    5*iIter+4);
    //
    //
    TH1D *h1 = (TH1D*)g1->GetHistogram(); h1->SetName(elName);
    TH1D *h2 = (TH1D*)g2->GetHistogram(); h2->SetName(piName);
    TH1D *h3 = (TH1D*)g3->GetHistogram(); h3->SetName(kaName);
    TH1D *h4 = (TH1D*)g4->GetHistogram(); h4->SetName(prName);
    //
    //
    histLineShapesCArr[5*iIter+0] = (TH1D*)htot;
    histLineShapesCArr[5*iIter+1] = (TH1D*)h1;
    histLineShapesCArr[5*iIter+2] = (TH1D*)h2;
    histLineShapesCArr[5*iIter+3] = (TH1D*)h3;
    histLineShapesCArr[5*iIter+4] = (TH1D*)h4;
    //
    //
    // !!!!!!!! this part is buggy 5 dimensional gauss is not properly streamed to TClonesArray
    TF1 *t1 = new TF1(elName,fitFunctionGenGaus,dEdxMin,dEdxMax);
    TF1 *t2 = new TF1(piName,fitFunctionGenGaus,dEdxMin,dEdxMax);
    TF1 *t3 = new TF1(kaName,fitFunctionGenGaus,dEdxMin,dEdxMax);
    TF1 *t4 = new TF1(prName,fitFunctionGenGaus,dEdxMin,dEdxMax);
    Double_t *par1 = g1->GetParameters();
    Double_t *par2 = g2->GetParameters();
    Double_t *par3 = g3->GetParameters();
    Double_t *par4 = g4->GetParameters();
    t1->SetParameters(par1[0],par1[1],par1[2],par1[3],par1[4]);
    t2->SetParameters(par2[0],par2[1],par2[2],par2[3],par2[4]);
    t3->SetParameters(par3[0],par3[1],par3[2],par3[3],par3[4]);
    t4->SetParameters(par4[0],par4[1],par4[2],par4[3],par4[4]);
    funcLineShapesCArr[5*iIter+0] = (TF1*)fitTot;
    funcLineShapesCArr[5*iIter+1] = (TF1*)t1;
    funcLineShapesCArr[5*iIter+2] = (TF1*)t2;
    funcLineShapesCArr[5*iIter+3] = (TF1*)t3;
    funcLineShapesCArr[5*iIter+4] = (TF1*)t4;

    /*

    TFile file("LineShapes_0_1_100_0.20_2.20_-0.10_30.00.root")
    TH1D *h = (TH1D*)histLineShapesCArr->FindObject("particle_2_bin_5_bin_5_bin_2_bin_2")
    TF1  *fC = (TF1*)funcLineShapesCArr->FindObject("Cparticle_2_bin_5_bin_5_bin_2_bin_2")
    TF1  *f = (TF1*)funcLineShapesObjArr->FindObject("particle_2_bin_5_bin_5_bin_2_bin_2")

    cout << fC->GetParameter(0) << "  " << f->GetParameter(0) << endl;
    cout << fC->GetParameter(1) << "  " << f->GetParameter(1) << endl;
    cout << fC->GetParameter(2) << "  " << f->GetParameter(2) << endl;
    cout << fC->GetParameter(3) << "  " << f->GetParameter(3) << endl;
    cout << fC->GetParameter(4) << "  " << f->GetParameter(4) << endl;
    cout << " -------------------------------------- " << endl;

    TH1D *k = histLineShapesCArr->FindObject("h_Total_bin_5_bin_5_bin_2_bin_2")
    k->Draw();
    f->Draw("same")
    h->Draw("same")


    h->GetMaximum()
    f->GetMaximum()
    xmin=0
    xmax=100
    TAxis *axis = h->GetXaxis();
    int bmin = axis->FindBin(xmin); //in your case xmin=-1.5
    int bmax = axis->FindBin(xmax); //in your case xmax=0.8
    double integral = h->Integral(bmin,bmax);
    integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/
    axis->GetBinWidth(bmin);
    integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax))/
    axis->GetBinWidth(bmax);
    integral/f->Integral(0,100)

    */

}
