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

#include "ROOT/RDataFrame.hxx"

#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::cout;
using std::setw;
using ROOT::RDataFrame;


// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************* Modification region Start ********************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************
// ******************************************************************************************************************************

// dumping bools
Int_t fVariableSkewnes        = 0;      // 0 for all particles have variable skew, 1: only pions have variable skew
Double_t applySmoothAfter     = 0.5;
Bool_t normalisedCleanSamples = kFALSE;
Bool_t dump                   = kTRUE;  // dump additional plots, not necessary when running for all fits
Bool_t exWRTMean              = kTRUE;   // find
// Initial Settings
Bool_t KSestimation           = kTRUE;   // make KS estimation from samples --> if false use gauss fits
Bool_t AutoK_FixedS           = kTRUE;    // if true only kurtosis is being estimated, output skewness values should be absolute
// otherwise should be <input>*estimatedSkewness
Int_t AutoK_FixedS_Ksetting   = 0;        // kurtosis setting --> 0; pi,K,pr auto, 1; pi auto, K,pr=piKurtosis, 2; pi auto, K,pr = 2
Bool_t fixedK                 = kTRUE;   // fix to auto found (KSestimation=kTRUE) or gauss pars (KSestimation=kFALSE) --> keep true
Bool_t fixedS                 = kTRUE;   // fix to auto found (KSestimation=kTRUE) or gauss pars (KSestimation=kFALSE) --> keep true
Bool_t automaticKS            = kFALSE;   // if true kurtosis and skewness set automatically
Bool_t takeKurtosisFromTable  = kFALSE;
Bool_t specialCareForLargeEta = kFALSE;   // apply only gauss fit for large eta --> >0.6

Int_t nBinsKSestimate         = 60;      // number of bins to be used in the 0.2-1.4 window for KS estimate ; 200 or 50
Bool_t useSafeRegion          = kTRUE;   // use either safe regions or clean samples for KS estimation
Double_t sysAmp               = 0.002;
Double_t sysMean              = 0.002;
Int_t assump                  = 1;        // choose a value 0,1: 0--> use only expecteds 1--> use clean + expected
Int_t fDeuteronsExist=1; // check if deutreon exist.

Double_t piKurtosisScan       = 1.;    // factor to be multiplied with the automatically found kurtosis
Double_t kaKurtosisScan       = 1.;    // factor to be multiplied with the automatically found kurtosis
Double_t prKurtosisScan       = 1.;    // factor to be multiplied with the automatically found kurtosis

Double_t piSkewnessScan       = 0.7;
Double_t kaSkewnessScan       = 0.5;
Double_t prSkewnessScan       = 0.5;

Double_t rebinFactor          = 2.;     // !!!!! dikkat
Int_t nBinsInLookUpTable      =1000;
TString cleanPionSmapleName   = "h2DCleanPi_KineCut";  // h2CleanPiTOF   or   h2CleanPiTight  or  h2CleanPion
Double_t fitWinMeanLowNsigma  = 3.;
Double_t fitWinMeanUpNsigma   = 2.;
//
Double_t sampleFitWindowDown = 2.5;
Double_t sampleFitWindowUp   = 3.;
Double_t electronCleanWindow = 8.;
Double_t skewborder[] = {1.,0.};

TString kurtosisTableFile     = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics_cRows_80_16EtaBin_mombin20MeV/Fits/KurtosisTable/KurtosisTable.root";


const Int_t nParticle  = 7;
Double_t skewnessFixFit[nParticle]   = {0., 0., 0., 0., 0., 0., 0.};
Double_t skewnessFixClean[nParticle] = {0., 0., 0., 0., 0., 0., 0.};

Double_t kurtosisFixFit[nParticle]   = {2., 2., 2., 2., 2., 0., 0.};
Double_t kurtosisFixClean[nParticle] = {2., 2., 2., 2., 2., 0., 0.};

Double_t kurtosisFix[nParticle]      = {2., 2., 2., 2., 2., 0., 0.};
Double_t skewnessFix[nParticle]      = {0., 0., 0., 0., 0., 0., 0.};

TString smoothFunction = "pol2";   // outlier removal fit function to be used in TLinearFitter
Double_t windowSetting = 4.;       // in window arrangment for mean "distToClosestParticle/windowSetting"
Double_t smoothSeg     = 5.;       // outlier removal range setting for the TLinearFitter (nSlice/smoothSeg)
Double_t smoothLevel   = 9.;       // Increases the iterations in smoothing (10 is too much) for TGraphSmooth
Double_t analysisRange = 5.;       // Maximum p value used in the analysis

Float_t elFitSafeMin = 0.24,   elFitSafeMax = 0.38;
Float_t piFitSafeMin = 0.40,   piFitSafeMax = 0.55;
Float_t kaFitSafeMin = 0.34,   kaFitSafeMax = 0.4;
Float_t prFitSafeMin = 0.56,   prFitSafeMax = 0.7;

Float_t elCleanSafeMin = 0.22, elCleanSafeMax = 0.4;
Float_t piCleanSafeMin = 0.7,  piCleanSafeMax = 0.8;
Float_t kaCleanSafeMin = 0.8,  kaCleanSafeMax = 1.;
Float_t prCleanSafeMin = 1.3,  prCleanSafeMax = 1.4;

// For free Kurtosis and skewness
Double_t skewMin       = 0.;
Double_t skewMax       = 1.5;
Double_t kurtosisMin   = 1.7;
Double_t kurtosisMax   = 2.2;

const Int_t   dEdxMin         = 20;
const Int_t   dEdxMax         = 1020;
const Int_t   fnEtaBins       = 16;   // MC: 8, Data:16
const Float_t fEtaRangeDown   = -0.8;
const Float_t fEtaRangeUp     = 0.8;
const Int_t   fnMomBins       = 150;
const Float_t fMomRangeDown   = 0.1;
const Float_t fMomRangeUp     = 3.1;
const Int_t   fnCentBins      = 9;
Float_t xCentBins[] = {0, 5,  10,  20, 30, 40, 50, 60, 70, 80};
Float_t xEtaBins[] = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

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
Double_t elMeanMin ,piMeanMin ,kaMeanMin ,prMeanMin, deMeanMin;
Double_t elMeanMax ,piMeanMax ,kaMeanMax ,prMeanMax, deMeanMax;
Double_t elAmpMin  ,piAmpMin  ,kaAmpMin  ,prAmpMin,  deAmpMin;
Double_t elAmpMax  ,piAmpMax  ,kaAmpMax  ,prAmpMax,  deAmpMax;
Double_t elSigmaMin,piSigmaMin,kaSigmaMin,prSigmaMin,deSigmaMin;
Double_t elSigmaMax,piSigmaMax,kaSigmaMax,prSigmaMax,deSigmaMax;

TVectorD *elExMean=0x0, *elExSigma=0x0;
TVectorD *piExMean=0x0, *piExSigma=0x0;
TVectorD *kaExMean=0x0, *kaExSigma=0x0;
TVectorD *prExMean=0x0, *prExSigma=0x0;
TVectorD *deExMean=0x0, *deExSigma=0x0;

// Defien max position in each expected particle mean bin
Double_t maxBin,maxBin0,maxBin1,maxBin2,maxBin3,maxBin4;

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
// TTreeSRedirector *lineShapesLookUp = 0;
TTreeSRedirector *sampleFile    = 0;
TGraphSmooth *gs;

//   Some objects
TObjArray    * Fit;
TObjArray    * Clean;
TObjArray    * ScaledClean;
TObjArray    * ScaledExWrtClean;
TObjArray    * ScaledEx;
TObjArray    * OutlierSmooth;
TObjArray    * Windows;

// helper graphs
TGraphErrors *grMeanCleanWrtExScaled[nParticle];
TGraphErrors *grSigmaCleanWrtExScaled[nParticle];

TGraphErrors *grFitAmpSmooth[nParticle];
TGraphErrors *grFitMeanSmooth[nParticle];
TGraphErrors *grFitSigmaSmooth[nParticle];
//
TGraphErrors *grFitAmp[nParticle];
TGraphErrors *grFitMean[nParticle];
TGraphErrors *grFitSigma[nParticle];
//
TGraphErrors *grMeanModExScaled[nParticle];
TGraphErrors *grMeanExScaled[nParticle];
TGraphErrors *grSigmaExScaled[nParticle];
//
TGraphErrors *grCleanAmp[nParticle];
TGraphErrors *grCleanMean[nParticle];
TGraphErrors *grCleanSigma[nParticle];
//
TGraphErrors *grCleanAmpScaled[nParticle];
TGraphErrors *grCleanMeanScaled[nParticle];
TGraphErrors *grCleanSigmaScaled[nParticle];
//
TGraphErrors *grFitAmpOutlierSmooth[nParticle];
TGraphErrors *grFitSigmaOutlierSmooth[nParticle];
//
TGraphErrors *grCleanMeanSmooth[nParticle];
//
TGraphErrors *grCleanMeanScaledSmooth[nParticle];
TGraphErrors *grCleanSigmaScaledSmooth[nParticle];
TGraphErrors *grCleanAmpScaledSmooth[nParticle];
//
TGraphErrors *grFreeKSAmp[nParticle];
TGraphErrors *grFreeKSMean[nParticle];
TGraphErrors *grFreeKSSigma[nParticle];
TGraphErrors *grFreeKSKurtosis[nParticle];
TGraphErrors *grFreeKSSkew[nParticle];
//
TGraphErrors *grLookUpCleanSkew[4];
TGraphErrors *grLookUpCleanKurtosis[4];
TGraphErrors *grLookUpFitSkew[4];
TGraphErrors *grLookUpFitKurtosis[4];

TF1 *fitLookUpCleanSkew[4];
TF1 *fitLookUpCleanKurtosis[4];
TF1 *fitLookUpFitSkew[4];
TF1 *fitLookUpFitKurtosis[4];

TH2D ** hCleanSamples=NULL;
TH2D ** hExpectedMean=NULL;
TH2D ** hExpectedSigma=NULL;
TH2D  * h2D=NULL;
TFile * inputFile=NULL;
TFile * splinesFile=NULL;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TString fitFunctionGenGaus = "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
// TString fitFunctionGenGaus = "[0]*exp(-(abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/1.414213))";
Float_t fEtaDown, fCentDown;
Float_t fEtaUp, fCentUp;
Int_t particleBin,centBin,etaBin,momBin,signBin;
TString centEtaStr = "";
TString corrStr = "";
Int_t dEdxCorr = -1;


enum ParticleType {
  kEl = 0,
  kPi = 1,
  kKa = 2,
  kPr = 3,
  kDe = 4,
  kKaTOF = 5,
  kPrTOF = 6,
} pType;
TString particleStr[]={"el", "pi", "ka", "pr", "de", "kaTOF", "prTOF"};

TString expectedsFile = "/Volumes/Extern/data/expectedTree.root";

// Main Functions for iterative fitting procedure
void          SmoothAmplitudes(TString fileFit, TString fileSample, TString fileOut, const Int_t nSlice, Int_t iIter);
void          GetExandSampleParameters(Int_t islice, TFile *sFile, Double_t *arrMean, Double_t *arrMeanModif, Double_t *arrSigma, Double_t *cleanParams);
void          GetHelperGraphs(Int_t iIter, TFile *ressFile);
void          GetCleanExParams(TH1D *hClean, Double_t *arrSigma, Double_t *arrCleanSigma, Double_t *arrCleanMean, Double_t *arrCleanMax, ParticleType pSpecy);

void          EstimateKS(Int_t iks, TH2D *h2DKS);
void          FitAllSamples(const Int_t nSlice, const Double_t ptMin, const Double_t ptMax, TString sampleFileName);
void          IterativeFitting(Int_t iIter, const Int_t nSlice, const Double_t ptMin, const Double_t ptMax, TString fileIn, TString readSmooth, TString readSamp,TString lineShapesFile);
void          ProduceGraphs(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grCleanTmp, TGraph **grCleanSmoothTmp);
void          ApplyOutlierSmooth(const Int_t nSlice, TTree *t, TObjArray *arr, TObjArray *barr, TGraphErrors **grTmp, TGraph **grSmoothTmp);
void          ApplyInitialFits(Int_t islice, Double_t pt, TH1D *h1D, Double_t *arrSigma, Double_t *arrMean);
void          SetStyleFor1DSlice(TH1D *h);
void          ReadHistograms();
void          GetKurtosisFromFile();


// Helper functions for severel computations
void          Fit1DSliceForKSEstimate(Int_t iks, TH1D *h1DKS, Double_t *pars, Double_t *arrMean, Double_t *arrSigma);
void          FitParticleSampleForKSEstimate(Int_t iks, TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy);
TMatrixD     *FitParticleSample(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy, Double_t pt);
void          FitParticleSampleFreeKS(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy, Double_t pt);
TF1          *FitParamInRange(TGraphErrors* gr, TString funcType, Double_t min, Double_t max);
TH1D         *FuncToHist(TF1 *f, TH1D * h, Int_t islice, Double_t *arrMean, Double_t *arrSigma);
TH1D         *GraphShade(TGraphErrors *grFitParam, Double_t percent);
TGraphErrors *RemoveOutliers(TGraphErrors *grsmooth, Int_t fitWindow);
TGraphErrors *ApplyScalingWrtExpected(const Int_t nSlice, TString exp, TString fit, const Double_t safeP, TObjArray * arr, TTree *t);
TGraphErrors *ApplyScalingWrtSample(const Int_t nSlice, TString sample, TString fit, const Double_t safeP, TObjArray * arr, TObjArray * arrSample, TTree *tSample, TGraph **grScaledSmooth,Int_t k);
TGraphErrors *ProduceWindowGraphs(Int_t nSlice,Int_t index, TString str, TTree *tree);

void          ModifyExpectedMeanAndApplyVariableSkewness(Double_t *arrMean, Double_t *arrSigma, Double_t *arrMeanModif, Double_t pt);

Double_t      ComputeGaussIntegral(Double_t parAmp,Double_t parMean,Double_t parSigma,Double_t parKurtosis, Double_t parSkew);
TH1D         *GetClean1DSlice(TH2D *h2Clean, TString parName, Int_t ptDownBin);
Int_t         GetMaxIndex(TGraphErrors *grIndex);
Int_t         GetPtIndex(TGraphErrors *grIndex,Double_t p);
Double_t      GetClosestParticleMean(Double_t sig, Double_t ref, Double_t gr1, Double_t gr2, Double_t gr3);
TGraphErrors *GraphSmooth(TGraphErrors * gr);
TGraphErrors *HistToGraphErrors(TH1D * h);
TCanvas      *GetFitResCanvas(TObjArray *arr, TObjArray *arrclean, TObjArray *arrscaledclean, TObjArray *arrscaledex);
TH2D         *TH2Subset(TH2D *h,Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2);
void          MakeAnimation(Int_t iIter, TH1D* h1ToGif, TH2D* h2Danimation, Double_t pt, Double_t ptMin, Double_t ptMax);

void          CreateSplinesFromExpectedsTree(TString expectedsFile);


// Set functions for the parameters
void          SetTotalFitParameters(Int_t iIter, Int_t islice, Int_t nSlice, TF1* total, Double_t *arrMean, Double_t *arrSigma, Double_t *arrMeanWindow, Double_t *cleanParams, TH1D* h1D, Double_t pt);
void          CheckIfMeanIsZero(TF1 *total, Double_t *arrMean);
void          SetParNamesOfTotalFit(TF1 *total);
void          SetParamsForKSEstimate(Int_t iks, TF1 *total, TH1D *h1DKS, Double_t *arrMean, Double_t *arrSigma);
void          SetFixedKSparameterForTotalFit(TF1 *total, Double_t pt);
void          SetFixedKSparameterForTotalFitForPions(TF1 *total);
void          SetFixedKSparameterForCleanSampleFit(TF1 *asyGaus, ParticleType pSpecy, Double_t pt);
void          SetFixedKSparameterForIndividualFits(TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4, TF1 *g5);
void          PrepareLineShapes(TH1D *h,TF1 *ftot,TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4, TF1 *g5, Int_t centbin, Int_t etabin, Int_t mombin, Int_t signbin, Int_t iIter);
void          ResetParametersOfEachFunction(TF1 *g1,TF1 *g2,TF1 *g3,TF1 *g4,TF1 *g5,TF1 *total,TH1D *h1D,const Double_t pt);
void          SetParams1stIteration(Double_t pt, TH1D *h1D, Double_t *arrMean, Double_t *arrSigma,Double_t *arrMeanWindow, Double_t *cleanParams);


// -------------------------------------------------------------------------------------------------------
void RealData_PIDIterativeFitting_RUN2(TString corr="", Int_t sign=0, TString fileData="", TString splinesData="", const Int_t nSlice=3, Double_t ptMin=0.2, Double_t ptMax=2.2, Double_t etaMin=0., Double_t etaMax=0.1,Double_t centMin=0., Double_t centMax=5., Int_t maxIter=8, Double_t piSkewnessScanIn=0.7, Double_t kaSkewnessScanIn=0.5, Double_t prSkewnessScanIn=0.5, Double_t piKurtosisScanIn=1., Double_t kaKurtosisScanIn=1., Double_t prKurtosisScanIn=1.)
{

  //
  // Apply gaussian fits to the ptot slices --> Main Function
  /*

  meld /home/marsland/Desktop/ubuntu_desktop/workdir/gsi_marsland_EbyeRatios/RealData_PIDIterativeFitting_RUN2.C /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_PIDIterativeFitting_RUN2.C

  cd  /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN2/LHC15o_pass1/LHC15o_pass1_NoSelection_20082019/Fits/test
  aliroot -l -b
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_PIDIterativeFitting_RUN2.C+
  RealData_PIDIterativeFitting_RUN2("", 0,"rootFiles/Hists.root", "rootFiles/Splines_Syst0_PlotVS_ptot_CutON_p.root", 100,0.2,2.2,  0.2,0.3,  0,5,     3,     0.7,0.5,0.5,1,1,1);  2> err.log

  TFile f("Samples_0_100_0.20_2.20_0.20_0.root")
  CleanSamples->SetMarkerStyle(20);
  CleanSamples->Draw("deMean:p","","")
  CleanSamples->Draw("kaMean:p","","same")
  CleanSamples->Draw("elMean:p","","same")
  CleanSamples->Draw("piMean:p","","same")
  CleanSamples->Draw("prMean:p","","same")
  CleanSamples->SetMarkerColor(kRed);
  CleanSamples->Draw("deCleanMean:p","","same")
  CleanSamples->Draw("kaCleanMean:p","","same")
  CleanSamples->Draw("elCleanMean:p","","same")
  CleanSamples->Draw("piCleanMean:p","","same")
  CleanSamples->Draw("prCleanMean:p","","same")


  aliroot -l -b
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_PIDIterativeFitting_RUN2.C+
  TString splines = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN2/LHC15o_pass1/LHC15o_pass1_NoSelection_20082019//Fits/test/rootFiles/Splines_Syst0_PlotVS_ptot_CutON_p.root"

  RealData_PIDIterativeFitting_RUN2("Corr", 0,"../HistsCorr.root", splines, 100,0.2,2.2,  0.2,0.3,  0,5,     3,     0.7,0.5,0.5,1,1,1);  2> err.log




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

  if (!KSestimation && !automaticKS) {

    skewnessFixFit[kEl] = 0.;                      kurtosisFixFit[kEl] = 2.;
    skewnessFixFit[kPi] = piSkewnessScanIn;        kurtosisFixFit[kPi] = piKurtosisScanIn*2.;
    skewnessFixFit[kKa] = kaSkewnessScanIn;        kurtosisFixFit[kKa] = kaKurtosisScanIn*2.;
    skewnessFixFit[kPr] = prSkewnessScanIn;        kurtosisFixFit[kPr] = prKurtosisScanIn*2.;
    skewnessFixFit[kDe] = 0.;                      kurtosisFixFit[kDe] = 2.;

    skewnessFixClean[kEl] = 0.;                    kurtosisFixClean[kEl] = 2.;
    skewnessFixClean[kPi] = piSkewnessScanIn;      kurtosisFixClean[kPi] = piKurtosisScanIn*2.;
    skewnessFixClean[kKa] = kaSkewnessScanIn;      kurtosisFixClean[kKa] = kaKurtosisScanIn*2.;
    skewnessFixClean[kPr] = prSkewnessScanIn;      kurtosisFixClean[kPr] = prKurtosisScanIn*2.;
    skewnessFixClean[kDe] = 0.;                    kurtosisFixClean[kDe] = 2.;

    skewnessFix[kEl] = 0.;                         kurtosisFix[kEl] = 2.;
    skewnessFix[kPi] = piSkewnessScanIn;           kurtosisFix[kPi] = piKurtosisScanIn*2.;
    skewnessFix[kKa] = kaSkewnessScanIn;           kurtosisFix[kKa] = kaKurtosisScanIn*2.;
    skewnessFix[kPr] = prSkewnessScanIn;           kurtosisFix[kPr] = prKurtosisScanIn*2.;
    skewnessFix[kDe] = 0.;                         kurtosisFix[kDe] = 2.;

  }
  //
  // Get eta and cent info
  fEtaDown     = etaMin;
  fCentDown    = centMin;
  fEtaUp       = etaMax;
  fCentUp      = centMax;
  Double_t etaBinCenter = (etaMin + etaMax)/2.;
  centEtaStr = Form("cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f",centMin,centMax,etaMin,etaMax);
  corrStr = corr;
  //
  //
  etaBin    = fhEta -> FindBin(fEtaDown - 0.0000001)  + 1;
  centBin   = fhCent-> FindBin(fCentDown - 0.0000001) + 1;
  signBin   = sign+1;

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "    << endl;
  cout << " file Path      = " << fileData                << endl;
  cout << " eta            = " << fEtaDown                << endl;
  cout << " cent           = " << fCentDown               << endl;
  cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "    << endl;
  cout << " KSestimation   = " << KSestimation            << endl;
  cout << " AutoK_FixedS   = " << AutoK_FixedS            << endl;
  cout << " automaticKS    = " << automaticKS             << endl;
  cout << " CareLargeEta   = " << specialCareForLargeEta  << endl;
  cout << " useSafeRegion  = " << useSafeRegion           << endl;
  cout << " nBinsKSestimate= " << nBinsKSestimate         << endl;
  cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "    << endl;
  cout << " fixedK         = " << fixedK            << endl;
  cout << " fixedS         = " << fixedS            << endl;
  cout << " piSkewness     = " << piSkewnessScan    << endl;
  cout << " kaSkewness     = " << kaSkewnessScan    << endl;
  cout << " prSkewness     = " << prSkewnessScan    << endl;
  cout << " piKurtScan     = " << piKurtosisScan    << endl;
  cout << " kaKurtScan     = " << kaKurtosisScan    << endl;
  cout << " prKurtScan     = " << prKurtosisScan    << endl;
  cout << " ================================== "    << endl;
  cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "    << endl;

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // Create an debug file
  debugFile = new TTreeSRedirector(Form("debugFile_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown),"recreate");
  inputFile = TFile::Open(fileData);
  splinesFile = TFile::Open(splinesData);
  ReadHistograms();
  //
  //  Estimate kurtosis and Skewness for kaon, proton and pion by fitting 100 slices in between [0.4,1.4]GeV/c
  cout << "Get The Clone of h2D " << h2D->GetNbinsX() << "   " << h2D->GetName() << endl;
  TH2D *h2DKS = (TH2D*)h2D->Clone();
  if (KSestimation){
    //
    for (Int_t iks = 0; iks<2; iks++){
      cout << " =========== KS estimate  ============  Iteration = " << iks << endl;
      EstimateKS(iks,h2DKS);
    }
    //
    // semi auto KS estimation skewness is given from outside
    if (AutoK_FixedS) {

      // set skewwness
      skewnessFix[kPi] = piSkewnessScanIn;
      skewnessFix[kKa] = kaSkewnessScanIn;
      skewnessFix[kPr] = prSkewnessScanIn;
      skewnessFix[kDe] = 0.;
      // set kurtosis
      if (AutoK_FixedS_Ksetting==0){
        cout << " piKurtosisScanIn = " << piKurtosisScanIn <<"    " << kurtosisFixFit[kPi] << endl;
        cout << " kaKurtosisScanIn = " << kaKurtosisScanIn <<"    " << kurtosisFixClean[kKa]<< endl;
        cout << " prKurtosisScanIn = " << prKurtosisScanIn <<"    " << kurtosisFixFit[kPr]<< endl;
        //
        skewnessFix[kEl] = piSkewnessScanIn;
        //
        kurtosisFix[kEl] = kurtosisFixClean[kEl];
        kurtosisFix[kPi] = kurtosisFixFit[kPi]*piKurtosisScanIn;
        kurtosisFix[kKa] = kurtosisFixFit[kKa]*kaKurtosisScanIn;
        kurtosisFix[kPr] = kurtosisFixFit[kPr]*prKurtosisScanIn;
        kurtosisFix[kDe] = 2.;
      }
      if (AutoK_FixedS_Ksetting==1){
        kurtosisFix[kPi] = kurtosisFixFit[kPi]*piKurtosisScanIn;
        kurtosisFix[kKa] = kurtosisFix[kPi];
        kurtosisFix[kPr] = kurtosisFix[kPi];
        kurtosisFix[kEl] = 2.;
        kurtosisFix[kDe] = 2.;
      }
      if (AutoK_FixedS_Ksetting==2){
        kurtosisFix[kPi] = kurtosisFixFit[kPi]*piKurtosisScanIn;
        kurtosisFix[kKa] = 2.;
        kurtosisFix[kPr] = 2.;
        kurtosisFix[kEl] = 2.;
        kurtosisFix[kDe] = 2.;
      }
    }

  }
  //
  // First fit the samples
  cout << " =========== first make sample file and clean them up =========== " << endl;
  TString fileSample  = Form("Samples_%d_%d_%3.2f_%3.2f_%3.2f_%d.root",fSign,nSlice,ptMin,ptMax,fEtaDown,0);      // care only eta dependence
  FitAllSamples(nSlice,ptMin,ptMax,fileSample);
  //
  cout << " ================================== "    << endl;
  cout << "skewnessFixElectron = " <<  skewnessFix[kEl] << endl;
  cout << "skewnessFixPion     = " <<  skewnessFix[kPi] << endl;
  cout << "skewnessFixKaon     = " <<  skewnessFix[kKa] << endl;
  cout << "skewnessFixProton   = " <<  skewnessFix[kPr] << endl;
  cout << "skewnessFixDeuteron = " <<  skewnessFix[kDe] << endl;
  cout << " ================================== "    << endl;
  cout << "kurtosisFixElectron = " <<  kurtosisFix[kEl] << endl;
  cout << "kurtosisFixPion     = " <<  kurtosisFix[kPi] << endl;
  cout << "kurtosisFixKaon     = " <<  kurtosisFix[kKa] << endl;
  cout << "kurtosisFixProton   = " <<  kurtosisFix[kPr] << endl;
  cout << "kurtosisFixDeuteron = " <<  kurtosisFix[kDe] << endl;
  cout << " ================================== "    << endl;
  //
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
  // //
  // // for large eta bin do not apply auto ks estimation ???
  // if ( specialCareForLargeEta && (TMath::Abs(etaBinCenter)>=0.6) && KSestimation )
  // {
  //   for (Int_t i=0;i<nParticle;i++)
  //   {
  //     skewnessFix[i] = 0.; kurtosisFix[i] = 2.;
  //   }
  // }
  //
  // Check the fit params just in case kurt and skewness 0
  for (Int_t i=0;i<nParticle;i++) {
    if ((skewnessFix[i]<=-0.1) || (kurtosisFix[i]<1.))
    {
      skewnessFix[i] = 0.;
      kurtosisFix[i] = 2.;
    }
  }
  //
  // iteratively apply fits
  TStopwatch timer; timer.Start();
  for (Int_t iIter = 0; iIter<10; iIter++){
    if(iIter>maxIter) break;
    TString fileFit    = Form("Trees_Iter%d_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",iIter,fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown);
    TString lineShapesFile = Form("LineShapes_%d_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",iIter,fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown);
    TString fileOut    = Form("Results_Iter%d_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",iIter,fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown);
    Int_t previousIter = (iIter==6) ? 1 : iIter-1;
    //     Int_t previousIter = iIter-1;
    TString readSmooth = Form("Results_Iter%d_%d_%d_%3.2f_%3.2f_%3.2f_%3.2f.root",previousIter,fSign,nSlice,ptMin,ptMax,fEtaDown,fCentDown);
    cout << " =========== Iterative fit procedure is being started  ============  Iteration = " << iIter << endl;
    IterativeFitting(iIter,nSlice,ptMin,ptMax,fileFit,readSmooth,fileSample,lineShapesFile);
    cout << " =========== Iteration = " << iIter << "   has finished go to iteration = " << iIter+1 << endl;
    if (nSlice>50) SmoothAmplitudes(fileFit,fileSample,fileOut,nSlice,iIter);  // ?????
  }
  timer.Stop(); timer.Print();

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
  // lineShapesLookUp = new TTreeSRedirector(lineShapesFile,"recreate");

  // OutPut Arrays
  fitResArr    . SetOwner(kTRUE);
  fitResiduals . SetOwner(kTRUE);
  histLineShapesCArr.SetOwner(kTRUE);
  funcLineShapesCArr.SetOwner(kTRUE);
  funcLineShapesObjArr.SetOwner(kTRUE);

  // loop over slices
  Double_t elIntegralFit[allBins], elFitAmp[allBins], elFitMean[allBins], elFitSigma[allBins], elFitMax[allBins];
  Double_t piIntegralFit[allBins], piFitAmp[allBins], piFitMean[allBins], piFitSigma[allBins], piFitMax[allBins];
  Double_t kaIntegralFit[allBins], kaFitAmp[allBins], kaFitMean[allBins], kaFitSigma[allBins], kaFitMax[allBins];
  Double_t prIntegralFit[allBins], prFitAmp[allBins], prFitMean[allBins], prFitSigma[allBins], prFitMax[allBins];
  Double_t deIntegralFit[allBins], deFitAmp[allBins], deFitMean[allBins], deFitSigma[allBins], deFitMax[allBins];

  Double_t elFitKurtosis[allBins], elFitSkew[allBins], elChi2Fit[allBins];
  Double_t piFitKurtosis[allBins], piFitSkew[allBins], piChi2Fit[allBins];
  Double_t kaFitKurtosis[allBins], kaFitSkew[allBins], kaChi2Fit[allBins];
  Double_t prFitKurtosis[allBins], prFitSkew[allBins], prChi2Fit[allBins];
  Double_t deFitKurtosis[allBins], deFitSkew[allBins], deChi2Fit[allBins];

  Int_t ptDownBin = 0;
  Int_t ptUpBin   = 0;
  Float_t ptStep = 0.;
  Float_t ptDown = 0.;
  Float_t ptUp   = 0.;
  Float_t pt     = 0.;
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
    //
    //
    Double_t * arrSigma      = new Double_t[10];
    Double_t * arrMean       = new Double_t[10];
    Double_t * arrMeanModif  = new Double_t[10];
    Double_t * arrMeanWindow = new Double_t[10];
    Double_t * cleanParams   = new Double_t[40];

    for (Int_t ipar=0; ipar<10; ipar++)
    {
      arrMean[ipar]       = 0;
      arrMeanModif[ipar]  = 0;
      arrSigma[ipar]      = 0;
      arrMeanWindow[ipar] = 0;
    }
    for (Int_t ipar=0; ipar<40; ipar++) cleanParams[ipar]  = 0;

    //
    GetExandSampleParameters(islice,sampFile,arrMean,arrMeanModif,arrSigma,cleanParams);
    if (dump) ApplyInitialFits(islice,pt,h1D,arrSigma,arrMean);

    // Fit part
    maxBin = h1D->GetMaximum();

    elMeanMin = TMath::Max(30.,(arrMean[kEl]-2.*arrSigma[kEl]));                                   // electron
    piMeanMin = TMath::Max(30.,(arrMean[kPi]-fitWinMeanLowNsigma*arrSigma[kPi]));                  // pion
    kaMeanMin = TMath::Max(30.,(arrMean[kKa]-fitWinMeanLowNsigma*arrSigma[kKa]));                  // kaon
    prMeanMin = TMath::Max(30.,(arrMean[kPr]-fitWinMeanLowNsigma*arrSigma[kPr]));                  // proton
    deMeanMin = TMath::Max(30.,(arrMean[kDe]-fitWinMeanLowNsigma*arrSigma[kDe]));

    elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kEl]+2.*arrSigma[kEl]));
    piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kPi]+fitWinMeanUpNsigma*arrSigma[kPi]));
    kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kKa]+fitWinMeanUpNsigma*arrSigma[kKa]));
    prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kPr]+fitWinMeanUpNsigma*arrSigma[kPr]));
    deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kDe]+fitWinMeanUpNsigma*arrSigma[kDe]));

    if (pt<0.3) {
      prMeanMin = TMath::Max(30.              ,(arrMean[kPr]-arrSigma[kPr]*3.));
      prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kPr]+arrSigma[kPr]*3.));
    }

    // Parameters: 0 --> Amplitude, 1 --> Mean, 2 --> Sigma, 3 --> kurtosis, 4 --> skewness
    TF1 *g1 = new TF1("g1",fitFunctionGenGaus,elMeanMin,elMeanMax);        // g1->SetNpx(nBinsInLookUpTable);
    TF1 *g2 = new TF1("g2",fitFunctionGenGaus,piMeanMin,piMeanMax);        // g2->SetNpx(nBinsInLookUpTable);
    TF1 *g3 = new TF1("g3",fitFunctionGenGaus,kaMeanMin,kaMeanMax);        // g3->SetNpx(nBinsInLookUpTable);
    TF1 *g4 = new TF1("g4",fitFunctionGenGaus,prMeanMin,prMeanMax);        // g4->SetNpx(nBinsInLookUpTable);
    TF1 *g5 = new TF1("g5",fitFunctionGenGaus,deMeanMin,deMeanMax);        // g5->SetNpx(nBinsInLookUpTable);

    g1->SetParameters(maxBin,arrMean[kEl],arrSigma[kEl],2.,0.);
    g2->SetParameters(maxBin,arrMean[kPi],arrSigma[kPi],2.,0.);
    g3->SetParameters(maxBin,arrMean[kKa],arrSigma[kKa],2.,0.);
    g4->SetParameters(maxBin,arrMean[kPr],arrSigma[kPr],2.,0.);
    g5->SetParameters(maxBin,arrMean[kDe],arrSigma[kDe],2.,0.);

    g1->SetLineColor(4);
    g2->SetLineColor(kOrange-1);
    g3->SetLineColor(kGreen+2);
    g4->SetLineColor(kMagenta+1);
    g5->SetLineColor(kBlack);

    // restrict the parameters of gaus functions of individual particles
    g1->SetParLimits(0,0.,maxBin);
    g2->SetParLimits(0,0.,maxBin);
    g3->SetParLimits(0,0.,maxBin);
    g4->SetParLimits(0,0.,maxBin);
    g5->SetParLimits(0,0.,maxBin);

    g1->SetParLimits(1,elMeanMin,elMeanMax);
    g2->SetParLimits(1,piMeanMin,piMeanMax);
    g3->SetParLimits(1,kaMeanMin,kaMeanMax);
    g4->SetParLimits(1,prMeanMin,prMeanMax);
    g5->SetParLimits(1,deMeanMin,deMeanMax);

    g1->SetParLimits(2,arrSigma[0]/5.,arrSigma[0]*2.);
    g2->SetParLimits(2,arrSigma[1]/5.,arrSigma[1]*2.);
    g3->SetParLimits(2,arrSigma[2]/5.,arrSigma[2]*2.);
    g4->SetParLimits(2,arrSigma[3]/5.,arrSigma[3]*2.);
    g5->SetParLimits(2,arrSigma[4]/5.,arrSigma[4]*2.);

    g1->SetParLimits(3,kurtosisMin,kurtosisMax);
    g2->SetParLimits(3,kurtosisMin,kurtosisMax);
    g3->SetParLimits(3,kurtosisMin,kurtosisMax);
    g4->SetParLimits(3,kurtosisMin,kurtosisMax);
    g5->SetParLimits(3,kurtosisMin,kurtosisMax);

    g1->SetParLimits(4,skewMin,skewMax);
    g2->SetParLimits(4,skewMin,skewMax);
    g3->SetParLimits(4,skewMin,skewMax);
    g4->SetParLimits(4,skewMin,skewMax);
    g5->SetParLimits(4,skewMin,skewMax);

    if (fixedK || fixedS) SetFixedKSparameterForIndividualFits(g1,g2,g3,g4,g5);

    cout << " ============================================ " << endl;
    cout << " ElMean  = "   << arrMean[0]          ;
    cout << " ElSigma = "   << arrSigma[0]          << endl;
    cout << " PiMean  = "   << arrMean[1]          ;
    cout << " PiSigma = "   << arrSigma[1]          << endl;
    cout << " KaMean  = "   << arrMean[2]          ;
    cout << " KaSigma = "   << arrSigma[2]          << endl;
    cout << " PrMean  = "   << arrMean[3]          ;
    cout << " PrSigma = "   << arrSigma[3]          << endl;
    cout << " DeMean  = "   << arrMean[4]          ;
    cout << " DeSigma = "   << arrSigma[4]          << endl;

    // Clone histogram for total fit
    TH1D *hClone = (TH1D*)h1D->Clone();
    SetStyleFor1DSlice(hClone);
    SetStyleFor1DSlice(h1D);

    // Fit h1D for individual particles
    if (arrMean[0]>=10) h1D->Fit(g1,"QNR");
    if (arrMean[1]>=10) h1D->Fit(g2,"QNR+");
    if (arrMean[2]>=10) h1D->Fit(g3,"QNR+");
    if (arrMean[3]>=10) h1D->Fit(g4,"QNR+");
    if (arrMean[4]>=10) h1D->Fit(g5,"QNR+");

    // Get Fit parameters from the individual fits and set for the total fit
    Double_t par[25] = {0};
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[5]);
    g3->GetParameters(&par[10]);
    g4->GetParameters(&par[15]);
    g5->GetParameters(&par[20]);

    TF1 *total = new TF1("total","g1+g2+g3+g4+g5",dEdxMin,dEdxMax);
    total->SetLineWidth(3); total->SetLineColor(2); total->SetParameters(par);  total->SetNpx(nBinsInLookUpTable);

    // Some setters for the total fit
    SetTotalFitParameters(iIter,islice,nSlice,total,arrMean,arrSigma,arrMeanWindow,cleanParams,h1D,pt);
    SetParNamesOfTotalFit(total);
    CheckIfMeanIsZero(total,arrMean);

    // Apply total fit
    // gStyle->SetOptFit(1111);
    h1D->Fit(total,"QMR+");

    // new hists for KS test
    TH1D *hFitKS    = (TH1D*)h1D->Clone();
    TH1D *hTotalFit = FuncToHist(total,hFitKS,islice,arrMean,arrSigma);
    hTotalFit->SetName(Form("FitHist_%d",islice)); hFitKS ->SetName(Form("HistKS_%d",islice));
    hTotalFit->SetMarkerColor(kRed);               hFitKS ->SetMarkerColor(kBlue);
    hTotalFit->SetLineColor(kRed);                 hFitKS ->SetLineColor(kBlue);
    Double_t binWidthTotal = hTotalFit->GetXaxis()->GetBinWidth(50);
    hTotalFit->Scale(1./binWidthTotal);

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
    ResetParametersOfEachFunction(g1,g2,g3,g4,g5,total,h1D,pt);
    // //
    // // Create Line Shapes
    // PrepareLineShapes(hLineShape,total,g1,g2,g3,g4,g5,centBin,etaBin,momBin,signBin,islice);

    // Dump histograms into arrays
    h1D       -> SetName(Form("FinalFit_%d_%4.3f"  ,islice,pt));
    hResidual -> SetName(Form("Residual_%d_%4.3f"  ,islice,pt));
    if (dump) { debugFile->GetFile()->cd(); h1D -> Write(Form("FinalFit_%d_%4.3f",islice,pt)); }

    fitResArr[islice]    = (TH1D*)h1D;
    fitResiduals[islice] = (TH1D*)hResidual;
    //
    // make an animation file
    //////////////////////////////////////////////////////////////////////////////
    TH2D * h2Danimation = (TH2D*)h2D->Clone(); h2Danimation->SetName("h2Danimation");
    TH1D * h1ToGif      = (TH1D*)h1D->Clone();
    MakeAnimation(iIter, h1ToGif, h2Danimation, pt, ptMin, ptMax);
    //////////////////////////////////////////////////////////////////////////////
    //
    // Dump some info into TTreeStream

    elFitMax[islice]   = g1->GetMaximumX();
    piFitMax[islice]   = g2->GetMaximumX();
    kaFitMax[islice]   = g3->GetMaximumX();
    prFitMax[islice]   = g4->GetMaximumX();
    deFitMax[islice]   = g5->GetMaximumX();

    elFitAmp[islice]   = g1->GetParameter(0);
    piFitAmp[islice]   = g2->GetParameter(0);
    kaFitAmp[islice]   = g3->GetParameter(0);
    prFitAmp[islice]   = g4->GetParameter(0);
    deFitAmp[islice]   = g5->GetParameter(0);

    elFitMean[islice]  = g1->GetParameter(1);
    piFitMean[islice]  = g2->GetParameter(1);
    kaFitMean[islice]  = g3->GetParameter(1);
    prFitMean[islice]  = g4->GetParameter(1);
    deFitMean[islice]  = g5->GetParameter(1);

    elFitSigma[islice] = g1->GetParameter(2);
    piFitSigma[islice] = g2->GetParameter(2);
    kaFitSigma[islice] = g3->GetParameter(2);
    prFitSigma[islice] = g4->GetParameter(2);
    deFitSigma[islice] = g5->GetParameter(2);

    elFitKurtosis[islice]  = g1->GetParameter(3);
    piFitKurtosis[islice]  = g2->GetParameter(3);
    kaFitKurtosis[islice]  = g3->GetParameter(3);
    prFitKurtosis[islice]  = g4->GetParameter(3);
    deFitKurtosis[islice]  = g5->GetParameter(3);

    elFitSkew[islice]  = g1->GetParameter(4);
    piFitSkew[islice]  = g2->GetParameter(4);
    kaFitSkew[islice]  = g3->GetParameter(4);
    prFitSkew[islice]  = g4->GetParameter(4);
    deFitSkew[islice]  = g5->GetParameter(4);

    elChi2Fit[islice]  = g1->GetChisquare();
    piChi2Fit[islice]  = g2->GetChisquare();
    kaChi2Fit[islice]  = g3->GetChisquare();
    prChi2Fit[islice]  = g4->GetChisquare();
    deChi2Fit[islice]  = g5->GetChisquare();

    elIntegralFit[islice] = ComputeGaussIntegral(elFitAmp[islice],elFitMean[islice],elFitSigma[islice],g1->GetParameter(3),g1->GetParameter(4));
    piIntegralFit[islice] = ComputeGaussIntegral(piFitAmp[islice],piFitMean[islice],piFitSigma[islice],g2->GetParameter(3),g2->GetParameter(4));
    kaIntegralFit[islice] = ComputeGaussIntegral(kaFitAmp[islice],kaFitMean[islice],kaFitSigma[islice],g3->GetParameter(3),g3->GetParameter(4));
    prIntegralFit[islice] = ComputeGaussIntegral(prFitAmp[islice],prFitMean[islice],prFitSigma[islice],g4->GetParameter(3),g4->GetParameter(4));
    deIntegralFit[islice] = ComputeGaussIntegral(deFitAmp[islice],deFitMean[islice],deFitSigma[islice],g5->GetParameter(3),g5->GetParameter(4));

      // parSample[7] = asyGaus->GetMaximumX();

    Double_t totalChi2    = total->GetChisquare();
    Double_t normNDF      = total->GetNDF();
    Double_t totChi2      = (normNDF<1) ? 1 : totalChi2/normNDF;
    Double_t pValue       = total->GetProb();
    Double_t ksValue      = hFitKS->KolmogorovTest(hTotalFit);
    //     Double_t chi2Value    = hFitKS->Chi2Test(hTotalFit);

    // Get the fit param errors
    for (Int_t i=0;i<total->GetNpar();i++) parErrors[i] = total->GetParError(i);

    cout << " ElSkew = "   << skewnessFix[0]         ;
    cout << " ElKurt = "   << kurtosisFix[0]  << endl;
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

    if (corrStr=="") dEdxCorr = 0;
    if (corrStr=="Corr") dEdxCorr = 1;

    // Fit results dump ttrees
    fitFile->GetFile()->cd();
    *fitFile << "treeId" <<

    "corr="      << dEdxCorr                <<
    "sl="        << islice                <<
    "sign="      << fSign                  <<
    "it="        << iIter                 <<
    "p="      << pt                    <<
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

    *fitFile << "treeId" <<

    "elMax="    << elFitMax[islice]      <<
    "piMax="    << piFitMax[islice]      <<
    "kaMax="    << kaFitMax[islice]      <<
    "prMax="    << prFitMax[islice]      <<
    "deMax="    << deFitMax[islice]      <<

    "elA="    << elFitAmp[islice]      <<
    "piA="    << piFitAmp[islice]      <<
    "kaA="    << kaFitAmp[islice]      <<
    "prA="    << prFitAmp[islice]      <<
    "deA="    << deFitAmp[islice]      <<

    "elM="    << elFitMean[islice]     <<
    "piM="    << piFitMean[islice]     <<
    "kaM="    << kaFitMean[islice]     <<
    "prM="    << prFitMean[islice]     <<
    "deM="    << deFitMean[islice]     <<

    "elSi="   << elFitSigma[islice]    <<
    "piSi="   << piFitSigma[islice]    <<
    "kaSi="   << kaFitSigma[islice]    <<
    "prSi="   << prFitSigma[islice]    <<
    "deSi="   << deFitSigma[islice]    ;

    *fitFile << "treeId" <<

    "elK="    << elFitKurtosis[islice]     <<
    "piK="    << piFitKurtosis[islice]     <<
    "kaK="    << kaFitKurtosis[islice]     <<
    "prK="    << prFitKurtosis[islice]     <<
    "deK="    << deFitKurtosis[islice]     <<

    "elSk="   << elFitSkew[islice]     <<
    "piSk="   << piFitSkew[islice]     <<
    "kaSk="   << kaFitSkew[islice]     <<
    "prSk="   << prFitSkew[islice]     <<
    "deSk="   << deFitSkew[islice]     <<

    "elAmpErr="      << parErrors[0]      <<
    "piAmpErr="      << parErrors[5]      <<
    "kaAmpErr="      << parErrors[10]     <<
    "prAmpErr="      << parErrors[15]     <<
    "deAmpErr="      << parErrors[20]     ;

    *fitFile << "treeId" <<

    "elMeanErr="     << parErrors[1]      <<
    "piMeanErr="     << parErrors[6]      <<
    "kaMeanErr="     << parErrors[11]      <<
    "prMeanErr="     << parErrors[16]      <<
    "deMeanErr="     << parErrors[21]      <<

    "elSigmaErr="    << parErrors[2]      <<
    "piSigmaErr="    << parErrors[7]      <<
    "kaSigmaErr="    << parErrors[12]      <<
    "prSigmaErr="    << parErrors[17]      <<
    "deSigmaErr="    << parErrors[22]      <<

    "elKurtosisErr=" << parErrors[3]      <<
    "piKurtosisErr=" << parErrors[8]      <<
    "kaKurtosisErr=" << parErrors[13]      <<
    "prKurtosisErr=" << parErrors[18]      <<
    "deKurtosisErr=" << parErrors[23]      <<

    "elSkewErr="     << parErrors[4]      <<
    "piSkewErr="     << parErrors[9]      <<
    "kaSkewErr="     << parErrors[14]      <<
    "prSkewErr="     << parErrors[19]      <<
    "deSkewErr="     << parErrors[24]      <<

    "\n";

    // Fit results dump ttrees
    *fitFile << "FitResults" <<

    "elMeanMod="    << arrMeanModif[0]         <<
    "piMeanMod="    << arrMeanModif[1]         <<
    "kaMeanMod="    << arrMeanModif[2]         <<
    "prMeanMod="    << arrMeanModif[3]         <<
    "deMeanMod="    << arrMeanModif[4]         <<
    "kaTOFMeanMod=" << arrMeanModif[5]         <<
    "prTOFMeanMod=" << arrMeanModif[6]         <<


    "elMeanSpline="       << arrMean[0]         <<
    "piMeanSpline="       << arrMean[1]         <<
    "kaMeanSpline="       << arrMean[2]         <<
    "prMeanSpline="       << arrMean[3]         <<
    "deMeanSpline="       << arrMean[4]         <<
    "kaTOFMeanSpline="    << arrMean[5]         <<
    "prTOFMeanSpline="    << arrMean[6]         <<


    "elSigmaSpline="      << arrSigma[0]        <<
    "piSigmaSpline="      << arrSigma[1]        <<
    "kaSigmaSpline="      << arrSigma[2]        <<
    "prSigmaSpline="      << arrSigma[3]        <<
    "deSigmaSpline="      << arrSigma[4]        <<
    "kaTOFSigmaSpline="   << arrSigma[5]        <<
    "prTOFSigmaSpline="   << arrSigma[6]        <<


    "elWindow="     << arrMeanWindow[0]   <<
    "piWindow="     << arrMeanWindow[1]   <<
    "kaWindow="     << arrMeanWindow[2]   <<
    "prWindow="     << arrMeanWindow[3]   <<
    "deWindow="     << arrMeanWindow[4]   ;


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

    "elFitMax="     << elFitMax[islice] <<
    "piFitMax="     << piFitMax[islice] <<
    "kaFitMax="     << kaFitMax[islice] <<
    "prFitMax="     << prFitMax[islice] <<
    "deFitMax="     << deFitMax[islice] <<

    "elFitInt="     << elIntegralFit[islice] <<
    "piFitInt="     << piIntegralFit[islice] <<
    "kaFitInt="     << kaIntegralFit[islice] <<
    "prFitInt="     << prIntegralFit[islice] <<
    "deFitInt="     << deIntegralFit[islice] <<

    "elFitMean="    << elFitMean[islice]  <<
    "piFitMean="    << piFitMean[islice]  <<
    "kaFitMean="    << kaFitMean[islice]  <<
    "prFitMean="    << prFitMean[islice]  <<
    "deFitMean="    << deFitMean[islice]  ;

    *fitFile << "FitResults" <<

    "elFitSigma="   << elFitSigma[islice] <<
    "piFitSigma="   << piFitSigma[islice] <<
    "kaFitSigma="   << kaFitSigma[islice] <<
    "prFitSigma="   << prFitSigma[islice] <<
    "deFitSigma="   << deFitSigma[islice] <<

    "elFitAmp="     << elFitAmp[islice]   <<
    "piFitAmp="     << piFitAmp[islice]   <<
    "kaFitAmp="     << kaFitAmp[islice]   <<
    "prFitAmp="     << prFitAmp[islice]   <<
    "deFitAmp="     << deFitAmp[islice]   <<

    "elKurtosis="    << elFitKurtosis[islice]     <<
    "piKurtosis="    << piFitKurtosis[islice]     <<
    "kaKurtosis="    << kaFitKurtosis[islice]     <<
    "prKurtosis="    << prFitKurtosis[islice]     <<
    "deKurtosis="    << deFitKurtosis[islice]     ;

    *fitFile << "FitResults" <<

    "elFitSkew="    << elFitSkew[islice]  <<
    "piFitSkew="    << piFitSkew[islice]  <<
    "kaFitSkew="    << kaFitSkew[islice]  <<
    "prFitSkew="    << prFitSkew[islice]  <<
    "deFitSkew="    << deFitSkew[islice]  <<

    "elChi2Fit="    << elChi2Fit[islice]  <<
    "piChi2Fit="    << piChi2Fit[islice]  <<
    "kaChi2Fit="    << kaChi2Fit[islice]  <<
    "prChi2Fit="    << prChi2Fit[islice]  <<
    "deChi2Fit="    << deChi2Fit[islice]  <<

    "elAmpErr="      << parErrors[0]      <<
    "piAmpErr="      << parErrors[5]      <<
    "kaAmpErr="      << parErrors[10]      <<
    "prAmpErr="      << parErrors[15]      <<
    "deAmpErr="      << parErrors[20]      ;

    *fitFile << "FitResults" <<

    "elMeanErr="     << parErrors[1]      <<
    "piMeanErr="     << parErrors[6]      <<
    "kaMeanErr="     << parErrors[11]      <<
    "prMeanErr="     << parErrors[16]      <<
    "deMeanErr="     << parErrors[21]      <<

    "elSigmaErr="    << parErrors[2]      <<
    "piSigmaErr="    << parErrors[7]      <<
    "kaSigmaErr="    << parErrors[12]      <<
    "prSigmaErr="    << parErrors[17]      <<
    "deSigmaErr="    << parErrors[22]      <<

    "elKurtosisErr=" << parErrors[3]      <<
    "piKurtosisErr=" << parErrors[8]      <<
    "kaKurtosisErr=" << parErrors[13]      <<
    "prKurtosisErr=" << parErrors[18]      <<
    "deKurtosisErr=" << parErrors[23]      <<

    "elSkewErr="     << parErrors[4]      <<
    "piSkewErr="     << parErrors[9]      <<
    "kaSkewErr="     << parErrors[14]      <<
    "prSkewErr="     << parErrors[19]      <<
    "deSkewErr="     << parErrors[24]      <<

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
    delete h2Danimation;
    delete h1ToGif;

  }


  histsFile -> GetFile()->cd();
  if (dump) {
    fitResArr . Write("fitResArr",TObject::kSingleKey);
    fitResiduals . Write("fitResiduals",TObject::kSingleKey);
  } else {
    if (iIter ==0 || iIter >=3) fitResArr . Write("fitResArr",TObject::kSingleKey);
    if (iIter >=3)  fitResiduals . Write("fitResiduals",TObject::kSingleKey);
  }
  // //
  // // Fill line shapes
  // lineShapesLookUp->GetFile()->cd();
  // histLineShapesCArr . Write("histLineShapesCArr",TObject::kSingleKey);
  // histLineShapesCArr.Clear("C");
  // funcLineShapesCArr . Write("funcLineShapesCArr",TObject::kSingleKey);
  // funcLineShapesCArr.Clear("C");
  // funcLineShapesObjArr . Write("funcLineShapesObjArr",TObject::kSingleKey);
  // funcLineShapesObjArr.Delete();

  // delete arrays
  fitResArr    . Clear("C");
  fitResiduals . Clear("C");

  delete histsFile;
  delete fitFile;
  // delete lineShapesLookUp;

  if (ressFile!=0x0) ressFile->Close();
  if (sampFile!=0x0) sampFile->Close();

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

    //
    //   set initials params using PIDresponse
    if (iIter>=0) SetParams1stIteration(pt,h1D,arrMean,arrSigma,arrMeanWindow,cleanParams);

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=1) {

      Double_t elAssump=0., piAssump=0., kaAssump=0., prAssump=0.;
      if (assump==0){    // use only PID response
        elAssump = grMeanModExScaled[kEl]       -> GetY()[islice];
        piAssump = grCleanMeanScaledSmooth[kPi] -> GetY()[islice];
        kaAssump = grMeanModExScaled[kKa]       -> GetY()[islice];
        prAssump = grMeanModExScaled[kPr]       -> GetY()[islice];
      } else if (assump==1){ // use clean samples
        elAssump = (pt<0.6) ? grCleanMeanScaledSmooth[kEl] -> GetY()[islice]: grMeanExScaled[kEl] -> GetY()[islice];
        piAssump =            grCleanMeanScaledSmooth[kPi] -> GetY()[islice];
        kaAssump = (pt>0.5) ? grCleanMeanScaledSmooth[kKa] -> GetY()[islice]: grMeanExScaled[kKa] -> GetY()[islice];
        prAssump = (pt>0.7) ? grCleanMeanScaledSmooth[kPr] -> GetY()[islice]: grMeanExScaled[kPr] -> GetY()[islice];
      }
      // Recalculate the windows with the clean samples
      // ***************  Electron - Kaon ***************
      if ((pt>=0.45 && pt<=0.65) || (pt>=0.8 && pt<=1.2) ){
        elAmpMax = grFitAmpSmooth[kEl]->GetY()[islice] + grFitAmpSmooth[kEl]->GetY()[islice]*0.05;
      }
      if ( pt>=1 && pt<=2.8 ){
        prAmpMin   = grCleanAmpScaledSmooth[kPrTOF]->GetY()[islice] - grCleanAmpScaledSmooth[kPrTOF]->GetY()[islice]*0.1;
        prAmpMax   = grCleanAmpScaledSmooth[kPrTOF]->GetY()[islice] + grCleanAmpScaledSmooth[kPrTOF]->GetY()[islice]*0.1;
      }
      //
      if ( (pt>=0.45 && pt<=0.65) ){
        elMeanMin  = grCleanMeanScaledSmooth[kEl]->GetY()[islice] - grCleanMeanScaledSmooth[kEl]  ->GetY()[islice]*0.005;
        elMeanMax  = grCleanMeanScaledSmooth[kEl]->GetY()[islice] + grCleanMeanScaledSmooth[kEl]  ->GetY()[islice]*0.005;
        kaMeanMin  = grMeanExScaled[kKa]->GetY()[islice] - grMeanExScaled[kKa]  ->GetY()[islice]*0.03;
        kaMeanMax  = grMeanExScaled[kKa]->GetY()[islice] + grMeanExScaled[kKa]  ->GetY()[islice]*0.03;
      }
      //
      //
      if ( (pt>=0.8) ){
        piAmpMin = grFitAmpSmooth[kPi]->GetY()[islice] - grFitAmpSmooth[kPi]->GetY()[islice]*0.05;
        piAmpMax = grFitAmpSmooth[kPi]->GetY()[islice] + grFitAmpSmooth[kPi]->GetY()[islice]*0.005;
      }
      //
      // ***************  deuteron ***************
      if (fDeuteronsExist) {
        deMeanMin  = grMeanExScaled[kDe]->GetY()[islice]  - grMeanExScaled[kDe]  ->GetY()[islice]*0.005;
        deMeanMax  = grMeanExScaled[kDe]->GetY()[islice]  + grMeanExScaled[kDe]  ->GetY()[islice]*0.005;
        deSigmaMin = grSigmaExScaled[kDe]->GetY()[islice] - grSigmaExScaled[kDe] ->GetY()[islice]*0.005;
        deSigmaMax = grSigmaExScaled[kDe]->GetY()[islice] + grSigmaExScaled[kDe] ->GetY()[islice]*0.005;
        deAmpMin   = grCleanAmpScaledSmooth[kDe]->GetY()[islice] - grCleanAmpScaledSmooth[kDe]->GetY()[islice]*0.1;
        deAmpMax   = grCleanAmpScaledSmooth[kDe]->GetY()[islice] + grCleanAmpScaledSmooth[kDe]->GetY()[islice]*0.1;
      } else{
        deMeanMin  = 900.;
        deMeanMax  = 1000.;
        deSigmaMin = 0.;
        deSigmaMax = 1.;
        deAmpMin   = 0.;
        deAmpMax   = 1.;
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=2) {
      //
      if ((pt>=0.45 && pt<=0.65) || (pt>=0.8 && pt<=1.2) ){
        elMeanMin = grMeanCleanWrtExScaled[kEl]->GetY()[islice] - grMeanCleanWrtExScaled[kEl]->GetY()[islice]*0.005;
        elMeanMax = grMeanCleanWrtExScaled[kEl]->GetY()[islice] + grMeanCleanWrtExScaled[kEl]->GetY()[islice]*0.005;
        elSigmaMin = grSigmaCleanWrtExScaled[kEl]->GetY()[islice] - grSigmaCleanWrtExScaled[kEl] ->GetY()[islice]*0.01;
        elSigmaMax = grSigmaCleanWrtExScaled[kEl]->GetY()[islice] + grSigmaCleanWrtExScaled[kEl] ->GetY()[islice]*0.01;
      }
      // ***************  Electron - Kaon ***************
      if ( pt>=0.65 ){
        elMeanMin = grMeanCleanWrtExScaled[kEl]->GetY()[islice] - grMeanCleanWrtExScaled[kEl]->GetY()[islice]*0.005;
        elMeanMax = grMeanCleanWrtExScaled[kEl]->GetY()[islice] + grMeanCleanWrtExScaled[kEl]->GetY()[islice]*0.005;
        elSigmaMin = grSigmaCleanWrtExScaled[kEl]->GetY()[islice] - grSigmaCleanWrtExScaled[kEl]->GetY()[islice]*0.01;
        elSigmaMax = grSigmaCleanWrtExScaled[kEl]->GetY()[islice] + grSigmaCleanWrtExScaled[kEl]->GetY()[islice]*0.01;
        //
        piMeanMin = grCleanMeanScaled[kPi]->GetY()[islice] - grCleanMeanScaled[kPi]->GetY()[islice]*0.005;
        piMeanMax = grCleanMeanScaled[kPi]->GetY()[islice] + grCleanMeanScaled[kPi]->GetY()[islice]*0.005;
        piSigmaMin = grCleanSigmaScaled[kPi]->GetY()[islice] - grCleanSigmaScaled[kPi] ->GetY()[islice]*0.005;
        piSigmaMax = grCleanSigmaScaled[kPi]->GetY()[islice] + grCleanSigmaScaled[kPi] ->GetY()[islice]*0.005;

      }
        if ( pt>=0.8 ){
        piMeanMin = grCleanMeanScaledSmooth[kPi]->GetY()[islice] - grCleanMeanScaledSmooth[kPi]->GetY()[islice]*0.005;
        piMeanMax = grCleanMeanScaledSmooth[kPi]->GetY()[islice] + grCleanMeanScaledSmooth[kPi]->GetY()[islice]*0.005;
        piSigmaMin = grCleanSigmaScaledSmooth[kPi]->GetY()[islice] - grCleanSigmaScaledSmooth[kPi] ->GetY()[islice]*0.005;
        piSigmaMax = grCleanSigmaScaledSmooth[kPi]->GetY()[islice] + grCleanSigmaScaledSmooth[kPi] ->GetY()[islice]*0.005;
        piAmpMin  = grFitAmpSmooth[kPi]->GetY()[islice] - grFitAmpSmooth[kPi]->GetY()[islice]*0.01;
        piAmpMax  = grFitAmpSmooth[kPi]->GetY()[islice] + grFitAmpSmooth[kPi]->GetY()[islice]*0.01;
      }
      //
      if ( pt>=0.75 ){
        kaMeanMin = grCleanMeanScaled[kKa]->GetY()[islice] - grCleanMeanScaled[kKa]->GetY()[islice]*0.01;
        kaMeanMax = grCleanMeanScaled[kKa]->GetY()[islice] + grCleanMeanScaled[kKa]->GetY()[islice]*0.01;
        kaSigmaMin = grCleanSigmaScaled[kKa]->GetY()[islice] - grCleanSigmaScaled[kKa] ->GetY()[islice]*0.03;
        kaSigmaMax = grCleanSigmaScaled[kKa]->GetY()[islice] + grCleanSigmaScaled[kKa] ->GetY()[islice]*0.03;
      }

      if ( pt>=2.2 ){
        kaMeanMin = grCleanMeanScaled[kKa]->GetY()[islice] - grCleanMeanScaled[kKa]->GetY()[islice]*0.05;
        kaMeanMax = grCleanMeanScaled[kKa]->GetY()[islice] + grCleanMeanScaled[kKa]->GetY()[islice]*0.05;
        kaSigmaMin = grCleanSigmaScaled[kKa]->GetY()[islice] - grCleanSigmaScaled[kKa] ->GetY()[islice]*0.1;
        kaSigmaMax = grCleanSigmaScaled[kKa]->GetY()[islice] + grCleanSigmaScaled[kKa] ->GetY()[islice]; // contamination is so much that sigma can not grow further
      }
      //
      if ( pt>=0.8 ){
        prMeanMin  = grCleanMeanScaled[kPr]->GetY()[islice] - grCleanMeanScaled[kPr] ->GetY()[islice]*0.005;
        prMeanMax  = grCleanMeanScaled[kPr]->GetY()[islice] + grCleanMeanScaled[kPr] ->GetY()[islice]*0.005;
        prSigmaMin = grCleanSigmaScaled[kPr]->GetY()[islice] - grCleanSigmaScaled[kPr] ->GetY()[islice]*0.005;
        prSigmaMax = grCleanSigmaScaled[kPr]->GetY()[islice] + grCleanSigmaScaled[kPr] ->GetY()[islice]*0.005;
      }
    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=3) {

      if ((pt>=0.45 && pt<=0.65) || (pt>=0.8 && pt<=1.2) ){
        elAmpMin = grFitAmpSmooth[kEl]->GetY()[islice] - grFitAmpSmooth[kEl]->GetY()[islice]*0.01;
        elAmpMax = grFitAmpSmooth[kEl]->GetY()[islice] + grFitAmpSmooth[kEl]->GetY()[islice]*0.01;
      }

      if ( pt>=0.65 ){
        piAmpMin   = grFitAmpSmooth[kPi]->GetY()[islice] - grFitAmpSmooth[kPi]->GetY()[islice]*0.005;
        piAmpMax   = grFitAmpSmooth[kPi]->GetY()[islice] + grFitAmpSmooth[kPi]->GetY()[islice]*0.005;
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=4) {

      if ( pt>=0.4 ){
        elMeanMin  = grFitMeanSmooth[kEl]->GetY()[islice]  - grFitMeanSmooth[kEl]->GetY()[islice]*0.005;
        elMeanMax  = grFitMeanSmooth[kEl]->GetY()[islice]  + grFitMeanSmooth[kEl]->GetY()[islice]*0.005;
        elSigmaMin = grFitSigmaSmooth[kEl]->GetY()[islice] - grFitSigmaSmooth[kEl]->GetY()[islice]*0.01;
        elSigmaMax = grFitSigmaSmooth[kEl]->GetY()[islice] + grFitSigmaSmooth[kEl]->GetY()[islice]*0.01;
      }

      if ( pt>=0.8 ){
        piMeanMin  = grFitMeanSmooth[kPi]->GetY()[islice] - grFitMeanSmooth[kPi]->GetY()[islice]*0.005;
        piMeanMax  = grFitMeanSmooth[kPi]->GetY()[islice] + grFitMeanSmooth[kPi]->GetY()[islice]*0.005;
        piSigmaMin = grFitSigmaSmooth[kPi]->GetY()[islice] - grFitSigmaSmooth[kPi] ->GetY()[islice]*0.005;
        piSigmaMax = grFitSigmaSmooth[kPi]->GetY()[islice] + grFitSigmaSmooth[kPi] ->GetY()[islice]*0.005;
        //
        prAmpMin  = grFitAmpSmooth[kPr]->GetY()[islice] - grFitAmpSmooth[kPr]->GetY()[islice]*0.05;
        prAmpMax  = grFitAmpSmooth[kPr]->GetY()[islice] + grFitAmpSmooth[kPr]->GetY()[islice]*0.05;
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=5) {

      if ( pt>=0.6 ){
        prAmpMin  = grFitAmpSmooth[kPr]->GetY()[islice] - grFitAmpSmooth[kPr]->GetY()[islice]*0.005;
        prAmpMax  = grFitAmpSmooth[kPr]->GetY()[islice] + grFitAmpSmooth[kPr]->GetY()[islice]*0.005;
        kaAmpMin  = grFitAmpSmooth[kKa]->GetY()[islice] - grFitAmpSmooth[kKa]->GetY()[islice]*0.005;
        kaAmpMax  = grFitAmpSmooth[kKa]->GetY()[islice] + grFitAmpSmooth[kKa]->GetY()[islice]*0.005;
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (iIter>=6) {

      if ( pt>=0.6 ){
        kaAmpMin  = grFitAmpSmooth[kKa]->GetY()[islice] - grFitAmpSmooth[kKa]->GetY()[islice]*0.001;
        kaAmpMax  = grFitAmpSmooth[kKa]->GetY()[islice] + grFitAmpSmooth[kKa]->GetY()[islice]*0.001;
        //
        prAmpMin  = grFitAmpSmooth[kPr]->GetY()[islice] - grFitAmpSmooth[kPr]->GetY()[islice]*0.001;
        prAmpMax  = grFitAmpSmooth[kPr]->GetY()[islice] + grFitAmpSmooth[kPr]->GetY()[islice]*0.001;
      }

    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  }  // close settings for 20MeV bin width
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //   Almost no restriction in low momentum
  if (pt<0.3){

    Double_t tmpmaxBin  = h1D->GetMaximum();

    elAmpMin = 1e-4;  elAmpMax = tmpmaxBin;
    piAmpMin = 1e-4;  piAmpMax = tmpmaxBin;
    kaAmpMin = 1e-4;  kaAmpMax = tmpmaxBin;
    prAmpMin = 1e-4;  prAmpMax = tmpmaxBin;
    deAmpMin = 1e-4;  deAmpMax = tmpmaxBin;

    elMeanMin = TMath::Max(30.,(arrMean[kEl]-1.*arrSigma[kEl]));
    piMeanMin = TMath::Max(30.,(arrMean[kPi]-1.*arrSigma[kPi]));
    kaMeanMin = TMath::Max(30.,(arrMean[kKa]-4.*arrSigma[kKa]));

    elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kEl]+1.*arrSigma[kEl]));
    piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kPi]+1.*arrSigma[kPi]));
    kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kKa]+4.*arrSigma[kKa]));

    // pr and de is special for low momentum
    prMeanMin = TMath::Max(30.,(arrMean[3]-4.*arrSigma[3]));
    prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+4.*arrSigma[3]));

    deMeanMin = TMath::Max(30.,(arrMean[4]-4.*arrSigma[4]));
    deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[4]+4.*arrSigma[4]));

    elSigmaMin = arrSigma[0]/5.;  elSigmaMax = arrSigma[0]*3.;
    piSigmaMin = arrSigma[1]/5.;  piSigmaMax = arrSigma[1]*3.;
    kaSigmaMin = arrSigma[2]/5.;  kaSigmaMax = arrSigma[2]*3.;
    prSigmaMin = arrSigma[3]/5.;  prSigmaMax = arrSigma[3]*3.;
    deSigmaMin = arrSigma[4]/5.;  deSigmaMax = arrSigma[4]*3.;

  }

  //   Let protons free below 0.7 gev
  if (pt<0.7){

    // ------- proton -------
    prAmpMin = 1e-4;
    prAmpMax = h1D->GetMaximum();
    prMeanMin = TMath::Max(30.,(arrMean[kPr]-3.*arrSigma[kPr]));
    prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kPr]+3.*arrSigma[kPr]));
    prSigmaMin = arrSigma[3]/5.;
    prSigmaMax = arrSigma[3]*3.;

  }

  //   Let deuterons be free below 0.7 gev
  if (pt<1){

    // ------- proton -------
    deAmpMin = 1e-4;
    deAmpMax = h1D->GetMaximum();
    deMeanMin = TMath::Max(150.,(arrMean[kDe]-3.*arrSigma[kDe]));
    deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[kDe]+3.*arrSigma[kDe]));
    deSigmaMin = arrSigma[4]/5.;
    deSigmaMax = arrSigma[4]*3.;

  }

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //   Finally set the params

  // set the Kurtosis and Skewness as free, fix them if needed (incase modify fixedK and fixedS)
  total->SetParLimits(4,skewMin,skewMax);    // electron
  total->SetParLimits(9,skewMin,skewMax);    // pion
  total->SetParLimits(14,skewMin,skewMax);   // kaon
  total->SetParLimits(19,skewMin,skewMax);   // proton
  total->SetParLimits(24,skewMin,skewMax);   // deuteron

  total->SetParLimits(3,kurtosisMin,kurtosisMax);   // electron
  total->SetParLimits(8,kurtosisMin,kurtosisMax);   // pion
  total->SetParLimits(13,kurtosisMin,kurtosisMax);  // kaon
  total->SetParLimits(18,kurtosisMin,kurtosisMax);  // proton
  total->SetParLimits(23,kurtosisMin,kurtosisMax);  // deuteron

  // if needed fix Kurtosis and Skewness
  if ( (fixedK || fixedS) ) SetFixedKSparameterForTotalFit(total,pt);

  // ++++++++++++++++++++++++ more freedom for some particles +++++++++++++++++++

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

}
// -------------------------------------------------------------------------------------------------------
void SmoothAmplitudes(TString fileFit, TString fileSample, TString fileOut, const Int_t nSlice, Int_t iIter)
{

  //
  // Analyse first results and smooth the amplitude graphs
  /*
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/tests/trash
  aliroot -l
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/RealData_PIDIterativeFitting.C+
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

  TObjArray *GrArrScaledExWrtClean  = new TObjArray(2*nArrEntriesClean);  GrArrScaledExWrtClean   -> SetOwner(kTRUE);
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
  //
  TGraphErrors **grScaledExWrtClean  = new TGraphErrors *[nArrEntries];
  TGraph       **grScaledExWrtCleanSmooth = new TGraph *[nArrEntries];

  //   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //   Produce smooth and normal graphs
  ProduceGraphs(nSlice, treeFit  , GrArrFit  , branchArrFit  , grFit  , grFitSmooth);
  ProduceGraphs(nSlice, treeClean, GrArrClean, branchArrClean, grClean, grCleanSmooth);

  //   Apply outlier removal to fits
  ApplyOutlierSmooth(nSlice, treeFit, GrArrOutlierSmooth, branchArrFit, grOutlierCut, grOutlierSmooth);

  //   Apply Scalings for expected mean and sigma
  Double_t scalingPointKa2 = 0.65; // clean and incspectra has stat only in
  Double_t scalingPointPr2 = 1.21; // clean and incspectra has stat only in
  //
  Double_t scalingPointEl  = 0.35;
  Double_t scalingPointPi  = 0.55;
  Double_t scalingPointKa  = 0.41;
  Double_t scalingPointPr  = 0.75;
  Double_t scalingPointDe  = 1.2;

  grScaledEx[0] = ApplyScalingWrtExpected(nSlice, "elMeanSpline"   ,"elFitMean",scalingPointEl,GrArrFit,treeFit);
  grScaledEx[1] = ApplyScalingWrtExpected(nSlice, "piMeanSpline"   ,"piFitMean",scalingPointPi,GrArrFit,treeFit);
  grScaledEx[2] = ApplyScalingWrtExpected(nSlice, "kaMeanSpline"   ,"kaFitMean",scalingPointKa,GrArrFit,treeFit);
  grScaledEx[3] = ApplyScalingWrtExpected(nSlice, "prMeanSpline"   ,"prFitMean",scalingPointPr,GrArrFit,treeFit);
  grScaledEx[4] = ApplyScalingWrtExpected(nSlice, "deMeanSpline"   ,"deFitMean",scalingPointDe,GrArrFit,treeFit);
  grScaledEx[5] = ApplyScalingWrtExpected(nSlice, "kaTOFMeanSpline","kaFitMean",scalingPointKa,GrArrFit,treeFit);
  grScaledEx[6] = ApplyScalingWrtExpected(nSlice, "prTOFMeanSpline","prFitMean",scalingPointPr,GrArrFit,treeFit);


  grScaledEx[7] = ApplyScalingWrtExpected(nSlice, "elSigmaSpline"    ,"elFitSigma",scalingPointEl,GrArrFit,treeFit);
  grScaledEx[8] = ApplyScalingWrtExpected(nSlice, "piSigmaSpline"    ,"piFitSigma",scalingPointPi,GrArrFit,treeFit);
  grScaledEx[9] = ApplyScalingWrtExpected(nSlice, "kaSigmaSpline"    ,"kaFitSigma",scalingPointKa,GrArrFit,treeFit);
  grScaledEx[10] = ApplyScalingWrtExpected(nSlice, "prSigmaSpline"   ,"prFitSigma",scalingPointPr,GrArrFit,treeFit);
  grScaledEx[11] = ApplyScalingWrtExpected(nSlice, "deSigmaSpline"   ,"deFitSigma",scalingPointDe,GrArrFit,treeFit);
  grScaledEx[12] = ApplyScalingWrtExpected(nSlice, "kaTOFSigmaSpline","kaFitSigma",scalingPointKa,GrArrFit,treeFit);
  grScaledEx[13] = ApplyScalingWrtExpected(nSlice, "prTOFSigmaSpline","prFitSigma",scalingPointPr,GrArrFit,treeFit);


  grScaledEx[14] = ApplyScalingWrtExpected(nSlice, "elMeanMod"   ,"elFitMean",scalingPointEl,GrArrFit,treeFit);
  grScaledEx[15] = ApplyScalingWrtExpected(nSlice, "piMeanMod"   ,"piFitMean",scalingPointPi,GrArrFit,treeFit);
  grScaledEx[16] = ApplyScalingWrtExpected(nSlice, "kaMeanMod"   ,"kaFitMean",scalingPointKa,GrArrFit,treeFit);
  grScaledEx[17] = ApplyScalingWrtExpected(nSlice, "prMeanMod"   ,"prFitMean",scalingPointPr,GrArrFit,treeFit);
  grScaledEx[18] = ApplyScalingWrtExpected(nSlice, "kaTOFMeanMod","kaFitMean",scalingPointKa,GrArrFit,treeFit);
  grScaledEx[19] = ApplyScalingWrtExpected(nSlice, "prTOFMeanMod","prFitMean",scalingPointPr,GrArrFit,treeFit);

  for (Int_t i=0; i<20; i++) GrArrScaledEx->AddAt(grScaledEx[i],i);

  //   Apply Scalings clean data mean and sigma KAON and PROTON AMPLITUDE !!!!!!
  grScaledClean[0]  = ApplyScalingWrtSample(nSlice,"elCleanAmp"   ,"elFitAmp"  ,scalingPointEl ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,0);
  grScaledClean[1]  = ApplyScalingWrtSample(nSlice,"piCleanAmp"   ,"piFitAmp"  ,scalingPointPi ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,1);
  grScaledClean[2]  = ApplyScalingWrtSample(nSlice,"kaCleanAmp"   ,"kaFitAmp"  ,scalingPointKa2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,2);
  grScaledClean[3]  = ApplyScalingWrtSample(nSlice,"prCleanAmp"   ,"prFitAmp"  ,scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,3);
  grScaledClean[4]  = ApplyScalingWrtSample(nSlice,"deCleanAmp"   ,"deFitAmp"  ,scalingPointDe ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,4);
  grScaledClean[5]  = ApplyScalingWrtSample(nSlice,"kaTOFCleanAmp","kaFitAmp"  ,scalingPointKa2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,5);
  grScaledClean[6]  = ApplyScalingWrtSample(nSlice,"prTOFCleanAmp","prFitAmp"  ,scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,6);

  grScaledClean[7]  = ApplyScalingWrtSample(nSlice,"elCleanMean"    ,"elFitMean",scalingPointEl ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,7);
  grScaledClean[8]  = ApplyScalingWrtSample(nSlice,"piCleanMean"    ,"piFitMean",scalingPointPi ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,8);
  grScaledClean[9]  = ApplyScalingWrtSample(nSlice,"kaCleanMean"    ,"kaFitMean",scalingPointKa2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,9);
  grScaledClean[10]  = ApplyScalingWrtSample(nSlice,"prCleanMean"   ,"prFitMean",scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,10);
  grScaledClean[11]  = ApplyScalingWrtSample(nSlice,"deCleanMean"   ,"deFitMean",scalingPointDe ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,11);
  grScaledClean[12]  = ApplyScalingWrtSample(nSlice,"kaTOFCleanMean","kaFitMean",scalingPointKa2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,12);
  grScaledClean[13]  = ApplyScalingWrtSample(nSlice,"prTOFCleanMean","prFitMean",scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,13);

  grScaledClean[14] = ApplyScalingWrtSample(nSlice,"elCleanSigma"   ,"elFitSigma",scalingPointEl ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,14);
  grScaledClean[15] = ApplyScalingWrtSample(nSlice,"piCleanSigma"   ,"piFitSigma",scalingPointPi ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,15);
  grScaledClean[16] = ApplyScalingWrtSample(nSlice,"kaCleanSigma"   ,"kaFitSigma",scalingPointKa2 ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,16);
  grScaledClean[17] = ApplyScalingWrtSample(nSlice,"prCleanSigma"   ,"prFitSigma",scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,17);
  grScaledClean[18] = ApplyScalingWrtSample(nSlice,"deCleanSigma"   ,"deFitSigma",scalingPointDe ,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,18);
  grScaledClean[19] = ApplyScalingWrtSample(nSlice,"kaTOFCleanSigma","kaFitSigma",scalingPointKa2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,19);
  grScaledClean[20] = ApplyScalingWrtSample(nSlice,"prTOFCleanSigma","prFitSigma",scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledCleanSmooth,20);


  for (Int_t i=0; i<21; i++) {
    GrArrScaledClean->AddAt(grScaledClean[i],2*i);
    GrArrScaledClean->AddAt(grScaledCleanSmooth[i],2*i+1);
  }
  //   Apply Scalings splines mean and sigma wrt clean samples  !!!!!!
  grScaledExWrtClean[0]  = ApplyScalingWrtSample(nSlice,"elMeanSpline"   ,"elCleanMean"   ,0.6            ,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,0);
  grScaledExWrtClean[1]  = ApplyScalingWrtSample(nSlice,"piMeanSpline"   ,"piCleanMean"   ,0.7            ,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,1);
  grScaledExWrtClean[2]  = ApplyScalingWrtSample(nSlice,"kaMeanSpline"   ,"kaCleanMean"   ,scalingPointKa2,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,2);
  grScaledExWrtClean[3]  = ApplyScalingWrtSample(nSlice,"prMeanSpline"   ,"prCleanMean"   ,scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,3);
  grScaledExWrtClean[4]  = ApplyScalingWrtSample(nSlice,"deMeanSpline"   ,"deCleanMean"   ,scalingPointDe ,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,4);
  grScaledExWrtClean[5]  = ApplyScalingWrtSample(nSlice,"kaTOFMeanSpline","kaTOFCleanMean",scalingPointKa2,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,5);
  grScaledExWrtClean[6]  = ApplyScalingWrtSample(nSlice,"prTOFMeanSpline","prTOFCleanMean",scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,6);

  grScaledExWrtClean[7]  = ApplyScalingWrtSample(nSlice,"elSigmaSpline"   ,"elCleanSigma"   ,0.6            ,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,7);
  grScaledExWrtClean[8]  = ApplyScalingWrtSample(nSlice,"piSigmaSpline"   ,"piCleanSigma"   ,0.7            ,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,8);
  grScaledExWrtClean[9]  = ApplyScalingWrtSample(nSlice,"kaSigmaSpline"   ,"kaCleanSigma"   ,scalingPointKa2,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,9);
  grScaledExWrtClean[10] = ApplyScalingWrtSample(nSlice,"prSigmaSpline"   ,"prCleanSigma"   ,scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,10);
  grScaledExWrtClean[11] = ApplyScalingWrtSample(nSlice,"deSigmaSpline"   ,"deCleanSigma"   ,scalingPointDe ,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,11);
  grScaledExWrtClean[12] = ApplyScalingWrtSample(nSlice,"kaTOFSigmaSpline","kaTOFCleanSigma",scalingPointKa2,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,12);
  grScaledExWrtClean[13] = ApplyScalingWrtSample(nSlice,"prTOFSigmaSpline","prTOFCleanSigma",scalingPointPr2,GrArrFit,GrArrClean,treeClean,grScaledExWrtCleanSmooth,13);


  for (Int_t i=0; i<14; i++) {
    GrArrScaledExWrtClean->AddAt(grScaledExWrtClean[i],2*i);
    GrArrScaledExWrtClean->AddAt(grScaledExWrtCleanSmooth[i],2*i+1);
  }

  //   Make windows around mean of particles
  grWindowsAroundMean[0]  = ProduceWindowGraphs(nSlice, 1  ,"(elMeanSpline-elWindow):p",treeFit);
  grWindowsAroundMean[1]  = ProduceWindowGraphs(nSlice, 2  ,"(elMeanSpline+elWindow):p",treeFit);
  grWindowsAroundMean[2]  = ProduceWindowGraphs(nSlice, 3  ,"elMeanSpline:p"           ,treeFit);
  grWindowsAroundMean[3]  = ProduceWindowGraphs(nSlice, 4  ,"(piMeanSpline-piWindow):p",treeFit);
  grWindowsAroundMean[4]  = ProduceWindowGraphs(nSlice, 5  ,"(piMeanSpline+piWindow):p",treeFit);
  grWindowsAroundMean[5]  = ProduceWindowGraphs(nSlice, 6  ,"piMeanSpline:p"           ,treeFit);
  grWindowsAroundMean[6]  = ProduceWindowGraphs(nSlice, 7  ,"(kaMeanSpline-kaWindow):p",treeFit);
  grWindowsAroundMean[7]  = ProduceWindowGraphs(nSlice, 8  ,"(kaMeanSpline+kaWindow):p",treeFit);
  grWindowsAroundMean[8]  = ProduceWindowGraphs(nSlice, 9  ,"kaMeanSpline:p"           ,treeFit);
  grWindowsAroundMean[9]  = ProduceWindowGraphs(nSlice, 10 ,"(prMeanSpline-prWindow):p",treeFit);
  grWindowsAroundMean[10] = ProduceWindowGraphs(nSlice, 11 ,"(prMeanSpline+prWindow):p",treeFit);
  grWindowsAroundMean[11] = ProduceWindowGraphs(nSlice, 12 ,"prMeanSpline:p"           ,treeFit);
  grWindowsAroundMean[12] = ProduceWindowGraphs(nSlice, 13 ,"(deMeanSpline-deWindow):p",treeFit);
  grWindowsAroundMean[13] = ProduceWindowGraphs(nSlice, 14 ,"(deMeanSpline+deWindow):p",treeFit);
  grWindowsAroundMean[14] = ProduceWindowGraphs(nSlice, 15 ,"deMeanSpline:p"           ,treeFit);
  for (Int_t i=0; i<15; i++) GrArrWindows->AddAt(grWindowsAroundMean[i],i);


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
  TCanvas *canResults = GetFitResCanvas(GrArrFit,GrArrClean,GrArrScaledClean,GrArrScaledEx);
  canResults->SaveAs(Form("summary_iter%d_%d_eta%2.1f_cent%3.2f.png",iIter,fSign,fEtaDown,fCentDown));

  //   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cout << " ----- Dump results to File ----- " <<  endl;

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
  GrArrScaledExWrtClean -> Write("ScaledExWrtClean" ,TObject::kSingleKey);
  GrArrWindows       -> Write("Windows"     ,TObject::kSingleKey);
  GrArrOutlierSmooth -> Write("OulierSmooth",TObject::kSingleKey);

  //   delete arrays
  GrArrFit           -> Delete();
  GrArrClean         -> Delete();
  GrArrScaledEx      -> Delete();
  GrArrScaledClean   -> Delete();
  GrArrScaledExWrtClean -> Delete();
  GrArrOutlierSmooth -> Delete();
  delete GrArrFit;
  delete GrArrClean;
  delete GrArrScaledEx;
  delete GrArrScaledClean;
  delete GrArrScaledExWrtClean;
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

  cout << " ----- Dumping is done ----- " <<  endl;

}
// -------------------------------------------------------------------------------------------------------
void GetExandSampleParameters(Int_t islice, TFile *sFile, Double_t *arrMean, Double_t *arrMeanModif, Double_t *arrSigma, Double_t *cleanParams)
{

  //
  // Read Sample file
  //

  TTree *treeSample = (TTree*)sFile->Get("CleanSamples");

  Double_t elMean =0., elSigma = 0., elCleanMean = 0., elCleanSigma = 0.;
  Double_t piMean =0., piSigma = 0., piCleanMean = 0., piCleanSigma = 0.;
  Double_t kaMean =0., kaSigma = 0., kaCleanMean = 0., kaCleanSigma = 0.;
  Double_t prMean =0., prSigma = 0., prCleanMean = 0., prCleanSigma = 0.;
  Double_t deMean =0., deSigma = 0., deCleanMean = 0., deCleanSigma = 0.;
  Double_t kaTOFMean =0., kaTOFSigma = 0., kaTOFCleanMean = 0., kaTOFCleanSigma = 0.;
  Double_t prTOFMean =0., prTOFSigma = 0., prTOFCleanMean = 0., prTOFCleanSigma = 0.;


  Double_t elCleanSkew =0., elCleanKurtosis = 0., elCleanAmp = 0., elMeanMod=0.;
  Double_t piCleanSkew =0., piCleanKurtosis = 0., piCleanAmp = 0., piMeanMod=0.;
  Double_t kaCleanSkew =0., kaCleanKurtosis = 0., kaCleanAmp = 0., kaMeanMod=0.;
  Double_t prCleanSkew =0., prCleanKurtosis = 0., prCleanAmp = 0., prMeanMod=0.;
  Double_t deCleanSkew =0., deCleanKurtosis = 0., deCleanAmp = 0., deMeanMod=0.;
  Double_t kaTOFCleanSkew =0., kaTOFCleanKurtosis = 0., kaTOFCleanAmp = 0., kaTOFMeanMod=0.;
  Double_t prTOFCleanSkew =0., prTOFCleanKurtosis = 0., prTOFCleanAmp = 0., prTOFMeanMod=0.;


  // get Clean Mean and expected mean-sigma from Clean Tree
  treeSample->SetBranchAddress("elMeanSpline" ,&elMean);
  treeSample->SetBranchAddress("piMeanSpline" ,&piMean);
  treeSample->SetBranchAddress("kaMeanSpline" ,&kaMean);
  treeSample->SetBranchAddress("prMeanSpline" ,&prMean);
  treeSample->SetBranchAddress("deMeanSpline" ,&deMean);
  treeSample->SetBranchAddress("kaTOFMeanSpline" ,&kaTOFMean);
  treeSample->SetBranchAddress("prTOFMeanSpline" ,&prTOFMean);


  treeSample->SetBranchAddress("elMeanMod" ,&elMeanMod);
  treeSample->SetBranchAddress("piMeanMod" ,&piMeanMod);
  treeSample->SetBranchAddress("kaMeanMod" ,&kaMeanMod);
  treeSample->SetBranchAddress("prMeanMod" ,&prMeanMod);
  treeSample->SetBranchAddress("deMeanMod" ,&deMeanMod);
  treeSample->SetBranchAddress("kaTOFMeanMod" ,&kaTOFMeanMod);
  treeSample->SetBranchAddress("prTOFMeanMod" ,&prTOFMeanMod);

  treeSample->SetBranchAddress("elSigmaSpline" ,&elSigma);
  treeSample->SetBranchAddress("piSigmaSpline" ,&piSigma);
  treeSample->SetBranchAddress("kaSigmaSpline" ,&kaSigma);
  treeSample->SetBranchAddress("prSigmaSpline" ,&prSigma);
  treeSample->SetBranchAddress("deSigmaSpline" ,&deSigma);
  treeSample->SetBranchAddress("kaTOFSigmaSpline" ,&kaTOFSigma);
  treeSample->SetBranchAddress("prTOFSigmaSpline" ,&prTOFSigma);


  treeSample->SetBranchAddress("elCleanMean" ,&elCleanMean);
  treeSample->SetBranchAddress("piCleanMean" ,&piCleanMean);
  treeSample->SetBranchAddress("kaCleanMean" ,&kaCleanMean);
  treeSample->SetBranchAddress("prCleanMean" ,&prCleanMean);
  treeSample->SetBranchAddress("deCleanMean" ,&deCleanMean);
  treeSample->SetBranchAddress("kaTOFCleanMean" ,&kaTOFCleanMean);
  treeSample->SetBranchAddress("prTOFCleanMean" ,&prTOFCleanMean);


  treeSample->SetBranchAddress("elCleanSigma" ,&elCleanSigma);
  treeSample->SetBranchAddress("piCleanSigma" ,&piCleanSigma);
  treeSample->SetBranchAddress("kaCleanSigma" ,&kaCleanSigma);
  treeSample->SetBranchAddress("prCleanSigma" ,&prCleanSigma);
  treeSample->SetBranchAddress("deCleanSigma" ,&deCleanSigma);
  treeSample->SetBranchAddress("kaTOFCleanSigma" ,&kaTOFCleanSigma);
  treeSample->SetBranchAddress("prTOFCleanSigma" ,&prTOFCleanSigma);


  treeSample->SetBranchAddress("elCleanSkew" ,&elCleanSkew);
  treeSample->SetBranchAddress("piCleanSkew" ,&piCleanSkew);
  treeSample->SetBranchAddress("kaCleanSkew" ,&kaCleanSkew);
  treeSample->SetBranchAddress("prCleanSkew" ,&prCleanSkew);
  treeSample->SetBranchAddress("deCleanSkew" ,&deCleanSkew);
  treeSample->SetBranchAddress("kaTOFCleanSkew" ,&kaTOFCleanSkew);
  treeSample->SetBranchAddress("prTOFCleanSkew" ,&prTOFCleanSkew);


  treeSample->SetBranchAddress("elCleanKurtosis" ,&elCleanKurtosis);
  treeSample->SetBranchAddress("piCleanKurtosis" ,&piCleanKurtosis);
  treeSample->SetBranchAddress("kaCleanKurtosis" ,&kaCleanKurtosis);
  treeSample->SetBranchAddress("prCleanKurtosis" ,&prCleanKurtosis);
  treeSample->SetBranchAddress("deCleanKurtosis" ,&deCleanKurtosis);
  treeSample->SetBranchAddress("kaTOFCleanKurtosis" ,&kaTOFCleanKurtosis);
  treeSample->SetBranchAddress("prTOFCleanKurtosis" ,&prTOFCleanKurtosis);


  treeSample->SetBranchAddress("elCleanAmp" ,&elCleanAmp);
  treeSample->SetBranchAddress("piCleanAmp" ,&piCleanAmp);
  treeSample->SetBranchAddress("kaCleanAmp" ,&kaCleanAmp);
  treeSample->SetBranchAddress("prCleanAmp" ,&prCleanAmp);
  treeSample->SetBranchAddress("deCleanAmp" ,&deCleanAmp);
  treeSample->SetBranchAddress("kaTOFCleanAmp" ,&kaTOFCleanAmp);
  treeSample->SetBranchAddress("prTOFCleanAmp" ,&prTOFCleanAmp);


  treeSample->GetEntry(islice);
  cleanParams[0] = elCleanMean;
  cleanParams[1] = piCleanMean;
  cleanParams[2] = kaCleanMean;
  cleanParams[3] = prCleanMean;
  cleanParams[4] = deCleanMean;

  cleanParams[5] = elCleanSigma;
  cleanParams[6] = piCleanSigma;
  cleanParams[7] = kaCleanSigma;
  cleanParams[8] = prCleanSigma;
  cleanParams[9] = deCleanSigma;

  cleanParams[10] = elCleanSkew;
  cleanParams[11] = piCleanSkew;
  cleanParams[12] = kaCleanSkew;
  cleanParams[13] = prCleanSkew;
  cleanParams[14] = deCleanSkew;

  cleanParams[15] = elCleanKurtosis;
  cleanParams[16] = piCleanKurtosis;
  cleanParams[17] = kaCleanKurtosis;
  cleanParams[18] = prCleanKurtosis;
  cleanParams[19] = deCleanKurtosis;

  cleanParams[20] = elCleanAmp;
  cleanParams[21] = piCleanAmp;
  cleanParams[22] = kaCleanAmp;
  cleanParams[23] = prCleanAmp;
  cleanParams[24] = deCleanAmp;

  cleanParams[25] = kaTOFCleanMean;
  cleanParams[26] = kaTOFCleanSigma;
  cleanParams[27] = kaTOFCleanSkew;
  cleanParams[28] = kaTOFCleanKurtosis;
  cleanParams[29] = kaTOFCleanAmp;

  cleanParams[30] = prTOFCleanMean;
  cleanParams[31] = prTOFCleanSigma;
  cleanParams[32] = prTOFCleanSkew;
  cleanParams[33] = prTOFCleanKurtosis;
  cleanParams[34] = prTOFCleanAmp;

  arrMean[0] = elMean;     arrSigma[0] = elSigma;     arrMeanModif[0] = elMeanMod;
  arrMean[1] = piMean;     arrSigma[1] = piSigma;     arrMeanModif[1] = piMeanMod;
  arrMean[2] = kaMean;     arrSigma[2] = kaSigma;     arrMeanModif[2] = kaMeanMod;
  arrMean[3] = prMean;     arrSigma[3] = prSigma;     arrMeanModif[3] = prMeanMod;
  arrMean[4] = deMean;     arrSigma[4] = deSigma;     arrMeanModif[4] = deMeanMod;
  arrMean[5] = kaTOFMean;  arrSigma[5] = kaTOFSigma;  arrMeanModif[5] = kaTOFMeanMod;
  arrMean[6] = prTOFMean;  arrSigma[6] = prTOFSigma;  arrMeanModif[6] = prTOFMeanMod;


  // // If there is a missing PIDres point the take the next slice's parameters
  // for (Int_t i = islice; i<islice+10; i++) {
  //   treeSample->GetEntry(i);
  //   if (elMean!=0 && piMean!=0 && kaMean!=0 && prMean!=0 && deMean!=0 && elSigma!=0 && piSigma!=0 && kaSigma!=0 && prSigma!=0 && deSigma!=0){
  //     arrMean[0] = elMean;  arrSigma[0] = elSigma;
  //     arrMean[1] = piMean;  arrSigma[1] = piSigma;
  //     arrMean[2] = kaMean;  arrSigma[2] = kaSigma;
  //     arrMean[3] = prMean;  arrSigma[3] = prSigma;
  //     arrMean[4] = deMean;  arrSigma[4] = deSigma;
  //     break;
  //   }
  // }
  //
  // // If there is a missing PIDres point and next parameters does not help use previous ones
  // if (arrMean[0]==0||arrMean[1]==0||arrMean[2]==0||arrMean[3]==0||arrMean[4]==0||arrSigma[0]==0||arrSigma[1]==0||arrSigma[2]==0||arrMean[3]==0||arrMean[4]==0){
  //   for (Int_t i = islice; i>0; i--) {
  //     treeSample->GetEntry(i);
  //     if (arrMean[0]==0 || arrSigma[0]==0) { arrMean[0] = elMean;  arrSigma[0] = elSigma;}
  //     if (arrMean[1]==0 || arrSigma[1]==0) { arrMean[1] = piMean;  arrSigma[1] = piSigma;}
  //     if (arrMean[2]==0 || arrSigma[2]==0) { arrMean[2] = kaMean;  arrSigma[2] = kaSigma;}
  //     if (arrMean[3]==0 || arrSigma[3]==0) { arrMean[3] = prMean;  arrSigma[3] = prSigma;}
  //     if (arrMean[4]==0 || arrSigma[4]==0) { arrMean[4] = deMean;  arrSigma[4] = deSigma;}
  //     if (elMean!=0 && piMean!=0 && kaMean!=0 && prMean!=0 && elSigma!=0 && piSigma!=0 && kaSigma!=0 && prSigma!=0) break;
  //   }
  // }


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
  Double_t pMin = 0.2;
  Double_t pMax = 1.4;

  TH1D hCleanKurtosis[nParticle];
  TH1D hCleanSkewness[nParticle];
  TH1D hFitKurtosis[nParticle];
  TH1D hFitSkewness[nParticle];

  // arrays to create lookup tables
  Double_t axisPt[nBinsKSestimate];
  Double_t elLookUpFitSkew[nBinsKSestimate];     Double_t elLookUpCleanSkew[nBinsKSestimate];
  Double_t piLookUpFitSkew[nBinsKSestimate];      Double_t piLookUpCleanSkew[nBinsKSestimate];
  Double_t kaLookUpFitSkew[nBinsKSestimate];      Double_t kaLookUpCleanSkew[nBinsKSestimate];
  Double_t prLookUpFitSkew[nBinsKSestimate];      Double_t prLookUpCleanSkew[nBinsKSestimate];
  Double_t elLookUpFitKurtosis[nBinsKSestimate];  Double_t elLookUpCleanKurtosis[nBinsKSestimate];
  Double_t piLookUpFitKurtosis[nBinsKSestimate];  Double_t piLookUpCleanKurtosis[nBinsKSestimate];
  Double_t kaLookUpFitKurtosis[nBinsKSestimate];  Double_t kaLookUpCleanKurtosis[nBinsKSestimate];
  Double_t prLookUpFitKurtosis[nBinsKSestimate];  Double_t prLookUpCleanKurtosis[nBinsKSestimate];

  for (Int_t i = 0; i<nParticle; i++){
    hCleanKurtosis[i] = TH1D(Form("hCleanKurtosis_%s",particleStr[i].Data()),Form("hCleanKurtosis_%s",particleStr[i].Data()),100,1.5,3.5);
    hFitKurtosis[i]   = TH1D(Form("hFitKurtosis_%s"  ,particleStr[i].Data()),Form("hFitKurtosis_%s"  ,particleStr[i].Data()),100,1.5,3.5);
    hCleanSkewness[i] = TH1D(Form("hCleanSkewness_%s",particleStr[i].Data()),Form("hCleanSkewness_%s",particleStr[i].Data()),100,0.,2.);
    hFitSkewness[i]   = TH1D(Form("hFitSkewness_%s"  ,particleStr[i].Data()),Form("hFitSkewness_%s"  ,particleStr[i].Data()),100,0.,2.);
  }

  Int_t ptDownBin = 0;
  Int_t ptDownBinExpected = 0;
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
    ptDownBinExpected = (hExpectedMean[0]->GetXaxis()->FindBin(ptDown+0.001));                                // TODO
    ptDownBin = (h2D->GetXaxis()->FindBin(ptDown+0.001));                                // TODO
    ptUpBin   = (h2D->GetXaxis()->FindBin(ptUp+0.001));

    // initialise 1D histograms tobe fitted
    TH1D * h1CleanSamples[nParticle];
    TH1D * h1ExpectedMean[nParticle];
    TH1D * h1ExpectedSigma[nParticle];

    // Real particle parameters
    Double_t * arrMean       = new Double_t[10];
    Double_t * arrSigma      = new Double_t[10];
    Double_t * arrCleanSigma = new Double_t[10];
    Double_t * arrCleanMean  = new Double_t[10];
    Double_t * cleanParams   = new Double_t[40];

    // Clean Particle parameters
    Double_t * elParClean    = new Double_t[10];
    Double_t * piParClean    = new Double_t[10];
    Double_t * kaParClean    = new Double_t[10];
    Double_t * prParClean    = new Double_t[10];
    Double_t * deParClean    = new Double_t[10];
    Double_t * kaTOFParClean    = new Double_t[10];
    Double_t * prTOFParClean    = new Double_t[10];


    // Total fit parameters
    Double_t * pars          = new Double_t[31]; for (Int_t j=0; j<31; j++) pars[j] = 0;

    for (Int_t ipar=0; ipar<40; ipar++) cleanParams[ipar]  = 0;
    for (Int_t ipar=0; ipar<10; ipar++)
    {
      arrMean[ipar]       = 0;   piParClean[ipar]    = 0;
      arrSigma[ipar]      = 0;   kaParClean[ipar]    = 0;
      arrCleanMean[ipar]  = 0;   prParClean[ipar]    = 0;
      arrCleanSigma[ipar] = 0;   elParClean[ipar]    = 0;
    }

    // Get 1d slice
    TH1D * h1DKS = h2DKS->ProjectionY(Form("h1D_KS_%d",islice),ptDownBin,ptDownBin);
    h1DKS->Rebin(rebinFactor); h1DKS->Scale(1./rebinFactor);
    Double_t binWidth = h1DKS->GetXaxis()->GetBinWidth(50);
    h1DKS->Scale(1./binWidth);

    // Get Clean Samples 1d slice
    for (Int_t iParticle=0;iParticle<nParticle; iParticle++){
      // if (iParticle==2){   // TODO
      //   h1CleanSamples[iParticle] = GetClean1DSlice(hCleanSamples[iParticle], Form("h1D%s_slice_%f" ,particleStr[iParticle].Data()),Double_t(ptDownBinExpected));
      // } else {
      //   h1CleanSamples[iParticle] = GetClean1DSlice(hCleanSamples[iParticle], Form("h1D%s_slice_%f" ,particleStr[iParticle].Data()),Double_t(ptDownBin));
      // }
      h1CleanSamples[iParticle] = GetClean1DSlice(hCleanSamples[iParticle], Form("h1D%s_slice_%d" ,particleStr[iParticle].Data(), ptDownBin), ptDownBin);
      h1ExpectedMean[iParticle]  = GetClean1DSlice(hExpectedMean[iParticle], Form("h1D%s_ExMean_slice_%d",particleStr[iParticle].Data(), ptDownBin), ptDownBinExpected);
      h1ExpectedSigma[iParticle] = GetClean1DSlice(hExpectedSigma[iParticle],Form("h1D%s_ExSigma_slice_%d",particleStr[iParticle].Data(), ptDownBin),ptDownBinExpected);
    }

    // Get Expected Mean and Sigma
    for (Int_t iParticle=0;iParticle<nParticle; iParticle++){
      arrMean[iParticle]  = (exWRTMean) ? h1ExpectedMean[iParticle]->GetMean(): h1ExpectedMean[iParticle]->GetBinCenter(h1ExpectedMean[iParticle]->GetMaximumBin());
      arrSigma[iParticle] = (exWRTMean) ? h1ExpectedSigma[iParticle]->GetMean(): h1ExpectedSigma[iParticle]->GetBinCenter(h1ExpectedSigma[iParticle]->GetMaximumBin());
    }

    //  Make the fits
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[0] , elParClean, arrMean, arrSigma, kEl);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[1] , piParClean, arrMean, arrSigma, kPi);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[2] , kaParClean, arrMean, arrSigma, kKa);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[3] , prParClean, arrMean, arrSigma, kPr);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[4] , deParClean, arrMean, arrSigma, kDe);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[5] , kaTOFParClean, arrMean, arrSigma, kKaTOF);
    FitParticleSampleForKSEstimate(iks, h1CleanSamples[6] , prTOFParClean, arrMean, arrSigma, kPrTOF);
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

      axisPt[islice]=pt;
      elLookUpFitSkew[islice]=pars[4];        elLookUpCleanSkew[islice]=elParClean[4];
      piLookUpFitSkew[islice]=pars[9];        piLookUpCleanSkew[islice]=piParClean[4];
      kaLookUpFitSkew[islice]=pars[14];       kaLookUpCleanSkew[islice]=kaParClean[4];
      prLookUpFitSkew[islice]=pars[19];       prLookUpCleanSkew[islice]=prParClean[4];
      elLookUpFitKurtosis[islice]=pars[3];    elLookUpCleanKurtosis[islice]=elParClean[3];
      piLookUpFitKurtosis[islice]=pars[8];    piLookUpCleanKurtosis[islice]=piParClean[3];
      kaLookUpFitKurtosis[islice]=pars[13];   kaLookUpCleanKurtosis[islice]=kaParClean[3];
      prLookUpFitKurtosis[islice]=pars[18];   prLookUpCleanKurtosis[islice]=prParClean[3];

      // Safe region Kurtosis
      if (pt>=elCleanSafeMin && pt<=elCleanSafeMax) hCleanKurtosis[kEl].Fill(elParClean[3]);
      if (pt>=piCleanSafeMin && pt<=piCleanSafeMax) hCleanKurtosis[kPi].Fill(piParClean[3]);
      if (pt>=kaCleanSafeMin && pt<=kaCleanSafeMax) hCleanKurtosis[kKa].Fill(kaParClean[3]);
      if (pt>=prCleanSafeMin && pt<=prCleanSafeMax) hCleanKurtosis[kPr].Fill(prParClean[3]);

      if (pt>=elFitSafeMin && pt<=elFitSafeMax) hFitKurtosis[kEl].Fill(pars[3]);
      if (pt>=piFitSafeMin && pt<=piFitSafeMax) hFitKurtosis[kPi].Fill(pars[8]);
      if (pt>=kaFitSafeMin && pt<=kaFitSafeMax) hFitKurtosis[kKa].Fill(pars[13]);
      if (pt>=prFitSafeMin && pt<=prFitSafeMax) hFitKurtosis[kPr].Fill(pars[18]);

    } else {

      // Safe region Skewness
      if (pt>=elCleanSafeMin && pt<=elCleanSafeMax) hCleanSkewness[kEl].Fill(elParClean[4]);
      if (pt>=piCleanSafeMin && pt<=piCleanSafeMax) hCleanSkewness[kPi].Fill(piParClean[4]);
      if (pt>=kaCleanSafeMin && pt<=kaCleanSafeMax) hCleanSkewness[kKa].Fill(kaParClean[4]);
      if (pt>=prCleanSafeMin && pt<=prCleanSafeMax) hCleanSkewness[kPr].Fill(prParClean[4]);

      if (pt>=elFitSafeMin && pt<=elFitSafeMax) hFitSkewness[kEl].Fill(pars[4]);
      if (pt>=piFitSafeMin && pt<=piFitSafeMax) hFitSkewness[kPi].Fill(pars[9]);
      if (pt>=kaFitSafeMin && pt<=kaFitSafeMax) hFitSkewness[kKa].Fill(pars[14]);
      if (pt>=prFitSafeMin && pt<=prFitSafeMax) hFitSkewness[kPr].Fill(pars[19]);

    }

    // dump all info to tree
    debugFile -> GetFile()->cd();
    *debugFile << "ks" <<

    "iks="          << iks                <<
    "p="            << pt                 <<
    "eta="          << fEtaDown             <<
    "cent="         << fCentDown            <<

    "elMeanSpline="       << arrMean[0]         <<
    "piMeanSpline="       << arrMean[1]         <<
    "kaMeanSpline="       << arrMean[2]         <<
    "prMeanSpline="       << arrMean[3]         <<
    "deMeanSpline="       << arrMean[4]         <<
    "kaTOFMeanSpline="    << arrMean[5]         <<
    "prTOFMeanSpline="    << arrMean[6]         <<


    "elSigmaSpline="      << arrSigma[0]        <<
    "piSigmaSpline="      << arrSigma[1]        <<
    "kaSigmaSpline="      << arrSigma[2]        <<
    "prSigmaSpline="      << arrSigma[3]        <<
    "deSigmaSpline="      << arrSigma[4]        <<
    "kaTOFSigmaSpline="   << arrSigma[5]        <<
    "prTOFSigmaSpline="   << arrSigma[6]        ;

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
    "prCleanChi2="     << prParClean[6]    <<

    "kaTOFCleanAmp="      << kaTOFParClean[0]    <<
    "kaTOFCleanMean="     << kaTOFParClean[1]    <<
    "kaTOFCleanSigma="    << kaTOFParClean[2]    <<
    "kaTOFCleanKurtosis=" << kaTOFParClean[3]    <<
    "kaTOFCleanSkew="     << kaTOFParClean[4]    <<
    "kaTOFCleanInt="      << kaTOFParClean[5]    <<
    "kaTOFCleanChi2="     << kaTOFParClean[6]    <<

    "prTOFCleanAmp="      << prTOFParClean[0]    <<
    "prTOFCleanMean="     << prTOFParClean[1]    <<
    "prTOFCleanSigma="    << prTOFParClean[2]    <<
    "prTOFCleanKurtosis=" << prTOFParClean[3]    <<
    "prTOFCleanSkew="     << prTOFParClean[4]    <<
    "prTOFCleanInt="      << prTOFParClean[5]    <<
    "prTOFCleanChi2="     << prTOFParClean[6]    ;

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

    "totChi2="         << pars[25]       <<

    "elFitMax="        << pars[26]       <<
    "piFitMax="        << pars[27]       <<
    "kaFitMax="        << pars[28]       <<
    "prFitMax="        << pars[29]       <<
    "deFitMax="        << pars[30]       <<

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
    delete [] deParClean;
    delete [] kaTOFParClean;
    delete [] prTOFParClean;
    delete [] pars;

  }



  // Set the kurtosis and skewness to be used in further iterations
  if (iks==0){

    kurtosisFixClean[kEl] = hCleanKurtosis[kEl].GetMean();
    kurtosisFixClean[kPi] = hCleanKurtosis[kPi].GetMean();
    kurtosisFixClean[kKa] = hCleanKurtosis[kKa].GetMean();
    kurtosisFixClean[kPr] = hCleanKurtosis[kPr].GetMean();

    kurtosisFixFit[kEl]   = hFitKurtosis[kEl].GetMean();
    kurtosisFixFit[kPi]   = hFitKurtosis[kPi].GetMean();
    kurtosisFixFit[kKa]   = hFitKurtosis[kKa].GetMean();
    kurtosisFixFit[kPr]   = hFitKurtosis[kPr].GetMean();

  }

  // dump final fixed KS into debug file
  debugFile -> GetFile()->cd();
  *debugFile << "fixedKS" <<

  "iks="             << iks                     <<
  "eta="             << fEtaDown                <<
  "cent="            << fCentDown               <<

  "FitSkewEl="       << skewnessFixFit[kEl]     <<
  "FitSkewPi="       << skewnessFixFit[kPi]     <<
  "FitSkewKa="       << skewnessFixFit[kKa]     <<
  "FitSkewPr="       << skewnessFixFit[kPr]     <<

  "CleanSkewEl="     << skewnessFixClean[kEl]   <<
  "CleanSkewPi="     << skewnessFixClean[kPi]   <<
  "CleanSkewKa="     << skewnessFixClean[kKa]   <<
  "CleanSkewPr="     << skewnessFixClean[kPr]   <<

  "FitKurtEl="       << kurtosisFixFit[kEl]     <<
  "FitKurtPi="       << kurtosisFixFit[kPi]     <<
  "FitKurtKa="       << kurtosisFixFit[kKa]     <<
  "FitKurtPr="       << kurtosisFixFit[kPr]     <<

  "CleanKurtEl="     << kurtosisFixClean[kEl]   <<
  "CleanKurtPi="     << kurtosisFixClean[kPi]   <<
  "CleanKurtKa="     << kurtosisFixClean[kKa]   <<
  "CleanKurtPr="     << kurtosisFixClean[kPr]   <<

  "\n";

  debugFile -> GetFile()->cd();
  for (Int_t i = 1; i<4; i++){
    if (iks==0) hFitKurtosis[i]  .Write(Form("Fit_Kurtosis_%s"  ,particleStr[i].Data()));
    if (iks==1) hFitSkewness[i]  .Write(Form("Fit_Skewness_%s"  ,particleStr[i].Data()));
    if (iks==0) hCleanKurtosis[i].Write(Form("Clean_Kurtosis_%s",particleStr[i].Data()));
    if (iks==1) hCleanSkewness[i].Write(Form("Clean_Skewness_%s",particleStr[i].Data()));
  }

  if (iks==1) {

    // decide whether to let automaticly found KS or by eye
    if (automaticKS) {
      skewnessFixFit[kEl]   = hFitSkewness[kEl].GetMean()*piSkewnessScan;  // multiplied with piSkewnessScan since skewness decrease with p
      skewnessFixFit[kPi]   = hFitSkewness[kPi].GetMean()*piSkewnessScan;  // multiplied with piSkewnessScan since skewness decrease with p
      skewnessFixFit[kKa]   = hFitSkewness[kKa].GetMean()*kaSkewnessScan;
      skewnessFixFit[kPr]   = hFitSkewness[kPr].GetMean()*prSkewnessScan;

      skewnessFixClean[kEl] = hCleanSkewness[kEl].GetMean()*piSkewnessScan;
      skewnessFixClean[kPi] = hCleanSkewness[kPi].GetMean()*piSkewnessScan;
      skewnessFixClean[kKa] = hCleanSkewness[kKa].GetMean()*kaSkewnessScan;
      skewnessFixClean[kPr] = hCleanSkewness[kPr].GetMean()*prSkewnessScan;
    } else {
      skewnessFixFit[kEl]   = TMath::Min(hFitSkewness[kEl].GetMean(),piSkewnessScan);
      skewnessFixFit[kPi]   = TMath::Min(hFitSkewness[kPi].GetMean(),piSkewnessScan);
      skewnessFixFit[kKa]   = TMath::Min(hFitSkewness[kKa].GetMean(),kaSkewnessScan);
      skewnessFixFit[kPr]   = TMath::Min(hFitSkewness[kPr].GetMean(),prSkewnessScan);

      skewnessFixClean[kEl] = TMath::Min(hCleanSkewness[kEl].GetMean(),piSkewnessScan);
      skewnessFixClean[kPi] = TMath::Min(hCleanSkewness[kPi].GetMean(),piSkewnessScan);
      skewnessFixClean[kKa] = TMath::Min(hCleanSkewness[kKa].GetMean(),kaSkewnessScan);
      skewnessFixClean[kPr] = TMath::Min(hCleanSkewness[kPr].GetMean(),prSkewnessScan);
    }

  }

  // Set the final kurtosis and skewness params --> take the estimation form clean samples
  if (iks==1){

    if (useSafeRegion){
      skewnessFix[kEl] = skewnessFixClean[kEl];
      skewnessFix[kPi] = skewnessFixFit[kPi];
      skewnessFix[kKa] = skewnessFixClean[kKa];
      skewnessFix[kPr] = skewnessFixFit[kPr];
      kurtosisFix[kEl] = kurtosisFixClean[kEl]*piKurtosisScan;
      kurtosisFix[kPi] = kurtosisFixFit[kPi]*piKurtosisScan;
      kurtosisFix[kKa] = kurtosisFixClean[kKa]*kaKurtosisScan;
      kurtosisFix[kPr] = kurtosisFixFit[kPr]*prKurtosisScan;
    } else {
      skewnessFix[kEl] = skewnessFixClean[kEl];
      skewnessFix[kPi] = skewnessFixClean[kPi];
      skewnessFix[kKa] = skewnessFixClean[kKa];
      skewnessFix[kPr] = skewnessFixClean[kPr];
      kurtosisFix[kEl] = kurtosisFixClean[kEl];
      kurtosisFix[kPi] = kurtosisFixClean[kPi];
      kurtosisFix[kKa] = kurtosisFixClean[kKa];
      kurtosisFix[kPr] = kurtosisFixClean[kPr];
    }

  }

  if (iks==0){

    fitLookUpCleanSkew[kEl] = new TF1(Form("%s_fCleanSkew",particleStr[kEl].Data()),"pol1",0.25,0.6);
    fitLookUpCleanSkew[kPi] = new TF1(Form("%s_fCleanSkew",particleStr[kPi].Data()),"pol1",0.4,1.4);
    fitLookUpCleanSkew[kKa] = new TF1(Form("%s_fCleanSkew",particleStr[kKa].Data()),"pol1",0.5,1.);
    fitLookUpCleanSkew[kPr] = new TF1(Form("%s_fCleanSkew",particleStr[kPr].Data()),"pol1",0.7,1.1);

    fitLookUpFitSkew[kEl] = new TF1(Form("%s_fFitSkew",particleStr[kEl].Data()),"pol1",0.2,0.4);
    fitLookUpFitSkew[kPi] = new TF1(Form("%s_fFitSkew",particleStr[kPi].Data()),"pol1",0.41,0.6);
    fitLookUpFitSkew[kKa] = new TF1(Form("%s_fFitSkew",particleStr[kKa].Data()),"pol1",0.2,0.4);
    fitLookUpFitSkew[kPr] = new TF1(Form("%s_fFitSkew",particleStr[kPr].Data()),"pol1",0.4,0.75);

    grLookUpCleanSkew[kEl] = new TGraphErrors(nBinsKSestimate,axisPt,elLookUpCleanSkew,0,0);
    grLookUpCleanSkew[kPi] = new TGraphErrors(nBinsKSestimate,axisPt,piLookUpCleanSkew,0,0);
    grLookUpCleanSkew[kKa] = new TGraphErrors(nBinsKSestimate,axisPt,kaLookUpCleanSkew,0,0);
    grLookUpCleanSkew[kPr] = new TGraphErrors(nBinsKSestimate,axisPt,prLookUpCleanSkew,0,0);
    //
    grLookUpFitSkew[kEl] = new TGraphErrors(nBinsKSestimate,axisPt,elLookUpFitSkew,0,0);
    grLookUpFitSkew[kPi] = new TGraphErrors(nBinsKSestimate,axisPt,piLookUpFitSkew,0,0);
    grLookUpFitSkew[kKa] = new TGraphErrors(nBinsKSestimate,axisPt,kaLookUpFitSkew,0,0);
    grLookUpFitSkew[kPr] = new TGraphErrors(nBinsKSestimate,axisPt,prLookUpFitSkew,0,0);
    //
    grLookUpFitKurtosis[kEl] = new TGraphErrors(nBinsKSestimate,axisPt,elLookUpFitKurtosis,0,0);
    grLookUpFitKurtosis[kPi] = new TGraphErrors(nBinsKSestimate,axisPt,piLookUpFitKurtosis,0,0);
    grLookUpFitKurtosis[kKa] = new TGraphErrors(nBinsKSestimate,axisPt,kaLookUpFitKurtosis,0,0);
    grLookUpFitKurtosis[kPr] = new TGraphErrors(nBinsKSestimate,axisPt,prLookUpFitKurtosis,0,0);
    //
    grLookUpCleanKurtosis[kEl] = new TGraphErrors(nBinsKSestimate,axisPt,elLookUpCleanKurtosis,0,0);
    grLookUpCleanKurtosis[kPi] = new TGraphErrors(nBinsKSestimate,axisPt,piLookUpCleanKurtosis,0,0);
    grLookUpCleanKurtosis[kKa] = new TGraphErrors(nBinsKSestimate,axisPt,kaLookUpCleanKurtosis,0,0);
    grLookUpCleanKurtosis[kPr] = new TGraphErrors(nBinsKSestimate,axisPt,prLookUpCleanKurtosis,0,0);
    //
    // Apply Fits
    grLookUpCleanSkew[kEl]->Fit(fitLookUpCleanSkew[kEl],"QMR+");
    grLookUpCleanSkew[kPi]->Fit(fitLookUpCleanSkew[kPi],"QMR+");
    grLookUpCleanSkew[kKa]->Fit(fitLookUpCleanSkew[kKa],"QMR+");
    grLookUpCleanSkew[kPr]->Fit(fitLookUpCleanSkew[kPr],"QMR+");
    grLookUpFitSkew[kEl]->Fit(fitLookUpFitSkew[kEl],"QMR+");
    grLookUpFitSkew[kPi]->Fit(fitLookUpFitSkew[kPi],"QMR+");
    grLookUpFitSkew[kKa]->Fit(fitLookUpFitSkew[kKa],"QMR+");
    grLookUpFitSkew[kPr]->Fit(fitLookUpFitSkew[kPr],"QMR+");

    debugFile -> GetFile()->cd();
    for (Int_t i=0;i<4;i++){
      grLookUpCleanSkew[i]    ->SetMarkerColor(kBlack);  grLookUpCleanSkew[i]    ->SetLineColor(kBlack);   grLookUpCleanSkew[i]    ->SetMarkerStyle(20);
      grLookUpFitSkew[i]      ->SetMarkerColor(kRed);    grLookUpFitSkew[i]      ->SetLineColor(kRed);     grLookUpFitSkew[i]      ->SetMarkerStyle(24);
      grLookUpCleanKurtosis[i]->SetMarkerColor(kBlack);  grLookUpCleanKurtosis[i]->SetLineColor(kBlack);   grLookUpCleanKurtosis[i]->SetMarkerStyle(20);
      grLookUpFitKurtosis[i]  ->SetMarkerColor(kRed);    grLookUpFitKurtosis[i]  ->SetLineColor(kRed);     grLookUpFitKurtosis[i]  ->SetMarkerStyle(24);
      //
      grLookUpCleanSkew[i]    ->Write(Form("%s_grCleanSkew",particleStr[i].Data()));
      grLookUpFitSkew[i]      ->Write(Form("%s_grFitSkew",particleStr[i].Data()));
      grLookUpCleanKurtosis[i]->Write(Form("%s_grCleanKurtosis",particleStr[i].Data()));
      grLookUpFitKurtosis[i]  ->Write(Form("%s_grFitKurtosis",particleStr[i].Data()));

      // fitLookUpCleanSkew[i]->SetRange(0.,3.);
      // fitLookUpFitSkew[i]->SetRange(0.,3.);
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
  Int_t ptDownBinExpected = 0;
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
    ptDownBinExpected = (hExpectedMean[0]->GetXaxis()->FindBin(ptDown+0.001));                                // TODO
    ptDownBin = (h2D->GetXaxis()->FindBin(ptDown+0.001));                                // TODO
    ptUpBin   = (h2D->GetXaxis()->FindBin(ptUp+0.001));

    // Storage for the histograms
    TObjArray * cleanSampArr = new TObjArray(10);  cleanSampArr -> SetOwner(kTRUE);
    cleanSampArr->SetName(Form("cleanSamples_%d_%4.3f",islice,pt));

    TObjArray * cleanSampArrFreeKS = new TObjArray(10);  cleanSampArrFreeKS -> SetOwner(kTRUE);
    cleanSampArrFreeKS->SetName(Form("cleanSamplesFreeKS_%d_%4.3f",islice,pt));

    TH1D * h1CleanSamples[nParticle];
    TH1D * h1CleanSamplesFreeKS[nParticle];
    TH1D * h1ExpectedMean[nParticle];
    TH1D * h1ExpectedSigma[nParticle];

    // Real particle parameters
    Double_t * arrSigma      = new Double_t[10];
    Double_t * arrMean       = new Double_t[10];
    Double_t * arrMeanModif  = new Double_t[10];
    Double_t * arrCleanSigma = new Double_t[10];
    Double_t * arrCleanMean  = new Double_t[10];
    Double_t * arrCleanMax   = new Double_t[10];

    // Clean Particle parameters
    Double_t * elParClean    = new Double_t[10];
    Double_t * piParClean    = new Double_t[10];
    Double_t * kaParClean    = new Double_t[10];
    Double_t * prParClean    = new Double_t[10];
    Double_t * deParClean    = new Double_t[10];
    Double_t * kaTOFParClean    = new Double_t[10];
    Double_t * prTOFParClean    = new Double_t[10];

    Double_t * elParCleanFreeKS    = new Double_t[10];
    Double_t * piParCleanFreeKS    = new Double_t[10];
    Double_t * kaParCleanFreeKS    = new Double_t[10];
    Double_t * prParCleanFreeKS    = new Double_t[10];
    Double_t * deParCleanFreeKS    = new Double_t[10];
    Double_t * kaTOFParCleanFreeKS    = new Double_t[10];
    Double_t * prTOFParCleanFreeKS    = new Double_t[10];

    for (Int_t ipar=0; ipar<10; ipar++)
    {
      arrMean[ipar]      = 0;   elParClean[ipar] = 0;   piParCleanFreeKS[ipar] = 0;
      arrSigma[ipar]     = 0;   piParClean[ipar] = 0;   kaParCleanFreeKS[ipar] = 0;
      arrCleanMean[ipar] = 0;   kaParClean[ipar] = 0;   prParCleanFreeKS[ipar] = 0;
      arrCleanMax[ipar]  = 0;   prParClean[ipar] = 0;   elParCleanFreeKS[ipar] = 0;
      arrCleanSigma[ipar]= 0;   deParClean[ipar] = 0;   deParCleanFreeKS[ipar] = 0;
      arrMeanModif[ipar] = 0;
      //
      kaTOFParClean[ipar] = 0;
      prTOFParClean[ipar] = 0;
      kaTOFParCleanFreeKS[ipar] = 0;
      prTOFParCleanFreeKS[ipar] = 0;
    }

    // Get Clean Samples 1d slice
    for (Int_t iParticle=0;iParticle<nParticle; iParticle++){
      // if (iParticle==2 || iParticle==5){   // TODO
      //   h1CleanSamples[iParticle] = GetClean1DSlice(hCleanSamples[iParticle], Form("h1D%s_slice_%f" ,particleStr[iParticle].Data()),Double_t(ptDownBinExpected));
      // } else {
      //   h1CleanSamples[iParticle] = GetClean1DSlice(hCleanSamples[iParticle], Form("h1D%s_slice_%f" ,particleStr[iParticle].Data()),Double_t(ptDownBin));
      // }
      h1CleanSamples[iParticle] = GetClean1DSlice(hCleanSamples[iParticle], Form("h1D%s_slice_%d" ,particleStr[iParticle].Data(), ptDownBin), ptDownBin);
      h1ExpectedMean[iParticle]  = GetClean1DSlice(hExpectedMean[iParticle], Form("h1D%s_ExMean_slice_%d",particleStr[iParticle].Data(), ptDownBin), ptDownBinExpected);
      h1ExpectedSigma[iParticle] = GetClean1DSlice(hExpectedSigma[iParticle],Form("h1D%s_ExSigma_slice_%d",particleStr[iParticle].Data(), ptDownBin), ptDownBinExpected);
      h1CleanSamplesFreeKS[iParticle] = (TH1D*)h1CleanSamples[iParticle]->Clone();
    }

    // Get Expected Mean and Sigma
    for (Int_t iParticle=0;iParticle<nParticle; iParticle++){
      arrMean[iParticle]  = (exWRTMean) ? h1ExpectedMean[iParticle]->GetMean():  h1ExpectedMean[iParticle] ->GetBinCenter(h1ExpectedMean[iParticle]->GetMaximumBin());
      arrSigma[iParticle] = (exWRTMean) ? h1ExpectedSigma[iParticle]->GetMean(): h1ExpectedSigma[iParticle]->GetBinCenter(h1ExpectedSigma[iParticle]->GetMaximumBin());
    }
    //
    // Modify mean positions due to change in the skewness and kurtosis
    ModifyExpectedMeanAndApplyVariableSkewness(arrMean, arrSigma, arrMeanModif, pt);

    //  Make the fits
    TMatrixD *elCleanMatrix    = FitParticleSample(h1CleanSamples[0], elParClean, arrMean, arrSigma, kEl, pt);
    TMatrixD *piCleanMatrix    = FitParticleSample(h1CleanSamples[1], piParClean, arrMean, arrSigma, kPi, pt);
    TMatrixD *kaCleanMatrix    = FitParticleSample(h1CleanSamples[2], kaParClean, arrMean, arrSigma, kKa, pt);
    TMatrixD *prCleanMatrix    = FitParticleSample(h1CleanSamples[3], prParClean, arrMean, arrSigma, kPr, pt);
    TMatrixD *deCleanMatrix    = FitParticleSample(h1CleanSamples[4], deParClean, arrMean, arrSigma, kDe, pt);
    TMatrixD *kaTOFCleanMatrix = FitParticleSample(h1CleanSamples[5], kaTOFParClean, arrMean, arrSigma, kKaTOF, pt);
    TMatrixD *prTOFCleanMatrix = FitParticleSample(h1CleanSamples[6], prTOFParClean, arrMean, arrSigma, kPrTOF, pt);


    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[0], elParCleanFreeKS, arrMean, arrSigma, kEl, pt);
    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[1], piParCleanFreeKS, arrMean, arrSigma, kPi, pt);
    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[2], kaParCleanFreeKS, arrMean, arrSigma, kKa, pt);
    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[3], prParCleanFreeKS, arrMean, arrSigma, kPr, pt);
    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[4], deParCleanFreeKS, arrMean, arrSigma, kDe, pt);
    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[5], kaTOFParCleanFreeKS, arrMean, arrSigma, kKaTOF, pt);
    FitParticleSampleFreeKS(h1CleanSamplesFreeKS[6], prTOFParCleanFreeKS, arrMean, arrSigma, kPrTOF, pt);


    GetCleanExParams(h1CleanSamples[0],arrSigma,arrCleanSigma,arrCleanMean,arrCleanMax,kEl);
    GetCleanExParams(h1CleanSamples[1],arrSigma,arrCleanSigma,arrCleanMean,arrCleanMax,kPi);
    GetCleanExParams(h1CleanSamples[2],arrSigma,arrCleanSigma,arrCleanMean,arrCleanMax,kKa);
    GetCleanExParams(h1CleanSamples[3],arrSigma,arrCleanSigma,arrCleanMean,arrCleanMax,kPr);
    GetCleanExParams(h1CleanSamples[4],arrSigma,arrCleanSigma,arrCleanMean,arrCleanMax,kDe);
    GetCleanExParams(h1CleanSamples[5],arrSigma,arrCleanSigma,arrCleanMean,arrCleanMax,kKaTOF);
    GetCleanExParams(h1CleanSamples[6],arrSigma,arrCleanSigma,arrCleanMean,arrCleanMax,kPrTOF);


    // Dump histograms into arrays
    h1CleanSamples[0] -> SetName(Form("CleanElectron_%d_%4.3f" ,islice,pt));
    h1CleanSamples[1] -> SetName(Form("CleanPion_%d_%4.3f"     ,islice,pt));
    h1CleanSamples[2] -> SetName(Form("CleanKaon_%d_%4.3f"     ,islice,pt));
    h1CleanSamples[3] -> SetName(Form("CleanProton_%d_%4.3f"   ,islice,pt));
    h1CleanSamples[4] -> SetName(Form("CleanDeuteron_%d_%4.3f" ,islice,pt));
    h1CleanSamples[5] -> SetName(Form("CleanKaonTOF_%d_%4.3f"  ,islice,pt));
    h1CleanSamples[6] -> SetName(Form("CleanProtonTOF_%d_%4.3f",islice,pt));


    h1CleanSamplesFreeKS[0] -> SetName(Form("CleanElectronFreeKS_%d_%4.3f" ,islice,pt));
    h1CleanSamplesFreeKS[1] -> SetName(Form("CleanPionFreeKS_%d_%4.3f"     ,islice,pt));
    h1CleanSamplesFreeKS[2] -> SetName(Form("CleanKaonFreeKS_%d_%4.3f"     ,islice,pt));
    h1CleanSamplesFreeKS[3] -> SetName(Form("CleanProtonFreeKS_%d_%4.3f"   ,islice,pt));
    h1CleanSamplesFreeKS[4] -> SetName(Form("CleanDeuteronFreeKS_%d_%4.3f" ,islice,pt));
    h1CleanSamplesFreeKS[5] -> SetName(Form("CleanKaonTOFFreeKS_%d_%4.3f" ,islice,pt));
    h1CleanSamplesFreeKS[6] -> SetName(Form("CleanProtonTOFFreeKS_%d_%4.3f" ,islice,pt));


    cleanSampArr -> AddAt(h1CleanSamples[0],0);
    cleanSampArr -> AddAt(h1CleanSamples[1],1);
    cleanSampArr -> AddAt(h1CleanSamples[2],2);
    cleanSampArr -> AddAt(h1CleanSamples[3],3);
    cleanSampArr -> AddAt(h1CleanSamples[4],4);
    cleanSampArr -> AddAt(h1CleanSamples[5],5);
    cleanSampArr -> AddAt(h1CleanSamples[6],6);


    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[0],0);
    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[1],1);
    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[2],2);
    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[3],3);
    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[4],4);
    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[5],5);
    cleanSampArrFreeKS -> AddAt(h1CleanSamplesFreeKS[6],6);


    cleanResArr      . AddAt(cleanSampArr,islice);
    cleanResArrFreeKS. AddAt(cleanSampArrFreeKS,islice);

    if (corrStr=="") dEdxCorr = 0;
    if (corrStr=="Corr") dEdxCorr = 1;

    // Clean Samples dump tree
    sampleFile->GetFile()->cd();
    *sampleFile << "CleanSamples" <<

    "corr="      << dEdxCorr                <<
    "slice="        << islice             <<
    "p="            << pt                 <<
    "eta="          << fEtaDown           <<
    "cent="         << fCentDown          <<

    "elMeanMod="    << arrMeanModif[0]         <<
    "piMeanMod="    << arrMeanModif[1]         <<
    "kaMeanMod="    << arrMeanModif[2]         <<
    "prMeanMod="    << arrMeanModif[3]         <<
    "deMeanMod="    << arrMeanModif[4]         <<
    "kaTOFMeanMod=" << arrMeanModif[5]         <<
    "prTOFMeanMod=" << arrMeanModif[6]         <<


    "elMeanSpline="       << arrMean[0]         <<
    "piMeanSpline="       << arrMean[1]         <<
    "kaMeanSpline="       << arrMean[2]         <<
    "prMeanSpline="       << arrMean[3]         <<
    "deMeanSpline="       << arrMean[4]         <<
    "kaTOFMeanSpline="    << arrMean[5]         <<
    "prTOFMeanSpline="    << arrMean[6]         <<


    "elSigmaSpline="      << arrSigma[0]        <<
    "piSigmaSpline="      << arrSigma[1]        <<
    "kaSigmaSpline="      << arrSigma[2]        <<
    "prSigmaSpline="      << arrSigma[3]        <<
    "deSigmaSpline="      << arrSigma[4]        <<
    "kaTOFSigmaSpline="   << arrSigma[5]        <<
    "prTOFSigmaSpline="   << arrSigma[6]        ;


    *sampleFile << "CleanSamples" <<

    "elExCMean="    << arrCleanMean[0]    <<
    "piExCMean="    << arrCleanMean[1]    <<
    "kaExCMean="    << arrCleanMean[2]    <<
    "prExCMean="    << arrCleanMean[3]    <<
    "deExCMean="    << arrCleanMean[4]    <<
    "kaTOFExCMean=" << arrCleanMean[5]    <<
    "prTOFExCMean=" << arrCleanMean[6]    <<



    "elExCSigma="   << arrCleanSigma[0]   <<
    "piExCSigma="   << arrCleanSigma[1]   <<
    "kaExCSigma="   << arrCleanSigma[2]   <<
    "prExCSigma="   << arrCleanSigma[3]   <<
    "deExCSigma="   << arrCleanSigma[4]   <<
    "kaTOFExCSigma="<< arrCleanSigma[5]   <<
    "prTOFExCSigma="<< arrCleanSigma[6]   ;


    *sampleFile << "CleanSamples" <<

    "elCleanAmp="      << elParClean[0]    <<
    "elCleanMean="     << elParClean[1]    <<
    "elCleanSigma="    << elParClean[2]    <<
    "elCleanKurtosis=" << elParClean[3]    <<
    "elCleanSkew="     << elParClean[4]    <<
    "elCleanInt="      << elParClean[5]    <<
    "elCleanChi2="     << elParClean[6]    <<
    "elCleanMax="      << elParClean[7]    <<
    "elCleanMatrix.="  << elCleanMatrix    <<

    "piCleanAmp="      << piParClean[0]    <<
    "piCleanMean="     << piParClean[1]    <<
    "piCleanSigma="    << piParClean[2]    <<
    "piCleanKurtosis=" << piParClean[3]    <<
    "piCleanSkew="     << piParClean[4]    <<
    "piCleanInt="      << piParClean[5]    <<
    "piCleanChi2="     << piParClean[6]    <<
    "piCleanMax="      << piParClean[7]    <<
    "piCleanMatrix.="  << piCleanMatrix    <<

    "kaCleanAmp="      << kaParClean[0]    <<
    "kaCleanMean="     << kaParClean[1]    <<
    "kaCleanSigma="    << kaParClean[2]    <<
    "kaCleanKurtosis=" << kaParClean[3]    <<
    "kaCleanSkew="     << kaParClean[4]    <<
    "kaCleanInt="      << kaParClean[5]    <<
    "kaCleanChi2="     << kaParClean[6]    <<
    "kaCleanMax="      << kaParClean[7]    <<
    "kaCleanMatrix.="  << kaCleanMatrix    ;

    *sampleFile << "CleanSamples" <<

    "prCleanAmp="      << prParClean[0]    <<
    "prCleanMean="     << prParClean[1]    <<
    "prCleanSigma="    << prParClean[2]    <<
    "prCleanKurtosis=" << prParClean[3]    <<
    "prCleanSkew="     << prParClean[4]    <<
    "prCleanInt="      << prParClean[5]    <<
    "prCleanChi2="     << prParClean[6]    <<
    "prCleanMax="      << prParClean[7]    <<
    "prCleanMatrix.="  << prCleanMatrix    <<

    "deCleanAmp="      << deParClean[0]    <<
    "deCleanMean="     << deParClean[1]    <<
    "deCleanSigma="    << deParClean[2]    <<
    "deCleanKurtosis=" << deParClean[3]    <<
    "deCleanSkew="     << deParClean[4]    <<
    "deCleanInt="      << deParClean[5]    <<
    "deCleanChi2="     << deParClean[6]    <<
    "deCleanMax="      << deParClean[7]    <<
    "deCleanMatrix.="  << deCleanMatrix    <<

    "kaTOFCleanAmp="      << kaTOFParClean[0]    <<
    "kaTOFCleanMean="     << kaTOFParClean[1]    <<
    "kaTOFCleanSigma="    << kaTOFParClean[2]    <<
    "kaTOFCleanKurtosis=" << kaTOFParClean[3]    <<
    "kaTOFCleanSkew="     << kaTOFParClean[4]    <<
    "kaTOFCleanInt="      << kaTOFParClean[5]    <<
    "kaTOFCleanChi2="     << kaTOFParClean[6]    <<
    "kaTOFCleanMax="      << kaTOFParClean[7]    <<
    "kaTOFCleanMatrix.="  << kaTOFCleanMatrix    <<

    "prTOFCleanAmp="      << prTOFParClean[0]    <<
    "prTOFCleanMean="     << prTOFParClean[1]    <<
    "prTOFCleanSigma="    << prTOFParClean[2]    <<
    "prTOFCleanKurtosis=" << prTOFParClean[3]    <<
    "prTOFCleanSkew="     << prTOFParClean[4]    <<
    "prTOFCleanInt="      << prTOFParClean[5]    <<
    "prTOFCleanChi2="     << prTOFParClean[6]    <<
    "prTOFCleanMax="      << prTOFParClean[7]    <<
    "prTOFCleanMatrix.="  << prTOFCleanMatrix    ;

    *sampleFile << "CleanSamples" <<

    "elCleanAmpFreeKS="      << elParCleanFreeKS[0]    <<
    "elCleanMeanFreeKS="     << elParCleanFreeKS[1]    <<
    "elCleanSigmaFreeKS="    << elParCleanFreeKS[2]    <<
    "elCleanKurtosisFreeKS=" << elParCleanFreeKS[3]    <<
    "elCleanSkewFreeKS="     << elParCleanFreeKS[4]    <<
    "elCleanIntFreeKS="      << elParCleanFreeKS[5]    <<
    "elCleanChi2FreeKS="     << elParCleanFreeKS[6]    <<
    "elCleanMaxFreeKS="      << elParCleanFreeKS[7]    <<

    "piCleanAmpFreeKS="      << piParCleanFreeKS[0]    <<
    "piCleanMeanFreeKS="     << piParCleanFreeKS[1]    <<
    "piCleanSigmaFreeKS="    << piParCleanFreeKS[2]    <<
    "piCleanKurtosisFreeKS=" << piParCleanFreeKS[3]    <<
    "piCleanSkewFreeKS="     << piParCleanFreeKS[4]    <<
    "piCleanIntFreeKS="      << piParCleanFreeKS[5]    <<
    "piCleanChi2FreeKS="     << piParCleanFreeKS[6]    <<
    "piCleanMaxFreeKS="      << piParCleanFreeKS[7]    ;

    *sampleFile << "CleanSamples" <<

    "kaCleanAmpFreeKS="      << kaParCleanFreeKS[0]    <<
    "kaCleanMeanFreeKS="     << kaParCleanFreeKS[1]    <<
    "kaCleanSigmaFreeKS="    << kaParCleanFreeKS[2]    <<
    "kaCleanKurtosisFreeKS=" << kaParCleanFreeKS[3]    <<
    "kaCleanSkewFreeKS="     << kaParCleanFreeKS[4]    <<
    "kaCleanIntFreeKS="      << kaParCleanFreeKS[5]    <<
    "kaCleanChi2FreeKS="     << kaParCleanFreeKS[6]    <<
    "kaCleanMaxFreeKS="      << kaParCleanFreeKS[7]    <<

    "prCleanAmpFreeKS="      << prParCleanFreeKS[0]    <<
    "prCleanMeanFreeKS="     << prParCleanFreeKS[1]    <<
    "prCleanSigmaFreeKS="    << prParCleanFreeKS[2]    <<
    "prCleanKurtosisFreeKS=" << prParCleanFreeKS[3]    <<
    "prCleanSkewFreeKS="     << prParCleanFreeKS[4]    <<
    "prCleanIntFreeKS="      << prParCleanFreeKS[5]    <<
    "prCleanChi2FreeKS="     << prParCleanFreeKS[6]    <<
    "prCleanMaxFreeKS="      << prParCleanFreeKS[7]    <<

    "deCleanAmpFreeKS="      << deParCleanFreeKS[0]    <<
    "deCleanMeanFreeKS="     << deParCleanFreeKS[1]    <<
    "deCleanSigmaFreeKS="    << deParCleanFreeKS[2]    <<
    "deCleanKurtosisFreeKS=" << deParCleanFreeKS[3]    <<
    "deCleanSkewFreeKS="     << deParCleanFreeKS[4]    <<
    "deCleanIntFreeKS="      << deParCleanFreeKS[5]    <<
    "deCleanChi2FreeKS="     << deParCleanFreeKS[6]    <<
    "deCleanMaxFreeKS="      << deParCleanFreeKS[7]    <<

    "kaTOFCleanAmpFreeKS="      << kaTOFParCleanFreeKS[0]    <<
    "kaTOFCleanMeanFreeKS="     << kaTOFParCleanFreeKS[1]    <<
    "kaTOFCleanSigmaFreeKS="    << kaTOFParCleanFreeKS[2]    <<
    "kaTOFCleanKurtosisFreeKS=" << kaTOFParCleanFreeKS[3]    <<
    "kaTOFCleanSkewFreeKS="     << kaTOFParCleanFreeKS[4]    <<
    "kaTOFCleanIntFreeKS="      << kaTOFParCleanFreeKS[5]    <<
    "kaTOFCleanChi2FreeKS="     << kaTOFParCleanFreeKS[6]    <<
    "kaTOFCleanMaxFreeKS="      << kaTOFParCleanFreeKS[7]    <<

    "prTOFCleanAmpFreeKS="      << prTOFParCleanFreeKS[0]    <<
    "prTOFCleanMeanFreeKS="     << prTOFParCleanFreeKS[1]    <<
    "prTOFCleanSigmaFreeKS="    << prTOFParCleanFreeKS[2]    <<
    "prTOFCleanKurtosisFreeKS=" << prTOFParCleanFreeKS[3]    <<
    "prTOFCleanSkewFreeKS="     << prTOFParCleanFreeKS[4]    <<
    "prTOFCleanIntFreeKS="      << prTOFParCleanFreeKS[5]    <<
    "prTOFCleanChi2FreeKS="     << prTOFParCleanFreeKS[6]    <<
    "prTOFCleanMaxFreeKS="      << prTOFParCleanFreeKS[7]    <<

    "\n";

    delete [] arrSigma;
    delete [] arrMean;
    delete [] arrCleanSigma;
    delete [] arrCleanMean;

    delete [] piParClean;
    delete [] elParClean;
    delete [] kaParClean;
    delete [] prParClean;
    delete [] deParClean;
    delete [] kaTOFParClean;
    delete [] prTOFParClean;

    delete [] piParCleanFreeKS;
    delete [] elParCleanFreeKS;
    delete [] kaParCleanFreeKS;
    delete [] prParCleanFreeKS;
    delete [] deParCleanFreeKS;
    delete [] kaTOFParCleanFreeKS;
    delete [] prTOFParCleanFreeKS;


    delete elCleanMatrix;
    delete piCleanMatrix;
    delete kaCleanMatrix;
    delete prCleanMatrix;
    delete deCleanMatrix;
    delete kaTOFCleanMatrix;
    delete prTOFCleanMatrix;


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
  if(arrMean[4]<=10){
    total->FixParameter(20,0);
    total->FixParameter(21,0);
    total->FixParameter(22,0);
    total->FixParameter(23,0);
    total->FixParameter(24,0);
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

  total->SetParName(20,"deAmplitude");
  total->SetParName(21,"deMean");
  total->SetParName(22,"deSigma");
  total->SetParName(23,"deKurtosis");
  total->SetParName(24,"deSkewness");

}
// --------------------------------------------------------------------------------------------
void SetFixedKSparameterForTotalFit(TF1 *total, Double_t pt)
{

  //
  // Set variable skewness
  if (fVariableSkewnes==0){
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[0] )  skewnessFix[kPi] = skewborder[0];
    if ( fitLookUpCleanSkew[kKa]->Eval(pt)>skewborder[0] )  skewnessFix[kKa] = skewborder[0];
    if ( fitLookUpCleanSkew[kPr]->Eval(pt)>skewborder[0] )  skewnessFix[kPr] = skewborder[0];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<=skewborder[1] ) skewnessFix[kPi] = skewborder[1];
    if ( fitLookUpCleanSkew[kKa]->Eval(pt)<=skewborder[1] ) skewnessFix[kKa] = skewborder[1];
    if ( fitLookUpCleanSkew[kPr]->Eval(pt)<=skewborder[1] ) skewnessFix[kPr] = skewborder[1];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[1]  ) skewnessFix[kPi] = fitLookUpCleanSkew[kPi]->Eval(pt);
    if ( fitLookUpCleanSkew[kKa]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kKa]->Eval(pt)>skewborder[1]  ) skewnessFix[kKa] = fitLookUpCleanSkew[kKa]->Eval(pt);
    if ( fitLookUpCleanSkew[kPr]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kPr]->Eval(pt)>skewborder[1]  ) skewnessFix[kPr] = fitLookUpCleanSkew[kPr]->Eval(pt);
  }
  if (fVariableSkewnes==1){
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[0] )  skewnessFix[kPi] = skewborder[0];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<=skewborder[1] ) skewnessFix[kPi] = skewborder[1];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[1]  ) skewnessFix[kPi] = fitLookUpCleanSkew[kPi]->Eval(pt);
    //
    skewnessFix[kKa] = kaSkewnessScan;
    skewnessFix[kPr] = prSkewnessScan;
    skewnessFix[kEl] = 0;
    skewnessFix[kDe] = 0;
  }


  // Set kurtosis
  if (fixedK){
    total->FixParameter(3, kurtosisFix[kEl]);
    total->FixParameter(8, kurtosisFix[kPi]);
    total->FixParameter(13,kurtosisFix[kKa]);
    total->FixParameter(18,kurtosisFix[kPr]);
    total->FixParameter(23,kurtosisFix[kDe]);
  }

  // Set skewness
  if (fixedS) {
    total->FixParameter(4, skewnessFix[kEl]);
    total->FixParameter(9, skewnessFix[kPi]);
    total->FixParameter(14,skewnessFix[kKa]);
    total->FixParameter(19,skewnessFix[kPr]);
    total->FixParameter(24,skewnessFix[kDe]);
  }


}
// --------------------------------------------------------------------------------------------
void SetFixedKSparameterForTotalFitForPions(TF1 *total)
{

  // Set kurtosis
  total->FixParameter(3, kurtosisFix[kEl]);

  total->SetParLimits(8,kurtosisMin,kurtosisMax);

  total->FixParameter(13,kurtosisFix[kKa]    );
  total->FixParameter(18,kurtosisFix[kPr]  );
  total->FixParameter(23,kurtosisFix[kDe]);

  // Set skewness
  total->FixParameter(4, skewnessFix[kEl]);

  total->SetParLimits(9,skewMin,skewMax);

  total->FixParameter(14,skewnessFix[kKa]    );
  total->FixParameter(19,skewnessFix[kPr]  );
  total->FixParameter(24,skewnessFix[kDe]);

}
// --------------------------------------------------------------------------------------------
void SetFixedKSparameterForIndividualFits(TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4, TF1 *g5)
{

  // Set kurtosis
  if (fixedK) {
    g1->FixParameter(3,kurtosisFix[kEl]);
    g2->FixParameter(3,kurtosisFix[kPi]);
    g3->FixParameter(3,kurtosisFix[kKa]);
    g4->FixParameter(3,kurtosisFix[kPr]);
    g5->FixParameter(3,kurtosisFix[kDe]);
  }

  // Set Skewness
  if (fixedS) {
    g1->FixParameter(4,skewnessFix[kEl]);
    g2->FixParameter(4,skewnessFix[kPi]);
    g3->FixParameter(4,skewnessFix[kKa]);
    g4->FixParameter(4,skewnessFix[kPr]);
    g5->FixParameter(4,skewnessFix[kDe]);
  }

}
// -------------------------------------------------------------------------------------------------------
void SetFixedKSparameterForCleanSampleFit(TF1 *asyGaus, ParticleType pSpecy, Double_t pt)
{


  //
  // Set variable skewness
  if (fVariableSkewnes==0){
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[0] ) skewnessFix[kPi] = skewborder[0];
    if ( fitLookUpCleanSkew[kKa]->Eval(pt)>skewborder[0] ) skewnessFix[kKa] = skewborder[0];
    if ( fitLookUpCleanSkew[kPr]->Eval(pt)>skewborder[0] ) skewnessFix[kPr] = skewborder[0];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<=skewborder[1] ) skewnessFix[kPi] = skewborder[1];
    if ( fitLookUpCleanSkew[kKa]->Eval(pt)<=skewborder[1] ) skewnessFix[kKa] = skewborder[1];
    if ( fitLookUpCleanSkew[kPr]->Eval(pt)<=skewborder[1] ) skewnessFix[kPr] = skewborder[1];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[1]  ) skewnessFix[kPi] = fitLookUpCleanSkew[kPi]->Eval(pt);
    if ( fitLookUpCleanSkew[kKa]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kKa]->Eval(pt)>skewborder[1]  ) skewnessFix[kKa] = fitLookUpCleanSkew[kKa]->Eval(pt);
    if ( fitLookUpCleanSkew[kPr]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kPr]->Eval(pt)>skewborder[1]  ) skewnessFix[kPr] = fitLookUpCleanSkew[kPr]->Eval(pt);
  }

  if (fVariableSkewnes==1){
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[0] ) skewnessFix[kPi] = skewborder[0];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<=skewborder[1] ) skewnessFix[kPi] = skewborder[1];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[1]  ) skewnessFix[kPi] = fitLookUpCleanSkew[kPi]->Eval(pt);
    //
    skewnessFix[kKa] = kaSkewnessScan;
    skewnessFix[kPr] = prSkewnessScan;
    skewnessFix[kEl] = 0;
    skewnessFix[kDe] = 0;
  }

  // Fix kurtoris and skewness for a given particle specy
  switch (pSpecy)
  {
    case kEl:
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kEl]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kEl]);
    }
    break;
    case kPi:
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kPi]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kPi]);
    }
    break;
    case kKa:
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kKa]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kKa]);
    }
    break;
    case kPr:
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kPr]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kPr]);
    }
    break;
    case kDe:
    {
      if (fixedK) asyGaus->FixParameter(3,kurtosisFix[kDe]);
      if (fixedS) asyGaus->FixParameter(4,skewnessFix[kDe]);
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
    ScaledExWrtClean = (TObjArray*)ressFile->Get("ScaledExWrtClean");
    ScaledEx      = (TObjArray*)ressFile->Get("ScaledEx");
    OutlierSmooth = (TObjArray*)ressFile->Get("OulierSmooth");
    Windows       = (TObjArray*)ressFile->Get("Windows");
    //
    // previous fir params
    for (Int_t iPart=0; iPart<nParticle; iPart++){

      grMeanCleanWrtExScaled[iPart]  = (TGraphErrors*)ScaledExWrtClean->FindObject(Form("%sMeanSplineCleanWrtExScaled",particleStr[iPart].Data()));
      grSigmaCleanWrtExScaled[iPart] = (TGraphErrors*)ScaledExWrtClean->FindObject(Form("%sSigmaSplineCleanWrtExScaled",particleStr[iPart].Data()));

      grFitAmp[iPart]   = (TGraphErrors*)Fit ->FindObject(Form("%sFitAmp",particleStr[iPart].Data()));
      grFitMean[iPart]  = (TGraphErrors*)Fit ->FindObject(Form("%sFitMean",particleStr[iPart].Data()));
      grFitSigma[iPart] = (TGraphErrors*)Fit ->FindObject(Form("%sFitSigma",particleStr[iPart].Data()));
      //
      grCleanAmp[iPart]   = (TGraphErrors*)Clean->FindObject(Form("%sCleanAmp",particleStr[iPart].Data()));
      grCleanMean[iPart]  = (TGraphErrors*)Clean->FindObject(Form("%sCleanMean",particleStr[iPart].Data()));
      grCleanSigma[iPart] = (TGraphErrors*)Clean->FindObject(Form("%sCleanSigma",particleStr[iPart].Data()));
      //
      grCleanAmpScaled[iPart]   = (TGraphErrors*)ScaledClean->FindObject(Form("%sCleanAmpScaled",particleStr[iPart].Data()));
      grCleanMeanScaled[iPart]  = (TGraphErrors*)ScaledClean->FindObject(Form("%sCleanMeanScaled",particleStr[iPart].Data()));
      grCleanSigmaScaled[iPart] = (TGraphErrors*)ScaledClean->FindObject(Form("%sCleanSigmaScaled",particleStr[iPart].Data()));
      //
      grFitAmpSmooth[iPart]   = (TGraphErrors*)Fit->FindObject(Form("%sFitAmpSmooth",particleStr[iPart].Data()));
      grFitMeanSmooth[iPart]  = (TGraphErrors*)Fit->FindObject(Form("%sFitMeanSmooth",particleStr[iPart].Data()));
      grFitSigmaSmooth[iPart] = (TGraphErrors*)Fit->FindObject(Form("%sFitSigmaSmooth",particleStr[iPart].Data()));
      //
      grMeanModExScaled[iPart] = (TGraphErrors*)ScaledEx->FindObject(Form("%sMeanModExScaled",particleStr[iPart].Data()));
      grMeanExScaled[iPart]    = (TGraphErrors*)ScaledEx->FindObject(Form("%sMeanSplineExScaled",particleStr[iPart].Data()));
      grSigmaExScaled[iPart]   = (TGraphErrors*)ScaledEx->FindObject(Form("%sSigmaSplineExScaled",particleStr[iPart].Data()));
      //
      grFitAmpOutlierSmooth[iPart]   = (TGraphErrors*)OutlierSmooth->FindObject(Form("%sFitAmpOutlierSmooth",particleStr[iPart].Data()));
      grFitSigmaOutlierSmooth[iPart] = (TGraphErrors*)OutlierSmooth->FindObject(Form("%sFitSigmaOutlierSmooth",particleStr[iPart].Data()));
      //
      grCleanAmpScaledSmooth[iPart]   = (TGraphErrors*)ScaledClean->FindObject(Form("%sCleanAmpScaledSmooth",particleStr[iPart].Data()));
      grCleanMeanScaledSmooth[iPart]  = (TGraphErrors*)ScaledClean->FindObject(Form("%sCleanMeanScaledSmooth",particleStr[iPart].Data()));
      grCleanSigmaScaledSmooth[iPart] = (TGraphErrors*)ScaledClean->FindObject(Form("%sCleanSigmaScaledSmooth",particleStr[iPart].Data()));
      //
      grCleanMeanSmooth[iPart]= (TGraphErrors*)Clean->FindObject(Form("%sCleanMeanSmooth",particleStr[iPart].Data()));
      //
      grFreeKSAmp[iPart]      = (TGraphErrors*)Clean->FindObject(Form("%sCleanAmpFreeKS",particleStr[iPart].Data()));
      grFreeKSMean[iPart]     = (TGraphErrors*)Clean->FindObject(Form("%sCleanMeanFreeKS",particleStr[iPart].Data()));
      grFreeKSSigma[iPart]    = (TGraphErrors*)Clean->FindObject(Form("%sCleanSigmaFreeKS",particleStr[iPart].Data()));
      grFreeKSSkew[iPart]     = (TGraphErrors*)Clean->FindObject(Form("%sCleanSkewFreeKS",particleStr[iPart].Data()));
      grFreeKSKurtosis[iPart] = (TGraphErrors*)Clean->FindObject(Form("%sCleanKurtosisFreeKS",particleStr[iPart].Data()));

    }

  }

}
// -------------------------------------------------------------------------------------------------------
void  SetParamsForKSEstimate(Int_t iks, TF1 *total, TH1D *h1DKS, Double_t *arrMean, Double_t *arrSigma)
{

  //
  //    Set mean amp and sigma par limits for KS estimate
  //

  //   Part to make setting for fit params
  maxBin  = h1DKS->GetMaximum();
  maxBin0 = TMath::Min(h1DKS->GetBinContent(h1DKS->FindBin(arrMean[0])),maxBin);
  maxBin1 = TMath::Min(h1DKS->GetBinContent(h1DKS->FindBin(arrMean[1])),maxBin);
  maxBin2 = TMath::Min(h1DKS->GetBinContent(h1DKS->FindBin(arrMean[2])),maxBin);
  maxBin3 = TMath::Min(h1DKS->GetBinContent(h1DKS->FindBin(arrMean[3])),maxBin);

  Double_t fitWinMaxBin      = 0.2;                      // only used for 0th iteration   // TO FIX

  // Restrictions on the total fit
  elAmpMin = 1e-4;
  elAmpMax = maxBin/10.;
  piAmpMin = 1e-4;
  piAmpMax = maxBin + maxBin*fitWinMaxBin;
  kaAmpMin = 1e-4;
  kaAmpMax = maxBin + maxBin*fitWinMaxBin;
  prAmpMin = 1e-4;
  prAmpMax = maxBin + maxBin*fitWinMaxBin;

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
  (iks==0) ? total->SetParLimits(3,kurtosisMin,kurtosisMax) : total->FixParameter(3,kurtosisFix[kEl]);
  total->SetParLimits(4,skewMin,skewMax);

  total->SetParLimits(5,piAmpMin,piAmpMax);
  total->SetParLimits(6,piMeanMin,piMeanMax);
  total->SetParLimits(7,piSigmaMin,piSigmaMax);
  (iks==0) ? total->SetParLimits(8,kurtosisMin,kurtosisMax) : total->FixParameter(8,kurtosisFix[kPi]);
  total->SetParLimits(9,skewMin,skewMax);

  total->SetParLimits(10,kaAmpMin,kaAmpMax);
  total->SetParLimits(11,kaMeanMin,kaMeanMax);
  total->SetParLimits(12,kaSigmaMin,kaSigmaMax);
  (iks==0) ? total->SetParLimits(13,kurtosisMin,kurtosisMax) : total->FixParameter(13,kurtosisFix[kKa]);
  total->SetParLimits(14,skewMin,skewMax);

  total->SetParLimits(15,prAmpMin,prAmpMax);
  total->SetParLimits(16,prMeanMin,prMeanMax);
  total->SetParLimits(17,prSigmaMin,prSigmaMax);
  (iks==0) ? total->SetParLimits(18,kurtosisMin,kurtosisMax) : total->FixParameter(18,kurtosisFix[kPr]);
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
  meantmparr[4] = arrMean[4];

  //   Part to make setting for fit params
  maxBin  = h1D->GetBinContent(h1D->GetMaximumBin());
  maxBin0 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[0])),maxBin);
  maxBin1 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[1])),maxBin);
  maxBin2 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[2])),maxBin);
  maxBin3 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[3])),maxBin);
  maxBin4 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[4])),maxBin);

  elAmpMin = 1e-4;
  elAmpMax = maxBin*1.2;
  piAmpMin = 1e-4;
  piAmpMax = maxBin*1.2;
  kaAmpMin = 1e-4;
  kaAmpMax = maxBin2*1.2;
  prAmpMin = 1e-4;
  prAmpMax = maxBin3*1.2;
  deAmpMin = 1e-4;
  deAmpMax = maxBin4*1.2;

  elMeanMin = TMath::Max(30.,(meantmparr[0]-arrMeanWindow[0]));
  piMeanMin = TMath::Max(30.,(meantmparr[1]-arrMeanWindow[1]));
  kaMeanMin = TMath::Max(30.,(meantmparr[2]-arrMeanWindow[2]));
  prMeanMin = TMath::Max(30.,(meantmparr[3]-arrMeanWindow[3]));
  deMeanMin = TMath::Max(30.,(meantmparr[4]-arrMeanWindow[4]));

  elMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[0]+arrMeanWindow[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[1]+arrMeanWindow[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[2]+arrMeanWindow[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[3]+arrMeanWindow[3]));
  deMeanMax = TMath::Min((Double_t)dEdxMax,(meantmparr[4]+arrMeanWindow[4]));

  elSigmaMin = (KSestimation) ? arrSigma[0]/2. : arrSigma[0]+2;
  piSigmaMin = (KSestimation) ? arrSigma[1]/2. : arrSigma[1]+2;
  kaSigmaMin = (KSestimation) ? arrSigma[2]/2. : arrSigma[2]+1;
  prSigmaMin = (KSestimation) ? arrSigma[3]/2. : arrSigma[3]+1;
  deSigmaMin = (KSestimation) ? arrSigma[4]/2. : arrSigma[4]+1;

  elSigmaMax = arrSigma[0]*2.;
  piSigmaMax = arrSigma[1]*2.;
  kaSigmaMax = arrSigma[2]*2.;
  prSigmaMax = arrSigma[3]*2.;
  deSigmaMax = arrSigma[4]*2.;

}
// -------------------------------------------------------------------------------------------------------
TF1  *FitParamInRange(TGraphErrors* gr, TString funcType, Double_t min, Double_t max)
{

  //
  // Fit proton TGraph with pol2
  //

  TF1 *f = new TF1("f",funcType,min,max);
  gr->Fit(f,"QNR+");
  f->SetRange(0.2,fMomRangeUp);
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
TH1D *GetClean1DSlice(TH2D *h2Clean, TString parName, Int_t ptDownBin)
{

  //
  //  get 1D momentum slice
  //

  TH1D *h1Clean   = (TH1D*)h2Clean->ProjectionY(Form(parName,ptDownBin),ptDownBin,ptDownBin);
  h1Clean->Rebin(rebinFactor); h1Clean->Scale(1./rebinFactor);
  Double_t binWidth = h1Clean->GetXaxis()->GetBinWidth(50);
  h1Clean->Scale(1./binWidth);
  // Double_t mean = h1Clean->GetMean();
  // Double_t rms  = h1Clean->GetRMS();
  // Double_t rangeDown = TMath::Max((Double_t)dEdxMin,mean-5.*rms);
  // Double_t rangeUp   = TMath::Min((Double_t)dEdxMax,mean+5.*rms);
  // h1Clean->GetXaxis()->SetRangeUser(rangeDown,rangeUp);

  return h1Clean;
}
// -------------------------------------------------------------------------------------------------------
void FitParticleSampleForKSEstimate(Int_t iks, TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy)
{

  //
  // Fit Clean sample Slice
  //
  // helper numbers for setparams
  TVirtualFitter::SetMaxIterations(100000);
  Double_t maxBinClean = hClean->GetMaximum();
  if ((maxBinClean<4.)) {
    parSample[0] = 1.;      // avoid floating point exception for muons
    parSample[1] = arrMean[pSpecy];
    parSample[2] = arrSigma[pSpecy];
    parSample[3] = 2.;
    parSample[4] = 0.;
    parSample[5] = 1.;
    parSample[6] = 1.;
    cout << "particle = " << pSpecy << "  is yok  " << endl;
    if (pSpecy==4) fDeuteronsExist=0;
    return;
  }
  Double_t mean = arrMean[pSpecy];
  Double_t rms  = arrSigma[pSpecy];
  //
  // Mand Maximum position
  TH1D *htmp = (TH1D*)hClean->Clone();
  htmp->GetXaxis()->SetRangeUser(htmp->GetBinCenter(htmp->FindBin(mean-3*rms)),htmp->GetBinCenter(htmp->FindBin(mean+3*rms)));
  Double_t meantmp = htmp->GetBinCenter(htmp->GetMaximumBin());
  delete htmp;
  //
  // fit function
  TF1 *asyGaus=NULL;
  if (pSpecy==kEl){
    asyGaus = new TF1("asyGaus",fitFunctionGenGaus,meantmp-electronCleanWindow,meantmp+electronCleanWindow);
  } else {
    asyGaus = new TF1("asyGaus",fitFunctionGenGaus,meantmp-sampleFitWindowDown*rms,meantmp+sampleFitWindowUp*rms);
  }
  asyGaus->SetParNames("Abundance","Mean","Sigma","Kurtosis","Skewness");
  asyGaus->SetParLimits(0,1.,maxBinClean);
  asyGaus->SetParLimits(1,meantmp-2*rms,meantmp+rms);
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
  // For deutetorns estimate should be done via clean samples
  maxBin = h1DKS->GetBinContent(h1DKS->GetMaximumBin());

  elMeanMin = TMath::Max(30.,(arrMean[0]-2.*arrSigma[0]));                     // electron
  piMeanMin = TMath::Max(30.,(arrMean[1]-fitWinMeanLowNsigma*arrSigma[1]));    // pion
  kaMeanMin = TMath::Max(30.,(arrMean[2]-fitWinMeanLowNsigma*arrSigma[2]));    // kaon
  prMeanMin = TMath::Max(30.,(arrMean[3]-fitWinMeanLowNsigma*arrSigma[3]));    // proton

  elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+2.*arrSigma[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+fitWinMeanUpNsigma*arrSigma[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+fitWinMeanUpNsigma*arrSigma[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+fitWinMeanUpNsigma*arrSigma[3]));

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

  (iks==0)? g1->SetParLimits(3,kurtosisMin,kurtosisMax) : g1->FixParameter(3,kurtosisFix[kEl]);
  (iks==0)? g2->SetParLimits(3,kurtosisMin,kurtosisMax) : g2->FixParameter(3,kurtosisFix[kPi]);
  (iks==0)? g3->SetParLimits(3,kurtosisMin,kurtosisMax) : g3->FixParameter(3,kurtosisFix[kKa]);
  (iks==0)? g4->SetParLimits(3,kurtosisMin,kurtosisMax) : g4->FixParameter(3,kurtosisFix[kPr]);

  g1->SetParLimits(4,skewMin,skewMax);
  g2->SetParLimits(4,skewMin,skewMax);
  g3->SetParLimits(4,skewMin,skewMax);
  g4->SetParLimits(4,skewMin,skewMax);

  //
  if (arrMean[0]>=10) h1DKS->Fit(g1,"QNR");
  if (arrMean[1]>=10) h1DKS->Fit(g2,"QNR+");
  if (arrMean[2]>=10) h1DKS->Fit(g3,"QNR+");
  if (arrMean[3]>=10) h1DKS->Fit(g4,"QNR+");

  // Get Fit parameters from the individual fits and set for the total fit
  const Int_t nParameters = 20;
  Double_t par[nParameters] = {0};
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
  Double_t totalChi2 = total->GetChisquare();
  Double_t normNDF   = total->GetNDF();
  Double_t totChi2   = (normNDF<1) ? 1 : totalChi2/normNDF;

  for (Int_t i=0;i<nParameters;i++) pars[i] = total->GetParameter(i);
  pars[nParameters] = totChi2;
  pars[nParameters+1] = g1->GetMaximumX();
  pars[nParameters+2] = g2->GetMaximumX();
  pars[nParameters+3] = g3->GetMaximumX();
  pars[nParameters+4] = g4->GetMaximumX();

  delete total;
  delete g1;
  delete g2;
  delete g3;
  delete g4;
}
// -------------------------------------------------------------------------------------------------------
TMatrixD *FitParticleSample(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy, Double_t pt)
{

  //
  // Fit Clean sample slice
  //
  TVirtualFitter::SetMaxIterations(100000);
  Double_t mean = arrMean[pSpecy];
  Double_t rms  = arrSigma[pSpecy];
  Double_t maxBinClean = hClean->GetMaximum();
  //
  // Mand Maximum position
  TH1D *htmp = (TH1D*)hClean->Clone();
  htmp->GetXaxis()->SetRangeUser(htmp->GetBinCenter(htmp->FindBin(mean-3*rms)),htmp->GetBinCenter(htmp->FindBin(mean+3*rms)));
  Double_t meantmp = htmp->GetBinCenter(htmp->GetMaximumBin());
  delete htmp;
  //
  TLine* linetmp = new TLine(meantmp, 0., meantmp, maxBinClean);
  linetmp->SetLineColor(kGreen); linetmp->SetLineWidth(2);
  hClean->GetListOfFunctions()->Add(linetmp);
  //
  TLine* line = new TLine(arrMean[pSpecy], 0., arrMean[pSpecy], maxBinClean);
  line->SetLineColor(kRed); line->SetLineWidth(2);
  hClean->GetListOfFunctions()->Add(line);
  //
  Double_t integral = hClean->Integral(hClean->FindBin(meantmp-sampleFitWindowDown*rms),hClean->FindBin(meantmp+sampleFitWindowUp*rms));

  Bool_t histCheck = ( (mean>20 && mean<1020) && (rms>0 && rms<200) && (maxBinClean>=1 && maxBinClean<1e+30) && (integral>=5 && integral<1e+30));
  if (!histCheck) {   // there is too few entries, fit would fail
    parSample[0] = 1.;      // avoid floating point exception for muons
    parSample[1] = arrMean[pSpecy];
    parSample[2] = arrSigma[pSpecy];
    parSample[3] = 2.;
    parSample[4] = 0.;
    parSample[5] = 1.;
    parSample[6] = 1.;
    return new TMatrixD(5,5);
  }

  // prepare normalised samples
  if (normalisedCleanSamples && integral > 1e-4) {
    hClean->Scale(1./integral);
    maxBinClean /= integral;
  }

  // fit function
  TF1 *asyGaus=NULL;
  if (pSpecy==kEl){
    asyGaus = new TF1("asyGaus",fitFunctionGenGaus,meantmp-electronCleanWindow,meantmp+electronCleanWindow);
  } else {
    asyGaus = new TF1("asyGaus",fitFunctionGenGaus,meantmp-sampleFitWindowDown*rms,meantmp+sampleFitWindowUp*rms);
  }
  asyGaus->SetParNames("Amplitude","Mean","Sigma","Kurtosis","Skewness");
  // asyGaus->SetNpx(nBinsInLookUpTable);
  asyGaus ->SetLineColor(2);
  asyGaus ->SetLineWidth(2);

  // Set Parameters
  asyGaus->SetParLimits(0,1.,maxBinClean*2);
  asyGaus->SetParLimits(1,meantmp-2.*rms,meantmp+rms);
  asyGaus->SetParLimits(2,rms/2. ,3.*rms);
  asyGaus->SetParLimits(3,kurtosisMin,kurtosisMax);
  asyGaus->SetParLimits(4,skewMin,skewMax);
  if ((fixedK || fixedS)) SetFixedKSparameterForCleanSampleFit(asyGaus,pSpecy,pt);

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
  parSample[7] = asyGaus->GetMaximumX();

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
  //
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
void FitParticleSampleFreeKS(TH1D *hClean, Double_t *parSample, Double_t *arrMean, Double_t *arrSigma, ParticleType pSpecy, Double_t pt)
{

  //
  // Fit Clean sample slice
  //

  TVirtualFitter::SetMaxIterations(100000);
  Double_t mean = arrMean[pSpecy];
  Double_t rms  = arrSigma[pSpecy];
  Double_t maxBinClean = hClean->GetMaximum();
  //
  // Mand Maximum position
  TH1D *htmp = (TH1D*)hClean->Clone();
  htmp->GetXaxis()->SetRangeUser(htmp->GetBinCenter(htmp->FindBin(mean-3*rms)),htmp->GetBinCenter(htmp->FindBin(mean+3*rms)));
  Double_t meantmp = htmp->GetBinCenter(htmp->GetMaximumBin());
  delete htmp;
  //
  TLine* linetmp = new TLine(meantmp, 0., meantmp, maxBinClean);
  linetmp->SetLineColor(kGreen); linetmp->SetLineWidth(2);
  hClean->GetListOfFunctions()->Add(linetmp);
  //
  TLine* line = new TLine(mean, 0., mean, maxBinClean);
  line->SetLineColor(kRed); line->SetLineWidth(2);
  hClean->GetListOfFunctions()->Add(line);
  //
  Double_t integral = hClean->Integral(hClean->FindBin(meantmp-sampleFitWindowDown*rms),hClean->FindBin(meantmp+sampleFitWindowUp*rms));

  Bool_t histCheck = ( (mean>20 && mean<1020) && (rms>0 && rms<200) && (maxBinClean>=1 && maxBinClean<1e+30) && (integral>=5 && integral<1e+30));
  if (!histCheck) {   // there is too few entries, fit would fail
    parSample[0] = 1.;
    parSample[1] = arrMean[pSpecy];
    parSample[2] = arrSigma[pSpecy];
    parSample[3] = 2.;
    parSample[4] = 0.;
    parSample[5] = 1.;
    parSample[6] = 1.;
    return;
  }

  // fit function
  TF1 *asyGaus=NULL;
  if (pSpecy==kEl){
    asyGaus = new TF1("asyGaus",fitFunctionGenGaus,meantmp-electronCleanWindow,meantmp+electronCleanWindow);
  } else {
    asyGaus = new TF1("asyGaus",fitFunctionGenGaus,meantmp-sampleFitWindowDown*rms,meantmp+sampleFitWindowUp*rms);
  }
  asyGaus->SetParNames("Amplitude","Mean","Sigma","Kurtosis","Skewness");
  // asyGaus->SetNpx(nBinsInLookUpTable);
  asyGaus ->SetLineColor(kBlue);
  asyGaus ->SetLineWidth(2);

  // Set Parameters
  asyGaus->SetParLimits(0,1.,maxBinClean);
  asyGaus->SetParLimits(1,meantmp-2*rms,meantmp+2*rms);
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
  parSample[7] = asyGaus->GetMaximumX();
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
void GetCleanExParams(TH1D *hClean, Double_t *arrSigma, Double_t *arrCleanSigma, Double_t *arrCleanMean, Double_t *arrCleanMax, ParticleType pSpecy)
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
  gaus->SetParLimits(1,mean-2.*rms,mean+rms);
  gaus->SetParLimits(2,rms/2. ,3.*rms);

  // retrieve fit parameters
  hClean->Fit(gaus,"QNMR");
  arrCleanMean[pSpecy]  = gaus->GetParameter(1);
  arrCleanSigma[pSpecy] = gaus->GetParameter(2);
  arrCleanMax[pSpecy]   = gaus->GetMaximumX();

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
void ResetParametersOfEachFunction(TF1 *g1,TF1 *g2,TF1 *g3,TF1 *g4,TF1 *g5,TF1 *total,TH1D *h1D,const Double_t pt)
{

  //
  // Reset the parameters of individual particle after the total fit
  //
  TF1 *tot = (TF1*)total->Clone();

  // Electron
  g1->FixParameter(0,total->GetParameter(0));
  g1->FixParameter(1,total->GetParameter(1));
  g1->FixParameter(2,total->GetParameter(2));
  g1->FixParameter(3,total->GetParameter(3));
  g1->FixParameter(4,total->GetParameter(4));

  // Pion
  g2->FixParameter(0,total->GetParameter(5));
  g2->FixParameter(1,total->GetParameter(6));
  g2->FixParameter(2,total->GetParameter(7));
  g2->FixParameter(3,total->GetParameter(8));
  g2->FixParameter(4,total->GetParameter(9));

  // Kaon
  g3->FixParameter(0,total->GetParameter(10));
  g3->FixParameter(1,total->GetParameter(11));
  g3->FixParameter(2,total->GetParameter(12));
  g3->FixParameter(3,total->GetParameter(13));
  g3->FixParameter(4,total->GetParameter(14));

  // Proton
  g4->FixParameter(0,total->GetParameter(15));
  g4->FixParameter(1,total->GetParameter(16));
  g4->FixParameter(2,total->GetParameter(17));
  g4->FixParameter(3,total->GetParameter(18));
  g4->FixParameter(4,total->GetParameter(19));

  // deuteron
  g5->FixParameter(0,total->GetParameter(20));
  g5->FixParameter(1,total->GetParameter(21));
  g5->FixParameter(2,total->GetParameter(22));
  g5->FixParameter(3,total->GetParameter(23));
  g5->FixParameter(4,total->GetParameter(24));

  // Set colors
  g1->SetLineColor(4);
  g2->SetLineColor(kOrange-1);
  g3->SetLineColor(kGreen+2);
  g4->SetLineColor(kMagenta+1);
  g5->SetLineColor(kBlack);

  g1->SetFillColor(4);
  g2->SetFillColor(kOrange-1);
  g3->SetFillColor(kGreen+2);
  g4->SetFillColor(kMagenta+1);
  g5->SetFillColor(kBlack);

  g1->SetFillStyle(3004); g1->Draw("fchist");
  g2->SetFillStyle(3005); g2->Draw("fchist");
  g3->SetFillStyle(3006); g3->Draw("fchist");
  g4->SetFillStyle(3007); g4->Draw("fchist");
  g5->SetFillStyle(3008); g5->Draw("fchist");

  // Set ranges
  g1->SetRange(0,h1D->GetNbinsX());
  g2->SetRange(0,h1D->GetNbinsX());
  g3->SetRange(0,h1D->GetNbinsX());
  g4->SetRange(0,h1D->GetNbinsX());
  g5->SetRange(0,h1D->GetNbinsX());

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
  leg->AddEntry(g5   ,"d"   ,"F");

  // Attach functions and TLegend to the histogram
  h1D->GetListOfFunctions()->Add(tot);
  h1D->GetListOfFunctions()->Add(g1);
  h1D->GetListOfFunctions()->Add(g2);
  h1D->GetListOfFunctions()->Add(g3);
  h1D->GetListOfFunctions()->Add(g4);
  h1D->GetListOfFunctions()->Add(g5);
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

    Double_t bassLevel = smoothLevel;
    if (branchName.Contains("deFitAmp"))   bassLevel=5.;
    if (branchName.Contains("elFitAmp"))   bassLevel=5.;
    if (branchName.Contains("kaFitAmp"))   bassLevel=8.;
    if (branchName.Contains("kaFitSigma")) bassLevel=9.;



    //
    // create the graph object
    t->Draw(drawStr,"","goff");
    grTmp[objIndex] = new TGraphErrors(nSlice,t->GetV2(),t->GetV1());
    grTmp[objIndex]->SetName(branchName);
    grTmp[objIndex]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    grTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grTmp[objIndex]->SetLineColor(kBlack);
    grTmp[objIndex]->SetMarkerStyle(7);
    grTmp[objIndex]->SetMarkerColor(kBlack);
    //
    //
    grSmoothTmp[objIndex] = (TGraph*)grTmp[objIndex]->Clone();
    grSmoothTmp[objIndex]->SetName(branchName.Append("Smooth"));
    grSmoothTmp[objIndex]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    grSmoothTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grSmoothTmp[objIndex]->SetLineColor(kRed);
    grSmoothTmp[objIndex]->SetMarkerStyle(7);
    grSmoothTmp[objIndex]->SetMarkerColor(kRed);
    //
    //
    Int_t nPoints = 0;
    if (branchName.Contains("el")) nPoints = t->Draw(drawStr,Form("p>%f",0.),"goff");
    else nPoints = t->Draw(drawStr,Form("p>%f",applySmoothAfter),"goff");
    TGraphErrors *grtmp = new TGraphErrors(nPoints,t->GetV2(),t->GetV1());
    gs = new TGraphSmooth("normal");
    TGraph *grtmpsmooth = gs->SmoothSuper((TGraph*)grtmp,"",bassLevel);
    //
    for (Int_t i=0; i<grSmoothTmp[objIndex]->GetN();i++){
      for (Int_t j=0; j<grtmpsmooth->GetN();j++){

        //
        // pointer treatment
        if ( grSmoothTmp[objIndex]->GetX()[i] == grtmpsmooth->GetX()[j] ) {
          *(grSmoothTmp[objIndex]->GetY()+i) = grtmpsmooth->GetY()[j];
        }

      }
    }
    // //
    // // Smooth Grapgh
    // gs = new TGraphSmooth("normal");
    // if ( branchName.Contains("Kurtosis") || branchName.Contains("Skew") ) {
    //   grSmoothTmp[objIndex] = (TGraph*)grTmp[objIndex]->Clone();
    // } else {
    //   grSmoothTmp[objIndex] = gs->SmoothSuper((TGraph*)grTmp[objIndex],"",bassLevel);
    // }
    // grSmoothTmp[objIndex]->SetName(branchName.Append("Smooth"));
    // grSmoothTmp[objIndex]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    // grSmoothTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    // grSmoothTmp[objIndex]->SetLineColor(kRed);
    // grSmoothTmp[objIndex]->SetMarkerStyle(7);
    // grSmoothTmp[objIndex]->SetMarkerColor(kRed);
    //
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

    if (branchName.Contains("deFitAmp")) continue;
    if (branchName.Contains("elFitAmp")) smoothSeg=10.;
    if (branchName.Contains("kaFitSigma")) smoothSeg=10.;

    // create the graph object
    t->Draw(drawStr,"","goff");
    grTmp[objIndex] = new TGraphErrors(nSlice,t->GetV2(),t->GetV1());
    grTmp[objIndex] = RemoveOutliers(grTmp[objIndex],TMath::Nint((Double_t)nSlice/smoothSeg));
    grTmp[objIndex] = RemoveOutliers(grTmp[objIndex],TMath::Nint((Double_t)nSlice/smoothSeg));
    grTmp[objIndex]->SetName(branchName.Append("Outlier"));
    grTmp[objIndex]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    grTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grTmp[objIndex]->SetLineColor(kBlack);
    grTmp[objIndex]->SetMarkerStyle(7);
    grTmp[objIndex]->SetMarkerColor(kBlack);
    //
    grSmoothTmp[objIndex] = (TGraph*)grTmp[objIndex]->Clone();
    grSmoothTmp[objIndex]->SetName(branchName.Append("Smooth"));
    grSmoothTmp[objIndex]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    grSmoothTmp[objIndex]->GetYaxis()->SetTitle(branchName);
    grSmoothTmp[objIndex]->SetLineColor(kRed);
    grSmoothTmp[objIndex]->SetMarkerStyle(7);
    grSmoothTmp[objIndex]->SetMarkerColor(kRed);
    cout << grSmoothTmp[objIndex]->GetName() << endl;
    //
    //
    Int_t nPoints = 0;
    if (branchName.Contains("el")) nPoints = t->Draw(drawStr,Form("p>%f",0.),"goff");
    else nPoints = t->Draw(drawStr,Form("p>%f",applySmoothAfter),"goff");
    TGraphErrors *grtmp = new TGraphErrors(nPoints,t->GetV2(),t->GetV1());
    gs = new TGraphSmooth("normal");
    TGraph *grtmpsmooth = gs->SmoothSuper((TGraph*)grtmp,"",9.);
    //
    for (Int_t i=0; i<grSmoothTmp[objIndex]->GetN();i++){
      for (Int_t j=0; j<grtmpsmooth->GetN();j++){
        if ( grSmoothTmp[objIndex]->GetX()[i] == grtmpsmooth->GetX()[j] ) {
          // pointer treatment
          *(grSmoothTmp[objIndex]->GetY()+i) = grtmpsmooth->GetY()[j];
        }
      }
    }
    //
    // Add graphs to ObjArrays
    arr -> AddAt(grTmp[objIndex]      ,2*objIndex);
    arr -> AddAt(grSmoothTmp[objIndex],2*objIndex+1);

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
  Double_t scaleFactor = grexp->GetY()[index]/grfit->GetY()[index];  // calculate scaling factor

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

  //   get a new string for newname assignment
  TString newname  = sample;
  TGraphErrors *grsample  = NULL;
  TGraphErrors *grfit     = NULL;
  //
  //   get the graphs to be used in scaling
  if (newname.Contains("Clean")) grsample  = (TGraphErrors*)arrSample -> FindObject(sample);
  else grsample  = (TGraphErrors*)arr -> FindObject(sample);
  //
  if (newname.Contains("Clean")) grfit = (TGraphErrors*)arr -> FindObject(fit);
  else grfit = (TGraphErrors*)arrSample -> FindObject(fit);

  //   Find the index corresponding to safeP value
  Int_t index = GetPtIndex(grsample,safeP);

  //   calculate scaling factor  samplemean/fitmean
  Double_t scaleFactor = (grsample->GetY()[index]<2) ? 1: grfit->GetY()[index]/grsample->GetY()[index];

  //   draw scaled graph and return it
  sample.Append("*%f:p");
  TString drawStr  = Form(sample,scaleFactor);
  if (newname.Contains("Clean")) newname.Append("Scaled");
  else newname.Append("CleanWrtExScaled");

  //   create scaled graph
  tSample->Draw(drawStr,"","goff");
  TGraphErrors *grScaled = new TGraphErrors(nSlice,tSample->GetV2(),tSample->GetV1());
  grScaled->SetName(newname);
  grScaled->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  grScaled->GetYaxis()->SetTitle(newname);
  grScaled->SetLineColor(kGreen);
  grScaled->SetMarkerStyle(7);
  grScaled->SetMarkerColor(kGreen);
  //
  grScaledSmooth[k] = (TGraph*)grScaled->Clone();
  grScaledSmooth[k]->SetName(newname.Append("Smooth"));
  grScaledSmooth[k]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  grScaledSmooth[k]->GetYaxis()->SetTitle(newname);
  grScaledSmooth[k]->SetLineColor(kRed);
  grScaledSmooth[k]->SetMarkerStyle(7);
  grScaledSmooth[k]->SetMarkerColor(kRed);
  //
  //
  Int_t nPoints = 0;
  if (newname.Contains("el")) nPoints = tSample->Draw(drawStr,Form("p>%f",0.),"goff");
  else nPoints = tSample->Draw(drawStr,Form("p>%f",applySmoothAfter),"goff");
  TGraphErrors *grtmp = new TGraphErrors(nPoints,tSample->GetV2(),tSample->GetV1());
  gs = new TGraphSmooth("normal");
  TGraph *grtmpsmooth = gs->SmoothSuper((TGraph*)grtmp,"",4);
  //
  for (Int_t i=0; i<grScaledSmooth[k]->GetN();i++){
    for (Int_t j=0; j<grtmpsmooth->GetN();j++){
      if ( grScaledSmooth[k]->GetX()[i] == grtmpsmooth->GetX()[j] ) {
        // pointer treatment
        *(grScaledSmooth[k]->GetY()+i) = grtmpsmooth->GetY()[j];
      }
    }
  }

  return grScaled;
}
// -------------------------------------------------------------------------------------------------------
TCanvas *GetFitResCanvas(TObjArray *arr, TObjArray *arrclean, TObjArray *arrscaledclean, TObjArray *arrscaledex)
{

  //
  //   Fit Results Canvas
  //

  TGraphErrors *sigma[nParticle];
  TGraphErrors *sigmaSmooth[nParticle];
  TGraphErrors *mean[nParticle];
  TGraphErrors *meanSmooth[nParticle];
  TGraphErrors *amp[nParticle];
  TGraphErrors *ampSmooth[nParticle];

  TGraphErrors *cleansigma[nParticle];
  TGraphErrors *cleansigmascaled[nParticle];
  TGraphErrors *sigmaexscaled[nParticle];

  TGraphErrors *cleanmean[nParticle];
  TGraphErrors *cleanmeanscaled[nParticle];
  TGraphErrors *meanexscaled[nParticle];

  TGraphErrors *exmeanMod[nParticle];

  for (Int_t i=0; i<nParticle-2; i++)
  {

    sigma[i]=0x0;  sigmaSmooth[i]=0x0;
    mean[i]=0x0;   meanSmooth[i]=0x0;
    amp[i]=0x0;    ampSmooth[i]=0x0;
    cleansigma[i]=0x0;  cleansigmascaled[i]=0x0;  sigmaexscaled[i]=0x0;
    cleanmean[i]=0x0;   cleanmeanscaled[i]=0x0;   meanexscaled[i]=0x0;
    exmeanMod[i]=0x0;

    sigma[i]       = (TGraphErrors*)arr->FindObject(Form("%sFitSigma",particleStr[i].Data()));       sigma[i]->SetLineColor(kBlack);
    mean[i]        = (TGraphErrors*)arr->FindObject(Form("%sFitMean",particleStr[i].Data()));        mean[i]->SetLineColor(kBlack);
    amp[i]         = (TGraphErrors*)arr->FindObject(Form("%sFitAmp",particleStr[i].Data()));         amp[i]->SetLineColor(kBlack);

    sigmaSmooth[i] = (TGraphErrors*)arr->FindObject(Form("%sFitSigmaSmooth",particleStr[i].Data())); sigmaSmooth[i]->SetLineColor(kRed);
    meanSmooth[i]  = (TGraphErrors*)arr->FindObject(Form("%sFitMeanSmooth",particleStr[i].Data()));  meanSmooth[i]->SetLineColor(kRed);
    ampSmooth[i]   = (TGraphErrors*)arr->FindObject(Form("%sFitAmpSmooth",particleStr[i].Data()));   ampSmooth[i]->SetLineColor(kRed);

    cleansigma[i]  = (TGraphErrors*)arrclean->FindObject(Form("%sCleanSigma",particleStr[i].Data()));
    cleansigma[i]->SetLineColor(kGreen+2);  cleansigma[i]->SetMarkerColor(kGreen+2);
    cleanmean[i]   = (TGraphErrors*)arrclean->FindObject(Form("%sCleanMean",particleStr[i].Data()));
    cleanmean[i]->SetLineColor(kGreen+2);   cleanmean[i]->SetMarkerColor(kGreen+2);

    cleansigmascaled[i] = (TGraphErrors*)arrscaledclean->FindObject(Form("%sCleanSigmaScaled",particleStr[i].Data()));
    if (cleansigmascaled[i]) { cleansigmascaled[i]->SetLineColor(kBlue+2); cleansigmascaled[i]->SetMarkerColor(kBlue+2); }
    cleanmeanscaled[i]  = (TGraphErrors*)arrscaledclean->FindObject(Form("%sCleanMeanScaled",particleStr[i].Data()));
    if (cleanmeanscaled[i]) { cleanmeanscaled[i]->SetLineColor(kBlue+2);  cleanmeanscaled[i]->SetMarkerColor(kBlue+2); }

    sigmaexscaled[i]    = (TGraphErrors*)arrscaledex->FindObject(Form("%sSigmaSplineExScaled",particleStr[i].Data()));
    if (sigmaexscaled[i])  { sigmaexscaled[i]->SetLineColor(kMagenta); sigmaexscaled[i]->SetMarkerColor(kMagenta); }
    meanexscaled[i]     = (TGraphErrors*)arrscaledex->FindObject(Form("%sMeanSplineExScaled",particleStr[i].Data()));
    if (meanexscaled[i])   { meanexscaled[i]->SetLineColor(kMagenta); meanexscaled[i]->SetMarkerColor(kMagenta); }

    exmeanMod[i]   = (TGraphErrors*)arrscaledex->FindObject(Form("%sMeanModExScaled",particleStr[i].Data()));
    if (exmeanMod[i]) { exmeanMod[i]->SetLineColor(kOrange);  exmeanMod[i]->SetMarkerColor(kOrange); }


    sigma[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    mean[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
    amp[i]->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");

  }


  TCanvas *can = new TCanvas("can", "can", 6000, 4000);
  can->Divide(6,4);

  // ++++++++++++++++++  Sigmas +++++++++++++++++
  can->cd(1);
  can->GetPad(1)->SetGrid();
  sigma[kEl]->GetYaxis()->SetTitle("electron sigma");
  sigma[kEl]->Draw("alp");
  sigmaSmooth[kEl]->Draw("lp");
  if (cleansigmascaled[kEl]) cleansigmascaled[kEl]->Draw("lp");
  if (sigmaexscaled[kEl]) sigmaexscaled[kEl]->Draw("lp");
  if (cleansigma[kEl]) cleansigma[kEl]->Draw("lp");

  can->cd(2);
  can->GetPad(2)->SetGrid();
  sigma[kPi]->GetYaxis()->SetTitle("pion sigma");
  sigma[kPi]->Draw("alp");
  sigmaSmooth[kPi]->Draw("lp");
  if (cleansigmascaled[kPi]) cleansigmascaled[kPi]->Draw("lp");
  if (sigmaexscaled[kPi]) sigmaexscaled[kPi]->Draw("lp");
  if (cleansigma[kPi]) cleansigma[kPi]->Draw("lp");

  can->cd(3);
  can->GetPad(3)->SetGrid();
  sigma[kKa]->GetYaxis()->SetTitle("kaon sigma");
  sigma[kKa]->Draw("alp");
  sigmaSmooth[kKa]->Draw("lp");
  if (cleansigmascaled[kKa]) cleansigmascaled[kKa]->Draw("lp");
  if (sigmaexscaled[kKa]) sigmaexscaled[kKa]->Draw("lp");
  if (cleansigma[kKa]) cleansigma[kKa]->Draw("lp");

  can->cd(4);
  can->GetPad(4)->SetGrid();
  sigma[kPr]->GetYaxis()->SetTitle("proton sigma");
  sigma[kPr]->Draw("alp");
  sigmaSmooth[kPr]->Draw("lp");
  if (cleansigmascaled[kPr]) cleansigmascaled[kPr]->Draw("lp");
  if (sigmaexscaled[kPr]) sigmaexscaled[kPr]->Draw("lp");
  if (cleansigma[kPr]) cleansigma[kPr]->Draw("lp");

  can->cd(5);
  can->GetPad(5)->SetGrid();
  sigma[kDe]->GetYaxis()->SetTitle("deuteron sigma");
  sigma[kDe]->Draw("alp");
  sigmaSmooth[kDe]->Draw("lp");
  if (cleansigmascaled[kDe]) cleansigmascaled[kDe]->Draw("lp");
  if (sigmaexscaled[kDe]) sigmaexscaled[kDe]->Draw("lp");
  if (cleansigma[kDe]) cleansigma[kDe]->Draw("lp");

  can->cd(6);
  can->GetPad(6)->SetGrid();
  TGraphErrors *chi2          = (TGraphErrors*)arr->FindObject("totChi2");
  chi2->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  chi2->GetYaxis()->SetTitle("Chi2");
  chi2->GetYaxis()->SetTitleOffset(1.6);
  chi2->SetMarkerStyle(7);
  chi2->SetMarkerColor(kRed);
  chi2->SetLineColor(kRed);
  chi2->Draw("alp");

  // ++++++++++++++++++  Means +++++++++++++++++
  can->cd(7);
  can->GetPad(7)->SetGrid();
  mean[kEl]->GetYaxis()->SetTitle("electron mean");
  mean[kEl]->Draw("alp");
  meanSmooth[kEl]->Draw("lp");
  if (cleanmeanscaled[kEl]) cleanmeanscaled[kEl]->Draw("lp");
  if (meanexscaled[kEl]) meanexscaled[kEl]->Draw("lp");
  if (cleanmean[kEl]) cleanmean[kEl]->Draw("lp");
  // if (exmeanMod[kEl]) exmeanMod[kEl]->Draw("lp");

  can->cd(8);
  can->GetPad(8)->SetGrid();
  mean[kPi]->GetYaxis()->SetTitle("pion mean");
  mean[kPi]->Draw("alp");
  meanSmooth[kPi]->Draw("lp");
  if (cleanmeanscaled[kPi]) cleanmeanscaled[kPi]->Draw("lp");
  if (meanexscaled[kPi]) meanexscaled[kPi]->Draw("lp");
  if (cleanmean[kPi]) cleanmean[kPi]->Draw("lp");
  // if (exmeanMod[kPi]) exmeanMod[kPi]->Draw("lp");

  can->cd(9);
  can->GetPad(9)->SetGrid();
  mean[kKa]->GetYaxis()->SetTitle("kaon mean");
  mean[kKa]->Draw("alp");
  meanSmooth[kKa]->Draw("lp");
  if (cleanmeanscaled[kKa]) cleanmeanscaled[kKa]->Draw("lp");
  if (meanexscaled[kKa]) meanexscaled[kKa]->Draw("lp");
  if (cleanmean[kKa]) cleanmean[kKa]->Draw("lp");
  // if (exmeanMod[kKa]) exmeanMod[kKa]->Draw("lp");

  can->cd(10);
  can->GetPad(10)->SetGrid();
  mean[kPr]->GetYaxis()->SetTitle("proton mean");
  mean[kPr]->Draw("alp");
  meanSmooth[kPr]->Draw("lp");
  if (cleanmeanscaled[kPr]) cleanmeanscaled[kPr]->Draw("lp");
  if (meanexscaled[kPr]) meanexscaled[kPr]->Draw("lp");
  if (cleanmean[kPr]) cleanmean[kPr]->Draw("lp");
  // if (exmeanMod[kPr]) exmeanMod[kPr]->Draw("lp");

  can->cd(11);
  can->GetPad(11)->SetGrid();
  mean[kDe]->GetYaxis()->SetTitle("deuteron mean");
  mean[kDe]->Draw("alp");
  meanSmooth[kDe]->Draw("lp");
  if (cleanmeanscaled[kDe]) cleanmeanscaled[kDe]->Draw("lp");
  if (meanexscaled[kDe]) meanexscaled[kDe]->Draw("lp");
  if (cleanmean[kDe]) cleanmean[kDe]->Draw("lp");
  // if (exmeanMod[kDe]) exmeanMod[kDe]->Draw("lp");


  can->cd(12);
  can->GetPad(12)->SetGrid();
  TGraphErrors *maxRes          = (TGraphErrors*)arr->FindObject("maxRes");
  maxRes->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  maxRes->GetYaxis()->SetTitle("Maximum Residual per Slice");
  maxRes->GetYaxis()->SetTitleOffset(1.6);
  maxRes->SetMarkerStyle(7);
  maxRes->SetMarkerColor(kBlue);
  maxRes->SetLineColor(kBlue);
  maxRes->Draw("alp");

  // ++++++++++++++++++  Amplitudes +++++++++++++++++
  can->cd(13);
  can->GetPad(13)->SetGrid();
  amp[kEl]->GetYaxis()->SetTitle("electron amplitude");
  amp[kEl]->Draw("alp");
  ampSmooth[kEl]->Draw("lp");

  can->cd(14);
  can->GetPad(14)->SetGrid();
  amp[kPi]->GetYaxis()->SetTitle("pion amplitude");
  amp[kPi]->Draw("alp");
  ampSmooth[kPi]->Draw("lp");

  can->cd(15);
  can->GetPad(15)->SetGrid();
  amp[kKa]->GetYaxis()->SetTitle("kaon amplitude");
  amp[kKa]->Draw("alp");
  ampSmooth[kKa]->Draw("lp");

  can->cd(16);
  can->GetPad(16)->SetGrid();
  amp[kPr]->GetYaxis()->SetTitle("proton amplitude");
  amp[kPr]->Draw("alp");
  ampSmooth[kPr]->Draw("lp");

  can->cd(17);
  can->GetPad(17)->SetGrid();
  amp[kDe]->GetYaxis()->SetTitle("deuteron amplitude");
  amp[kDe]->Draw("alp");
  ampSmooth[kDe]->Draw("lp");

  can->cd(18);
  can->GetPad(18)->SetGrid();
  TGraphErrors *maxPer = (TGraphErrors*)arr->FindObject("maxPer");
  maxPer->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  maxPer->GetYaxis()->SetTitle("|data-fit|/data");
  maxPer->GetYaxis()->SetTitleOffset(1.6);
  maxPer->SetMarkerStyle(7);
  maxPer->SetMarkerColor(kGreen);
  maxPer->SetLineColor(kGreen);
  maxPer->Draw("alp");

  // ++++++++++++++++++  Yields +++++++++++++++++
  can->cd(19);
  can->GetPad(19)->SetGrid();
  TGraphErrors *elint       = (TGraphErrors*)arr->FindObject("elFitInt");
  elint->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  elint->GetYaxis()->SetTitle("electron yield");
  elint->GetYaxis()->SetTitleOffset(1.6);
  elint->Draw("alp");

  can->cd(20);
  can->GetPad(20)->SetGrid();
  TGraphErrors *piint       = (TGraphErrors*)arr->FindObject("piFitInt");
  piint->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  piint->GetYaxis()->SetTitle("pion yield");
  piint->GetYaxis()->SetTitleOffset(1.6);
  piint->Draw("alp");

  can->cd(21);
  can->GetPad(21)->SetGrid();
  TGraphErrors *kaint       = (TGraphErrors*)arr->FindObject("kaFitInt");
  kaint->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  kaint->GetYaxis()->SetTitle("kaon yield");
  kaint->GetYaxis()->SetTitleOffset(1.6);
  kaint->Draw("alp");

  can->cd(22);
  can->GetPad(22)->SetGrid();
  TGraphErrors *print       = (TGraphErrors*)arr->FindObject("prFitInt");
  print->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  print->GetYaxis()->SetTitle("proton yield");
  print->GetYaxis()->SetTitleOffset(1.6);
  print->Draw("alp");

  can->cd(23);
  can->GetPad(23)->SetGrid();
  TGraphErrors *deint       = (TGraphErrors*)arr->FindObject("deFitInt");
  deint->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  deint->GetYaxis()->SetTitle("deuteron yield");
  deint->GetYaxis()->SetTitleOffset(1.6);
  deint->Draw("alp");

  can->cd(24);
  can->GetPad(24)->SetGrid();
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

  maxBin = h1D->GetMaximum();

  maxBin0 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[0])),maxBin);
  maxBin1 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[1])),maxBin);
  maxBin2 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[2])),maxBin);
  maxBin3 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[3])),maxBin);
  maxBin4 = TMath::Min(h1D->GetBinContent(h1D->FindBin(arrMean[4])),maxBin);

  TH1D * h1FitWindows  = (TH1D*)h1D->Clone(); h1FitWindows->SetName(Form("FitWindows_%d",islice));
  SetStyleFor1DSlice(h1FitWindows);
  TH1D * h1InitialFits = (TH1D*)h1D->Clone(); h1InitialFits->SetName(Form("InitialFits_%d",islice));
  SetStyleFor1DSlice(h1InitialFits);

  // Fix the mean position between +-1 sigma
  elMeanMin = TMath::Max(30.,(arrMean[0]-fitWinMeanLowNsigma*arrSigma[0]));            // electron
  piMeanMin = TMath::Max(30.,(arrMean[1]-fitWinMeanLowNsigma*arrSigma[1]));            // pion
  kaMeanMin = TMath::Max(30.,(arrMean[2]-fitWinMeanLowNsigma*arrSigma[2]));            // kaon
  prMeanMin = TMath::Max(30.,(arrMean[3]-fitWinMeanLowNsigma*arrSigma[3]));            // proton
  deMeanMin = TMath::Max(30.,(arrMean[4]-fitWinMeanLowNsigma*arrSigma[4]));            // deuteron

  elMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[0]+fitWinMeanUpNsigma*arrSigma[0]));
  piMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[1]+fitWinMeanUpNsigma*arrSigma[1]));
  kaMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[2]+fitWinMeanUpNsigma*arrSigma[2]));
  prMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[3]+fitWinMeanUpNsigma*arrSigma[3]));
  deMeanMax = TMath::Min((Double_t)dEdxMax,(arrMean[4]+fitWinMeanUpNsigma*arrSigma[4]));

  TF1 *g1b = new TF1("g1b",fitFunctionGenGaus,elMeanMin,elMeanMax);
  TF1 *g2b = new TF1("g2b",fitFunctionGenGaus,piMeanMin,piMeanMax);
  TF1 *g3b = new TF1("g3b",fitFunctionGenGaus,kaMeanMin,kaMeanMax);
  TF1 *g4b = new TF1("g4b",fitFunctionGenGaus,prMeanMin,prMeanMax);
  TF1 *g5b = new TF1("g5b",fitFunctionGenGaus,deMeanMin,deMeanMax);

  // g1b->SetNpx(nBinsInLookUpTable);
  // g2b->SetNpx(nBinsInLookUpTable);
  // g3b->SetNpx(nBinsInLookUpTable);
  // g4b->SetNpx(nBinsInLookUpTable);

  g1b->SetLineColor(4);         // Blue
  g2b->SetLineColor(kOrange-1);
  g3b->SetLineColor(kGreen+2);
  g4b->SetLineColor(kMagenta+1);
  g5b->SetLineColor(kBlack);

  g1b->SetParLimits(0,0.,maxBin);
  g2b->SetParLimits(0,0.,maxBin);
  g3b->SetParLimits(0,0.,maxBin);
  g4b->SetParLimits(0,0.,maxBin);
  g5b->SetParLimits(0,0.,maxBin);

  g1b->SetParLimits(1,elMeanMin,elMeanMax);
  g2b->SetParLimits(1,piMeanMin,piMeanMax);
  g3b->SetParLimits(1,kaMeanMin,kaMeanMax);
  g4b->SetParLimits(1,prMeanMin,prMeanMax);
  g5b->SetParLimits(1,deMeanMin,deMeanMax);

  Double_t fitWinSigma = 0.5;
  g1b->SetParLimits(2,arrSigma[0]-fitWinSigma*arrSigma[0],arrSigma[0]+fitWinSigma*arrSigma[0]);
  g2b->SetParLimits(2,arrSigma[1]-fitWinSigma*arrSigma[1],arrSigma[1]+fitWinSigma*arrSigma[1]);
  g3b->SetParLimits(2,arrSigma[2]-fitWinSigma*arrSigma[2],arrSigma[2]+fitWinSigma*arrSigma[2]);
  g4b->SetParLimits(2,arrSigma[3]-fitWinSigma*arrSigma[3],arrSigma[3]+fitWinSigma*arrSigma[3]);
  g5b->SetParLimits(2,arrSigma[4]-fitWinSigma*arrSigma[4],arrSigma[4]+fitWinSigma*arrSigma[4]);

  g1b->FixParameter(3,0.);
  g2b->FixParameter(3,0.);
  g3b->FixParameter(3,0.);
  g4b->FixParameter(3,0.);
  g5b->FixParameter(3,0.);


  if (arrMean[0]>0) h1FitWindows->Fit(g1b,"QR");
  if (arrMean[1]>0) h1FitWindows->Fit(g2b,"QR+");
  if (arrMean[2]>0) h1FitWindows->Fit(g3b,"QR+");
  if (arrMean[3]>0) h1FitWindows->Fit(g4b,"QR+");
  if (arrMean[4]>0) h1FitWindows->Fit(g5b,"QR+");

  ///////////////////////////////////////////////////////////////////////////////////
  // add TLegend
  TLegend *leg = new TLegend(0.75, 0.6, 0.85, 0.85);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);  // White
  //   leg->SetHeader(Form("p=%4.3f",pt));
  TLegend *legClone = (TLegend*)leg->Clone();

  // Plot lines for electrons
  TLine* lineg1 = new TLine(arrMean[0], 0., arrMean[0], maxBin);
  lineg1     ->SetLineColor(4);
  lineg1     ->SetLineWidth(2);
  TLine* lineg1up = new TLine(elMeanMax, 0., elMeanMax, maxBin);
  lineg1up   ->SetLineColor(4);
  lineg1up   ->SetLineWidth(1);
  lineg1up   ->SetLineStyle(2);
  TLine* lineg1down = new TLine(elMeanMin, 0., elMeanMin, maxBin);
  lineg1down ->SetLineColor(4);
  lineg1down ->SetLineWidth(1);
  lineg1down ->SetLineStyle(2);
  if (arrMean[0]>0 && arrMean[0]<dEdxMax) {
    h1FitWindows->GetListOfFunctions()->Add(lineg1);
    h1FitWindows->GetListOfFunctions()->Add(lineg1up);
    h1FitWindows->GetListOfFunctions()->Add(lineg1down);
  }

  // Plot lines for Pions
  TLine* lineg2 = new TLine(arrMean[1], 0., arrMean[1], maxBin);
  lineg2      ->SetLineColor(kOrange-1);
  lineg2      ->SetLineWidth(2);
  TLine* lineg2up = new TLine(piMeanMax, 0., piMeanMax, maxBin);
  lineg2up    ->SetLineColor(kOrange-1);
  lineg2up    ->SetLineWidth(1);
  lineg2up    ->SetLineStyle(2);
  TLine* lineg2down = new TLine(piMeanMin, 0., piMeanMin, maxBin);
  lineg2down  ->SetLineColor(kOrange-1);
  lineg2down  ->SetLineWidth(1);
  lineg2down  ->SetLineStyle(2);
  if (arrMean[1]>0 && arrMean[1]<dEdxMax) {
    h1FitWindows->GetListOfFunctions()->Add(lineg2);
    h1FitWindows->GetListOfFunctions()->Add(lineg2up);
    h1FitWindows->GetListOfFunctions()->Add(lineg2down);
  }

  // Plot lines for Kaon
  TLine* lineg3 = new TLine(arrMean[2], 0., arrMean[2], maxBin);
  lineg3      ->SetLineColor(kGreen+2);
  lineg3      ->SetLineWidth(2);
  TLine* lineg3up = new TLine(kaMeanMax, 0., kaMeanMax, maxBin);
  lineg3up    ->SetLineColor(kGreen+2);
  lineg3up    ->SetLineWidth(1);
  lineg3up    ->SetLineStyle(2);
  TLine* lineg3down = new TLine(kaMeanMin, 0., kaMeanMin, maxBin);
  lineg3down  ->SetLineColor(kGreen+2);
  lineg3down  ->SetLineWidth(1);
  lineg3down  ->SetLineStyle(2);
  if (arrMean[2]>0 && arrMean[2]<dEdxMax) {
    h1FitWindows->GetListOfFunctions()->Add(lineg3);
    h1FitWindows->GetListOfFunctions()->Add(lineg3up);
    h1FitWindows->GetListOfFunctions()->Add(lineg3down);
  }

  // Plot lines for Proton
  TLine* lineg4 = new TLine(arrMean[3], 0., arrMean[3], maxBin);
  lineg4     ->SetLineColor(kMagenta+1);
  lineg4     ->SetLineWidth(2);
  TLine* lineg4up = new TLine(prMeanMax, 0., prMeanMax, maxBin);
  lineg4up   ->SetLineColor(kMagenta+1);
  lineg4up   ->SetLineWidth(1);
  lineg4up   ->SetLineStyle(2);
  TLine* lineg4down = new TLine(prMeanMin, 0., prMeanMin, maxBin);
  lineg4down ->SetLineColor(kMagenta+1);
  lineg4down ->SetLineWidth(1);
  lineg4down ->SetLineStyle(2);
  if (arrMean[3]>0 && arrMean[3]<dEdxMax) {
    h1FitWindows->GetListOfFunctions()->Add(lineg4);
    h1FitWindows->GetListOfFunctions()->Add(lineg4up);
    h1FitWindows->GetListOfFunctions()->Add(lineg4down);
  }

  // Plot lines for deuteron
  TLine* lineg5 = new TLine(arrMean[4], 0., arrMean[4], maxBin);
  lineg5     ->SetLineColor(kBlack);
  lineg5     ->SetLineWidth(2);
  TLine* lineg5up = new TLine(deMeanMax, 0., deMeanMax, maxBin);
  lineg5up   ->SetLineColor(kBlack);
  lineg5up   ->SetLineWidth(1);
  lineg5up   ->SetLineStyle(2);
  TLine* lineg5down = new TLine(deMeanMin, 0., deMeanMin, maxBin);
  lineg5down ->SetLineColor(kBlack);
  lineg5down ->SetLineWidth(1);
  lineg5down ->SetLineStyle(2);
  if (arrMean[4]>0 && arrMean[4]<dEdxMax) {
    h1FitWindows->GetListOfFunctions()->Add(lineg5);
    h1FitWindows->GetListOfFunctions()->Add(lineg5up);
    h1FitWindows->GetListOfFunctions()->Add(lineg5down);
  }

  leg->AddEntry(h1FitWindows,"Data","LEP");
  leg->AddEntry(g1b,"e","L");
  leg->AddEntry(g2b,"#pi","L");
  leg->AddEntry(g3b,"K","L");
  leg->AddEntry(g4b,"p","L");
  leg->AddEntry(g5b,"d","L");

  h1FitWindows->GetListOfFunctions()->Add(g1b);
  h1FitWindows->GetListOfFunctions()->Add(g2b);
  h1FitWindows->GetListOfFunctions()->Add(g3b);
  h1FitWindows->GetListOfFunctions()->Add(g4b);
  h1FitWindows->GetListOfFunctions()->Add(g5b);
  h1FitWindows->GetListOfFunctions()->Add(leg);

  ///////////////////////////////////////////////////////////////////////////////////

  // Initial fits
  TF1 * g1bClone = (TF1*)g1b->Clone();
  TF1 * g2bClone = (TF1*)g2b->Clone();
  TF1 * g3bClone = (TF1*)g3b->Clone();
  TF1 * g4bClone = (TF1*)g4b->Clone();
  TF1 * g5bClone = (TF1*)g5b->Clone();
  //   TF1 *g1bClone = new TF1(); g1b->Copy(*g1bClone);
  //   TF1 *g2bClone = new TF1(); g2b->Copy(*g2bClone);
  //   TF1 *g3bClone = new TF1(); g3b->Copy(*g3bClone);
  //   TF1 *g4bClone = new TF1(); g4b->Copy(*g4bClone);

  g1bClone->SetRange(dEdxMin,dEdxMax);
  g2bClone->SetRange(dEdxMin,dEdxMax);
  g3bClone->SetRange(dEdxMin,dEdxMax);
  g4bClone->SetRange(dEdxMin,dEdxMax);
  g5bClone->SetRange(dEdxMin,dEdxMax);

  legClone->AddEntry(h1InitialFits,"Data","LEP");
  legClone->AddEntry(g1bClone,"e","L");
  legClone->AddEntry(g2bClone,"#pi","L");
  legClone->AddEntry(g3bClone,"K","L");
  legClone->AddEntry(g4bClone,"p","L");
  legClone->AddEntry(g5bClone,"d","L");

  h1InitialFits->GetListOfFunctions()->Add(g1bClone);
  h1InitialFits->GetListOfFunctions()->Add(g2bClone);
  h1InitialFits->GetListOfFunctions()->Add(g3bClone);
  h1InitialFits->GetListOfFunctions()->Add(g4bClone);
  h1InitialFits->GetListOfFunctions()->Add(g5bClone);
  h1InitialFits->GetListOfFunctions()->Add(legClone);

  ///////////////////////////////////////////////////////////////////////////////////

  debugFile->GetFile()->cd();
  h1FitWindows  -> Write(Form("Windows_%d_%4.3f",islice,pt));
  h1InitialFits -> Write(Form("InitialFits_%d_%4.3f",islice,pt));

  delete g1bClone;
  delete g2bClone;
  delete g3bClone;
  delete g4bClone;
  delete g5bClone;
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
  for (Int_t iPart=0; iPart<nParticle; iPart++)
  {

    grMeanCleanWrtExScaled[iPart]=0x0;
    grSigmaCleanWrtExScaled[iPart]=0x0;

    grFitAmpSmooth[iPart]=0x0;   grFitAmp[iPart]=0x0;   grCleanAmp[iPart]=0x0;    grCleanAmpScaled[iPart]=0x0;
    grFitMeanSmooth[iPart]=0x0;  grFitMean[iPart]=0x0;  grCleanMean[iPart]=0x0;   grCleanMeanScaled[iPart]=0x0;
    grFitSigmaSmooth[iPart]=0x0; grFitSigma[iPart]=0x0; grCleanSigma[iPart]=0x0;  grCleanSigmaScaled[iPart]=0x0;
    //
    grMeanExScaled[iPart]=0x0;   grFitAmpOutlierSmooth[iPart]=0x0;
    grSigmaExScaled[iPart]=0x0;  grFitSigmaOutlierSmooth[iPart]=0x0;
    grMeanModExScaled[iPart]=0x0;
    //
    grCleanMeanSmooth[iPart]=0x0;
    grCleanAmpScaledSmooth[iPart]=0x0;
    grCleanMeanScaledSmooth[iPart]=0x0;
    grCleanSigmaScaledSmooth[iPart]=0x0;
    //
    grFreeKSAmp[iPart]=0x0;
    grFreeKSMean[iPart]=0x0;
    grFreeKSSigma[iPart]=0x0;
    grFreeKSSkew[iPart]=0x0;
    grFreeKSKurtosis[iPart]=0x0;
  }
  //
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
  cout << " main hist name = " << Form("h2Dall%s_%s",corrStr.Data(),centEtaStr.Data()) << endl;
  //
  if (fSign==0)       h2D = (TH2D *)inputFile->Get(Form("h2Dall%s_%s",corrStr.Data(),centEtaStr.Data()));
  else if (fSign==1 ) h2D = (TH2D *)inputFile->Get(Form("h2Dpos%s_%s",corrStr.Data(),centEtaStr.Data()));
  else if (fSign==-1) h2D = (TH2D *)inputFile->Get(Form("h2Dneg%s_%s",corrStr.Data(),centEtaStr.Data()));
  //
  // Clean Samples
  hCleanSamples[0]  = (TH2D *)inputFile->Get(Form("h2DCleanEl_KineCut%s_%s"           ,corrStr.Data(),centEtaStr.Data()));
  hCleanSamples[1]  = (TH2D *)inputFile->Get(Form("%s%s_%s",cleanPionSmapleName.Data(),corrStr.Data(),centEtaStr.Data()));
  hCleanSamples[2]  = (TH2D *)inputFile->Get(Form("h2DCleanKa_TOFTRD%s_%s"            ,corrStr.Data(),centEtaStr.Data()));
  hCleanSamples[3]  = (TH2D *)inputFile->Get(Form("h2DCleanPr%s_%s"                   ,corrStr.Data(),centEtaStr.Data()));
  hCleanSamples[4]  = (TH2D *)inputFile->Get(Form("h2DCleanDe%s_%s"                   ,corrStr.Data(),centEtaStr.Data()));
  //
  if (fSign==0)       hCleanSamples[5] = (TH2D *)inputFile->Get(Form("h2DallKaTOF%s_%s"   ,corrStr.Data(),centEtaStr.Data()));
  else if (fSign==1 ) hCleanSamples[5] = (TH2D *)inputFile->Get(Form("h2DallKaTOFPos%s_%s",corrStr.Data(),centEtaStr.Data()));
  else if (fSign==-1) hCleanSamples[5] = (TH2D *)inputFile->Get(Form("h2DallKaTOFNeg%s_%s",corrStr.Data(),centEtaStr.Data()));
  if (fSign==0)       hCleanSamples[6] = (TH2D *)inputFile->Get(Form("h2DallPrTOF%s_%s"   ,corrStr.Data(),centEtaStr.Data()));
  else if (fSign==1 ) hCleanSamples[6] = (TH2D *)inputFile->Get(Form("h2DallPrTOFPos%s_%s",corrStr.Data(),centEtaStr.Data()));
  else if (fSign==-1) hCleanSamples[6] = (TH2D *)inputFile->Get(Form("h2DallPrTOFNeg%s_%s",corrStr.Data(),centEtaStr.Data()));
  //
  CreateSplinesFromExpectedsTree(expectedsFile);
  //
  hExpectedMean[5]  = (TH2D*)hExpectedMean[2]->Clone();  hExpectedMean[5]->SetName(Form("h2Expected_5_%s",centEtaStr.Data()));
  hExpectedMean[6]  = (TH2D*)hExpectedMean[3]->Clone();  hExpectedMean[6]->SetName(Form("h2Expected_6_%s",centEtaStr.Data()));
  hExpectedSigma[5]  = (TH2D*)hExpectedSigma[2]->Clone();  hExpectedSigma[5]->SetName(Form("h2ExpectedSigma_5_%s",centEtaStr.Data()));
  hExpectedSigma[6]  = (TH2D*)hExpectedSigma[3]->Clone();  hExpectedSigma[6]->SetName(Form("h2ExpectedSigma_6_%s",centEtaStr.Data()));
  //
  // Get the Profiles of expecteds
  TProfile *hProfExMean[nParticle], *hProfExSigma[nParticle];
  for (Int_t i=0;i<nParticle;i++){
    hProfExMean[i]=NULL; hProfExSigma[i]=NULL;
    hProfExMean[i] = hExpectedMean[i]->ProfileX(); hProfExSigma[i] = hExpectedSigma[i]->ProfileX();
    hProfExMean[i]->SetMarkerColor(kBlack);        hProfExSigma[i]->SetMarkerColor(kBlack);
    hProfExMean[i]->SetLineColor(kBlack);          hProfExSigma[i]->SetLineColor(kBlack);
  }
  // for proper error calculations
  cout << " =========== Sumw2 ============ " << endl;
  h2D -> Sumw2();
  for (Int_t i=0;i<nParticle;i++) hCleanSamples[i]->Sumw2();

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  cout << " =========== Write histograms into the debug file ============ " << endl;
  debugFile -> GetFile()->cd();
  h2D               -> Write("h2Dall");
  hCleanSamples[0]  -> Write("h2Electron");
  hCleanSamples[1]  -> Write("h2Pion");
  hCleanSamples[2]  -> Write("h2Kaon");
  hCleanSamples[3]  -> Write("h2Proton");
  hCleanSamples[4]  -> Write("h2Deuteron");
  hCleanSamples[5]  -> Write("h2DKaonTOF");
  hCleanSamples[6]  -> Write("h2DProtonTOF");
  //
  if (hExpectedMean[0]) hExpectedMean[0]  -> Write("h2ExpectedEl");
  if (hExpectedMean[1]) hExpectedMean[1]  -> Write("h2ExpectedPi");
  if (hExpectedMean[2]) hExpectedMean[2]  -> Write("h2ExpectedKa");
  if (hExpectedMean[3]) hExpectedMean[3]  -> Write("h2ExpectedPr");
  if (hExpectedMean[4]) hExpectedMean[4]  -> Write("h2ExpectedDe");
  if (hExpectedMean[5]) hExpectedMean[5]  -> Write("h2ExpectedKaTOF");
  if (hExpectedMean[6]) hExpectedMean[6]  -> Write("h2ExpectedPrTOF");
  if (hExpectedSigma[0]) hExpectedSigma[0] -> Write("h2ExpectedSigmaEl");
  if (hExpectedSigma[1]) hExpectedSigma[1] -> Write("h2ExpectedSigmaPi");
  if (hExpectedSigma[2]) hExpectedSigma[2] -> Write("h2ExpectedSigmaKa");
  if (hExpectedSigma[3]) hExpectedSigma[3] -> Write("h2ExpectedSigmaPr");
  if (hExpectedSigma[4]) hExpectedSigma[4] -> Write("h2ExpectedSigmaDe");
  if (hExpectedSigma[5]) hExpectedSigma[5] -> Write("h2ExpectedSigmaKaTOF");
  if (hExpectedSigma[6]) hExpectedSigma[6] -> Write("h2ExpectedSigmaPrTOF");
  //
  if (hProfExMean[0]) hProfExMean[0]  -> Write("hProfExpectedEl");
  if (hProfExMean[1]) hProfExMean[1]  -> Write("hProfExpectedPi");
  if (hProfExMean[2]) hProfExMean[2]  -> Write("hProfExpectedKa");
  if (hProfExMean[3]) hProfExMean[3]  -> Write("hProfExpectedPr");
  if (hProfExMean[4]) hProfExMean[4]  -> Write("hProfExpectedDe");
  if (hProfExMean[5]) hProfExMean[5]  -> Write("hProfExpectedKaTOF");
  if (hProfExMean[6]) hProfExMean[6]  -> Write("hProfExpectedPrTOF");
  if (hProfExSigma[0]) hProfExSigma[0] -> Write("hProfExpectedSigmaEl");
  if (hProfExSigma[1]) hProfExSigma[1] -> Write("hProfExpectedSigmaPi");
  if (hProfExSigma[2]) hProfExSigma[2] -> Write("hProfExpectedSigmaKa");
  if (hProfExSigma[3]) hProfExSigma[3] -> Write("hProfExpectedSigmaPr");
  if (hProfExSigma[4]) hProfExSigma[4] -> Write("hProfExpectedSigmaDe");
  if (hProfExSigma[5]) hProfExSigma[5] -> Write("hProfExpectedSigmaKaTOF");
  if (hProfExSigma[6]) hProfExSigma[6] -> Write("hProfExpectedSigmaPrTOF");
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  cout << " elKurt   = " << kurtosisFix[kEl] << endl;
  cout << " piKurt   = " << kurtosisFix[kPi]     << endl;
  cout << " kaKurt   = " << kurtosisFix[kKa]     << endl;
  cout << " prKurt   = " << kurtosisFix[kPr]   << endl;
  cout << " ============================================================================" << endl;
  cout << " ============================== New info is =================================" << endl;

  // Read kurtosis tabel tree
  TFile fKurtosisTable(kurtosisTableFile);
  TTree * treeKurtosisTable = (TTree*)fKurtosisTable.Get("treeId");

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
      kurtosisFix[kEl] = elKurtosisTable;
      kurtosisFix[kPi] = piKurtosisTable;
      kurtosisFix[kKa] = kaKurtosisTable;
      kurtosisFix[kPr] = prKurtosisTable;
      cout << " sign     = " << signTable     << endl;
      cout << " etaDown  = " << etaTable  << endl;
      cout << " centDown = " << centTable << endl;
      cout << " elKurt   = " << kurtosisFix[kEl] << endl;
      cout << " piKurt   = " << kurtosisFix[kPi]     << endl;
      cout << " kaKurt   = " << kurtosisFix[kKa]     << endl;
      cout << " prKurt   = " << kurtosisFix[kPr]   << endl;
      break;
    }
  }

}
// -------------------------------------------------------------------------------------------------------
void PrepareLineShapes(TH1D *h,TF1 *ftot,TF1 *g1, TF1 *g2, TF1 *g3, TF1 *g4, TF1 *g5, Int_t centbin, Int_t etabin, Int_t mombin, Int_t signbin, Int_t iIter)
{

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

}

TH2D * TH2Subset(TH2D *h,Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2)
{
  Double_t x1 = h->GetXaxis()->GetBinLowEdge(binx1);
  Double_t x2 = h->GetXaxis()->GetBinLowEdge(binx2)+h->GetXaxis()->GetBinWidth(binx2);
  Double_t y1 = h->GetYaxis()->GetBinLowEdge(biny1);
  Double_t y2 = h->GetYaxis()->GetBinLowEdge(biny2)+h->GetYaxis()->GetBinWidth(biny2);
  Int_t    n1 = binx1+binx2+1;
  Int_t    n2 = biny1+biny2+1;

  Int_t nbinsx = h->GetNbinsX();
  Int_t nbinsy = h->GetNbinsY();

  TH2D *hs = new TH2D(h->GetName(), h->GetTitle(), n1, x1, x2, n2, y1, y2);

  Double_t content, x, y;

  for (Int_t i=1; i<=nbinsx; i++) {
    for (Int_t j=1; j<=nbinsx; j++) {
      content = h->GetBinContent(i,j);
      x = h->GetXaxis()->GetBinCenter(i);
      y = h->GetYaxis()->GetBinCenter(j);
      hs->Fill(x, y, content);
    }
  }

  return hs;
}

void MakeAnimation(Int_t iIter, TH1D* h1ToGif, TH2D* h2Danimation, Double_t pt, Double_t ptMin, Double_t ptMax){
  TLine* lineSlice = new TLine(pt, 20., pt, 400);
  h2Danimation->GetListOfFunctions()->Add(lineSlice);
  TCanvas gifC;
  gifC.SetCanvasSize(900,1200);
  gifC.Divide(1,2);
  gifC.cd(1);
  gifC.GetPad(1)->Clear();
  gifC.GetPad(1)->SetLogz();
  gifC.GetPad(1)->SetPad(.2, .75, .8, .995);
  h2Danimation->SetTitle("");
  h2Danimation->GetXaxis()->SetRangeUser(ptMin,ptMax);
  h2Danimation->GetXaxis()->SetTitle("p(GeV/c)");
  h2Danimation->GetYaxis()->SetRangeUser(20,400);
  h2Danimation->GetYaxis()->SetTitle("TPC dE/dx Signal");
  h2Danimation->Draw("colz");
  gifC.cd(2);
  if (pt<1.) {
    h1ToGif->GetXaxis()->SetRangeUser(0,350);
  } else {
    h1ToGif->GetXaxis()->SetRangeUser(0,180);
  }
  gifC.GetPad(2)->SetLogy();
  gifC.GetPad(2)->SetPad(.005, .005, .990, .75);
  h1ToGif->SetTitle("");
  h1ToGif->GetXaxis()->SetTitle("TPC dE/dx Signal");
  h1ToGif->Draw();
  if (iIter == 3) gifC.Print(Form("Iteration3_%s.gif+100", centEtaStr.Data()));   // delay in between plots is 120*10 ms
  if (iIter == 4) gifC.Print(Form("Iteration4_%s.gif+100", centEtaStr.Data()));   // delay in between plots is 120*10 ms
  if (iIter == 5) gifC.Print(Form("Iteration5_%s.gif+100", centEtaStr.Data()));   // delay in between plots is 120*10 ms
  if (iIter == 6) gifC.Print(Form("Iteration6_%s.gif+100", centEtaStr.Data()));   // delay in between plots is 120*10 ms
  if (iIter == 7) gifC.Print(Form("Iteration7_%s.gif+100", centEtaStr.Data()));   // delay in between plots is 120*10 ms
  if (iIter == 8) gifC.Print(Form("Iteration8_%s.gif+100", centEtaStr.Data()));   // delay in between plots is 120*10 ms
  h2Danimation->GetListOfFunctions()->Remove(lineSlice);
}

void  ModifyExpectedMeanAndApplyVariableSkewness(Double_t *arrMean, Double_t *arrSigma, Double_t *arrMeanModif, Double_t pt)
{


  if (fVariableSkewnes==0){
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[0] ) skewnessFix[kPi] = skewborder[0];
    if ( fitLookUpCleanSkew[kKa]->Eval(pt)>skewborder[0] ) skewnessFix[kKa] = skewborder[0];
    if ( fitLookUpCleanSkew[kPr]->Eval(pt)>skewborder[0] ) skewnessFix[kPr] = skewborder[0];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<=skewborder[1] ) skewnessFix[kPi] = skewborder[1];
    if ( fitLookUpCleanSkew[kKa]->Eval(pt)<=skewborder[1] ) skewnessFix[kKa] = skewborder[1];
    if ( fitLookUpCleanSkew[kPr]->Eval(pt)<=skewborder[1] ) skewnessFix[kPr] = skewborder[1];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[1]  ) skewnessFix[kPi] = fitLookUpCleanSkew[kPi]->Eval(pt);
    if ( fitLookUpCleanSkew[kKa]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kKa]->Eval(pt)>skewborder[1]  ) skewnessFix[kKa] = fitLookUpCleanSkew[kKa]->Eval(pt);
    if ( fitLookUpCleanSkew[kPr]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kPr]->Eval(pt)>skewborder[1]  ) skewnessFix[kPr] = fitLookUpCleanSkew[kPr]->Eval(pt);
  }

  if (fVariableSkewnes==1){
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[0] ) skewnessFix[kPi] = skewborder[0];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<=skewborder[1] ) skewnessFix[kPi] = skewborder[1];
    if ( fitLookUpCleanSkew[kPi]->Eval(pt)<skewborder[0] && fitLookUpCleanSkew[kPi]->Eval(pt)>skewborder[1]  ) skewnessFix[kPi] = fitLookUpCleanSkew[kPi]->Eval(pt);
    //
    skewnessFix[kKa] = kaSkewnessScan;
    skewnessFix[kPr] = prSkewnessScan;
    skewnessFix[kEl] = 0;
    skewnessFix[kDe] = 0;
  }


  for (Int_t pSpecy=0;pSpecy<nParticle;pSpecy++){
    TF1 *gengaus = new TF1("gengaus",fitFunctionGenGaus,arrMean[pSpecy]-6*arrSigma[pSpecy],arrMean[pSpecy]+6*arrSigma[pSpecy]);
    gengaus->SetParameters(1000,arrMean[pSpecy],arrSigma[pSpecy],2,0);
    Double_t max1 = gengaus->GetMaximumX();
    gengaus->SetParameter(4,skewnessFix[pSpecy]);
    Double_t max2 = gengaus->GetMaximumX();
    arrMeanModif[pSpecy] = arrMean[pSpecy]-(max2-max1);
    delete gengaus;
  }

}

void CreateSplinesFromExpectedsTree(TString expectedsFile)
{
  auto rdfExpecteds = RDataFrame("expecteds", expectedsFile);

  std::vector<TString> particles = {"El", "Pi", "Ka", "Pr", "De"};
  TString centEtaStr = Form("cent_%3.2f_%3.2f_Eta_%3.2f_%3.2f", fCentDown, fCentUp, fEtaDown, fEtaUp);
  TString centEtaFilter = Form("cent >= %f && cent < %f && eta >= %f && eta < %f", fCentDown, fCentUp, fEtaDown, fEtaUp);

  for (UInt_t iPart = 0; iPart < particles.size(); iPart++) {
    hExpectedMean[iPart] = (TH2D*) rdfExpecteds.Filter(centEtaFilter.Data())
                           .Filter(Form("dEdx%s > 0", particles.at(iPart).Data()))
                           .Histo2D({Form("h2Expected_%d_%s", iPart, centEtaStr.Data()), Form("Expected mean %s;ptot;dEdx", particles.at(iPart).Data()),
                                    150, 0.1, 3., 1000, 25., 1000.}, "ptot", Form("dEdx%s", particles.at(iPart).Data())).GetPtr()->Clone();
    hExpectedSigma[iPart] = (TH2D*) rdfExpecteds.Filter(centEtaFilter.Data())
                           .Filter(Form("dEdx%s > 0", particles.at(iPart).Data()))
                           .Histo2D({Form("h2Expected_%d_%s", iPart, centEtaStr.Data()), Form("Expected sigma %s;ptot;dEdx", particles.at(iPart).Data()),
                                    150, 0.1, 3., 1000, 2., 60.}, "ptot", Form("sigma%s", particles.at(iPart).Data())).GetPtr()->Clone();
  }
}