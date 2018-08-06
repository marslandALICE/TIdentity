
//Author: Anar Rustamov
//identity method for third moments in case of 4 particles.

#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <TF1.h>
#include <TH2F.h>
#include <math.h>
#include "TRandom3.h"
#include "TRandom1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include <algorithm>
#include <math.h>
#include <numeric>
#include "TStyle.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLorentzVector.h"

using namespace std;

Double_t RapPion;
Double_t RapKaon;
Double_t RapElectron;

TFile *fitHists = NULL;
Double_t step_M;
const Int_t my_size = 2;

Double_t electronArray[my_size];
Double_t pionArray[my_size];
Double_t kaonArray[my_size];
Double_t protonArray[my_size];

TH1D *hElectron;
TH1D *hPion;
TH1D *hKaon;
TH1D *hProton;

Int_t wrongCount = 0;
Float_t meanKaon;
Float_t meanProton;
Float_t meanPion;
Float_t meanElectron;

Int_t numAllEvents;
Int_t numAllCutEvents;

TF1 *funProton, *funKaon, *funPion, *funElectron, *funSum;
TF1 *funProtonAV, *funKaonAV, *funPionAV, *funElectronAV;

Double_t evtNumber;
Double_t WK_sum, WP_sum, WPi_sum, WEl_sum, WMult_sum;
Int_t multEv;

vector<int> sumMult;

vector<double> W2P_sum_vec;

vector<double> W3P_sum_vec;
vector<double> W3Pi_sum_vec;
vector<double> W3K_sum_vec;

vector<double> W2K_sum_vec;
vector<double> W2Pi_sum_vec;
vector<double> WPK_sum_vec;
vector<double> WPPi_sum_vec;
vector<double> WPiK_sum_vec;
vector<double> WPiPrK_sum_vec;

vector<double> WPr2Pi_sum_vec;
vector<double> WPr2K_sum_vec;
vector<double> WPi2K_sum_vec;

vector<double> WPi2Pr_sum_vec;
vector<double> WK2Pr_sum_vec;
vector<double> WK2Pi_sum_vec;

Double_t twopiroot;

Double_t sigmaProtonF;
Double_t sigmaKaonF;
Double_t sigmaPionF;

TFile *outFile = NULL;

Double_t recP2_av, recPi2_av, recEl2_av, recK2_av, recPPi_av, recPK_av, recPEl_av, recPiK_av, recKEl_av, recPiEl_av;
Double_t wP_P, wK_P, wPi_P, wP_P2, wK_P2, wPi_P2, wP_P3, wK_P3, wPi_P3, wP_K, wK_K, wPi_K, wP_K2, wK_K2, wPi_K2;
Double_t wP_K3, wK_K3, wPi_K3, wP_Pi, wK_Pi, wPi_Pi, wP_Pi2, wK_Pi2, wPi_Pi2, wP_Pi3, wK_Pi3, wPi_Pi3;

///mixed
Double_t wPK_P, wPPi_P, wPiK_P, wPK_K, wPPi_K, wPiK_K, wPK_Pi, wPPi_Pi, wPiK_Pi, wP2Pi_P, wP2Pi_Pi, wP2Pi_K, wPi2P_P, wPi2P_Pi;
Double_t wPi2P_K, wP2K_P, wP2K_Pi, wP2K_K, wK2P_P, wK2P_Pi, wK2P_K, wPi2K_P, wPi2K_Pi, wPi2K_K, wK2Pi_P, wK2Pi_Pi, wK2Pi_K, wPiPrK_Pr, wPiPrK_K, wPiPrK_Pi;
//
TF1 *funP_P = NULL, *funK_P = NULL, *funPi_P = NULL, *funP_P2 = NULL, *funK_P2 = NULL, *funPi_P2 = NULL, *funP_P3 = NULL, *funK_P3 = NULL, *funPi_P3 = NULL;
TF1 *funP_K = NULL, *funK_K = NULL, *funPi_K = NULL, *funP_K2 = NULL, *funK_K2 = NULL, *funPi_K2 = NULL, *funP_K3 = NULL, *funK_K3 = NULL, *funPi_K3 = NULL;
TF1 *funP_Pi = NULL, *funK_Pi = NULL, *funPi_Pi = NULL, *funP_Pi2 = NULL, *funK_Pi2 = NULL, *funPi_Pi2 = NULL, *funP_Pi3 = NULL, *funK_Pi3 = NULL, *funPi_Pi3 = NULL;
TF1 *funPK_P = NULL, *funPPi_P = NULL, *funPiK_P = NULL, *funPK_K = NULL, *funPPi_K = NULL, *funPiK_K = NULL, *funPK_Pi = NULL, *funPPi_Pi = NULL, *funPiK_Pi = NULL;
TF1 *funP2Pi_P = NULL, *funP2Pi_Pi = NULL, *funP2Pi_K = NULL, *funPi2P_P = NULL, *funPi2P_Pi = NULL, *funPi2P_K = NULL, *funP2K_P = NULL, *funP2K_Pi = NULL, *funP2K_K = NULL;
TF1 *funK2P_P = NULL, *funK2P_Pi = NULL, *funK2P_K = NULL, *funPi2K_P = NULL, *funPi2K_Pi = NULL, *funPi2K_K = NULL, *funK2Pi_P = NULL, *funK2Pi_Pi = NULL;
TF1 *funK2Pi_K = NULL, *funPiPrK_P = NULL, *funPiPrK_K = NULL, *funPiPrK_Pi = NULL;
////////////////////////////

Double_t fMin = 0., fMax = 50.;

Int_t WmomB, WptB, WphiB;
Double_t WmomX, WmomY, WmomZ;

Int_t size_size = 3000;

Float_t wmean[3][3];
Float_t wmean2[3][3];
Float_t wmean3[3][3];
Float_t wmix[3][3];

vector<double> WPM_sum_vec;
vector<double> WKM_sum_vec;
vector<double> WPiM_sum_vec;
vector<double> WElM_sum_vec;

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
// ======= Modification Part =============================================================================
const Int_t fnParticleBins = 4;
TString treeIdentity = "tracks";
TString lookUpCloneArrayName = "funcLineShapesCArr";
const Int_t nBinsLineShape = 1000;
Int_t fnTestEntries = 0;
Bool_t lookUpTableForLine = kFALSE;
Int_t lookUpTableLineMode = 0; // 0 for hist and 1 for func
//
Double_t fTreeVariablesArray[5];
const Int_t nBranches = 1;
TString branchNames[nBranches] = {"cent"};

// Some Global declarations
Double_t nEvents = 0;
TFile *fLineShapesLookUpTable = NULL;
TClonesArray *cloneArrFunc = NULL;
TH1D **hLineShape;
TF1 **fLineShape;
TH1D *hDedxDebug;
static TH1D *hDedxDebugLineShapes[fnParticleBins];
enum momentType
{
  kEl = 0,
  kPi = 1,
  kKa = 2,
  kPr = 3,
  kDe = 4,
};

void InitializeObjects()
{
  cout << " ================================================================================= " << endl;
  cout << " ================================================================================= " << endl;
  //
  outFile = new TFile("TIdentity_Moments_NetParticles.root", "recreate");
  hDedxDebug = new TH1D("hDedxDebug", "hDedxDebug", 1500, 0., 50.);
  //
  // Initialize pointers to lookup table
  fLineShape = new TF1 *[fnParticleBins];
  for (Int_t ipart = 0; ipart < fnParticleBins; ipart++)
  {
    fLineShape[ipart] = NULL;
  }
  hLineShape = new TH1D *[fnParticleBins];
  for (Int_t ipart = 0; ipart < fnParticleBins; ipart++)
  {
    hLineShape[ipart] = NULL;
  }
}

void ReadFitParamsFromLineShapes(TString paramTreeName)
{
  cout << " --- In ReadFitParamsFromLineShapes --- " << endl;
  fLineShapesLookUpTable = new TFile(paramTreeName);
  cloneArrFunc = (TClonesArray *)fLineShapesLookUpTable->Get(lookUpCloneArrayName);
  if (!cloneArrFunc)
    cout << " ReadFitParamsFromLineShapes.Error: cloneArrFunc is empty " << endl;
  for (Int_t ipart = 0; ipart < fnParticleBins; ipart++)
  {
    TString objName = Form("particle_%d", ipart);
    fLineShape[ipart] = (TF1 *)cloneArrFunc->FindObject(objName);
    fLineShape[ipart]->SetName(objName);
    fLineShape[ipart]->SetNpx(nBinsLineShape);
    hLineShape[ipart] = (TH1D *)fLineShape[ipart]->GetHistogram();
    hLineShape[ipart]->SetName(objName);
  }

  hElectron = hLineShape[0];
  hPion = hLineShape[1];
  hKaon = hLineShape[2];
  hProton = hLineShape[3];
  funElectron = fLineShape[0];
  funPion = fLineShape[1];
  funKaon = fLineShape[2];
  funProton = fLineShape[3];
}
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================

void setpars(TF1 *fun, Double_t a, Double_t b, Double_t c)
{
  fun->SetParameter(0, a);
  fun->SetParameter(1, b);
  fun->SetParameter(2, c);
}

Double_t getW2(Int_t mix, Int_t a, Int_t b, Int_t c, Int_t i, Int_t l)
{

  return (wmix[mix][i] * wmean[c][l] - wmean[a][i] * wmean[b][i] * wmean[c][l]);
}

Double_t myIntegral(TF1 *fun, Double_t min = 0., Double_t max = 50.)
{

  Double_t xx, sum = 0;
  Double_t step = 50. / (2 * size_size);
  for (Int_t i = 1; i < 2 * size_size + 1; i++)
  {
    xx = 0. + step * i;
    sum += fun->Eval(xx);
  }
  return sum * step;
}

Double_t setFunParam(Double_t x, Int_t i)
{
  Int_t my_bin_num = (Int_t)((x - fMin) / step_M);
  Double_t elValue = electronArray[my_bin_num];
  Double_t piValue = pionArray[my_bin_num];
  Double_t kValue = kaonArray[my_bin_num];
  Double_t prValue = protonArray[my_bin_num];

  if (i == 0)
  {
    return elValue;
  }
  else if (i == 1)
    return piValue;
  else if (i == 2)
    return kValue;
  else if (i == 3)
  {
    return prValue;
  }
  else
    return (elValue + piValue + kValue + prValue);
}

void resetArrays()
{
  for (Int_t m = 0; m < my_size; m++)
  {
    electronArray[m] = 0.;
    pionArray[m] = 0.;
    kaonArray[m] = 0.;
    protonArray[m] = 0.;
  }
}

Double_t getFunctions(Double_t *xx, Double_t *par)
{
  Double_t x = xx[0];
  Int_t j = (Int_t)par[0];
  Int_t i = (Int_t)par[1];
  Int_t k = (Int_t)par[2];
  Double_t val[3];
  val[0] = funProton->Eval(x);
  val[1] = funKaon->Eval(x);
  val[2] = funPion->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  Double_t relVal[3];
  if (sumVal < 1e-15)
    return 0.;
  Double_t power = k;
  for (Int_t m = 0; m < 3; m++)
  {
    relVal[m] = val[m] / sumVal;
    relVal[m] = TMath::Power(relVal[m], power);
  }
  return relVal[j] * val[i];
}

Double_t getFunctionsMix(Double_t *xx, Double_t *par)
{
  Double_t x = xx[0];
  Int_t j = (Int_t)par[0];
  Int_t k = (Int_t)par[1];

  Double_t val[3];
  val[0] = funProton->Eval(x);
  val[1] = funKaon->Eval(x);
  val[2] = funPion->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  if (sumVal < 1e-15)
    return 0.;

  Double_t myVal[3];

  myVal[0] = val[0] * val[1] / sumVal / sumVal;
  myVal[1] = val[0] * val[2] / sumVal / sumVal;
  myVal[2] = val[1] * val[2] / sumVal / sumVal;

  return myVal[j] * val[k];
}

Double_t getFunctionsMix2(Double_t *xx, Double_t *par)
{

  Double_t x = xx[0];
  Int_t m = (Int_t)par[0];
  Int_t j = (Int_t)par[1];
  Int_t k = (Int_t)par[2];

  Double_t val[3];
  val[0] = funProton->Eval(x);
  val[1] = funKaon->Eval(x);
  val[2] = funPion->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  if (sumVal < 1e-15)
    return 0.;

  Double_t myVal[3];

  myVal[0] = val[m] * val[m] * val[0] / sumVal / sumVal / sumVal;
  myVal[1] = val[m] * val[m] * val[1] / sumVal / sumVal / sumVal;
  myVal[2] = val[m] * val[m] * val[2] / sumVal / sumVal / sumVal;

  return myVal[j] * val[k];
}

Double_t getFunctionsMix3(Double_t *xx, Double_t *par)
{
  Double_t x = xx[0];
  Int_t k = (Int_t)par[0];
  Double_t val[3];
  val[0] = funProton->Eval(x);
  val[1] = funKaon->Eval(x);
  val[2] = funPion->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  if (sumVal < 1e-15)
    return 0.;
  Double_t myVal;
  myVal = val[0] * val[1] * val[2] / sumVal / sumVal / sumVal;
  return myVal * val[k];
}

void myDelete(TF1 *fun)
{
  if (fun)
  {
    delete fun;
    fun = NULL;
  }
}

void initFunctions()
{

  funP_P = new TF1("funP_P", getFunctions, fMin, fMax, 3);
  setpars(funP_P, 0, 0, 1);
  funK_P = new TF1("funK_P", getFunctions, fMin, fMax, 3);
  setpars(funK_P, 1, 0, 1);
  funPi_P = new TF1("funPi_P", getFunctions, fMin, fMax, 3);
  setpars(funPi_P, 2, 0, 1);
  funP_P2 = new TF1("funP_P2", getFunctions, fMin, fMax, 3);
  setpars(funP_P2, 0, 0, 2);
  funK_P2 = new TF1("funK_P2", getFunctions, fMin, fMax, 3);
  setpars(funK_P2, 1, 0, 2);
  funPi_P2 = new TF1("funPi_P2", getFunctions, fMin, fMax, 3);
  setpars(funPi_P2, 2, 0, 2);
  funP_P3 = new TF1("funP_P3", getFunctions, fMin, fMax, 3);
  setpars(funP_P3, 0, 0, 3);
  funK_P3 = new TF1("funK_P3", getFunctions, fMin, fMax, 3);
  setpars(funK_P3, 1, 0, 3);
  funPi_P3 = new TF1("funPi_P3", getFunctions, fMin, fMax, 3);
  setpars(funPi_P3, 2, 0, 3);
  funP_K = new TF1("funP_K", getFunctions, fMin, fMax, 3);
  setpars(funP_K, 0, 1, 1);
  funK_K = new TF1("funK_K", getFunctions, fMin, fMax, 3);
  setpars(funK_K, 1, 1, 1);
  funPi_K = new TF1("funPi_K", getFunctions, fMin, fMax, 3);
  setpars(funPi_K, 2, 1, 1);
  funP_K2 = new TF1("funP_K2", getFunctions, fMin, fMax, 3);
  setpars(funP_K2, 0, 1, 2);
  funK_K2 = new TF1("funK_K2", getFunctions, fMin, fMax, 3);
  setpars(funK_K2, 1, 1, 2);
  funPi_K2 = new TF1("funPi_K2", getFunctions, fMin, fMax, 3);
  setpars(funPi_K2, 2, 1, 2);
  funP_K3 = new TF1("funP_K3", getFunctions, fMin, fMax, 3);
  setpars(funP_K3, 0, 1, 3);
  funK_K3 = new TF1("funK_K3", getFunctions, fMin, fMax, 3);
  setpars(funK_K3, 1, 1, 3);
  funPi_K3 = new TF1("funPi_K3", getFunctions, fMin, fMax, 3);
  setpars(funPi_K3, 2, 1, 3);
  funP_Pi = new TF1("funP_Pi", getFunctions, fMin, fMax, 3);
  setpars(funP_Pi, 0, 2, 1);
  funK_Pi = new TF1("funK_Pi", getFunctions, fMin, fMax, 3);
  setpars(funK_Pi, 1, 2, 1);
  funPi_Pi = new TF1("funPi_Pi", getFunctions, fMin, fMax, 3);
  setpars(funPi_Pi, 2, 2, 1);
  funP_Pi2 = new TF1("funP_Pi2", getFunctions, fMin, fMax, 3);
  setpars(funP_Pi2, 0, 2, 2);
  funK_Pi2 = new TF1("funK_Pi2", getFunctions, fMin, fMax, 3);
  setpars(funK_Pi2, 1, 2, 2);
  funPi_Pi2 = new TF1("funPi_Pi2", getFunctions, fMin, fMax, 3);
  setpars(funPi_Pi2, 2, 2, 2);
  funP_Pi3 = new TF1("funP_Pi3", getFunctions, fMin, fMax, 3);
  setpars(funP_Pi3, 0, 2, 3);
  funK_Pi3 = new TF1("funK_Pi3", getFunctions, fMin, fMax, 3);
  setpars(funK_Pi3, 1, 2, 3);
  funPi_Pi3 = new TF1("funPi_Pi3", getFunctions, fMin, fMax, 3);
  setpars(funPi_Pi3, 2, 2, 3);
  ////////////////// mixed ones
  funPK_P = new TF1("funPK_P", getFunctionsMix, fMin, fMax, 3);
  setpars(funPK_P, 0, 0, -10);
  funPPi_P = new TF1("funPPi_P", getFunctionsMix, fMin, fMax, 3);
  setpars(funPPi_P, 1, 0, -10);
  funPiK_P = new TF1("funPiK_P", getFunctionsMix, fMin, fMax, 3);
  setpars(funPiK_P, 2, 0, -10);
  funPK_K = new TF1("funPK_K", getFunctionsMix, fMin, fMax, 3);
  setpars(funPK_K, 0, 1, -10);
  funPPi_K = new TF1("funPPi_K", getFunctionsMix, fMin, fMax, 3);
  setpars(funPPi_K, 1, 1, -10);
  funPiK_K = new TF1("funPiK_K", getFunctionsMix, fMin, fMax, 3);
  setpars(funPiK_K, 2, 1, -10);
  funPK_Pi = new TF1("funPK_Pi", getFunctionsMix, fMin, fMax, 3);
  setpars(funPK_Pi, 0, 2, -10);
  funPPi_Pi = new TF1("funPPi_Pi", getFunctionsMix, fMin, fMax, 3);
  setpars(funPPi_Pi, 1, 2, -10);
  funPiK_Pi = new TF1("funPiK_Pi", getFunctionsMix, fMin, fMax, 3);
  setpars(funPiK_Pi, 2, 2, -10);
  funP2Pi_P = new TF1("funP2Pi_P", getFunctionsMix2, fMin, fMax, 3);
  setpars(funP2Pi_P, 0, 2, 0);
  funP2Pi_Pi = new TF1("funP2Pi_Pi", getFunctionsMix2, fMin, fMax, 3);
  setpars(funP2Pi_Pi, 0, 2, 2);
  funP2Pi_K = new TF1("funP2Pi_K", getFunctionsMix2, fMin, fMax, 3);
  setpars(funP2Pi_K, 0, 2, 1);
  funPi2P_P = new TF1("funPi2P_P", getFunctionsMix2, fMin, fMax, 3);
  setpars(funPi2P_P, 2, 0, 0);
  funPi2P_Pi = new TF1("funPi2P_Pi", getFunctionsMix2, fMin, fMax, 3);
  setpars(funPi2P_Pi, 2, 0, 2);
  funPi2P_K = new TF1("funPi2P_K", getFunctionsMix2, fMin, fMax, 3);
  setpars(funPi2P_K, 2, 0, 1);
  funP2K_P = new TF1("funP2K_P", getFunctionsMix2, fMin, fMax, 3);
  setpars(funP2K_P, 0, 1, 0);
  funP2K_Pi = new TF1("funP2K_Pi", getFunctionsMix2, fMin, fMax, 3);
  setpars(funP2K_Pi, 0, 1, 2);
  funP2K_K = new TF1("funP2K_K", getFunctionsMix2, fMin, fMax, 3);
  setpars(funP2K_K, 0, 1, 1);
  funK2P_P = new TF1("funK2P_P", getFunctionsMix2, fMin, fMax, 3);
  setpars(funK2P_P, 1, 0, 0);
  funK2P_Pi = new TF1("funK2P_Pi", getFunctionsMix2, fMin, fMax, 3);
  setpars(funK2P_Pi, 1, 0, 2);
  funK2P_K = new TF1("funK2P_K", getFunctionsMix2, fMin, fMax, 3);
  setpars(funK2P_K, 1, 0, 1);
  funPi2K_P = new TF1("funPi2K_P", getFunctionsMix2, fMin, fMax, 3);
  setpars(funPi2K_P, 2, 1, 0);
  funPi2K_Pi = new TF1("funPi2K_Pi", getFunctionsMix2, fMin, fMax, 3);
  setpars(funPi2K_Pi, 2, 1, 2);
  funPi2K_K = new TF1("funPi2K_K", getFunctionsMix2, fMin, fMax, 3);
  setpars(funPi2K_K, 2, 1, 1);
  funK2Pi_P = new TF1("funK2Pi_P", getFunctionsMix2, fMin, fMax, 3);
  setpars(funK2Pi_P, 1, 2, 0);
  funK2Pi_Pi = new TF1("funK2Pi_Pi", getFunctionsMix2, fMin, fMax, 3);
  setpars(funK2Pi_Pi, 1, 2, 2);
  funK2Pi_K = new TF1("funK2Pi_K", getFunctionsMix2, fMin, fMax, 3);
  setpars(funK2Pi_K, 1, 2, 1);
  funPiPrK_P = new TF1("funPiPrK_P", getFunctionsMix3, fMin, fMax, 3);
  setpars(funPiPrK_P, 0, -100, -100);
  funPiPrK_K = new TF1("funPiPrK_K", getFunctionsMix3, fMin, fMax, 3);
  setpars(funPiPrK_K, 1, -100, -100);
  funPiPrK_Pi = new TF1("funPiPrK_Pi", getFunctionsMix3, fMin, fMax, 3);
  setpars(funPiPrK_Pi, 2, -100, -100);
}

void calcIntegrals(TF1 *funP, TF1 *funK, TF1 *funPi, TF1 *funEl, Float_t *normme)
{
  Double_t nkNorm = normme[2];  /// * 0.01;
  Double_t npNorm = normme[0];  // * 0.01;
  Double_t npiNorm = normme[1]; /// * 0.01;
  Double_t nelNorm = normme[3]; //// * 0.01;

  if (nkNorm == 0)
    nkNorm = 1;
  if (npNorm == 0)
    npNorm = 1;
  if (npiNorm == 0)
    npiNorm = 1;
  if (nelNorm == 0)
    nelNorm = 1;

  wP_P += myIntegral(funP_P, fMin, fMax) / npNorm;
  wK_P += myIntegral(funK_P, fMin, fMax) / npNorm;
  wPi_P += myIntegral(funPi_P, fMin, fMax) / npNorm;
  wP_P2 += myIntegral(funP_P2, fMin, fMax) / npNorm;
  wK_P2 += myIntegral(funK_P2, fMin, fMax) / npNorm;
  wPi_P2 += myIntegral(funPi_P2, fMin, fMax) / npNorm;

  wP_P3 += myIntegral(funP_P3, fMin, fMax) / npNorm;
  wK_P3 += myIntegral(funK_P3, fMin, fMax) / npNorm;
  wPi_P3 += myIntegral(funPi_P3, fMin, fMax) / npNorm;

  wP_K += myIntegral(funP_K, fMin, fMax) / nkNorm;
  wK_K += myIntegral(funK_K, fMin, fMax) / nkNorm;
  wPi_K += myIntegral(funPi_K, fMin, fMax) / nkNorm;
  wP_K2 += myIntegral(funP_K2, fMin, fMax) / nkNorm;
  wK_K2 += myIntegral(funK_K2, fMin, fMax) / nkNorm;
  wPi_K2 += myIntegral(funPi_K2, fMin, fMax) / nkNorm;

  wP_K3 += myIntegral(funP_K3, fMin, fMax) / nkNorm;
  wK_K3 += myIntegral(funK_K3, fMin, fMax) / nkNorm;
  wPi_K3 += myIntegral(funPi_K3, fMin, fMax) / nkNorm;

  wP_Pi += myIntegral(funP_Pi, fMin, fMax) / npiNorm;
  wK_Pi += myIntegral(funK_Pi, fMin, fMax) / npiNorm;
  wPi_Pi += myIntegral(funPi_Pi, fMin, fMax) / npiNorm;
  wP_Pi2 += myIntegral(funP_Pi2, fMin, fMax) / npiNorm;
  wK_Pi2 += myIntegral(funK_Pi2, fMin, fMax) / npiNorm;
  wPi_Pi2 += myIntegral(funPi_Pi2, fMin, fMax) / npiNorm;

  wP_Pi3 += myIntegral(funP_Pi3, fMin, fMax) / npiNorm;
  wK_Pi3 += myIntegral(funK_Pi3, fMin, fMax) / npiNorm;
  wPi_Pi3 += myIntegral(funPi_Pi3, fMin, fMax) / npiNorm;

  ///mixed
  wPK_P += myIntegral(funPK_P, fMin, fMax) / npNorm;
  wPPi_P += myIntegral(funPPi_P, fMin, fMax) / npNorm;
  wPiK_P += myIntegral(funPiK_P, fMin, fMax) / npNorm;

  wPK_K += myIntegral(funPK_K, fMin, fMax) / nkNorm;
  wPPi_K += myIntegral(funPPi_K, fMin, fMax) / nkNorm;
  wPiK_K += myIntegral(funPiK_K, fMin, fMax) / nkNorm;

  wPK_Pi += myIntegral(funPK_Pi, fMin, fMax) / npiNorm;
  wPPi_Pi += myIntegral(funPPi_Pi, fMin, fMax) / npiNorm;
  wPiK_Pi += myIntegral(funPiK_Pi, fMin, fMax) / npiNorm;

  ///////////////////////////////////////////////////////////////////////

  wP2Pi_P += myIntegral(funP2Pi_P, fMin, fMax) / npNorm;
  wP2Pi_Pi += myIntegral(funP2Pi_Pi, fMin, fMax) / npiNorm;
  wP2Pi_K += myIntegral(funP2Pi_K, fMin, fMax) / nkNorm;

  wPi2P_P += myIntegral(funPi2P_P, fMin, fMax) / npNorm;
  wPi2P_Pi += myIntegral(funPi2P_Pi, fMin, fMax) / npiNorm;
  wPi2P_K += myIntegral(funPi2P_K, fMin, fMax) / nkNorm;

  wP2K_P += myIntegral(funP2K_P, fMin, fMax) / npNorm;
  wP2K_Pi += myIntegral(funP2K_Pi, fMin, fMax) / npiNorm;
  wP2K_K += myIntegral(funP2K_K, fMin, fMax) / nkNorm;

  wK2P_P += myIntegral(funK2P_P, fMin, fMax) / npNorm;
  wK2P_Pi += myIntegral(funK2P_Pi, fMin, fMax) / npiNorm;
  wK2P_K += myIntegral(funK2P_K, fMin, fMax) / nkNorm;

  wPi2K_P += myIntegral(funPi2K_P, fMin, fMax) / npNorm;
  wPi2K_Pi += myIntegral(funPi2K_Pi, fMin, fMax) / npiNorm;
  wPi2K_K += myIntegral(funPi2K_K, fMin, fMax) / nkNorm;

  wK2Pi_P += myIntegral(funK2Pi_P, fMin, fMax) / npNorm;
  wK2Pi_Pi += myIntegral(funK2Pi_Pi, fMin, fMax) / npiNorm;
  wK2Pi_K += myIntegral(funK2Pi_K, fMin, fMax) / nkNorm;

  wPiPrK_Pr += myIntegral(funPiPrK_P, fMin, fMax) / npNorm;
  wPiPrK_K += myIntegral(funPiPrK_K, fMin, fMax) / nkNorm;
  wPiPrK_Pi += myIntegral(funPiPrK_Pi, fMin, fMax) / npiNorm;
}

/////////////////////////////

template <class T>
Float_t Sum(vector<T> &v)
{

  return accumulate(v.begin(), v.end(), 0.0);
}

template <class T>
void setMultipl(T np, T nk, T npi)
{
  funProton->SetParameter(0, np / (twopiroot * sigmaProtonF));
  funKaon->SetParameter(0, nk / (twopiroot * sigmaKaonF));
  funPion->SetParameter(0, npi / (twopiroot * sigmaPionF));
}

template <class T>
void setMultiplAV(T np, T nk, T npi)
{
  funProtonAV->SetParameter(0, np / (twopiroot * sigmaProtonF));
  funKaonAV->SetParameter(0, nk / (twopiroot * sigmaKaonF));
  funPionAV->SetParameter(0, npi / (twopiroot * sigmaPionF));
}

void resetValues()
{
  WP_sum = WK_sum = WPi_sum = WEl_sum = WMult_sum = 0.;
}

void addParticles(Float_t mVal, Int_t &count)
{
  Double_t kValue, prValue, piValue, elValue;
  Double_t WP1, WK1, WPi1, WEl1;
  prValue = funProton->Eval(mVal);
  kValue = funKaon->Eval(mVal);
  piValue = funPion->Eval(mVal);
  elValue = funElectron->Eval(mVal);
  // elValue = 0.;

  Double_t sumValue = prValue + kValue + piValue + elValue;
  if (sumValue == 0)
  {
    cout << "why:   " << mVal << endl;
    wrongCount++;
    WP1 = WK1 = WPi1 = WEl1 = 0;
  }
  else
  {
    WP1 = prValue / sumValue;
    WK1 = kValue / sumValue;
    WPi1 = piValue / sumValue;
    WEl1 = elValue / sumValue;
    // WEl1 = 0; //elValue / sumValue;

    count++;

    WP_sum += WP1;
    WK_sum += WK1;
    WPi_sum += WPi1;
    WEl_sum += WEl1;
    WMult_sum += 1;
  }
}

int identity4Particle()
{

  Char_t eName[255], piName[255], kName[255], prName[255];
  step_M = (fMax - fMin) / (my_size - 1);

  TFile *outputFile = new TFile("aaa.root", "recreate");
  cout << "burada " << endl;
  TFile *dataFile = new TFile("DataTree_Net.root");
  InitializeObjects();
  ReadFitParamsFromLineShapes("LineShapes_Net.root");

  wrongCount = 0;
  TTree *myTree = (TTree *)dataFile->Get("tracks");
  cout << "Funtions are inited" << endl;
  Int_t prevEvt = -1000;
  cout << "set address" << endl;
  //
  // Int_t myBin[4];
  // Double_t myDeDx;
  // Int_t evtNum;
  // myTree->SetBranchAddress("myBin", myBin);
  // myTree->SetBranchAddress("myDeDx", &myDeDx);
  // myTree->SetBranchAddress("evtNum", &evtNum);
  ULong64_t fEventNum;
  Float_t fDEdx;
  UInt_t fCutBit;
  Int_t fSign;
  myTree->SetBranchAddress("gid", &fEventNum);
  myTree->SetBranchAddress("dEdx", &fDEdx);
  myTree->SetBranchAddress("sign", &fSign);
  myTree->SetBranchAddress("cutBit", &fCutBit);
  cout << "setted address" << endl;
  resetValues();
  //
  Int_t count = 0;
  meanKaon = 0;
  meanProton = 0;
  meanPion = 0;
  meanElectron = 0;

  numAllEvents = 0;
  numAllCutEvents = 0;

  Int_t runEvt = (Int_t)myTree->GetEntries();
  cout << "total number " << runEvt << endl;
  Int_t runStat1 = 0;
  Int_t runStat2 = runEvt;

  cout << "running from " << runStat1 << " to " << runStat2 << endl;
  Int_t countVeto = 0;
  Int_t prevEvtVeto = -1;
  Int_t prevEvtPos = -1;
  Bool_t isBinFitted = kFALSE;

  myTree->GetEntry(runStat2 - 1);
  Int_t remEvent = fEventNum;

  for (Int_t i = runStat1; i < runStat2; i++)
  {
    if (i % 500000 == 0)
      cout << "track " << i << " of " << runStat2 - runStat1 << endl;
    myTree->GetEntry(i);
    if ((Int_t)fEventNum == remEvent)
      continue;
    if ((Int_t)fEventNum != prevEvtVeto && prevEvtVeto > 0)
      countVeto++;

    prevEvtVeto = fEventNum;
    prevEvtPos = fEventNum;

    if (fDEdx < fMin || fDEdx > fMax) continue;
    if (!meanKaon)   meanKaon = myIntegral(funKaon); //// * 100.; ////my_par[8]*100.;
    if (!meanProton) meanProton = myIntegral(funProton); //// * 100.; ///my_par[9]*100.;
    if (!meanPion)   meanPion = myIntegral(funPion); /// * 100.; ///my_par[7]*100.;

    numAllEvents += multEv;

    if ((Int_t)fEventNum == prevEvt)
    {
      addParticles(fDEdx, count);
    }
    else
    {
      if (count != 0)
      {
        if (WMult_sum)
        {
          sumMult.push_back(count);
          W2P_sum_vec.push_back(WP_sum * WP_sum);

          W3P_sum_vec.push_back(TMath::Power(WP_sum, 3));
          W3Pi_sum_vec.push_back(TMath::Power(WPi_sum, 3));
          W3K_sum_vec.push_back(TMath::Power(WK_sum, 3));
          W2K_sum_vec.push_back(WK_sum * WK_sum);
          W2Pi_sum_vec.push_back(WPi_sum * WPi_sum);
          //W2El_sum_vec.push_back(WEl_sum*WEl_sum);
          WPM_sum_vec.push_back(WP_sum);
          WKM_sum_vec.push_back(WK_sum);
          WPiM_sum_vec.push_back(WPi_sum);
          // WElM_sum_vec.push_back(WEl_sum);
          WPK_sum_vec.push_back(WP_sum * WK_sum);
          WPPi_sum_vec.push_back(WP_sum * WPi_sum);
          WPiK_sum_vec.push_back(WPi_sum * WK_sum);

          WPiPrK_sum_vec.push_back(WPi_sum * WP_sum * WK_sum);

          WPr2Pi_sum_vec.push_back(TMath::Power(WP_sum, 2) * WPi_sum);
          WPr2K_sum_vec.push_back(TMath::Power(WP_sum, 2) * WK_sum);
          WPi2K_sum_vec.push_back(TMath::Power(WPi_sum, 2) * WK_sum);

          WPi2Pr_sum_vec.push_back(WPi_sum * WPi_sum * WP_sum);
          WK2Pr_sum_vec.push_back(WK_sum * WK_sum * WP_sum);
          WK2Pi_sum_vec.push_back(WK_sum * WK_sum * WPi_sum);
        }
      }
      resetValues();
      count = 0;
      addParticles(fDEdx, count);
      numAllCutEvents++;
    }
    prevEvt = fEventNum;
  } //end event

  Double_t nEvents = countVeto;
  Double_t countVetoBack = countVeto;

  nEvents = sumMult.size();
  countVeto = sumMult.size();

  cout << "# of Events == " << nEvents << endl;
  Double_t corrFactor = multEv / nEvents;
  multEv = countVeto;
  corrFactor = 1;
  cout << "corrFactor == " << corrFactor << endl;

  ////////////////////

  Double_t meanMult = Sum<int>(sumMult) / nEvents;
  cout << "mean mult   " << meanMult << endl;

  Double_t proton_aver = 0;
  Double_t kaon_aver = 0;
  Double_t pion_aver = 0;
  Double_t electron_aver = 0;

  kaon_aver += meanKaon;
  proton_aver += meanProton;
  pion_aver += meanPion;
  electron_aver += meanElectron;

  Double_t proton_aver_fun = proton_aver * corrFactor;
  Double_t kaon_aver_fun = kaon_aver * corrFactor;
  Double_t pion_aver_fun = pion_aver * corrFactor;
  Double_t electron_aver_fun = electron_aver * corrFactor;

  cout << "***********" << endl;
  cout << "***********" << endl;
  cout << "***********" << endl;
  cout << "mean multiplicities : ONLY TRUE FOR ALL STATISTICS !" << endl;

  cout << "proton " << proton_aver * corrFactor << endl;
  cout << "pkaon  " << kaon_aver * corrFactor << endl;
  cout << "pion   " << pion_aver * corrFactor << endl;
  cout << "electron   " << electron_aver * corrFactor << endl;

  proton_aver *= corrFactor;
  kaon_aver *= corrFactor;
  pion_aver *= corrFactor;
  electron_aver *= corrFactor;

  cout << "***********" << endl;
  cout << "***********" << endl;
  cout << "***********" << endl;

  //Double_t averMult = Sum<double>(WMult_sum_vec)/nEvents;
  //Double_t averMult2 = Sum<double>(W2Mult_sum_vec)/nEvents;

  Double_t W2P_aver = Sum<double>(W2P_sum_vec) / nEvents;
  Double_t W2K_aver = Sum<double>(W2K_sum_vec) / nEvents;
  Double_t W2Pi_aver = Sum<double>(W2Pi_sum_vec) / nEvents;

  Double_t WPM_aver = Sum<double>(WPM_sum_vec) / nEvents;
  Double_t WKM_aver = Sum<double>(WKM_sum_vec) / nEvents;
  Double_t WPiM_aver = Sum<double>(WPiM_sum_vec) / nEvents;
  Double_t WElM_aver = Sum<double>(WElM_sum_vec) / nEvents;

  Double_t W3P_aver = Sum<double>(W3P_sum_vec) / nEvents;
  Double_t W3Pi_aver = Sum<double>(W3Pi_sum_vec) / nEvents;
  Double_t W3K_aver = Sum<double>(W3K_sum_vec) / nEvents;
  Double_t WPiPrK_aver = Sum<double>(WPiPrK_sum_vec) / nEvents;
  //    Double_t WPi2P_aver     = Sum<double>(WPi2P_sum_vec)/nEvents;

  Double_t WPr2Pi_aver = Sum<double>(WPr2Pi_sum_vec) / nEvents;
  Double_t WPr2K_aver = Sum<double>(WPr2K_sum_vec) / nEvents;
  Double_t WPi2K_aver = Sum<double>(WPi2K_sum_vec) / nEvents;

  Double_t WPi2Pr_aver = Sum<double>(WPi2Pr_sum_vec) / nEvents;
  Double_t WK2Pr_aver = Sum<double>(WK2Pr_sum_vec) / nEvents;
  Double_t WK2Pi_aver = Sum<double>(WK2Pi_sum_vec) / nEvents;

  Double_t WPK_aver = Sum<double>(WPK_sum_vec) / nEvents;
  Double_t WPPi_aver = Sum<double>(WPPi_sum_vec) / nEvents;
  Double_t WPiK_aver = Sum<double>(WPiK_sum_vec) / nEvents;

  Double_t pr_aver = 0.;  //Sum<int>(P_mult_vec)/nEvents;
  Double_t k_aver = 0.;   ///Sum<int>(K_mult_vec)/nEvents;
  Double_t pi_aver = 0.;  ///Sum<int>(Pi_mult_vec)/nEvents;
  Double_t pr2_aver = 0.; ////Sum<int>(P2_mult_vec)/nEvents;
  Double_t k2_aver = 0.;  ////Sum<int>(K2_mult_vec)/nEvents;
  Double_t pi2_aver = 0.; ///Sum<int>(Pi2_mult_vec)/nEvents;

  Double_t prpi_aver = 0.; ///Sum<int>(PPi_mult_vec)/nEvents;
  Double_t prk_aver = 0.;  ///Sum<int>(PK_mult_vec)/nEvents;
  Double_t pik_aver = 0.;  ////Sum<int>(PiK_mult_vec)/nEvents;

  ///////////////////////////////////////////////////////////////////////
  Double_t pr3_aver = 0.;   //Sum<int>(P3_mult_vec)/nEvents;
  Double_t k3_aver = 0.;    ///Sum<int>(K3_mult_vec)/nEvents;
  Double_t pi3_aver = 0.;   //Sum<int>(Pi3_mult_vec)/nEvents;
  Double_t piprk_aver = 0.; ///Sum<int>(PiPrK_mult_vec)/nEvents;

  Double_t pr2pi_aver = 0.; ////Sum<int>(Pr2Pi_mult_vec)/nEvents;
  Double_t pr2k_aver = 0.;  ///Sum<int>(Pr2K_mult_vec)/nEvents;
  Double_t pi2k_aver = 0.;  ///Sum<int>(Pi2K_mult_vec)/nEvents;

  Double_t pi2pr_aver = 0.; ///Sum<int>(Pi2Pr_mult_vec)/nEvents;
  Double_t k2pr_aver = 0.;  ///Sum<int>(K2Pr_mult_vec)/nEvents;
  Double_t k2pi_aver = 0.;  ///Sum<int>(K2Pi_mult_vec)/nEvents;

  cout << "W2P_aver " << W2P_aver << endl;
  cout << "W2K_aver " << W2K_aver << endl;
  cout << "W2Pi_aver " << W2Pi_aver << endl;
  cout << "WPK_aver " << WPK_aver << endl;
  cout << "WPPi_aver " << WPPi_aver << endl;
  cout << "W2PiK_aver " << WPiK_aver << endl;

  cout << "mean from identities " << endl;
  cout << "WPM_aver " << WPM_aver << " -> " << meanProton << endl;
  cout << "WKM_aver " << WKM_aver << " -> " << meanKaon << endl;
  cout << "WPiM_aver " << WPiM_aver << endl;
  cout << "WElM_aver " << WElM_aver << endl;

  /*
   proton_aver   = WPM_aver;
   kaon_aver     = WKM_aver;
   pion_aver     = WPiM_aver;
   electron_aver = WElM_aver;
   */

  cout << "start to calculate integrals " << endl;
  wP_P = 0., wK_P = 0., wPi_P = 0., wP_P2 = 0., wK_P2 = 0., wPi_P2 = 0., wP_P3 = 0., wK_P3 = 0., wPi_P3 = 0., wP_K = 0., wK_K = 0., wPi_K = 0., wP_K2 = 0., wK_K2 = 0.;
  wPi_K2 = 0., wP_K3 = 0., wK_K3 = 0., wPi_K3 = 0., wP_Pi = 0., wK_Pi = 0., wPi_Pi = 0., wP_Pi2 = 0., wK_Pi2 = 0., wPi_Pi2 = 0., wP_Pi3 = 0., wK_Pi3 = 0., wPi_Pi3 = 0.;
  ///mixed
  wPK_P = 0., wPPi_P = 0., wPiK_P = 0., wPK_K = 0., wPPi_K = 0., wPiK_K = 0., wPK_Pi = 0., wPPi_Pi = 0., wPiK_Pi = 0., wP2Pi_P = 0., wP2Pi_Pi = 0., wP2Pi_K = 0.;
  wPi2P_P = 0., wPi2P_Pi = 0., wPi2P_K = 0., wP2K_P = 0., wP2K_Pi = 0., wP2K_K = 0., wK2P_P = 0., wK2P_Pi = 0., wK2P_K = 0., wPi2K_P = 0., wPi2K_Pi = 0., wPi2K_K = 0.;
  wK2Pi_P = 0., wK2Pi_Pi = 0., wK2Pi_K = 0., wPiPrK_Pr = 0., wPiPrK_K = 0., wPiPrK_Pi = 0.;
  Int_t count_pr = 0;
  //Float_t outInt[27] = {0.};
  Float_t normme[] = {(Float_t)proton_aver, (Float_t)pion_aver, (Float_t)kaon_aver, (Float_t)electron_aver};
  ///////////

  initFunctions();
  calcIntegrals(funProton, funKaon, funPion, funElectron, normme);

  // for (Int_t kk = 0; kk < fitEnt_M; kk++)    // ??????????????
  // {
  //   count_pr++;
  //   fitTree_M->GetEntry(kk);
  //   calcIntegrals(funProton, funKaon, funPion, funElectron, normme);
  //   if (count_pr % 200 == 0)
  //     cout << "still calculating integral " << endl;
  // }

  pr_aver = proton_aver;
  pi_aver = pion_aver;
  k_aver = kaon_aver;

  cout << "test " << wK_P << "  " << wK_K << "  " << wK_Pi << endl;
  cout << "1 " << W2K_aver << endl;
  cout << "2 " << pr_aver * (wK_P2 - wK_P * wK_P) << endl;
  cout << "3 " << pi_aver * (wK_Pi2 - wK_Pi * wK_Pi) << endl;
  cout << "4 " << k_aver * (wK_K2 - wK_K * wK_K) << endl;

  TMatrixD A2(6, 6);
  A2(0, 0) = wP_P * wP_P;
  A2(0, 1) = wP_Pi * wP_Pi;
  A2(0, 2) = wP_K * wP_K;
  A2(0, 3) = 2. * wP_P * wP_Pi;
  A2(0, 4) = 2. * wP_P * wP_K;
  A2(0, 5) = 2. * wP_Pi * wP_K;

  A2(1, 0) = wPi_P * wPi_P;
  A2(1, 1) = wPi_Pi * wPi_Pi;
  A2(1, 2) = wPi_K * wPi_K;
  A2(1, 3) = 2. * wPi_P * wPi_Pi;
  A2(1, 4) = 2. * wPi_P * wPi_K;
  A2(1, 5) = 2. * wPi_Pi * wPi_K;

  A2(2, 0) = wK_P * wK_P;
  A2(2, 1) = wK_Pi * wK_Pi;
  A2(2, 2) = wK_K * wK_K;
  A2(2, 3) = 2. * wK_P * wK_Pi;
  A2(2, 4) = 2. * wK_P * wK_K;
  A2(2, 5) = 2. * wK_Pi * wK_K;

  A2(3, 0) = wP_P * wPi_P;
  A2(3, 1) = wP_Pi * wPi_Pi;
  A2(3, 2) = wP_K * wPi_K;
  A2(3, 3) = wP_P * wPi_Pi + wP_Pi * wPi_P;
  A2(3, 4) = wP_P * wPi_K + wP_K * wPi_P;
  A2(3, 5) = wP_Pi * wPi_K + wP_K * wPi_Pi;

  A2(4, 0) = wP_P * wK_P;
  A2(4, 1) = wP_Pi * wK_Pi;
  A2(4, 2) = wP_K * wK_K;
  A2(4, 3) = wP_P * wK_Pi + wP_Pi * wK_P;
  A2(4, 4) = wP_P * wK_K + wP_K * wK_P;
  A2(4, 5) = wP_Pi * wK_K + wP_K * wK_Pi;

  A2(5, 0) = wPi_P * wK_P;
  A2(5, 1) = wPi_Pi * wK_Pi;
  A2(5, 2) = wPi_K * wK_K;
  A2(5, 3) = wPi_P * wK_Pi + wPi_Pi * wK_P;
  A2(5, 4) = wPi_P * wK_K + wPi_K * wK_P;
  A2(5, 5) = wPi_Pi * wK_K + wPi_K * wK_Pi;

  proton_aver = WPM_aver;
  kaon_aver = WKM_aver;
  pion_aver = WPiM_aver;
  electron_aver = WElM_aver;
  //    }

  pr_aver = proton_aver;
  pi_aver = pion_aver;
  k_aver = kaon_aver;

  Double_t B2[6];
  B2[0] = W2P_aver -
          pr_aver * (wP_P2 - wP_P * wP_P) -
          pi_aver * (wP_Pi2 - wP_Pi * wP_Pi) -
          k_aver * (wP_K2 - wP_K * wP_K);

  B2[1] = W2Pi_aver -
          pr_aver * (wPi_P2 - wPi_P * wPi_P) -
          pi_aver * (wPi_Pi2 - wPi_Pi * wPi_Pi) -
          k_aver * (wPi_K2 - wPi_K * wPi_K);

  B2[2] = W2K_aver -
          pr_aver * (wK_P2 - wK_P * wK_P) -
          pi_aver * (wK_Pi2 - wK_Pi * wK_Pi) -
          k_aver * (wK_K2 - wK_K * wK_K);

  B2[3] = WPPi_aver -
          pr_aver * (wPPi_P - wP_P * wPi_P) -
          pi_aver * (wPPi_Pi - wP_Pi * wPi_Pi) -
          k_aver * (wPPi_K - wP_K * wPi_K);

  B2[4] = WPK_aver -
          pr_aver * (wPK_P - wP_P * wK_P) -
          pi_aver * (wPK_Pi - wP_Pi * wK_Pi) -
          k_aver * (wPK_K - wP_K * wK_K);

  B2[5] = WPiK_aver -
          pr_aver * (wPiK_P - wPi_P * wK_P) -
          pi_aver * (wPiK_Pi - wPi_Pi * wK_Pi) -
          k_aver * (wPiK_K - wPi_K * wK_K);

  cout << "A2(0,0) " << A2(0, 0) << "  " << B2[0] << endl;
  A2.Print();

  TMatrixD invA2 = A2.Invert();
  //   invA.Print();
  Double_t recP2_av = 0.;
  Double_t recPi2_av = 0.;
  Double_t recK2_av = 0.;
  Double_t recPPi_av = 0.;
  Double_t recPK_av = 0.;
  Double_t recPiK_av = 0.;

  //recP2_av = recPi2_av = recK2_av  = 0.;
  //recPPi_av = recPK_av = recPiK_av = 0.;
  for (Int_t tt = 0; tt < 6; tt++)
  {
    recP2_av += invA2(0, tt) * B2[tt];
    recPi2_av += invA2(1, tt) * B2[tt];
    recK2_av += invA2(2, tt) * B2[tt];
    //cout<< "kaon test "<<invA2(2,tt)<<"  "<<B2[tt]<<endl;
    recPPi_av += invA2(3, tt) * B2[tt];
    recPK_av += invA2(4, tt) * B2[tt];
    recPiK_av += invA2(5, tt) * B2[tt];
  }

  cout << "print second moments " << endl;
  cout << "pr2  " << recP2_av  << endl;
  cout << "pi2  " << recPi2_av << endl;
  cout << "k2   " << recK2_av  << endl;
  cout << "prpi " << recPPi_av << endl;
  cout << "prk  " << recPK_av  << endl;
  cout << "pik  " << recPiK_av << endl;

  pr2_aver = recP2_av;
  pi2_aver = recPi2_av;
  k2_aver = recK2_av;
  prpi_aver = recPPi_av;
  prk_aver = recPK_av;
  pik_aver = recPiK_av;

  TMatrixD A(10, 10);
  TMatrixD BB(10, 9);

  wmean[0][0] = wP_P;
  wmean[0][1] = wP_Pi;
  wmean[0][2] = wP_K;
  wmean[1][0] = wPi_P;
  wmean[1][1] = wPi_Pi;
  wmean[1][2] = wPi_K;
  wmean[2][0] = wK_P;
  wmean[2][1] = wK_Pi;
  wmean[2][2] = wK_K;

  wmean2[0][0] = wP_P2;
  wmean2[0][1] = wP_Pi2;
  wmean2[0][2] = wP_K2;

  wmean2[1][0] = wPi_P2;
  wmean2[1][1] = wPi_Pi2;
  wmean2[1][2] = wPi_K2;

  wmean2[2][0] = wK_P2;
  wmean2[2][1] = wK_Pi2;
  wmean2[2][2] = wK_K2;

  wmean3[0][0] = wP_P3;
  wmean3[0][1] = wP_Pi3;
  wmean3[0][2] = wP_K3;
  wmean3[1][0] = wPi_P3;
  wmean3[1][1] = wPi_Pi3;
  wmean3[1][2] = wPi_K3;
  wmean3[2][0] = wK_P3;
  wmean3[2][1] = wK_Pi3;
  wmean3[2][2] = wK_K3;

  wmix[0][0] = wPPi_P;
  wmix[0][1] = wPPi_Pi;
  wmix[0][2] = wPPi_K;

  wmix[1][0] = wPK_P;
  wmix[1][1] = wPK_Pi;
  wmix[1][2] = wPK_K;

  wmix[2][0] = wPiK_P;
  wmix[2][1] = wPiK_Pi;
  wmix[2][2] = wPiK_K;

  Float_t wmix2A[3][3];

  wmix2A[0][0] = wP2Pi_P;
  wmix2A[0][1] = wP2Pi_Pi;
  wmix2A[0][2] = wP2Pi_K;

  wmix2A[1][0] = wP2K_P;
  wmix2A[1][1] = wP2K_Pi;
  wmix2A[1][2] = wP2K_K;

  wmix2A[2][0] = wPi2K_P;
  wmix2A[2][1] = wPi2K_Pi;
  wmix2A[2][2] = wPi2K_K;

  Float_t wmix2B[3][3];

  wmix2B[0][0] = wPi2P_P;
  wmix2B[0][1] = wPi2P_Pi;
  wmix2B[0][2] = wPi2P_K;

  wmix2B[1][0] = wK2P_P;
  wmix2B[1][1] = wK2P_Pi;
  wmix2B[1][2] = wK2P_K;

  wmix2B[2][0] = wK2Pi_P;
  wmix2B[2][1] = wK2Pi_Pi;
  wmix2B[2][2] = wK2Pi_K;

  Float_t wmeanPrPiK[3];

  wmeanPrPiK[0] = wPiPrK_Pr;
  wmeanPrPiK[1] = wPiPrK_Pi;
  wmeanPrPiK[2] = wPiPrK_K;

  Int_t npart = 3;
  Int_t nn = 0;
  for (Int_t p = 0; p < npart; p++)
  {
    nn = 0;
    for (Int_t i = 0; i < npart; i++)
    {
      A(p, nn) = TMath::Power(wmean[p][i], 3);
      nn++;
    }

    for (Int_t i = 0; i < (npart - 1); i++)
      for (Int_t l = i + 1; l < npart; l++)
      {
        A(p, nn) = 3. * TMath::Power(wmean[p][i], 2) * wmean[p][l];
        nn++;
      }

    for (Int_t i = 0; i < (npart - 1); i++)
      for (Int_t l = i + 1; l < npart; l++)
      {
        A(p, nn) = 3. * TMath::Power(wmean[p][l], 2) * wmean[p][i];
        nn++;
      }

    A(p, nn) = 6. * wmean[p][0] * wmean[p][1] * wmean[p][2];
    nn++;
  }

  //cout<<"nn == "<<nn<<endl;
  for (Int_t p = 0; p < npart; p++)
  {
    nn = 0;
    for (Int_t i = 0; i < npart; i++)
    {
      BB(p, nn) = 3. * (wmean2[p][i] * wmean[p][i] - TMath::Power(wmean[p][i], 3));
      nn++;
    }
    for (Int_t i = 0; i < (npart - 1); i++)
      for (Int_t l = i + 1; l < npart; l++)
      {
        BB(p, nn) = 3. * (wmean2[p][i] * wmean[p][l] + wmean2[p][l] * wmean[p][i] - TMath::Power(wmean[p][i], 2) * wmean[p][l] - TMath::Power(wmean[p][l], 2) * wmean[p][i]);
        nn++;
      }
    for (Int_t i = 0; i < npart; i++)
    {
      BB(p, nn) = 2. * TMath::Power(wmean[p][i], 3) + wmean3[p][i] - 3. * wmean2[p][i] * wmean[p][i];
      nn++;
    }
  }
  //cout<<"nn == "<<nn<<endl;

  //////////////////////
  //////////////////////
  //////////////////////

  nn = 0;
  for (Int_t i = 0; i < 3; i++)
  {
    A(3, nn) = wmean[0][i] * wmean[1][i] * wmean[2][i];
    nn++;
  }

  for (Int_t i = 0; i < (npart - 1); i++)
    for (Int_t l = i + 1; l < npart; l++)
    {
      A(3, nn) = wmean[0][i] * wmean[1][i] * wmean[2][l] + wmean[0][i] * wmean[1][l] * wmean[2][i] + wmean[0][l] * wmean[1][i] * wmean[2][i];
      nn++;
    }

  for (Int_t i = 0; i < (npart - 1); i++)
    for (Int_t l = i + 1; l < npart; l++)
    {
      A(3, nn) = wmean[0][l] * wmean[1][l] * wmean[2][i] + wmean[0][l] * wmean[1][i] * wmean[2][l] + wmean[0][i] * wmean[1][l] * wmean[2][l];
      nn++;
    }

  A(3, nn) = wmean[0][0] * wmean[1][1] * wmean[2][2] + wmean[0][1] * wmean[1][2] * wmean[2][0] + wmean[0][2] * wmean[1][0] * wmean[2][1] +
             wmean[0][0] * wmean[1][2] * wmean[2][1] + wmean[0][1] * wmean[1][0] * wmean[2][2] + wmean[0][2] * wmean[1][1] * wmean[2][0];

  nn = 0;
  Int_t p = 3;
  for (Int_t i = 0; i < npart; i++)
  {
    BB(3, nn) = wmix[0][i] * wmean[2][i] + wmix[1][i] * wmean[1][i] + wmix[2][i] * wmean[0][i] - 3. * wmean[0][i] * wmean[1][i] * wmean[2][i];
    nn++;
  }

  for (Int_t i = 0; i < (npart - 1); i++)
    for (Int_t l = i + 1; l < npart; l++)
    {
      BB(3, nn) = getW2(0, 0, 1, 2, i, l) + getW2(0, 0, 1, 2, l, i) + getW2(1, 0, 2, 1, i, l) + getW2(1, 0, 2, 1, l, i) +
                  getW2(2, 1, 2, 0, i, l) + getW2(2, 1, 2, 0, l, i);
      nn++;
    }

  for (Int_t i = 0; i < npart; i++)
  {
    BB(p, nn) = wmeanPrPiK[i] + 2. * wmean[0][i] * wmean[1][i] * wmean[2][i] - wmix[0][i] * wmean[2][i] - wmix[1][i] * wmean[1][i] - wmix[2][i] * wmean[0][i];
    nn++;
  }

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////

  Int_t mm = 4;
  for (Int_t p = 0; p < (npart - 1); p++)
    for (Int_t q = p + 1; q < npart; q++)
    {
      nn = 0;
      for (Int_t i = 0; i < 3; i++)
      {
        A(mm, nn) = TMath::Power(wmean[p][i], 2) * wmean[q][i];
        nn++;
        // cout<<"mm "<<mm<<endl;
      }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          A(mm, nn) = TMath::Power(wmean[p][i], 2) * wmean[q][l] + 2. * wmean[p][i] * wmean[p][l] * wmean[q][i];
          nn++;
        }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          A(mm, nn) = TMath::Power(wmean[p][l], 2) * wmean[q][i] + 2. * wmean[p][i] * wmean[p][l] * wmean[q][l];
          nn++;
        }

      A(mm, nn) = 2. * (wmean[p][0] * wmean[p][1] * wmean[q][2] + wmean[p][0] * wmean[p][2] * wmean[q][1] + wmean[p][1] * wmean[p][2] * wmean[q][0]);

      nn = 0;
      for (Int_t i = 0; i < 3; i++)
      {

        BB(mm, nn) = wmean2[p][i] * wmean[q][i] - 3. * TMath::Power(wmean[p][i], 2) * wmean[q][i] + 2. * wmix[mm - 4][i] * wmean[p][i];
        nn++;
      }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          BB(mm, nn) = (2. * getW2(mm - 4, p, q, p, i, l) + 2. * getW2(mm - 4, p, q, p, l, i) + wmean2[p][i] * wmean[q][l] + wmean2[p][l] * wmean[q][i] -
                        TMath::Power(wmean[p][i], 2) * wmean[q][l] - TMath::Power(wmean[p][l], 2) * wmean[q][i]);
          nn++;
        }

      for (Int_t i = 0; i < 3; i++)
      {
        BB(mm, nn) = 2. * TMath::Power(wmean[p][i], 2) * wmean[q][i] + wmix2A[mm - 4][i] - wmean2[p][i] * wmean[q][i] - 2. * wmix[mm - 4][i] * wmean[p][i];
        nn++;
      }
      mm++;
    }
  ////////////////////////////////////////////////////

  mm = 7;
  for (Int_t q = 0; q < (npart - 1); q++)
    for (Int_t p = q + 1; p < npart; p++)
    {
      nn = 0;
      for (Int_t i = 0; i < 3; i++)
      {
        A(mm, nn) = TMath::Power(wmean[p][i], 2) * wmean[q][i];
        nn++;
        // cout<<"mm "<<mm<<endl;
      }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          A(mm, nn) = TMath::Power(wmean[p][i], 2) * wmean[q][l] + 2. * wmean[p][i] * wmean[p][l] * wmean[q][i];
          nn++;
        }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          A(mm, nn) = TMath::Power(wmean[p][l], 2) * wmean[q][i] + 2. * wmean[p][i] * wmean[p][l] * wmean[q][l];
          nn++;
        }

      A(mm, nn) = 2. * (wmean[p][0] * wmean[p][1] * wmean[q][2] + wmean[p][0] * wmean[p][2] * wmean[q][1] + wmean[p][1] * wmean[p][2] * wmean[q][0]);

      nn = 0;
      for (Int_t i = 0; i < 3; i++)
      {

        BB(mm, nn) = wmean2[p][i] * wmean[q][i] - 3. * TMath::Power(wmean[p][i], 2) * wmean[q][i] + 2. * wmix[mm - 7][i] * wmean[p][i];
        nn++;
      }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          // cout << "nn == " << mm << "  " << nn << endl;
          BB(mm, nn) = (2. * getW2(mm - 7, p, q, p, i, l) + 2. * getW2(mm - 7, p, q, p, l, i) + wmean2[p][i] * wmean[q][l] + wmean2[p][l] * wmean[q][i] - TMath::Power(wmean[p][i], 2) * wmean[q][l] - TMath::Power(wmean[p][l], 2) * wmean[q][i]);
          nn++;
        }

      for (Int_t i = 0; i < 3; i++)
      {
        BB(mm, nn) = 2. * TMath::Power(wmean[p][i], 2) * wmean[q][i] + wmix2B[mm - 7][i] - wmean2[p][i] * wmean[q][i] - 2. * wmix[mm - 7][i] * wmean[p][i];
        nn++;
      }
      mm++;
    }

  Double_t B[10] = {W3P_aver, W3Pi_aver, W3K_aver, WPiPrK_aver, WPr2Pi_aver,
                    WPr2K_aver, WPi2K_aver, WPi2Pr_aver, WK2Pr_aver, WK2Pi_aver};

  for (Int_t i = 0; i < 10; i++)
  {
    B[i] -= (BB(i, 0) * pr2_aver + BB(i, 1) * pi2_aver +
             BB(i, 2) * k2_aver + BB(i, 3) * prpi_aver +
             BB(i, 4) * prk_aver + BB(i, 5) * pik_aver +
             BB(i, 6) * pr_aver + BB(i, 7) * pi_aver + BB(i, 8) * k_aver);
  }

  TString mom3Names[] = {"NP3", "NPi3", "NK3", "Pr2Pi", "Pr2K",
                         "Pi2K", "Pi2Pr", "K2Pr", "K2Pi", "PiPrK"};

  TMatrixD invA = A.Invert();
  Double_t rec3mom[10] = {0};
  cout << "print third moments " << endl;
  for (Int_t i = 0; i < 10; i++)
  {
    for (Int_t tt = 0; tt < 10; tt++)
    {
      rec3mom[i] += invA(i, tt) * B[tt];
    }
      // N^3+3N^2+N --> poisson 3rd moment
    cout << mom3Names[i].Data() << "   :  "  << rec3mom[i] << endl;
  }

  //cout<<"rec3mom 0 "<<rec3mom[0]<<endl;
  Double_t skew = rec3mom[0] - 3. * pr2_aver * proton_aver + 2. * TMath::Power(proton_aver, 3);

  // cout<<"third variance "<<skew/proton_aver<<endl;
  cout << " =========================== " << endl;
  cout << "ratio to poisson 3rd moment  " << (TMath::Power(proton_aver, 3) + 3. * TMath::Power(proton_aver, 2) + proton_aver) / rec3mom[0] << endl;
  cout << "ratio to poisson 2nd moment  " << (TMath::Power(proton_aver, 2) + proton_aver) / pr2_aver << endl;
  cout << " =========================== " << endl;

  skew /= (pr2_aver - proton_aver * proton_aver);
  cout << "skew === " << skew << endl;


  outputFile->cd();
  outputFile->Close();

  Char_t outFileName[255];

  return 0;
}

int main(int argc, char *argv[])
{
  cout << "***********************************************" << endl;
  cout << "***********************************************" << endl;
  cout << "******************** IDENTITY METHOD **********" << endl;
  cout << "***********************************************" << endl;
  cout << "***********************************************" << endl;
  cout << "      " << endl;
  cout << "      " << endl;
  cout << "      " << endl;
  TROOT IdentityMethod("IdentityMethod", "compiled identity method");
  return identity4Particle();
}
