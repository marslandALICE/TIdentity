#include "TIdentity2D.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TFile.h"
#include "TVectorF.h"
#include "TSystem.h"
#include "TObject.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include <iomanip>
#include "iostream"
#include "string"

using namespace std;
using std::cout;
using std::setw;

// =======================================================================================================
// Helper Functions
void      InitializeObjects();
void      ReadFitParamsFromLineShapes(TString paramTreeName);
void      RetrieveMoments(TIdentity2D *tidenObj, TVectorF *vecMom, TVectorF *vecInt);
Double_t  EvalFitValue(Int_t particle, Double_t x);
// =======================================================================================================
//
// ======= Modification Part =============================================================================
const Int_t fnParticleBins      = 4;
TString treeIdentity            = "tracks";
const Int_t nBinsLineShape      = 2000;
Bool_t      fTestMode           = kFALSE;
Bool_t      lookUpTableForLine  = kFALSE;
Int_t       lookUpTableLineMode = 0;
// =======================================================================================================
//
// Inputs
Char_t  inputfileNameDataTree[255];     //file name of tree
Char_t  inputfileNameLineShapes[255];   // file name for fit function
TString fileNameDataTree = "";
TString fileNameLineShapes = "";
// =======================================================================================================
//
Double_t nEvents = 0;
Double_t nnorm   = 1.;
const Int_t nMoments = 15;
TVectorF *fIntegrals;
TVectorF *fMmoments;
TFile *fLineShapesLookUpTable = NULL;
TClonesArray *cloneArrFunc=NULL;
TH1D **hLineShape;
TF1 **fLineShape;
enum particles{electron=0, pion=1, kaon=2, proton=3};
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kElEl=4,kPiPi=5,kKaKa=6,kPrPr=7,
  kElPi=8,kElKa=9,kElPr=10,kPiKa=11,kPiPr=12,kKaPr=13,};
//
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
//
int main(int argc, char *argv[])
{
  //  Arguments: $dataTree $lineShapes
  cout << " main.Info: NUMBER OF ARGUMENTS "<<argc<<endl;
  if(argc == 3)
  {
    sprintf(inputfileNameDataTree,"%s",argv[1]);
    sprintf(inputfileNameLineShapes,"%s",argv[2]);
    cout<<" main.Info: read file names from input "<<endl;
  }
  else
  {
    cout<<" main.Error: wrong input list"<<endl;
  }
  //
  InitializeObjects();
  fileNameDataTree   = inputfileNameDataTree;
  fileNameLineShapes = inputfileNameLineShapes;
  //
  // Initialize objects and get the bin information
  TROOT IdentityMethod("IdentityMethod","compiled identity method");
  ReadFitParamsFromLineShapes(fileNameLineShapes);
  //
  // Create the TIdentity2D object and start analysis
  TIdentity2D *iden4 = new TIdentity2D(4);      // Set the number of particles to 4
  iden4 -> SetFileName(fileNameDataTree);
  iden4 -> SetFunctionPointers(EvalFitValue);
  iden4 -> SetLimits(0.,1020.,10.); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
  iden4 -> SetUseSign(0);  // pass input sign value to TIdentity module
  Long_t nEntries;
  iden4 -> GetTree(nEntries,treeIdentity);
  iden4 -> Reset();
  //
  // track by track loop --> read all track info  and add tracks to the iden4 object
  if (fTestMode) nEntries = 10000000;
  for( Int_t i = 0; i < nEntries; i++ )
  {
    if( !iden4 ->  GetEntry(i) ) continue;
    iden4 -> AddEntry();
  }
  iden4 -> Finalize();
  iden4 -> AddIntegrals(0); // real sign information passed for the check with real data tree
  iden4 -> CalcMoments();
  RetrieveMoments(iden4,fMmoments,fIntegrals);
  delete iden4;
  return 1;
}
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
void ReadFitParamsFromLineShapes(TString paramTreeName)
{

  fLineShapesLookUpTable = new TFile(paramTreeName);
  cloneArrFunc   = (TClonesArray*)fLineShapesLookUpTable->Get("funcLineShapesCArr");
  if (!cloneArrFunc) cout << " Error:: cloneArrFunc is empty " << endl;
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    TString objName = Form("particle_%d",ipart);
    fLineShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
    fLineShape[ipart]->SetName(objName);
    hLineShape[ipart] = (TH1D*)fLineShape[ipart]->GetHistogram();
    hLineShape[ipart]->SetName(objName);
  }

}
// =======================================================================================================
Double_t EvalFitValue(Int_t particle, Double_t x)
{

  Int_t bin = hLineShape[particle]->FindBin(x);
  if (lookUpTableLineMode==0) return hLineShape[particle]->GetBinContent(bin);
  if (lookUpTableLineMode==1) return fLineShape[particle]->Eval(x);

}
// =======================================================================================================
void InitializeObjects()
{

  cout << " ================================================================================= " << endl;
  cout << " InitializeObjects.Info: treeIdentity          = " << treeIdentity                << endl;
  cout << " InitializeObjects.Info: data Tree             = " << inputfileNameDataTree       << endl;
  cout << " InitializeObjects.Info: Line Shapes           = " << inputfileNameLineShapes     << endl;
  cout << " ================================================================================= " << endl;
  //
  fMmoments  = new TVectorF(nMoments);
  fIntegrals = new TVectorF(nMoments);
  for(Int_t i=0;i<nMoments; i++){
    (*fMmoments)[i]=0.;
    (*fIntegrals)[i]=0.;
  }
  //
  hLineShape = new TH1D *[fnParticleBins];
  fLineShape = new TF1 *[fnParticleBins];
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++) hLineShape[ipart] = NULL;
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++) fLineShape[ipart] = NULL;

}
// =======================================================================================================
void RetrieveMoments(TIdentity2D *tidenObj, TVectorF *vecMom, TVectorF *vecInt)
{

  (*vecMom)[kEl] = tidenObj -> GetMean(electron);
  (*vecMom)[kPi] = tidenObj -> GetMean(pion);
  (*vecMom)[kKa] = tidenObj -> GetMean(kaon);
  (*vecMom)[kPr] = tidenObj -> GetMean(proton);
  //
  // Second Moments
  (*vecMom)[kElEl] = tidenObj -> GetSecondMoment(electron);
  (*vecMom)[kPiPi] = tidenObj -> GetSecondMoment(pion);
  (*vecMom)[kKaKa] = tidenObj -> GetSecondMoment(kaon);
  (*vecMom)[kPrPr] = tidenObj -> GetSecondMoment(proton);
  //
  // Mixed Moments
  (*vecMom)[kElPi] = tidenObj -> GetMixedMoment(electron,pion);
  (*vecMom)[kElKa] = tidenObj -> GetMixedMoment(electron,kaon);
  (*vecMom)[kElPr] = tidenObj -> GetMixedMoment(electron,proton);
  (*vecMom)[kPiKa] = tidenObj -> GetMixedMoment(pion,kaon);
  (*vecMom)[kPiPr] = tidenObj -> GetMixedMoment(pion,proton);
  (*vecMom)[kKaPr] = tidenObj -> GetMixedMoment(kaon,proton);
  //
  //Integrals:
  (*vecInt)[kEl] = tidenObj -> GetMeanI(electron);
  (*vecInt)[kPi] = tidenObj -> GetMeanI(pion);
  (*vecInt)[kKa] = tidenObj -> GetMeanI(kaon);
  (*vecInt)[kPr] = tidenObj -> GetMeanI(proton);
  //
  // Printing
  nnorm     = (*vecMom)[kPi]/(*vecInt)[kPi];
  nEvents   = tidenObj -> GetNEvents();
  cout << " =============================== Summary of Moments =============================== "<<endl;
  cout << " events      : "<< nEvents << endl;
  cout << " ================================================================================== "<<endl;
  cout << " electron    : "<< (*vecMom)[kEl]   <<" int: "<< (*vecInt)[kEl]*nnorm << "  ratio: " << (*vecMom)[kEl]/((*vecInt)[kEl]*nnorm) << endl;
  cout << " pion        : "<< (*vecMom)[kPi]   <<" int: "<< (*vecInt)[kPi]*nnorm << "  ratio: " << (*vecMom)[kPi]/((*vecInt)[kPi]*nnorm) << endl;
  cout << " kaon        : "<< (*vecMom)[kKa]   <<" int: "<< (*vecInt)[kKa]*nnorm << "  ratio: " << (*vecMom)[kKa]/((*vecInt)[kKa]*nnorm) << endl;
  cout << " proton      : "<< (*vecMom)[kPr]   <<" int: "<< (*vecInt)[kPr]*nnorm << "  ratio: " << (*vecMom)[kPr]/((*vecInt)[kPr]*nnorm) << endl;
  cout << " electron2   : "<< (*vecMom)[kElEl]  <<endl;
  cout << " pion2       : "<< (*vecMom)[kPiPi]  <<endl;
  cout << " kaon2       : "<< (*vecMom)[kKaKa]  <<endl;
  cout << " proton2     : "<< (*vecMom)[kPrPr]  <<endl;

}
