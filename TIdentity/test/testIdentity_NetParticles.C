#include "TIdentity2D.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
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
void      RetrieveMoments(TIdentity2D *tidenObj);
Double_t  EvalFitValue(Int_t particle, Double_t x);
// =======================================================================================================
//
// ======= Modification Part =============================================================================
Int_t       fnSubsamples         = 25;
const Int_t fnParticleBins       = 10;
TString     treeIdentity         = "tracks";
TString     lookUpCloneArrayName = "funcLineShapesCArr";
const Int_t nBinsLineShape       = 1000;
Int_t       fnTestEntries        = 0;
Bool_t      lookUpTableForLine   = kFALSE;
Int_t       lookUpTableLineMode  = 0;  // 0 for hist and 1 for func
//
// fixed tree branches --> [0]=event; [1]=dEdx; [2]=sign; [3]=cutBit; ||||  [4]=cent; [5]=subsampleindex;
Double_t fTreeVariablesArray[6];
const Int_t nBranches = 2;
TString branchNames[nBranches]={"cent","subsampleindex"};
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
Int_t fUsedSign;
const Int_t nMoments = 19;
TFile *fLineShapesLookUpTable = NULL;
TClonesArray *cloneArrFunc=NULL;
TH1D **hLineShape;
TF1 **fLineShape;
TFile *outFile;
TH1D *hDedxDebug;
static TH1D *hDedxDebugLineShapes[fnParticleBins];
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kDe=4,};
//
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
//
int main(int argc, char *argv[])
{
  //  Arguments: $dataTree $lineShapes
  cout << " main.Info: NUMBER OF ARGUMENTS "<<argc<<endl;
  if(argc == 4)
  {
    sprintf(inputfileNameDataTree,"%s",argv[1]);
    sprintf(inputfileNameLineShapes,"%s",argv[2]);
    fUsedSign       = atoi(argv[3]);
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
  //
  TIdentity2D *iden4 = new TIdentity2D(fnParticleBins);      // Set the number of particles to 4
  iden4 -> SetFileName(fileNameDataTree);
  iden4 -> SetBranchNames(nBranches,branchNames);
  iden4 -> SetFunctionPointers(EvalFitValue);
  iden4 -> SetLimits(-50.,50.,3000.,30.,50); // --> (dEdxMin,dEdxMax,nBinsUsed in dEdx), if slice histograms are scaled wrt binwidth, then binwidth=1
  iden4 -> SetUseSign(fUsedSign);  // pass input sign value to TIdentity module
  Long_t nEntries;
  iden4 -> GetTree(nEntries,treeIdentity);
  iden4 -> Reset();
  //
  // track by track loop --> read all track info  and add tracks to the iden4 object
  if (fnTestEntries>0) nEntries = fnTestEntries;
  Int_t countEntry = 0;
  for( Int_t i = 0; i < nEntries; i++ )
  {
    //
    // Read the entries and corresponding bins
    if( !iden4 ->  GetEntry(i) ) continue;
    // reads identity tree and retrives mybin[] info
    iden4 -> GetBins(nBranches, fTreeVariablesArray);
    // print some info
    if(i%200000 == 0) {
      cout << " main.Info: track " << i << " of " << nEntries << "   " << fTreeVariablesArray[2] << "   " << fTreeVariablesArray[5] << endl;
    }
    //
    hDedxDebug->Fill(fTreeVariablesArray[1]);
    iden4 -> AddEntry();
    countEntry++;
  }

  cout << "main.Info: Total number of tracks processed = " << countEntry << endl;
  iden4 -> Finalize();
  iden4  -> AddIntegrals(fUsedSign); // real sign information passed for the check with real data tree
  //
  iden4 -> CalcMoments();
  RetrieveMoments(iden4);
  delete iden4;

  outFile->cd();
  hDedxDebug->Write();
  for (Int_t i=0;i<fnParticleBins;i++){
    hDedxDebugLineShapes[i] = (TH1D*)hLineShape[i]->Clone();
    hDedxDebugLineShapes[i]->SetName(Form("LineShape_part_%d",i));
    hDedxDebugLineShapes[i]->Write();
  }
  outFile -> Close();
  delete outFile;
  return 1;
}
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
void ReadFitParamsFromLineShapes(TString paramTreeName)
{

  cout << " --- In ReadFitParamsFromLineShapes --- " << endl;
  fLineShapesLookUpTable = new TFile(paramTreeName);
  cloneArrFunc   = (TClonesArray*)fLineShapesLookUpTable->Get(lookUpCloneArrayName);
  if (!cloneArrFunc) cout << " ReadFitParamsFromLineShapes.Error: cloneArrFunc is empty " << endl;
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    TString objName = Form("particle_%d",ipart);
    fLineShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
    fLineShape[ipart]->SetName(objName);
    fLineShape[ipart]->SetNpx(nBinsLineShape);
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
  cout << " InitializeObjects.Info: Input sign            = " << fUsedSign                   << endl;
  cout << " InitializeObjects.Info: treeIdentity          = " << treeIdentity                << endl;
  cout << " InitializeObjects.Info: data Tree             = " << inputfileNameDataTree       << endl;
  cout << " InitializeObjects.Info: Line Shapes           = " << inputfileNameLineShapes     << endl;
  cout << " ================================================================================= " << endl;
  //
  outFile = new TFile("TIdentity_Moments_NetParticles.root","recreate");
  hDedxDebug = new TH1D("hDedxDebug","hDedxDebug",3000,-50.,50.);
  //
  // Initialize pointers to lookup table
  fLineShape = new TF1 *[fnParticleBins];
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    fLineShape[ipart] = NULL;
  }
  hLineShape = new TH1D *[fnParticleBins];
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    hLineShape[ipart] = NULL;
  }

}
// =======================================================================================================
void RetrieveMoments(TIdentity2D *tidenObj)
{

  Double_t nnorm = tidenObj->GetMean(kPi)/tidenObj->GetMeanI(kPi);
  cout << " Integrals " << endl;
  for (Int_t i=0;i<fnParticleBins;i++) cout << i << " = " <<  tidenObj -> GetMean(i)/(tidenObj->GetMeanI(i)*nnorm) << endl;
  cout << " First Moments " << endl;
  for (Int_t i=0;i<fnParticleBins;i++) cout << i << " = " <<  tidenObj -> GetMean(i) << endl;
  cout << " Second Moments " << endl;
  for (Int_t i=0;i<fnParticleBins;i++) cout << i << " = " <<  tidenObj -> GetSecondMoment(i) << " " << (tidenObj -> GetMean(i))*(tidenObj -> GetMean(i))+tidenObj -> GetMean(i) << endl;
  cout << " Mixed Moments " << endl;
  for (Int_t i=0;i<fnParticleBins;i++) {
    for (Int_t j=0;j<fnParticleBins;j++) {
      if (i<j) cout << i << "," << j << " = " <<  tidenObj -> GetMixedMoment(i,j) << endl;
    }
  }

}
