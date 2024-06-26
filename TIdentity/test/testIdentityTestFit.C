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
void      RetrieveMoments(TIdentity2D *tidenObj, TVectorF *vecMom, TVectorF *vecInt);
Double_t  EvalFitValue(Int_t particle, Double_t x);
// =======================================================================================================
//
// ======= Modification Part =============================================================================
Int_t       fnSubsamples         = 25;
const Int_t fnParticleBins       = 5;
Int_t       fnSignBins           = 3;
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
Double_t nnorm   = 1.;
Int_t fUsedSign;
Int_t fSubsampleIndex;
Int_t fSignBin;
const Int_t nMoments = 19;
TVectorF *fIntegrals;
TVectorF *fMmoments;
TTree *momTree=NULL;
TFile *fLineShapesLookUpTable = NULL;
TClonesArray *cloneArrFunc=NULL;
TH1D ***hLineShape;
TF1 ***fLineShape;
TFile *outFile;
TH1D *hDedxDebug[2];
static TH1D *hDedxDebugLineShapes[2][fnParticleBins];
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kDe=4,
  kElEl=5,kPiPi=6,kKaKa=7,kPrPr=8,kDeDe=9,
  kElPi=10,kElKa=11,kElPr=12,kElDe=13,
  kPiKa=14,kPiPr=15,kPiDe=16,
  kKaPr=17,kKaDe=18,};
TString momNames[nMoments] = {"#LT e #GT","#LT #pi #GT","#LT K #GT","#LT p #GT","#LT d #GT",
  "#LT e^{2} #GT","#LT #pi^{2} #GT","#LT K^{2} #GT","#LT p^{2} #GT","#LT d^{2} #GT",
  "#LT e,#pi #GT","#LT e,K #GT","#LT e,p #GT","#LT e,d #GT",
  "#LT #pi,K #GT","#LT #pi,p #GT","#LT #pi,d #GT",
  "#LT K,p #GT","#LT K,d #GT"};

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
  for (Int_t iss = 0; iss<fnSubsamples; iss++){

    fSubsampleIndex = iss;
    TIdentity2D *iden4 = new TIdentity2D(fnParticleBins);      // Set the number of particles to 4
    iden4 -> SetFileName(fileNameDataTree);
    iden4 -> SetBranchNames(nBranches,branchNames);
    iden4 -> SetFunctionPointers(EvalFitValue);
    iden4 -> SetLimits(0.,50.,1500.,30.,50); // --> (dEdxMin,dEdxMax,nBinsUsed in dEdx), if slice histograms are scaled wrt binwidth, then binwidth=1
    iden4 -> SetUseSign(fUsedSign);  // pass input sign value to TIdentity module
    Long_t nEntries;
    iden4 -> GetTree(nEntries,treeIdentity);
    iden4 -> Reset();
    //
    // track by track loop --> read all track info  and add tracks to the iden4 object
    if (fnTestEntries>0) nEntries = fnTestEntries;
    Int_t fUsedBins[fnSignBins]={0};
    Int_t countEntry = 0;
    for( Int_t i = 0; i < nEntries; i++ )
    {
      //
      // Read the entries and corresponding bins
      if( !iden4 ->  GetEntry(i) ) continue;
      iden4 -> GetBins(nBranches, fTreeVariablesArray);    // reads identity tree and retrives mybin[] info
      if(i%20000000 == 0) {
        cout << " main.Info: track " << i << " of " << nEntries;
        cout << " -- used sign  =  " << fSubsampleIndex << "  ----  " << fSignBin << "   " << fTreeVariablesArray[2] << "   " << fTreeVariablesArray[5] << endl;
      }
      //
      if (fTreeVariablesArray[2]==-1) { fSignBin=0; hDedxDebug[0]->Fill(fTreeVariablesArray[1]); }
      if (fTreeVariablesArray[2]== 1) { fSignBin=1; hDedxDebug[1]->Fill(fTreeVariablesArray[1]); }
      //
      //
      if(fTreeVariablesArray[5]!=fSubsampleIndex) continue;
      //
      // tag the bins an read only required bins
      if (fUsedSign==0) { fUsedBins[fSignBin] = 1; iden4 -> AddEntry(); countEntry++; }
      else { fUsedBins[fSignBin] = 1; iden4 -> AddEntry(); countEntry++; }

    }

    cout << "main.Info: Total number of tracks processed = " << countEntry << endl;
    iden4 -> Finalize();
    //
    // Calculate 2. order moments only for full range
    for(Int_t isign = 0; isign < fnSignBins-1; isign++) {
      if(fUsedBins[isign] != 1) continue;
      cout << " main.Info: used sign Bin  =  " << isign << endl;
      fSignBin = isign;     // to be used in retrival of obj from the lookup table
      iden4  -> AddIntegrals(fUsedSign); // real sign information passed for the check with real data tree
    }
    //
    iden4 -> CalcMoments();
    RetrieveMoments(iden4,fMmoments,fIntegrals);
    momTree -> Fill();
    delete iden4;
  }
  outFile->cd();
  momTree -> Write();
  for (Int_t isign=0;isign<2;isign++){
    hDedxDebug[isign]->Write();
    for (Int_t i=0;i<fnParticleBins;i++){
      hDedxDebugLineShapes[isign][i] = (TH1D*)hLineShape[i][isign]->Clone();
      hDedxDebugLineShapes[isign][i]->SetName(Form("LineShape_sign_%d_part_%d",isign,i));
      hDedxDebugLineShapes[isign][i]->Write();
    }
  }

  outFile -> Close();
  delete outFile; //yeni eklave etdim.
  // delete iden4;
  return 1;
}
// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
void ReadFitParamsFromLineShapes(TString paramTreeName)
{

  fLineShapesLookUpTable = new TFile(paramTreeName);
  cloneArrFunc   = (TClonesArray*)fLineShapesLookUpTable->Get(lookUpCloneArrayName);
  if (!cloneArrFunc) cout << " ReadFitParamsFromLineShapes.Error: cloneArrFunc is empty " << endl;
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    for (Int_t isign = 0; isign<fnSignBins; isign++){
      TString objName = Form("particle_%d_bin_%d",ipart,isign);
      fLineShape[ipart][isign] = (TF1*)cloneArrFunc->FindObject(objName);
      fLineShape[ipart][isign]->SetName(objName);
      fLineShape[ipart][isign]->SetNpx(nBinsLineShape);
      hLineShape[ipart][isign] = (TH1D*)fLineShape[ipart][isign]->GetHistogram();
      hLineShape[ipart][isign]->SetName(objName);
    }
  }

}
// =======================================================================================================
Double_t EvalFitValue(Int_t particle, Double_t x)
{

  Int_t bin = hLineShape[particle][fSignBin]->FindBin(x);
  if (lookUpTableLineMode==0) return hLineShape[particle][fSignBin]->GetBinContent(bin);
  if (lookUpTableLineMode==1) return fLineShape[particle][fSignBin]->Eval(x);

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
  outFile = new TFile("TIdentity_Moments.root","recreate");
  for (Int_t i=0; i<2; i++) hDedxDebug[i] = new TH1D(Form("hDedxDebug_%d",i),Form("hDedxDebug_%d",i),1500,0.,50.);
  //
  fMmoments  = new TVectorF(nMoments);
  fIntegrals = new TVectorF(nMoments);
  for(Int_t i=0;i<nMoments; i++){
    (*fMmoments)[i]=0.;
    (*fIntegrals)[i]=0.;
  }
  //
  // Initialize pointers to lookup table
  fLineShape = new TF1 **[fnParticleBins];
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    fLineShape[ipart] = new TF1*[fnSignBins];
    for (Int_t isign = 0; isign<fnSignBins; isign++){
      fLineShape[ipart][isign] = NULL;
    }
  }
  hLineShape = new TH1D **[fnParticleBins];
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    hLineShape[ipart] = new TH1D*[fnSignBins];
    for (Int_t isign = 0; isign<fnSignBins; isign++){
      hLineShape[ipart][isign] = NULL;
    }
  }
  //
  // initialize output tree
  momTree = new TTree("momTree","momTree");
  momTree -> Branch("sign",&fUsedSign);
  momTree -> Branch("subsampleindex",&fSubsampleIndex);
  momTree -> Branch("moment",&fMmoments);

}
// =======================================================================================================
void RetrieveMoments(TIdentity2D *tidenObj, TVectorF *vecMom, TVectorF *vecInt)
{

  (*vecMom)[kEl] = tidenObj -> GetMean(kEl);
  (*vecMom)[kPi] = tidenObj -> GetMean(kPi);
  (*vecMom)[kKa] = tidenObj -> GetMean(kKa);
  (*vecMom)[kPr] = tidenObj -> GetMean(kPr);
  (*vecMom)[kDe] = tidenObj -> GetMean(kDe);
  //
  // Second Moments
  (*vecMom)[kElEl] = tidenObj -> GetSecondMoment(kEl);
  (*vecMom)[kPiPi] = tidenObj -> GetSecondMoment(kPi);
  (*vecMom)[kKaKa] = tidenObj -> GetSecondMoment(kKa);
  (*vecMom)[kPrPr] = tidenObj -> GetSecondMoment(kPr);
  (*vecMom)[kDeDe] = tidenObj -> GetSecondMoment(kDe);
  //
  // Mixed Moments
  (*vecMom)[kElPi] = tidenObj -> GetMixedMoment(kEl,kPi);
  (*vecMom)[kElKa] = tidenObj -> GetMixedMoment(kEl,kKa);
  (*vecMom)[kElPr] = tidenObj -> GetMixedMoment(kEl,kPr);
  (*vecMom)[kElDe] = tidenObj -> GetMixedMoment(kEl,kDe);

  (*vecMom)[kPiKa] = tidenObj -> GetMixedMoment(kPi,kKa);
  (*vecMom)[kPiPr] = tidenObj -> GetMixedMoment(kPi,kPr);
  (*vecMom)[kPiDe] = tidenObj -> GetMixedMoment(kPi,kDe);

  (*vecMom)[kKaPr] = tidenObj -> GetMixedMoment(kKa,kPr);
  (*vecMom)[kKaDe] = tidenObj -> GetMixedMoment(kKa,kDe);

  //
  //Integrals:
  (*vecInt)[kEl] = tidenObj -> GetMeanI(kEl);
  (*vecInt)[kPi] = tidenObj -> GetMeanI(kPi);
  (*vecInt)[kKa] = tidenObj -> GetMeanI(kKa);
  (*vecInt)[kPr] = tidenObj -> GetMeanI(kPr);
  (*vecInt)[kDe] = tidenObj -> GetMeanI(kDe);

  for (Int_t i=0;i<nMoments;i++) (*vecMom)[i] = (*vecMom)[i]*fnSubsamples;
  //
  // Printing
  nnorm     = (*vecMom)[kPi]/(*vecInt)[kPi];
  nEvents   = tidenObj -> GetNEvents();
  cout << " =============================== Summary of Moments =============================== "<<endl;
  cout << " #Particles  : " << fnParticleBins << endl;
  cout << " events      : " << nEvents+1 << endl;
  cout << " ================================================================================== "<<endl;
  cout << " electron    : "<< (*vecMom)[kEl]   <<" int: "<< (*vecInt)[kEl] << "  ratio: " << (*vecMom)[kEl]/((*vecInt)[kEl]) << endl;
  cout << " pion        : "<< (*vecMom)[kPi]   <<" int: "<< (*vecInt)[kPi] << "  ratio: " << (*vecMom)[kPi]/((*vecInt)[kPi]) << endl;
  cout << " kaon        : "<< (*vecMom)[kKa]   <<" int: "<< (*vecInt)[kKa] << "  ratio: " << (*vecMom)[kKa]/((*vecInt)[kKa]) << endl;
  cout << " proton      : "<< (*vecMom)[kPr]   <<" int: "<< (*vecInt)[kPr] << "  ratio: " << (*vecMom)[kPr]/((*vecInt)[kPr]) << endl;
  cout << " ================================================================================== "<<endl;
  cout << " electron    : "<< (*vecMom)[kEl]   <<" int: "<< (*vecInt)[kEl]*nnorm << "  ratio: " << (*vecMom)[kEl]/((*vecInt)[kEl]*nnorm) << endl;
  cout << " pion        : "<< (*vecMom)[kPi]   <<" int: "<< (*vecInt)[kPi]*nnorm << "  ratio: " << (*vecMom)[kPi]/((*vecInt)[kPi]*nnorm) << endl;
  cout << " kaon        : "<< (*vecMom)[kKa]   <<" int: "<< (*vecInt)[kKa]*nnorm << "  ratio: " << (*vecMom)[kKa]/((*vecInt)[kKa]*nnorm) << endl;
  cout << " proton      : "<< (*vecMom)[kPr]   <<" int: "<< (*vecInt)[kPr]*nnorm << "  ratio: " << (*vecMom)[kPr]/((*vecInt)[kPr]*nnorm) << endl;
  cout << " electron2   : "<< (*vecMom)[kElEl]  <<endl;
  cout << " pion2       : "<< (*vecMom)[kPiPi]  <<endl;
  cout << " kaon2       : "<< (*vecMom)[kKaKa]  <<endl;
  cout << " proton2     : "<< (*vecMom)[kPrPr]  <<endl;

  TFile *fGen = new TFile("ToyMC_Moments.root");
  TH1D *hGen   = (TH1D*)fGen->Get("hGen"); hGen->SetLineColor(kBlack);
  TH1D *hRec   = (TH1D*)hGen->Clone();   hRec->SetName("hRec"); hRec->SetLineColor(kRed+1);
  TH1D *hRatio = (TH1D*)hGen->Clone();   hRatio->SetName("hRatio"); hRatio->GetYaxis()->SetTitle("gen/rec");
  for (Int_t i=1;i<nMoments+1;i++) {
    cout << momNames[i-1] << "  " << (*vecMom)[i-1] << endl;
    hRec->SetBinContent(i,(*vecMom)[i-1]);
  }

  // Read generated histogram
  for (Int_t i=1;i<nMoments+1;i++) {
    Double_t gen = hGen->GetBinContent(i);
    Double_t rec = hRec->GetBinContent(i);
    Double_t ratio = hGen->GetBinContent(i)/hRec->GetBinContent(i);
    cout << std::setw(10) << momNames[i-1] << " -->  gen:   ";
    cout << std::setw(10) << gen           << " -->  ";
    cout << std::setw(10) << rec           << "    ///  ";
    cout << std::setw(10) << ratio << endl;
    hRatio->SetBinContent(i,ratio);
  }

  outFile->cd();
  hRec->Write();
  hGen->Write();
  hRatio->Write();

}
