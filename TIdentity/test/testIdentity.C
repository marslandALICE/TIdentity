#include "TIdentity2D.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
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
TString   PrintNumInBinary(UInt_t num);
Bool_t    ApplyTreeSelection(UInt_t cut, Int_t syst);
Double_t  fitFunctionGenGaus(Double_t *x, Double_t *par);
void      InitializeObjects();
void      ReadFitParamsFromLineShapes(TString paramTreeName);
void      ReadFitParamsFromTree(TString paramTreeName, Int_t fitIter);
Double_t  EvalFitValue(Int_t particle, Double_t x);
void      RetrieveMoments(TIdentity2D *tidenObj, TVectorF *vecMom, TVectorF *vecInt);
void      PrintInitInfo();

//
// =======================================================================================================
//
// enums
enum momentType {
  kEl=0,
  kPi=1,
  kKa=2,
  kPr=3,
  kElEl=4,
  kPiPi=5,
  kKaKa=6,
  kPrPr=7,
  kElPi=8,
  kElKa=9,
  kElPr=10,
  kPiKa=11,
  kPiPr=12,
  kKaPr=13,
};
TString momNames[14] = {"kEl","kPi","kKa","kPr",
  "kElEl","kPiPi","kKaKa","kPrPr",
  "kElPi","kElKa","kElPr",
  "kPiKa","kPiPr",
  "kKaPr"};
enum trackCutBit {
  kCRows60=0,
  kCRows80=1,
  kCRows100=2,
  kChi2TPC3=3,
  kChi2TPC4=4,
  kChi2TPC5=5,
  kDCAXYSmall=6,
  kDCAXY=7,
  kDCAXYLarge=8,
  kVZSmall=9,
  kVZ=10,
  kVZLarge=11,
  kEventVertexZSmall=12,
  kEventVertexZ=13,
  kEventVertexZLarge=14,
  kClusterRequirementITS=15,
  kNewITSCut=16,
};
//
//
// ======= Modification Part =============================================================================
const Int_t fnSignBins     = 2;
const Int_t fnEtaBins      = 16;  // MC: 8, Data:16
const Int_t fnCentBins     = 9;
const Int_t fnParticleBins = 4;
const Int_t fnMomBins      = 150;
//
//
// Look up table related
const Int_t nBinsLineShape      = 2000;
Int_t       fnTestEntries       = 0;
Bool_t      lookUpTableForLine  = kFALSE;
Int_t       lookUpTableLineMode = 0;     // 0 for TH1D, 1 for TF1
//
//
const Float_t fEtaRangeDown = -0.8;
const Float_t fEtaRangeUp   = 0.8;
const Float_t fMomRangeDown = 0.2;
const Float_t fMomRangeUp   = 3.2;
Float_t xCentBins[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
//
// Default braches are --> ULong64_t ("gid"); Float_t ("dEdx"), Int_t ("sign"), UInt_t ("cutBit")
// the rest of them can be given as an array
//
// fixed tree branches --> [0]=event; [1]=dEdx; [2]=sign; [3]=cutBit; ||||  [4]=eta; [5]=cent; [6]=ptot; [7]=cRows; [8]=tpcchi2;
Double_t fTreeVariablesArray[9];
const Int_t nBranches = 5;
TString branchNames[nBranches]={"eta","cent","ptot","cRows","chi2TPC"};
//
//
Int_t   fNthFitIteration = 6;
TString treeLineShapes   = "treeId";
TString treeIdentity     = "fIdenTree";   // for data "fIdenTree",    for MC "fIdenTreeMC"  , new data "tracks"
//
//
// =======================================================================================================
//
// Inputs
Char_t  inputfileNameDataTree[255];     //file name of tree
Char_t  inputfileNameLineShapes[255];   // file name for fit functions
Int_t   fSubsample=-100;
Int_t   fSignInput=-100;
Float_t fCentInput=-100;
Float_t fEtaDownInput=-100;
Float_t fEtaUpInput=-100;
Float_t fpDownInput=-100;
Float_t fpUpInput=-100;
Int_t fSystematic=-100;
//
//
TString fileNameDataTree = "";
TString fileNameLineShapes = "";
Int_t fCentInputBin = 0;
Int_t fpDownBin     = 0;
Int_t fpUpBin       = 0;
Int_t fEtaDownBin   = 0;
Int_t fEtaUpBin     = 0;
//
//
Double_t nEvents = 0;
Double_t nnorm   = 1.;
Int_t fEtaBin, fCentBin, fMomBin, fSignBin;
Int_t fUsedBins[fnEtaBins][fnCentBins][fnMomBins][fnSignBins];
UInt_t fCutBit;
const Int_t nMoments = 14;
TVectorF *fIntegrals, *fIntegralsPos, *fIntegralsNeg;
TVectorF *fMmoments,  *fMmomentsPos,  *fMmomentsNeg;
//
// to be initialized
TH1D *fChi2    = NULL;
TH1D *fcRows   = NULL;
TH1D *fhPtot   = NULL;
TH1D *fhEta    = NULL;
TH1D *fhCent   = NULL;
static TH1D ******hLineShape;
static TF1 ******fLineShape;
TTree *momTree = NULL;
//
// members
TString fitFunctionGenGausStr = "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
TFile *fLineShapesLookUpTable = NULL;
TClonesArray *cloneArrHist=NULL;
TClonesArray *cloneArrFunc=NULL;
TTree *treeLookUp = NULL;
TFile *outFile = NULL;
TF1   *fgengaus = 0;
TStopwatch timer;
//
Double_t fAmpArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins][fnSignBins]; // fAmpArr[fEtaBin][fCentBin][fMomBin][particleType]
Double_t fMeanArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins][fnSignBins];
Double_t fSigmaArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins][fnSignBins];
Double_t fSkewArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins][fnSignBins];
Double_t fKurtosisArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBins][fnSignBins];
//
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
//
int main(int argc, char *argv[])
{
  //  Arguments: $dataTree $lineShapes $fSubsample $fSign $fCent $fpDown $fpUp $fEtaDown $fEtaUp $fSystematic
  cout << " main.Info: NUMBER OF ARGUMENTS "<<argc<<endl;
  if(argc == 11)
  {
    sprintf(inputfileNameDataTree,"%s",argv[1]);
    sprintf(inputfileNameLineShapes,"%s",argv[2]);
    fSubsample       = atoi(argv[3]);
    fSignInput       = atoi(argv[4]);
    fCentInput       = atof(argv[5]);
    fpDownInput      = atof(argv[6]);
    fpUpInput        = atof(argv[7]);
    fEtaDownInput    = atof(argv[8]);
    fEtaUpInput      = atof(argv[9]);
    fSystematic      = atoi(argv[10]);
    cout<<" main.Info: read file names from input "<<endl;
    TString outPutFileNAme = Form("TIMoments_sub%d_cent_%3.2f_mom_%3.2f_%3.2f_eta_%3.2f_%3.2f.root"
    ,fSubsample,fCentInput,fpDownInput,fpUpInput,fEtaDownInput,fEtaUpInput);
    outFile = new TFile(outPutFileNAme,"recreate");
    cout << " main.Info: write output into:    " << outPutFileNAme << endl;
  }
  else
  {
    cout<<" main.Error: wrong input list"<<endl;
  }
  //
  //
  InitializeObjects();
  fileNameDataTree   = inputfileNameDataTree;
  fileNameLineShapes = inputfileNameLineShapes;
  fpDownBin       = fhPtot -> FindBin(fpDownInput   + 0.0000001) - 1;
  fpUpBin         = fhPtot -> FindBin(fpUpInput     - 0.0000001) - 1;
  fEtaDownBin     = fhEta  -> FindBin(fEtaDownInput + 0.0000001) - 1;
  fEtaUpBin       = fhEta  -> FindBin(fEtaUpInput   - 0.0000001) - 1;
  fCentInputBin   = fhCent -> FindBin(fCentInput    + 0.0000001) - 1;
  Int_t etaRange[2] = {fEtaDownBin,fEtaUpBin};
  Int_t momRange[2] = {fpDownBin,fpUpBin};
  //
  // Initialize objects and get the bin information
  TROOT IdentityMethod("IdentityMethod","compiled identity method");
  PrintInitInfo();
  if (lookUpTableForLine) ReadFitParamsFromLineShapes(fileNameLineShapes);
  else  ReadFitParamsFromTree(fileNameLineShapes,fNthFitIteration);
  //
  // Create the TIdentity2D object and start analysis
  fgengaus = new TF1("fgengaus",fitFunctionGenGaus,0,1020,5);
  TIdentity2D *iden4 = new TIdentity2D(fnParticleBins);
  iden4 -> SetBranchNames(nBranches,branchNames);
  iden4 -> SetFileName(fileNameDataTree);
  iden4 -> SetFunctionPointers(EvalFitValue);
  iden4 -> SetLimits(20.,1020.,250.,1500.,10); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
  iden4 -> SetUseSign(fSignInput);  // pass input sign value to TIdentity module
  Long_t nEntries;
  iden4 -> GetTree(nEntries,treeIdentity);
  iden4 -> Reset();
  //
  // track by track loop --> read all track info  and add tracks to the iden4 object
  //
  if (fnTestEntries>0) nEntries = fnTestEntries;
  Int_t countEntry=0;
  for( Int_t i = 0; i < nEntries; i++ )
  {
    //
    if( !iden4 ->  GetEntry(i) ) continue;
    iden4      ->  GetBins(nBranches, fTreeVariablesArray);    // reads identity tree and retrives mybin[] info
    if(i%2000000 == 0) {
      cout << " main.Info: track " << i << " of " << nEntries;
      cout << " -- bin 0 =  " << fTreeVariablesArray[0];
      cout << " -- bin 1 =  " << fTreeVariablesArray[1];
      cout << " -- bin 2 =  " << fTreeVariablesArray[2];
      cout << " -- bin 3 =  " << fTreeVariablesArray[3] << endl;
    }
    //
    // Choose which kind of tree input is used
    if (treeIdentity=="fIdenTree"){
      // Dikkkkaaaaat this is for backward compatibility
      fEtaBin  = Int_t(fTreeVariablesArray[0]);
      fCentBin = Int_t(fTreeVariablesArray[1]);
      fMomBin  = Int_t(fTreeVariablesArray[2]);
      if (Int_t(fTreeVariablesArray[3])==-1) fSignBin=0;
      if (Int_t(fTreeVariablesArray[3])== 1) fSignBin=1;
    }
    if(treeIdentity=="tracks") {
      if (Int_t(fTreeVariablesArray[2])==-1) fSignBin=0;
      if (Int_t(fTreeVariablesArray[2])== 1) fSignBin=1;
      fCutBit  = (UInt_t)fTreeVariablesArray[3];
      fEtaBin  = fhEta  -> FindBin(fTreeVariablesArray[4] + 0.0000001) -1;
      fCentBin = fhCent -> FindBin(fTreeVariablesArray[5] + 0.0000001) -1;
      fMomBin  = fhPtot -> FindBin(fTreeVariablesArray[6] + 0.0000001) -1;
      if (!ApplyTreeSelection(fCutBit,fSystematic)) continue;
      fChi2->Fill(fTreeVariablesArray[8]);
      fcRows->Fill(fTreeVariablesArray[7]);
    }
    //
    //
    if( fCentBin != fCentInputBin ) continue;
    if( fMomBin < momRange[0] || fMomBin > momRange[1] ) continue;
    if( fEtaBin < etaRange[0] || fEtaBin > etaRange[1] ) continue;
    //
    if (fSignInput==0) { fUsedBins[fEtaBin][fCentBin][fMomBin][fSignBin] = 1; iden4 -> AddEntry(); countEntry++;}
    else { fUsedBins[fEtaBin][fCentBin][fMomBin][fSignBin] = 1; iden4 -> AddEntry(); countEntry++;}
    //
  }
  cout << "main.Info: Total number of tracks processed = " << countEntry << endl;
  iden4 -> Finalize();
  //
  // Calculate 2. order moments only for full range
  cout << " ==================================" << endl;
  cout << " main.Info: calculating integrals " <<endl;
  timer.Reset(); timer.Start();
  cout << " ==================================" << endl;
  for(Int_t i = 0; i < fnEtaBins; i++){
    for(Int_t j = 0; j < fnCentBins; j++){
      for(Int_t k = 0; k < fnMomBins; k++){
        for(Int_t l = 0; l < fnSignBins; l++) {

          // if (l!=(fSignInput+1))
          if(fUsedBins[i][j][k][l] != 1) continue;
          fEtaBin   =  i;
          fCentBin  =  j;
          fMomBin   =  k;
          fSignBin  =  l;
          iden4  -> AddIntegrals(fSignInput); // real sign information passed for the check with real data tree
        }
      }
    }
  }
  iden4 -> CalcMoments();
  //
  // Retrive Moments
  if (fSignInput== 0) RetrieveMoments(iden4,fMmoments,fIntegrals);
  if (fSignInput== 1) RetrieveMoments(iden4,fMmomentsPos,fIntegrals);
  if (fSignInput==-1) RetrieveMoments(iden4,fMmomentsNeg,fIntegrals);
  cout << "====================================" << endl;
  cout << " main.Info: calculation is finished " << endl;
  timer.Stop(); timer.Print();
  cout << "====================================" << endl;
  //
  // Fill output tree
  cout << " main.Info: fill tree " << endl;
  momTree -> Fill();
  //
  // Close file and clear memory
  cout << " main.Info: dump output " << endl;
  outFile -> cd();
  fChi2   -> Write("fChi2");
  fcRows  -> Write("fcRows");
  momTree -> Write();
  outFile -> Close();
  delete outFile; //yeni eklave etdim.
  delete iden4;
  return 1;
}
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
Bool_t ApplyTreeSelection(UInt_t cut, Int_t syst)
{
  /*
  syst: 0 -->  Reference
  1 -->  CRows60
  2 -->  CRows100
  3 -->  Chi2TPC3
  4 -->  Chi2TPC5
  5 -->  DCAXYSmall
  6 -->  DCAXYLarge
  7 -->  VZSmall
  8 -->  VZLarge
  9 -->  EventVertexZSmall
  10 --> EventVertexZLarge
  11 --> ClusterRequirementITS
  12 --> NewITSCut
  */

  const Int_t fnCutBins=7;
  Int_t fCutArr[fnCutBins]={0};
  if(syst==0) {  // 0 -->  Reference
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==1) {  // 1 -->  CRows60
    Int_t fCutArrTmp[fnCutBins] = {kCRows60,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==2) {  // 2 -->  CRows100
    Int_t fCutArrTmp[fnCutBins] = {kCRows100, kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==3) {  // 3 -->  Chi2TPC3
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC3, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==4) {  // 4 -->  Chi2TPC5
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC5, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==5) {  // 5 -->  kDCAXYSmall
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXYSmall, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==6) {  // 6 -->  kDCAXYLarge
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXYLarge, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==7) {  // 7 -->  kVZSmall
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZSmall, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==8) {  // 8 -->  kVZLarge
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZLarge, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==9) {  // 9 -->  kEventVertexZSmall
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZSmall, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==10) { // 10 -->  kEventVertexZLarge
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZLarge, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==11) { // 11 -->  kNewITSCut
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kNewITSCut, kCRows80};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==12) { // 12 -->  kClusterRequirementITS
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kCRows80};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  //
  //  Apply conditions
  for (Int_t i=0;i<fnCutBins;i++){
    if( ((cut >> fCutArr[i]) & 1) == 0 ) return kFALSE;
  }

  return kTRUE;

}
// -----------------------------------------------------------------------------------------
void ReadFitParamsFromLineShapes(TString paramTreeName)
{
  //
  // Read the fit parameters from the paramtree file "ParamTree.root" which comes from PID FITS
  //
  // Open the lookup table and initialize array
  fLineShapesLookUpTable = new TFile(paramTreeName);
  cloneArrFunc   = (TClonesArray*)fLineShapesLookUpTable->Get("funcLineShapes");
  cloneArrHist   = (TClonesArray*)fLineShapesLookUpTable->Get("histLineShapes");
  if (!cloneArrFunc) cout << " Error:: cloneArrFunc is empty " << endl;
  if (!cloneArrHist) cout << " Error:: cloneArrFunc is empty " << endl;
  Int_t centBinRange[2] = {TMath::Max(fCentInputBin-1,0), TMath::Min(fCentInputBin+1,fnCentBins)};
  Int_t etaBinRange[2]  = {TMath::Max(fEtaDownBin-1,0), TMath::Min(fEtaUpBin+1,fnEtaBins)};
  Int_t momBinRange[2]  = {TMath::Max(fpDownBin-1,0), TMath::Min(fpUpBin+1,fnMomBins)};
  cout << "==================================" << endl;
  cout << "   Reading the file is started " << endl;
  cout << "==================================" << endl;
  cout << "cent Bin window = " << centBinRange[0] << "  " << centBinRange[1] << endl;
  cout << "eta  Bin window = " << etaBinRange[0]  << "  " << etaBinRange[1] << endl;
  cout << "mom  Bin window = " << momBinRange[0]  << "  " << momBinRange[1] << endl;
  cout << "==================================" << endl;
  //
  //
  timer.Start();
  Int_t counter = 0;
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    for (Int_t icent = centBinRange[0]; icent< centBinRange[1]; icent++){
      for (Int_t ieta = etaBinRange[0]; ieta<etaBinRange[1]; ieta++){
        for (Int_t imom = momBinRange[0]; imom<momBinRange[1]; imom++){
          for (Int_t isign = 0; isign<fnSignBins; isign++){
            //
            //
            TString objName = Form("particle_%d_bin_%d_bin_%d_bin_%d_bin_%d",ipart,icent,ieta,imom,isign);
            TString tmpName = Form("tmp_%d_bin_%d_bin_%d_bin_%d_bin_%d",ipart,icent,ieta,imom,isign);
            //
            // recreate the function while reading !!! there is a bug in TClonesArray
            TF1 *tmp = (TF1*)cloneArrFunc->FindObject(objName); tmp->SetName(tmpName);
            Double_t *par = tmp->GetParameters();
            fLineShape[ipart][icent][ieta][imom][isign] = new TF1(objName,fitFunctionGenGausStr,20.,1020.);
            fLineShape[ipart][icent][ieta][imom][isign]->SetParameters(par[0],par[1],par[2],par[3],par[4]);
            fLineShape[ipart][icent][ieta][imom][isign]->SetNpx(nBinsLineShape);
            //
            //
            TH1D *htemp = NULL;
            htemp = (TH1D*)fLineShape[ipart][icent][ieta][imom][isign]->GetHistogram();
            if (!htemp) {
              cout << "testIdentity::Error: " << objName << " failed during fitting " << endl;
              htemp = (TH1D*)fLineShape[ipart][icent][ieta][(Int_t)TMath::Max(0,imom-1)][isign]->GetHistogram();
            }
            if (!htemp) htemp = (TH1D*)fLineShape[ipart][icent][ieta][(Int_t)TMath::Min(fnMomBins,imom+1)][isign]->GetHistogram();
            //
            //
            htemp->SetName(objName);
            hLineShape[ipart][icent][ieta][imom][isign] = htemp;
            if (counter%2000==0 && counter>0){
              cout << counter << " par  = " << ipart << " ------ input function = " << tmp->GetMaximum();
              cout << " new function = "    << fLineShape[ipart][icent][ieta][imom][isign]->GetMaximum();
              cout << " final histogram = " << hLineShape[ipart][icent][ieta][imom][isign]->GetMaximum() << endl;
            }
            counter++;
          }
        }
      }
    }
  }
  cout << "==================================" << endl;
  cout << "   Read both histogram and function lookuptable  " << endl;
  timer.Stop(); timer.Print();
  cout << "==================================" << endl;
}
// -----------------------------------------------------------------------------------------
void ReadFitParamsFromTree(TString paramTreeName, Int_t fitIter)
{
  //
  // Read the fit parameters from the paramtree file "ParamTree.root" which comes from PID FITS
  //
  // Open the lookup table and initialize array
  fLineShapesLookUpTable = new TFile(paramTreeName);
  cout << " ReadFitParamsFromTree.Info: Read line shapes from ttree " << endl;
  treeLookUp = (TTree*)fLineShapesLookUpTable->Get(treeLineShapes);
  Int_t nent = treeLookUp -> GetEntries();
  cout << nent <<  " ReadFitParamsFromTree.Info: tree is fine go ahead " << endl;

  Int_t sign       = 0;
  Int_t it         = 0;
  Int_t sl         = 0;
  Int_t signBin = -100;

  Int_t myBin[3]  = {0};  // 0; eta, 1;cent, 2;ptot

  Double_t elA  = 0., elM  = 0., elSi = 0., elK = 0., elSk = 0.;
  Double_t piA  = 0., piM  = 0., piSi = 0., piK = 0., piSk = 0.;
  Double_t kaA  = 0., kaM  = 0., kaSi = 0., kaK = 0., kaSk = 0.;
  Double_t prA  = 0., prM  = 0., prSi = 0., prK = 0., prSk = 0.;

  treeLookUp->SetBranchAddress("sign"   ,&sign);
  treeLookUp->SetBranchAddress("it"     ,&it);
  treeLookUp->SetBranchAddress("sl"     ,&sl);
  treeLookUp->SetBranchAddress("myBin"  ,myBin);

  treeLookUp->SetBranchAddress("elM" ,&elM);
  treeLookUp->SetBranchAddress("piM" ,&piM);
  treeLookUp->SetBranchAddress("kaM" ,&kaM);
  treeLookUp->SetBranchAddress("prM" ,&prM);

  treeLookUp->SetBranchAddress("elSi" ,&elSi);
  treeLookUp->SetBranchAddress("piSi" ,&piSi);
  treeLookUp->SetBranchAddress("kaSi" ,&kaSi);
  treeLookUp->SetBranchAddress("prSi" ,&prSi);

  treeLookUp->SetBranchAddress("elA" ,&elA);
  treeLookUp->SetBranchAddress("piA" ,&piA);
  treeLookUp->SetBranchAddress("kaA" ,&kaA);
  treeLookUp->SetBranchAddress("prA" ,&prA);

  treeLookUp->SetBranchAddress("elSk" ,&elSk);
  treeLookUp->SetBranchAddress("piSk" ,&piSk);
  treeLookUp->SetBranchAddress("kaSk" ,&kaSk);
  treeLookUp->SetBranchAddress("prSk" ,&prSk);

  treeLookUp->SetBranchAddress("elK" ,&elK);
  treeLookUp->SetBranchAddress("piK" ,&piK);
  treeLookUp->SetBranchAddress("kaK" ,&kaK);
  treeLookUp->SetBranchAddress("prK" ,&prK);

  for(Int_t i = 0; i < nent; ++i)
  {
    // myBin[0] --> Eta, myBin[1]--> Centrality, myBin[2]-->Momentum
    treeLookUp -> GetEntry(i);
    // signBin = sign+1;
    if (sign==-1) signBin=0;
    if (sign== 1) signBin=1;

    // Analyse only 4 iteration of the iterative fitting procedure
    if (it != fitIter) continue;
    fAmpArr[myBin[0]][myBin[1]][myBin[2]][kEl][signBin] = elA;
    fAmpArr[myBin[0]][myBin[1]][myBin[2]][kPi][signBin] = piA;
    fAmpArr[myBin[0]][myBin[1]][myBin[2]][kKa][signBin] = kaA;
    fAmpArr[myBin[0]][myBin[1]][myBin[2]][kPr][signBin] = prA;

    fMeanArr[myBin[0]][myBin[1]][myBin[2]][kEl][signBin] = elM;
    fMeanArr[myBin[0]][myBin[1]][myBin[2]][kPi][signBin] = piM;
    fMeanArr[myBin[0]][myBin[1]][myBin[2]][kKa][signBin] = kaM;
    fMeanArr[myBin[0]][myBin[1]][myBin[2]][kPr][signBin] = prM;

    fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kEl][signBin] = elSi;
    fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kPi][signBin] = piSi;
    fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kKa][signBin] = kaSi;
    fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kPr][signBin] = prSi;

    fSkewArr[myBin[0]][myBin[1]][myBin[2]][kEl][signBin] = elSk;
    fSkewArr[myBin[0]][myBin[1]][myBin[2]][kPi][signBin] = piSk;
    fSkewArr[myBin[0]][myBin[1]][myBin[2]][kKa][signBin] = kaSk;
    fSkewArr[myBin[0]][myBin[1]][myBin[2]][kPr][signBin] = prSk;

    fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kEl][signBin] = elK;
    fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kPi][signBin] = piK;
    fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kKa][signBin] = kaK;
    fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kPr][signBin] = prK;

  }


}
// -----------------------------------------------------------------------------------------
Double_t EvalFitValue(Int_t particle, Double_t x)
{

  if (lookUpTableForLine){
    Int_t bin = hLineShape[particle][fCentBin][fEtaBin][fMomBin][fSignBin]->FindBin(x);
    if (lookUpTableLineMode==0){
      return hLineShape[particle][fCentBin][fEtaBin][fMomBin][fSignBin]->GetBinContent(bin);
    }
    if (lookUpTableLineMode==1){
      return fLineShape[particle][fCentBin][fEtaBin][fMomBin][fSignBin]->Eval(x);
    }

  } else {
    fgengaus -> SetParameter(0,fAmpArr[fEtaBin][fCentBin][fMomBin][particle][fSignBin]);
    fgengaus -> SetParameter(1,fMeanArr[fEtaBin][fCentBin][fMomBin][particle][fSignBin]);
    fgengaus -> SetParameter(2,fSigmaArr[fEtaBin][fCentBin][fMomBin][particle][fSignBin]);
    fgengaus -> SetParameter(3,fKurtosisArr[fEtaBin][fCentBin][fMomBin][particle][fSignBin]);
    fgengaus -> SetParameter(4,fSkewArr[fEtaBin][fCentBin][fMomBin][particle][fSignBin]);
    return fgengaus->Eval(x);
  }

}
// -----------------------------------------------------------------------------------------
void InitializeObjects()
{
  //
  fChi2   = new TH1D("fChi2" ,"fChi2",200 ,0, 10 );
  fcRows  = new TH1D("fcRows","fcRows",200 ,0, 200 );
  fhEta   = new TH1D("fhEta" ,"Eta Bins"       ,fnEtaBins ,fEtaRangeDown, fEtaRangeUp );
  fhPtot  = new TH1D("fhPtot","Momentum Bins"  ,fnMomBins ,fMomRangeDown, fMomRangeUp );
  fhCent  = new TH1D("fhCent","Centrality Bins",fnCentBins ,xCentBins );
  //
  // initialize lookup arrays
  hLineShape = new TH1D *****[fnParticleBins];
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    hLineShape[ipart] = new TH1D****[fnCentBins];
    for (Int_t icent = 0; icent<fnCentBins; icent++){
      hLineShape[ipart][icent] = new TH1D***[fnEtaBins];
      for (Int_t ieta = 0; ieta<fnEtaBins; ieta++){
        hLineShape[ipart][icent][ieta] = new TH1D**[fnMomBins];
        for (Int_t imom = 0; imom<fnMomBins; imom++){
          hLineShape[ipart][icent][ieta][imom] = new TH1D*[fnSignBins];
          for (Int_t isign = 0; isign<fnSignBins; isign++){
            hLineShape[ipart][icent][ieta][imom][isign] = NULL;
          }
        }
      }
    }
  }
  //
  //
  fLineShape = new TF1 *****[fnParticleBins];
  for (Int_t ipart = 0; ipart<fnParticleBins; ipart++){
    fLineShape[ipart] = new TF1****[fnCentBins];
    for (Int_t icent = 0; icent<fnCentBins; icent++){
      fLineShape[ipart][icent] = new TF1***[fnEtaBins];
      for (Int_t ieta = 0; ieta<fnEtaBins; ieta++){
        fLineShape[ipart][icent][ieta] = new TF1**[fnMomBins];
        for (Int_t imom = 0; imom<fnMomBins; imom++){
          fLineShape[ipart][icent][ieta][imom] = new TF1*[fnSignBins];
          for (Int_t isign = 0; isign<fnSignBins; isign++){
            fLineShape[ipart][icent][ieta][imom][isign] = NULL;
          }
        }
      }
    }
  }
  //
  //
  for( Int_t i = 0; i < fnEtaBins; i++ ){
    for( Int_t j = 0; j < fnCentBins; j++ ){
      for( Int_t k = 0; k < fnMomBins; k++ ){
        for(Int_t l = 0; l < fnSignBins; l++)
        {
          fUsedBins[i][j][k][l] = -1;
        }
      }
    }
  }
  //
  // initialize counters
  fIntegrals    = new TVectorF(fnParticleBins);
  fIntegralsPos = new TVectorF(fnParticleBins);
  fIntegralsNeg = new TVectorF(fnParticleBins);
  fMmoments     = new TVectorF(nMoments);
  fMmomentsPos  = new TVectorF(nMoments);
  fMmomentsNeg  = new TVectorF(nMoments);
  for(Int_t i=0;i<nMoments; i++){
    if (i<fnParticleBins){
      (*fIntegrals)[i]=0.;
      (*fIntegralsPos)[i]=0.;
      (*fIntegralsNeg)[i]=0.;
    }
    (*fMmoments)[i]=0.;
    (*fMmomentsPos)[i]=0.;
    (*fMmomentsNeg)[i]=0.;
  }
  //
  // initialize output tree
  //
  momTree = new TTree("momTree","momTree");
  momTree -> Branch("nEvents",&nEvents);
  momTree -> Branch("sign",&fSignInput);
  momTree -> Branch("nnorm",&nnorm);
  momTree -> Branch("pDown",&fpDownInput);
  momTree -> Branch("pUp",&fpUpInput);
  momTree -> Branch("etaDown",&fEtaDownInput);
  momTree -> Branch("etaUp",&fEtaUpInput);
  momTree -> Branch("cent",&fCentInput);
  momTree -> Branch("centBin",&fCentInputBin);
  momTree -> Branch("fIntegrals",&fIntegrals);
  momTree -> Branch("fIntegralsPos",&fIntegralsPos);
  momTree -> Branch("fIntegralsNeg",&fIntegralsNeg);
  momTree -> Branch("fMmoments",&fMmoments);
  momTree -> Branch("fMmomentsPos",&fMmomentsPos);
  momTree -> Branch("fMmomentsNeg",&fMmomentsNeg);

}
// -----------------------------------------------------------------------------------------
TString PrintNumInBinary(UInt_t num)
{
  TString bin="";
  Int_t numberOfBits = sizeof(UInt_t)*8;
  for (Int_t i=numberOfBits-1; i>=0; i--) {
    Bool_t isBitSet = (num & (1<<i));
    if (isBitSet) {
      bin+="1";
    } else {
      bin+="0";
    }
  }
  return bin;
}
// -----------------------------------------------------------------------------------------
Double_t fitFunctionGenGaus(Double_t *x, Double_t *par)
{
  //
  // Generalised gauss function --> "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
  // Skew-normal distribution   --> "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]/TMath::Sqrt(2)))";
  // par[0] --> Amplitude
  // par[1] --> Mean
  // par[2] --> Sigma
  // par[3] --> Kurtosis
  // par[4] --> Skewness
  //

  Double_t fun = par[0]*exp(-TMath::Power((TMath::Abs(x[0]-par[1])/par[2]),par[3]))
  *(1+TMath::Erf(par[4]*(x[0]-par[1])/par[2]/TMath::Sqrt(2)));
  return fun;
}
// -----------------------------------------------------------------------------------------
void RetrieveMoments(TIdentity2D *tidenObj, TVectorF *vecMom, TVectorF *vecInt)
{

  (*vecMom)[kEl] = tidenObj -> GetMean(kEl);
  (*vecMom)[kPi] = tidenObj -> GetMean(kPi);
  (*vecMom)[kKa] = tidenObj -> GetMean(kKa);
  (*vecMom)[kPr] = tidenObj -> GetMean(kPr);
  //
  // Second Moments
  (*vecMom)[kElEl] = tidenObj -> GetSecondMoment(kEl);
  (*vecMom)[kPiPi] = tidenObj -> GetSecondMoment(kPi);
  (*vecMom)[kKaKa] = tidenObj -> GetSecondMoment(kKa);
  (*vecMom)[kPrPr] = tidenObj -> GetSecondMoment(kPr);
  //
  // Mixed Moments
  (*vecMom)[kElPi] = tidenObj -> GetMixedMoment(kEl,kPi);
  (*vecMom)[kElKa] = tidenObj -> GetMixedMoment(kEl,kKa);
  (*vecMom)[kElPr] = tidenObj -> GetMixedMoment(kEl,kPr);
  (*vecMom)[kPiKa] = tidenObj -> GetMixedMoment(kPi,kKa);
  (*vecMom)[kPiPr] = tidenObj -> GetMixedMoment(kPi,kPr);
  (*vecMom)[kKaPr] = tidenObj -> GetMixedMoment(kKa,kPr);
  //
  //Integrals:
  (*vecInt)[kEl] = tidenObj -> GetMeanI(kEl);
  (*vecInt)[kPi] = tidenObj -> GetMeanI(kPi);
  (*vecInt)[kKa] = tidenObj -> GetMeanI(kKa);
  (*vecInt)[kPr] = tidenObj -> GetMeanI(kPr);
  //
  // Printing
  nnorm     = (*vecMom)[kPi]/(*vecInt)[kPi];
  nEvents   = tidenObj -> GetNEvents();
  cout << " =============================== Summary of Moments =============================== "<<endl;
  cout << " events      : "<< nEvents << endl;
  cout << " sign        : "<< fSignInput << endl;
  cout << " ================================================================================== "<<endl;
  cout << " electron    : "<< (*vecMom)[kEl]   <<" int: "<< (*vecInt)[kEl]*nnorm << "  ratio: " << (*vecMom)[kEl]/((*vecInt)[kEl]*nnorm) << endl;
  cout << " pion        : "<< (*vecMom)[kPi]   <<" int: "<< (*vecInt)[kPi]*nnorm << "  ratio: " << (*vecMom)[kPi]/((*vecInt)[kPi]*nnorm) << endl;
  cout << " kaon        : "<< (*vecMom)[kKa]   <<" int: "<< (*vecInt)[kKa]*nnorm << "  ratio: " << (*vecMom)[kKa]/((*vecInt)[kKa]*nnorm) << endl;
  cout << " proton      : "<< (*vecMom)[kPr]   <<" int: "<< (*vecInt)[kPr]*nnorm << "  ratio: " << (*vecMom)[kPr]/((*vecInt)[kPr]*nnorm) << endl;
  cout << " electron2   : "<< (*vecMom)[kElEl]  <<endl;
  cout << " pion2       : "<< (*vecMom)[kPiPi]  <<endl;
  cout << " kaon2       : "<< (*vecMom)[kKaKa]  <<endl;
  cout << " proton2     : "<< (*vecMom)[kPrPr]  <<endl;
  //
  //
  for (Int_t i=1;i<nMoments+1;i++) cout << momNames[i-1] << "  " << (*vecMom)[i-1] << endl;



}
// -----------------------------------------------------------------------------------------
void PrintInitInfo()
{

  //
  cout << " ================================================================================= " << endl;
  cout << " InitializeObjects.Info: treeLineShapes        = " << treeLineShapes       << endl;
  cout << " InitializeObjects.Info: treeIdentity          = " << treeIdentity       << endl;
  cout << " ================================================================================= " << endl;
  cout << " InitializeObjects.Info: Inputs: " << endl;
  cout << " InitializeObjects.Info: data Tree             = " << fileNameDataTree       << endl;
  cout << " InitializeObjects.Info: Line Shapes           = " << fileNameLineShapes     << endl;
  cout << " InitializeObjects.Info: Centrality            = " << fCentInput      << endl;
  cout << " InitializeObjects.Info: Subsample index       = " << fSubsample      << endl;
  cout << " InitializeObjects.Info: Eta Range             = " << fEtaDownInput       << " - " <<  fEtaUpInput << endl;
  cout << " InitializeObjects.Info: Momentum Range        = " << fpDownInput         << " - " <<  fpUpInput << endl;
  cout << " InitializeObjects.Info: Fit iteration         = " << fNthFitIteration    << endl;
  cout << " InitializeObjects.Info: Charge                = " << fSignInput         << endl;
  cout << " ================================================================================= " << endl;
  cout << " InitializeObjects.Info: Centrality (Bin)      = " << fCentInputBin      << endl;
  cout << " InitializeObjects.Info: Eta Range (Bin)       = " << fEtaDownBin    << " - " <<  fEtaUpBin << endl;
  cout << " InitializeObjects.Info: Momentum Range (Bin)  = " << fpDownBin    << " - " <<  fpUpBin << endl;
  cout << " ================================================================================= " << endl;
  //
}


/*

TFile f("/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/RUN1/PHD/Systematics_cRows_80_16EtaBin_mombin20MeV/ParamTrees/LineShapes_ClonesArray.root")
// TFile f("LineShapes_ClonesArray.root")
cloneArrFunc   = (TObjArray*)f->Get("funcLineShapes");
cloneArrHist   = (TClonesArray*)f->Get("histLineShapes");

hPi  = (TH1D*)cloneArrHist->FindObject("particle_1_bin_0_bin_0_bin_0_bin_0");
fPi  = (TF1*)cloneArrFunc->FindObject("particle_1_bin_0_bin_0_bin_0_bin_0");
hPi->Draw();
fPi->Draw("same");

TF1 *uu = new TF1("uu","[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))",20.,1020.);
Double_t *par = fPi->GetParameters();
uu->SetParameters(par[0],par[1],par[2],par[3],par[4]);
uu->SetNpx(4000);
uu->SetLineColor(kBlack);
uu->Draw("same");

TH1D *k = (TH1D*)uu->GetHistogram()->Clone();
k->SetLineColor(kMagenta);
k->Draw("same");


*/
