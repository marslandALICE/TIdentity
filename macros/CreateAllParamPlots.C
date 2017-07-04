#include <TFile.h>
#include <TH3.h>
#include "THn.h"
#include "TCut.h"
#include "TCutG.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <THnSparse.h>
#include "TMath.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TTreeStream.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLinearFitter.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "AliXRDPROOFtoolkit.h"
#include "TStatToolkit.h"
#include "AliMathBase.h"
#include <fstream>
#include <iostream>

using namespace std;

// 
// Make TH2D from the TTree chain which will be used for the pt slices
// 

void InitInitials();
void MergeFitResultsTree(TString fitTreeFile, TString paramFile);
void CheckResults(TString paramFile);

    
// ======= Modification part =======  
Bool_t fpp     = kFALSE;
const Int_t ncentbins = 10; Float_t xCentBins[ncentbins] = {0, 5,  10,  20, 30, 40, 50, 60, 70, 80}; 
// const Int_t ncentbins = 19; Float_t xCentBins[ncentbins] = {0, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80}; 
Int_t firstIter = 5;
Int_t nIter     = 8;
Int_t nCentbins = ncentbins-1;
 
Double_t ptMin      = 0.2;
Double_t ptMax      = 3.2;
Double_t ptNbins    = 150;   // for 20MeV slice 150, for 10MeV slice 300
Double_t etaMin     = -0.8;
Double_t etaMax     = 0.8;
Int_t nEtabins      = 32;

// ======= Modification part =======  

TH1D  *hEta, *hCent, *hMom;


void CreateAllParamPlots(Int_t actionType, TString pidTreeFile, TString paramFile){
  //
  // Produce all hists form MC and Real data
  /*
     
    cd /hera/alice/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_mombin20MeV/Fits/FitResults_pi0.8_ka0.5_pr0.4
    aliroot -l -b
    .L ~/PHD/macros/marsland_EbyeRatios/CreateAllParamPlots.C+
    CreateAllParamPlots(0,"AllPIDtrees_sign0.root","ParamTree_sign0.root")     
    CreateAllParamPlots(1,"","ParamTree_sign0.root")    
    
  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllParamPlots.C .");
  
  if (actionType==0){
    MergeFitResultsTree(pidTreeFile,paramFile);
  } else if (actionType==1){
    CheckResults(paramFile);  
  } 
   
}
//____________________________________________________________________________________________________________
void MergeFitResultsTree(TString fitTreeFile, TString paramFile)
{
  
  //
  // Obtain ttree with sequentially ordered events
  // fitTreeFile   : merged tree of fits form PIDIterativeFitting.C
  /*
  
  cd /hera/alice/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_mombin20MeV/ParamTrees
  aliroot -l 
  .L ~/PHD/macros/marsland_EbyeRatios/CreateAllParamPlots.C+ 
  MergeFitResultsTree("AllPIDtrees_sign0.root","ParamTree_sign0.root")
   
  */
  //
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllParamPlots.C .");

  // Initialise Histograms for the myBin[] varibale
  InitInitials();
  
  TTree *tree;
  if (fitTreeFile.Contains(".root")){
    TFile *f = TFile::Open(fitTreeFile);
    tree = (TTree*)f->Get("IdMethodInput");  
  } else {
    TChain *chain = AliXRDPROOFtoolkit::MakeChainRandom(fitTreeFile,"IdMethodInput",0,500000,0);
    chain -> SetCacheSize(10000000000000);
    tree = (TTree*)chain;
  }
 
  Double_t nEntries = tree->GetEntries();
  
  Int_t myBin[3]      = {0.,0.,0.};
  Double_t totChi2    = 0.;
  Double_t ndf        = 0.; 
  Int_t sign          = 0;
  Int_t sl            = 0;   
  Int_t it            = 0;
  Double_t piSScan   = 0.;
  Double_t kaSScan   = 0.;
  Double_t prSScan   = 0.;
  Double_t piKScan   = 0.;
  Double_t kaKScan   = 0.;
  Double_t prKScan   = 0.;

  Int_t fSign         = 0;
  Int_t iter          = 0;
  Int_t slice         = 0;    
  Double_t p          = 0.;
  Float_t eta         = 0.;
  Float_t cent        = 0;
  Double_t totalChi2  = 0.;
  Double_t normNDF    = 0.; 
  Double_t piSkewScan = 0.;
  Double_t kaSkewScan = 0.;
  Double_t prSkewScan = 0.;
  Double_t piKurtScan = 0.;
  Double_t kaKurtScan = 0.;
  Double_t prKurtScan = 0.;

  
  Double_t elAmp      = 0., elA = 0.;
  Double_t piAmp      = 0., piA = 0.;
  Double_t kaAmp      = 0., kaA = 0.;
  Double_t prAmp      = 0., prA = 0.;
        
  Double_t elMean     = 0., elM = 0.;
  Double_t piMean     = 0., piM = 0.;
  Double_t kaMean     = 0., kaM = 0.;
  Double_t prMean     = 0., prM = 0.;
        
  Double_t elSigma    = 0., elSi = 0.;
  Double_t piSigma    = 0., piSi = 0.;
  Double_t kaSigma    = 0., kaSi = 0.;
  Double_t prSigma    = 0., prSi = 0.;
  
  Double_t elSkew     = 0., elSk = 0.;
  Double_t piSkew     = 0., piSk = 0.;
  Double_t kaSkew     = 0., kaSk = 0.;
  Double_t prSkew     = 0., prSk = 0.;
               
  Double_t elKurtosis = 0., elK = 0.;
  Double_t piKurtosis = 0., piK = 0.;
  Double_t kaKurtosis = 0., kaK = 0.;
  Double_t prKurtosis = 0., prK = 0.;
  
  Double_t elAmpErr      = 0., elAErr = 0.;
  Double_t piAmpErr      = 0., piAErr = 0.;
  Double_t kaAmpErr      = 0., kaAErr = 0.;
  Double_t prAmpErr      = 0., prAErr = 0.;
        
  Double_t elMeanErr     = 0., elMErr = 0.;
  Double_t piMeanErr     = 0., piMErr = 0.;
  Double_t kaMeanErr     = 0., kaMErr = 0.;
  Double_t prMeanErr     = 0., prMErr = 0.;
        
  Double_t elSigmaErr    = 0., elSiErr = 0.;
  Double_t piSigmaErr    = 0., piSiErr = 0.;
  Double_t kaSigmaErr    = 0., kaSiErr = 0.;
  Double_t prSigmaErr    = 0., prSiErr = 0.;
  
  Double_t elSkewErr     = 0., elSkErr = 0.;
  Double_t piSkewErr     = 0., piSkErr = 0.;
  Double_t kaSkewErr     = 0., kaSkErr = 0.;
  Double_t prSkewErr     = 0., prSkErr = 0.;
               
  Double_t elKurtosisErr = 0., elKErr = 0.;
  Double_t piKurtosisErr = 0., piKErr = 0.;
  Double_t kaKurtosisErr = 0., kaKErr = 0.;
  Double_t prKurtosisErr = 0., prKErr = 0.;
  
  
  Double_t elInt = 0., elI = 0.;
  Double_t piInt = 0., piI = 0.;
  Double_t kaInt = 0., kaI = 0.;
  Double_t prInt = 0., prI = 0.;
       
  // kinematic variables
  tree->SetBranchAddress("fSign"     ,&fSign);
  tree->SetBranchAddress("iter"      ,&iter);
  tree->SetBranchAddress("slice"     ,&slice);
  tree->SetBranchAddress("p"         ,&p);
  tree->SetBranchAddress("eta"       ,&eta);
  tree->SetBranchAddress("cent"      ,&cent);
  tree->SetBranchAddress("totalChi2" ,&totalChi2);
  tree->SetBranchAddress("normNDF"   ,&normNDF);
  
  // kurtosis and skewness scan paramemters
  tree->SetBranchAddress("piSkewScan"  ,&piSkewScan);
  tree->SetBranchAddress("kaSkewScan"  ,&kaSkewScan);
  tree->SetBranchAddress("prSkewScan"  ,&prSkewScan);
  tree->SetBranchAddress("piKurtScan"  ,&piKurtScan);
  tree->SetBranchAddress("kaKurtScan"  ,&kaKurtScan);
  tree->SetBranchAddress("prKurtScan"  ,&prKurtScan);

  // Fit Params
  tree->SetBranchAddress("elMean" ,&elMean);
  tree->SetBranchAddress("piMean" ,&piMean);
  tree->SetBranchAddress("kaMean" ,&kaMean);
  tree->SetBranchAddress("prMean" ,&prMean);
  
  tree->SetBranchAddress("elSigma" ,&elSigma);
  tree->SetBranchAddress("piSigma" ,&piSigma);
  tree->SetBranchAddress("kaSigma" ,&kaSigma);
  tree->SetBranchAddress("prSigma" ,&prSigma);
  
  tree->SetBranchAddress("elAmp" ,&elAmp);
  tree->SetBranchAddress("piAmp" ,&piAmp);
  tree->SetBranchAddress("kaAmp" ,&kaAmp);
  tree->SetBranchAddress("prAmp" ,&prAmp);
  
  tree->SetBranchAddress("elSkew" ,&elSkew);
  tree->SetBranchAddress("piSkew" ,&piSkew);
  tree->SetBranchAddress("kaSkew" ,&kaSkew);
  tree->SetBranchAddress("prSkew" ,&prSkew);
  
  tree->SetBranchAddress("elKurtosis" ,&elKurtosis);
  tree->SetBranchAddress("piKurtosis" ,&piKurtosis);
  tree->SetBranchAddress("kaKurtosis" ,&kaKurtosis);
  tree->SetBranchAddress("prKurtosis" ,&prKurtosis);
  
  tree->SetBranchAddress("elMeanErr" ,&elMeanErr);
  tree->SetBranchAddress("piMeanErr" ,&piMeanErr);
  tree->SetBranchAddress("kaMeanErr" ,&kaMeanErr);
  tree->SetBranchAddress("prMeanErr" ,&prMeanErr);
  
  tree->SetBranchAddress("elSigmaErr" ,&elSigmaErr);
  tree->SetBranchAddress("piSigmaErr" ,&piSigmaErr);
  tree->SetBranchAddress("kaSigmaErr" ,&kaSigmaErr);
  tree->SetBranchAddress("prSigmaErr" ,&prSigmaErr);
  
  tree->SetBranchAddress("elAmpErr" ,&elAmpErr);
  tree->SetBranchAddress("piAmpErr" ,&piAmpErr);
  tree->SetBranchAddress("kaAmpErr" ,&kaAmpErr);
  tree->SetBranchAddress("prAmpErr" ,&prAmpErr);
  
  tree->SetBranchAddress("elSkewErr" ,&elSkewErr);
  tree->SetBranchAddress("piSkewErr" ,&piSkewErr);
  tree->SetBranchAddress("kaSkewErr" ,&kaSkewErr);
  tree->SetBranchAddress("prSkewErr" ,&prSkewErr);
  
  tree->SetBranchAddress("elKurtosisErr" ,&elKurtosisErr);
  tree->SetBranchAddress("piKurtosisErr" ,&piKurtosisErr);
  tree->SetBranchAddress("kaKurtosisErr" ,&kaKurtosisErr);
  tree->SetBranchAddress("prKurtosisErr" ,&prKurtosisErr);
  
  tree->SetBranchAddress("elInt" ,&elInt);
  tree->SetBranchAddress("piInt" ,&piInt);
  tree->SetBranchAddress("kaInt" ,&kaInt);
  tree->SetBranchAddress("prInt" ,&prInt);

  // New Tree
  TFile f(paramFile,"recreate");
  TTree *treeId = new TTree("treeId","IdentityInput FitParamTree");
  
  // kinematic variables
  treeId->Branch("sign"    ,&sign     ,"sign/I");
  treeId->Branch("it"      ,&it       ,"it/I");
  treeId->Branch("sl"      ,&sl       ,"sl/I");
  treeId->Branch("myBin"   ,myBin     ,"myBin[3]/I");
  treeId->Branch("totChi2" ,&totChi2  ,"totChi2/D");
  treeId->Branch("ndf"     ,&ndf      ,"ndf/D");
  
  // kurtosis and skewness scan paramemters
  treeId->Branch("piSScan"   ,&piSScan    ,"piSScan/D");
  treeId->Branch("kaSScan"   ,&kaSScan    ,"kaSScan/D");
  treeId->Branch("prSScan"   ,&prSScan    ,"prSScan/D");
  treeId->Branch("piKScan"   ,&piKScan    ,"piKScan/D");
  treeId->Branch("kaKScan"   ,&kaKScan    ,"kaKScan/D");
  treeId->Branch("prKScan"   ,&prKScan    ,"prKScan/D");
  
  // Fit Parameters
  treeId->Branch("elM"     ,&elM      ,"elM/D");
  treeId->Branch("piM"     ,&piM      ,"piM/D");
  treeId->Branch("kaM"     ,&kaM      ,"kaM/D");
  treeId->Branch("prM"     ,&prM      ,"prM/D");
  
  treeId->Branch("elSi"    ,&elSi     ,"elSi/D");
  treeId->Branch("piSi"    ,&piSi     ,"piSi/D");
  treeId->Branch("kaSi"    ,&kaSi     ,"kaSi/D");
  treeId->Branch("prSi"    ,&prSi     ,"prSi/D");
  
  treeId->Branch("elA"     ,&elA      ,"elA/D");
  treeId->Branch("piA"     ,&piA      ,"piA/D");
  treeId->Branch("kaA"     ,&kaA      ,"kaA/D");
  treeId->Branch("prA"     ,&prA      ,"prA/D");
  
  treeId->Branch("elSk"    ,&elSk     ,"elSk/D");
  treeId->Branch("piSk"    ,&piSk     ,"piSk/D");
  treeId->Branch("kaSk"    ,&kaSk     ,"kaSk/D");
  treeId->Branch("prSk"    ,&prSk     ,"prSk/D");
  
  treeId->Branch("elK"     ,&elK      ,"elK/D");
  treeId->Branch("piK"     ,&piK      ,"piK/D");
  treeId->Branch("kaK"     ,&kaK      ,"kaK/D");
  treeId->Branch("prK"     ,&prK      ,"prK/D");
  
   // Fit Parameter Errors
  treeId->Branch("elMErr"     ,&elMErr      ,"elMErr/D");
  treeId->Branch("piMErr"     ,&piMErr      ,"piMErr/D");
  treeId->Branch("kaMErr"     ,&kaMErr      ,"kaMErr/D");
  treeId->Branch("prMErr"     ,&prMErr      ,"prMErr/D");
  
  treeId->Branch("elSiErr"    ,&elSiErr     ,"elSiErr/D");
  treeId->Branch("piSiErr"    ,&piSiErr     ,"piSiErr/D");
  treeId->Branch("kaSiErr"    ,&kaSiErr     ,"kaSiErr/D");
  treeId->Branch("prSiErr"    ,&prSiErr     ,"prSiErr/D");
  
  treeId->Branch("elAErr"     ,&elAErr      ,"elAErr/D");
  treeId->Branch("piAErr"     ,&piAErr      ,"piAErr/D");
  treeId->Branch("kaAErr"     ,&kaAErr      ,"kaAErr/D");
  treeId->Branch("prAErr"     ,&prAErr      ,"prAErr/D");
  
  treeId->Branch("elSkErr"    ,&elSkErr     ,"elSkErr/D");
  treeId->Branch("piSkErr"    ,&piSkErr     ,"piSkErr/D");
  treeId->Branch("kaSkErr"    ,&kaSkErr     ,"kaSkErr/D");
  treeId->Branch("prSkErr"    ,&prSkErr     ,"prSkErr/D");
  
  treeId->Branch("elKErr"     ,&elKErr      ,"elKErr/D");
  treeId->Branch("piKErr"     ,&piKErr      ,"piKErr/D");
  treeId->Branch("kaKErr"     ,&kaKErr      ,"kaKErr/D");
  treeId->Branch("prKErr"     ,&prKErr      ,"prKErr/D");
  
  treeId->Branch("elI"     ,&elI      ,"elI/D");
  treeId->Branch("piI"     ,&piI      ,"piI/D");
  treeId->Branch("kaI"     ,&kaI      ,"kaI/D");
  treeId->Branch("prI"     ,&prI      ,"prI/D");
  
  cout << nEntries << "  number of entries will be processed " << endl; 
  for (Int_t itrack=0; itrack<nEntries; ++itrack) {
   
    tree->GetEntry(itrack);
    Int_t tmpIter = iter;
    
    // Parameter boundary conditions --> Scan of the previous and next slices on case there is a fit fail
    Bool_t elParamCheck = ((elMean>20 && elMean<1020) && (elSigma>0 && elSigma<200) && (elAmp>=1 && elAmp<1e+30) && (elSkew>=0 && elSkew<3) && (elKurtosis>1.5 && elKurtosis<4.));
    Bool_t piParamCheck = ((piMean>20 && piMean<1020) && (piSigma>0 && piSigma<200) && (piAmp>=1 && piAmp<1e+30) && (piSkew>=0 && piSkew<3) && (piKurtosis>1.5 && piKurtosis<4.));
    Bool_t kaParamCheck = ((kaMean>20 && kaMean<1020) && (kaSigma>0 && kaSigma<200) && (kaAmp>=1 && kaAmp<1e+30) && (kaSkew>=0 && kaSkew<3) && (kaKurtosis>1.5 && kaKurtosis<4.));
    Bool_t prParamCheck = ((prMean>20 && prMean<1020) && (prSigma>0 && prSigma<200) && (prAmp>=1 && prAmp<1e+30) && (prSkew>=0 && prSkew<3) && (prKurtosis>1.5 && prKurtosis<4.));
   
    elM=elMean;  elSi=elSigma;  elA=elAmp;  elSk=elSkew;  elK=elKurtosis;  elI=elInt;
    piM=piMean;  piSi=piSigma;  piA=piAmp;  piSk=piSkew;  piK=piKurtosis;  piI=piInt;
    kaM=kaMean;  kaSi=kaSigma;  kaA=kaAmp;  kaSk=kaSkew;  kaK=kaKurtosis;  kaI=kaInt;
    prM=prMean;  prSi=prSigma;  prA=prAmp;  prSk=prSkew;  prK=prKurtosis;  prI=prInt;
    
    elMErr=elMeanErr;  elSiErr=elSigmaErr;  elAErr=elAmpErr;  elSkErr=elSkewErr;  elKErr=elKurtosisErr;  
    piMErr=piMeanErr;  piSiErr=piSigmaErr;  piAErr=piAmpErr;  piSkErr=piSkewErr;  piKErr=piKurtosisErr; 
    kaMErr=kaMeanErr;  kaSiErr=kaSigmaErr;  kaAErr=kaAmpErr;  kaSkErr=kaSkewErr;  kaKErr=kaKurtosisErr; 
    prMErr=prMeanErr;  prSiErr=prSigmaErr;  prAErr=prAmpErr;  prSkErr=prSkewErr;  prKErr=prKurtosisErr; 
  
    // avoid from failed fits lower limit
    if ( !elParamCheck ) {elM=0.;elSigma=1.;elA=0.;elSk=0.;elK=2.;elI=1;}
    if ( !piParamCheck ) {piM=0.;piSigma=1.;piA=0.;piSk=0.;piK=2.;piI=1;}
    if ( !kaParamCheck ) {kaM=0.;kaSigma=1.;kaA=0.;kaSk=0.;kaK=2.;kaI=1;}
    if ( !prParamCheck ) {prM=0.;prSigma=1.;prA=0.;prSk=0.;prK=2.;prI=1;}

    // Set eta momentum cent etc. and fill the treeId
    myBin[0] = hEta ->FindBin(eta+0.0001)-1;
    myBin[1] = (fpp) ? -1 : hCent->FindBin(cent+0.0001)-1;
    myBin[2] = hMom->FindBin(p+0.0001)-1;
    totChi2=totalChi2; ndf=normNDF; it=iter; sl=slice; sign=fSign;
    piSScan=piSkewScan; kaSScan=kaSkewScan; prSScan=prSkewScan;
    piKScan=piKurtScan; kaKScan=kaKurtScan; prKScan=prKurtScan;
    treeId->Fill();
    
  }
  
  f.cd();
  treeId->Write();
  hEta->Write();
  hCent->Write();
  hMom->Write();

}
//____________________________________________________________________________________________________________
void CheckResults(TString paramFile){
  
  //
  // Look at all fit params for the smoothness of parameters vs momentum
  //
  
  /*
  
  cd /hera/alice/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_mombin20MeV/ParamTrees
  aliroot -l 
  .L ~/PHD/macros/marsland_EbyeRatios/CreateAllParamPlots.C+ 
  CheckResults("ParamTree.root")
   
  */
  //
  
  TTreeSRedirector *streamer = new TTreeSRedirector("ParameterControl.root","recreate");
  
  TFile *fResControl    = TFile::Open(paramFile);
  TTree *treeResControl = (TTree*)fResControl->Get("treeId");  
  treeResControl->SetCacheSize(1000000000);
  // loop over particle types
  for (Int_t iType=0; iType<3; iType++){
    
    TString parType; 
    if (iType==0) parType = "pi";  
    if (iType==1) parType = "ka";  
    if (iType==2) parType = "pr";
  
  // loop over all iterations 
    for (Int_t iIter=firstIter;iIter<nIter;iIter++){
    
      cout << "particle type = " << parType <<  "      iter = " << iIter << endl;
    
      TGraphErrors *grAmp[nEtabins][nCentbins];
      TGraphErrors *grMean[nEtabins][nCentbins];
      TGraphErrors *grSigma[nEtabins][nCentbins];
      TGraphErrors *grKurt[nEtabins][nCentbins];
      TGraphErrors *grSkew[nEtabins][nCentbins];
      TGraphErrors *grInt[nEtabins][nCentbins];
      TGraphErrors *grChi2[nEtabins][nCentbins];
    
      for (Int_t i=0; i<nEtabins; i++)
      {
        for (Int_t j=0; j<nCentbins; j++){
          grAmp[i][j]  =NULL; 
          grMean[i][j] =NULL; 
          grSigma[i][j]=NULL; 
          grKurt[i][j] =NULL;
          grSkew[i][j] =NULL;  
          grInt[i][j]  =NULL;
          grChi2[i][j] =NULL;
        }
      }
    
    // get the graphs for each eta and cent bin
      for (Int_t iEta=0; iEta<nEtabins; iEta++)
      {
        for (Int_t iCent=0; iCent<nCentbins; iCent++)
        {
          //TString mainCut   = Form("it==%d && myBin[1]==%d && myBin[0]==%d && myBin[2]<=140",iIter,iCent,iEta);
          TString mainCut   = Form("it==%d && myBin[1]==%d && myBin[0]==%d",iIter,iCent,iEta);
          TString drawAmp   = Form("%sA:myBin[2]" ,parType.Data());
          TString drawMean  = Form("%sM:myBin[2]" ,parType.Data());
          TString drawSigma = Form("%sSi:myBin[2]",parType.Data());
          TString drawKurt  = Form("%sK:myBin[2]" ,parType.Data());
          TString drawSkew  = Form("%sSk:myBin[2]",parType.Data());
          TString drawInt   = Form("%sI:myBin[2]" ,parType.Data());
          TString drawChi2  = "totChi2/ndf:myBin[2]";
      
          Int_t markerStyle = (iEta<nEtabins/2) ? 20 : 24;
          Int_t markerColor = (iEta<nEtabins/2) ? iEta+1 : (nEtabins/2)-iEta%(nEtabins/2);
          grSigma[iEta][iCent] = TStatToolkit::MakeGraphErrors(treeResControl,drawSigma.Data(),mainCut.Data(),markerStyle,markerColor,0.5);
          grMean[iEta][iCent]  = TStatToolkit::MakeGraphErrors(treeResControl,drawMean.Data() ,mainCut.Data(),markerStyle,markerColor,0.5);
          grAmp[iEta][iCent]   = TStatToolkit::MakeGraphErrors(treeResControl,drawAmp.Data()  ,mainCut.Data(),markerStyle,markerColor,0.5);
          grKurt[iEta][iCent]  = TStatToolkit::MakeGraphErrors(treeResControl,drawKurt.Data() ,mainCut.Data(),markerStyle,markerColor,0.5);
          grSkew[iEta][iCent]  = TStatToolkit::MakeGraphErrors(treeResControl,drawSkew.Data() ,mainCut.Data(),markerStyle,markerColor,0.5);
          grInt[iEta][iCent]   = TStatToolkit::MakeGraphErrors(treeResControl,drawInt.Data()  ,mainCut.Data(),markerStyle,markerColor,0.5);
          grChi2[iEta][iCent]  = TStatToolkit::MakeGraphErrors(treeResControl,drawChi2.Data() ,mainCut.Data(),markerStyle,markerColor,0.5);

          grKurt[iEta][iCent]->GetYaxis()->SetRangeUser(1.7,2.1);
          grSkew[iEta][iCent]->GetYaxis()->SetRangeUser(0.,1.1);

          grSigma[iEta][iCent]->GetXaxis()->SetTitle("p (GeV/c)");  grSigma[iEta][iCent]->GetYaxis()->SetTitle("sigma");
          grMean[iEta][iCent] ->GetXaxis()->SetTitle("p (GeV/c)");  grMean[iEta][iCent]->GetYaxis()->SetTitle("mean");
          grAmp[iEta][iCent]  ->GetXaxis()->SetTitle("p (GeV/c)");  grAmp[iEta][iCent]->GetYaxis()->SetTitle("amplitude");
          grKurt[iEta][iCent] ->GetXaxis()->SetTitle("p (GeV/c)");  grKurt[iEta][iCent]->GetYaxis()->SetTitle("kurtosis");
          grSkew[iEta][iCent] ->GetXaxis()->SetTitle("p (GeV/c)");  grSkew[iEta][iCent]->GetYaxis()->SetTitle("skewness");
          grInt[iEta][iCent]  ->GetXaxis()->SetTitle("p (GeV/c)");  grInt[iEta][iCent]->GetYaxis()->SetTitle("Integral");
          grChi2[iEta][iCent] ->GetXaxis()->SetTitle("p (GeV/c)");  grChi2[iEta][iCent]->GetYaxis()->SetTitle("#chi^{2}/NDF");
          
        } // end of cent loop
      } // end of eta loop
    
      TObjArray paramsPerCentArr(nCentbins); paramsPerCentArr.SetOwner(kTRUE);
      // plots the canvases
      for (Int_t iCent=0;iCent<nCentbins;iCent++){
      
        TCanvas *paramCan = new TCanvas(Form("%s_params_iter%d_cent%d",parType.Data(),iIter,iCent), "All params", 1200, 800);   
        paramCan->Divide(3,3);
    
      // Sigmas
        paramCan->cd(1);
        grSigma[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grSigma[i][iCent]->Draw("p");
        TLegend *lSigma = new TLegend(0.25, 0.55, 0.85, 0.85);
        lSigma->SetTextFont(62); lSigma->SetTextSize(0.03);  lSigma->SetFillColor(0);  lSigma-> SetNColumns(2);
        lSigma->AddEntry(grSigma[0][iCent], " Eta=[-0.8,-0.7]","LP"); lSigma->AddEntry(grSigma[15][iCent]," Eta=[ 0.7, 0.8]","LP");
        lSigma->AddEntry(grSigma[1][iCent], " Eta=[-0.7,-0.6]","LP"); lSigma->AddEntry(grSigma[14][iCent]," Eta=[ 0.6, 0.7]","LP");
        lSigma->AddEntry(grSigma[2][iCent], " Eta=[-0.6,-0.5]","LP"); lSigma->AddEntry(grSigma[13][iCent]," Eta=[ 0.5, 0.6]","LP");
        lSigma->AddEntry(grSigma[3][iCent], " Eta=[-0.5,-0.4]","LP"); lSigma->AddEntry(grSigma[12][iCent]," Eta=[ 0.4, 0.5]","LP");
        lSigma->AddEntry(grSigma[4][iCent], " Eta=[-0.4,-0.3]","LP"); lSigma->AddEntry(grSigma[11][iCent]," Eta=[ 0.3, 0.4]","LP");
        lSigma->AddEntry(grSigma[5][iCent], " Eta=[-0.3,-0.2]","LP"); lSigma->AddEntry(grSigma[10][iCent]," Eta=[ 0.2, 0.3]","LP");
        lSigma->AddEntry(grSigma[6][iCent], " Eta=[-0.2,-0.1]","LP"); lSigma->AddEntry(grSigma[9][iCent], " Eta=[ 0.1, 0.2]","LP");
        lSigma->AddEntry(grSigma[7][iCent], " Eta=[-0.1, 0  ]","LP"); lSigma->AddEntry(grSigma[8][iCent], " Eta=[ 0  , 0.1]" ,"LP");
        lSigma->Draw("same");
      
        // Means
        paramCan->cd(2);
        grMean[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grMean[i][iCent]->Draw("p");
      
        // Amps
        paramCan->cd(3);
        grAmp[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grAmp[i][iCent]->Draw("p");
      
        // Kurtosis
        paramCan->cd(4);
        grKurt[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grKurt[i][iCent]->Draw("p");
      
        // skewness
        paramCan->cd(5);
        grSkew[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grSkew[i][iCent]->Draw("p");
      
        // Integral
        paramCan->cd(6);
        grInt[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grInt[i][iCent]->Draw("p");
        
        // Chi2/ndf
        paramCan->cd(8);
        grChi2[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grChi2[i][iCent]->Draw("p");
      
        // fill the graphs in tobjarray
        paramsPerCentArr.AddAt(paramCan,iCent);   
        if (iIter>=nIter-2) paramCan->SaveAs(Form("%s_params_iter%d_cent%d.png",parType.Data(),iIter,iCent));
       
      }
      
      streamer->GetFile()->cd();
      paramsPerCentArr.Write(Form("%s_FitParams_iter%d",parType.Data(),iIter) ,TObject::kSingleKey);  
    
    } // end of iter loop
  } // end of particle type loop
  
  delete streamer;
}
//____________________________________________________________________________________________________________
void InitInitials(){
  
  // 
  // Initialise histograms to be used for binning
  //
  
  // Create the histograms to be used in the binning of eta, cent and momentum
  hEta  =  new TH1D("hEta" ,"Eta Bins" ,nEtabins   ,etaMin, etaMax );
  hMom  =  new TH1D("hMom" ,"Mom Bins" ,ptNbins    ,ptMin,  ptMax );
  hCent =  new TH1D("hCent","Cent Bins",ncentbins-1, xCentBins );
  hEta  -> FillRandom("gaus",100000);
  hCent -> FillRandom("gaus",100000);
  hMom  -> FillRandom("gaus",100000);
               
}


