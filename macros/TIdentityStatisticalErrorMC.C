#include <TFile.h>
#include "TBranch.h"
#include <iostream>
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TTreeStream.h"
#include "TFrame.h"
#include "TROOT.h"
#include "TKey.h"
#include "TVirtualFitter.h"
#include "TClass.h"
#include "TF1.h"
#include "TTree.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TStatToolkit.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "AliXRDPROOFtoolkit.h"
#include "TChain.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::cout;
using std::setw;


void GraphShade();
void ErrorBandTH1();
void AnalyseHistsInFile(TString file);
void FitHistogram(TH1D *hClean);
void AnalyseWs(TString wFile, TString wSliceFile);
void ProcessErrorCalculationVScent(TString momentsFile,Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp);
void ProcessErrorCalculationVSmomentum(TString momentsFile, Double_t centDown, Double_t centUp);
void ProcessErrorCalculationVSeta(TString momentsFile,Double_t centDown, Double_t centUp, Double_t pDown, Double_t pUp);
void MakeRatioPlots(TString mfile,Int_t scanType, Double_t pDown, Double_t pUp, Double_t etaDown, Double_t etaUp);
void RatioScan(TString etaScanFile,TString momScanFile );
void SetYTitleName(TString branchName);
void MultiplydNch(TString momFile, Int_t flag);
void TIdentityStatisticalErrorMCNew(TString momentsFile, Int_t dataType, Int_t nSubsample, Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp);
Double_t GetNuDyn(Double_t a, Double_t b, Double_t a2, Double_t b2, Double_t ab);
Double_t GetNetParCumulant2(Double_t a, Double_t b, Double_t a2, Double_t b2, Double_t ab);


void FastPlotVScent(Int_t style, Int_t col, TString momentsFile, Int_t dataType, TString var, Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp, Int_t same);
void FastPlotVSeta(Int_t style, Int_t col, TString momentsFile, Int_t dataType, TString var, Double_t centDown, Double_t pDown, Double_t pUp, Int_t same);
void FastPlotVSmom(Int_t style,Int_t col,TString momentsFile,Int_t dataType,TString var,Double_t centDown,Int_t same);

TGraphErrors * RatioPlot(TTreeSRedirector *ratiofile, TGraphErrors *g1, TGraphErrors *g2, Int_t scanType, Double_t pDown, Double_t pUp, Double_t etaDown, Double_t etaUp);
Double_t       CalculateStatErrors(TGraphErrors *grS);
TCanvas      * TwoScales(TTree *t);
TCanvas      * OverlayPads(TTree *t);

// ******************************* Modification region ******************************************
Bool_t lowStat       = kFALSE;    // USed only for W analysis
Bool_t pp            = kFALSE;
Bool_t fastGen       = kFALSE;
TTreeSRedirector *MCssMomets=NULL;
TTreeSRedirector *statResults=NULL;
TTreeSRedirector *ratioresults=NULL;
TTreeSRedirector *correctionLookup=NULL;
TString yTitleLatex;                  // title name should be wrtiiten in latex format

// Double_t eff[4] = {1.,0.67, 0.47, 0.64};
Double_t eff[4] = {0.6, 0.6, 0.6, 0.6};
Double_t pDownArray[4]   = {0.2, 0.2, 0.3, 0.3};
Double_t pUpArray[4]     = {1.5, 1.8, 1.5, 1.8};
Double_t etaDownArray[2] = {-0.5, -0.8};
Double_t etaUpArray[2]   = { 0.5,  0.8};
Double_t centDownArray[9]= {0.,5. ,10.,20.,30.,40.,50.,60.,70.};
Double_t centUpArray[9]  = {5.,10.,20.,30.,40.,50.,60.,70.,80.};

Double_t dnchArray[9]    = {1601,1294,966,649,426,261,149,76,35};  // eta <0.5
Double_t dnchArrayErr[9] = {60,49,37,23,15,9,6,4,2};
Double_t nPartArr[9]     = {382.8, 329.7, 260.5, 186.4, 128.9, 85.0, 52.8, 30.0, 15.8};
Double_t nPartArrErr[9]  = {3.1, 4.6, 4.4, 3.9, 3.3, 2.6, 2.0, 1.3, 0.6};
Double_t dnchArrayHIJING[9]    = {1725.8, 1414.1, 1064.7, 725.9, 477.9, 297.5, 171.5, 89.1, 41.1 };  
Double_t dnchArrayHIJINGErr[9] = {60,49,37,23,15,9,6,4,2};
Double_t dnchArrayAMPT[9]      = {1592.5, 1291.8, 962.6, 646.4, 422.1, 260.9, 148.6, 76.11, 34.35};  // eta <0.5
Double_t dnchArrayAMPTErr[9]   = {60,49,37,23,15,9,6,4,2};
// ******************************* Modification region ******************************************

enum momentType {kPi=0,kKa=1,kPr=2,kPiPi=3,kKaKa=4,kPrPr=5,kPiKa=6,kPiPr=7,kKaPr=8,};
enum momentTypeUnlike {
    kPiPosPiNeg=0,
    kPiPosKaNeg=1,
    kPiPosPrNeg=2,
    kKaPosPiNeg=3,
    kKaPosKaNeg=4,
    kKaPosPrNeg=5,
    kPrPosPiNeg=6,
    kPrPosKaNeg=7,
    kPrPosPrNeg=8,
};
enum nudynType {
    kNudyn_PiKa=0,
    kNudyn_PiPr=1,
    kNudyn_KaPr=2,
    
    kNudyn_PiPosKaPos=3,
    kNudyn_PiPosPrPos=4,
    kNudyn_KaPosPrPos=5,
    
    kNudyn_PiNegKaNeg=6,
    kNudyn_PiNegPrNeg=7,
    kNudyn_KaNegPrNeg=8,
    
    kNudyn_PiPosPiNeg=9,
    kNudyn_PiPosKaNeg=10,
    kNudyn_PiPosPrNeg=11,
    kNudyn_KaPosPiNeg=12,
    kNudyn_KaPosKaNeg=13,
    kNudyn_KaPosPrNeg=14,
    kNudyn_PrPosPiNeg=15,
    kNudyn_PrPosKaNeg=16,
    kNudyn_PrPosPrNeg=17,
};
enum netCumType {
    kNetPi2=0,
    kNetKa2=1,
    kNetPr2=2,
};



void TIdentityStatisticalErrorMC(TString momentsFile, Int_t dataType, Int_t nSubsample, Double_t etaDown=-0.8, Double_t etaUp=0.8, Double_t pDown=0.2, Double_t pUp=3.){

  //
  // Calculate statistical errors for MCgen and MCrec form "MC_recgen_moments.root" using subsample method
  // 
                                                                                       
  // To RUN
  /*                                        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/MC_genVSrec/SSanalysis/EtaMomentumScan/cRows80/centBinWidth_10_all
  tpcdev
  aliroot -l -b
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  TIdentityStatisticalErrorMC("RecGen_moments.root",0,21)
  
  TIdentityStatisticalErrorMC("RecGen_moments.root",0,1)

  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C .");

  
  // Output file
  MCssMomets = new TTreeSRedirector(Form("MomentsTree_DT_%d_eta_%2.1f-%2.1f_mom_%2.1f-%2.1f.root",dataType,etaDown,etaUp,pDown,pUp),"recreate");

  // Get the trees 
  TTree *tree;
//   TString treeName = "mcMoments";
  TString treeName = (dataType==0) ? "mcRec": "mcGen";

  if (momentsFile.Contains(".root")){
    TFile *f = TFile::Open(momentsFile);
    cout << " get the ttree " << endl;
    tree = (TTree*)f->Get(treeName);
  } else {
    TChain *chain;
    cout << " make the chain " << endl;
    chain = AliXRDPROOFtoolkit::MakeChainRandom(momentsFile,treeName,0,500000,0);  
    chain -> SetCacheSize(100000000000000000);
    tree = (TTree*)chain;
  }
  if (!tree) return;
  
  // Second array to hold centDown
  Int_t nPoints = sizeof(centDownArray)/sizeof(centDownArray[0]);   
  cout << " nPoints == " << nPoints << endl;   
  
     // loop over centrality bins
  for(Int_t subsample=0; subsample<nSubsample;subsample++){
    cout << " =========== subsample = " << subsample << " =========== " << endl;
    for (Int_t icent=0;icent<nPoints;icent++){
        
    // Moment arrays for each centrality bin
      Float_t pTruth = 0;
      Float_t pT=0;
      Double_t firstMoments[4];
      Double_t secondMoments[4][4];
      Double_t nudynArr[4][4];
      Double_t nudynTerm1[4][4];
      Double_t nudynTerm2[4][4];
      Double_t nudynTerm3[4][4];
      
      Double_t nudynTerm1a[4][4];
      Double_t nudynTerm1b[4][4];
      Double_t nudynTerm2a[4][4];
      Double_t nudynTerm2b[4][4];
      
      Double_t nudynDyn[4][4];
      Double_t nudynStat[4][4];
      
      for (Int_t i=0; i<4;i++){
        firstMoments[i]=0.;
        eff[i]=1.;
        for (Int_t j=0;j<4;j++){
          secondMoments[i][j]=0.;
          nudynArr[i][j]=0.;
          nudynTerm1[i][j]=0.;
          nudynTerm2[i][j]=0.;
          nudynTerm3[i][j]=0.;
          nudynTerm1a[i][j]=0.;
          nudynTerm1b[i][j]=0.;
          nudynTerm2a[i][j]=0.;
          nudynTerm2b[i][j]=0.; 
          
          nudynDyn[i][j]=0.; 
          nudynStat[i][j]=0.; 
          
        }
      }
      TCut centCut   = Form("centDown==%f && isample==%d && dataType==%d",centDownArray[icent],subsample,dataType);
      TCut etaMomCut = Form("etaDown==%f && etaUp==%f && pDown==%f && pUp==%f"   ,etaDown, etaUp, pDown, pUp);

      // get list of the branches to loop over them
      TObjArray *branchArr   =  (TObjArray*)tree  ->GetListOfBranches(); 
      Int_t nEntries = branchArr->GetEntriesFast();
      for (Int_t iBranch=0; iBranch<nEntries; ++iBranch){
  
        TString branchName = branchArr->At(iBranch)->GetName();
        if ( branchName.Contains("cent")      || branchName.Contains("eta")      || branchName.Contains("ptot") ) continue;
        if ( branchName.Contains("de")        || branchName.Contains("el")       || branchName.Contains("mu") ) continue;
        if ( branchName.Contains("isample")   || branchName.Contains("dataType") || branchName.Contains("event") ) continue;
        if ( branchName.Contains("pDown")     || branchName.Contains("pUp")      || branchName.Contains("pBin") ) continue;
        if ( branchName.Contains("mctrackcount") || branchName.Contains("trackcount") ) continue;
        
    // plot the graph hist for each moment and get mean
        tree->Draw(branchName,centCut && etaMomCut,"goff");
        TH1D *htemp = (TH1D*)tree->GetHistogram()->Clone();
    
    // get first moments
        if (branchName.Contains("pTruth"))  pTruth           = htemp->GetMean();
        if (branchName.Contains("pT"))      pT               = htemp->GetMean();
        if (branchName.Contains("piCount")) firstMoments[1]  = htemp->GetMean();
        if (branchName.Contains("kaCount")) firstMoments[2]  = htemp->GetMean();
        if (branchName.Contains("prCount")) firstMoments[3]  = htemp->GetMean();
        if (branchName.Contains("pipi")) secondMoments[1][1] = htemp->GetMean();
        if (branchName.Contains("kaka")) secondMoments[2][2] = htemp->GetMean();
        if (branchName.Contains("prpr")) secondMoments[3][3] = htemp->GetMean();
        if (branchName.Contains("pika")) secondMoments[1][2] = htemp->GetMean();
        if (branchName.Contains("pipr")) secondMoments[1][3] = htemp->GetMean();
        if (branchName.Contains("kapr")) secondMoments[2][3] = htemp->GetMean();
        delete htemp; 
      }
  
    // make the matrix symmetric
      for (Int_t i=0; i<4;i++){
        for (Int_t j=0;j<4;j++){
          if (secondMoments[i][j]!=0.) secondMoments[j][i]=secondMoments[i][j];
        }
      }
  
    // Calculate nudyn
      for (Int_t inudyn=0; inudyn<4;inudyn++){
        for (Int_t jnudyn=0;jnudyn<4;jnudyn++){
          if ((firstMoments[inudyn]!=0) && (firstMoments[jnudyn]!=0)){
            nudynArr[inudyn][jnudyn] =  secondMoments[inudyn][inudyn]/firstMoments[inudyn]/firstMoments[inudyn];
            nudynArr[inudyn][jnudyn] += secondMoments[jnudyn][jnudyn]/firstMoments[jnudyn]/firstMoments[jnudyn];
            nudynArr[inudyn][jnudyn] -= 2.*secondMoments[inudyn][jnudyn]/firstMoments[inudyn]/firstMoments[jnudyn];
            nudynArr[inudyn][jnudyn] -= (1./firstMoments[inudyn]) + (1./firstMoments[jnudyn]);  // TO FIX
     
            nudynTerm1[inudyn][jnudyn]= secondMoments[inudyn][inudyn]/firstMoments[inudyn]/firstMoments[inudyn]-(1./firstMoments[inudyn]);
            nudynTerm2[inudyn][jnudyn]= secondMoments[jnudyn][jnudyn]/firstMoments[jnudyn]/firstMoments[jnudyn]-(1./firstMoments[jnudyn]);
            nudynTerm3[inudyn][jnudyn]= 2.*secondMoments[inudyn][jnudyn]/firstMoments[inudyn]/firstMoments[jnudyn];
          
            nudynTerm1b[inudyn][jnudyn]=1./firstMoments[inudyn];
            nudynTerm1a[inudyn][jnudyn]=secondMoments[inudyn][inudyn]/firstMoments[inudyn]/firstMoments[inudyn];
            nudynTerm2b[inudyn][jnudyn]=1./firstMoments[jnudyn];
            nudynTerm2a[inudyn][jnudyn]=secondMoments[jnudyn][jnudyn]/firstMoments[jnudyn]/firstMoments[jnudyn]; 
          
            nudynDyn[inudyn][jnudyn] = secondMoments[inudyn][inudyn]/firstMoments[inudyn]/firstMoments[inudyn];
            nudynDyn[inudyn][jnudyn] += secondMoments[jnudyn][jnudyn]/firstMoments[jnudyn]/firstMoments[jnudyn];
            nudynDyn[inudyn][jnudyn] -= 2.*secondMoments[inudyn][jnudyn]/firstMoments[inudyn]/firstMoments[jnudyn];
            nudynStat[inudyn][jnudyn] = (1./firstMoments[inudyn] + 1./firstMoments[jnudyn]);
          
          } 
        }
      }

      Double_t pBin   =(pDown+pUp)/2.;
      Double_t etaBin =(etaDown+etaUp)/2.;
      Double_t centBin=(centDownArray[icent]+centUpArray[icent])/2.;  
      cout << "etaRange = " << etaDown << " - "     << etaUp << "    pRange = " << pDown << " - " << pUp << " ------------ ";
      cout << "etaBin = " << etaBin << "   pBin = " << pBin << "    centBin = " << centBin << endl;
       
      // dump everything into ttree   
      MCssMomets -> GetFile()->cd();
      *MCssMomets << "mcMoments" <<
          
          "dataType="     << dataType            <<
          "subsample="    << subsample           <<
          "pDown="        << pDown               << 
          "pUp="          << pUp                 << 
          "pBin="         << pBin                << 
          "etaDown="      << etaDown             << 
          "etaUp="        << etaUp               << 
          "etaBin="       << etaBin              << 
          "centDown="     << centDownArray[icent]<< 
          "centUp="       << centUpArray[icent]  << 
          "centBin="      << centBin             << 
          "pTruth="       << pTruth              << 
          "pT="           << pT                  <<
          
          "pipiNuDyn="    << nudynArr[1][1]      <<
          "kakaNuDyn="    << nudynArr[2][2]      <<
          "prprNuDyn="    << nudynArr[3][3]      << 
          "pikaNuDyn="    << nudynArr[1][2]      <<
          "piprNuDyn="    << nudynArr[1][3]      <<
          "kaprNuDyn="    << nudynArr[2][3]      <<
        
          "pi="        << firstMoments[1]     << 
          "ka="        << firstMoments[2]     << 
          "pr="        << firstMoments[3]     << 
          "pipi="      << secondMoments[1][1] << 
          "kaka="      << secondMoments[2][2] << 
          "prpr="      << secondMoments[3][3] << 
          "pika="      << secondMoments[1][2] << 
          "pipr="      << secondMoments[1][3] << 
          "kapr="      << secondMoments[2][3] << 
           
//           "pipiNuDynDyn="    << nudynDyn[1][1]      <<
//           "kakaNuDynDyn="    << nudynDyn[2][2]      <<
//           "prprNuDynDyn="    << nudynDyn[3][3]      << 
//           "pikaNuDynDyn="    << nudynDyn[1][2]      <<
//           "piprNuDynDyn="    << nudynDyn[1][3]      <<
//           "kaprNuDynDyn="    << nudynDyn[2][3]      <<
          
//           "pipiNuStat="    << nudynStat[1][1]      <<
//           "kakaNuStat="    << nudynStat[2][2]      <<
//           "prprNuStat="    << nudynStat[3][3]      << 
//           "pikaNuStat="    << nudynStat[1][2]      <<
//           "piprNuStat="    << nudynStat[1][3]      <<
//           "kaprNuStat="    << nudynStat[2][3]      <<
          
//           "pipiNuTerm1="    << nudynTerm1[1][1]      <<
//           "kakaNuTerm1="    << nudynTerm1[2][2]      <<
//           "prprNuTerm1="    << nudynTerm1[3][3]      << 
//           "pikaNuTerm1="    << nudynTerm1[1][2]      <<
//           "piprNuTerm1="    << nudynTerm1[1][3]      <<
//           "kaprNuTerm1="    << nudynTerm1[2][3]      <<
          
//           "pipiNuTerm1a="   << nudynTerm1a[1][1]      <<
//           "kakaNuTerm1a="   << nudynTerm1a[2][2]      <<
//           "prprNuTerm1a="   << nudynTerm1a[3][3]      << 
//           "pikaNuTerm1a="   << nudynTerm1a[1][2]      <<
//           "piprNuTerm1a="   << nudynTerm1a[1][3]      <<
//           "kaprNuTerm1a="   << nudynTerm1a[2][3]      <<
          
//           "pipiNuTerm1b="   << nudynTerm1b[1][1]      <<
//           "kakaNuTerm1b="   << nudynTerm1b[2][2]      <<
//           "prprNuTerm1b="   << nudynTerm1b[3][3]      << 
//           "pikaNuTerm1b="   << nudynTerm1b[1][2]      <<
//           "piprNuTerm1b="   << nudynTerm1b[1][3]      <<
//           "kaprNuTerm1b="   << nudynTerm1b[2][3]      ;
      
//       *MCssMomets << "mcMoments" <<
        
//           "pipiNuTerm2="    << nudynTerm2[1][1]      <<
//           "kakaNuTerm2="    << nudynTerm2[2][2]      <<
//           "prprNuTerm2="    << nudynTerm2[3][3]      << 
//           "pikaNuTerm2="    << nudynTerm2[1][2]      <<
//           "piprNuTerm2="    << nudynTerm2[1][3]      <<
//           "kaprNuTerm2="    << nudynTerm2[2][3]      <<
          
//           "pipiNuTerm2a="   << nudynTerm2a[1][1]      <<
//           "kakaNuTerm2a="   << nudynTerm2a[2][2]      <<
//           "prprNuTerm2a="   << nudynTerm2a[3][3]      << 
//           "pikaNuTerm2a="   << nudynTerm2a[1][2]      <<
//           "piprNuTerm2a="   << nudynTerm2a[1][3]      <<
//           "kaprNuTerm2a="   << nudynTerm2a[2][3]      <<
          
//           "pipiNuTerm2b="   << nudynTerm2b[1][1]      <<
//           "kakaNuTerm2b="   << nudynTerm2b[2][2]      <<
//           "prprNuTerm2b="   << nudynTerm2b[3][3]      << 
//           "pikaNuTerm2b="   << nudynTerm2b[1][2]      <<
//           "piprNuTerm2b="   << nudynTerm2b[1][3]      <<
//           "kaprNuTerm2b="   << nudynTerm2b[2][3]      <<
        
//           "pipiNuTerm3="    << nudynTerm3[1][1]      <<
//           "kakaNuTerm3="    << nudynTerm3[2][2]      <<
//           "prprNuTerm3="    << nudynTerm3[3][3]      << 
//           "pikaNuTerm3="    << nudynTerm3[1][2]      <<
//           "piprNuTerm3="    << nudynTerm3[1][3]      <<
//           "kaprNuTerm3="    << nudynTerm3[2][3]      <<
          
          "\n"; 
        
    } // end of centrality loop
  } // end of subsample loop
  delete MCssMomets;
  
}
// =====================================================================================================
void TIdentityStatisticalErrorMCNew(TString momentsFile, Int_t dataType, Int_t nSubsample, Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp){

  //
  // Calculate statistical errors for MCgen and MCrec form "MC_recgen_moments.root" using subsample method
  // 
                                                                                       
  // To RUN
  /*     
        
  TString input = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN1/11a10abis_ResRhoPhiOff_cRows_80_10EtaBin_mombin20MeV/mergedPeriods/AnalysisResults.list"    
        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN1/11a10abis_ResRhoPhiOff_cRows_80_10EtaBin_mombin20MeV
  aliroot -l -b
  TString input = "AnalysisResults.root"
  //   TString input = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/RUN1/11a10abis_ResRhoPhiOff_cRows_80_10EtaBin_mombin20MeV/mergedPeriods/AnalysisResults.list"    
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  TIdentityStatisticalErrorMCNew(input,1,1,-1,1,0.2,1.5)
  TIdentityStatisticalErrorMCNew(input,1,1,-0.8,0.8,0.2,1.5)
  TIdentityStatisticalErrorMCNew(input,1,1,-0.6,0.6,0.2,1.5)
  TIdentityStatisticalErrorMCNew(input,1,1,-0.4,0.4,0.2,1.5)
  TIdentityStatisticalErrorMCNew(input,1,1,-0.2,0.2,0.2,1.5)
  
  TFile ("MomentsTree_DT_1_eta_-1.0-1.0_mom_0.2-1.5.root")
  mcMoments->SetMarkerStyle(20); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("momentPos.fElements[0]+momentNeg.fElements[0]:centBin","","")
  mcMoments->SetMarkerStyle(24); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("netpar.fElements[0]:centBin","","same")
  mcMoments->SetMarkerStyle(21); mcMoments->SetMarkerColor(kGreen+1); 
  mcMoments->Draw("noResmomentPos.fElements[0]+noResmomentNeg.fElements[0]:centBin","","same")
  mcMoments->SetMarkerStyle(25); mcMoments->SetMarkerColor(kGreen+1); 
  mcMoments->Draw("noResnetpar.fElements[0]:centBin","","same")

  new TCanvas
  mcMoments->SetMarkerStyle(20); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("momentPos.fElements[1]+momentNeg.fElements[1]:centBin","","")
  mcMoments->SetMarkerStyle(24); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("netpar.fElements[1]:centBin","","same")
  mcMoments->SetMarkerStyle(21); mcMoments->SetMarkerColor(kGreen+1); 
  mcMoments->Draw("noResmomentPos.fElements[1]+noResmomentNeg.fElements[1]:centBin","","same")
  mcMoments->SetMarkerStyle(25); mcMoments->SetMarkerColor(kGreen+1); 
  mcMoments->Draw("noResnetpar.fElements[1]:centBin","","same")

  new TCanvas
  mcMoments->SetMarkerStyle(20); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("momentPos.fElements[2]+momentNeg.fElements[2]:centBin","","")
  mcMoments->SetMarkerStyle(24); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("netpar.fElements[2]:centBin","","same")
  mcMoments->SetMarkerStyle(21); mcMoments->SetMarkerColor(kGreen+1); 
  mcMoments->Draw("noResmomentPos.fElements[2]+noResmomentNeg.fElements[2]:centBin","","same")
  mcMoments->SetMarkerStyle(25); mcMoments->SetMarkerColor(kGreen+1); 
  mcMoments->Draw("noResnetpar.fElements[2]:centBin","","same")

  
  new TCanvas
  mcMoments->SetMarkerStyle(20); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("nudyn.fElements[0]:centBin","","")
  mcMoments->SetMarkerStyle(24); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("noResnudyn.fElements[0]:centBin","","same")
 
  new TCanvas
  mcMoments->SetMarkerStyle(20); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("nudyn.fElements[1]:centBin","","")
  mcMoments->SetMarkerStyle(24); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("noResnudyn.fElements[1]:centBin","","same")

  new TCanvas
  mcMoments->SetMarkerStyle(20); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("nudyn.fElements[2]:centBin","","")
  mcMoments->SetMarkerStyle(24); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("noResnudyn.fElements[2]:centBin","","same")

  new TCanvas
  mcMoments->SetMarkerStyle(20); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("nudyn.fElements[9]:centBin","","")
  mcMoments->SetMarkerStyle(24); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("noResnudyn.fElements[9]:centBin","","same")

  new TCanvas
  mcMoments->SetMarkerStyle(20); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("nudyn.fElements[13]:centBin","","")
  mcMoments->SetMarkerStyle(24); mcMoments->SetMarkerColor(kRed+1); 
  mcMoments->Draw("noResnudyn.fElements[13]:centBin","","same")

  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C .");

  
  // Output file
  MCssMomets = new TTreeSRedirector(Form("MomentsTree_DT_%d_eta_%2.1f-%2.1f_mom_%2.1f-%2.1f.root",dataType,etaDown,etaUp,pDown,pUp),"recreate");

  // Get the trees 
  TTree *tree;
//   TString treeName = "mcMoments";
  TString treeName = (dataType==0) ? "mcRec": "mcGen";

  if (momentsFile.Contains(".root")){
    TFile *f = TFile::Open(momentsFile);
    cout << " get the ttree " << treeName << endl;
    tree = (TTree*)f->Get(treeName);
  } else {
    TChain *chain;
    cout << " make the chain " << treeName << endl;
    chain = AliXRDPROOFtoolkit::MakeChainRandom(momentsFile,treeName,0,500000,0);  
    chain -> SetCacheSize(100000000000000000);
    tree = (TTree*)chain;
  }
  if (!tree) return;
  
  // Second array to hold centDown
  Int_t nPoints = sizeof(centDownArray)/sizeof(centDownArray[0]);   
  cout << " nPoints == " << nPoints << endl;   
  const Int_t nMoments = 9;

  // loop over centrality bins
  for(Int_t subsample=0; subsample<nSubsample;subsample++){
    cout << " =========== subsample = " << subsample << " =========== " << endl;
    for (Int_t icent=0;icent<nPoints;icent++){
        
      TVectorF moment(nMoments);
      TVectorF momentPos(nMoments);
      TVectorF momentNeg(nMoments);
      TVectorF momentCross(nMoments);
      TVectorF noResmoment(nMoments);
      TVectorF noResmomentPos(nMoments);
      TVectorF noResmomentNeg(nMoments);
      TVectorF noResmomentCross(nMoments);
      // initialize counters 
      for(Int_t i=0;i<nMoments; i++){  
          moment[i]=0.; 
          momentPos[i]=0.;  
          momentNeg[i]=0.; 
          momentCross[i]=0.;
          noResmoment[i]=0.; 
          noResmomentPos[i]=0.;  
          noResmomentNeg[i]=0.; 
          noResmomentCross[i]=0.;
      }
      
      // get the moments from tree   
      TCut centCut   = Form("abs(centDown-%f)<0.01 && isample==%d && dataType==%d",centDownArray[icent],subsample,dataType);
      TCut etaMomCut = Form("abs(etaUp-%f)<0.001 && abs(pDown-%f)<0.001 && abs(pUp-%f)<0.001", etaUp, pDown, pUp);

      // get list of the branches to loop over them
      TObjArray *branchArr   =  (TObjArray*)tree  ->GetListOfBranches(); 
      Int_t nEntries = branchArr->GetEntriesFast();
      for (Int_t ielement=0; ielement<nMoments; ++ielement){
  
        // plot the graph hist for each moment and get mean
        tree->Draw(Form("moment.fElements[%d]",ielement),centCut && etaMomCut,"goff");
        TH1D *hmoment = (TH1D*)tree->GetHistogram()->Clone();
        moment[ielement]=hmoment->GetMean();
        tree->Draw(Form("momentPos.fElements[%d]",ielement),centCut && etaMomCut,"goff");
        TH1D *hmomentPos = (TH1D*)tree->GetHistogram()->Clone();
        momentPos[ielement]=hmomentPos->GetMean();
        tree->Draw(Form("momentNeg.fElements[%d]",ielement),centCut && etaMomCut,"goff");
        TH1D *hmomentNeg = (TH1D*)tree->GetHistogram()->Clone();
        momentNeg[ielement]=hmomentNeg->GetMean();
        tree->Draw(Form("momentCross.fElements[%d]",ielement),centCut && etaMomCut,"goff");
        TH1D *hmomentCross = (TH1D*)tree->GetHistogram()->Clone();
        momentCross[ielement]=hmomentCross->GetMean();
       
        tree->Draw(Form("noResmoment.fElements[%d]",ielement),centCut && etaMomCut,"goff");
        TH1D *hnoResmoment = (TH1D*)tree->GetHistogram()->Clone();
        noResmoment[ielement]=hnoResmoment->GetMean();
        tree->Draw(Form("noResmomentPos.fElements[%d]",ielement),centCut && etaMomCut,"goff");
        TH1D *hnoResmomentPos = (TH1D*)tree->GetHistogram()->Clone();
        noResmomentPos[ielement]=hnoResmomentPos->GetMean();
        tree->Draw(Form("noResmomentNeg.fElements[%d]",ielement),centCut && etaMomCut,"goff");
        TH1D *hnoResmomentNeg = (TH1D*)tree->GetHistogram()->Clone();
        noResmomentNeg[ielement]=hnoResmomentNeg->GetMean();
        tree->Draw(Form("noResmomentCross.fElements[%d]",ielement),centCut && etaMomCut,"goff");
        TH1D *hnoResmomentCross = (TH1D*)tree->GetHistogram()->Clone();
        noResmomentCross[ielement]=hnoResmomentCross->GetMean();
        
        delete hmoment; 
        delete hmomentPos; 
        delete hmomentNeg; 
        delete hmomentCross; 
        delete hnoResmoment; 
        delete hnoResmomentPos; 
        delete hnoResmomentNeg; 
        delete hnoResmomentCross; 

      }
      
      // Calculate nudyn
      TVectorF nudyn(2*nMoments); for(Int_t i=0;i<2*nMoments; i++) nudyn[i]=0.;
      TVectorF noResnudyn(2*nMoments); for(Int_t i=0;i<2*nMoments; i++) nudyn[i]=0.;
      nudyn[kNudyn_PiKa]       = GetNuDyn(moment[kPi],moment[kKa],moment[kPiPi],moment[kKaKa],moment[kPiKa]);
      nudyn[kNudyn_PiPr]       = GetNuDyn(moment[kPi],moment[kPr],moment[kPiPi],moment[kPrPr],moment[kPiPr]);
      nudyn[kNudyn_KaPr]       = GetNuDyn(moment[kKa],moment[kPr],moment[kKaKa],moment[kPrPr],moment[kKaPr]);
      nudyn[kNudyn_PiPosKaPos] = GetNuDyn(momentPos[kPi],momentPos[kKa],momentPos[kPiPi],momentPos[kKaKa],momentPos[kPiKa]);
      nudyn[kNudyn_PiPosPrPos] = GetNuDyn(momentPos[kPi],momentPos[kPr],momentPos[kPiPi],momentPos[kPrPr],momentPos[kPiPr]);
      nudyn[kNudyn_KaPosPrPos] = GetNuDyn(momentPos[kKa],momentPos[kPr],momentPos[kKaKa],momentPos[kPrPr],momentPos[kKaPr]);
      nudyn[kNudyn_PiNegKaNeg] = GetNuDyn(momentNeg[kPi],momentNeg[kKa],momentNeg[kPiPi],momentNeg[kKaKa],momentNeg[kPiKa]);
      nudyn[kNudyn_PiNegPrNeg] = GetNuDyn(momentNeg[kPi],momentNeg[kPr],momentNeg[kPiPi],momentNeg[kPrPr],momentNeg[kPiPr]);
      nudyn[kNudyn_KaNegPrNeg] = GetNuDyn(momentNeg[kKa],momentNeg[kPr],momentNeg[kKaKa],momentNeg[kPrPr],momentNeg[kKaPr]);
      nudyn[kNudyn_PiPosPiNeg] = GetNuDyn(momentPos[kPi],momentNeg[kPi],momentPos[kPiPi],momentNeg[kPiPi],momentCross[kPiPosPiNeg]);
      nudyn[kNudyn_PiPosKaNeg] = GetNuDyn(momentPos[kPi],momentNeg[kKa],momentPos[kPiPi],momentNeg[kKaKa],momentCross[kPiPosKaNeg]);
      nudyn[kNudyn_PiPosPrNeg] = GetNuDyn(momentPos[kPi],momentNeg[kPr],momentPos[kPiPi],momentNeg[kPrPr],momentCross[kPiPosPrNeg]);
      nudyn[kNudyn_KaPosPiNeg] = GetNuDyn(momentPos[kKa],momentNeg[kPi],momentPos[kKaKa],momentNeg[kPiPi],momentCross[kKaPosPiNeg]);
      nudyn[kNudyn_KaPosKaNeg] = GetNuDyn(momentPos[kKa],momentNeg[kKa],momentPos[kKaKa],momentNeg[kKaKa],momentCross[kKaPosKaNeg]);
      nudyn[kNudyn_KaPosPrNeg] = GetNuDyn(momentPos[kKa],momentNeg[kPr],momentPos[kKaKa],momentNeg[kPrPr],momentCross[kKaPosPrNeg]);
      nudyn[kNudyn_PrPosPiNeg] = GetNuDyn(momentPos[kPr],momentNeg[kPi],momentPos[kPrPr],momentNeg[kPiPi],momentCross[kPrPosPiNeg]);
      nudyn[kNudyn_PrPosKaNeg] = GetNuDyn(momentPos[kPr],momentNeg[kKa],momentPos[kPrPr],momentNeg[kKaKa],momentCross[kPrPosKaNeg]);
      nudyn[kNudyn_PrPosPrNeg] = GetNuDyn(momentPos[kPr],momentNeg[kPr],momentPos[kPrPr],momentNeg[kPrPr],momentCross[kPrPosPrNeg]);
      
      noResnudyn[kNudyn_PiKa]       = GetNuDyn(noResmoment[kPi],noResmoment[kKa],noResmoment[kPiPi],noResmoment[kKaKa],noResmoment[kPiKa]);
      noResnudyn[kNudyn_PiPr]       = GetNuDyn(noResmoment[kPi],noResmoment[kPr],noResmoment[kPiPi],noResmoment[kPrPr],noResmoment[kPiPr]);
      noResnudyn[kNudyn_KaPr]       = GetNuDyn(noResmoment[kKa],noResmoment[kPr],noResmoment[kKaKa],noResmoment[kPrPr],noResmoment[kKaPr]);
      noResnudyn[kNudyn_PiPosKaPos] = GetNuDyn(noResmomentPos[kPi],noResmomentPos[kKa],noResmomentPos[kPiPi],noResmomentPos[kKaKa],noResmomentPos[kPiKa]);
      noResnudyn[kNudyn_PiPosPrPos] = GetNuDyn(noResmomentPos[kPi],noResmomentPos[kPr],noResmomentPos[kPiPi],noResmomentPos[kPrPr],noResmomentPos[kPiPr]);
      noResnudyn[kNudyn_KaPosPrPos] = GetNuDyn(noResmomentPos[kKa],noResmomentPos[kPr],noResmomentPos[kKaKa],noResmomentPos[kPrPr],noResmomentPos[kKaPr]);
      noResnudyn[kNudyn_PiNegKaNeg] = GetNuDyn(noResmomentNeg[kPi],noResmomentNeg[kKa],noResmomentNeg[kPiPi],noResmomentNeg[kKaKa],noResmomentNeg[kPiKa]);
      noResnudyn[kNudyn_PiNegPrNeg] = GetNuDyn(noResmomentNeg[kPi],noResmomentNeg[kPr],noResmomentNeg[kPiPi],noResmomentNeg[kPrPr],noResmomentNeg[kPiPr]);
      noResnudyn[kNudyn_KaNegPrNeg] = GetNuDyn(noResmomentNeg[kKa],noResmomentNeg[kPr],noResmomentNeg[kKaKa],noResmomentNeg[kPrPr],noResmomentNeg[kKaPr]);
      noResnudyn[kNudyn_PiPosPiNeg] = GetNuDyn(noResmomentPos[kPi],noResmomentNeg[kPi],noResmomentPos[kPiPi],noResmomentNeg[kPiPi],noResmomentCross[kPiPosPiNeg]);
      noResnudyn[kNudyn_PiPosKaNeg] = GetNuDyn(noResmomentPos[kPi],noResmomentNeg[kKa],noResmomentPos[kPiPi],noResmomentNeg[kKaKa],noResmomentCross[kPiPosKaNeg]);
      noResnudyn[kNudyn_PiPosPrNeg] = GetNuDyn(noResmomentPos[kPi],noResmomentNeg[kPr],noResmomentPos[kPiPi],noResmomentNeg[kPrPr],noResmomentCross[kPiPosPrNeg]);
      noResnudyn[kNudyn_KaPosPiNeg] = GetNuDyn(noResmomentPos[kKa],noResmomentNeg[kPi],noResmomentPos[kKaKa],noResmomentNeg[kPiPi],noResmomentCross[kKaPosPiNeg]);
      noResnudyn[kNudyn_KaPosKaNeg] = GetNuDyn(noResmomentPos[kKa],noResmomentNeg[kKa],noResmomentPos[kKaKa],noResmomentNeg[kKaKa],noResmomentCross[kKaPosKaNeg]);
      noResnudyn[kNudyn_KaPosPrNeg] = GetNuDyn(noResmomentPos[kKa],noResmomentNeg[kPr],noResmomentPos[kKaKa],noResmomentNeg[kPrPr],noResmomentCross[kKaPosPrNeg]);
      noResnudyn[kNudyn_PrPosPiNeg] = GetNuDyn(noResmomentPos[kPr],noResmomentNeg[kPi],noResmomentPos[kPrPr],noResmomentNeg[kPiPi],noResmomentCross[kPrPosPiNeg]);
      noResnudyn[kNudyn_PrPosKaNeg] = GetNuDyn(noResmomentPos[kPr],noResmomentNeg[kKa],noResmomentPos[kPrPr],noResmomentNeg[kKaKa],noResmomentCross[kPrPosKaNeg]);
      noResnudyn[kNudyn_PrPosPrNeg] = GetNuDyn(noResmomentPos[kPr],noResmomentNeg[kPr],noResmomentPos[kPrPr],noResmomentNeg[kPrPr],noResmomentCross[kPrPosPrNeg]);

      //calculate net particles
      TVectorF netpar(2*nMoments); for(Int_t i=0;i<2*nMoments; i++) netpar[i]=0.;
      TVectorF noResnetpar(2*nMoments); for(Int_t i=0;i<2*nMoments; i++) noResnetpar[i]=0.;
      netpar[kNetPi2] = GetNetParCumulant2(momentPos[kPi],momentNeg[kPi],momentPos[kPiPi],momentNeg[kPiPi],momentCross[kPiPosPiNeg]);
      netpar[kNetKa2] = GetNetParCumulant2(momentPos[kKa],momentNeg[kKa],momentPos[kKaKa],momentNeg[kKaKa],momentCross[kKaPosKaNeg]);
      netpar[kNetPr2] = GetNetParCumulant2(momentPos[kPr],momentNeg[kPr],momentPos[kPrPr],momentNeg[kPrPr],momentCross[kPrPosPrNeg]);
      
      noResnetpar[kNetPi2] = GetNetParCumulant2(noResmomentPos[kPi],noResmomentNeg[kPi],noResmomentPos[kPiPi],noResmomentNeg[kPiPi],noResmomentCross[kPiPosPiNeg]);
      noResnetpar[kNetKa2] = GetNetParCumulant2(noResmomentPos[kKa],noResmomentNeg[kKa],noResmomentPos[kKaKa],noResmomentNeg[kKaKa],noResmomentCross[kKaPosKaNeg]);
      noResnetpar[kNetPr2] = GetNetParCumulant2(noResmomentPos[kPr],noResmomentNeg[kPr],noResmomentPos[kPrPr],noResmomentNeg[kPrPr],noResmomentCross[kPrPosPrNeg]);

      Double_t pBin   =(pDown+pUp)/2.;
      Double_t etaBin =(etaDown+etaUp)/2.;
      Double_t centBin=(centDownArray[icent]+centUpArray[icent])/2.;  
      cout << "etaRange = " << etaDown << " - "     << etaUp << "    pRange = " << pDown << " - " << pUp << " ------------ ";
      cout << "etaBin = " << etaBin << "   pBin = " << pBin << "    centBin = " << centBin << endl;
       
      // dump everything into ttree   
      MCssMomets -> GetFile()->cd();
      *MCssMomets << "mcMoments" <<
          
          "dataType="     << dataType            <<
          "subsample="    << subsample           <<
          "pDown="        << pDown               << 
          "pUp="          << pUp                 << 
          "pBin="         << pBin                << 
          "etaDown="      << etaDown             << 
          "etaUp="        << etaUp               << 
          "etaBin="       << etaBin              << 
          "centDown="     << centDownArray[icent]<< 
          "centUp="       << centUpArray[icent]  << 
          "centBin="      << centBin             << 
          
          "netpar.="       << &netpar <<             // second moments for particle+antiparticle
          "noResnetpar.="  << &noResnetpar <<             // second moments for particle+antiparticle
          "nudyn.="        << &nudyn <<             // second moments for particle+antiparticle
          "noResnudyn.="   << &noResnudyn <<             // second moments for particle+antiparticle
          "moment.="       << &moment <<             // second moments for particle+antiparticle
          "momentPos.="    << &momentPos <<          // second moment of positive particles
          "momentNeg.="    << &momentNeg <<          // second moment of negative particles
          "momentCross.="  << &momentCross <<        // second moment of unlikesign particles
          "noResmoment.="      << &noResmoment <<             // second moments for particle+antiparticle
          "noResmomentPos.="   << &noResmomentPos <<          // second moment of positive particles
          "noResmomentNeg.="   << &noResmomentNeg <<          // second moment of negative particles
          "noResmomentCross.=" << &noResmomentCross <<        // second moment of unlikesign particles
           
          
          "\n"; 
        
    } // end of centrality loop
  } // end of subsample loop
  delete MCssMomets;
  
}
// =====================================================================================================
void ProcessErrorCalculationVScent(TString momentsFile,Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp){

  //
  // Calculate statistical errors using subsample method
  // First merge all moment outputs
  // hadd Tree_TImoments_All.root cent_#/#/TI#.root 
  // 
                                                                                       
  // To RUN
  /*                                        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/test
  tpcdev
  aliroot -l
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  ProcessErrorCalculationVScent("MomentsTree.root",-0.8,0.8,0.2,1.5)
  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C .");

  
  // Output file --> Eta_-0.8_0.8_Mom_0.2_1.5
  TString statFileName = Form("Eta_%2.1f_%2.1f_Mom_%2.1f_%2.1f.root",etaDown,etaUp,pDown,pUp); 
  statResults = new TTreeSRedirector(statFileName,"recreate");

  // Get the trees
  TFile *fMoments = TFile::Open(momentsFile);  
  TTree *treeM    = (TTree*)fMoments -> Get("mcMoments"); 
  TTree *treeSS   = (TTree*)treeM->Clone();
  
  // GetList of branches for each tree
  TObjArray *branchM = (TObjArray*)treeM -> GetListOfBranches();
  Int_t nArrEntries  = branchM->GetEntriesFast();
  cout << " nArrEntries = " << nArrEntries << endl;
  TObjArray *momentArr    = new TObjArray(nArrEntries);  momentArr -> SetOwner(kTRUE); 
  TObjArray *subsamArr    = new TObjArray(nArrEntries);  subsamArr -> SetOwner(kTRUE); 

  TGraphErrors *grMoments[nArrEntries]; for (Int_t i=0;i<nArrEntries;i++) grMoments[i]=NULL;
  
  // Create all "moment vs cent" TGraphErrors
  for (Int_t idatatype=0; idatatype<2; idatatype++) {
    if (fastGen && idatatype==0) continue;
    Int_t objIndex = 0;
    for (Int_t i=0; i<nArrEntries; ++i){
        
      // Get the brach name 
      TString branchName = branchM->At(i)->GetName();
      TString drawMoment = branchM->At(i)->GetName();
      TString drawSubSam = branchM->At(i)->GetName();
   
      drawMoment.Append(":centBin");
      drawSubSam.Append(":subsample");
      TCut allStat   = Form("subsample==0 && dataType==%d",idatatype);
      TCut subsamCut = "subsample>0 && centBin==%f && dataType==%d";
      TCut etaMomCut = Form("etaDown==%f && etaUp==%f && pDown==%f && pUp==%f"   ,etaDown, etaUp, pDown, pUp);

      // look at only amplitude Mean and Sigma
      if ( (branchName.Contains("cent")  || branchName.Contains("subsample") || branchName.Contains("dataType")) ) continue;
      if ( (branchName.Contains("Up")    || branchName.Contains("Down")      || branchName.Contains("Bin")) ) continue;
      if ( (branchName.Contains("pT")) ) continue;
      if ( (branchName.Contains("pT")) ) continue;
      if ( (branchName.Contains("DynDyn")) ) continue;
      if ( (branchName.Contains("Stat")) ) continue;

    
      // create the graph object
      treeM->Draw(drawMoment,allStat && etaMomCut,"goff");
      Int_t N = treeM->GetSelectedRows();
    
      // Fix the ordering of centralitues in graph
      Double_t x[N];
      Double_t y[N];
      Int_t index[N];
      Double_t *tmpx = treeM->GetV2();  
      Double_t *tmpy = treeM->GetV1(); 
      TMath::Sort(N,tmpx,index);
      for (Int_t k=0; k<N;k++){ x[k] = tmpx[index[k]]; y[k] = tmpy[index[k]]; }
      
      SetYTitleName(branchName);
      grMoments[objIndex] = new TGraphErrors(N,x,y); 
      grMoments[objIndex]->SetDrawOption("a3");
      grMoments[objIndex]->SetName(branchName);
      grMoments[objIndex]->GetXaxis()->SetTitle("Centrality (%)");
      grMoments[objIndex]->GetYaxis()->SetTitle(yTitleLatex);
      grMoments[objIndex]->SetFillColor(kBlue);
      grMoments[objIndex]->SetMarkerStyle(21);
      grMoments[objIndex]->SetMarkerColor(kBlue);

      // Get the subsample
      TObjArray *subsamCentArr = new TObjArray(N);  subsamCentArr -> SetOwner(kTRUE); 
      TGraphErrors *grSS[N];
      for (Int_t is = 0; is<N; is++){
            
        treeSS->Draw(drawSubSam,Form(subsamCut,x[is],idatatype) && etaMomCut,"goff");
        grSS[is] = new TGraphErrors(treeSS->GetSelectedRows(),treeSS->GetV2(),treeSS->GetV1());
        grSS[is] -> SetMarkerStyle(20);
        TString sName = branchName;
        subsamCentArr -> SetName(sName.Append("_SS"));
        grSS[is]      -> SetName(sName.Append(Form("_cent%d",Int_t(x[is]))));      
        subsamCentArr -> AddAt(grSS[is],is);   
      }
      
      // Calculate the statistical errors
      Int_t ngrPoints = (sizeof x)/(sizeof x[0]);
      for (Int_t ierr = 0; ierr<ngrPoints; ierr++) {
        if (!grMoments[objIndex]) continue;
        grMoments[objIndex]->SetPointError(ierr,0,CalculateStatErrors(grSS[ierr]));
      }
      
      // Add graphs to ObjArrays
      momentArr -> AddAt(grMoments[objIndex],objIndex);
      subsamArr -> AddAt(subsamCentArr,objIndex);
      cout << "dataType = " << idatatype << "    moment " << objIndex << " is " << branchName << endl;
      objIndex++;
    
    }
    statResults->GetFile()->cd();
    momentArr->Write(Form("Moments_%d",idatatype),TObject::kSingleKey);
    subsamArr->Write(Form("SubSamples_%d",idatatype),TObject::kSingleKey);
  }
  
  delete statResults;
  fMoments->Close();
  
  if(!fastGen) MakeRatioPlots(statFileName,10,pDown,pUp,etaDown,etaUp);
 
}
// =====================================================================================================
void ProcessErrorCalculationVSmomentum(TString momentsFile,Double_t centDown, Double_t centUp){

  //
  // Calculate statistical errors using subsample method
  // First merge all moment outputs
  // hadd Tree_TImoments_All.root cent_#/#/TI#.root 
  // 
                                                                                       
  // To RUN
  /*                                        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/test
  tpcdev
  aliroot -l -b
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  ProcessErrorCalculationVSmomentum("rootFiles/MomScan_MomentsTree.root",0,5)
  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C .");
  
  // Output file
  TString statFileName = Form("MomScan_cent_%d_%d.root",Int_t(centDown),Int_t(centUp)); 
  statResults = new TTreeSRedirector(statFileName,"recreate");

  // Get the trees
  TFile *fMoments = TFile::Open(momentsFile);  
  TTree *treeM    = (TTree*)fMoments -> Get("mcMoments"); 
  TTree *treeSS   = (TTree*)treeM->Clone();
  
  // GetList of branches for each tree
  TObjArray *branchM = (TObjArray*)treeM -> GetListOfBranches();
  Int_t nArrEntries  = branchM->GetEntriesFast();
  cout << " nArrEntries = " << nArrEntries << endl;

//   Output TObjArrays
  TObjArray *momentArr    = new TObjArray(nArrEntries);  momentArr -> SetOwner(kTRUE); 
  TObjArray *subsamArr    = new TObjArray(nArrEntries);  subsamArr -> SetOwner(kTRUE); 

  TGraphErrors *grMoments[nArrEntries]; for (Int_t i=0;i<nArrEntries;i++) grMoments[i]=NULL;
  
  // Create all "moment vs cent" TGraphErrors
//   Int_t nCentBins = (sizeof centDownArray)/(sizeof centDownArray[0]);
  for (Int_t idatatype=0; idatatype<2; idatatype++) {
    Int_t objIndex = 0;
    for (Int_t i=0; i<nArrEntries; ++i){
        
    // Get the brach name 
      TString branchName = branchM->At(i)->GetName();
      TString drawMoment = branchM->At(i)->GetName();
      TString drawSubSam = branchM->At(i)->GetName();
   
      drawMoment.Append(":pBin");
      drawSubSam.Append(":subsample");
      
      TCut allStat   = Form("subsample==0&&dataType==%d",idatatype);
      TCut fixedCut  = Form("etaDown==-0.8&&etaUp==0.8&&(pUp-pDown)<0.5&&centDown==%f&&centUp==%f",centDown,centUp);
      TCut subsamCut = "subsample>0 && subsample<26 && pBin==%f && dataType==%d";

    // look at only amplitude Mean and Sigma
      if ( (branchName.Contains("cent")  || branchName.Contains("subsample") || branchName.Contains("dataType")) ) continue;
      if ( (branchName.Contains("Up")    || branchName.Contains("Down")      || branchName.Contains("Bin")) ) continue;
      if ( (branchName.Contains("DynDyn")|| branchName.Contains("Stat")      || branchName.Contains("Term")) ) continue;
      if ( (branchName.Contains("pipi")  || branchName.Contains("kaka")      || branchName.Contains("prpr")) ) continue;
      if ( (branchName.Contains("pT")) ) continue;

      
    // create the graph object
      treeM->Draw(drawMoment,allStat && fixedCut,"goff");
      Int_t N = treeM->GetSelectedRows();
    
    // Fix the ordering of centralitues in graph
      Double_t x[N]; Double_t y[N]; Int_t index[N];
      Double_t *tmpx = treeM->GetV2();  
      Double_t *tmpy = treeM->GetV1(); 
      TMath::Sort(N,tmpx,index);
      for (Int_t k=0; k<N;k++){ x[k] = tmpx[index[k]]; y[k] = tmpy[index[k]]; }

      SetYTitleName(branchName);
      grMoments[objIndex] = new TGraphErrors(N,x,y); 
      grMoments[objIndex]->SetDrawOption("a3");
      grMoments[objIndex]->SetName(branchName);
      grMoments[objIndex]->GetXaxis()->SetTitle("p (GeV/#it{c})");
      grMoments[objIndex]->GetYaxis()->SetTitle(yTitleLatex);
      grMoments[objIndex]->SetFillColor(kBlue);
      grMoments[objIndex]->SetMarkerStyle(21);
      grMoments[objIndex]->SetMarkerColor(kBlue);

    
    // Get the subsample
      TObjArray *subsamCentArr = new TObjArray(N);  subsamCentArr -> SetOwner(kTRUE); 
      TGraphErrors *grSS[N];
      for (Int_t is = 0; is<N; is++){   
        treeSS->Draw(drawSubSam,fixedCut && Form(subsamCut,x[is],idatatype),"goff");
        grSS[is] = new TGraphErrors(treeSS->GetSelectedRows(),treeSS->GetV2(),treeSS->GetV1());
        grSS[is] -> SetMarkerStyle(20);
        TString sName = branchName;
        subsamCentArr -> SetName(sName.Append("_SS"));
        grSS[is]      -> SetName(sName.Append(Form("_pBin%4.2f",x[is])));      
        subsamCentArr -> AddAt(grSS[is],is);   
      }
    
    // Calculate the statistical errors
      Int_t ngrPoints = (sizeof x)/(sizeof x[0]);
      for (Int_t ierr = 0; ierr<ngrPoints; ierr++) {
        if (!grMoments[objIndex]) continue;
        grMoments[objIndex]->SetPointError(ierr,0,CalculateStatErrors(grSS[ierr]));
      }
    
    // Add graphs to ObjArrays
      momentArr -> AddAt(grMoments[objIndex],objIndex);
      subsamArr -> AddAt(subsamCentArr,objIndex);

      cout << "dataType = " << idatatype << "    moment " << objIndex << " is " << branchName << endl;
      objIndex++;
    
    }
    statResults->GetFile()->cd();
    momentArr->Write(Form("Moments_%d",idatatype),TObject::kSingleKey);
    subsamArr->Write(Form("SubSamples_%d",idatatype),TObject::kSingleKey);
  }
  
  delete statResults;
  fMoments->Close();
  
  MakeRatioPlots(statFileName,0,0.2,3.,-0.8,0.8);
//   gSystem->Exec(Form("rm %s",statFileName.Data()));
 
}
// =====================================================================================================
void ProcessErrorCalculationVSeta(TString momentsFile,Double_t centDown, Double_t centUp, Double_t pDown, Double_t pUp){
                                                                                    
  /*                                        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/test
  tpcdev
  aliroot -l -b
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  ProcessErrorCalculationVSeta("rootFiles/EtaScan_MomentsTree.root",0,5,0.2,1.5)
  */
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C .");
  
  // Output file
  TString statFileName = Form("EtaScan_cent_%d_%d_p_%2.1f_%2.1f.root",Int_t(centDown),Int_t(centUp),pDown,pUp); 
  statResults = new TTreeSRedirector(statFileName,"recreate");

  // Get the trees
  TFile *fMoments = TFile::Open(momentsFile);  
  TTree *treeM    = (TTree*)fMoments -> Get("mcMoments"); 
  TTree *treeSS   = (TTree*)treeM->Clone();
  
  // GetList of branches for each tree
  TObjArray *branchM = (TObjArray*)treeM -> GetListOfBranches();
  Int_t nArrEntries  = branchM->GetEntriesFast();
  cout << " nArrEntries = " << nArrEntries << endl;

//   Output TObjArrays
  TObjArray *momentArr    = new TObjArray(nArrEntries);  momentArr -> SetOwner(kTRUE); 
  TObjArray *subsamArr    = new TObjArray(nArrEntries);  subsamArr -> SetOwner(kTRUE); 

  TGraphErrors *grMoments[nArrEntries]; for (Int_t i=0;i<nArrEntries;i++) grMoments[i]=NULL;
  
  // Create all "moment vs cent" TGraphErrors
  for (Int_t idatatype=0; idatatype<2; idatatype++) {
    Int_t objIndex = 0;
    for (Int_t i=0; i<nArrEntries; ++i){
        
    // Get the brach name 
      TString branchName = branchM->At(i)->GetName();
      TString drawMoment = branchM->At(i)->GetName();
      TString drawSubSam = branchM->At(i)->GetName();
   
      drawMoment.Append(":etaBin");
      drawSubSam.Append(":subsample");
      
      TCut allStat   = Form("subsample==0&&dataType==%d",idatatype);
      TCut fixedCut  = Form("pDown==%f && pUp==%f && (etaUp-etaDown)<0.4 && centDown==%f && centUp==%f",pDown,pUp,centDown,centUp);
      TCut subsamCut = "subsample>0 && subsample<26 && etaBin==%f && dataType==%d";

    // look at only amplitude Mean and Sigma
      if ( (branchName.Contains("cent")  || branchName.Contains("subsample") || branchName.Contains("dataType")) ) continue;
      if ( (branchName.Contains("Up")    || branchName.Contains("Down")      || branchName.Contains("Bin")) ) continue;
      if ( (branchName.Contains("DynDyn")|| branchName.Contains("Stat")      || branchName.Contains("Term")) ) continue;
      if ( (branchName.Contains("pipi")  || branchName.Contains("kaka")      || branchName.Contains("prpr")) ) continue;
      if ( (branchName.Contains("pT")) ) continue;

    // create the graph object
      treeM->Draw(drawMoment,allStat && fixedCut,"goff");
      Int_t N = treeM->GetSelectedRows();
    
    // Fix the ordering of centralitues in graph
      Double_t x[N]; Double_t y[N]; Int_t index[N];
      Double_t *tmpx = treeM->GetV2();  
      Double_t *tmpy = treeM->GetV1(); 
      TMath::Sort(N,tmpx,index);
      for (Int_t k=0; k<N;k++){ x[k] = tmpx[index[k]]; y[k] = tmpy[index[k]]; }

      SetYTitleName(branchName);
      grMoments[objIndex] = new TGraphErrors(N,x,y); 
      grMoments[objIndex]->SetDrawOption("a3");
      grMoments[objIndex]->SetName(branchName);
      grMoments[objIndex]->GetXaxis()->SetTitle("#Delta#eta");
      grMoments[objIndex]->GetYaxis()->SetTitle(yTitleLatex);
      grMoments[objIndex]->SetFillColor(kBlue);
      grMoments[objIndex]->SetMarkerStyle(21);
      grMoments[objIndex]->SetMarkerColor(kBlue);

    
    // Get the subsample
      TObjArray *subsamCentArr = new TObjArray(N);  subsamCentArr -> SetOwner(kTRUE); 
      TGraphErrors *grSS[N];
      for (Int_t is = 0; is<N; is++){   
        treeSS->Draw(drawSubSam,fixedCut && Form(subsamCut,x[is],idatatype),"goff");
        grSS[is] = new TGraphErrors(treeSS->GetSelectedRows(),treeSS->GetV2(),treeSS->GetV1());
        grSS[is] -> SetMarkerStyle(20);
        TString sName = branchName;
        subsamCentArr -> SetName(sName.Append("_SS"));
        grSS[is]      -> SetName(sName.Append(Form("_etaBin%4.2f",x[is])));      
        subsamCentArr -> AddAt(grSS[is],is);   
      }
    
    // Calculate the statistical errors
      Int_t ngrPoints = (sizeof x)/(sizeof x[0]);
      for (Int_t ierr = 0; ierr<ngrPoints; ierr++) {
        if (!grMoments[objIndex]) continue;
        grMoments[objIndex]->SetPointError(ierr,0,CalculateStatErrors(grSS[ierr]));
      }
    
    // Add graphs to ObjArrays
      momentArr -> AddAt(grMoments[objIndex],objIndex);
      subsamArr -> AddAt(subsamCentArr,objIndex);

      cout << "dataType = " << idatatype << "    moment " << objIndex << " is " << branchName << endl;
      objIndex++;
    
    }
    statResults->GetFile()->cd();
    momentArr->Write(Form("Moments_%d",idatatype),TObject::kSingleKey);
    subsamArr->Write(Form("SubSamples_%d",idatatype),TObject::kSingleKey);
  }
  
  delete statResults;
  fMoments->Close();
  
  MakeRatioPlots(statFileName,1,0.2,3.,-0.8,0.8);
//   gSystem->Exec(Form("rm %s",statFileName.Data()));
 
}
// =====================================================================================================
void MakeRatioPlots(TString mfile,Int_t scanType, Double_t pDown, Double_t pUp, Double_t etaDown, Double_t etaUp){
    
  //   
  //  Process the moments file and produce ratio plots with proper error propagation
  //   
   
  gROOT->SetBatch(kTRUE);
  TFile *fmcfile     = TFile::Open(mfile); 
  
  // Get the centrality info from the filename
  TObjArray *objArr2  = mfile.Tokenize("_");
  Int_t centBin     = atoi((objArr2->At(2))->GetName());
  cout << " filename = " << mfile << "      cent = " << centBin << endl;
 
  TString tmpname = mfile;
  tmpname.Prepend("Ratios_");
  ratioresults = new TTreeSRedirector(tmpname,"recreate");
  
  TObjArray * mcarr    = (TObjArray*)fmcfile->Get("Moments_0");
  TObjArray * mcgenarr = (TObjArray*)fmcfile->Get("Moments_1");
  TObjArray * ssarr    = (TObjArray*)fmcfile->Get("SubSamples_0");
  TObjArray * ssgenarr = (TObjArray*)fmcfile->Get("SubSamples_1");
  Int_t nplots = mcarr->GetEntriesFast();
  
  TObjArray * ratioarr = new TObjArray(nplots);  ratioarr -> SetOwner(kTRUE);  
  // loop ober all objects in TObjArray
  for (Int_t iplot=0; iplot<nplots; iplot++){
  
    TGraphErrors *grMC    = (TGraphErrors*)mcarr->At(iplot);
    TGraphErrors *grMCgen = (TGraphErrors*)mcgenarr->At(iplot);
    grMC   ->SetMarkerStyle(20);      grMC   ->SetLineWidth(2);  
    grMC   ->SetMarkerColor(kRed+1);  grMC   ->SetLineColor(kRed+1); 
    grMC   ->SetMarkerSize(1.6);
    grMCgen->SetMarkerStyle(21);      grMCgen->SetLineWidth(2);  
    grMCgen->SetMarkerColor(kBlue+1); grMCgen->SetLineColor(kBlue+1);
    grMCgen->SetMarkerSize(1.4);
    grMCgen->GetYaxis()->SetTitle(grMCgen->GetName());
  
    // Make the ratio plot
    TGraphErrors *grMCratio = RatioPlot(ratioresults,grMC,grMCgen,scanType,pDown,pUp,etaDown,etaUp);
    TString ytitle = grMC->GetName();
    cout << iplot << "    " << ytitle << endl;
    grMCratio->SetName(Form("MCratio_%s",ytitle.Data())); 
    ratioarr->AddAt(grMCratio,iplot);
    
    // Dump ratio plots to look up table
    if (correctionLookup){
      TString ratioGraphName = grMC->GetName();
      if( (ratioGraphName== "pi") || (ratioGraphName =="ka") || (ratioGraphName == "pr") ) {
        correctionLookup->GetFile()->cd();
        grMCratio->GetYaxis()->SetTitle("rec/gen");
        if ( scanType==0 ) {
          grMCratio->GetXaxis()->SetTitle("p (GeV/#it{c})");
          grMCratio->Write(Form("mom_%s_%d",ratioGraphName.Data(),centBin));
        }
        else if (scanType==1) {
          grMCratio->GetXaxis()->SetTitle("#Delta#eta");
          grMCratio->Write(Form("eta_%s_%d",ratioGraphName.Data(),centBin));
        }
      } 
    }
  }
  
  ratioresults->GetFile()->cd();
  ratioarr -> Write("Ratios",TObject::kSingleKey);
  mcarr    -> Write("Moments_rec",TObject::kSingleKey);
  mcgenarr -> Write("Moments_gen",TObject::kSingleKey);
  ssarr    -> Write("SubSamples_rec",TObject::kSingleKey);
  ssgenarr -> Write("SubSamples_gen",TObject::kSingleKey);
  delete ratioresults;
  
} 
// =====================================================================================================
void FastPlotVScent(Int_t style, Int_t col, TString momentsFile, Int_t dataType, TString var, Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp, Int_t same){
 
  /*
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  FastPlotVSeta(20,1,"rootFiles/EtaScan_MomentsTree.root",0,"pi",-0.8,0.8,0.2,2.5,0)
  */
  
  TString treeName = "mcMoments";
  TFile *f = TFile::Open(momentsFile);
  TTree *tree = (TTree*)f->Get(treeName);
  // define cuts
  TString drawStat = Form("%s:centBin",var.Data());
  TString mainCut  = "subsample==0 &&";
  TString extraCut = Form("etaDown==%f&&etaUp==%f&&pDown==%f&&pUp==%f&&dataType==%d",etaDown,etaUp,pDown,pUp,dataType);
  mainCut.Append(extraCut);
  // draw the graph
  if (same) {
    TStatToolkit::MakeGraphErrors(tree,drawStat.Data(),mainCut.Data(),style,col,1.)->Draw("p");
  } else {
    TStatToolkit::MakeGraphErrors(tree,drawStat.Data(),mainCut.Data(),style,col,1.)->Draw("ap");
  }
}
// =====================================================================================================
void FastPlotVSeta(Int_t style, Int_t col, TString momentsFile, Int_t dataType, TString var, Double_t centDown, Double_t pDown, Double_t pUp, Int_t same){
 
  /*
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  FastPlotVSeta(20,1,"rootFiles/EtaScan_MomentsTree.root",0,"pi",0,0.2,2.5,0)
  */
  
  TString treeName = "mcMoments";
  TFile *f = TFile::Open(momentsFile);
  TTree *tree = (TTree*)f->Get(treeName);
  // define cuts
  TString drawStat = Form("%s:etaBin",var.Data());
  TString mainCut  = "subsample==0 && etaBin!=0 &&";
  TString extraCut = Form("centDown==%f && pDown==%f && pUp==%f && dataType==%d",centDown,pDown,pUp,dataType);
  mainCut.Append(extraCut);
  // draw the graph
  if (same) {
    TStatToolkit::MakeGraphErrors(tree,drawStat.Data(),mainCut.Data(),style,col,1.)->Draw("p");
  } else {
    TStatToolkit::MakeGraphErrors(tree,drawStat.Data(),mainCut.Data(),style,col,1.)->Draw("ap");
  }
}
// =====================================================================================================
void FastPlotVSmom(Int_t style, Int_t col, TString momentsFile, Int_t dataType, TString var, Double_t centDown, Int_t same){
  
  /*
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  FastPlotVSmom(20,1,"rootFiles/MomScan_MomentsTree.root",1,"pi",70,0)
  FastPlotVSmom(20,2,"rootFiles/MomScan_MomentsTree.root",0,"pi",70,1)
  */
  
  TString treeName = "mcMoments";
  TFile *f = TFile::Open(momentsFile);
  TTree *tree = (TTree*)f->Get(treeName);
  // define cuts
  TString drawStat = Form("%s:pBin",var.Data());
  TString mainCut  = "(pUp-pDown)<0.5 && subsample==0 && etaBin==0 &&";
  TString extraCut = Form("centDown==%f && dataType==%d",centDown,dataType);
  mainCut.Append(extraCut);
  // draw the graph
  if (same) {
    TStatToolkit::MakeGraphErrors(tree,drawStat.Data(),mainCut.Data(),style,col,1.)->Draw("p");
  } else {
    TStatToolkit::MakeGraphErrors(tree,drawStat.Data(),mainCut.Data(),style,col,1.)->Draw("ap");
  }
}
// =====================================================================================================
void RatioScan(TString etaScanFile,TString momScanFile ){

  // To RUN
  /*                                        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/MC_genVSrec/SSanalysis/EtaMomentumScan/cRows80/centBinWidth_10_all/Results
  tpcdev
  aliroot -l -b
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  RatioScan("rootFiles/EtaScan_MomentsTree.root","rootFiles/MomScan_MomentsTree.root")
  */
  //
  
  correctionLookup = new TTreeSRedirector("CorrectionLookup.root","recreate");
  
//   Int_t nMomentumBins = (sizeof pDownArray)/(sizeof pDownArray[0]);
  Int_t nEtaBins  = (sizeof etaDownArray)/(sizeof etaDownArray[0]);
  Int_t nCentBins = (sizeof centDownArray)/(sizeof centDownArray[0]);

  cout << " ================= Cent Scan for EffvsMOM Lookuptable ================= " << endl;
  for (Int_t icent=0; icent<nCentBins; icent++)
  {
    cout << "    etaRange    = " << -0.8 << " - " << 0.8 ;
    cout << "    centRange   = " << centDownArray[icent]     << " - " << centUpArray[icent]     << endl;
    ProcessErrorCalculationVSmomentum(momScanFile,centDownArray[icent],centUpArray[icent]);
  }
  
  cout << " ================= Cent Scan for EffvsETA Lookuptable ================= " << endl;
  for (Int_t icent=0; icent<nCentBins; icent++)
  {
    cout << "    p Range     = " << 0.2 << " - " << 3. ;
    cout << "    centRange   = " << centDownArray[icent]     << " - " << centUpArray[icent]     << endl;
    ProcessErrorCalculationVSeta(etaScanFile,centDownArray[icent],centUpArray[icent],0.2,2.);
  }

  cout << " ================= All ratios ================= " << endl;
  for (Int_t ip=0; ip<10; ip++)
    for (Int_t ieta=0; ieta<nEtaBins; ieta++)
  {
    cout << "    etaRange = " << etaDownArray[ieta] << " - " << etaUpArray[ieta] ;
    cout << "    pRange   = " << pDownArray[ip]     << " - " << pUpArray[ip]     << endl;
    ProcessErrorCalculationVScent(etaScanFile,etaDownArray[ieta],etaUpArray[ieta],pDownArray[ip],pUpArray[ip]);
  }
  
}
// =====================================================================================================
Double_t CalculateStatErrors(TGraphErrors *grS){

  //
  // calculate staistical error of each point of the moments
  //

  Double_t meanMS = TMath::Abs(grS->GetMean(2));
  Double_t n      = grS->GetN();
  Double_t sum = 0;
    
  for (Int_t i = 0; i < n; i++) {
    Double_t y= TMath::Abs(grS->GetY()[i]);
    sum = sum+(y-meanMS)*(y-meanMS);
  }
  
  Double_t err = (n!=0) ? TMath::Sqrt(sum/(n*(n-1))) : 0.;
  
  if (err) return err;
  else return 0.;
  
}
// =====================================================================================================
void AnalyseWs(TString wFile, TString wSliceFile){

  //
  // Calculate statistical errors using subsample method
  /*
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/IdDataTrees_300000_Results_OK2
  aliroot -l
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  AnalyseWs("streamer_All.root","streamer_Slice.root")
  */
  //
  
  TTreeSRedirector *wstream = new TTreeSRedirector("PlotsWs.root","recreate");
  
  // Get the tree
  TFile *fWs     = TFile::Open(wFile); 
  TFile *fSliceWs = TFile::Open(wSliceFile);       

  TTree *tree       = (TTree*)fWs -> Get("debugW");
  TTree *treeSlice  = (TTree*)fSliceWs -> Get("debugW");
  TTree *treeSliceW = (TTree*)fSliceWs -> Get("debugw");
  
  TCanvas *wCan     = TwoScales(treeSliceW);
  TCanvas *wOverlay = OverlayPads(treeSliceW);


  Double_t nTreeEntries = (lowStat) ? 10000000. : tree->GetEntries();

  // Debug histos
  TH2D * h2dEdx = new TH2D("h2dEdx","h2dEdx",300,0,300,2000,0,1000);
  TH1D * h1dEdx = new TH1D("h1dEdx","h1dEdx",2000,0,1000);
  TH1D * hWel  = new TH1D("hWel"   ,"hWel"  ,4000,0,2000);
  TH1D * hWpi  = new TH1D("hWpi"   ,"hWpi"  ,4000,0,2000);
  TH1D * hWka  = new TH1D("hWka"   ,"hWka"  ,4000,0,2000);
  TH1D * hWpr  = new TH1D("hWpr"   ,"hWpr"  ,4000,0,2000);
  
  TH2D * h2SlicedEdx = new TH2D("h2SlicedEdx","h2SlicedEdx",300,0,300,2000,0,1000);
  TH1D * h1SlicedEdx = new TH1D("h1SlicedEdx","h1SlicedEdx",3000,0,1000);
  TH1D * hSliceWel  = new TH1D("hSliceWel"   ,"hSliceWel"  ,3000,0,1000);
  TH1D * hSliceWpi  = new TH1D("hSliceWpi"   ,"hSliceWpi"  ,3000,0,1000);
  TH1D * hSliceWka  = new TH1D("hSliceWka"   ,"hSliceWka"  ,3000,0,1000);
  TH1D * hSliceWpr  = new TH1D("hSliceWpr"   ,"hSliceWpr"  ,3000,0,1000);

  // Projections
  tree ->Project("h2dEdx"  ,"dEdx:p" ,""        ,"",nTreeEntries);
  tree ->Project("h1dEdx"  ,"dEdx"   ,""        ,"",nTreeEntries);
  tree ->Project("hWel"    ,"W"      ,"pType==0","",nTreeEntries);
  tree ->Project("hWpi"    ,"W"      ,"pType==1","",nTreeEntries);
  tree ->Project("hWka"    ,"W"      ,"pType==3","",nTreeEntries);
  tree ->Project("hWpr"    ,"W"      ,"pType==2","",nTreeEntries);
  
  treeSlice ->Project("h2SlicedEdx"  ,"dEdx:p" ,""        ,"",nTreeEntries);
  treeSlice ->Project("h1SlicedEdx"  ,"dEdx"   ,""        ,"",nTreeEntries);
  treeSlice ->Project("hSliceWel"    ,"W"      ,"pType==0","",nTreeEntries);
  treeSlice ->Project("hSliceWpi"    ,"W"      ,"pType==1","",nTreeEntries);
  treeSlice ->Project("hSliceWka"    ,"W"      ,"pType==3","",nTreeEntries);
  treeSlice ->Project("hSliceWpr"    ,"W"      ,"pType==2","",nTreeEntries);

  wstream->GetFile()->cd();
  wCan     ->Write("w-dEdx");
  wOverlay ->Write("wOverlay-dEdx");

  // one slice at p bin = 19
  hSliceWel->Write("hWel");
  hSliceWpi->Write("hWpi");
  hSliceWka->Write("hWka");
  hSliceWpr->Write("hWpr");
  h1SlicedEdx->Write("h1dEdx");
  h2SlicedEdx->Write("h2dEdx");
  // all eta and p bins
  hWel->Write("hWel_All");
  hWpi->Write("hWpi_All");
  hWka->Write("hWka_All");
  hWpr->Write("hWpr_All");
  h1dEdx->Write("h1dEdx_All");
  h2dEdx->Write("h2dEdx_All");
    
  delete wstream;
}
// =====================================================================================================
void AnalyseHistsInFile(TString file){
  // 
  // delete TH1D object from file
  /*
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/IdDataTrees_MC_060514/
  aliroot -l 
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
  AnalyseHistsInFile("TIdenResults_0_0_GAUS.root"); 2> /dev/null
  AnalyseHistsInFile("TIdenResults_0_0_GENGAUS.root"); 2> /dev/null
  AnalyseHistsInFile("TIdenResults_0_0_GENGAUS_auto.root"); 2> /dev/null
  AnalyseHistsInFile("DataTreeMC_50000_0_5.root"); 2> /dev/null
  */
  //   
  
  TFile * f = TFile::Open(file);   
  TH1D *h[100];
  TH1D *hNoRebin[100];
  TString name[100];
  TString nameNoRebin[100];
  
  
  TIter next(f->GetListOfKeys());
  TKey *key;
  
  Int_t i=0;
  while ((key = (TKey*)next())) {
    
    TString histname(key->GetName());
     
    // loop over only on th1 objects
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if ( !cl->InheritsFrom("TH1") )       continue;
    if ( !histname.Contains("histo") )    continue;
    if ( histname.Contains("histoBin") )  continue;
    if ( histname.Contains("histodEdx") ) continue;
    if ( histname.Contains("histoDe") )   continue;
    if ( histname.Contains("histoMu") )   continue;
    if ( histname.Contains("histoAll") )  continue;

    
    hNoRebin[i] = (TH1D*)key->ReadObj();
    h[i] = (TH1D*)key->ReadObj();
    cout << histname << "  binWidth is " << h[i]->GetXaxis()->GetBinWidth(50) << endl;
    if ( !file.Contains("stream") ) h[i]->Rebin(3);
    name[i]=histname;
    nameNoRebin[i]=histname;
    nameNoRebin[i].Append("_noRebin");
    
    // Fit the histogram
    FitHistogram(h[i]);    
    i++;
  }
  
  TString outputfile;
  if (file.Contains("_GAUS."))         outputfile = "Whists_GAUS.root";
  if (file.Contains("_GENGAUS."))      outputfile = "Whists_GENGAUS.root";
  if (file.Contains("_GENGAUS_auto.")) outputfile = "Whists_GENGAUS_auto.root";
  if (file.Contains("DataTreeMC"))     outputfile = "Whists_RealMC.root";
  if (file.Contains("stream"))         outputfile = "Whists_From_Streamer.root";

  TTreeSRedirector * pcstream = new TTreeSRedirector(outputfile,"recreate");
  
  pcstream->GetFile()->cd();
  for(Int_t j=0; j<i; j++) hNoRebin[j] -> Write(nameNoRebin[j]);
  for(Int_t j=0; j<i; j++) h[j]        -> Write(name[j]);
    
  
  
  delete pcstream;
  f->Close();
}
// =====================================================================================================
void FitHistogram(TH1D *hClean){
 
  //
  // Fit histogram with gengaus
  //

//   gStyle->SetOptStat(220002210); // --> (Kurtosis,Skewness,integral,overflows,underflows,Rms,Mean,entries,name) --> (KSiRMe)
  gStyle->SetOptFit(1111);
  TString hname = hClean->GetName();
  hClean->SetLineColor(kBlack);
  
  TVirtualFitter::SetMaxIterations(1000000);
  
    // helper numbers for setparams
  Double_t maxBin          = hClean->GetBinContent(hClean->GetMaximumBin());  if ((maxBin<5.)) return;   
  Double_t mean            = hClean->GetMean();
  Double_t rms             = hClean->GetRMS();
  Double_t sampleFitWindow = 4.; 
        
    // fit function
  TString fitFunctionGenGaus = "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
  TF1 *asyGaus = new TF1("asyGaus",fitFunctionGenGaus,mean-sampleFitWindow*rms,mean+4.2*rms); 
  asyGaus->SetParNames("fAmplitude","fMean","fSigma","fKurtosis","fSkewness"); 
  asyGaus->SetParameters(maxBin,mean,rms,2.,1.); 
  asyGaus->SetParLimits(0,maxBin/2.   ,maxBin);
  asyGaus->SetParLimits(1,mean-2.*rms ,mean+2.*rms);
  asyGaus->SetParLimits(2,rms         ,3.*rms);
  asyGaus->SetParLimits(3,1.5         ,3.);
  asyGaus->SetParLimits(4,0.0         ,3.);
  asyGaus->SetLineWidth(3);
  asyGaus->SetLineColor(kBlack);
  
  if ( !hname.Contains("W") ) {
    hClean->SetLineColor(kRed);
    asyGaus->SetLineColor(kRed);
  }
 
  hClean->GetXaxis()->SetRangeUser(TMath::Max(0.,mean-5.*rms),mean+5.*rms);
    
  // retrieve fit parameters
  hClean->Fit(asyGaus,"QMR");  
  delete asyGaus;
}
// =====================================================================================================
void GraphShade() {
  TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);

  c1->SetGrid();
  c1->DrawFrame(0,0,2.2,12);
   
  const Int_t n = 20;
  Double_t x[n], y[n],ymin[n], ymax[n];
  Int_t i;
  for (i=0;i<n;i++) {
    x[i] = 0.1+i*0.1;
    ymax[i] = 10*sin(x[i]+0.2);
    ymin[i] = 8*sin(x[i]+0.1);
    y[i] = 9*sin(x[i]+0.15);
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
}
// =====================================================================================================
void ErrorBandTH1(){
  TH1F *h = new TH1F("h","",10,0,10);

  h->SetLineColor(kRed);

  for(int i=0; i<9; i++){
    h->Fill(i,i*0.5);
  }

  h->DrawCopy("hist"); 
  h->SetFillColor(kBlue);
  h->SetFillStyle(3018);
  h->Draw("e2same");
 
}
// =====================================================================================================
TCanvas *TwoScales(TTree *t){
  //
  //example of macro illustrating how to superimpose two histograms
  //with different scales in the "same" pad.
  //  
  
  gStyle->SetOptStat(kFALSE);
  Double_t nEntries = (lowStat) ? 10000000. : t->GetEntries() ;
  
  // dEdx histo for the right TGaxis
  TH1D * hdEdx = new TH1D("hdEdx","hdEdx",1000,0,400);
  t->Project("hdEdx","dEdx","eta==1","",nEntries);
  Double_t maxBin = hdEdx->GetBinContent(hdEdx->GetMaximumBin()); 
  hdEdx->Scale(1./maxBin);
  hdEdx->SetLineWidth(2);
  hdEdx->SetLineColor(kMagenta);

  // w distributions of each particle
  TH2D * hwEl = new TH2D("hwEl","hwEl",1000,0,400,100,0,1.2);
  TH2D * hwPi = new TH2D("hwPi","hwPi",1000,0,400,100,0,1.2);
  TH2D * hwKa = new TH2D("hwKa","hwKa",1000,0,400,100,0,1.2);
  TH2D * hwPr = new TH2D("hwPr","hwPr",1000,0,400,100,0,1.2);

  t->Project("hwEl","w:dEdx","pType==0 && eta==1","",nEntries);
  t->Project("hwPi","w:dEdx","pType==1 && eta==1","",nEntries);
  t->Project("hwKa","w:dEdx","pType==3 && eta==1","",nEntries);
  t->Project("hwPr","w:dEdx","pType==2 && eta==1","",nEntries);

      
  hwEl->SetMarkerStyle(20);
  hwPi->SetMarkerStyle(20);
  hwKa->SetMarkerStyle(20);
  hwPr->SetMarkerStyle(20);

  hwEl->SetMarkerColor(kBlack);
  hwPi->SetMarkerColor(kGreen);
  hwKa->SetMarkerColor(kRed);
  hwPr->SetMarkerColor(kBlue);


  TCanvas *c2 = new TCanvas("c2","hists with different scales",600,400);
  c2->SetGrid();
  
//   hdEdx->GetYaxis()->SetLabelColor(kMagenta);
  hdEdx->Draw(); 
  hdEdx->GetXaxis()->SetTitle("TPC dEdx Signal");
  hdEdx->GetYaxis()->SetTitle("identity variable #omega (%)");
//   hdEdx->GetYaxis()->SetTitleColor(kMagenta);

  c2->Update();  

//   Float_t rightmax = 1.1*hwPi->GetMaximum();
  Float_t rightmax = 1.05*1.;
  Float_t scale    = gPad->GetUymax()/rightmax;
  cout << rightmax << "   " << scale << "   " << gPad->GetUymax() <<  endl;
  
  hwEl->Scale(scale);
  hwPi->Scale(scale);
  hwKa->Scale(scale);
  hwPr->Scale(scale);
  
  hwEl->Draw("same"); 
  hwPi->Draw("same"); 
  hwKa->Draw("same"); 
  hwPr->Draw("same"); 
   
   //draw an axis on the right side
//   TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(),0,rightmax,510,"+L");
//   axis->SetLineColor(kBlack);
//   axis->SetLabelColor(kBlack);
//   axis->SetTitle("identity variable w [probability]");
//   axis->Draw();
  
  TLegend *leg = new TLegend(0.7, 0.6, 0.9, 0.9);
  leg->SetHeader("#omega : identity variable");
  leg->SetFillColor(0);
  leg->AddEntry(hwEl,"Electron","p");
  leg->AddEntry(hwPi,"Pion"    ,"p");
  leg->AddEntry(hwKa,"Kaon"    ,"p");
  leg->AddEntry(hwPr,"Proton"  ,"p");
  leg->Draw();
  
  return c2;
  
}
// =====================================================================================================
TCanvas *OverlayPads(TTree *t){
  //
  //example of macro illustrating how to superimpose two histograms
  //with different scales in the "same" pad.
  //  
  
  gStyle->SetOptStat(kFALSE);
  Double_t nEntries = (lowStat) ? 1000000. : t->GetEntries() ;
  
  // -----------------------------------------------------
  gROOT->Reset();
  TCanvas *c1 = new TCanvas("c1","gerrors2",200,10,700,500);
  TPad *pad = new TPad("pad","",0,0,1,1);
  pad->SetFillColor(0);
  pad->SetGrid();
  pad->Draw();
  pad->cd();

      // draw a frame to define the range
  TH1F *hr = c1->DrawFrame(0.,0,400.,1.2);
  hr->SetXTitle("TPC dEdx Signal");
  hr->SetYTitle("w");
  pad->GetFrame()->SetFillColor(0);
  pad->GetFrame()->SetBorderSize(12);

  // w distributions of each particle
  t->Draw("w:dEdx","pType==0 && eta==1","goff",nEntries);
  TGraphErrors *hwEl = new TGraphErrors(t->GetSelectedRows(),t->GetV2(),t->GetV1());
  hwEl->SetName("hwEl");
  hwEl->SetMarkerStyle(7);
  hwEl->SetMarkerColor(kBlack);

  t->Draw("w:dEdx","pType==1 && eta==1","goff",nEntries);
  TGraphErrors *hwPi = new TGraphErrors(t->GetSelectedRows(),t->GetV2(),t->GetV1());
  hwPi->SetName("hwPi");
  hwPi->SetMarkerStyle(7);
  hwPi->SetMarkerColor(kGreen);

  t->Draw("w:dEdx","pType==3 && eta==1","goff",nEntries);
  TGraphErrors *hwKa = new TGraphErrors(t->GetSelectedRows(),t->GetV2(),t->GetV1());
  hwKa->SetName("hwKa");
  hwKa->SetMarkerStyle(7);
  hwKa->SetMarkerColor(kRed);

  t->Draw("w:dEdx","pType==2 && eta==1","goff",nEntries);
  TGraphErrors *hwPr = new TGraphErrors(t->GetSelectedRows(),t->GetV2(),t->GetV1());
  hwPr->SetName("hwPr");
  hwPr->SetMarkerStyle(7);
  hwPr->SetMarkerColor(kBlue);
 
  hwEl->Draw("p"); 
  hwPi->Draw("p"); 
  hwKa->Draw("p"); 
  hwPr->Draw("p"); 
  
  // Create second pad
  c1->cd();
  TPad *overlay = new TPad("overlay","",0,0,1,1);
  overlay->SetFillStyle(4000);
  overlay->SetFillColor(0);
  overlay->SetFrameFillStyle(4000);
  overlay->Draw();
  overlay->cd();
  
  // dEdx histo for the right TGaxis
  TH1D * hdEdx = new TH1D("hdEdx","hdEdx",1000,0,400);
  t->Project("hdEdx","dEdx","eta==1","",nEntries);
  Double_t maxBin = hdEdx->GetBinContent(hdEdx->GetMaximumBin()); 
  hdEdx->Scale(1./maxBin);
  hdEdx->SetLineWidth(2);
  hdEdx->SetLineColor(kMagenta);
  hdEdx->GetYaxis()->SetRangeUser(0.,1.2);
  
  
  Double_t xmin = pad->GetUxmin();
  Double_t ymin = 0;
  Double_t xmax = pad->GetUxmax();
  Double_t ymax = 1.2;
  TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
  hframe->GetXaxis()->SetLabelOffset(99);
  hframe->GetYaxis()->SetLabelOffset(99);
  
  hdEdx->Draw();
  
  //Draw an axis on the right side
  TGaxis *axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
  axis->SetLineColor(kMagenta);
  axis->SetLabelColor(kMagenta);
  axis->Draw();
  
  // Draw the TLegend
  TLegend *leg = new TLegend(0.7, 0.6, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->AddEntry(hwEl,"Electron","p");
  leg->AddEntry(hwPi,"Pion"    ,"p");
  leg->AddEntry(hwKa,"Kaon"    ,"p");
  leg->AddEntry(hwPr,"Proton"  ,"p");
  leg->Draw();
  
  return c1;
  
}
// =====================================================================================================
TGraphErrors * RatioPlot(TTreeSRedirector *ratiofile, TGraphErrors *g1, TGraphErrors *g2, Int_t scanType, Double_t pDown, Double_t pUp, Double_t etaDown, Double_t etaUp) {
  
  const Int_t N = g1->GetN();
  Double_t x[N];
  Double_t y[N];
  Double_t errx[N];
  Double_t erry[N];
  
  for (Int_t i=0; i<N; i++){
    if ((g2->GetY()[i]!=0.)) {
      y[i]=g1->GetY()[i]/g2->GetY()[i];
    } else {
      y[i]=0.;
    }
    x[i]=g1->GetX()[i];
    errx[i]=0.;
    erry[i]=0.;
  }
  
  TLegend *leg = new TLegend(0.2,0.7,0.6,0.85); 
  leg->SetTextFont(62);  leg->SetTextSize(0.03); leg->SetFillColor(0);  leg->SetBorderSize(0);
  TLine* line = new TLine(0., 1.,80., 1.);
  line -> SetLineColor(1); line -> SetLineWidth(4); line -> SetLineStyle(2);
       
  TGraphErrors *grRatio = new TGraphErrors(N,x,y,errx,erry);  
  grRatio->SetMinimum(-100);
  
    // Calculate y axis ranges 
  Double_t max1 = TMath::MaxElement(g1->GetN(),g1->GetY()); 
  Double_t max2 = TMath::MaxElement(g2->GetN(),g2->GetY()); 
  Double_t maxx = TMath::Max(max1,max2)*1.1;
  Double_t min1 = TMath::MinElement(g1->GetN(),g1->GetY()); 
  Double_t min2 = TMath::MinElement(g2->GetN(),g2->GetY()); 
  Double_t minn = TMath::Min(min1,min2)*0.8;
  g1->GetYaxis()->SetRangeUser(minn,maxx);
  g1->GetYaxis()->SetTitleOffset(1.25);
  leg->AddEntry(g1,"MC reconstructed","lep");
  leg->AddEntry(g2,"MC generated","lep");
  TLegendEntry *header;
  if (scanType>1){
    //     leg->SetHeader(Form("%2.1f<#eta<%2.1f, %2.1f<#it{p}<%2.1f GeV/#it{c}",etaDown,etaUp,pDown,pUp));
    leg->SetHeader(Form("|#eta|<%2.1f, %2.1f<#it{p}<%2.1f GeV/#it{c}",etaUp,pDown,pUp));
    header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    header->SetTextSize(0.045);
  }
   
  // plot ratio and graphs together in one canvas 
  if (!gROOT->GetListOfCanvases()->FindObject("cRatios")) delete gROOT->GetListOfCanvases()->FindObject("cRatios");
  TCanvas *cRatios = new TCanvas("cRatios","cRatios",600,700);  cRatios->SetBottomMargin(20.); cRatios->SetLeftMargin(20.);
  cRatios->Clear(); cRatios->Update(); 
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);  pad1->SetBottomMargin(0); pad1->SetLeftMargin(0.14); pad1->SetTicks(1,1);
  pad1->Draw(); 
  pad1->cd();
  g1->Draw("ap");
  g2->Draw("p");
  leg->Draw("same");
  
  cRatios->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
  pad2->SetTopMargin(0.); pad2->SetLeftMargin(0.14); pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();
  pad2->SetTicks(1,1);
  grRatio->SetMarkerStyle(20); 
  grRatio->SetMarkerSize(1.3); 
  grRatio->SetMarkerColor(kBlack);
  grRatio->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle()); grRatio->GetXaxis()->SetTitleOffset(1.);
  grRatio->GetYaxis()->SetTitle("rec/gen");
  grRatio->GetYaxis()->SetTitleSize(0.13);
  grRatio->GetYaxis()->SetRangeUser(0.3,1.7);
  
  // Calculate y axis ranges 
//   Double_t max = TMath::MaxElement(grRatio->GetN(),grRatio->GetY())*1.5; 
//   Double_t min = TMath::MinElement(grRatio->GetN(),grRatio->GetY())*0.5; 
//   grRatio->GetYaxis()->SetRangeUser(min,max);
//   grRatio->GetYaxis()->SetRangeUser(0.,2.);
  
  // Do the error propagation
  for (Int_t ierrpro=0; ierrpro<g1->GetN(); ierrpro++){
    Double_t xx    = g1   ->GetY()[ierrpro];
    Double_t errxx = g1   ->GetErrorY(ierrpro);
    Double_t yy    = g2->GetY()[ierrpro];
    Double_t erryy = g2->GetErrorY(ierrpro);
    Double_t errRatio=0;
    if (yy!=0) errRatio = TMath::Sqrt((errxx*errxx)/(yy*yy)+(xx*xx*erryy*erryy)/(yy*yy*yy*yy));
    grRatio->SetPointError(ierrpro,0,errRatio); 
  }
  
  TGraphErrors *gr = (TGraphErrors*)grRatio->Clone(); // graph to return
  grRatio->GetYaxis()->SetLabelSize(0.1);
  grRatio->GetXaxis()->SetLabelSize(0.1);
  grRatio->GetXaxis()->SetTitleSize(0.1);
  grRatio->GetYaxis()->SetTitleSize(0.1);
  grRatio->GetYaxis()->SetTitleOffset(0.5);
  
  grRatio->Draw("ap");
  line->Draw("same");
  
  TString canvasname = g1->GetName();
  cRatios->SaveAs(Form("canvas_%s.pdf",canvasname.Data()));
  ratiofile->GetFile()->cd();
  cRatios->Write(Form("canvas_%s",canvasname.Data()));
  delete cRatios;
  return gr;
}
// =====================================================================================================
void SetYTitleName(TString branchName){
  
  //
  // Set the title name 
  //
  
  if (branchName=="pikaNuDynDyn") yTitleLatex="#nu [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="piprNuDynDyn") yTitleLatex="#nu [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="kaprNuDynDyn") yTitleLatex="#nu [p+#bar{p},K^{+}+K^{-}]";
  if (branchName=="pikaNuStat") yTitleLatex="#nu_{stat} [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="piprNuStat") yTitleLatex="#nu_{stat} [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="kaprNuStat") yTitleLatex="#nu_{stat} [p+#bar{p},K^{+}+K^{-}]";
  if (branchName=="pikaNuDyn") yTitleLatex="#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="piprNuDyn") yTitleLatex="#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="kaprNuDyn") yTitleLatex="#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]";
  if (branchName=="pipiNuDyn") yTitleLatex="#nu_{dyn} [#pi^{+}+#pi^{-},#pi^{+}+#pi^{-}]";
  if (branchName=="kakaNuDyn") yTitleLatex="#nu_{dyn} [K^{+}+K^{-},K^{+}+K^{-}]";
  if (branchName=="prprNuDyn") yTitleLatex="#nu_{dyn} [p+#bar{p},p+#bar{p}]";
  if (branchName=="pika") yTitleLatex="#LT#pi^{+}+#pi^{-},K^{+}+K^{-}#GT";
  if (branchName=="pipr") yTitleLatex="#LT#pi^{+}+#pi^{-},p+#bar{p}#GT";
  if (branchName=="kapr") yTitleLatex="#LTK^{+}+K^{-},p+#bar{p}#GT";
  if (branchName=="pipi") yTitleLatex="#LT#pi^{+}+#pi^{-},#pi^{+}+#pi^{-}#GT";
  if (branchName=="kaka") yTitleLatex="#LTK^{+}+K^{-},K^{+}+K^{-}#GT";
  if (branchName=="prpr") yTitleLatex="#LTp+#bar{p},p+#bar{p}#GT";
  if (branchName=="pi") yTitleLatex="#LT#pi^{+}+#pi^{-}#GT";
  if (branchName=="ka") yTitleLatex="#LTK^{+}+K^{-}#GT";
  if (branchName=="pr") yTitleLatex="#LTp+#bar{p}#GT";
  
  if (branchName=="deltapika") yTitleLatex="#Delta [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="deltapipr") yTitleLatex="#Delta [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="deltakapr") yTitleLatex="#Delta [p+#bar{p},K^{+}+K^{-}]";
  
  if (branchName=="sigmapika") yTitleLatex="#Sigma [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="sigmapipr") yTitleLatex="#Sigma [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="sigmakapr") yTitleLatex="#Sigma [p+#bar{p},K^{+}+K^{-}]";
 
  if (branchName=="phipika") yTitleLatex="#Phi [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="phipipr") yTitleLatex="#Phi [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="phikapr") yTitleLatex="#Phi [p+#bar{p},K^{+}+K^{-}]";
  
  if (branchName=="wpi") yTitleLatex="#omega_{#pi}";
  if (branchName=="wka") yTitleLatex="#omega_{K}";
  if (branchName=="wpr") yTitleLatex="#omega_{p}";
  
  if (branchName=="nudynpika") yTitleLatex="#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="nudynpipr") yTitleLatex="#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="nudynkapr") yTitleLatex="#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]";
 
  if (branchName=="cDeltapika") yTitleLatex="C_{#Delta} [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="cDeltapipr") yTitleLatex="C_{#Delta} [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="cDeltakapr") yTitleLatex="C_{#Delta} [p+#bar{p},K^{+}+K^{-}]";
 
  if (branchName=="cSigmapika") yTitleLatex="C_{#Sigma} [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="cSigmapipr") yTitleLatex="C_{#Sigma} [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="cSigmakapr") yTitleLatex="C_{#Sigma} [p+#bar{p},K^{+}+K^{-}]";
           
}
// =====================================================================================================
void MultiplydNch(TString momFile, Int_t flag){

  //
  // correct the centrality axis of the moments
  //
  /*
   
   cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_HIJING_TightCuts/test
   aliroot -l 
   .L ~/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C+
   MultiplydNch("Ratios_CentScan_eta_-0.8-0.8_mom_0.2-1.5.root", 0)
   MultiplydNch("Ratios_CentScan_eta_-0.8-0.8_mom_0.2-1.5.root", 1)
   MultiplydNch("Ratios_CentScan_eta_-0.8-0.8_mom_0.2-1.5.root", 2)
   MultiplydNch("Ratios_CentScan_eta_-0.8-0.8_mom_0.2-1.5.root", 3)
   MultiplydNch("Ratios_CentScan_eta_-0.8-0.8_mom_0.2-1.5.root", 4)
   MultiplydNch("Ratios_CentScan_eta_-0.8-0.8_mom_0.2-1.5.root", 5)
   MultiplydNch("Ratios_CentScan_eta_-0.8-0.8_mom_0.2-1.5.root", 6)

  
  */
  
  TTreeSRedirector *scaling = new TTreeSRedirector(Form("Nudyn_%d_Scaling.root",flag),"recreate");
  
  TFile *f = TFile::Open(momFile);
  TObjArray *arr = (TObjArray*)f->Get("Moments_gen");
  TGraphErrors *grdn = (TGraphErrors*)arr->FindObject("pi");
  TGraphErrors *grpi = (TGraphErrors*)arr->FindObject("pi");
  TGraphErrors *grka = (TGraphErrors*)arr->FindObject("ka");
  TGraphErrors *grpr = (TGraphErrors*)arr->FindObject("pr");

  TString nuName;
  for (Int_t igr=0; igr<3;igr++){
    
    if(igr==0) nuName = "pikaNuDyn";
    if(igr==1) nuName = "piprNuDyn";
    if(igr==2) nuName = "kaprNuDyn";

    TGraphErrors *grNu = (TGraphErrors*)arr->FindObject(nuName);
    Int_t NN = grNu->GetN(); 
    Double_t yErr[NN];
    Double_t y[NN];
    Double_t x[NN];
    Double_t errY;
    Double_t errYdnch;
    
    // Prepare data points and do error propagation
    for (Int_t i=0;i<NN;i++){
      
      // scaling wrt to number of pions 
      if (flag==0) {
	y[i] = grNu->GetY()[i]*grdn->GetY()[i];
	x[i] = grNu->GetX()[i];
	errY     = grNu->GetEY()[i];
	errYdnch = grdn->GetEY()[i];
	cout << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << grdn->GetY()[i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(grdn->GetY()[i]*grdn->GetY()[i])*(errY*errY) );
      }
      
      // scaling wrt to Npart
      if (flag==1) {
	y[i] = grNu->GetY()[i]*nPartArr[NN-1-i];
	x[i] = grNu->GetX()[i];
	errY     = grNu->GetEY()[i];
	errYdnch = nPartArrErr[NN-1-i];
	cout << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << nPartArr[NN-1-i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(nPartArr[NN-1-i]*nPartArr[NN-1-i])*(errY*errY) );
      }
      
      // scaling wrt to dnch HIJING
      if (flag==2) {
	y[i] = grNu->GetY()[i]*dnchArrayHIJING[NN-1-i];
	x[i] = grNu->GetX()[i];
	errY     = grNu->GetEY()[i];
	errYdnch = dnchArrayHIJINGErr[NN-1-i];
	cout << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << dnchArrayHIJING[NN-1-i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(dnchArrayHIJING[NN-1-i]*dnchArrayHIJING[NN-1-i])*(errY*errY) );
      }
      
      // scaling wrt to dnch AMPT 
      if (flag==3) {
	y[i] = grNu->GetY()[i]*dnchArrayAMPT[NN-1-i];
	x[i] = grNu->GetX()[i];
	errY     = grNu->GetEY()[i];
	errYdnch = dnchArrayAMPTErr[NN-1-i];
	cout << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << dnchArrayAMPT[NN-1-i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(dnchArrayAMPT[NN-1-i]*dnchArrayAMPT[NN-1-i])*(errY*errY) );
      }
     
      // scaling wrt to anar
      if (flag==4) {
	if (igr==0) y[i] = grNu->GetY()[i]/(1./grpi->GetY()[i]+1./grka->GetY()[i]);
	if (igr==1) y[i] = grNu->GetY()[i]/(1./grpi->GetY()[i]+1./grpr->GetY()[i]);
	if (igr==2) y[i] = grNu->GetY()[i]/(1./grpr->GetY()[i]+1./grka->GetY()[i]);
	x[i]     = grNu->GetX()[i];
	cout << i << "  " << grNu->GetY()[i] << "  " << 1./grpi->GetY()[i]+1./grka->GetY()[i] << "  " << y[i] <<  "  " << grpi->GetY()[i] << "  " << grka->GetY()[i] << "  " << grpr->GetY()[i] << endl;
	yErr[i] = 0.;
      }
      
      // scaling wrt to dnch HIJING
      if (flag==5) {
	y[i] = grNu->GetY()[i]*dnchArrayHIJING[NN-1-i];
	x[i] = dnchArrayHIJING[NN-1-i];
	errY     = grNu->GetEY()[i];
	errYdnch = dnchArrayHIJINGErr[NN-1-i];
	cout << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << dnchArrayHIJING[NN-1-i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(dnchArrayHIJING[NN-1-i]*dnchArrayHIJING[NN-1-i])*(errY*errY) );
      }
      
      // scaling wrt to dnch AMPT 
      if (flag==6) {
	y[i] = grNu->GetY()[i]*dnchArrayAMPT[NN-1-i];
	x[i] = dnchArrayAMPT[NN-1-i];
	errY     = grNu->GetEY()[i];
	errYdnch = dnchArrayAMPTErr[NN-1-i];
	cout << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << dnchArrayAMPT[NN-1-i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(dnchArrayAMPT[NN-1-i]*dnchArrayAMPT[NN-1-i])*(errY*errY) );
      }
    
   
     
    }
    
    // Final plots 
    TGraphErrors *gr = new TGraphErrors(NN,x,y,grdn->GetEX(),yErr);
    gr->SetDrawOption("a3");
    gr->SetMarkerSize(2);
    gr->SetName(grNu->GetName());
    gr->GetXaxis()->SetTitle("centrality (%)");
    if (flag==0 || flag==2 || flag==3) {
      if (nuName=="pikaNuDyn") gr->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]xdN_{ch}/d#eta");
      if (nuName=="piprNuDyn") gr->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]xdN_{ch}/d#eta");
      if (nuName=="kaprNuDyn") gr->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]xdN_{ch}/d#eta");
    }
    if (flag==1) {
      if (nuName=="pikaNuDyn") gr->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]xN_{part}");
      if (nuName=="piprNuDyn") gr->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]xN_{part}");
      if (nuName=="kaprNuDyn") gr->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]xN_{part}");
    }
     if (flag==4) {
      if (nuName=="pikaNuDyn") gr->GetYaxis()->SetTitle(" #frac{#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]}{1/<#pi>+1/<K>} ");
      if (nuName=="piprNuDyn") gr->GetYaxis()->SetTitle(" #frac{#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]}{1/<#pi>+1/<p>} ");
      if (nuName=="kaprNuDyn") gr->GetYaxis()->SetTitle(" #frac{#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]}{1/<p>+1/<K>} ");
    }
    if (flag==5 || flag==6) gr->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    
    gr->SetFillColor(kBlack);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlack);
    
    scaling->GetFile()->cd();
    gr->Write();
    delete gr;
    
    
  }
  
  delete scaling;
  
} 
// =====================================================================================================
Double_t GetNuDyn(Double_t a, Double_t b, Double_t a2, Double_t b2, Double_t ab )
{
    
    if (a&&b) {
        return a2/(a*a)+b2/(b*b)-2*ab/(a*b)-(1/a+1/b);
    } else {
        cout << "ERROR: Problem with the first moments" << endl;
        return 0.;
    }
    
}
// =====================================================================================================
Double_t GetNetParCumulant2(Double_t a, Double_t b, Double_t a2, Double_t b2, Double_t ab )
{
    
    if (a&&b) {
        return (a2-a*a)+(b2-b*b)-2*(ab-a*b);
    } else {
        cout << "ERROR: Problem with the first moments" << endl;
        return 0.;
    }
    
}






















