// #include <TFile.h>
// #include "TBranch.h"
// #include <iostream>
// #include "TMath.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "TTreeStream.h"
// #include "TFrame.h"
// #include "TROOT.h"
// #include "TKey.h"
// #include "TVirtualFitter.h"
// #include "TClass.h"
// #include "TF1.h"
// #include "TTree.h"
// #include "TStyle.h"
// #include "TGaxis.h"
// #include "TLegend.h"
// #include "TCut.h"
// #include "TCanvas.h"
// #include "TGraphErrors.h"
// #include "TStatToolkit.h"
// #include "TLegend.h"
// #include "TLegendEntry.h"
// #include "AliXRDPROOFtoolkit.h"
// #include "TChain.h"
// #include <fstream>
// #include <iostream>
// #include <iomanip>
using namespace std;
using std::cout;
using std::setw;


void FitHistogram(TH1D *hClean);
void ProcessErrorCalculationVScent(TString momentsFile,Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp);
void ProcessErrorCalculationVSeta(TString momentsFile,Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp);
void MakeRatioPlots(TString mfile,Int_t scanType, Double_t pDown, Double_t pUp, Double_t etaDown, Double_t etaUp);
void RatioScan(TString etaScanFile,TString momScanFile );
void SetYTitleName(TString branchName);


void FastPlotVScent(Int_t style, Int_t col, TString momentsFile, Int_t dataType, TString var, Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp, Int_t same);

TGraphErrors * RatioPlot(TTreeSRedirector *ratiofile, TGraphErrors *g1, TGraphErrors *g2, Int_t scanType, Double_t pDown, Double_t pUp, Double_t etaDown, Double_t etaUp);
Double_t       CalculateStatErrors(TGraphErrors *grS);


TTreeSRedirector *MCssMomets=NULL;
TTreeSRedirector *statResults=NULL;
TTreeSRedirector *ratioresults=NULL;
TTreeSRedirector *correctionLookup=NULL;
TString yTitleLatex;                  // title name should be wrtiiten in latex format


// ******************************* Modification region ******************************************

Int_t nSubsample     = 16;
Bool_t lowStat       = kFALSE;    // USed only for W analysis
Bool_t fpp           = kFALSE;
Bool_t fFastGen      = kFALSE;
Bool_t fNewTreeStyle = kTRUE;
Double_t pDownArray[4]    = {0.2, 0.2, 0.3, 0.3};
Double_t pUpArray[4]      = {1.5, 1.8, 1.5, 1.8};
Double_t etaDownArray[2]  = {-0.5, -0.8};
Double_t etaUpArray[2]    = { 0.5,  0.8};
Double_t centDownArray[9] = {0.,5. ,10.,20.,30.,40.,50.,60.,70.};
Double_t centUpArray[9]   = {5.,10.,20.,30.,40.,50.,60.,70.,80.};

// **********************************************************************************************

void TIdentitySubSampleMC(TString momentsFile, Int_t sampNo, Int_t dataType, Double_t pDown, Double_t pUp, Double_t etaDown, Double_t etaUp ){

  //
  // Calculate statistical errors for MCgen and MCrec form "MC_recgen_moments.root" using subsample method
  // 
                                                                                       
  // To RUN
  /*       
   
   cd /hera/alice/marsland/pFluct/files/analysis/Data/PbPb/StronglyIntensive/HIJING
   aliroot -l -b
   TString h = "/hera/alice/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_Fastgen_HIJING_dnchdeta/mergedPeriods/AnalysisResults.root"
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
   TIdentitySubSampleMC(h,0,1,0.2,1.5,-0.8,0.8)
   
   cd /hera/alice/marsland/pFluct/files/analysis/Data/PbPb/StronglyIntensive/AMPTsmONrsON
   aliroot -l -b
   TString a1 = "/hera/alice/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_Fastgen_AMPT_LHC13f3c_StringMelting_ON_Rescattering_ON_dnchdeta/mergedPeriods/AnalysisResults.root"
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
   TIdentitySubSampleMC(a1,0,1,0.2,1.5,-0.8,0.8) 
   
   cd /hera/alice/marsland/pFluct/files/analysis/Data/PbPb/StronglyIntensive/AMPTsmOFFrsON
   aliroot -l -b
   TString a2 = "/hera/alice/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_Fastgen_AMPT_LHC13f3b_StringMelting_OFF_Rescattering_ON_dnchdeta/mergedPeriods/AnalysisResults.root"
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
   TIdentitySubSampleMC(a2,0,1,0.2,1.5,-0.8,0.8) 
   
   cd /hera/alice/marsland/pFluct/files/analysis/Data/PbPb/StronglyIntensive/AMPTsmONrsOFF
   aliroot -l -b
   TString a3 = "/hera/alice/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_Fastgen_AMPT_LHC13f3a_StringMelting_ON_Rescattering_OFF_dnchdeta/mergedPeriods/AnalysisResults.root"   
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
   TIdentitySubSampleMC(a3,0,1,0.2,1.5,-0.8,0.8) 

   cd /hera/alice/marsland/pFluct/files/analysis/Data/PbPb/StronglyIntensive/HIJING_Closure
   aliroot -l -b
   TString closureFile = "/hera/alice/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_Closure_Test_HIJING/mergedPeriods/marsland_MCmoments.root"   
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
   TIdentitySubSampleMC(closureFile,0,1,0.2,1.5,-0.8,0.8) 
   TIdentitySubSampleMC(closureFile,0,0,0.2,1.5,-0.8,0.8) 
   
   TIdentitySubSampleMC(closureFile,0,1,0.6,1.5,-0.8,0.8) 
   TIdentitySubSampleMC(closureFile,0,0,0.6,1.5,-0.8,0.8) 

   TIdentitySubSampleMC(closureFile,0,1,0.6,1.8,-0.8,0.8) 
   TIdentitySubSampleMC(closureFile,0,0,0.6,1.8,-0.8,0.8) 
   
   
   cd /hera/alice/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_Closure_HIJING_EffMAtrix_TOFcutON/EfficiencyCheck
   aliroot -l -b
   TString newTOFon = "/hera/alice/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_Closure_HIJING_EffMAtrix_TOFcutON/mergedPeriods/AnalysisResults.root"
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
   TIdentitySubSampleMC(newTOFon,0,1,0.2,1.5,-0.8,0.8) 
   TIdentitySubSampleMC(newTOFon,0,0,0.2,1.5,-0.8,0.8) 
   TIdentitySubSampleMC(newTOFon,0,1,0.6,1.5,-0.8,0.8) 
   TIdentitySubSampleMC(newTOFon,0,0,0.6,1.5,-0.8,0.8) 
   TIdentitySubSampleMC(newTOFon,0,1,0.6,1.8,-0.8,0.8) 
   TIdentitySubSampleMC(newTOFon,0,0,0.6,1.8,-0.8,0.8) 
   

   
   TIdentitySubSampleMC(momFile,1,0,0.2,1.5,-0.8,0.8) 
  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C .");
  MCssMomets = new TTreeSRedirector(Form("MomentsTree_DT_samp%d_gen%d_eta_%2.1f-%2.1f_mom_%2.1f-%2.1f.root",sampNo,dataType,etaDown,etaUp,pDown,pUp),"recreate");

  // Get the trees 
  TTree *tree;
  TString treeName = (dataType==0) ? "mcRec": "mcGen";

  if (momentsFile.Contains(".root")){
    TFile *f = TFile::Open(momentsFile);
    cout << " get the ttree " << endl;
    tree = (TTree*)f->Get(treeName);
  } else {
    TChain *chain;
    cout << " make the chain " << endl;
    tree = (TTree*)AliXRDPROOFtoolkit::MakeChainRandom(momentsFile,treeName,0,500000,0);  
  }
  if (!tree) return;
  
  // Second array to hold centDown
  Int_t nPoints = sizeof(centDownArray)/sizeof(centDownArray[0]);   
  cout << " nPoints == " << nPoints << endl;   
  
  // loop over centrality bins
  cout << " =========== subsample = " << sampNo << " =========== " << endl;
  for (Int_t icent=0;icent<nPoints;icent++){
    
    // Moment arrays for each centrality bin
    Double_t trCount=0.;
    Double_t firstMoments[4];
    Double_t secondMoments[4][4];
    Double_t nudynArr[4][4];
    Double_t nudynDyn[4][4];
    Double_t nudynStat[4][4];
    for (Int_t i=0; i<4;i++){
      firstMoments[i]=0.;
      for (Int_t j=0;j<4;j++){
	secondMoments[i][j]=0.;
	nudynArr[i][j]=0.;
	nudynDyn[i][j]=0.; 
	nudynStat[i][j]=0.; 
      }
    }
    TString percentSign = "%";
    TCut centCut, etaMomCut;
    if (fFastGen) {
      centCut   = Form("impPar<14.5&&centDown==%f && dataType==%d",centDownArray[icent],dataType);
      etaMomCut = Form("impPar<14.5&&abs(pDown-%f)<0.001 && abs(pUp-%f)<0.001 && abs(etaUp-%f)<0.001",pDown,pUp,etaUp);
    } else {
      centCut   = Form("centDown==%f && dataType==%d",centDownArray[icent],dataType);
      etaMomCut = Form("abs(pDown-%f)<0.001 && abs(pUp-%f)<0.001 && abs(etaUp-%f)<0.001",pDown,pUp,etaUp);
    }
    TCut sampCut   = (sampNo==0) ? "isample==0" : Form("Entry$%s%d==%d",percentSign.Data(),nSubsample-1,sampNo-1);
    
    // get list of the branches to loop over them
    TObjArray *branchArr   =  (TObjArray*)tree->GetListOfBranches(); 
    Int_t nEntries = branchArr->GetEntriesFast();
    for (Int_t iBranch=0; iBranch<nEntries; ++iBranch){
      
      TString branchName = branchArr->At(iBranch)->GetName();
      if ( branchName.Contains("cent")      || branchName.Contains("eta")      || branchName.Contains("ptot") ) continue;
      if ( branchName.Contains("isample")   || branchName.Contains("dataType") || branchName.Contains("event") ) continue;
      if ( branchName.Contains("pDown")     || branchName.Contains("pUp")      || branchName.Contains("pBin") ) continue;
      
      // plot the graph hist for each moment and get mean

      TH1D *htemp[11];
      for (Int_t i=0;i<11;i++) htemp[i]=NULL;
      if (fNewTreeStyle) {
	
	if (branchName.Contains("trCount")) 
	{
	  tree->Draw(branchName,centCut && etaMomCut && sampCut,"goff");
	  htemp[0] = (TH1D*)tree->GetHistogram()->Clone();
	  trCount = htemp[0]->GetMean();
	}
	if (branchName.Contains("moment.")) 
	{
	  // 	  kPi=0,kKa=1,kPr=2,kPiPi=3,kKaKa=4,kPrPr=5,kPiKa=6,kPiPr=7,kKaPr=8
	  tree->Draw("moment.fElements[0]",centCut && etaMomCut && sampCut,"goff");
	  htemp[1] = (TH1D*)tree->GetHistogram()->Clone();
	  firstMoments[1]  = htemp[1]->GetMean();
	  
	  tree->Draw("moment.fElements[1]",centCut && etaMomCut && sampCut,"goff");
	  htemp[2] = (TH1D*)tree->GetHistogram()->Clone();
	  firstMoments[2]  = htemp[2]->GetMean();
	
	  tree->Draw("moment.fElements[2]",centCut && etaMomCut && sampCut,"goff");
	  htemp[3] = (TH1D*)tree->GetHistogram()->Clone();
	  firstMoments[3]  = htemp[3]->GetMean();
	
	  tree->Draw("moment.fElements[3]",centCut && etaMomCut && sampCut,"goff");
	  htemp[4] = (TH1D*)tree->GetHistogram()->Clone();
	  secondMoments[1][1] = htemp[4]->GetMean();
	  
	  tree->Draw("moment.fElements[4]",centCut && etaMomCut && sampCut,"goff");
	  htemp[5] = (TH1D*)tree->GetHistogram()->Clone();
	  secondMoments[2][2] = htemp[5]->GetMean();
	  
	  tree->Draw("moment.fElements[5]",centCut && etaMomCut && sampCut,"goff");
	  htemp[6] = (TH1D*)tree->GetHistogram()->Clone();
	  secondMoments[3][3] = htemp[6]->GetMean();
	
	  tree->Draw("moment.fElements[6]",centCut && etaMomCut && sampCut,"goff");
	  htemp[7] = (TH1D*)tree->GetHistogram()->Clone();
	  secondMoments[1][2] = htemp[7]->GetMean();
	  
	  tree->Draw("moment.fElements[7]",centCut && etaMomCut && sampCut,"goff");
	  htemp[8] = (TH1D*)tree->GetHistogram()->Clone();
	  secondMoments[1][3] = htemp[8]->GetMean();
	  
	  tree->Draw("moment.fElements[8]",centCut && etaMomCut && sampCut,"goff");
	  htemp[9] = (TH1D*)tree->GetHistogram()->Clone();
	  secondMoments[2][3] = htemp[9]->GetMean();
	    
	}
	  
      } else {
	tree->Draw(branchName,centCut && etaMomCut && sampCut,"goff");
	htemp[10] = (TH1D*)tree->GetHistogram()->Clone();
	// get first moments
	if (branchName.Contains("piCount")) firstMoments[1]  = htemp[10]->GetMean();
	if (branchName.Contains("kaCount")) firstMoments[2]  = htemp[10]->GetMean();
	if (branchName.Contains("prCount")) firstMoments[3]  = htemp[10]->GetMean();
	if (branchName.Contains("pipi")) secondMoments[1][1] = htemp[10]->GetMean();
	if (branchName.Contains("kaka")) secondMoments[2][2] = htemp[10]->GetMean();
	if (branchName.Contains("prpr")) secondMoments[3][3] = htemp[10]->GetMean();
	if (branchName.Contains("pika")) secondMoments[1][2] = htemp[10]->GetMean();
	if (branchName.Contains("pipr")) secondMoments[1][3] = htemp[10]->GetMean();
	if (branchName.Contains("kapr")) secondMoments[2][3] = htemp[10]->GetMean();
	if (branchName.Contains("trCount")) trCount = htemp[10]->GetMean();
      }
      
      for (Int_t i=0;i<11;i++) { if(htemp[i]) delete htemp[i];} 
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
    
    
    // Strongly intensive quantities
    
    Double_t pi1 = firstMoments[1];
    Double_t ka1 = firstMoments[2];
    Double_t pr1 = firstMoments[3];
    Double_t pi2 = secondMoments[1][1];
    Double_t ka2 = secondMoments[2][2];
    Double_t pr2 = secondMoments[3][3];
    Double_t pika = secondMoments[1][2];
    Double_t pipr = secondMoments[1][3];
    Double_t kapr = secondMoments[2][3];
    
    Double_t wpi = (pi2-pi1*pi1)/pi1;
    Double_t wka = (ka2-ka1*ka1)/ka1;
    Double_t wpr = (pr2-pr1*pr1)/pr1;
    
    Double_t cDeltapika = pi1-ka1;
    Double_t cDeltapipr = pi1-pr1;
    Double_t cDeltakapr = ka1-pr1;
    
    Double_t cSigmapika = pi1+ka1;
    Double_t cSigmapipr = pi1+pr1;
    Double_t cSigmakapr = ka1+pr1;
    
    Double_t deltapika = (1./cDeltapika)*(pi1*wka - ka1*wpi);
    Double_t deltapipr = (1./cDeltapipr)*(pi1*wpr - pr1*wpi);
    Double_t deltakapr = (1./cDeltakapr)*(ka1*wpr - pr1*wka);
    
    
    Double_t sigmapika = (1./cSigmapika)*(pi1*wka + ka1*wpi - 2*(pika-pi1*ka1));
    Double_t sigmapipr = (1./cSigmapipr)*(pi1*wpr + pr1*wpi - 2*(pipr-pi1*pr1));
    Double_t sigmakapr = (1./cSigmakapr)*(ka1*wpr + pr1*wka - 2*(kapr-ka1*pr1));
    
    Double_t phipika = (TMath::Sqrt(pi1*ka1)/(pi1+ka1))*(TMath::Sqrt(sigmapika)-1);
    Double_t phipipr = (TMath::Sqrt(pi1*pr1)/(pi1+pr1))*(TMath::Sqrt(sigmapipr)-1);
    Double_t phikapr = (TMath::Sqrt(ka1*pr1)/(ka1+pr1))*(TMath::Sqrt(sigmakapr)-1);

    
    
    // dump everything into ttree   
    MCssMomets -> GetFile()->cd();
    *MCssMomets << "mcMoments" <<
    
    "dataType="     << dataType            <<
    "subsample="    << sampNo               <<
    "pDown="        << pDown               << 
    "pUp="          << pUp                 << 
    "pBin="         << pBin                << 
    "etaDown="      << etaDown             << 
    "etaUp="        << etaUp               << 
    "etaBin="       << etaBin              << 
    "centDown="     << centDownArray[icent]<< 
    "centUp="       << centUpArray[icent]  << 
    "centBin="      << centBin             << 
    "trCount="      << trCount             <<
    
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
    
    "pikaNuDynDyn="    << nudynDyn[1][2]      <<
    "piprNuDynDyn="    << nudynDyn[1][3]      <<
    "kaprNuDynDyn="    << nudynDyn[2][3]      <<
    
    "pikaNuStat="    << nudynStat[1][2]      <<
    "piprNuStat="    << nudynStat[1][3]      <<
    "kaprNuStat="    << nudynStat[2][3]      ;
    
    *MCssMomets << "mcMoments" <<
    
    // Variance etc
    "wpi="       <<    wpi        <<
    "wka="       <<    wka        <<
    "wpr="       <<    wpr        <<
    
    "cDeltapika=" <<   cDeltapika <<
    "cDeltapipr=" <<   cDeltapipr <<
    "cDeltakapr=" <<   cDeltakapr <<
    
    "cSigmapika=" <<   cSigmapika <<
    "cSigmapipr=" <<   cSigmapipr <<
    "cSigmakapr=" <<   cSigmakapr <<
    
    "nudynpika="    << nudynArr[1][2]      <<
    "nudynpipr="    << nudynArr[1][3]      <<
    "nudynkapr="    << nudynArr[2][3]      <<
  
    // strongly intensive quantities
    "deltapika=" <<    deltapika  << 
    "deltapipr=" <<    deltapipr  << 
    "deltakapr=" <<    deltakapr  << 
    
    "sigmapika=" <<    sigmapika  << 
    "sigmapipr=" <<    sigmapipr  << 
    "sigmakapr=" <<    sigmakapr  << 
    
    "phipika="   <<    phipika  << 
    "phipipr="   <<    phipipr  << 
    "phikapr="   <<    phikapr  << 
    
    
    "\n"; 
    
  } // end of centrality loop
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
  cd /hera/alice/marsland/pFluct/files/analysis/IdentityFiles/test
  tpcdev
  aliroot -l -b
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
  ProcessErrorCalculationVScent("MomentsTree.root",-0.8,0.8,0.2,2)
  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C .");

  
  // Output file
  TString statFileName = Form("CentScan_eta_%2.1f-%2.1f_mom_%2.1f-%2.1f.root",etaDown,etaUp,pDown,pUp); 
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
   
      drawMoment.Append(":centBin");
      drawSubSam.Append(":subsample");
//       TCut allStat   = Form("subsample==0 && dataType==%d && %s!=0",idatatype,branchName.Data());
      TCut allStat   = Form("subsample==0 && dataType==%d",idatatype);
      TCut subsamCut = "subsample>0 && subsample<26 && centBin==%f && dataType==%d";
      TCut etaMomCut = Form("etaDown==%f && etaUp==%f && pDown==%f && pUp==%f"   ,etaDown, etaUp, pDown, pUp);

    // look at only amplitude Mean and Sigma
      if ( (branchName.Contains("cent")  || branchName.Contains("subsample") || branchName.Contains("dataType")) ) continue;
      if ( (branchName.Contains("Up")    || branchName.Contains("Down")      || branchName.Contains("Bin")) ) continue;
      if ( (branchName.Contains("pT")) ) continue;

    
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
      grMoments[objIndex]->GetXaxis()->SetTitle("centrality [%]");
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
  
  MakeRatioPlots(statFileName,10,pDown,pUp,etaDown,etaUp);
  gSystem->Exec(Form("rm %s",statFileName.Data()));
 
}
// =====================================================================================================
void ProcessErrorCalculationVSeta(TString momentsFile,Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp){
  
  //
  // Calculate statistical errors using subsample method
  // First merge all moment outputs
  // hadd Tree_TImoments_All.root cent_#/#/TI#.root 
  // 
                                                                                       
  // To RUN
  /*                                        
  cd /hera/alice/marsland/pFluct/files/analysis/IdentityFiles/test
  tpcdev
  aliroot -l -b
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
  ProcessErrorCalculationVScent("MomentsTree.root",-0.8,0.8,0.2,2)
  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalErrorMC.C .");

  
  // Output file
  TString statFileName = Form("CentScan_eta_%2.1f-%2.1f_mom_%2.1f-%2.1f.root",etaDown,etaUp,pDown,pUp); 
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
   
      drawMoment.Append(":abs(etaDown)+etaUp");
      drawSubSam.Append(":subsample");
//       TCut allStat   = Form("subsample==0 && dataType==%d && %s!=0",idatatype,branchName.Data());
      TCut allStat   = Form("subsample==0 && dataType==%d",idatatype);
      TCut subsamCut = "subsample>0 && subsample<26 && centBin==%f && dataType==%d";
      TCut etaMomCut = Form("abs(etaDown-%f)<0.001 && abs(etaUp-%f)<0.001 && abs(pDown-%f)<0.001 && abs(pUp-%f)<0.001",etaDown, etaUp, pDown, pUp);

    // look at only amplitude Mean and Sigma
      if ( (branchName.Contains("cent")  || branchName.Contains("subsample") || branchName.Contains("dataType")) ) continue;
      if ( (branchName.Contains("Up")    || branchName.Contains("Down")      || branchName.Contains("Bin")) ) continue;
      if ( (branchName.Contains("pT")) ) continue;

    
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
      grMoments[objIndex]->GetXaxis()->SetTitle("centrality [%]");
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
  
  MakeRatioPlots(statFileName,10,pDown,pUp,etaDown,etaUp);
  gSystem->Exec(Form("rm %s",statFileName.Data()));
 
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
    grMC   ->SetMarkerStyle(20);  
    grMC   ->SetMarkerColor(kRed);   
    grMCgen->SetMarkerStyle(21);  
    grMCgen->SetMarkerColor(kBlue);
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
          grMCratio->GetXaxis()->SetTitle("p (GeV/c)");
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
  
  
} 
// =====================================================================================================
void FastPlotVScent(Int_t style, Int_t col, TString momentsFile, Int_t dataType, TString var, Double_t etaDown, Double_t etaUp, Double_t pDown, Double_t pUp, Int_t same){
 
  /*
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
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
void RatioScan(TString etaScanFile,TString momScanFile ){

  // To RUN
  /*                                        
  cd /hera/alice/marsland/pFluct/files/analysis/IdentityFiles/MC_genVSrec/SSanalysis/EtaMomentumScan/cRows80/centBinWidth_10_all/Results
  tpcdev
  aliroot -l -b
  .L ~/PHD/macros/marsland_EbyeRatios/TIdentitySubSampleMC.C+
  RatioScan("rootFiles/EtaScan_MomentsTree.root","rootFiles/MomScan_MomentsTree.root")
  */
  //
  
  correctionLookup = new TTreeSRedirector("CorrectionLookup.root","recreate");
  
//   Int_t nMomentumBins = (sizeof pDownArray)/(sizeof pDownArray[0]);
  Int_t nEtaBins      = (sizeof etaDownArray)/(sizeof etaDownArray[0]);
  Int_t nCentBins     = (sizeof centDownArray)/(sizeof centDownArray[0]);
  
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

  Double_t tmpSum = 0.;
  Int_t nonzeroCount = 0;
  for (Int_t i = 0; i < grS->GetN(); i++) {
    if (grS->GetY()[i]<0.000001) continue;
    nonzeroCount++;
    tmpSum+=TMath::Abs(grS->GetY()[i]);
  }
  if (nonzeroCount==0) {
    cout << "aaaaaaaa  " << grS->GetName() << endl;
    return 0.;
  }
  Double_t meanMS = tmpSum/Double_t(nonzeroCount);
  
  Double_t n      = grS->GetN();
  Double_t sum = 0;
  for (Int_t i = 0; i < n; i++) {
    if (grS->GetY()[i]<0.000001) continue; // ?????
    Double_t y= TMath::Abs(grS->GetY()[i]);
    sum = sum+(y-meanMS)*(y-meanMS);
  }
  Double_t err = (nonzeroCount!=0) ? TMath::Sqrt(sum/(n*(n-1))) : 0.;
 
  
  
  
  
//   Double_t meanMS = TMath::Abs(grS->GetMean(2));
//   Double_t n      = grS->GetN();
//   Double_t sum = 0;
//     
//   for (Int_t i = 0; i < n; i++) {
//     if (grS->GetY()[i]<0.000001) continue; // ?????
//     Double_t y= TMath::Abs(grS->GetY()[i]);
//     sum = sum+(y-meanMS)*(y-meanMS);
//   }
//   Double_t err = (n!=0) ? TMath::Sqrt(sum/(n*(n-1))) : 0.;
  
  if (err) return err;
  else return 0.;
  
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
TGraphErrors * RatioPlot(TTreeSRedirector *ratiofile, TGraphErrors *g1, TGraphErrors *g2, Int_t scanType, Double_t pDown, Double_t pUp, Double_t etaDown, Double_t etaUp) {
  
  const Int_t N = g1->GetN();
  Double_t x[N];
  Double_t y[N];
  Double_t errx[N];
  Double_t erry[N];
  
  for (Int_t i=0; i<N; i++){
    if ((g2->GetY()[i]>0.)) {
      y[i]=g1->GetY()[i]/g2->GetY()[i];
    } else {
      y[i]=0.;
    }
    x[i]=g1->GetX()[i];
    errx[i]=0.;
    erry[i]=0.;
  }
  
  TLine* line = new TLine(0., 1.,80., 1.);
  line -> SetLineColor(1); line -> SetLineWidth(4); line -> SetLineStyle(2);
       
  TGraphErrors *grRatio = new TGraphErrors(N,x,y,errx,erry);  
  grRatio->SetMinimum(-100);
  
  TLegend *leg = new TLegend(0.2,0.7,0.6,0.8);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);  // White
  
  // plot ratio and graphs together in one canvas  
  TCanvas *c1 = new TCanvas("c1","example",600,700);
  c1->SetBottomMargin(20.); c1->SetLeftMargin(20.);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0.); pad1->SetLeftMargin(0.14); 
  pad1->Draw();
  pad1->cd();
  pad1->SetTicks(1,1);
  
  // Calculate y axis ranges 
  Double_t max1 = TMath::MaxElement(g1->GetN(),g1->GetY()); 
  Double_t max2 = TMath::MaxElement(g2->GetN(),g2->GetY()); 
  Double_t maxx = TMath::Max(max1,max2)*1.3;
  Double_t min1 = TMath::MinElement(g1->GetN(),g1->GetY()); 
  Double_t min2 = TMath::MinElement(g2->GetN(),g2->GetY()); 
  Double_t minn = TMath::Min(min1,min2)*0.7;

  g1->Draw("alp");
  g1->GetYaxis()->SetRangeUser(minn,maxx);
  g1->GetYaxis()->SetTitleOffset(1.45);
  g2->Draw("lp");
  leg->AddEntry(g1,"MC reconstruced","lep");
  leg->AddEntry(g2,"MC generated","lep");
  TLegendEntry *header;
  if (scanType>1){
    leg->SetHeader(Form("%2.1f<#eta<%2.1f, %2.1f<p<%2.1f",etaDown,etaUp,pDown,pUp));
    header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    header->SetTextSize(0.045);
  }
  leg->Draw();
  
  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
  pad2->SetTopMargin(0.); pad2->SetLeftMargin(0.14); pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();
  pad2->SetTicks(1,1);
  grRatio->SetMarkerStyle(22); 
  grRatio->SetMarkerColor(kBlack);
  grRatio->GetXaxis()->SetTitle(g1->GetXaxis()->GetTitle()); grRatio->GetXaxis()->SetTitleOffset(1.);
  grRatio->GetYaxis()->SetTitle("rec/gen");
   
  // Calculate y axis ranges 
  Double_t max = TMath::MaxElement(grRatio->GetN(),grRatio->GetY())*1.3; 
  Double_t min = TMath::MinElement(grRatio->GetN(),grRatio->GetY())*0.7; 
  grRatio->GetYaxis()->SetRangeUser(min,max);
  
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
  c1->cd();  
  
  TString canvasname = g1->GetName();
  ratiofile->GetFile()->cd();
  c1->Write(Form("canvas_%s",canvasname.Data()));
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
  if (branchName=="pika") yTitleLatex="#LT#pi,K#GT";
  if (branchName=="pipr") yTitleLatex="#LT#pi,p#GT";
  if (branchName=="kapr") yTitleLatex="#LTK,p#GT";
  if (branchName=="pipi") yTitleLatex="#LT#pi,#pi#GT";
  if (branchName=="kaka") yTitleLatex="#LTK,K#GT";
  if (branchName=="prpr") yTitleLatex="#LTp,p#GT";
  if (branchName=="pi") yTitleLatex="#LT#pi#GT";
  if (branchName=="ka") yTitleLatex="#LTK#GT";
  if (branchName=="pr") yTitleLatex="#LTp#GT";
          
}
