#include <TFile.h>
#include "TBranch.h"
#include <iostream>
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeStream.h"
#include "TFrame.h"
#include "TROOT.h"
#include "TKey.h"
#include "TVirtualFitter.h"
#include "TPaveLabel.h"
#include "TClass.h"
#include "TF1.h"
#include "TTree.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "AliXRDPROOFtoolkit.h"
#include "TChain.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::cout;
using std::setw;

void     FinalSyst();
void     GraphShade();
void     ErrorBandTH1();
void     AnalyseHistsInFile(TString file);
void     FitHistogram(TH1D *hClean);
void     AnalyseWs(TString wFile, TString wSliceFile);
void     AnalyseWsApprovedPlots(TString wFile, TString wSliceFile);
void     SetYTitleName(TString branchName);
void     CompareToReference(TString refFile, TString xFile, TString nuName);
void     CompareToReferenceSyst(TString refFile, TString xFile, Int_t iterREF, Int_t iterX, TString outName);
void     CompareDataVsMC(TString dataFile, TString hijingFile, TString amptFile);
void     CompareIdentityVsTraditional(TString dataFile, TString tradFile);


TCanvas  *TwoScales(TTree *t);
TCanvas  *OverlayPads(TTree *t);
TGraphErrors * CorrectCentBinOfMoments(TGraphErrors *grCorr);
TGraphErrors * MakeNpartMoments(TGraphErrors *grCorr);
TGraphErrors * MakeScaledMoments(TGraphErrors *grCorr, TGraphErrors *grpi);
TGraphErrors * DiffGraph(TGraphErrors *g1, TGraphErrors *g2);
TGraphErrors * RatioGraph(TGraphErrors *g1, TGraphErrors *g2);

Double_t       CalculateStatErrors(TGraphErrors *grS);
void           MultiplydNch(TString datafile,Int_t flag);
void           PlotNpartScalings(TString data, TString hijing, TString ampt);
void           PlotdnchScalings(TString data, TString hijing, TString ampt);
void           MultiplydNchStatAndSys(TString grStatName, TString grSystName);

// const Int_t nCentbins = 9; Float_t centArray[nCentbins] = {0,5,10,20,30,40,50,60,70};
// const Int_t ncentbins = 10; Float_t xCentBins[ncentbins] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
// const Int_t nCentbins = 18; Float_t centArray[nCentbins] = {0, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75}; 
// const Int_t ncentbins = 19; Float_t xCentBins[ncentbins] = {0, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80};

Double_t centArr[9]  = {2.5,7.5,15,25,35,45,55,65 ,75 };
Double_t centArr2[9] = {75 ,65 ,55,45,35,25,15,7.5,2.5};
// Double_t centArr[18]  = {1.25, 3.75, 6.25, 8.75, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5};
// Double_t centArr2[18] = {77.5 ,72.5 ,67.5, 62.5, 57.5, 52.5, 47.5, 42.5, 37.5, 32.5, 27.5, 22.5, 17.5, 12.5, 8.75, 6.25, 3.75, 1/25};


// nCharge and nPart arrays from publication
Double_t dnchArray[9]    = {1601 ,1294 ,966 ,649 ,426 ,261 ,149 ,76 ,35};  // eta <0.5
Double_t dnchArrayErr[9] = {60   ,49   ,37  ,23  ,15  ,9   ,6   ,4  ,2};
Double_t nPartArr[9]     = {382.8, 329.7, 260.5, 186.4, 128.9, 85.0, 52.8, 30.0, 15.8};
Double_t nPartArrErr[9]  = {3.1, 4.6, 4.4, 3.9, 3.3, 2.6, 2.0, 1.3, 0.6};

Bool_t savePDF     = kFALSE;
Bool_t noMostPerip = kFALSE;
Bool_t lowStat     = kFALSE;
Bool_t noErrorBar  = kTRUE;
Bool_t pp          = kFALSE;
Bool_t transparent = kFALSE;

TLatex * text;

/* 
 
 aliroot -l 
 TH2F *hist = new TH2F("fHist","fHist", 100,0.,28., 100,-0.005,0.005);
 hist->GetYaxis()->SetTitleOffset(1.5);  
 hist->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
//  hist->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
//  hist->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
 hist->GetXaxis()->SetTitle("subsample index"); 
 hist->Draw();
 //TString data="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/TIdenResults/TIdenResults_EtaMomScan_2015_03_14_13_10/Eta_-0.8_0.8_Mom_0.2_1.5/TIdenMoments.root"
 //TString data="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/TIdenResults/TIdenResults_EtaMomScan_le30centbinsDownscaled_GOOD/Eta_-0.8_0.8_Mom_0.2_1.5/TIdenMoments.root" 
 TString data="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics_cRows_80_16EtaBin_mombin20MeV/TIdenResults/TIdenResults_piS0.7_kaS0.5_prS0.5_pikaprKauto/Eta_-0.8_0.8_Mom_0.2_1.5/TIdenMoments.root" 
 TFile f(data)
 momTree->SetMarkerStyle(20)
 momTree->Draw("pikaNuDyn:subsample+1","cent==0","same")
 momTree->SetMarkerColor(kRed)
 momTree->Draw("pikaNuDyn:subsample+1","cent==0&&subsample==0","same")
 
 ****************************************************************************
 
 aliroot -l 
 gStyle->SetOptStat(0);
 TString smallCent = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_18CentBin_16EtaBin_mombin20MeV_TightCuts/TIdenResults/TIdenResults_sign0_piS0.7_kaS0.5_prS0.5_prK2_kapiKAutoKpi_allSS/Eta_-0.8_0.8_Mom_0.2_1.5/TIdenMoments.root"
 //TString smallCent = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_18CentBin_16EtaBin_mombin20MeV_TightCuts/TIdenResults/TIdenResults_EtaMomScan_2015_08_20_02_52/Eta_-0.8_0.8_Mom_0.2_1.5/TIdenMoments.root"
 TString largeCent = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/TIdenResults/TIdenResults_EtaMomScan_le30centbinsDownscaled_GOOD/Eta_-0.8_0.8_Mom_0.2_1.5/TIdenMoments.root"
 TFile *small = TFile::Open(smallCent);
 TFile *large = TFile::Open(largeCent);
 TTree *treesmall = (TTree*)small->Get("momTree");
 TTree *treelarge = (TTree*)large->Get("momTree");
 treelarge->SetMarkerColor(kBlack);
 treesmall->SetMarkerColor(kRed);
 treelarge->SetMarkerStyle(20);
 treesmall->SetMarkerStyle(20);
 treelarge->Draw("pikaNuDyn:cent+5","subsample>0","prof")
 treesmall->Draw("pikaNuDyn:cent+2.5","subsample>0 && cent<70","profsame")
 new TCanvas
 treelarge->Draw("prkaNuDyn:cent+5","subsample>0","prof")
 treesmall->Draw("prkaNuDyn:cent+2.5","subsample>0 && cent<70","profsame")
 new TCanvas
 treelarge->Draw("prpiNuDyn:cent+5","subsample>0","prof")
 treesmall->Draw("prpiNuDyn:cent+2.5","subsample>0 && cent<70","profsame")

 */

void TIdentityStatisticalError(TString momentsFile,Int_t iter=0)
{

  //
  // Calculate statistical errors using subsample method
  // First merge all moment outputs
  // hadd Tree_TImoments_All.root cent_#/#/TI#.root 
  // 
                                                                                       
  // To RUN
  /*                                        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/test
  aliroot -l 
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
  TIdentityStatisticalError("TIdenMoments.root",3)
  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C .");

  
  // Output file
  TTreeSRedirector *statResults = new TTreeSRedirector(Form("PlotsMomentsStatErr_iter%d.root",iter),"recreate");
  
  // Get the trees
  TFile *fMoments = TFile::Open(momentsFile);  
  TTree *treeM    = (TTree*)fMoments -> Get("momTree"); 
  TTree *treeSS   = (TTree*)treeM->Clone();
  
  // GetList of branches for each tree
  TObjArray *branchM = (TObjArray*)treeM -> GetListOfBranches();
  Int_t nArrEntries  = branchM->GetEntriesFast();
  cout << " nArrEntries = " << nArrEntries << endl;

//   Output TObjArrays
  TObjArray *momentArr                 = new TObjArray(nArrEntries);  momentArr -> SetOwner(kTRUE); 
  TObjArray *momentArrSystematic       = new TObjArray(nArrEntries);  momentArrSystematic -> SetOwner(kTRUE); 
  TObjArray *momentArrSystematicNpart  = new TObjArray(nArrEntries);  momentArrSystematicNpart -> SetOwner(kTRUE); 
  TObjArray *subsamArr                 = new TObjArray(nArrEntries);  subsamArr -> SetOwner(kTRUE); 
  TObjArray *momentNpartArr            = new TObjArray(nArrEntries);  momentNpartArr -> SetOwner(kTRUE); 
  TObjArray *momentScaledArr           = new TObjArray(nArrEntries);  momentScaledArr -> SetOwner(kTRUE); 
  TObjArray *momentScaledArrNpart      = new TObjArray(nArrEntries);  momentScaledArrNpart -> SetOwner(kTRUE); 
  
  TGraphErrors *grMoments[nArrEntries];
  TGraphErrors *grMomentsCorrected[nArrEntries];
  TGraphErrors *grMomentsNpart[nArrEntries];
  TGraphErrors *grMomentsScaledMult[nArrEntries];
  TGraphErrors *grMomentsScaledMultNpart[nArrEntries];

  TGraphAsymmErrors * grMomentsSystematic[nArrEntries];
  TGraphAsymmErrors * grMomentsSystematicNpart[nArrEntries];


  
  // Create all "moment vs cent" TGraphErrors
  Int_t objIndex = 0;
  for (Int_t i=0; i<nArrEntries; ++i){
        
    // Get the brach name 
    TString branchName = branchM->At(i)->GetName();
    TString drawMoment = branchM->At(i)->GetName();
    TString drawSubSam = branchM->At(i)->GetName();
    
    if (pp){
      drawMoment.Append(":centBin+90");
    } else {
      drawMoment.Append(":centBin");
    }
    drawSubSam.Append(":subsample");
    TCut allStat, subsamCut;
    if (noMostPerip) {
      allStat   = "subsample==0 && centBin<70";         
      subsamCut = "subsample>0 && subsample<26 && abs(centBin-%f)<0.001 && centBin<70"; 
    } else {
      allStat   = "subsample==0";         
      subsamCut = "subsample>0 && subsample<26 && abs(centBin-%f)<0.001"; 
    }
    TCut outsideCut   = Form("fitIter==%d",iter);
    
    // look at only amplitude Mean and Sigma
    if ( (branchName.Contains("cent") || branchName.Contains("subsample")) ) continue;
    if ( (branchName.Contains("pMin") || branchName.Contains("pMax")) ) continue;
    if ( (branchName.Contains("etaMin") || branchName.Contains("etaMax")) ) continue;
    if ( (branchName.Contains("nEvent") || branchName.Contains("fitIter")) ) continue;
    if ( (branchName.Contains("Int") || branchName.Contains("cal")) ) continue;
    if ( (branchName.Contains("pipiNuDyn") || branchName.Contains("kakaNuDyn")) ) continue;
    if ( (branchName.Contains("prprNuDyn") || branchName.Contains("el")) ) continue;
    if (  branchName.Contains("2") ) continue;
    if ( (branchName.Contains("etaDown") || branchName.Contains("etaUp")) ) continue;
    if ( (branchName.Contains("pDown") || branchName.Contains("pUp")) ) continue;

    // create the graph object
    treeM->Draw(drawMoment,allStat && outsideCut,"goff"); //TO FIX
    Int_t N = treeM->GetSelectedRows();
    
    // Fix the ordering of centralitues in graph
    Double_t x[N];
    Double_t y[N];
    Int_t index[N];
    Double_t *tmpx = treeM->GetV2();  
    Double_t *tmpy = treeM->GetV1(); 
    TMath::Sort(N,tmpx,index);
    for (Int_t k=0; k<N;k++){ x[k] = tmpx[index[k]]; y[k] = tmpy[index[k]]; }
//     for (Int_t k=0; k<N;k++){ x[k] = centArr[k]; y[k] = tmpy[index[k]]; }

        
    SetYTitleName(branchName);
    grMoments[objIndex] = new TGraphErrors(N,x,y); 
    grMoments[objIndex]->SetDrawOption("a3");
    grMoments[objIndex]->SetName(branchName);
    grMoments[objIndex]->GetXaxis()->SetTitle("Centrality (%)");
    grMoments[objIndex]->GetYaxis()->SetTitle(branchName);
    grMoments[objIndex]->SetFillColor(kBlack);
    grMoments[objIndex]->SetMarkerStyle(20);
    grMoments[objIndex]->SetMarkerColor(kBlack);
    
    // Get the subsample
    TObjArray *subsamCentArr = new TObjArray(N);  subsamCentArr -> SetOwner(kTRUE); 
    TGraphErrors *grSS[N];
    for (Int_t is = 0; is<N; is++){
            
      if (!pp) { 
        treeSS->Draw(drawSubSam,Form(subsamCut,x[is]) && outsideCut   ,"goff");    //TO FIX
      } else {
        treeSS->Draw(drawSubSam,Form(subsamCut,0),"goff");
      }
      grSS[is] = new TGraphErrors(treeSS->GetSelectedRows(),treeSS->GetV2(),treeSS->GetV1()); grSS[is] -> SetMarkerStyle(20);
      TString sName = branchName;
      subsamCentArr -> SetName(sName.Append("_SS"));
      grSS[is]      -> SetName(sName.Append(Form("_cent%d",Int_t(x[is]))));
      subsamCentArr -> AddAt(grSS[is],is);
      
    }
    cout << "number of cent bins = " << N << endl;
    
    // Calculate the statistical errors
    for (Int_t ierr = 0; ierr<N; ierr++) {
      if(!noErrorBar) {
	grMoments[objIndex]->SetPointError(ierr,0,CalculateStatErrors(grSS[ierr]));
      }else {
	grMoments[objIndex]->SetPointError(ierr,0,0);
      }

    }
    
    // Add graphs to ObjArrays
    grMomentsCorrected[objIndex] = CorrectCentBinOfMoments(grMoments[objIndex]);
    grMomentsNpart[objIndex]     = MakeNpartMoments(grMoments[objIndex]);
    
    // Add systematic errors
    grMomentsSystematic[objIndex] = new TGraphAsymmErrors(N,centArr2,y);
    grMomentsSystematic[objIndex]->SetName(branchName);
    for (Int_t ia=0; ia<N; ia++) {
      grMomentsSystematic[objIndex]->SetPointEYhigh(ia,TMath::Abs(grMoments[objIndex]->GetY()[ia]));
      grMomentsSystematic[objIndex]->SetPointEYlow(ia,TMath::Abs(grMoments[objIndex]->GetY()[ia]));
      grMomentsSystematic[objIndex]->SetPointEXhigh(ia,1);
      grMomentsSystematic[objIndex]->SetPointEXlow(ia,1);
    }
    grMomentsSystematic[objIndex]->SetMarkerStyle(20);
    grMomentsSystematic[objIndex]->SetMarkerColor(kBlack);
    grMomentsSystematic[objIndex]->SetFillStyle(3001);
    
    grMomentsSystematicNpart[objIndex] = new TGraphAsymmErrors(N,grMomentsNpart[objIndex]->GetX(),grMomentsNpart[objIndex]->GetY());
    grMomentsSystematicNpart[objIndex]->SetName(branchName);
    grMomentsSystematicNpart[objIndex]->SetFillStyle(3001);
    
    for (Int_t ia=0; ia<N; ia++) {
      grMomentsSystematicNpart[objIndex]->SetPointEYhigh(ia,TMath::Abs(grMomentsNpart[objIndex]->GetY()[ia]));
      grMomentsSystematicNpart[objIndex]->SetPointEYlow(ia,TMath::Abs(grMomentsNpart[objIndex]->GetY()[ia]));
      grMomentsSystematicNpart[objIndex]->SetPointEXhigh(ia,1);
      grMomentsSystematicNpart[objIndex]->SetPointEXlow(ia,1);
    }
    grMomentsSystematicNpart[objIndex]->SetMarkerStyle(20);
    grMomentsSystematicNpart[objIndex]->SetMarkerColor(kBlack);
    
    // Add graphs to array and write in the file
    momentNpartArr -> AddAt(grMomentsNpart[objIndex],objIndex);
    momentArr      -> AddAt(grMomentsCorrected[objIndex],objIndex);
    momentArrSystematic -> AddAt(grMomentsSystematic[objIndex],objIndex);
    momentArrSystematicNpart -> AddAt(grMomentsSystematicNpart[objIndex],objIndex);
    subsamArr      -> AddAt(subsamCentArr,objIndex);

    cout << " moment " << objIndex << " is " << branchName << endl;
    objIndex++;
    
  }
  
  // Get the scaled nudyns
  
  TGraphErrors *pi        = (TGraphErrors*)momentArr->FindObject("pi1");
  TGraphErrors *nudynprpi = (TGraphErrors*)momentArr->FindObject("prpiNuDyn");
  TGraphErrors *nudynprka = (TGraphErrors*)momentArr->FindObject("prkaNuDyn");
  TGraphErrors *nudynpika = (TGraphErrors*)momentArr->FindObject("pikaNuDyn");
    
  TGraphErrors *npartpi        = (TGraphErrors*)momentNpartArr->FindObject("pi1");
  TGraphErrors *npartnudynprpi = (TGraphErrors*)momentNpartArr->FindObject("prpiNuDyn");
  TGraphErrors *npartnudynprka = (TGraphErrors*)momentNpartArr->FindObject("prkaNuDyn");
  TGraphErrors *npartnudynpika = (TGraphErrors*)momentNpartArr->FindObject("pikaNuDyn");
    
  grMomentsScaledMultNpart[0] = MakeScaledMoments(npartnudynprpi, npartpi);
  grMomentsScaledMultNpart[1] = MakeScaledMoments(npartnudynprka, npartpi);
  grMomentsScaledMultNpart[2] = MakeScaledMoments(npartnudynpika, npartpi);
    
  grMomentsScaledMult[0] = MakeScaledMoments(nudynprpi, pi);
  grMomentsScaledMult[1] = MakeScaledMoments(nudynprka, pi);
  grMomentsScaledMult[2] = MakeScaledMoments(nudynpika, pi);


  for (Int_t i=0;i<3;i++){
    momentScaledArr -> AddAt(grMomentsScaledMult[i],i);
    momentScaledArrNpart -> AddAt(grMomentsScaledMultNpart[i],i);
  }
  
  statResults->GetFile()->cd();
  momentArr->Write("Moments",TObject::kSingleKey);
  momentNpartArr->Write("MomentsNpart",TObject::kSingleKey);
  momentArrSystematic->Write("momentSystematic",TObject::kSingleKey);
  momentArrSystematicNpart->Write("momentSystematicNpart",TObject::kSingleKey);
  momentScaledArr->Write("MomentsScaled",TObject::kSingleKey);
  momentScaledArrNpart->Write("MomentsScaledNpart",TObject::kSingleKey);
  subsamArr->Write("SubSamples",TObject::kSingleKey);

  fMoments->Close();
  delete statResults;
}
// =====================================================================================================
void CompareToReference(TString refFile, TString xFile, TString nuName)
{
  
  //
  // Compare refeence to calculation for the systematic uncertainty
  // 
                                                                                       
  // To RUN
  /*                                        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_cRows/TIdenResults/TIdenResults_EtaMomScan_minus1ksFixed_GAUSS
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_cRows/TIdenResults/TIdenResults_EtaMomScan_minus1ksFree_GAUSS

  aliroot -l 
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
  CompareToReference("PlotsMomentsStatErr_MC.root","PlotsMomentsStatErr_iter3.root","pikaNuDyn")
  */
  //
 
  TFile *f = TFile::Open(refFile);
  TFile *g = TFile::Open(xFile);
  TObjArray *refArr = (TObjArray*)f->Get("Moments_rec");
  TObjArray *xArr   = (TObjArray*)g->Get("Moments");
  
  
  TString nuName1; 
  if (nuName.Contains("pikaNuDyn")) nuName1="pikaNuDyn";
  if (nuName.Contains("piprNuDyn")) nuName1="prpiNuDyn";
  if (nuName.Contains("kaprNuDyn")) nuName1="prkaNuDyn";
  
  TGraphErrors *refgr = (TGraphErrors*)refArr->FindObject(nuName);
  TGraphErrors *xgr   = (TGraphErrors*)xArr  ->FindObject(nuName1);
  refgr->SetMarkerColor(kRed+1);    
  refgr->SetLineColor(kRed+1);            
  refgr->SetMarkerStyle(20); 
  
  Int_t N = refgr->GetN();
  
  Double_t *x1=new Double_t[N];
  Double_t *x2=new Double_t[N];
  Double_t *y1=new Double_t[N];
  Double_t *y2=new Double_t[N];
  Double_t *ex1=new Double_t[N];
  Double_t *ex2=new Double_t[N];
  Double_t *ey1=new Double_t[N];
  Double_t *ey2=new Double_t[N];
  Double_t *rx=new Double_t[N];
  Double_t *ry=new Double_t[N];
  Double_t *erx=new Double_t[N];
  Double_t *ery=new Double_t[N];

  Double_t *xdiff=new Double_t[N];
  Double_t *ydiff=new Double_t[N];
  
  for (Int_t i=0; i<N; i++){
    
//     cout << refgr->GetX()[N-1-i] << endl;
    x1[i] = centArr2[i];  xdiff[i] = centArr2[i]; x2[i] = centArr2[i]; rx[i] = centArr[i];  ex1[i] = 0; ex2[i] = 0; erx[i] = 0; 
    
    y1[i] = xgr  ->GetY()[N-1-i];
    y2[i] = refgr->GetY()[i];
    ry[i] = y2[i]/y1[i];
    ydiff[i] = TMath::Abs(y2[i]-y1[i]);
    
    ey1[i] = xgr  ->GetEY()[N-1-i];
    ey2[i] = refgr->GetEY()[i];
    
    ery[i] = TMath::Sqrt( (ey1[i]*ey1[i])/(y2[i]*y2[i]) + (y1[i]*y1[i]*ey2[i]*ey2[i])/(y2[i]*y2[i]*y2[i]*y2[i]) );

  }
   
  TGraphErrors *gr1     = new TGraphErrors(N,x1,y1,ex1,ey1);
//   TGraphErrors *gr2     = new TGraphErrors(N,x2,y2,ex2,ey2);
  TGraphErrors *grratio = new TGraphErrors(N,rx,ry,erx,ery);
  TGraphErrors *grdiff  = new TGraphErrors(N,xdiff,ydiff,0,0);

  
  grratio->SetMarkerStyle(20);      grratio->GetYaxis()->SetRangeUser(0,2); grratio->SetMarkerColor(kBlack);
  grdiff->SetMarkerStyle(20);       grdiff->SetMarkerColor(kBlack);
  gr1->SetMarkerColor(kGreen+2);    gr1->SetLineColor(kGreen+2);            gr1->SetMarkerStyle(21);                       
  grratio->GetXaxis()->SetTitle("Centrality (%)"); 
  grratio->GetYaxis()->SetTitle("(MC rec.)/(Id. Method)"); 

  grdiff->GetXaxis()->SetTitle("Centrality (%)"); 
  grdiff->GetYaxis()->SetTitle("(MC rec.)-(Id. Method)"); 

  
  // Draw Nudyn
  TCanvas *nu = new TCanvas("nu", "Alice Figure Template", 800, 600);   // Rectangular
  nu->cd();
  refgr->Draw("alp");
  gr1->Draw("lp");
  TLegend leg(0.20, 0.85, 0.8, 0.7);
  leg.SetTextFont(62); leg.SetTextSize(0.045); leg.SetFillColor(0); leg.SetBorderSize(0);
  leg.AddEntry(refgr," MC reconstructed","LPE");
  leg.AddEntry(gr1, " Identity Method","LPE");
  leg.Draw("same");
  
  // Draw ratio
  TCanvas *can = new TCanvas("can", "Alice Figure Template", 800, 600);   // Rectangular
  can->cd();
  can->SetGridy();
  TLine* line = new TLine(0., 1.,80., 1.);
  line -> SetLineColor(1); line -> SetLineWidth(4); line -> SetLineStyle(2);
  grratio->Draw("alp");
  line->Draw("same");
  
  if(savePDF) {
    nu ->SaveAs(Form("nudiff_%s.pdf",nuName.Data()));
    can->SaveAs(Form("nuratio_%s.pdf",nuName.Data()));
  }
  TFile *fOut = TFile::Open(Form("%s.root",nuName.Data()),"recreate");
  nu     ->Write(Form("nudiff_%s",nuName.Data()));
  can    ->Write(Form("nuratio_%s",nuName.Data()));
  grdiff ->Write(Form("grdiff%s",nuName.Data()));
  grratio->Write(Form("grratio%s",nuName.Data()));
  refgr  ->Write(Form("refgr%s",nuName.Data()));
  gr1    ->Write(Form("gr1%s",nuName.Data()));

  fOut->Close();

   
}
// =====================================================================================================
void CompareToReferenceSyst(TString refFile, TString xFile, Int_t iterREF, Int_t iterX, TString outName)
{
  
  //
  // Compare refeence to calculation for the systematic uncertainty
  // 
                                                                                       
  // To RUN
  /*                                        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics
  aliroot -l 
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
  TString ref="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts_PxPyPz//TIdenResults/TIdenResults_EtaMomScan_2015_09_10_00_10/Eta_-0.8_0.8_Mom_0.2_1.5/PlotsMomentsStatErr_iter6.root"
  TString x  ="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts_PxPyPz//TIdenResults/TIdenResults_EtaMomScan_2015_09_10_00_10/Eta_-0.8_0.8_Mom_0.2_1.5/PlotsMomentsStatErr_iter6.root"
  CompareToReferenceSyst(ref,x,6)
  */
  //
 
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C .");

  // Output file
  TTreeSRedirector *systResults = new TTreeSRedirector(Form("Systematic_%s_iterREF%d_iterX_%d.root",outName.Data(),iterREF,iterX),"recreate");

  // Get the trees
  TFile *fRef = TFile::Open(refFile);  
  TFile *fX   = TFile::Open(xFile);  
  TTree *treeRef  = (TTree*)fRef -> Get("momTree"); 
  TTree *treeX    = (TTree*)fX   -> Get("momTree"); 
  TGraphErrors *nudynprpiREF=NULL;
  TGraphErrors *nudynprkaREF=NULL;
  TGraphErrors *nudynpikaREF=NULL;
  TGraphErrors *nudynprpiX=NULL;
  TGraphErrors *nudynprkaX=NULL;
  TGraphErrors *nudynpikaX=NULL;
 
  // GetList of branches for each tree
  TObjArray *branchM = (TObjArray*)treeRef -> GetListOfBranches();
  Int_t nArrEntries  = branchM->GetEntriesFast();
  cout << " nArrEntries = " << nArrEntries << endl;

  // Create all "moment vs cent" TGraphErrors
  for (Int_t i=0; i<nArrEntries; ++i){
        
    // Get the brach name 
    TString branchName = branchM->At(i)->GetName();
    TString drawMoment = branchM->At(i)->GetName();
    TString drawMoment1 = Form("%s:cent",drawMoment.Data());
    TCut allStatREF = Form("subsample==0 && fitIter==%d && etaDown<-0.7",iterREF); 
    TCut allStatX   = Form("subsample==0 && fitIter==%d && etaDown<-0.7",iterX);       
    if ( !(branchName.Contains("prpiNuDyn") || branchName.Contains("prkaNuDyn") || branchName.Contains("pikaNuDyn")) ) continue;
    
    // create the graph object
    treeRef->Draw(drawMoment1,allStatREF,"goff"); 
    treeX->Draw(drawMoment1,allStatX,"goff"); 
    Int_t N = treeRef->GetSelectedRows();
    
    if (branchName.Contains("prpiNuDyn")) {
      nudynprpiREF=new TGraphErrors(N,treeRef->GetV2(),treeRef->GetV1());   nudynprpiREF->GetXaxis()->SetTitle("Centrality (%)"); 
      nudynprpiREF->GetYaxis()->SetTitle(branchName);                       nudynprpiREF->SetMarkerStyle(20);
    } 
    if (branchName.Contains("prkaNuDyn")) {
      nudynprkaREF=new TGraphErrors(N,treeRef->GetV2(),treeRef->GetV1());   nudynprkaREF->GetXaxis()->SetTitle("Centrality (%)"); 
      nudynprkaREF->GetYaxis()->SetTitle(branchName);                       nudynprkaREF->SetMarkerStyle(20);
    } 
    if (branchName.Contains("pikaNuDyn")) {
      nudynpikaREF=new TGraphErrors(N,treeRef->GetV2(),treeRef->GetV1());   nudynpikaREF->GetXaxis()->SetTitle("Centrality (%)"); 
      nudynpikaREF->GetYaxis()->SetTitle(branchName);                       nudynpikaREF->SetMarkerStyle(20);
    }
    if (branchName.Contains("prpiNuDyn")) {
      nudynprpiX=new TGraphErrors(N,treeX->GetV2(),treeX->GetV1());  nudynprpiX->GetXaxis()->SetTitle("Centrality (%)"); 
      nudynprpiX->GetYaxis()->SetTitle(branchName);                  nudynprpiX->SetMarkerStyle(20);
    }
    if (branchName.Contains("prkaNuDyn")) {
      nudynprkaX=new TGraphErrors(N,treeX->GetV2(),treeX->GetV1());  nudynprkaX->GetXaxis()->SetTitle("Centrality (%)"); 
      nudynprkaX->GetYaxis()->SetTitle(branchName);                  nudynprkaX->SetMarkerStyle(20);
    } 
    if (branchName.Contains("pikaNuDyn")) {
      nudynpikaX=new TGraphErrors(N,treeX->GetV2(),treeX->GetV1());  nudynpikaX->GetXaxis()->SetTitle("Centrality (%)"); 
      nudynpikaX->GetYaxis()->SetTitle(branchName);                  nudynpikaX->SetMarkerStyle(20);
    }
    
    // Add graphs to array and write in the file
    cout << " nudyn to be processed = " << branchName << endl;
    
  }
  TGraphErrors *nudynprpiRatio = RatioGraph(nudynprpiREF,nudynprpiX);
  TGraphErrors *nudynprkaRatio = RatioGraph(nudynprkaREF,nudynprkaX);
  TGraphErrors *nudynpikaRatio = RatioGraph(nudynpikaREF,nudynpikaX);

  TGraphErrors *nudynprpiDiff = DiffGraph(nudynprpiREF,nudynprpiX);
  TGraphErrors *nudynprkaDiff = DiffGraph(nudynprkaREF,nudynprkaX);
  TGraphErrors *nudynpikaDiff = DiffGraph(nudynpikaREF,nudynpikaX);
  
  nudynprpiRatio->GetXaxis()->SetTitle("Centrality (%)"); nudynprpiRatio->GetYaxis()->SetTitle("nudynprpiRatio"); nudynprpiRatio->SetMarkerStyle(20);
  nudynprkaRatio->GetXaxis()->SetTitle("Centrality (%)"); nudynprkaRatio->GetYaxis()->SetTitle("nudynprkaRatio"); nudynprkaRatio->SetMarkerStyle(20);
  nudynpikaRatio->GetXaxis()->SetTitle("Centrality (%)"); nudynpikaRatio->GetYaxis()->SetTitle("nudynpikaRatio"); nudynpikaRatio->SetMarkerStyle(20);
  nudynprpiDiff->GetXaxis()->SetTitle("Centrality (%)"); nudynprpiDiff->GetYaxis()->SetTitle("nudynprpiDiff"); nudynprpiDiff->SetMarkerStyle(20);
  nudynprkaDiff->GetXaxis()->SetTitle("Centrality (%)"); nudynprkaDiff->GetYaxis()->SetTitle("nudynprkaDiff"); nudynprkaDiff->SetMarkerStyle(20);
  nudynpikaDiff->GetXaxis()->SetTitle("Centrality (%)"); nudynpikaDiff->GetYaxis()->SetTitle("nudynpikaDiff"); nudynpikaDiff->SetMarkerStyle(20);

  systResults->GetFile()->cd();
  nudynprpiRatio->Write("nudynprpiRatio");
  nudynprkaRatio->Write("nudynprkaRatio");
  nudynpikaRatio->Write("nudynpikaRatio");
  nudynprpiDiff->Write("nudynprpiDiff");
  nudynprkaDiff->Write("nudynprkaDiff");
  nudynpikaDiff->Write("nudynpikaDiff");
  
  nudynprpiREF->Write("nudynprpiREF");
  nudynprkaREF->Write("nudynprkaREF");
  nudynpikaREF->Write("nudynpikaREF");
  nudynprpiX->Write("nudynprpiX");
  nudynprkaX->Write("nudynprkaX");
  nudynpikaX->Write("nudynpikaX");

  
  delete systResults;
   
}
// =====================================================================================================
void FinalSyst()
{
  
  //
  // Compare refeence to calculation for the systematic uncertainty
  // 
                                                                                       
  // To RUN
  /*                                        
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics
  aliroot -l 
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
  FinalSyst()
  */
  //
 
  TFile *file[11];
  TGraphErrors *prpi[11]={0};
  TGraphErrors *prka[11]={0};
  TGraphErrors *pika[11]={0};
  Double_t xprpi[9];
  Double_t yprpi[9];
  Double_t xprka[9];
  Double_t yprka[9];
  Double_t xpika[9];
  Double_t ypika[9];
  
  file[0]    = TFile::Open("Systematic_cRows_60_iter6.root");  
  file[1]    = TFile::Open("Systematic_cRows_100_iter6.root");  
  file[2]    = TFile::Open("Systematic_CL1_iter6.root");  
  file[3]    = TFile::Open("Systematic_smallDCAxy_iter6.root");  
  file[4]    = TFile::Open("Systematic_largeDCAxy_iter6.root");  
  file[5]    = TFile::Open("Systematic_chi2ndf5_iter6.root");  
  file[6]    = TFile::Open("Systematic_chi2ndf3_iter6.root");  
  file[7]    = TFile::Open("Systematic__iterREF5_iterX_7.root");  
  file[8]    = TFile::Open("Systematic__iterREF5_iterX_6.root");  
  file[9]    = TFile::Open("Systematic__iterREF5_iterX_4.root");  
  file[10]    = TFile::Open("Systematic__iterREF5_iterX_3.root");  

  for (Int_t i=0;i<11;i++) {
    prpi[i] = (TGraphErrors*)file[i]->Get("nudynprpiDiff");
    prka[i] = (TGraphErrors*)file[i]->Get("nudynprkaDiff");
    pika[i] = (TGraphErrors*)file[i]->Get("nudynpikaDiff");
  }
  
  Float_t systfactor=1.;
  for (Int_t k=0;k<9;k++){
    Double_t sumprpi=0.;
    Double_t sumprka=0.;
    Double_t sumpika=0.;
    for (Int_t i=0;i<11;i++) {
      if (i==10 || i==7 ||i==8) systfactor=10000000.;
      if (i==9) systfactor=1.;
      xprpi[k]=k; xprka[k]=k; xpika[k]=k;
      if (prpi[i]->GetY()[k]>-5 && prpi[i]->GetY()[k]<5 ) sumprpi+=(prpi[i]->GetY()[k]/systfactor)*(prpi[i]->GetY()[k]/systfactor);
      if (prka[i]->GetY()[k]>-5 && prka[i]->GetY()[k]<5 ) sumprka+=(prka[i]->GetY()[k]/systfactor)*(prka[i]->GetY()[k]/systfactor);
      if (pika[i]->GetY()[k]>-5 && pika[i]->GetY()[k]<5 ) sumpika+=(pika[i]->GetY()[k]/systfactor)*(pika[i]->GetY()[k]/systfactor);
      yprpi[k]=TMath::Sqrt(sumprpi);
      yprka[k]=TMath::Sqrt(sumprka);
      ypika[k]=TMath::Sqrt(sumpika);
      
    }
  }

  TGraphErrors *grfinalprpi = new TGraphErrors(9,xprpi,yprpi,0,0);
  TGraphErrors *grfinalprka = new TGraphErrors(9,xprka,yprka,0,0);
  TGraphErrors *grfinalpika = new TGraphErrors(9,xpika,ypika,0,0); 

  TFile *fOut = TFile::Open("FinalSyst.root","recreate");
  grfinalprpi->Write("prpiFinalSyst");
  grfinalprka->Write("prkaFinalSyst");
  grfinalpika->Write("pikaFinalSyst");
  fOut->Close();
   
}
// =====================================================================================================
void CompareDataVsMC(TString dataFile, TString hijingFile, TString amptFile)
{
  
  //
  // Compare data to montecarlo data t
  // 
                                                                                       
  // To RUN
  /*                                        
   aliroot -l 
   TString data = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/TIdenResults/TIdenResults_EtaMomScan_le30centbinsDownscaled_EtaScan_GOOD/Eta_-0.8_0.8_Mom_0.2_1.5/PlotsMomentsStatErr_iter6.root"
   TString ampt="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_AMPT_TightCuts/EfficiencyCheck/Ratios_CentScan_eta_-0.8-0.8_mom_0.2-1.5.root"
   TString hijing="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_HIJING_TightCuts/EfficiencyCheck/Ratios_CentScan_eta_-0.8-0.8_mom_0.2-1.5.root"

   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
    CompareDataVsMC(data,hijing,ampt)

  */
  //
  TGaxis::SetMaxDigits(3);
  gStyle->SetEndErrorSize(3);
  gStyle->SetErrorX(3);
 
  TGraphErrors *grHijing[3];
  TGraphErrors *grAmpt[3];
  TGraphErrors *grData[3];
  TGraphAsymmErrors *grDataSyst[3];
  TFile *finalSyst = TFile::Open("/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics/FinalSyst.root"); 
  TGraphErrors * prpiSystALICE = (TGraphErrors*)finalSyst->Get("prpiFinalSyst");
  TGraphErrors * prkaSystALICE = (TGraphErrors*)finalSyst->Get("prkaFinalSyst");
  TGraphErrors * pikaSystALICE = (TGraphErrors*)finalSyst->Get("pikaFinalSyst");

  TFile *fdata = TFile::Open(dataFile);
  TObjArray *dataArr = (TObjArray*)fdata->Get("Moments"); 
  grData[0] = (TGraphErrors*)dataArr->FindObject("pikaNuDyn");
  grData[1] = (TGraphErrors*)dataArr->FindObject("prpiNuDyn");
  grData[2] = (TGraphErrors*)dataArr->FindObject("prkaNuDyn");
  
  TObjArray *dataArrSyst = (TObjArray*)fdata->Get("momentSystematic"); 
  grDataSyst[0] = (TGraphAsymmErrors*)dataArrSyst->FindObject("pikaNuDyn");
  grDataSyst[1] = (TGraphAsymmErrors*)dataArrSyst->FindObject("prpiNuDyn");
  grDataSyst[2] = (TGraphAsymmErrors*)dataArrSyst->FindObject("prkaNuDyn");
  
  // Common legends
  TLegend legErr(0.3, 0.55, 0.45, 0.62);  legErr.SetBorderSize(0);
  legErr.SetTextFont(62);
  legErr.SetTextSize(0.035); //legErr.SetTextFont(62);
  legErr.SetFillColor(0);
  legErr.SetFillStyle(4000);
  legErr.AddEntry(grDataSyst[0],  " Systematic uncertainty"   ,"f");
  legErr.AddEntry(grData[0],      " Statistical uncertainty"  ,"l");
 
  TLegend legLogo(0.35, 0.45, 0.70, 0.55);  legLogo.SetBorderSize(0);
  legLogo.SetTextFont(62);
  legLogo.SetTextSize(0.035); //legLogo.SetTextFont(62);
  legLogo.SetFillColor(0);
  legLogo.SetFillStyle(4000);
  legLogo.AddEntry((TObject*)0, "ALICE Preliminary, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV"," ");
  legLogo.AddEntry((TObject*)0, "   |#eta|<0.8, 0.2<#it{p}<1.5 GeV/#it{c} ", " ");
 
  
  // error propagation
  for (Int_t i=0;i<9;i++){
    
    Double_t sys0 = TMath::Abs(pikaSystALICE->GetY()[8-i]);
    Double_t mc0  = TMath::Abs(grDataSyst[0]->GetY()[i])*0.2;
    Double_t sys1 = TMath::Abs(prpiSystALICE->GetY()[8-i]);
    Double_t mc1  = TMath::Abs(grDataSyst[1]->GetY()[i])*0.2;
    Double_t sys2 = TMath::Abs(prkaSystALICE->GetY()[8-i]);
    Double_t mc2  = TMath::Abs(grDataSyst[2]->GetY()[i])*0.2;
    
    grDataSyst[0]->SetPointEYhigh(i,TMath::Sqrt(sys0*sys0+mc0*mc0));
    grDataSyst[0]->SetPointEYlow(i ,TMath::Sqrt(sys0*sys0+mc0*mc0));
    grDataSyst[1]->SetPointEYhigh(i,TMath::Sqrt(sys1*sys1+mc1*mc1));
    grDataSyst[1]->SetPointEYlow(i ,TMath::Sqrt(sys1*sys1+mc1*mc1));
    grDataSyst[2]->SetPointEYhigh(i,TMath::Sqrt(sys2*sys2+mc2*mc2));
    grDataSyst[2]->SetPointEYlow(i ,TMath::Sqrt(sys2*sys2+mc2*mc2)); 
   
  }
  
  TFile *fhijing = TFile::Open(hijingFile);
  TObjArray *hijingArr = (TObjArray*)fhijing->Get("Moments_gen"); 
  grHijing[0] = (TGraphErrors*)hijingArr->FindObject("pikaNuDyn");
  grHijing[1] = (TGraphErrors*)hijingArr->FindObject("piprNuDyn");
  grHijing[2] = (TGraphErrors*)hijingArr->FindObject("kaprNuDyn");
  
  TFile *fampt = TFile::Open(amptFile);
  TObjArray *amptArr = (TObjArray*)fampt->Get("Moments_gen"); 
  grAmpt[0] = (TGraphErrors*)amptArr->FindObject("pikaNuDyn");
  grAmpt[1] = (TGraphErrors*)amptArr->FindObject("piprNuDyn");
  grAmpt[2] = (TGraphErrors*)amptArr->FindObject("kaprNuDyn");
  
  for (Int_t i=0; i<3; i++){
    grData[i]     ->SetMarkerStyle(20);   grData[i]     ->SetMarkerColor(kRed+1);    grData[i]     ->SetLineColor(kRed+1); 
    grDataSyst[i] ->SetMarkerStyle(20);   grDataSyst[i] ->SetMarkerColor(kRed+1);    grDataSyst[i] ->SetLineColor(kRed+1);  
    grAmpt[i]     ->SetMarkerStyle(21);   grAmpt[i]     ->SetMarkerColor(kGreen+2);  grAmpt[i]     ->SetLineColor(kGreen+2);  
    grHijing[i]   ->SetMarkerStyle(34);   grHijing[i]   ->SetMarkerColor(kBlack);    grHijing[i]   ->SetLineColor(kBlack); 
     
    grData[i]  ->SetMarkerSize(1.5);   grData[i]  ->SetLineWidth(2.); grData[i]->GetYaxis()->SetTitleOffset(1.3);
    grHijing[i]->SetMarkerSize(2.5);    grHijing[i]->SetLineWidth(2.); grHijing[i]->GetYaxis()->SetTitleOffset(1.3);
    grAmpt[i]  ->SetMarkerSize(1.7);   grAmpt[i]  ->SetLineWidth(2.); grAmpt[i]->GetYaxis()->SetTitleOffset(1.3);

    grData[i]->GetXaxis()->SetTitle("Centrality (%)"); 
    grHijing[i]->GetXaxis()->SetTitle("Centrality (%)"); 
    grAmpt[i]->GetXaxis()->SetTitle("Centrality (%)"); 
    grDataSyst[i]->SetFillStyle(0); 
    grData[i]->SetFillStyle(0);      
    
    grDataSyst[i]->GetXaxis()->SetTitle("Centrality (%)"); 
    grDataSyst[i]->GetYaxis()->SetTitleOffset(1.3);
    
    grDataSyst[i]->SetLineWidth(2); 
    grDataSyst[i]->SetLineWidth(2); 
    
    grData[i]     ->GetYaxis()->SetTitleOffset(1.);
    grDataSyst[i] ->GetYaxis()->SetTitleOffset(1.);   
    grAmpt[i]     ->GetYaxis()->SetTitleOffset(1.);  
    grHijing[i]   ->GetYaxis()->SetTitleOffset(1.);
   
    
    cout << TMath::MinElement(grData[i]->GetN(),grData[i]->GetY())       << "   " << TMath::MaxElement(grData[i]->GetN(),grData[i]->GetY()) << endl;
    cout << TMath::MinElement(grHijing[i]->GetN(),grHijing[i]->GetY())   << "   " << TMath::MaxElement(grHijing[i]->GetN(),grHijing[i]->GetY()) << endl;
    cout << TMath::MinElement(grAmpt[i]->GetN(),grAmpt[i]->GetY())       << "   " << TMath::MaxElement(grAmpt[i]->GetN(),grAmpt[i]->GetY()) << endl;

  }
  grData[0]     ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
  grData[1]     ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
  grData[2]     ->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
  
  grDataSyst[0] ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
  grDataSyst[1] ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
  grDataSyst[2] ->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
  
  grHijing[0]   ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
  grHijing[1]   ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
  grHijing[2]   ->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
  
  grAmpt[0]     ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
  grAmpt[1]     ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
  grAmpt[2]     ->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
  
  // Legend for model comparison
  TLegend leg(0.23, 0.67, 0.55, 0.86);
  leg.SetTextFont(62); leg.SetTextSize(0.035); leg.SetFillColor(0); leg.SetBorderSize(0); leg.SetFillStyle(4000);
  leg.AddEntry(grData[0]  , "ALICE Data, stat. errors ","LPE");
  leg.AddEntry(grDataSyst[0],  " Systematic uncertainty"   ,"f");
  leg.AddEntry(grHijing[0], "HIJING ","PE");
  leg.AddEntry(grAmpt[0]  , "AMPT ","PE");
  
  // For Centrality dependence plot
  TLegend legCent(0.25, 0.75, 0.45, 0.85);
  legCent.SetTextFont(62); legCent.SetTextSize(0.035); legCent.SetFillColor(0); legCent.SetBorderSize(0); legCent.SetFillStyle(4000);
  legCent.AddEntry(grData[0], "ALICE Data, stat. errors: ", "LPE");
  legCent.AddEntry(grDataSyst[0],  " Systematic uncertainty"   ,"f");
//   legCent.AddEntry(grData[0],      " Statistical uncertainty"  ,"l");
  
  Int_t nPoints = grData[0]->GetN();
  // pika comparison
  Int_t nMaxX = 82;
  TLine* line0 = new TLine(0., 0.,nMaxX, 0.); line0 -> SetLineStyle(2); line0 -> SetLineWidth(3);  
  TCanvas *pikaDiff = new TCanvas("pikaDiff", "pika Data vs MC comparison", 900, 600); 
  pikaDiff->cd(); 
  if (transparent){
    pikaDiff->SetFillStyle(4000);
    pikaDiff->SetFrameFillStyle(4000);
  }
  grHijing[0]  ->GetYaxis()->SetRangeUser(-0.004,TMath::MaxElement(nPoints,grHijing[0]->GetY())*1.2); 
  grHijing[0]  ->Draw("ap");
  grData[0]    ->Draw("p");
  grAmpt[0]    ->Draw("p");
  grDataSyst[0]->Draw("p5");
  line0->Draw("same");
  leg.DrawClone("same");
  legLogo.DrawClone("same");

  // pika comparison
  TCanvas *piprDiff = new TCanvas("piprDiff", "pipr Data vs MC comparison", 800, 600); 
  piprDiff->cd();
  if (transparent){
    piprDiff->SetFillStyle(4000);
    piprDiff->SetFrameFillStyle(4000);
  }
  grHijing[1]  ->GetYaxis()->SetRangeUser(-0.03,TMath::MaxElement(nPoints,grHijing[1]->GetY())*1.2); 
  grHijing[1]  ->Draw("ap");
  grData[1]    ->Draw("p");
  grAmpt[1]    ->Draw("p");
  grDataSyst[1]->Draw("p5");
  line0->Draw("same");
  leg.DrawClone("same");
  legLogo.DrawClone("same");
  
  // pika comparison
  TCanvas *kaprDiff = new TCanvas("kaprDiff", "kapr Data vs MC comparison", 800, 600);
  kaprDiff->cd(); 
  if (transparent){
    kaprDiff->SetFillStyle(4000);
    kaprDiff->SetFrameFillStyle(4000);
  }
  grHijing[2]  ->GetYaxis()->SetRangeUser(-0.01,TMath::MaxElement(nPoints,grHijing[2]->GetY())*1.2); 
  grHijing[2]  ->Draw("ap");
  grData[2]    ->Draw("p");
  grAmpt[2]    ->Draw("p");
  grDataSyst[2]->Draw("p5");
  line0->Draw("same");
  leg.DrawClone("same");
  legLogo.DrawClone("same");

 
  for (Int_t i=0; i<3; i++){ grData[i]->SetMarkerStyle(20); grData[i]->SetMarkerColor(kRed+1); grData[i]->SetLineColor(kRed+1); }
  
  // Centrality dependence plots
  TCanvas *pikaData = new TCanvas("pikaData", "pika Data", 800, 600); pikaData->cd(); 
  if (transparent){
    pikaData->SetFillStyle(4000);
    pikaData->SetFrameFillStyle(4000);
  }
  grDataSyst[0]->GetYaxis()->SetRangeUser(0.,0.045);
  grDataSyst[0]->Draw("ap5"); 
  grData[0]    ->Draw("p");
  legCent.DrawClone("same");
  legLogo.DrawClone("same");

  
  // pika comparison
  TCanvas *piprData = new TCanvas("piprData", "pipr Data", 800, 600); piprData->cd(); 
  if (transparent){
    piprData->SetFillStyle(4000);
    piprData->SetFrameFillStyle(4000);
  }
  grDataSyst[1]->GetYaxis()->SetRangeUser(-0.025,0.003); 
  grDataSyst[1]->Draw("ap5"); 
  grData[1]    ->Draw("p");
  legCent.DrawClone("same");
  legLogo.DrawClone("same");
  
  // kapr comparison
  TCanvas *kaprData = new TCanvas("kaprData", "kapr Data", 800, 600); kaprData->cd(); 
  if (transparent){
    kaprData->SetFillStyle(4000);
    kaprData->SetFrameFillStyle(4000);
  }
  grDataSyst[2]->GetYaxis()->SetRangeUser(-0.001,0.025); 
  grDataSyst[2]->Draw("ap5"); 
  grData[2]    ->Draw("p");
  legCent.DrawClone("same");
  legLogo.DrawClone("same");

  
  TFile *fOut = TFile::Open("FinalErrors.root","recreate");
  grDataSyst[0]->Write("pikasyst");
  grDataSyst[1]->Write("prpisyst");
  grDataSyst[2]->Write("prkasyst");
  grData[0]->Write("pikastat");
  grData[1]->Write("prpistat");
  grData[2]->Write("prkastat");
  fOut->Close();
  
}
// =====================================================================================================
void CompareIdentityVsTraditional(TString dataFile, TString tradFile)
{
  
  //
  // Compare data to montecarlo data t
  // 
                                                                                       
  // To RUN
  /*                                        
  aliroot -l 
  TString data = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/TIdenResults/TIdenResults_EtaMomScan_le30centbinsDownscaled_EtaScan_GOOD/Eta_-0.5_0.5_Mom_0.5_1.5/PlotsMomentsStatErr_iter6.root"
  TString trad="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/traditional/deepika_NuDynInNpart.root"
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
  CompareIdentityVsTraditional(data,trad)
  
  */
  //
  
  const Int_t nALICEtra = 8;
  
  TGraphErrors *grData[3];
  TGraphAsymmErrors *grDataSyst[3];
  TGraphErrors *grDataTra[3];
  TGraphAsymmErrors *grDataSystTra[3];

  
  TFile *fdata = TFile::Open(dataFile);
  TObjArray *dataArr = (TObjArray*)fdata->Get("MomentsNpart"); 
  grData[0] = (TGraphErrors*)dataArr->FindObject("pikaNuDyn");
  grData[1] = (TGraphErrors*)dataArr->FindObject("prpiNuDyn");
  grData[2] = (TGraphErrors*)dataArr->FindObject("prkaNuDyn");
  
  TObjArray *dataArrSyst = (TObjArray*)fdata->Get("momentSystematicNpart"); 
  grDataSyst[0] = (TGraphAsymmErrors*)dataArrSyst->FindObject("pikaNuDyn");
  grDataSyst[1] = (TGraphAsymmErrors*)dataArrSyst->FindObject("prpiNuDyn");
  grDataSyst[2] = (TGraphAsymmErrors*)dataArrSyst->FindObject("prkaNuDyn");
  
  TFile *fdataTrad = TFile::Open(tradFile);
  grDataTra[0] = (TGraphErrors*)fdataTrad->Get("gr_stat_bwc_P0515_UC_NkaNpi_NuDyn_etabin4");
  grDataTra[1] = (TGraphErrors*)fdataTrad->Get("gr_stat_bwc_P0515_UC_NprNpi_NuDyn_etabin4");
  grDataTra[2] = (TGraphErrors*)fdataTrad->Get("gr_stat_bwc_P0515_UC_NprNka_NuDyn_etabin4");
  grDataSystTra[0] = (TGraphAsymmErrors*)fdataTrad->Get("gr_syst_bwc_P0515_UC_NkaNpi_NuDyn_etabin4");
  grDataSystTra[1] = (TGraphAsymmErrors*)fdataTrad->Get("gr_syst_bwc_P0515_UC_NprNpi_NuDyn_etabin4");
  grDataSystTra[2] = (TGraphAsymmErrors*)fdataTrad->Get("gr_syst_bwc_P0515_UC_NprNka_NuDyn_etabin4");
  
  // some make up 
  for (Int_t i=0; i<3; i++){
    grData[i]      ->SetMarkerStyle(20);   grData[i]     ->SetMarkerColor(kBlack);     grData[i]     ->SetLineColor(kBlack); 
    grDataSyst[i]  ->SetMarkerStyle(20);   grDataSyst[i] ->SetMarkerColor(kBlack);     grDataSyst[i] ->SetLineColor(kBlack);  
    grDataTra[i]   ->SetMarkerStyle(21);   grDataTra[i]  ->SetMarkerColor(kRed+1);       grDataTra[i]  ->SetLineColor(kRed+1); 

    grData[i]     ->SetMarkerSize(1.5);
    grDataTra[i]  ->SetMarkerSize(1.5);

    grData[i]->GetXaxis()->SetTitle("<N_{part}>"); 
    grData[i]->GetYaxis()->SetTitleOffset(1.6);
    grData[i]->GetXaxis()->SetTitleSize(0.05); grData[i]->GetXaxis()->SetLabelSize(0.045); 
    grData[i]->GetYaxis()->SetTitleSize(0.05); grData[i]->GetYaxis()->SetLabelSize(0.045); 
    
    grDataSyst[i]->GetXaxis()->SetTitle("<N_{part}>"); 
    grDataSyst[i]->GetYaxis()->SetTitleOffset(1.6);
    grDataSyst[i]->GetXaxis()->SetTitleSize(0.05);         grDataSyst[i]->GetXaxis()->SetLabelSize(0.045); 
    grDataSyst[i]->GetYaxis()->SetTitleSize(0.05);         grDataSyst[i]->GetYaxis()->SetLabelSize(0.045); 
    
    grDataSystTra[i]->GetXaxis()->SetTitle("<N_{part}>"); 
    grDataSystTra[i]->GetYaxis()->SetTitleOffset(1.6);
    grDataSystTra[i]->GetXaxis()->SetTitleSize(0.05);         grDataSystTra[i]->GetXaxis()->SetLabelSize(0.045); 
    grDataSystTra[i]->GetYaxis()->SetTitleSize(0.05);         grDataSystTra[i]->GetYaxis()->SetLabelSize(0.045); 
    
    grDataTra[i]->GetXaxis()->SetTitle("<N_{part}>"); 
    grDataTra[i]->GetYaxis()->SetTitleOffset(1.6);
    grDataTra[i]->GetXaxis()->SetTitleSize(0.05);          grDataTra[i]->GetXaxis()->SetLabelSize(0.045); 
    grDataTra[i]->GetYaxis()->SetTitleSize(0.05);          grDataTra[i]->GetYaxis()->SetLabelSize(0.045); 
    
    grDataSystTra[i]->SetFillStyle(0); 
    grDataSyst[i]->SetFillStyle(0);      

  }
  grData[0]       ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
  grData[1]       ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
  grData[2]       ->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]");  
  grDataSyst[0]   ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
  grDataSyst[1]   ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
  grDataSyst[2]   ->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
  grDataTra[0]    ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
  grDataTra[1]    ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
  grDataTra[2]    ->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
  grDataSystTra[0]->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
  grDataSystTra[1]->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
  grDataSystTra[2]->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
  
 
  // prepare alegend for all plots
  TLegend leg(0.50, 0.85, 0.85, 0.7);
  leg.SetTextFont(62); leg.SetTextSize(0.045); leg.SetFillColor(0); leg.SetBorderSize(0);
  leg.AddEntry(grData[0],      " Identity Method ","LPE");
  leg.AddEntry(grDataTra[0]  , " Traditional Method ","LPE");

 
  // pika comparison
  TCanvas *pikaDiff = new TCanvas("pikaDiff", "pika TIdentity vs Traditional method comparison", 800, 600);   
  pikaDiff->cd();  pikaDiff->SetGridy();
  grDataTra[0] ->Draw("alp");      grDataTra[0]  ->GetYaxis()->SetRangeUser(0.,0.03); 
  grData[0]    ->Draw("lp5");
  grDataSystTra[0]->Draw("p3");
  grDataSyst[0]->Draw("p3");
  leg.Draw("same");
  
  // pika comparison
  TCanvas *piprDiff = new TCanvas("piprDiff", "pipr TIdentity vs Traditional method comparison", 800, 600);   
  piprDiff->cd(); piprDiff->SetGridy();
  grDataTra[1] ->Draw("alp");      grDataTra[1]  ->GetYaxis()->SetRangeUser(-0.004,0.013); 
  grData[1]    ->Draw("lp5");
  grDataSystTra[1]->Draw("p3");
  grDataSyst[1]->Draw("p3");
  leg.Draw("same");
  
  // pika comparison
  TCanvas *kaprDiff = new TCanvas("kaprDiff", "kapr TIdentity vs Traditional method comparison", 800, 600);   
  kaprDiff->cd(); kaprDiff->SetGridy();
  grDataTra[2] ->Draw("alp");      grDataTra[2]  ->GetYaxis()->SetRangeUser(-0.004,0.045); 
  grData[2]    ->Draw("lp5");
  grDataSystTra[2]->Draw("p3");
  grDataSyst[2]->Draw("p3");
  leg.Draw("same");
  
  pikaDiff->SaveAs("pika_IdenVsTra.pdf");
  piprDiff->SaveAs("pipr_IdenVsTra.pdf");
  kaprDiff->SaveAs("kapr_IdenVsTra.pdf");
   
}
// =====================================================================================================
Double_t CalculateStatErrors(TGraphErrors *grS)
{

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
  
  return TMath::Sqrt(sum/(n*(n-1)));
  
}
// =====================================================================================================
void AnalyseWs(TString wFile, TString wSliceFile)
{

  //
  // Calculate statistical errors using subsample method
  /*
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/IdDataTrees_300000_Results_OK2
  aliroot -l
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
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
void AnalyseHistsInFile(TString file)
{
  // 
  // delete TH1D object from file
  /*
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/IdDataTrees_MC_060514/
  aliroot -l 
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
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
void FitHistogram(TH1D *hClean)
{
 
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
    hClean->SetLineColor(kRed+1);
    asyGaus->SetLineColor(kRed+1);
  }
 
  hClean->GetXaxis()->SetRangeUser(TMath::Max(0.,mean-5.*rms),mean+5.*rms);
    
  // retrieve fit parameters
  hClean->Fit(asyGaus,"QMR");  
  delete asyGaus;
}
// =====================================================================================================
void GraphShade()
{
   TCanvas *c1 = new TCanvas("c1",
      "A Simple Graph Example",200,10,700,500);

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
void ErrorBandTH1()
{
  TH1F *h = new TH1F("h","",10,0,10);

  h->SetLineColor(kRed+1);

  for(int i=0; i<9; i++){
    h->Fill(i,i*0.5);
  }

  h->DrawCopy("hist"); 
  h->SetFillColor(kBlue);
  h->SetFillStyle(3018);
  h->Draw("e2same");
 
}
// =====================================================================================================
TCanvas *TwoScales(TTree *t)
{
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
  hwKa->SetMarkerColor(kRed+1);
  hwPr->SetMarkerColor(kBlue);


  TCanvas *c2 = new TCanvas("c2","hists with different scales",600,400);
  c2->SetGrid();
  
//   hdEdx->GetYaxis()->SetLabelColor(kMagenta);
  hdEdx->Draw(); 
  hdEdx->GetXaxis()->SetTitle("TPC dEdx Signal");
  hdEdx->GetYaxis()->SetTitle("identity variable w (%)");
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
  
  TLegend leg(0.7, 0.6, 0.9, 0.9);
//   leg.SetHeader("w : identity variable");
  leg.SetFillColor(0);
  leg.AddEntry(hwEl,"Electron","p");
  leg.AddEntry(hwPi,"Pion"    ,"p");
  leg.AddEntry(hwKa,"Kaon"    ,"p");
  leg.AddEntry(hwPr,"Proton"  ,"p");
  leg.Draw();
  
  return c2;
  
}
// =====================================================================================================
TCanvas *OverlayPads(TTree *t)
{
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
  hwKa->SetMarkerColor(kRed+1);

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
  TLegend leg(0.7, 0.6, 0.9, 0.9);
  leg.SetFillColor(0);
  leg.AddEntry(hwEl,"Electron","p");
  leg.AddEntry(hwPi,"Pion"    ,"p");
  leg.AddEntry(hwKa,"Kaon"    ,"p");
  leg.AddEntry(hwPr,"Proton"  ,"p");
  leg.Draw();
  
  return c1;
  
}
// =====================================================================================================
void SetYTitleName(TString branchName)
{
  
  //
  // Set the title name 
  //
  
  TString yTitleLatex;
  if (branchName=="pikaNuDyn") yTitleLatex="#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]";
  if (branchName=="prpiNuDyn") yTitleLatex="#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]";
  if (branchName=="prkaNuDyn") yTitleLatex="#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]";
  if (branchName=="pipiNuDyn") yTitleLatex="#nu_{dyn} [#pi^{+}+#pi^{-},#pi^{+}+#pi^{-}]";
  if (branchName=="kakaNuDyn") yTitleLatex="#nu_{dyn} [K^{+}+K^{-},K^{+}+K^{-}]";
  if (branchName=="prprNuDyn") yTitleLatex="#nu_{dyn} [p+#bar{p},p+#bar{p}]";
  if (branchName=="pika") yTitleLatex="#LT#pi,K#GT";
  if (branchName=="pipr") yTitleLatex="#LT#pi,p#GT";
  if (branchName=="kapr") yTitleLatex="#LTK,p#GT";
  if (branchName=="pipi") yTitleLatex="#LT#pi,#pi#GT";
  if (branchName=="kaka") yTitleLatex="#LTK,K#GT";
  if (branchName=="prpr") yTitleLatex="#LTp,pGT";
  if (branchName=="pi") yTitleLatex="#LT#pi#GT";
  if (branchName=="ka") yTitleLatex="#LTK#GT";
  if (branchName=="pr") yTitleLatex="#LTp#GT";
          
}
// =====================================================================================================
TGraphErrors * CorrectCentBinOfMoments(TGraphErrors *grCorr)
{

  //
  // correct the centrality axis of the moments
  //
  
  Int_t NN = grCorr->GetN();
  Double_t yCorr[NN], yECorr[NN];
  for (Int_t n=0;n<NN;n++){
    yCorr[n]  = grCorr->GetY()[NN-1-n];
    yECorr[n] = grCorr->GetEY()[NN-1-n];
  }
  TGraphErrors *gr = new TGraphErrors(NN,centArr,yCorr,0,yECorr);
  gr->SetDrawOption("a3");
  gr->SetName(grCorr->GetName());
  gr->GetXaxis()->SetTitle("Centrality (%)");
  gr->GetYaxis()->SetTitle(grCorr->GetName());
  gr->SetFillColor(kBlack);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlack);
  return gr;
  
} 
// =====================================================================================================
TGraphErrors * MakeNpartMoments(TGraphErrors *grCorr)
{
  
  //
  // correct the centrality axis of the moments
  //
  
  Int_t NN = grCorr->GetN();
  Double_t yCorr[NN], yECorr[NN];
  for (Int_t n=0;n<NN;n++){
    yCorr[n]  = grCorr->GetY()[NN-1-n];
    yECorr[n] = grCorr->GetEY()[NN-1-n];
  }
  TGraphErrors *gr = new TGraphErrors(NN,nPartArr,yCorr,0,yECorr);
  gr->SetDrawOption("a3");
  gr->SetName(grCorr->GetName());
  gr->GetXaxis()->SetTitle("N_{part}");
  gr->GetYaxis()->SetTitle(grCorr->GetName());
  gr->SetFillColor(kBlack);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlack);
  return gr;
  
} 
// =====================================================================================================
void MultiplydNchStatAndSys(TString grStatName, TString grSystName)
{
  
  //
  // correct the centrality axis of the moments
  //
  /*
   cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics
   aliroot -l 
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
   MultiplydNchStatAndSys("prpistat","prpisyst")
   MultiplydNchStatAndSys("prkastat","prkasyst")
   MultiplydNchStatAndSys("pikastat","pikasyst")

    
   */
  TFile *finalSyst = TFile::Open("/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics/FinalErrors.root"); 
  TGraphAsymmErrors * grSyst = (TGraphAsymmErrors*)finalSyst->Get(grSystName);
  TGraphErrors      * grStat = (TGraphErrors*)finalSyst->Get(grStatName);
  Int_t NN = grStat->GetN(); 
  
  TTreeSRedirector *scaling = new TTreeSRedirector(Form("Nudyn_dnch_Scaling_%s.root",grSystName.Data()),"recreate");
  
 
  // Statistical part
  Double_t yErrStat[NN];
  Double_t yStat[NN];
  Double_t xStat[NN];
  for (Int_t i=0;i<NN;i++){
    yStat[i] = grStat->GetY()[i]*dnchArray[i];
    xStat[i] = dnchArray[i];
    yErrStat[i] = TMath::Abs(dnchArray[i]*grStat->GetEY()[i]);
  }
  
  TGraphErrors *grStatScaled = new TGraphErrors(NN,xStat,yStat,0,yErrStat);
  grStatScaled->SetName(grStatName);
 
  
  //Systematic part
  Double_t yErrSyst[NN];
  Double_t ySyst[NN];
  Double_t xSyst[NN];
  for (Int_t i=0;i<NN;i++){
    Double_t errY     = grSyst->GetEYhigh()[NN-1-i];
    Double_t errYdnch = dnchArrayErr[i];
    ySyst[i] = grSyst->GetY()[NN-1-i]*dnchArray[i];
    xSyst[i] = dnchArray[i];
    cout << errY << "  " << errYdnch << "            "  << grSyst->GetY()[NN-1-i] << "  " <<  dnchArray[i] << endl; 
    //      yErrSyst[i] =  grSyst->GetY()[i] ;
    yErrSyst[i] = TMath::Sqrt( (grSyst->GetY()[NN-1-i]*grSyst->GetY()[NN-1-i])*(errYdnch*errYdnch)+(dnchArray[i]*dnchArray[i])*(errY*errY) );
    cout << yErrSyst[i] << endl;
  }
  
  Double_t xErrSyst[NN]; for (Int_t xx=0;xx<NN; xx++) xErrSyst[xx]=15.; 
  TGraphErrors *grSystScaled = new TGraphErrors(NN,xSyst,ySyst,xErrSyst,yErrSyst);
  grSystScaled->SetName(grSystName);
    
    
  grSystScaled->GetXaxis()->SetTitle("Centrality (%)");
  grSystScaled->SetMarkerSize(2); grSystScaled->GetXaxis()->SetTitle("dN_{ch}/d#eta"); 
  grSystScaled->SetFillColor(kBlack); grSystScaled->SetMarkerStyle(20); grSystScaled->SetMarkerColor(kBlack);
  if (grSystName.Contains("pika")) grSystScaled->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]#timesdN_{ch}/d#eta");
  if (grSystName.Contains("prpi")) grSystScaled->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]#timesdN_{ch}/d#eta");
  if (grSystName.Contains("prka")) grSystScaled->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]#timesdN_{ch}/d#eta");
  grStatScaled->GetXaxis()->SetTitle("Centrality (%)");
  grStatScaled->SetMarkerSize(2); grStatScaled->GetXaxis()->SetTitle("dN_{ch}/d#eta"); 
  grStatScaled->SetFillColor(kBlack); grStatScaled->SetMarkerStyle(20); grStatScaled->SetMarkerColor(kBlack);
  if (grStatName.Contains("pika")) grStatScaled->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]#timesdN_{ch}/d#eta");
  if (grStatName.Contains("prpi")) grStatScaled->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]#timesdN_{ch}/d#eta");
  if (grStatName.Contains("prka")) grStatScaled->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]#timesdN_{ch}/d#eta");
  
  scaling->GetFile()->cd();
  grStatScaled->Write();
  grSystScaled->Write();
  delete grStatScaled;
  delete grSystScaled;
  delete scaling;
  
  
} 
// =====================================================================================================
void MultiplydNch(TString datafile, Int_t flag)
{

  //
  // correct the centrality axis of the moments
  //
  /*
   cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/test
   aliroot -l 
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
   MultiplydNch("PlotsMomentsStatErr_iter6.root",0)
   MultiplydNch("PlotsMomentsStatErr_iter6.root",1)
   MultiplydNch("PlotsMomentsStatErr_iter6.root",2)
   MultiplydNch("PlotsMomentsStatErr_iter6.root",4)
   MultiplydNch("PlotsMomentsStatErr_iter6.root",5)


  */
  TFile *finalSyst = TFile::Open("/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics/FinalErrors.root"); 
  TGraphAsymmErrors * prpiSystALICE = (TGraphAsymmErrors*)finalSyst->Get("prpisyst");
  TGraphAsymmErrors * prkaSystALICE = (TGraphAsymmErrors*)finalSyst->Get("prkastat");
  TGraphAsymmErrors * pikaSystALICE = (TGraphAsymmErrors*)finalSyst->Get("pikasyst");

  TGraphErrors * prpiStatALICE = (TGraphErrors*)finalSyst->Get("prpistat");
  TGraphErrors * prkaStatALICE = (TGraphErrors*)finalSyst->Get("prkastat");
  TGraphErrors * pikaStatALICE = (TGraphErrors*)finalSyst->Get("pikastat");

  TTreeSRedirector *scaling = new TTreeSRedirector(Form("Nudyn_%d_Scaling.root",flag),"recreate");
  TFile *f       = TFile::Open(datafile);
  TObjArray *arr = (TObjArray*)f->Get("Moments");
  TGraphErrors *grdn = (TGraphErrors*)arr->FindObject("pi1");
  
  TGraphErrors *grpi = (TGraphErrors*)arr->FindObject("pi1");
  TGraphErrors *grka = (TGraphErrors*)arr->FindObject("ka1");
  TGraphErrors *grpr = (TGraphErrors*)arr->FindObject("pr1");

  TString nuName;
  TGraphErrors *tmpgr;
  for (Int_t igr=0; igr<3;igr++){
    
    if(igr==0) {nuName = "pikaNuDyn"; tmpgr=(TGraphErrors*)pikaSystALICE->Clone(); }
    if(igr==1) {nuName = "prpiNuDyn"; tmpgr=(TGraphErrors*)prpiSystALICE->Clone(); }
    if(igr==2) {nuName = "prkaNuDyn"; tmpgr=(TGraphErrors*)prkaSystALICE->Clone(); }
    
    TGraphErrors *grNu = (TGraphErrors*)arr->FindObject(nuName);
    Int_t NN = grNu->GetN(); 
    Double_t yErr[NN];
    Double_t y[NN];
    Double_t x[NN];
    Double_t errY;
    Double_t errYdnch;
    TString nuName = (TString)(grNu->GetName());
    
    for (Int_t i=0;i<NN;i++){
      
      errY = TMath::Sqrt( grNu->GetEY()[i]*grNu->GetEY()[i] + tmpgr->GetY()[i]*tmpgr->GetY()[i] );
      
      // scaling wrt to number of pions 
      if (flag==0) {
	y[i] = grNu->GetY()[i]*grdn->GetY()[i];
	x[i] = grNu->GetX()[i];
	errY     = grNu->GetEY()[i];
	errYdnch = grdn->GetEY()[i];
	cout << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << grdn->GetY()[i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(grdn->GetY()[i]*grdn->GetY()[i])*(errY*errY) );
      }
      
      // scaling wrt to nPart
      if (flag==1) {
	y[i] = grNu->GetY()[i]*nPartArr[i];
	x[i] = grNu->GetX()[i];
	errYdnch = nPartArrErr[i];
	cout << nuName << "  " << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << nPartArr[i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(nPartArr[i]*nPartArr[i])*(errY*errY) );
      }
      
      // scaling wrt to dnch
      if (flag==2) {
	y[i] = grNu->GetY()[i]*dnchArray[i];
	x[i] = grNu->GetX()[i];
	errYdnch = dnchArrayErr[i];
	cout << nuName << "  " << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << dnchArray[i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(dnchArray[i]*dnchArray[i])*(errY*errY) );
      }
      
      // Anar
      Double_t a=0.,b=0.,c=0.,aerr=0.,berr=0.,cerr=0.;
      if (flag==4) {
	a = grNu->GetY()[i];
// 	aerr = grNu->GetEY()[i];
	aerr = TMath::Sqrt( grNu->GetEY()[i]*grNu->GetEY()[i] + tmpgr->GetY()[i]*tmpgr->GetY()[i] );
	if (igr==0) {
	  y[i] = grNu->GetY()[i]/(1./grpi->GetY()[i]+1./grka->GetY()[i]); 
	  b=grpi->GetY()[i];c=grka->GetY()[i];
	  berr=grpi->GetEY()[i];cerr=grka->GetEY()[i];
	}
	if (igr==1) {
	  y[i] = grNu->GetY()[i]/(1./grpi->GetY()[i]+1./grpr->GetY()[i]); 
	  b=grpi->GetY()[i];c=grpr->GetY()[i];
	  berr=grpi->GetEY()[i];cerr=grpr->GetEY()[i];	  
	}
	if (igr==2) {
	  y[i] = grNu->GetY()[i]/(1./grpr->GetY()[i]+1./grka->GetY()[i]); 
	  b=grpr->GetY()[i];c=grka->GetY()[i];
	  berr=grpr->GetEY()[i];cerr=grka->GetEY()[i];	  
	}
	x[i]           = grNu->GetX()[i];
	cout << i << "  " << grNu->GetY()[i] << "  " << 1./grpi->GetY()[i]+1./grka->GetY()[i] << "  " << y[i] <<  "  " << grpi->GetY()[i] << "  " << grka->GetY()[i] << "  " << grpr->GetY()[i] << endl;
	yErr[i] = TMath::Sqrt(  (1./(1./b+1./c))*(1./(1./b+1./c))*aerr*aerr +
	                       ((b*b*a)/((b+c)*(b+c)))*((b*b*a)/((b+c)*(b+c)))*cerr*cerr +
	                       ((c*c*a)/((b+c)*(b+c)))*((c*c*a)/((b+c)*(b+c)))*berr*berr );
	cout << yErr[i] << endl;
      }  
      
      // scaling wrt to dnch
      if (flag==5) {
	y[i] = grNu->GetY()[i]*dnchArray[i];
	x[i] = dnchArray[i];
	errYdnch = dnchArrayErr[i];
	cout << nuName << "  " << i<< "   "  << errY << "  " << errYdnch << "  " << grNu->GetY()[i] << "  " << dnchArray[i] << endl; 
	yErr[i] = TMath::Sqrt( (grNu->GetY()[i]*grNu->GetY()[i])*(errYdnch*errYdnch)+(dnchArray[i]*dnchArray[i])*(errY*errY) );
      }
      
      
    }
    
    TGraphErrors *gr = new TGraphErrors(NN,x,y,0,yErr);
    gr->SetDrawOption("a3");
    gr->SetName(nuName);
    gr->GetXaxis()->SetTitle("Centrality (%)");
    gr->SetMarkerSize(2);
    if (flag==0 || flag==2 || flag==5 ) {
      if (nuName.Contains("pika")) gr->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]xdN_{ch}/d#eta");
      if (nuName.Contains("prpi")) gr->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]xdN_{ch}/d#eta");
      if (nuName.Contains("prka")) gr->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]xdN_{ch}/d#eta");
    }
    if (flag==1) {
      if (nuName.Contains("pika")) gr->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]xN_{part}");
      if (nuName.Contains("prpi")) gr->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]xN_{part}");
      if (nuName.Contains("prka")) gr->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]xN_{part}");
    }
    if (flag==4) {
      if (nuName=="pikaNuDyn") gr->GetYaxis()->SetTitle(" #frac{#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]}{1/<#pi>+1/<K>} ");
      if (nuName=="prpiNuDyn") gr->GetYaxis()->SetTitle(" #frac{#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]}{1/<#pi>+1/<p>} ");
      if (nuName=="prkaNuDyn") gr->GetYaxis()->SetTitle(" #frac{#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]}{1/<p>+1/<K>} ");
      gr->GetYaxis()->SetTitleOffset(1.5);
    }
    if (flag==5) gr->GetXaxis()->SetTitle("dN_{ch}/d#eta");
    
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
TGraphErrors * MakeScaledMoments(TGraphErrors *grCorr, TGraphErrors *grpi)
{
  
  //
  // correct the centrality axis of the moments
  //
  
  Int_t NN = grCorr->GetN();
  Double_t yCorr[NN], xCorr[NN];
  for (Int_t n=0;n<NN;n++){
    yCorr[n]  = grCorr->GetY()[n]*grpi->GetY()[n];
    xCorr[n]  = grCorr->GetX()[n];
  }
  TGraphErrors *gr = new TGraphErrors(NN,xCorr,yCorr,0,0);
  gr->SetDrawOption("a3");
  gr->SetName(grCorr->GetName());
  gr->GetXaxis()->SetTitle("Centrality (%)");
  gr->GetYaxis()->SetTitle(grCorr->GetName());
  gr->SetFillColor(kBlack);
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(kBlack);
  return gr;
  
} 
// =====================================================================================================
void PlotNpartScalings(TString data, TString hijing, TString ampt)
{
  
  //
  /*
   *   
   *   cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/test
   *   // npion
   *   aliroot -l 
   *   TString d = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/test/Nudyn_0_Scaling.root"
   *   TString a = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_AMPT_TightCuts/test/Nudyn_0_Scaling.root"
   *   TString h = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_HIJING_TightCuts/test/Nudyn_0_Scaling.root"
   *   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
   *   PlotNpartScalings(d,h,a)
   *   
   *   // npar
   *   aliroot -l 
   *   TString d = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/test/Nudyn_1_Scaling.root"
   *   TString a = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_AMPT_TightCuts/test/Nudyn_1_Scaling.root"
   *   TString h = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_HIJING_TightCuts/test/Nudyn_1_Scaling.root"
   *   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
   *   PlotNpartScalings(d,h,a)
   *  
   *   // nch
   *   aliroot -l 
   *   TString d = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/test/Nudyn_2_Scaling.root"
   *   TString a = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_AMPT_TightCuts/test/Nudyn_2_Scaling.root"
   *   TString h = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_HIJING_TightCuts/test/Nudyn_3_Scaling.root"
   *   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
   *   PlotNpartScalings(d,h,a)
   *   
   *   // anar
      aliroot -l 
      TString d = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/test/Nudyn_4_Scaling.root"
      TString a = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_AMPT_TightCuts/test/Nudyn_4_Scaling.root"
      TString h = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_HIJING_TightCuts/test/Nudyn_4_Scaling.root"
      .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
     PlotNpartScalings(d,h,a)
    
      // nudyn vs nch
      aliroot -l 
      TString d = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/test/Nudyn_5_Scaling.root"
      TString a = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_AMPT_TightCuts/test/Nudyn_5_Scaling.root"
      TString h = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_HIJING_TightCuts/test/Nudyn_6_Scaling.root"
      .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
      PlotNpartScalings(d,h,a)
     
     
   */
  
  TFile *fd = TFile::Open(data);
  TFile *fh = TFile::Open(hijing);
  TFile *fa = TFile::Open(ampt);
  
  TTreeSRedirector *scaling = new TTreeSRedirector("Scaling_Summary.root","recreate");
  
  
  TGraphErrors *data_pika = (TGraphErrors*)fd->Get("pikaNuDyn"); Double_t nMaxX = TMath::MaxElement(data_pika->GetN(),data_pika->GetX());
  data_pika->SetMarkerColor(kRed+1); data_pika->SetMarkerStyle(20); data_pika->SetLineColor(kRed+1); data_pika->SetMarkerSize(1.5); 
  TGraphErrors *data_pipr = (TGraphErrors*)fd->Get("prpiNuDyn"); 
  data_pipr->SetMarkerColor(kRed+1); data_pipr->SetMarkerStyle(20); data_pipr->SetLineColor(kRed+1); data_pipr->SetMarkerSize(1.5);
  TGraphErrors *data_kapr = (TGraphErrors*)fd->Get("prkaNuDyn"); 
  data_kapr->SetMarkerColor(kRed+1); data_kapr->SetMarkerStyle(20); data_kapr->SetLineColor(kRed+1); data_kapr->SetMarkerSize(1.5);
  
  TGraphErrors *hijing_pika = (TGraphErrors*)fh->Get("pikaNuDyn"); 
  hijing_pika->SetMarkerColor(kBlack); hijing_pika->SetMarkerStyle(1); hijing_pika->SetLineColor(kBlack); hijing_pika->SetLineWidth(3);
  TGraphErrors *hijing_pipr = (TGraphErrors*)fh->Get("piprNuDyn"); 
  hijing_pipr->SetMarkerColor(kBlack); hijing_pipr->SetMarkerStyle(1); hijing_pipr->SetLineColor(kBlack); hijing_pipr->SetLineWidth(3);
  TGraphErrors *hijing_kapr = (TGraphErrors*)fh->Get("kaprNuDyn"); 
  hijing_kapr->SetMarkerColor(kBlack); hijing_kapr->SetMarkerStyle(1); hijing_kapr->SetLineColor(kBlack); hijing_kapr->SetLineWidth(3);
  
  TGraphErrors *ampt_pika = (TGraphErrors*)fa->Get("pikaNuDyn"); 
  ampt_pika->SetMarkerColor(kGreen+2); ampt_pika->SetMarkerStyle(1); ampt_pika->SetLineColor(kGreen+2); ampt_pika->SetLineWidth(3);
  TGraphErrors *ampt_pipr = (TGraphErrors*)fa->Get("piprNuDyn"); 
  ampt_pipr->SetMarkerColor(kGreen+2); ampt_pipr->SetMarkerStyle(1); ampt_pipr->SetLineColor(kGreen+2); ampt_pipr->SetLineWidth(3);
  TGraphErrors *ampt_kapr = (TGraphErrors*)fa->Get("kaprNuDyn"); 
  ampt_kapr->SetMarkerColor(kGreen+2); ampt_kapr->SetMarkerStyle(1); ampt_kapr->SetLineColor(kGreen+2); ampt_kapr->SetLineWidth(3);
  
  TCanvas *cpika = new TCanvas("cpika","cpika",900,700); cpika->Clear(); cpika->Update(); cpika->SetTicks(1,1); cpika->cd(); 
  data_pika->Draw("ap"); data_pika->GetYaxis()->SetRangeUser(-2.,8.);
  hijing_pika->Draw("Xl");
  ampt_pika->Draw("Xl");
  TLine* linepika = new TLine(0., 0.,nMaxX, 0.); linepika -> SetLineColor(2); linepika -> SetLineStyle(2); linepika -> SetLineWidth(3);
  linepika->Draw("same");
  
  TLegend legpika(0.20, 0.85, 0.5, 0.7);
  legpika.SetTextFont(62); legpika.SetTextSize(0.045); legpika.SetFillColor(0); legpika.SetBorderSize(0);
  legpika.AddEntry(data_pika  ," ALICE Data ","LPE");
  legpika.AddEntry(hijing_pika," HIJING ","LP");
  legpika.AddEntry(ampt_pika  ," AMPT ","LP");
  legpika.Draw("same");
  
  
  TCanvas *cpipr = new TCanvas("cpipr","cpipr",900,700); cpipr->Clear(); cpipr->Update(); cpipr->SetTicks(1,1); cpipr->cd(); 
  data_pipr->Draw("ap"); data_pipr->GetYaxis()->SetRangeUser(-2.,8.);
  hijing_pipr->Draw("Xl");
  ampt_pipr->Draw("Xl");
  TLine* linepipr = new TLine(0., 0.,nMaxX, 0.); linepipr -> SetLineColor(2); linepipr -> SetLineStyle(2); linepipr -> SetLineWidth(3);
  linepipr->Draw("same");
  
  TLegend legpipr(0.20, 0.85, 0.5, 0.7);
  legpipr.SetTextFont(62); legpipr.SetTextSize(0.045); legpipr.SetFillColor(0); legpipr.SetBorderSize(0);
  legpipr.AddEntry(data_pipr  ," ALICE Data ","LPE");
  legpipr.AddEntry(hijing_pipr," HIJING ","LP");
  legpipr.AddEntry(ampt_pipr  ," AMPT ","LP");
  legpipr.Draw("same");
  
  
  TCanvas *ckapr = new TCanvas("ckapr","cpika",900,700); ckapr->Clear(); ckapr->Update(); ckapr->SetTicks(1,1); ckapr->cd(); 
  data_kapr->Draw("ap"); data_kapr->GetYaxis()->SetRangeUser(-2.,8.);
  hijing_kapr->Draw("Xl");
  ampt_kapr->Draw("Xl");
  TLine* linekapr = new TLine(0., 0.,nMaxX, 0.); linekapr -> SetLineColor(2); linekapr -> SetLineStyle(2); linekapr -> SetLineWidth(3);
  linekapr->Draw("same");
  
  TLegend legpkapr(0.20, 0.85, 0.5, 0.7);
  legpkapr.SetTextFont(62); legpkapr.SetTextSize(0.045); legpkapr.SetFillColor(0); legpkapr.SetBorderSize(0);
  legpkapr.AddEntry(data_kapr  ," ALICE Data ","LPE");
  legpkapr.AddEntry(hijing_kapr," HIJING ","LP");
  legpkapr.AddEntry(ampt_kapr  ," AMPT ","LP");
  legpkapr.Draw("same");
  
  
  
  scaling->GetFile()->cd();
  cpika->Write();
  cpipr->Write();
  ckapr->Write();
  
  delete scaling;
  
}
// =====================================================================================================
void PlotdnchScalings(TString data, TString hijing, TString ampt)
{
  
  //
  /*
   *   
  
    
      // nudyn vs nch
      aliroot -l 
      TString d = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_16EtaBin_mombin20MeV_TightCuts/test/Nudyn_5_Scaling.root"
      TString a = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_AMPT_TightCuts/test/Nudyn_6_Scaling.root"
      TString h = "/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_HIJING_TightCuts/test/Nudyn_5_Scaling.root"
      .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
      PlotdnchScalings(d,h,a)
     
     
   */
//   TGaxis::SetMaxDigits(3);
  gStyle->SetEndErrorSize(3);
 
  
  TFile *finalSyst = TFile::Open("/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Systematics/Nudyn_dnch_Scaling.root"); 
  TGraphAsymmErrors * prpiSystALICE = (TGraphAsymmErrors*)finalSyst->Get("prpisyst");
  TGraphAsymmErrors * prkaSystALICE = (TGraphAsymmErrors*)finalSyst->Get("prkasyst");
  TGraphAsymmErrors * pikaSystALICE = (TGraphAsymmErrors*)finalSyst->Get("pikasyst");
  TGraphErrors * prpiStatALICE = (TGraphErrors*)finalSyst->Get("prpistat");
  TGraphErrors * prkaStatALICE = (TGraphErrors*)finalSyst->Get("prkastat");
  TGraphErrors * pikaStatALICE = (TGraphErrors*)finalSyst->Get("pikastat");
  prpiSystALICE->SetFillStyle(0);
  prkaSystALICE->SetFillStyle(0);
  pikaSystALICE->SetFillStyle(0);
  prpiStatALICE->SetFillStyle(0);
  prkaStatALICE->SetFillStyle(0);
  pikaStatALICE->SetFillStyle(0);
   
  TLegend legLogo(0.35, 0.45, 0.70, 0.55);  legLogo.SetBorderSize(0);
  legLogo.SetTextFont(62);
  legLogo.SetTextSize(0.03);
  legLogo.SetFillColor(0); 
  legLogo.SetFillStyle(4000);
  legLogo.AddEntry((TObject*)0, "ALICE Preliminary, Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV"," ");
  legLogo.AddEntry((TObject*)0, "   |#eta|<0.8, 0.2<#it{p}<1.5 GeV/#it{c} ", " "); 
  

  pikaStatALICE->SetMarkerColor(kRed+1); pikaStatALICE->SetMarkerStyle(20); pikaStatALICE->SetLineColor(kRed+1); pikaStatALICE->SetMarkerSize(1.5); 
  prkaStatALICE->SetMarkerColor(kRed+1); prkaStatALICE->SetMarkerStyle(20); prkaStatALICE->SetLineColor(kRed+1); prkaStatALICE->SetMarkerSize(1.5); 
  prpiStatALICE->SetMarkerColor(kRed+1); prpiStatALICE->SetMarkerStyle(20); prpiStatALICE->SetLineColor(kRed+1); prpiStatALICE->SetMarkerSize(1.5); 
  
  pikaSystALICE->SetMarkerColor(kRed+1); pikaSystALICE->SetMarkerStyle(1); pikaSystALICE->SetLineColor(kRed+1); pikaSystALICE->SetMarkerSize(1.5); 
  prkaSystALICE->SetMarkerColor(kRed+1); prkaSystALICE->SetMarkerStyle(1); prkaSystALICE->SetLineColor(kRed+1); prkaSystALICE->SetMarkerSize(1.5); 
  prpiSystALICE->SetMarkerColor(kRed+1); prpiSystALICE->SetMarkerStyle(1); prpiSystALICE->SetLineColor(kRed+1); prpiSystALICE->SetMarkerSize(1.5); 

  TFile *fd = TFile::Open(data);
  TFile *fh = TFile::Open(hijing);
  TFile *fa = TFile::Open(ampt);
  
  TTreeSRedirector *scaling = new TTreeSRedirector("Scaling_Summary.root","recreate");
  
  
  TGraphErrors *data_pika = (TGraphErrors*)fd->Get("pikaNuDyn"); 
  data_pika->SetMarkerColor(kRed+1); data_pika->SetMarkerStyle(20); data_pika->SetLineColor(kRed+1); data_pika->SetMarkerSize(1.5); 
  TGraphErrors *data_pipr = (TGraphErrors*)fd->Get("prpiNuDyn"); 
  data_pipr->SetMarkerColor(kRed+1); data_pipr->SetMarkerStyle(20); data_pipr->SetLineColor(kRed+1); data_pipr->SetMarkerSize(1.5);
  TGraphErrors *data_kapr = (TGraphErrors*)fd->Get("prkaNuDyn"); 
  data_kapr->SetMarkerColor(kRed+1); data_kapr->SetMarkerStyle(20); data_kapr->SetLineColor(kRed+1); data_kapr->SetMarkerSize(1.5);
  
  TGraphErrors *hijing_pika = (TGraphErrors*)fh->Get("pikaNuDyn"); 
  hijing_pika->SetMarkerColor(kBlack); hijing_pika->SetMarkerStyle(1); hijing_pika->SetLineColor(kBlack); hijing_pika->SetLineWidth(3);
  TGraphErrors *hijing_pipr = (TGraphErrors*)fh->Get("piprNuDyn"); 
  hijing_pipr->SetMarkerColor(kBlack); hijing_pipr->SetMarkerStyle(1); hijing_pipr->SetLineColor(kBlack); hijing_pipr->SetLineWidth(3);
  TGraphErrors *hijing_kapr = (TGraphErrors*)fh->Get("kaprNuDyn"); 
  hijing_kapr->SetMarkerColor(kBlack); hijing_kapr->SetMarkerStyle(1); hijing_kapr->SetLineColor(kBlack); hijing_kapr->SetLineWidth(3);
  
  TGraphErrors *ampt_pika = (TGraphErrors*)fa->Get("pikaNuDyn"); 
  ampt_pika->SetMarkerColor(kGreen+2); ampt_pika->SetMarkerStyle(1); ampt_pika->SetLineColor(kGreen+2); ampt_pika->SetLineWidth(3);
  TGraphErrors *ampt_pipr = (TGraphErrors*)fa->Get("piprNuDyn"); 
  ampt_pipr->SetMarkerColor(kGreen+2); ampt_pipr->SetMarkerStyle(1); ampt_pipr->SetLineColor(kGreen+2); ampt_pipr->SetLineWidth(3);
  TGraphErrors *ampt_kapr = (TGraphErrors*)fa->Get("kaprNuDyn"); 
  ampt_kapr->SetMarkerColor(kGreen+2); ampt_kapr->SetMarkerStyle(1); ampt_kapr->SetLineColor(kGreen+2); ampt_kapr->SetLineWidth(3);
  
  Double_t nMaxX = TMath::MaxElement(hijing_pika->GetN(),hijing_pika->GetX());
  
  // Legend for data points
  TLegend legdata(0.23, 0.65, 0.55, 0.86); //legdata.SetNColumns(2);
  legdata.SetTextFont(62); legdata.SetTextSize(0.03); legdata.SetFillColor(0); legdata.SetBorderSize(0);
  legdata.SetFillStyle(4000);
  legdata.AddEntry(pikaStatALICE  , "ALICE Data, stat. errors ","LPE");
  legdata.AddEntry(pikaStatALICE,   "Systematic uncertainty",  "F");
  legdata.AddEntry(hijing_pika, "HIJING ","LP");
  legdata.AddEntry(ampt_pika  , "AMPT ","LP");
  TLine* line0 = new TLine(0., 0.,nMaxX, 0.); line0 -> SetLineStyle(2); line0 -> SetLineWidth(3);

  
  TCanvas *cpika = new TCanvas("cpika","cpika",900,700); cpika->Clear(); cpika->Update(); cpika->SetTicks(1,1); cpika->cd(); 
  if (transparent){
    cpika->SetFillStyle(4000);  
    cpika->SetFrameFillStyle(4000); 
  }
  pikaSystALICE->GetYaxis()->SetRangeUser(-0.5,3.5);
  pikaSystALICE->GetYaxis()->SetTitleOffset(1.);
  pikaSystALICE->Draw("a5"); pikaStatALICE->Draw("p"); 
  hijing_pika->Draw("Xl");
  ampt_pika->Draw("Xl");
  line0->Draw("same");
  legdata.Draw("same");
  legLogo.Draw("same");
  
 
  TCanvas *cpipr = new TCanvas("cpipr","cpipr",900,700); cpipr->Clear(); cpipr->Update(); cpipr->SetTicks(1,1); cpipr->cd();
  if (transparent){ 
    cpipr->SetFillStyle(4000);   
    cpipr->SetFrameFillStyle(4000); 
  }
  prpiSystALICE->GetYaxis()->SetRangeUser(-1.,5.);
  prpiSystALICE->GetYaxis()->SetTitleOffset(1.);
  prpiSystALICE->Draw("a5"); prpiStatALICE->Draw("p"); 
  hijing_pipr->Draw("Xl");
  ampt_pipr->Draw("Xl");
  line0->Draw("same");
  legdata.Draw("same");
  legLogo.Draw("same");

  
 
  TCanvas *ckapr = new TCanvas("ckapr","cpika",900,700); ckapr->Clear(); ckapr->Update(); ckapr->SetTicks(1,1); ckapr->cd(); 
  if (transparent){ 
    ckapr->SetFillStyle(4000);  
    ckapr->SetFrameFillStyle(4000); 
  }
  prkaSystALICE->GetYaxis()->SetRangeUser(-1.,8.);
  prkaSystALICE->GetYaxis()->SetTitleOffset(1.);
  prkaSystALICE->Draw("a5");  prkaStatALICE->Draw("p");
  hijing_kapr->Draw("Xl");
  ampt_kapr->Draw("Xl");
  line0->Draw("same");
  legdata.Draw("same");
  legLogo.Draw("same");

  scaling->GetFile()->cd();
  cpika->Write();
  cpipr->Write();
  ckapr->Write();
  
  delete scaling;
  
}
// =====================================================================================================
TGraphErrors * RatioGraph(TGraphErrors *g1, TGraphErrors *g2) 
{
  
  const Int_t N = g1->GetN();
  Double_t x[N];
  Double_t y[N]; 
  for (Int_t i=0; i<N; i++){    
    if ((g2->GetY()[i]!=0.)) {
      y[i]=TMath::Abs(g1->GetY()[i]/g2->GetY()[i]);
    } else {
      y[i]=0.;
    }
    x[i]=g1->GetX()[i];
  }
        
  TGraphErrors *grRatio = new TGraphErrors(N,x,y,0,0);  
  return grRatio;
}
// =====================================================================================================
TGraphErrors * DiffGraph(TGraphErrors *g1, TGraphErrors *g2) 
{
  
  const Int_t N = g1->GetN();
  Double_t x[N];
  Double_t y[N]; 
  for (Int_t i=0; i<N; i++){
    if ((g2->GetY()[i]!=0.)) {
      y[i]=TMath::Abs(g1->GetY()[i]-g2->GetY()[i]);
    } else {
      y[i]=0.;
    }
    x[i]=g1->GetX()[i];
  }
        
  TGraphErrors *grRatio = new TGraphErrors(N,x,y,0,0);  
  return grRatio;
}
// =====================================================================================================
void AnalyseWsApprovedPlots(TString wFile, TString wSliceFile)
{

  //
  // Calculate statistical errors using subsample method
  /*
  ##cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/IdDataTrees_300000_Results_OK2
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/IdentityFiles/Wanalysis
  aliroot -l
   .L /u/marsland/PHD/macros/marsland_EbyeRatios/TIdentityStatisticalError.C+
   TString inputFile1="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/PLOTS_AnaNote/rootfiles/streamer_all_300000.root"
   TString inputFile2="/lustre/nyx/alice/users/marsland/pFluct/files/analysis/PLOTS_AnaNote/rootfiles/streamer_slice_300000.root"
   AnalyseWsApprovedPlots(inputFile1,inputFile2)
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


// +++++++++++++++++++++++++ Plot Ratios +++++++++++++++++++++++++

/*
 *  TGraphErrors *grRatioAmpt[3];
 *  TGraphErrors *grRatioData[3];
 *  const Int_t nPoints=grData[0]->GetN();
 *  cout << "nPoints = " << nPoints << endl;
 *  Double_t yAmpt0[nPoints];
 *  Double_t yAmpt1[nPoints];
 *  Double_t yAmpt2[nPoints];
 *  Double_t yData0[nPoints];
 *  Double_t yData1[nPoints];
 *  Double_t yData2[nPoints];
 *  
 *  grRatioAmpt[0] = new TGraphErrors(nPoints,grData[0]->GetX(),yAmpt0,0,0);
 *  grRatioAmpt[1] = new TGraphErrors(nPoints,grData[0]->GetX(),yAmpt1,0,0);
 *  grRatioAmpt[2] = new TGraphErrors(nPoints,grData[0]->GetX(),yAmpt2,0,0);
 *  grRatioData[0] = new TGraphErrors(nPoints,grData[0]->GetX(),yData0,0,0);
 *  grRatioData[1] = new TGraphErrors(nPoints,grData[0]->GetX(),yData1,0,0);
 *  grRatioData[2] = new TGraphErrors(nPoints,grData[0]->GetX(),yData2,0,0);
 * 
 *  for (Int_t i=0;i<nPoints;i++){
 *    cout << " i = " << i << "      " << grHijing[0]->GetX()[i] << endl;
 *    cout << grHijing[0]->GetY()[i] <<"  " << grAmpt[0]->GetY()[i] << "  " << grData[0]->GetY()[nPoints-1-i] << endl;
 *    cout << grHijing[1]->GetY()[i] <<"  " << grAmpt[1]->GetY()[i] << "  " << grData[1]->GetY()[nPoints-1-i] << endl;
 *    cout << grHijing[2]->GetY()[i] <<"  " << grAmpt[2]->GetY()[i] << "  " << grData[2]->GetY()[i] << endl;
 *    cout << " ------------------ " << endl;
 *    grRatioAmpt[0] ->SetPoint(i, grHijing[0]->GetX()[i], grHijing[0]->GetY()[i]-grAmpt[0]->GetY()[i]);
 *    grRatioAmpt[1] ->SetPoint(i, grHijing[1]->GetX()[i], grHijing[1]->GetY()[i]-grAmpt[1]->GetY()[i]);
 *    grRatioAmpt[2] ->SetPoint(i, grHijing[2]->GetX()[i], grHijing[2]->GetY()[i]-grAmpt[2]->GetY()[i]);
 *    
 *    grRatioData[0] ->SetPoint(i, grHijing[0]->GetX()[i], grHijing[0]->GetY()[i]-grData[0]->GetY()[nPoints-1-i]);
 *    grRatioData[1] ->SetPoint(i, grHijing[1]->GetX()[i], grHijing[1]->GetY()[i]-grData[1]->GetY()[nPoints-1-i]);
 *    grRatioData[2] ->SetPoint(i, grHijing[2]->GetX()[i], grHijing[2]->GetY()[i]-grData[2]->GetY()[nPoints-1-i]);
 *  }
 *  
 *  
 * 
 *  
 *  for (Int_t i=0; i<3; i++){
 *    grRatioData[i]->SetMarkerStyle(20);   grRatioData[i]->SetMarkerColor(kBlack);     grRatioData[i]->SetLineColor(kBlack);  
 *    grRatioAmpt[i]->SetMarkerStyle(21);   grRatioAmpt[i]->SetMarkerColor(kGreen+2);   grRatioAmpt[i]->SetLineColor(kGreen+2);  
 *     
 *    grRatioData[i]  ->SetMarkerSize(1.5);
 *    grRatioAmpt[i]  ->SetMarkerSize(1.5);
 * 
 *    grRatioData[i]->GetXaxis()->SetTitle("Centrality (%)"); 
 *    grRatioAmpt[i]->GetXaxis()->SetTitle("Centrality (%)"); 
 * 
 *    grRatioData[i]->GetXaxis()->SetTitleSize(0.05); grRatioData[i]->GetXaxis()->SetLabelSize(0.045); 
 *    grRatioData[i]->GetYaxis()->SetTitleSize(0.05); grRatioData[i]->GetYaxis()->SetLabelSize(0.045); 
 *    
 *    grRatioAmpt[i]->GetXaxis()->SetTitleSize(0.05); grRatioAmpt[i]->GetXaxis()->SetLabelSize(0.045); 
 *    grRatioAmpt[i]->GetYaxis()->SetTitleSize(0.05); grRatioAmpt[i]->GetYaxis()->SetLabelSize(0.045); 
 * 
 *  }
 *  grRatioData[0]  ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
 *  grRatioData[1]  ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
 *  grRatioData[2]  ->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
 *   
 *  grRatioAmpt[0]  ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},K^{+}+K^{-}]"); 
 *  grRatioAmpt[1]  ->GetYaxis()->SetTitle("#nu_{dyn} [#pi^{+}+#pi^{-},p+#bar{p}]"); 
 *  grRatioAmpt[2]  ->GetYaxis()->SetTitle("#nu_{dyn} [p+#bar{p},K^{+}+K^{-}]"); 
 *  
 *  
 *  
 *  // prepare alegend for all plots
 *  TLegend *legRatio = new TLegend(0.20, 0.85, 0.5, 0.7);
 *  legRatio->SetTextFont(62); legRatio->SetTextSize(0.045); legRatio->SetFillColor(0);
 *  legRatio->AddEntry(grRatioData[0]  , " Data ","LPE");
 *  legRatio->AddEntry(grRatioAmpt[0]  , " Ampt ","LPE");
 *  
 *   // pika comparison
 *  TCanvas *pikaRatio = new TCanvas("pikaRatio", "pika comparison wrt HIJING", 800, 600);   
 *  pikaRatio->cd();  pikaRatio->SetGridy();
 *  grRatioData[0]->Draw("alp");
 *  grRatioAmpt[0]->Draw("lp");
 *  legRatio->Draw("same");
 *  
 *  // pika comparison
 *  TCanvas *piprRatio = new TCanvas("piprRatio", "pipr comparison wrt HIJING", 800, 600);   
 *  piprRatio->cd(); piprRatio->SetGridy();
 *  grRatioData[1]->Draw("alp");
 *  grRatioAmpt[1]->Draw("lp");
 *  legRatio->Draw("same");
 *  
 *  // pika comparison
 *  TCanvas *kaprRatio = new TCanvas("kaprRatio", "kapr comparison wrt HIJING", 800, 600);   
 *  kaprRatio->cd(); kaprRatio->SetGridy();
 *  grRatioData[2]->Draw("alp");
 *  grRatioAmpt[2]->Draw("lp");
 *  legRatio->Draw("same");
 *  
 *  pikaRatio->SaveAs("pikaRatio.pdf");
 *  piprRatio->SaveAs("piprRatio.pdf");
 *  kaprRatio->SaveAs("kaprRatio.pdf");
 */


