// #include <TFile.h>
// #include <TH3.h>
// #include "THn.h"
// #include "TCut.h"
// #include "TCutG.h"
// #include <TCanvas.h>
// #include <TStyle.h>
// #include <iostream>
// #include <THnSparse.h>
// #include "TMath.h"
// #include "TLine.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "TROOT.h"
// #include "TSystem.h"
// #include "TDirectory.h"
// #include "TTreeStream.h"
// #include "TFile.h"
// #include "TChain.h"
// #include "TTree.h"
// #include "TCanvas.h"
// #include "TGraphErrors.h"
// #include "TLinearFitter.h"
// #include "TF1.h"
// #include "TStopwatch.h"
// #include "TLegend.h"
// #include "AliXRDPROOFtoolkit.h"
// #include "TStatToolkit.h"
// #include "AliMathBase.h"
// #include <fstream>
// #include <iostream>

using namespace std;

// 
// Make TH2D from the TTree chain which will be used for the pt slices
// 

void InitInitials();
void RealDataHists(const TString dataTreeFile, const TString histFile, const Double_t etaLow, const Double_t etaUp, const Float_t centLow, const Float_t centUp);
void RealDataHistsFromThnSparse(const TString cleanFile, const TString histFile, const Double_t etaLow, const Double_t etaUp, const Float_t centLow, const Float_t centUp);
void IdenDataHists(const TString idenTreeFile, const TString dataTreeFile, const TString histFile, const Double_t etaLow, const Double_t etaUp, const Float_t centLow, const Float_t centUp);
void MCDataHists(const TString mcFile, const Double_t etaLow, const Double_t etaUp, const Float_t centLow, const Float_t centUp);
void MakeMainSubSample(TString idenTreeFile,TString treeName,Int_t centBin);
void MakeSubSamples(TString idenTreeFile,TString treeName,Int_t centBin,Int_t subsample);
void MergeFitResultsTree(TString infile);
void CheckResults(TString paramFile);
void GetDeDxTree(TString list);
void GetIdenTree(TString idenlist, TString datalist);


    
// ======= Modification part ======= 
Bool_t MCclosure      = kFALSE;
Bool_t pp             = kFALSE;
Bool_t test           = kFALSE;
Int_t nSubsample      = 20;  
Int_t rangeCleanPions = 0; // 100: full centrality is taken, 0: centrality bins are taken into account 

Double_t dEdxMin    = 20;
Double_t dEdxMax    = 1020;
Double_t dEdxNbins  = 1000;
Double_t ptMin      = 0.2;
Double_t ptMax      = 3.2;
Double_t ptNbins    = 150;   // mostly for real data analysis
Double_t etaMin     = -1;
Double_t etaMax     = 1;
Int_t nEtabins      = 20;

// Double_t dEdxMin    = 20;
// Double_t dEdxMax    = 1020;
// Double_t dEdxNbins  = 400;
// Double_t ptMin      = 0.2;
// Double_t ptMax      = 3.2;
// Double_t ptNbins    = 150;   // mostly for real data analysis
// Double_t etaMin     = -1;
// Double_t etaMax     = 1;
// Int_t nEtabins      = 10;

const Int_t nCentbins = 9;         Float_t centArray[nCentbins] = {0,5,10,20,30,40,50,60,70};
const Int_t nCentBinsPlusOne = 10; Float_t xCentBins[nCentBinsPlusOne] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
// const Int_t nCentbins = 18; Float_t centArray[nCentbins] = {0, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75}; 
// const Int_t nCentBinsPlusOne = 19; Float_t xCentBins[nCentBinsPlusOne] = {0, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80};

// ======= Modification part =======  
Int_t nIter           = 8;        // not used anymore
Double_t testEntries  = 1000000.;
TTree *tree=NULL, *armtree=NULL, *treeIden=NULL; 
TTreeSRedirector *outStream=0;
TH1D  *hEta=0x0, *hCent=0x0, *hMom=0x0;
TCutG *pionCutG=0, *antiProtonCutG=0, *protonCutG=0;
TH2D *h2Dall, *h2Dpos, *h2Dneg;
 


void CreateAllTIDENHists(Int_t tidenSwitch, Int_t centBin, Int_t subsample, const TString idenTreeName,const TString idenTreeFile, const TString dataTreeFile, const TString mcFile, const TString histFile, const Double_t etaLow, const Double_t etaUp, const Float_t centLow, const Float_t centUp)
{
  //
  // Produce all hists form MC and Real data
  
  /*

  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_Reference/Hists
  
  aliroot -l
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C
  TString d = "../mergedPeriods/marsland_DataTree.list"
  TString m = "../mergedPeriods/marsland_MCmoments.root"
  TString h = "../mergedPeriods/marsland_Hists.root"
  TString i = "../mergedPeriods/marsland_IdenTree.list"
  
  //CreateAllTIDENHists(0,0,0,"","",d,m,h,-0.8,-0.7,10,20)      // for hists
  CreateAllTIDENHists(3,0,0,"",i,d,m,h,-0.8,-0.7,10,20)      // for hists form identree
  CreateAllTIDENHists(1,0,0,"fIdenTree",i,"","","",0,0,0,0)    // for main subsample
  CreateAllTIDENHists(2,0,10,"fIdenTree",i,"","","",0,0,0,0)    // for subsamples
  
  */
  //
  
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C .");
  InitInitials();
  if (tidenSwitch==0)
  {
    TString outputFileName = Form("Hists_PbPb_eta_%3.2f_%3.2f_cent_%3.2f_%3.2f.root",etaLow,etaUp,centLow,centUp);
    // check if file exist
    TFile *ftemp = TFile::Open(outputFileName);
    if (!ftemp) { 
      std::cout << " yolun acik ola --> file does not exist  " << std::endl;
    } else if (ftemp->GetSize()>1000) { 
      std::cout << " yavassss --> file harbiden produced  " << std::endl;  
      return;
    }
    outStream = new TTreeSRedirector(outputFileName,"recreate");
    TStopwatch timer; timer.Start();
    if (!MCclosure) RealDataHists(dataTreeFile,histFile,etaLow,etaUp,centLow,centUp);
    MCDataHists(mcFile,etaLow,etaUp,centLow,centUp);
    timer.Stop(); timer.Print();
    delete outStream;
  } else if (tidenSwitch==1)
  {
    TStopwatch timer; timer.Start();
    MakeMainSubSample(idenTreeFile,idenTreeName,centBin);  
    timer.Stop(); timer.Print();
  } else if (tidenSwitch==2)
  {
    TStopwatch timer; timer.Start();
    MakeSubSamples(idenTreeFile,idenTreeName,centBin,subsample);  
    timer.Stop(); timer.Print();
  } else if (tidenSwitch==3)
  {
    TString outputFileName = Form("Hists_PbPb_eta_%3.2f_%3.2f_cent_%3.2f_%3.2f.root",etaLow,etaUp,centLow,centUp);
    // check if file exist
    TFile *ftemp = TFile::Open(outputFileName);
    if (!ftemp) { 
      std::cout << " yolun acik ola --> file does not exist  " << std::endl;
    } else if (ftemp->GetSize()>1000) { 
      std::cout << " yavassss --> file harbiden produced  " << std::endl;  
      return;
    }
    outStream = new TTreeSRedirector(outputFileName,"recreate");
    TStopwatch timer; timer.Start();
    IdenDataHists(idenTreeFile,dataTreeFile,histFile,etaLow,etaUp,centLow,centUp);
    timer.Stop(); timer.Print();
  }  else if (tidenSwitch==4){
    TString outputFileName = Form("Hists_PbPb_eta_%3.2f_%3.2f_cent_%3.2f_%3.2f.root",etaLow,etaUp,centLow,centUp);
    TFile *ftemp = TFile::Open(outputFileName);
    if (!ftemp) { // check if file exist
      std::cout << " yolun acik ola --> file does not exist  " << std::endl;
    } else if (ftemp->GetSize()>1000) { 
      std::cout << " yavassss --> file harbiden produced  " << std::endl;  
      return;
    }
    outStream = new TTreeSRedirector(outputFileName,"recreate");
    TStopwatch timer; timer.Start();
    RealDataHistsFromThnSparse(dataTreeFile,histFile,etaLow,etaUp,centLow,centUp);
    MCDataHists(mcFile,etaLow,etaUp,centLow,centUp);
    timer.Stop(); timer.Print();
    delete outStream;
  }
   
}
//____________________________________________________________________________________________________________
void RealDataHists(const TString dataTreeFile, const TString histFile, const Double_t etaLow, const Double_t etaUp, const Float_t centLow, const Float_t centUp)
{
  
  //
  // Produce expected, real data, clean sample histograms 
  //
   
  if(!hEta) InitInitials();
//   Float_t etaBin  = hEta ->FindBin(etaLow)-1;
//   Float_t centBin = (centLow>=0) ? hCent->FindBin(centLow)-1 : hCent->FindBin(centLow+2)-1 ;  // cent<0 --> pp data
  std::cout << " eta Range  = " << etaLow  << " - " << etaUp; 
  std::cout << " cent Range = " << centLow << " - " << centUp << std::endl;

  // Read TTrees
  GetDeDxTree(dataTreeFile);
  
  Double_t nTreeEntriesAll=tree->GetEntries();
  if (nTreeEntriesAll<10) { 
    std::cout << " === upss data tree is empty === " << std::endl; return;
  }
  Double_t nTreeEntries = (test) ? testEntries : nTreeEntriesAll;

  // ======= Get the copy tree ======
  std::cout << " tree Size is --> " << nTreeEntriesAll << std::endl;
  TFile * tmpFile = new TFile(Form("tmp_eta_%3.2f_%3.2f_cent_%3.2f_%3.2f.root",etaLow,etaUp,centLow,centUp),"recreate");
  TString etaCentString = Form("eta>=%f && eta<%f && cent>=%f && cent<%f",etaLow,etaUp,centLow,centUp);
  TString etaCutArmtree = Form("eta>=%f   && eta<%f" ,etaLow ,etaUp);
  Long64_t nAll = (test) ? (Long64_t)testEntries : 100000000000;
  if (tree   ->GetEntries()<10) { 
    std::cout << " === upss main tree is empty === " << std::endl; return;  
  }
  TTree * treeRestricted = tree->CopyTree(etaCentString.Data(),"",nAll,0);
  treeRestricted->Write();    tree->Delete();    
  
  // Tree for clean samples it is not filled for MC closure test
  TTree * treeRestrictedArm = armtree->CopyTree(etaCutArmtree.Data(),"",nAll,0);
  if (armtree->GetEntries()>10) { 
    std::cout << " === arpPod tree is OK === " << std::endl; 
    treeRestrictedArm->Write(); armtree->Delete();
  }
   
  // ================================

  // 2D dEdx plots
  std::cout << " ========= make 2D dEdx plots ========= " << std::endl; 
  TCut parPos     = "sign>0.";  
  TCut parNeg     = "sign<0.";
  h2Dall  = new TH2D("h2Dall","h2Dall"  ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  h2Dpos  = new TH2D("h2Dpos","h2Dpos"  ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  h2Dneg  = new TH2D("h2Dneg","h2Dneg"  ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TStopwatch timer; timer.Start();
  
  std::cout << " make 2D dEdx plot ----- total Entry = " << treeRestricted->GetEntries() << "   running over = " << nTreeEntries << std::endl;
  treeRestricted ->Project("h2Dall" ,"dEdx:ptot" , ""     ,"",   nTreeEntries);
  treeRestricted ->Project("h2Dpos" ,"dEdx:ptot" , parPos ,"",2.*nTreeEntries);
  treeRestricted ->Project("h2Dneg" ,"dEdx:ptot" , parNeg ,"",2.*nTreeEntries);
  timer.Stop(); timer.Print();
 
  std::cout << " ========= make PID response plots ========= " << std::endl; 
  timer.Start();
  // Read THnSparses 
  TFile fhist(histFile);
  TList * list  = (TList*)fhist.Get("cleanHists");
  THnSparse *fhnExpected = (THnSparse*)list->FindObject("fhnExpected");
  //   THnSparse *fhnCleanEl  = (THnSparse*)list->FindObject("fhnCleanEl");
  THnSparse *fhnCleanKa  = (THnSparse*)list->FindObject("fhnCleanKa");
  TH2D *hArmPod          = (TH2D*)list->FindObject("fHistArmPod");

  // Make projections form THnSparse
  fhnExpected->GetAxis(2)->SetRangeUser(centLow,centUp);   // centrality
  fhnExpected->GetAxis(3)->SetRangeUser(etaLow,etaUp);     // eta
  fhnExpected->GetAxis(5)->SetRangeUser(2.,60.);          // sigma
  fhnExpected->GetAxis(6)->SetRangeUser(25.,1000.);        // mean
  //   fhnCleanEl->GetAxis(1)->SetRangeUser(centLow,centUp);    // centrality
  //   fhnCleanEl->GetAxis(2)->SetRangeUser(etaLow,etaUp);      // eta
  fhnCleanKa->GetAxis(1)->SetRangeUser(centLow,centUp);    // centrality
  fhnCleanKa->GetAxis(2)->SetRangeUser(etaLow,etaUp);      // eta
  
  // get each particle PID response
  fhnExpected->GetAxis(0)->SetRangeUser(0,1);              // electron
  TH2D * h2ExpectedEl      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedEl->SetName("h2ExpectedEl");
  TH2D * h2ExpectedSigmaEl = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaEl->SetName("h2ExpectedSigmaEl");
  fhnExpected->GetAxis(0)->SetRangeUser(1,2);  // pion
  TH2D * h2ExpectedPi      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedPi->SetName("h2ExpectedPi");
  TH2D * h2ExpectedSigmaPi = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaPi->SetName("h2ExpectedSigmaPi");
  fhnExpected->GetAxis(0)->SetRangeUser(2,3);  // kaon
  TH2D * h2ExpectedKa      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedKa->SetName("h2ExpectedKa");
  TH2D * h2ExpectedSigmaKa = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaKa->SetName("h2ExpectedSigmaKa");
  fhnExpected->GetAxis(0)->SetRangeUser(3,4);  // proton
  TH2D * h2ExpectedPr      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedPr->SetName("h2ExpectedPr");
  TH2D * h2ExpectedSigmaPr = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaPr->SetName("h2ExpectedSigmaPr");
  
  std::cout << " make clean Kaon Electron " << std::endl; 
  // Clean Kaon and Electron
  //   TH2D * h2CleanEl = (TH2D*)fhnCleanEl->Projection(4,3); h2CleanEl->SetName("h2CleanEl");
  TH2D * h2CleanKa = (TH2D*)fhnCleanKa->Projection(4,3); h2CleanKa->SetName("h2CleanKa");
  timer.Stop(); timer.Print();
  
  // Clean Pion and Proton
  TH2D *h2CleanPi=NULL, *h2CleanPr=NULL;
  if (!MCclosure){
    Double_t narmTreeEntries = (test) ? testEntries : treeRestrictedArm->GetEntries();
    std::cout << " make clean Pion Proton ---- total Entry = " << treeRestrictedArm->GetEntries() << std::endl;
    timer.Start();
    h2CleanPi = new TH2D("h2CleanPi","h2CleanPi",ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
    h2CleanPr = new TH2D("h2CleanPr","h2CleanPr",ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
    
    TCut cleanCutPi = "abs(piTOFnSigma)<3. && abs(alfa)<0.5  && pionCutG";
    TCut cleanCutPr = "abs(prTOFnSigma)<3. && qt>0.07 && qt<0.11 && (antiProtonCutG||protonCutG)";
    treeRestrictedArm -> Project("h2CleanPi" ,"dEdx:ptot",cleanCutPi,"",narmTreeEntries);
    treeRestrictedArm -> Project("h2CleanPr" ,"dEdx:ptot",cleanCutPr,"",narmTreeEntries);
    timer.Stop(); timer.Print();
  }
  
  // Write hists to the streamer
  outStream->GetFile()->cd();
  h2Dall            -> Write("h2Dall"); 
  h2Dpos            -> Write("h2Dpos");
  h2Dneg            -> Write("h2Dneg");
  h2ExpectedEl      -> Write("h2ExpectedEl");
  h2ExpectedPi      -> Write("h2ExpectedPi");
  h2ExpectedKa      -> Write("h2ExpectedKa");
  h2ExpectedPr      -> Write("h2ExpectedPr");
  h2ExpectedSigmaEl -> Write("h2ExpectedSigmaEl");
  h2ExpectedSigmaPi -> Write("h2ExpectedSigmaPi");
  h2ExpectedSigmaKa -> Write("h2ExpectedSigmaKa");
  h2ExpectedSigmaPr -> Write("h2ExpectedSigmaPr");
  if (!MCclosure){
      //     h2CleanEl         -> Write("h2CleanElectron");
    h2CleanPi         -> Write("h2CleanPion");
    h2CleanKa         -> Write("h2CleanKaon");
    h2CleanPr         -> Write("h2CleanProton");
  }
  antiProtonCutG    -> Write();
  protonCutG        -> Write();
  pionCutG          -> Write();
  hArmPod           -> Write("hArmPod");
     
  delete tmpFile;
}
//____________________________________________________________________________________________________________
void RealDataHistsFromThnSparse(const TString cleanFile, const TString histFile, const Double_t etaLow, const Double_t etaUp, const Float_t centLow, const Float_t centUp)
{
  
  //
  // Produce expected, real data, clean sample histograms 
  //
  /*
  /lustre/nyx/alice/users/marsland/alice-tpc-notes/JIRA/ATO-123/code/SAMPAesds.C
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Syst_LooseCuts_cRows_80_16EtaBin_mombin20MeV_largeDCAxy/testMacros
  aliroot -l
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C+
  TString d = "../mergedPeriods/CleanSampleTree.root"
  TString m = "../mergedPeriods/marsland_MCmoments.root"
  TString h = "../mergedPeriods/marsland_Hists.root"
  TString i = "../mergedPeriods/marsland_IdenTree.list"
  CreateAllTIDENHists(4,0,0,"","",d,m,h,0,0.1,30,40)      
   
  */
  Double_t BSF = 0.001; // Bin securing factor for ThnSparses when using setrange user. Otherwise bin edge problems occur
  InitInitials();
  std::cout << " eta Range  = " << etaLow  << " - " << etaUp; 
  std::cout << " cent Range = " << centLow << " - " << centUp << std::endl;
    
  // ================================
  // Read THnSparses 
  std::cout << " ========= make histograms from THnSparse ========= " << std::endl; 
  TFile fhist(histFile);
  TList * list  = (TList*)fhist.Get("cleanHists");
  THnSparse *fhndEdx     = (THnSparse*)list->FindObject("hdEdx");
  THnSparse *fhnExpected = (THnSparse*)list->FindObject("hExpected");
  //   THnSparse *fhnCleanEl  = (THnSparse*)list->FindObject("hCleanEl");
  THnSparse *fhnCleanKa  = (THnSparse*)list->FindObject("hCleanKa");
  //   THnSparse *fhnCleanDe  = (THnSparse*)list->FindObject("hCleanDe");
  TH2D *hArmPod          = (TH2D*)list->FindObject("hArmPod");

  // ================================
  // prepare dEdx histograms for full data sample
  fhndEdx->GetAxis(1)->SetRangeUser(centLow+BSF,centUp-BSF);   // centrality
  fhndEdx->GetAxis(2)->SetRangeUser(etaLow+BSF,etaUp-BSF);     // eta
  h2Dall = (TH2D*)fhndEdx->Projection(4,3); h2Dall->SetName("h2Dall");
  fhndEdx->GetAxis(0)->SetRangeUser(0.+BSF,2-BSF);   // centrality
  h2Dpos = (TH2D*)fhndEdx->Projection(4,3); h2Dpos->SetName("h2Dpos");
  fhndEdx->GetAxis(0)->SetRangeUser(-1+BSF,0-BSF);   // centrality
  h2Dneg = (TH2D*)fhndEdx->Projection(4,3); h2Dneg->SetName("h2Dneg"); 
 
  // ================================
  // Make projections form THnSparse --> Expecteds and clean electrons and clean kaons 
  fhnExpected->GetAxis(2)->SetRangeUser(centLow+BSF,centUp-BSF);   // centrality
  fhnExpected->GetAxis(3)->SetRangeUser(etaLow+BSF,etaUp-BSF);     // eta
  fhnExpected->GetAxis(5)->SetRangeUser(2.+BSF,60.-BSF);           // sigma
  fhnExpected->GetAxis(6)->SetRangeUser(25.+BSF,1000.-BSF);        // mean
  //   fhnCleanEl->GetAxis(1)->SetRangeUser(centLow+BSF,centUp-BSF);    // centrality
  //   fhnCleanEl->GetAxis(2)->SetRangeUser(etaLow+BSF,etaUp-BSF);      // eta
  fhnCleanKa->GetAxis(1)->SetRangeUser(centLow+BSF,centUp-BSF);    // centrality
  fhnCleanKa->GetAxis(2)->SetRangeUser(etaLow+BSF,etaUp-BSF);      // eta
  //   fhnCleanDe->GetAxis(1)->SetRangeUser(centLow+BSF,centUp-BSF);    // centrality
  //   fhnCleanDe->GetAxis(2)->SetRangeUser(etaLow+BSF,etaUp-BSF);      // eta
  
  // ================================
  // get each particle PID response
  fhnExpected->GetAxis(0)->SetRangeUser(0+BSF,1-BSF);              // electron
  TH2D * h2ExpectedEl      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedEl->SetName("h2ExpectedEl");
  TH2D * h2ExpectedSigmaEl = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaEl->SetName("h2ExpectedSigmaEl");
  fhnExpected->GetAxis(0)->SetRangeUser(1+BSF,2-BSF);  // pion
  TH2D * h2ExpectedPi      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedPi->SetName("h2ExpectedPi");
  TH2D * h2ExpectedSigmaPi = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaPi->SetName("h2ExpectedSigmaPi");
  fhnExpected->GetAxis(0)->SetRangeUser(2+BSF,3-BSF);  // kaon
  TH2D * h2ExpectedKa      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedKa->SetName("h2ExpectedKa");
  TH2D * h2ExpectedSigmaKa = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaKa->SetName("h2ExpectedSigmaKa");
  fhnExpected->GetAxis(0)->SetRangeUser(3+BSF,4-BSF);  // proton
  TH2D * h2ExpectedPr      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedPr->SetName("h2ExpectedPr");
  TH2D * h2ExpectedSigmaPr = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaPr->SetName("h2ExpectedSigmaPr");
  fhnExpected->GetAxis(0)->SetRangeUser(4+BSF,5-BSF);  // proton
  TH2D * h2ExpectedDe      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedDe->SetName("h2ExpectedDe");
  TH2D * h2ExpectedSigmaDe = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaDe->SetName("h2ExpectedSigmaDe");
  
  // ================================
  // Clean Kaon and Electron
  std::cout << " make clean Kaon Electron " << std::endl; 
  TH2D * h2CleanKa = (TH2D*)fhnCleanKa->Projection(4,3); h2CleanKa->SetName("h2CleanKa");
  //   TH2D * h2CleanDe = (TH2D*)fhnCleanDe->Projection(4,3); h2CleanDe->SetName("h2CleanDe");
  //   TH2D * h2CleanEl = (TH2D*)fhnCleanEl->Projection(4,3); h2CleanEl->SetName("h2CleanEl");

  // ================================
  std::cout << " ========= make clean pions and protons from TTree ========= " << std::endl; 
  // Clean Pion and Proton
  TH2D *h2CleanPi=NULL, *h2CleanPr=NULL, *h2CleanPiTight=NULL, *h2CleanPiTOF=NULL;
  TCut etaCentCut = Form("eta>=%f && eta<%f && cent>=%f && cent<%f" ,etaLow ,etaUp,centLow-rangeCleanPions,centUp+rangeCleanPions);
  TCut piTOFcut   = "abs(piTOFnSigma)<2.";
  TCut prTOFcut   = "abs(prTOFnSigma)<2.";
  TCut piK0cut    = "piFromK0";
  TCut piPixedcut = "!v0haspixel";
  
  TTree *cleanPionTree=NULL, *cleanProtonTree=NULL;
  if (!MCclosure){
     // Read TTrees
    TFile fClean(cleanFile);
    cleanPionTree   = (TTree*)fClean.Get("pionTree");
    cleanProtonTree = (TTree*)fClean.Get("protonTree");
    Double_t narmTreeEntries = (test) ? testEntries : cleanPionTree->GetEntries();
    h2CleanPi      = new TH2D("h2CleanPi","h2CleanPi",ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
    h2CleanPiTight = new TH2D("h2CleanPiTight","h2CleanPiTight",ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
    h2CleanPiTOF   = new TH2D("h2CleanPiTOF","h2CleanPiTOF",ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
    h2CleanPr      = new TH2D("h2CleanPr","h2CleanPr",ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
    
    cleanPionTree   -> Project("h2CleanPi"      ,"dEdx:ptot",etaCentCut+piK0cut,"",narmTreeEntries);
    cleanPionTree   -> Project("h2CleanPiTight" ,"dEdx:ptot",etaCentCut+piK0cut+piPixedcut,"",narmTreeEntries);
    cleanPionTree   -> Project("h2CleanPiTOF"   ,"dEdx:ptot",etaCentCut+piTOFcut,"",narmTreeEntries);

    cleanProtonTree -> Project("h2CleanPr" ,"dEdx:ptot",etaCentCut+prTOFcut,"",narmTreeEntries);
  }
  
  // ================================
  // Write hists to the streamer
  outStream->GetFile()->cd();
  h2ExpectedEl      -> Write("h2ExpectedEl");
  h2ExpectedPi      -> Write("h2ExpectedPi");
  h2ExpectedKa      -> Write("h2ExpectedKa");
  h2ExpectedPr      -> Write("h2ExpectedPr");
  h2ExpectedDe      -> Write("h2ExpectedDe");
  h2ExpectedSigmaEl -> Write("h2ExpectedSigmaEl");
  h2ExpectedSigmaPi -> Write("h2ExpectedSigmaPi");
  h2ExpectedSigmaKa -> Write("h2ExpectedSigmaKa");
  h2ExpectedSigmaPr -> Write("h2ExpectedSigmaPr");
  h2ExpectedSigmaDe -> Write("h2ExpectedSigmaDe");
  if (!MCclosure){
    h2Dall            -> Write("h2Dall"); 
    h2Dpos            -> Write("h2Dpos");
    h2Dneg            -> Write("h2Dneg");
    //     h2CleanEl         -> Write("h2CleanElectron");
    //     if (h2CleanPiTight) h2CleanPiTight-> Write("h2CleanPionTight");
    //     if (h2CleanPi)      h2CleanPi     -> Write("h2CleanPion");
    //     if (h2CleanPiTOF)   h2CleanPiTOF  -> Write("h2CleanPionTOF");
    //     if (h2CleanKa)      h2CleanKa     -> Write("h2CleanKaon");
    //     if (h2CleanPr)      h2CleanPr     -> Write("h2CleanProton");
    //     if (h2CleanDe)      h2CleanDe     -> Write("h2CleanDeuteron");
    //     if (hArmPod)        hArmPod       -> Write("hArmPod");
    antiProtonCutG    -> Write();
    protonCutG        -> Write();
    pionCutG          -> Write();
    
  }
  
     
}
//____________________________________________________________________________________________________________
void IdenDataHists(const TString idenTreeFile, const TString dataTreeFile, const TString histFile, const Double_t etaLow, const Double_t etaUp, const Float_t centLow, const Float_t centUp)
{
  
  //
  // Produce expected, real data, clean sample histograms 
  //
  Int_t idenEtaBin  = hEta ->FindBin(etaLow)-1;
  Int_t idenCentBin = (centLow>=0) ? hCent->FindBin(centLow)-1 : hCent->FindBin(centLow+2)-1 ;  // cent<0 --> pp data
  
  if(!hEta) InitInitials();
  std::cout << " eta Range   = " << etaLow  << " - " << etaUp << std::endl; 
  std::cout << " cent Range  = " << centLow << " - " << centUp << std::endl;
  std::cout << " idenEtaBin  = " << idenEtaBin << std::endl; 
  std::cout << " idenCentBin = " << idenCentBin << std::endl;
  
  // Read TTrees
  GetIdenTree(idenTreeFile, dataTreeFile);
  
  Double_t nTreeEntriesAll=treeIden->GetEntries();
  if (nTreeEntriesAll<10) { std::cout << " === upss data tree is empty === " << std::endl; return; }
  Double_t nTreeEntries = (test) ? testEntries : nTreeEntriesAll;

  // ======= Get the copy tree ======
  std::cout << " tree Size is --> " << nTreeEntriesAll << std::endl;
  TFile * tmpFile = new TFile(Form("tmp_eta_%d_cent_%d.root",idenEtaBin,idenCentBin),"recreate");
  TString etaCentString = Form("myBin[0]==%d && myBin[1]==%d",idenEtaBin,idenCentBin);
  TString etaCutArmtree = Form("eta>=%f   && eta<%f" ,etaLow ,etaUp);
  Long64_t nAll = (test) ? (Long64_t)testEntries : 100000000000;
  if (treeIden   ->GetEntries()<10) { std::cout << " === upss main tree is empty === " << std::endl; return; }
  TTree * treeRestricted = treeIden->CopyTree(etaCentString.Data(),"",nAll,0);
  treeRestricted->Write();    treeIden->Delete();    
  // Tree for clean samples it is not filled for MC closure test
  TTree * treeRestrictedArm = armtree->CopyTree(etaCutArmtree.Data(),"",nAll,0);
  if (armtree->GetEntries()>10) { 
    std::cout << " === armPod tree is OK === " << std::endl; 
    treeRestrictedArm->Write(); armtree->Delete();
  }

  // ================================
  // 2D dEdx plots
  std::cout << " ========= make 2D dEdx plots ========= " << std::endl; 
  TCut parPos     = "sign>0.";  
  TCut parNeg     = "sign<0.";
  TStopwatch timer; timer.Start();
  nTreeEntries = treeRestricted->GetEntries();
  std::cout << " make 2D dEdx plot ----- total Entry = " << nTreeEntries << std::endl;
  
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // 2D hist projection
  //   h2Dall  = new TH2D("h2Dall","h2Dall"  ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  //   h2Dpos  = new TH2D("h2Dpos","h2Dpos"  ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  //   h2Dneg  = new TH2D("h2Dneg","h2Dneg"  ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  //   treeRestricted ->Project("h2Dall" ,"myDeDx:myBin[2]" , ""     ,"",   nTreeEntries);
  //   treeRestricted ->Project("h2Dpos" ,"myDeDx:myBin[2]" , parPos ,"",2.*nTreeEntries);
  //   treeRestricted ->Project("h2Dneg" ,"myDeDx:myBin[2]" , parNeg ,"",2.*nTreeEntries);
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // fill iden tree style
  //   Double_t myDeDx;
  //   Int_t sign;
  //   Int_t myBin[4];
  //   treeRestricted -> SetBranchAddress("sign",&sign);
  //   treeRestricted -> SetBranchAddress("myBin",myBin);
  //   treeRestricted -> SetBranchAddress("myDeDx",&myDeDx);
  //   for(Int_t i = 0; i < nTreeEntries; ++i)
  //   {
  //     treeRestricted -> GetEntry(i);
  //     h2Dall -> Fill(hMom -> GetBinCenter(myBin[2]), myDeDx);
  //     if (sign>0) h2Dpos -> Fill(hMom -> GetBinCenter(myBin[2]), myDeDx);
  //     if (sign<0) h2Dneg -> Fill(hMom -> GetBinCenter(myBin[2]), myDeDx);
  //   }
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Bin by bin create all 1D histograms
  
  std::cout << " =============== entries to be used in the end = " << nTreeEntries << "  ======================  " << std::endl;
  TObjArray hArrAll(ptNbins);  hArrAll.SetOwner(kTRUE); 
  TObjArray hArrPos(ptNbins);  hArrPos.SetOwner(kTRUE); 
  TObjArray hArrNeg(ptNbins);  hArrNeg.SetOwner(kTRUE); 
  for (Int_t imombin=0; imombin<ptNbins; imombin++){
    TH1F *hAll = new TH1F(Form("histAll_eta%d_cent%d_mom%d",idenEtaBin,idenCentBin,imombin),"hAll"  ,dEdxNbins,dEdxMin,dEdxMax);
    TH1F *hPos = new TH1F(Form("histPos_eta%d_cent%d_mom%d",idenEtaBin,idenCentBin,imombin),"hPos"  ,dEdxNbins,dEdxMin,dEdxMax);
    TH1F *hNeg = new TH1F(Form("histNeg_eta%d_cent%d_mom%d",idenEtaBin,idenCentBin,imombin),"hNeg"  ,dEdxNbins,dEdxMin,dEdxMax);
    TString cutAll = Form("myBin[0]==%d && myBin[1]==%d && myBin[2]==%d",idenEtaBin,idenCentBin,imombin);
    TString cutPos = Form("myBin[0]==%d && myBin[1]==%d && myBin[2]==%d && sign>0",idenEtaBin,idenCentBin,imombin);
    TString cutNeg = Form("myBin[0]==%d && myBin[1]==%d && myBin[2]==%d && sign<0",idenEtaBin,idenCentBin,imombin);
    treeRestricted ->Project(hAll->GetName() ,"myDeDx" , cutAll ,"",   nTreeEntries);
    treeRestricted ->Project(hPos->GetName() ,"myDeDx" , cutPos ,"",2.*nTreeEntries);
    treeRestricted ->Project(hNeg->GetName() ,"myDeDx" , cutNeg ,"",2.*nTreeEntries);
    hAll->SetLineColor(kBlack);
    hPos->SetLineColor(kBlue);
    hNeg->SetLineColor(kGreen+3);
    hAll->GetXaxis()->SetTitle("TPC d#it{E}/d#it{x} Signal (a.u.)");
    hPos->GetXaxis()->SetTitle("TPC d#it{E}/d#it{x} Signal (a.u.)");
    hNeg->GetXaxis()->SetTitle("TPC d#it{E}/d#it{x} Signal (a.u.)");
    hAll->GetYaxis()->SetTitle("entries");
    hPos->GetYaxis()->SetTitle("entries");
    hNeg->GetYaxis()->SetTitle("entries");

    hArrAll.AddAt(hAll,imombin);
    hArrPos.AddAt(hPos,imombin);
    hArrNeg.AddAt(hNeg,imombin);
    std::cout << "   EtaBin = " << idenEtaBin << "   CentBin = " << idenCentBin << "   Mombin = " << imombin << std::endl;
  }
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  timer.Stop(); timer.Print();
 
  std::cout << " ========= make PID response plots ========= " << std::endl; 
  timer.Start();
  // Read THnSparses 
  TFile fhist(histFile);
  TList * list  = (TList*)fhist.Get("cleanHists");
  THnSparse *fhnExpected = (THnSparse*)list->FindObject("fhnExpected");
  //   THnSparse *fhnCleanEl  = (THnSparse*)list->FindObject("fhnCleanEl");
  THnSparse *fhnCleanKa  = (THnSparse*)list->FindObject("fhnCleanKa");
  TH2D *hArmPod          = (TH2D*)list->FindObject("fHistArmPod");

  // Make projections form THnSparse
  fhnExpected->GetAxis(2)->SetRangeUser(centLow,centUp);   // centrality
  fhnExpected->GetAxis(3)->SetRangeUser(etaLow,etaUp);     // eta
  fhnExpected->GetAxis(5)->SetRangeUser(2.,60.);          // sigma
  fhnExpected->GetAxis(6)->SetRangeUser(25.,1000.);        // mean
  //   fhnCleanEl->GetAxis(1)->SetRangeUser(centLow,centUp);    // centrality
  //   fhnCleanEl->GetAxis(2)->SetRangeUser(etaLow,etaUp);      // eta
  fhnCleanKa->GetAxis(1)->SetRangeUser(centLow,centUp);    // centrality
  fhnCleanKa->GetAxis(2)->SetRangeUser(etaLow,etaUp);      // eta
  
  // get each particle PID response
  fhnExpected->GetAxis(0)->SetRangeUser(0,1);              // electron
  TH2D * h2ExpectedEl      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedEl->SetName("h2ExpectedEl");
  TH2D * h2ExpectedSigmaEl = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaEl->SetName("h2ExpectedSigmaEl");
  fhnExpected->GetAxis(0)->SetRangeUser(1,2);  // pion
  TH2D * h2ExpectedPi      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedPi->SetName("h2ExpectedPi");
  TH2D * h2ExpectedSigmaPi = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaPi->SetName("h2ExpectedSigmaPi");
  fhnExpected->GetAxis(0)->SetRangeUser(2,3);  // kaon
  TH2D * h2ExpectedKa      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedKa->SetName("h2ExpectedKa");
  TH2D * h2ExpectedSigmaKa = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaKa->SetName("h2ExpectedSigmaKa");
  fhnExpected->GetAxis(0)->SetRangeUser(3,4);  // proton
  TH2D * h2ExpectedPr      = (TH2D*)fhnExpected->Projection(6,4); h2ExpectedPr->SetName("h2ExpectedPr");
  TH2D * h2ExpectedSigmaPr = (TH2D*)fhnExpected->Projection(5,4); h2ExpectedSigmaPr->SetName("h2ExpectedSigmaPr");
  
  std::cout << " make clean Kaon Electron " << std::endl; 
  // Clean Kaon and Electron
  //   TH2D * h2CleanEl = (TH2D*)fhnCleanEl->Projection(4,3); h2CleanEl->SetName("h2CleanEl");
  TH2D * h2CleanKa = (TH2D*)fhnCleanKa->Projection(4,3); h2CleanKa->SetName("h2CleanKa");
  timer.Stop(); timer.Print();
  
  // Clean Pion and Proton
  TH2D *h2CleanPi=NULL, *h2CleanPr=NULL;
  if (!MCclosure){
    Double_t narmTreeEntries = (test) ? testEntries : treeRestrictedArm->GetEntries();
    std::cout << " make clean Pion Proton ---- total Entry = " << treeRestrictedArm->GetEntries() << std::endl;
    timer.Start();
    h2CleanPi = new TH2D("h2CleanPi","h2CleanPi",ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
    h2CleanPr = new TH2D("h2CleanPr","h2CleanPr",ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
    
    TCut cleanCutPi = "abs(piTOFnSigma)<3. && abs(alfa)<0.5  && pionCutG";
    TCut cleanCutPr = "abs(prTOFnSigma)<3. && qt>0.07 && qt<0.11 && (antiProtonCutG||protonCutG)";
    treeRestrictedArm -> Project("h2CleanPi" ,"dEdx:ptot",cleanCutPi,"",narmTreeEntries);
    treeRestrictedArm -> Project("h2CleanPr" ,"dEdx:ptot",cleanCutPr,"",narmTreeEntries);
    timer.Stop(); timer.Print();
  }
  
  // Write hists to the streamer
  outStream->GetFile()->cd();
//   h2Dall            -> Write("h2Dall"); 
//   h2Dpos            -> Write("h2Dpos");
//   h2Dneg            -> Write("h2Dneg");
  hArrAll           . Write("hArrAll",TObject::kSingleKey);  hArrAll.Delete();
  hArrPos           . Write("hArrPos",TObject::kSingleKey);  hArrPos.Delete();
  hArrNeg           . Write("hArrNeg",TObject::kSingleKey);  hArrNeg.Delete();
  h2ExpectedEl      -> Write("h2ExpectedEl");
  h2ExpectedPi      -> Write("h2ExpectedPi");
  h2ExpectedKa      -> Write("h2ExpectedKa");
  h2ExpectedPr      -> Write("h2ExpectedPr");
  h2ExpectedSigmaEl -> Write("h2ExpectedSigmaEl");
  h2ExpectedSigmaPi -> Write("h2ExpectedSigmaPi");
  h2ExpectedSigmaKa -> Write("h2ExpectedSigmaKa");
  h2ExpectedSigmaPr -> Write("h2ExpectedSigmaPr");
  if (!MCclosure){
      //     h2CleanEl         -> Write("h2CleanElectron");
    h2CleanPi         -> Write("h2CleanPion");
    h2CleanKa         -> Write("h2CleanKaon");
    h2CleanPr         -> Write("h2CleanProton");
  }
  antiProtonCutG    -> Write();
  protonCutG        -> Write();
  pionCutG          -> Write();
  hArmPod           -> Write("hArmPod");
  
  delete tmpFile;
}
//____________________________________________________________________________________________________________
void MCDataHists(const TString mcFile, const Double_t etaLow, const Double_t etaUp, const Float_t centLow, const Float_t centUp)
{
  
  //
  // Produce MC sample histograms 
  //
    
  // Eta and centrality cuts  
  TCut centCutMC = Form("cent>=%f && cent<%f" ,centLow,centUp);
  TCut etaCutMC  = Form("eta>=%f  && eta<%f"  ,etaLow ,etaUp);
  TCut parPos    = "sign>0.";
  TCut parNeg    = "sign<0.";
  TCut elCut = "el>0", piCut = "pi>0", kaCut = "ka>0", prCut = "pr>0", allCut = "dEdx>0";
  
  
  TTree *tree=NULL; 
  if (mcFile.Contains(".root")){
      TFile * fInputFile = TFile::Open(mcFile);
      tree               = (TTree*)fInputFile->Get("fTreeMC");
  } else {
      tree=(TTree*)AliXRDPROOFtoolkit::MakeChain(mcFile, "fTreeMC", 0, 1000000000,0,1);
  }

  if (tree->GetEntries()<10) { 
    std::cout << " === upss MC tree is empty === " << std::endl; return;
  }
  Double_t mcTreeEntries = (test) ? testEntries : tree->GetEntries();
  
   // Histograms
  TH2D * h2DMCall     = new TH2D("h2DMCall"    ,"h2DMCall"     ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2DMCpos     = new TH2D("h2DMCpos"    ,"h2DMCpos"     ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2DMCneg     = new TH2D("h2DMCneg"    ,"h2DMCneg"     ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCelectron = new TH2D("h2MCelectron","h2MCelectron" ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCpion     = new TH2D("h2MCpion"    ,"h2MCpion"     ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCkaon     = new TH2D("h2MCkaon"    ,"h2MCkaon"     ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCproton   = new TH2D("h2MCproton"  ,"h2MCproton"   ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  
  TH2D * h2MCelectronNeg = new TH2D("h2MCelectronNeg","h2MCelectronNeg" ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCpionNeg     = new TH2D("h2MCpionNeg"    ,"h2MCpionNeg"     ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCkaonNeg     = new TH2D("h2MCkaonNeg"    ,"h2MCkaonNeg"     ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCprotonNeg   = new TH2D("h2MCprotonNeg"  ,"h2MCprotonNeg"   ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
 
  TH2D * h2MCelectronPos = new TH2D("h2MCelectronPos","h2MCelectronPos" ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCpionPos     = new TH2D("h2MCpionPos"    ,"h2MCpionPos"     ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCkaonPos     = new TH2D("h2MCkaonPos"    ,"h2MCkaonPos"     ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
  TH2D * h2MCprotonPos   = new TH2D("h2MCprotonPos"  ,"h2MCprotonPos"   ,ptNbins,ptMin,ptMax,dEdxNbins,dEdxMin,dEdxMax);
 
  
  std::cout << " make MC sample plots ---- total Entry =  " << tree->GetEntries() << std::endl;
  // Projections
  tree->Project("h2DMCall"        ,"dEdx:ptot",allCut && etaCutMC && centCutMC,"",mcTreeEntries);
  tree->Project("h2DMCpos"        ,"dEdx:ptot",allCut && etaCutMC && centCutMC && parPos,"",mcTreeEntries);
  tree->Project("h2DMCneg"        ,"dEdx:ptot",allCut && etaCutMC && centCutMC && parNeg,"",mcTreeEntries);
  tree->Project("h2MCelectron"    ,"el:ptot"  ,elCut  && etaCutMC && centCutMC,"",mcTreeEntries);
  tree->Project("h2MCpion"        ,"pi:ptot"  ,piCut  && etaCutMC && centCutMC,"",mcTreeEntries);
  tree->Project("h2MCkaon"        ,"ka:ptot"  ,kaCut  && etaCutMC && centCutMC,"",mcTreeEntries);
  tree->Project("h2MCproton"      ,"pr:ptot"  ,prCut  && etaCutMC && centCutMC,"",mcTreeEntries);
  
  tree->Project("h2MCelectronNeg" ,"el:ptot"  ,elCut  && etaCutMC && centCutMC && parNeg,"",mcTreeEntries);
  tree->Project("h2MCpionNeg"     ,"pi:ptot"  ,piCut  && etaCutMC && centCutMC && parNeg,"",mcTreeEntries);
  tree->Project("h2MCkaonNeg"     ,"ka:ptot"  ,kaCut  && etaCutMC && centCutMC && parNeg,"",mcTreeEntries);
  tree->Project("h2MCprotonNeg"   ,"pr:ptot"  ,prCut  && etaCutMC && centCutMC && parNeg,"",mcTreeEntries);
 
  tree->Project("h2MCelectronPos" ,"el:ptot"  ,elCut  && etaCutMC && centCutMC && parPos,"",mcTreeEntries);
  tree->Project("h2MCpionPos"     ,"pi:ptot"  ,piCut  && etaCutMC && centCutMC && parPos,"",mcTreeEntries);
  tree->Project("h2MCkaonPos"     ,"ka:ptot"  ,kaCut  && etaCutMC && centCutMC && parPos,"",mcTreeEntries);
  tree->Project("h2MCprotonPos"   ,"pr:ptot"  ,prCut  && etaCutMC && centCutMC && parPos,"",mcTreeEntries);
 
  // Write hists to the streamer
  outStream->GetFile()->cd();
  h2DMCall     -> Write("h2Dall");
  h2DMCpos     -> Write("h2Dpos");
  h2DMCneg     -> Write("h2Dneg");
  h2MCelectron -> Write("h2MCElectron");
  h2MCpion     -> Write("h2MCPion");
  h2MCkaon     -> Write("h2MCKaon");
  h2MCproton   -> Write("h2MCProton");
  h2MCelectronPos -> Write("h2MCElectronPos");
  h2MCpionPos     -> Write("h2MCPionPos");
  h2MCkaonPos     -> Write("h2MCKaonPos");
  h2MCprotonPos   -> Write("h2MCProtonPos");
  h2MCelectronNeg -> Write("h2MCElectronNeg");
  h2MCpionNeg     -> Write("h2MCPionNeg");
  h2MCkaonNeg     -> Write("h2MCKaonNeg");
  h2MCprotonNeg   -> Write("h2MCProtonNeg");
  std::cout << " ================= " << std::endl; 
   
}
//____________________________________________________________________________________________________________
void MakeSubSamples(TString idenTreeFile,TString treeName,Int_t centBin,Int_t subsample)
{
  
  /*
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_Reference/SubSamples
  aliroot -l
  .L ~/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C 
  TString iden = "../mergedPeriods/marsland_IdenTree.list"
  MakeSubSamples(iden,"fIdenTree",3,5)  
  */
  
  // ================================
  // get the old tree
  TString ssFileName = Form("SubSample_cent%3.2f_ss%d.root",centArray[centBin],subsample);
  TTree * treeId;  
  if (idenTreeFile.Contains(".root")){
    TFile *oldfile = TFile::Open(idenTreeFile);
    treeId  = (TTree*)oldfile->Get(treeName); 
  } else {
    TChain *chain = AliXRDPROOFtoolkit::MakeChainRandom(idenTreeFile,treeName,0,500000,0);
    treeId = (TTree*)chain;
  }
  // ================================
  //
  // ================================
  if (treeId->GetEntries()<10) { std::cout << " === upss iden treeId is empty === " << std::endl; return; }
  Double_t sampTreeEntries = (test) ? testEntries : treeId->GetEntries()/(nSubsample-1);
  // ================================
  //
  // ================================
  // Create new file to put subsample treeId
  TFile * newfile = new TFile(ssFileName,"recreate");
  Long64_t firstEntry = sampTreeEntries*(subsample-1);  
  std::cout << " firstEntry = " << firstEntry << "  step = " << sampTreeEntries << std::endl;
  TTree * treeIdSample = treeId->CopyTree("1","1",sampTreeEntries,firstEntry);
  treeIdSample->Write(treeName);
  delete treeId;
  delete newfile;
                 
}
//____________________________________________________________________________________________________________
void MakeMainSubSample(TString idenTreeFile,TString treeName,Int_t centBin)
{
  
  /*
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C
  TString iden = "../mergedPeriods/marsland_IdenTree.list"
  MakeMainSubSample(iden,"fIdenTree",1)  
  */
  
  // ================================
  // get the old tree
  TTree *treeID;
  if (idenTreeFile.Contains(".root")){ 
    TFile *oldfile = TFile::Open(idenTreeFile); 
    treeID = (TTree*)oldfile->Get(treeName); 
  } else { 
    treeID = (TTree*)AliXRDPROOFtoolkit::MakeChain(idenTreeFile,treeName,0,-1);
  }
  // ================================
  Double_t sampTreeEntries = (test) ? testEntries : treeID->GetEntries();
  std::cout << " sampTreeEntries " << sampTreeEntries << "   ---   centrality = " << centArray[centBin]<< std::endl;
  // ================================
  //
  // ================================
  // Create new file to put subsample tree
  TFile * outfile = new TFile(Form("SubSample_cent%3.2f_ss0.root",centArray[centBin]),"recreate");
  TTree * treeIdSample = treeID->CopyTree(Form("myBin[1]==%d",centBin),"",sampTreeEntries*10.);
  treeIdSample->Write(treeName);
  delete treeID;
  delete outfile;
                 
}
//____________________________________________________________________________________________________________
void MergeFitResultsTree(TString infile)
{
  
  //
  // Obtain ttree with sequentially ordered events
  // infile   : merged tree of fits form PIDIterativeFitting.C
  /*
  
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/MC/MC_PbPb_cRows/ParamTrees
  aliroot -l 
  .L ~/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C+ 
  MergeFitResultsTree("AllPIDtrees.root")
   
  */
  //
  gSystem->Exec("cp /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C .");
    
  // Initialise Histograms for the myBin[] varibale
  if(!hEta) InitInitials();
  
  if (infile.Contains(".root")){
    TFile *f = TFile::Open(infile);
    tree = (TTree*)f->Get("IdMethodInput");  
  } else {
    TChain *chain = AliXRDPROOFtoolkit::MakeChainRandom(infile,"IdMethodInput",0,500000,0);
    //     chain -> SetCacheSize(100000000000000000);
    tree = (TTree*)chain;
  }
 
  Double_t nEntries = tree->GetEntries();
  
  Int_t myBin[3]      = {0,0,0};
  Double_t totChi2    = 0.;
  Double_t ndf        = 0.; 
  Int_t sl            = 0;   
  Int_t it            = 0;

  Int_t iter          = 0;
  Int_t slice         = 0;    
  Double_t p          = 0.;
  Float_t eta         = 0.;
  Float_t cent          = 0;
  Double_t totalChi2  = 0.;
  Double_t normNDF    = 0.; 
       
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
               
  Double_t elKurtosis = 0., elK = 0.;
  Double_t piKurtosis = 0., piK = 0.;
  Double_t kaKurtosis = 0., kaK = 0.;
  Double_t prKurtosis = 0., prK = 0.;
        
  Double_t elSkew     = 0., elSk = 0.;
  Double_t piSkew     = 0., piSk = 0.;
  Double_t kaSkew     = 0., kaSk = 0.;
  Double_t prSkew     = 0., prSk = 0.;
   
  tree->SetBranchAddress("iter"      ,&iter);
  tree->SetBranchAddress("slice"     ,&slice);
  tree->SetBranchAddress("p"         ,&p);
  tree->SetBranchAddress("eta"       ,&eta);
  tree->SetBranchAddress("cent"      ,&cent);
  tree->SetBranchAddress("totalChi2" ,&totalChi2);
  tree->SetBranchAddress("normNDF"   ,&normNDF);
  
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

  // New Tree
  TString paramtreefile = "ParamTree.root";
  TFile f(paramtreefile,"recreate");
  TTree *treeId = new TTree("treeId","IdentityInput FitParamTree");
  treeId->Branch("it"      ,&it       ,"it/I");
  treeId->Branch("sl"      ,&sl       ,"sl/I");
  treeId->Branch("myBin"   ,myBin     ,"myBin[3]/I");
  treeId->Branch("totChi2" ,&totChi2  ,"totChi2/D");
  treeId->Branch("ndf"     ,&ndf      ,"ndf/D");
  
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
  
  std::cout << nEntries << "  number of entries will be processed " << std::endl; 
  for (Int_t itrack=0; itrack<nEntries; ++itrack) {
   
    tree->GetEntry(itrack);
//     Int_t tmpIter = iter;
    
    // Parameter boundary conditions --> Scan of the previous and next slices on case there is a fit fail
    Bool_t elParamCheck = ((elMean>20 && elMean<1020) && (elSigma>0 && elSigma<200) && (elAmp>=1 && elAmp<1e+30) && (elSkew>=0 && elSkew<3) && (elKurtosis>1.5 && elKurtosis<4.));
    Bool_t piParamCheck = ((piMean>20 && piMean<1020) && (piSigma>0 && piSigma<200) && (piAmp>=1 && piAmp<1e+30) && (piSkew>=0 && piSkew<3) && (piKurtosis>1.5 && piKurtosis<4.));
    Bool_t kaParamCheck = ((kaMean>20 && kaMean<1020) && (kaSigma>0 && kaSigma<200) && (kaAmp>=1 && kaAmp<1e+30) && (kaSkew>=0 && kaSkew<3) && (kaKurtosis>1.5 && kaKurtosis<4.));
    Bool_t prParamCheck = ((prMean>20 && prMean<1020) && (prSigma>0 && prSigma<200) && (prAmp>=1 && prAmp<1e+30) && (prSkew>=0 && prSkew<3) && (prKurtosis>1.5 && prKurtosis<4.));
   
    /*
    // ================================= Set the fit Params =================================
    if (itrack%1000==0) //std::cout << " i = " << itrack << " ===== params will be set ==== " << std::endl;
    Int_t iterStep=5;
    
    // ============ Electron ============
    for (Int_t i = itrack; i<itrack+iterStep; i++) {
      tree->GetEntry(i);
      if (elParamCheck) { elM=elMean;  elSi=elSigma;  elA=elAmp;  elSk=elSkew;  elK=elKurtosis; break; }
    }
  
    if (!elParamCheck && tmpIter>20){
      for (Int_t i = itrack; i>0; i--) {
        tree->GetEntry(i);
        if (tmpIter!=iter) break;
        elM=elMean;  elSi=elSigma;  elA=elAmp;  elSk=elSkew;  elK=elKurtosis;
        if (elParamCheck) break;
      }
    }
    
    // ============ Pion ============
    for (Int_t i = itrack; i<itrack+iterStep; i++) {
      tree->GetEntry(i);
      if (piParamCheck) { piM=piMean;  piSi=piSigma;  piA=piAmp;  piSk=piSkew;  piK=piKurtosis; break; }
    }
  
    if (!piParamCheck && tmpIter>20){
      for (Int_t i = itrack; i>0; i--) {
        tree->GetEntry(i);
        if (tmpIter!=iter) break;
        piM=piMean;  piSi=piSigma;  piA=piAmp;  piSk=piSkew;  piK=piKurtosis;
        if (piParamCheck) break;
      }
    }
    
    // ============ Kaon ============
    for (Int_t i = itrack; i<itrack+iterStep; i++) {
      tree->GetEntry(i);
      if (kaParamCheck) { kaM=kaMean;  kaSi=kaSigma;  kaA=kaAmp;  kaSk=kaSkew;  kaK=kaKurtosis; break; }
    }
  
    if (!kaParamCheck && tmpIter>20){
      for (Int_t i = itrack; i>0; i--) {
        tree->GetEntry(i);
        if (tmpIter!=iter) break;
        kaM=kaMean;  kaSi=kaSigma;  kaA=kaAmp;  kaSk=kaSkew;  kaK=kaKurtosis;
        if (kaParamCheck) break;
      }
    }
    
    // ============ Proton ============
    for (Int_t i = itrack; i<itrack+iterStep; i++) {
      tree->GetEntry(i);
      if (prParamCheck) { prM=prMean;  prSi=prSigma;  prA=prAmp;  prSk=prSkew;  prK=prKurtosis; break; }
    }
  
    if (!prParamCheck && tmpIter>20){
      for (Int_t i = itrack; i>0; i--) {
        tree->GetEntry(i);
        if (tmpIter!=iter) break;
        prM=prMean;  prSi=prSigma;  prA=prAmp;  prSk=prSkew;  prK=prKurtosis;
        if (prParamCheck) break;
      }
    }
    // =====================================================================================
    */ 
     
    elM=elMean;  elSi=elSigma;  elA=elAmp;  elSk=elSkew;  elK=elKurtosis;
    piM=piMean;  piSi=piSigma;  piA=piAmp;  piSk=piSkew;  piK=piKurtosis;
    kaM=kaMean;  kaSi=kaSigma;  kaA=kaAmp;  kaSk=kaSkew;  kaK=kaKurtosis;
    prM=prMean;  prSi=prSigma;  prA=prAmp;  prSk=prSkew;  prK=prKurtosis;
    
    // avoid from failed fits lower limit
    if ( !elParamCheck ) {elM=0.;elSigma=1.;elA=0.;elSk=0.;elK=2.;}
    if ( !piParamCheck ) {piM=0.;piSigma=1.;piA=0.;piSk=0.;piK=2.;}
    if ( !kaParamCheck ) {kaM=0.;kaSigma=1.;kaA=0.;kaSk=0.;kaK=2.;}
    if ( !prParamCheck ) {prM=0.;prSigma=1.;prA=0.;prSk=0.;prK=2.;}

    // Set eta momentum cent etc. and fill the treeId
    myBin[0] = hEta ->FindBin(eta+0.0001)-1;
    myBin[1] = hCent->FindBin(cent+0.0001)-1;
    myBin[2] = hMom->FindBin(p+0.0001)-1;
    totChi2=totalChi2; ndf=normNDF; it = iter; sl = slice;
    treeId->Fill();
    
  }
  
  f.cd();
  treeId->Write();
  hEta->Write();
  hCent->Write();
  hMom->Write();

}
//____________________________________________________________________________________________________________
void CheckResults(TString paramFile)
{
  
  //
  // Look at all fit params for the smoothness of parameters vs momentum
  //
  
  /*
  
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/PbPb_cRows_80_mombin20MeV/ParamTrees
  aliroot -l 
  .L ~/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C+ 
  CheckResults("ParamTree_SemiAutoKS_properMomBinning.root")
   
  */
  //
  
  TTreeSRedirector *streamer = new TTreeSRedirector("ParameterControl.root","recreate");
  
  TFile *fResControl    = TFile::Open(paramFile);
  TTree *treeResControl = (TTree*)fResControl->Get("treeId");  
    
  // loop over particle types
  for (Int_t iType=0; iType<3; iType++){
    
    TString parType; 
    if (iType==0) parType = "pi";  
    if (iType==1) parType = "ka";  
    if (iType==2) parType = "pr";
  
  // loop over all iterations 
    for (Int_t iIter=0;iIter<nIter;iIter++){
    
      //std::cout << "particle type = " << parType <<  "      iter = " << iIter << std::endl;
    
      TGraphErrors *grAmp[nEtabins][nCentbins];
      TGraphErrors *grMean[nEtabins][nCentbins];
      TGraphErrors *grSigma[nEtabins][nCentbins];
      TGraphErrors *grKurt[nEtabins][nCentbins];
      TGraphErrors *grSkew[nEtabins][nCentbins];
      TGraphErrors *grChi2[nEtabins][nCentbins];
    
      for (Int_t i=0; i<nEtabins; i++)
      {
        for (Int_t j=0; j<nCentbins; j++){
          grAmp[i][j]=NULL; 
          grMean[i][j]=NULL; 
          grSigma[i][j]=NULL; 
          grKurt[i][j]=NULL; 
          grSkew[i][j]=NULL;
          grChi2[i][j]=NULL;
        }
      }
    
    // get the graphs for each eta and cent bin
      for (Int_t iEta=0; iEta<nEtabins; iEta++)
      {
        for (Int_t iCent=0; iCent<nCentbins; iCent++)
        {
        
          TString mainCut   = Form("it==%d && myBin[1]==%d && myBin[0]==%d && myBin[2]<=140",iIter,iCent,iEta);
          TString drawAmp   = Form("%sA:myBin[2]" ,parType.Data());
          TString drawMean  = Form("%sM:myBin[2]" ,parType.Data());
          TString drawSigma = Form("%sSi:myBin[2]",parType.Data());
          TString drawKurt  = Form("%sK:myBin[2]" ,parType.Data());
          TString drawSkew  = Form("%sSk:myBin[2]",parType.Data());
          TString drawChi2  = "totChi2/ndf:myBin[2]";
      
          grSigma[iEta][iCent] = TStatToolkit::MakeGraphErrors(treeResControl,drawSigma.Data(),mainCut.Data(),20+iEta+1,iEta+1,0.5);
          grMean[iEta][iCent]  = TStatToolkit::MakeGraphErrors(treeResControl,drawMean.Data() ,mainCut.Data(),20+iEta+1,iEta+1,0.7);
          grAmp[iEta][iCent]   = TStatToolkit::MakeGraphErrors(treeResControl,drawAmp.Data()  ,mainCut.Data(),20+iEta+1,iEta+1,0.7);
          grKurt[iEta][iCent]  = TStatToolkit::MakeGraphErrors(treeResControl,drawKurt.Data() ,mainCut.Data(),20+iEta+1,iEta+1,0.7);
          grSkew[iEta][iCent]  = TStatToolkit::MakeGraphErrors(treeResControl,drawSkew.Data() ,mainCut.Data(),20+iEta+1,iEta+1,0.7);
          grChi2[iEta][iCent]  = TStatToolkit::MakeGraphErrors(treeResControl,drawChi2.Data() ,mainCut.Data(),20+iEta+1,iEta+1,0.7);

        } // end of cent loop
      } // end of eta loop
    
      TObjArray paramsPerCentArr(nCentbins); paramsPerCentArr.SetOwner(kTRUE);
      // plots the canvases
      for (Int_t iCent=0;iCent<nCentbins;iCent++){
      
        TCanvas *paramCan = new TCanvas(Form("paramCan_%s_iter%d_cent%d",parType.Data(),iIter,iCent), "All params", 1200, 800);   
        paramCan->Divide(3,2);
    
      // Sigmas
        paramCan->cd(1);
        grSigma[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grSigma[i][iCent]->Draw("p");
        TLegend *lSigma = new TLegend(0.25, 0.55, 0.85, 0.85);
        lSigma->SetTextFont(62); lSigma->SetTextSize(0.03);  lSigma->SetFillColor(0);  lSigma-> SetNColumns(2);
        lSigma->AddEntry(grSigma[0][iCent]," Eta=[-0.8,-0.6]","LP");
        lSigma->AddEntry(grSigma[7][iCent]," Eta=[0.6,0.8]"  ,"LP");
        lSigma->AddEntry(grSigma[1][iCent]," Eta=[-0.6,-0.4]","LP");
        lSigma->AddEntry(grSigma[6][iCent]," Eta=[0.4,0.6]"  ,"LP");
        lSigma->AddEntry(grSigma[2][iCent]," Eta=[-0.4,-0.2]","LP");
        lSigma->AddEntry(grSigma[5][iCent]," Eta=[0.2,0.4]"  ,"LP");
        lSigma->AddEntry(grSigma[3][iCent]," Eta=[-0.2,0]"   ,"LP");
        lSigma->AddEntry(grSigma[4][iCent]," Eta=[0,0.2]"    ,"LP");
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
      
        // Skewness
        paramCan->cd(5);
        grSkew[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grSkew[i][iCent]->Draw("p");
      
        // Skewness
        paramCan->cd(6);
        grChi2[0][iCent]->Draw("ap");
        for (Int_t i=1;i<nEtabins;i++) grChi2[i][iCent]->Draw("p");
      
        // fill the graphs in tobjarray
        paramsPerCentArr.AddAt(paramCan,iCent);   
        if (iIter==6) paramsPerCentArr.SaveAs(Form("paramCan_%s_iter%d_cent%d.png",parType.Data(),iIter,iCent));
       
      }
      
      streamer->GetFile()->cd();
      paramsPerCentArr.Write(Form("FitParams_%s_iter%d",parType.Data(),iIter) ,TObject::kSingleKey);  
    
    } // end of iter loop
  } // end of particle type loop
  
  delete streamer;
}
//____________________________________________________________________________________________________________
void InitInitials()
{
  
  // 
  // Initialise histograms to be used for binning
  //
  
  std::cout <<  "pp =                 " << pp                 << std::endl;
  std::cout <<  "test =               " << test               << std::endl;
  std::cout <<  "testEntries =        " << testEntries        << std::endl;
  std::cout <<  "nSubsample =         " << nSubsample         << std::endl; 
  std::cout <<  "ptNbins =            " << ptNbins            << std::endl;
  std::cout <<  "centArray =          " << centArray[0]       << std::endl;
  std::cout <<  "nCentbins =          " << nCentbins          << std::endl;
  std::cout <<  "nEtabins =           " << nEtabins           << std::endl;
  std::cout <<  "nIter =              " << nIter              << std::endl;
  
  // Create the histograms to be used in the binning of eta, cent and momentum
  hEta  =  new TH1D("hEta" ,"Eta Bins"       ,nEtabins   ,etaMin, etaMax );
  hMom  =  new TH1D("hMom" ,"Mom Bins"       ,ptNbins    ,ptMin,  ptMax );
  hCent =  new TH1D("hCent","Centrality Bins",nCentBinsPlusOne-1,xCentBins );
  hEta  -> FillRandom("gaus",1000);
  hCent -> FillRandom("gaus",1000);
  hMom  -> FillRandom("gaus",1000);
  
  // Define graphical cuts from armenteros podolanski space
  TFile cutGFile("/lustre/nyx/alice/users/marsland/pFluct/ArmPodGraphicalCuts_Tight.root");
  pionCutG       = (TCutG*)cutGFile.Get("pionCutG");
  antiProtonCutG = (TCutG*)cutGFile.Get("antiProtonCutG");
  protonCutG     = (TCutG*)cutGFile.Get("protonCutG");
  antiProtonCutG -> SetName("antiProtonCutG");
  antiProtonCutG -> SetVarX("alfa");
  antiProtonCutG -> SetVarY("qt");
  protonCutG     -> SetName("protonCutG");
  protonCutG     -> SetVarX("alfa");
  protonCutG     -> SetVarY("qt");
  pionCutG       -> SetName("pionCutG");
  pionCutG       -> SetVarX("alfa");
  pionCutG       -> SetVarY("qt");
               
}
//____________________________________________________________________________________________________________
void GetDeDxTree(TString list)
{

  //
  // Get the dEdxTree (i.e identity tree with branches myBin, mydEdx, ...) either from a root file or from a chain 
  // "list" can be either a list of root files containing dEdxTree or a single file 
  // 
  
  if (list.Contains(".root")) {
    TFile *f = TFile::Open(list);
    tree      = (TTree*)f->Get("fTree");
    armtree   = (TTree*)f->Get("fArmPodTree");
    std::cout << "get dEdxTree from ---------> " << list << std::endl;
  } else {
    tree    = (TTree*)AliXRDPROOFtoolkit::MakeChain(list,"fTree", 0,-1);
    armtree = (TTree*)AliXRDPROOFtoolkit::MakeChain(list,"fArmPodTree", 0,-1);
  }
  
}
//____________________________________________________________________________________________________________
void GetIdenTree(TString idenlist, TString datalist)
{

  //
  // Get the dEdxTree (i.e identity tree with branches myBin, mydEdx, ...) either from a root file or from a chain 
  // "list" can be either a list of root files containing dEdxTree or a single file 
  // 
  
  if (idenlist.Contains(".root")) {
    TFile *fiden = TFile::Open(idenlist);
    TFile *fdata = TFile::Open(idenlist);
    if (MCclosure) {
      treeIden = (TTree*)fiden->Get("fIdenTreeMC");
    } else {
      treeIden = (TTree*)fiden->Get("fIdenTree");
    }
    armtree  = (TTree*)fdata->Get("fArmPodTree");
    std::cout << "get dEdxTree from the root file ---------> " << idenlist << std::endl;
  } else {
    std::cout << "get dEdxTree from the list file ---------> " << idenlist << std::endl;
    if (MCclosure) {
      treeIden = (TTree*)AliXRDPROOFtoolkit::MakeChain(idenlist,"fIdenTreeMC", 0,-1);
    } else {
      treeIden = (TTree*)AliXRDPROOFtoolkit::MakeChain(idenlist,"fIdenTree", 0,-1);
    }
    armtree = (TTree*)AliXRDPROOFtoolkit::MakeChain(datalist,"fArmPodTree", 0,-1);
  }
  
}
//____________________________________________________________________________________________________________
void CreateCleanSamples(TString cleanSamplesTreeList)
{
  
  //
  // Produce expected, real data, clean sample histograms 
  //
  /*
  cd /lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/PbPb/Real/Syst_LooseCuts_cRows_80_16EtaBin_mombin20MeV_smallDCAxy/mergedPeriods
  aliroot -l
  .L /u/marsland/PHD/macros/marsland_EbyeRatios/CreateAllTIDENHists.C+
  CreateCleanSamples("marsland_DataTree.list")      
   
  */
  InitInitials();
  TCut momDeDx    = "ptot<3.1 && ptot>0.1 && dEdx>30 && dEdx<200";
  //   TCut cleanCutPi = "abs(piTOFnSigma)<2.  && abs(alfa)<0.5  && pionCutG";
  //   TCut cleanCutPr = "abs(prTOFnSigma)<2.  && qt>0.07 && qt<0.11 && (antiProtonCutG||protonCutG)";
  TCut cleanCutPi = "abs(alfa)<0.5 && pionCutG";
  TCut cleanCutPr = "qt>0.07 && qt<0.11 && (antiProtonCutG||protonCutG)";
  // ================================
  //
  // ================================
  if (cleanSamplesTreeList.Contains(".root")){ TFile oldfile(cleanSamplesTreeList); tree = (TTree*)oldfile.Get("fArmPodTree"); 
  } else { tree = (TTree*)AliXRDPROOFtoolkit::MakeChain(cleanSamplesTreeList,"fArmPodTree",0,-1); }
  
  //   TTree * tree = (TTree*)AliXRDPROOFtoolkit::MakeChain(cleanSamplesTreeList,"fArmPodTree", 0,-1);
  Double_t allEntries = tree->GetEntries();
  std::cout << " all entries = " << allEntries << std::endl;
  // ================================
  //
  // ================================
  // prepare copy trees for clean pion and proton
  TFile * outfile = new TFile("CleanSampleTree.root","recreate");
  TTree * pionTree   = tree->CopyTree(cleanCutPi && momDeDx,"",allEntries*2,0);
  TTree * protonTree = tree->CopyTree(cleanCutPr && momDeDx,"",allEntries*2,0);
  std::cout << " clean pion   entries = " << pionTree  ->GetEntries() << std::endl;
  std::cout << " clean proton entries = " << protonTree->GetEntries() << std::endl;
  pionTree  ->Write("pionTree"); 
  protonTree->Write("protonTree"); 
  delete tree;
  delete outfile;
    
}

