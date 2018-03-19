#include "iostream"
#include "TIdentity2D.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorF.h"
#include "TChain.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"
#include <iomanip>
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "string"

using namespace std;
using std::cout;
using std::setw;

Float_t effPion, effProton, effKaon;
enum particles{electron, pion, proton, kaon};
enum momentType {kPi=0,kKa=1,kPr=2,kPiPi=3,kKaKa=4,kPrPr=5,kPiKa=6,kPiPr=7,kKaPr=8,kLa=9,kLaLa=10,kCh=11,kChCh=12};
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
    kLaPosLaNeg=9,
    kChPosChNeg=10,
};

// ======= Modification Part =============================================================================
const Int_t fnSignBins      = 3;
const Int_t fnEtaBins       = 16;///16;    MC: 8, Data:16
const Float_t fEtaRangeDown = -0.8;///16;
const Float_t fEtaRangeUp   = 0.8;///16;
const Int_t fnMomBins       = 150;
const Float_t fMomRangeDown = 0.2;
const Float_t fMomRangeUp   = 3.2;
const Int_t fnCentBins  = 9;
Float_t     xCentBins[] = {0, 5,  10,  20, 30, 40, 50, 60, 70, 80};

Int_t charge           = -100;
Bool_t fpp             = kFALSE;
TString treeName       = "fEventTree";
TString treeLineShapes = "treeId";
TString treeIdentity   = "fIdenTree";
TString treeIdentityMC = "fIdenTreeMC";
TString fitFunctionGenGausStr = "[0]*exp(-(abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/1.414213))";


Int_t nBinsLineShape   = 4000;
const Int_t nParticles = 4;   // keep it always 4
const Int_t nCent      = 9;
const Int_t nEta       = 16;
const Int_t nMom      = 65;
const Int_t nSign      = 3; 

Bool_t test = kFALSE;
Bool_t lookUpTableForLine = kTRUE;
Int_t lookUpTableLineMode = 0;     // 0 for TH1D, 1 for TF1


// =======================================================================================================
Double_t amp[fnEtaBins][fnCentBins][fnMomBins][nParticles][2]; // amp[fEtaBin][fCentBin][fMomBin][particleType], before kurtosis[9][12][1000][6];
Double_t mean[fnEtaBins][fnCentBins][fnMomBins][nParticles][2];
Double_t sigma[fnEtaBins][fnCentBins][fnMomBins][nParticles][2];
Double_t skew[fnEtaBins][fnCentBins][fnMomBins][nParticles][2];
Double_t kurtosis[fnEtaBins][fnCentBins][fnMomBins][nParticles][2];
Int_t fEtaBin, fCentBin, fMomBin;
Int_t fCharge;
TStopwatch timer; 
// Values to set from outside
TFile *fLineShapesLookUpTable = NULL;
TClonesArray *cloneArrHist=NULL;
TClonesArray *cloneArrFunc=NULL;
TTree *treeLookUp = NULL;
TString inDir("/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/");
TString outDir("./");
TTree *momTree = NULL;
TH1D *fhPtot   = NULL;
TH1D *fhEta    = NULL;
TH1D *fhCent   = NULL;
TH1D *hTempPt  = NULL;
TH1D *hTempEta = NULL;
static TH1D ******hLineShape;
static TF1 ******fLineShape;


// function and the tree object
TF1   *fgengaus = 0;
TTree *dataTree = 0;
Int_t ntest = 1;


Double_t nEvents = 0;
Double_t nnorm   = 0;
Int_t ptBin;
Char_t fileName[255];      //file name of tree
Char_t fitFileName[255];   // file name for fit functions
Int_t fcentInput;
Int_t subsample;
Int_t analyseIter;    
Float_t fEtaDown = -0.8;
Float_t fEtaUp   = 0.8;
Float_t fpDown   = 0.2;
Float_t fpUp     = 1.5;
Float_t fptDown  = 0.;
Float_t fptUp    = 15;

Int_t fcentInputBin = 0;
Int_t fpDownBin     = 0;
Int_t fpUpBin       = 200;
Int_t fetaDownBin   = 0;
Int_t fetaUpBin     = 32;  

TString fnameused = "";
TString fitnameused = "";
Int_t systematics = -1;
Int_t isIntegrated = 1;
Int_t rapInterval = 0;
Int_t momInterval = 0;
Int_t isSim = 0;
Double_t fCountMean;
Double_t fCountSecond;
Double_t fCountMix;

// Tree for the final results
Double_t prpiNuDyn=0., prkaNuDyn=0., pikaNuDyn=0., pipiNuDyn=0., kakaNuDyn=0., prprNuDyn=0.;
Double_t el1=0.,  pi1=0.,  ka1=0.,  pr1=0.,  de1=0.,  mu1=0.;
Double_t el1Int=0.,  pi1Int=0.,  ka1Int=0.,  pr1Int=0;
Double_t el2=0.,  pi2=0.,  ka2=0.,  pr2=0.,  de2=0.,  mu2=0.;
Double_t elpi=0., elka=0., elpr=0., elde=0., elmu=0., pika=0., pipr=0., pide=0.;
Double_t pimu=0., kapr=0., kade=0., kamu=0., prde=0., prmu=0., demu=0.;

//TString fitFunctionGenGaus  = "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
Double_t fitFunctionGenGaus(Double_t *x, Double_t *par){
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
TF1 *CheckFitParams(Int_t centdebug)
{
    
    //
    //   Have a look at the fit Params if everything is ok.
    //
    
    TString fname  = "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
    TF1 *g = new TF1("g",fname,20.,100.); g->SetNpx(1000); g->SetLineColor(kMagenta); g->SetLineWidth(4);
    g->FixParameter(0,amp[0][centdebug][19][pion][1]);
    g->FixParameter(1,mean[0][centdebug][19][pion][1]);
    g->FixParameter(2,sigma[0][centdebug][19][pion][1]);
    g->FixParameter(3,kurtosis[0][centdebug][19][pion][1]);
    g->FixParameter(4,skew[0][centdebug][19][pion][1]);
    
    return g;
    
}
// -----------------------------------------------------------------------------------------
void readFitParams(TString paramTreeName, Int_t analyseIter)
{
    //
    // Read the fit parameters from the paramtree file "ParamTree.root" which comes from PID FITS
    //
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
    
    // Open the lookup table and initialize array
    fLineShapesLookUpTable = new TFile(paramTreeName);
    if (lookUpTableForLine){
        //
        cloneArrFunc   = (TClonesArray*)fLineShapesLookUpTable->Get("funcLineShapes"); 
        cloneArrHist   = (TClonesArray*)fLineShapesLookUpTable->Get("histLineShapes"); 
        if (!cloneArrFunc) cout << " Error:: cloneArrFunc is empty " << endl;
        if (!cloneArrHist) cout << " Error:: cloneArrFunc is empty " << endl;
        Int_t centBinRange[2] = {TMath::Max(fcentInputBin-1,0), TMath::Min(fcentInputBin+1,nCent)};
        Int_t etaBinRange[2]  = {TMath::Max(fetaDownBin-1,0), TMath::Min(fetaUpBin+1,nEta)};
        Int_t momBinRange[2]  = {TMath::Max(fpDownBin-1,0), TMath::Min(fpUpBin+1,nMom)};
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
        for (Int_t ipart = 0; ipart<nParticles; ipart++){
            for (Int_t icent = centBinRange[0]; icent< centBinRange[1]; icent++){
                for (Int_t ieta = etaBinRange[0]; ieta<etaBinRange[1]; ieta++){
                    for (Int_t imom = momBinRange[0]; imom<momBinRange[1]; imom++){
                        for (Int_t isign = 0; isign<nSign; isign++){
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
                            if (!htemp) htemp = (TH1D*)fLineShape[ipart][icent][ieta][(Int_t)TMath::Min(nMom,imom+1)][isign]->GetHistogram();
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
        
    } else{
        //
        cout << " Read lane shapes from ttree " << endl;
        treeLookUp = (TTree*)fLineShapesLookUpTable->Get(treeLineShapes);
        Int_t nent = treeLookUp -> GetEntries();
        cout << nent <<  "   tree is fine go ahead " << endl;       
        
        Int_t sign       = 0;
        Int_t it         = 0;
        Int_t sl         = 0;
        Int_t signBin = -100;
        
        
        Int_t myBin[3]  = {0};  // 0; eta, 1;cent, 2;p
        
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
            //if(sign != charge ) continue;
            if(sign == 0) continue;
            
            signBin = sign;
            if(signBin == -1) signBin = 0;
            
            if (fpp) myBin[1]=0;
            // Analyse only 4 iteration of the iterative fitting procedure
            if (it != analyseIter) continue;
            amp[myBin[0]][myBin[1]][myBin[2]][electron][signBin]      = elA;
            amp[myBin[0]][myBin[1]][myBin[2]][pion][signBin]          = piA;
            amp[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]          = kaA;
            amp[myBin[0]][myBin[1]][myBin[2]][proton][signBin]        = prA;
            
            mean[myBin[0]][myBin[1]][myBin[2]][electron][signBin]     = elM;
            mean[myBin[0]][myBin[1]][myBin[2]][pion][signBin]         = piM;
            mean[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]         = kaM;
            mean[myBin[0]][myBin[1]][myBin[2]][proton][signBin]       = prM;
            
            sigma[myBin[0]][myBin[1]][myBin[2]][electron][signBin]    = elSi;
            sigma[myBin[0]][myBin[1]][myBin[2]][pion][signBin]        = piSi;
            sigma[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]        = kaSi;
            sigma[myBin[0]][myBin[1]][myBin[2]][proton][signBin]      = prSi;
            
            skew[myBin[0]][myBin[1]][myBin[2]][electron][signBin]     = elSk;
            skew[myBin[0]][myBin[1]][myBin[2]][pion][signBin]         = piSk;
            skew[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]         = kaSk;
            skew[myBin[0]][myBin[1]][myBin[2]][proton][signBin]       = prSk;
            
            kurtosis[myBin[0]][myBin[1]][myBin[2]][electron][signBin] = elK;
            kurtosis[myBin[0]][myBin[1]][myBin[2]][pion][signBin]     = piK;
            kurtosis[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]     = kaK;
            kurtosis[myBin[0]][myBin[1]][myBin[2]][proton][signBin]   = prK;
            
        }
    }
    
}
// -----------------------------------------------------------------------------------------
void InitializeObjects()
{
    
    //
    // initialize lookup arrays
    hLineShape = new TH1D *****[nParticles];
    for (Int_t ipart = 0; ipart<nParticles; ipart++){
        hLineShape[ipart] = new TH1D****[nCent];
        for (Int_t icent = 0; icent<nCent; icent++){
            hLineShape[ipart][icent] = new TH1D***[nEta];
            for (Int_t ieta = 0; ieta<nEta; ieta++){
                hLineShape[ipart][icent][ieta] = new TH1D**[nMom];
                for (Int_t imom = 0; imom<nMom; imom++){
                    hLineShape[ipart][icent][ieta][imom] = new TH1D*[nSign];
                    for (Int_t isign = 0; isign<nSign; isign++){
                        hLineShape[ipart][icent][ieta][imom][isign] = NULL;
                    }
                }
            }
        }
    }
    //
    // 
    fLineShape = new TF1 *****[nParticles];
    for (Int_t ipart = 0; ipart<nParticles; ipart++){
        fLineShape[ipart] = new TF1****[nCent];
        for (Int_t icent = 0; icent<nCent; icent++){
            fLineShape[ipart][icent] = new TF1***[nEta];
            for (Int_t ieta = 0; ieta<nEta; ieta++){
                fLineShape[ipart][icent][ieta] = new TF1**[nMom];
                for (Int_t imom = 0; imom<nMom; imom++){
                    fLineShape[ipart][icent][ieta][imom] = new TF1*[nSign];
                    for (Int_t isign = 0; isign<nSign; isign++){
                        fLineShape[ipart][icent][ieta][imom][isign] = NULL;
                    }
                }
            }
        }
    }
    //
    // initialize output tree
    // 
    momTree = new TTree("momTree","momTree");
    momTree -> Branch("fileName",&fnameused);
    momTree -> Branch("fitFileName",&fitnameused);
    momTree -> Branch("fCountMean",&fCountMean);
    momTree -> Branch("fCountSecond",&fCountSecond);
    momTree -> Branch("fCountMix",&fCountMix);
    momTree -> Branch("sign",&charge);
    momTree -> Branch("el1",&el1);
    momTree -> Branch("pi1",&pi1);
    momTree -> Branch("ka1",&ka1);
    momTree -> Branch("pr1",&pr1);
    momTree -> Branch("el2",&el2);
    momTree -> Branch("pi2",&pi2);
    momTree -> Branch("ka2",&ka2);
    momTree -> Branch("pr2",&pr2);
    momTree -> Branch("elpi",&elpi);
    momTree -> Branch("elka",&elka);
    momTree -> Branch("elpr",&elpr);
    momTree -> Branch("pika",&pika);
    momTree -> Branch("pipr",&pipr);
    momTree -> Branch("kapr",&kapr);
    momTree -> Branch("nEvents",&nEvents);
    momTree -> Branch("el1Int",&el1Int);
    momTree -> Branch("pi1Int",&pi1Int);
    momTree -> Branch("ka1Int",&ka1Int);
    momTree -> Branch("pr1Int",&pr1Int);
    momTree -> Branch("nnorm",&nnorm);
    momTree -> Branch("fpDown",&fpDown);
    momTree -> Branch("fpUp",&fpUp);
    momTree -> Branch("fpDownBin",&fpDownBin);
    momTree -> Branch("fpUpBin",&fpUpBin);
    momTree -> Branch("fptDown",&fptDown);
    momTree -> Branch("fptUp",&fptUp);
    momTree -> Branch("rapInt",&rapInterval);
    momTree -> Branch("fEtaDown",&fEtaDown);
    momTree -> Branch("fEtaUp",&fEtaUp);
    
    fhEta   = new TH1D("fhEta" ,"Eta Bins"       ,fnEtaBins ,fEtaRangeDown, fEtaRangeUp );
    fhPtot  = new TH1D("fhPtot","Momentum Bins"  ,fnMomBins ,fMomRangeDown, fMomRangeUp ); 
    fhCent  = new TH1D("fhCent","Centrality Bins",fnCentBins ,xCentBins );
    hTempPt = new TH1D("hTempPt" ,"hTempPt" ,26,  0.2, 1.5); 
    hTempEta= new TH1D("hTempEta","hTempEta",16, -0.8, 0.8);
    
    if(rapInterval == 0) { fEtaDown = -0.8; fEtaUp   =  0.8; }
    if(rapInterval == 1) { fEtaDown = -0.7; fEtaUp   =  0.7; }
    if(rapInterval == 2) { fEtaDown = -0.6; fEtaUp   =  0.6; }
    if(rapInterval == 3) { fEtaDown = -0.5; fEtaUp   =  0.5; }
    if(rapInterval == 4) { fEtaDown = -0.4; fEtaUp   =  0.4; }
    if(rapInterval == 5) { fEtaDown = -0.3; fEtaUp   =  0.3; }
    if(rapInterval == 6) { fEtaDown = -0.2; fEtaUp   =  0.2; }
    if(rapInterval == 7) { fEtaDown = -0.1; fEtaUp   =  0.1; }
    if(rapInterval == 8) { fEtaDown = -0.05;fEtaUp   =  0.05;}
    
    const Int_t nMoments = 13;
    TVectorF moments(nMoments);
    TVectorF momentsPos(nMoments);
    TVectorF momentsNeg(nMoments);
    TVectorF momentsCross(nMoments);
    
    // initialize counters 
    for(Int_t i=0;i<nMoments; i++){  
        moments[i]=0.; 
        momentsPos[i]=0.;  
        momentsNeg[i]=0.; 
        momentsCross[i]=0.;
    }
    
}
// -----------------------------------------------------------------------------------------
Double_t EvalFitValue(Int_t particle, Double_t x)
{
    
    if (lookUpTableForLine){
        Int_t bin = hLineShape[particle][fCentBin][fEtaBin][fMomBin][fCharge]->FindBin(x);
        if (lookUpTableLineMode==0){ 
            return hLineShape[particle][fCentBin][fEtaBin][fMomBin][fCharge]->GetBinContent(bin);
        }
        if (lookUpTableLineMode==1){ 
            return fLineShape[particle][fCentBin][fEtaBin][fMomBin][fCharge]->Eval(x);
        }
        
    } else {
        fgengaus -> SetParameter(0,amp[fEtaBin][fCentBin][fMomBin][particle][fCharge]);
        fgengaus -> SetParameter(1,mean[fEtaBin][fCentBin][fMomBin][particle][fCharge]);
        fgengaus -> SetParameter(2,sigma[fEtaBin][fCentBin][fMomBin][particle][fCharge]);
        fgengaus -> SetParameter(3,kurtosis[fEtaBin][fCentBin][fMomBin][particle][fCharge]);
        fgengaus -> SetParameter(4,skew[fEtaBin][fCentBin][fMomBin][particle][fCharge]);
        return fgengaus->Eval(x);
    }    
    
}
// -----------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    
    //     
    // 4 arguments should be given ./testIden <dEdxListName> <fitFileName> cent subsample
    //     
    
    cout<<"main.Info: NUMBER OF ARGUMENTS "<<argc<<endl;
    if(argc == 12)
    {
        sprintf(fileName,"%s",argv[1]);
        sprintf(fitFileName,"%s",argv[2]);
        fcentInput    = atoi(argv[3]);
        subsample     = atoi(argv[4]);
        analyseIter   = atoi(argv[5]);
        charge        = atoi(argv[6]);
        systematics   = atoi(argv[7]);
        isIntegrated  = atoi(argv[8]);
        rapInterval   = atoi(argv[9]);
        momInterval   = atoi(argv[10]);
        isSim         = atoi(argv[11]);
        cout<<"main.Info: read file names from input"<<endl;
    }
    else
    {
        cout<<"main.Error: wrong input list"<<endl;
    }
    
    fnameused   = fileName;
    fitnameused = fitFileName;
    InitializeObjects();
    
    Bool_t useMomCuts  = kTRUE; //kFalse if pt cut is used
    Float_t minMom = 0.6;
    Float_t maxMom = 1.5;
    Float_t minPt  = 0.;
    Float_t maxPt  = 15.;
    
    /////////// MOM  RANGE
    if(momInterval == 0 ) { minMom = 0.6; maxMom = 1.5; }
    if(momInterval == 1 ) { minMom = 0.6; maxMom = 1.2; }
    if(momInterval == 2 ) { minMom = 0.6; maxMom = 1.;  }
    if(momInterval == 3 ) { minMom = 0.6; maxMom = 0.8; }
    if(momInterval == 4 ) { minMom = 0.2; maxMom = 1.5; }
    /////////
    
    if(useMomCuts)
    {
        fpDown = minMom;
        fpUp   = maxMom;
    }
    
    fpDownBin     = fhPtot -> FindBin(fpDown + 0.0000001) - 1;
    fpUpBin       = fhPtot -> FindBin(fpUp - 0.0000001)   - 1;
    fetaDownBin   = fhEta  -> FindBin(fEtaDown + 0.0000001)  - 1;
    fetaUpBin     = fhEta  -> FindBin(fEtaUp   - 0.0000001)  - 1;  
    fcentInputBin = fhCent -> FindBin(fcentInput + 0.0000001) -1;
    
    Int_t etaRange[2] = {fetaDownBin,fetaUpBin};
    Int_t momRange[2] = {fpDownBin,fpUpBin};     
    
    TROOT IdentityMethod("IdentityMethod","compiled identity method");
    
    // full path to data tree list
    TString tmpFileName(fileName);
    tmpFileName.Prepend(inDir);
    
    // Full path to param tree file
    TString tmpfitFileName(fitFileName);
    TString paramTreeName(inDir);
    paramTreeName.Append(tmpfitFileName);
    
    cout << " ================================================================================= " << endl;
    cout << " main.Info: Inputs: " << endl;
    cout << " main.Info: dEdx List             = " << fileName       << endl;
    cout << " main.Info: Fit File              = " << fitFileName    << endl;
    cout << " main.Info: Centrality bin        = " << fcentInputBin      << endl;
    cout << " main.Info: Subsample index       = " << subsample      << endl;
    cout << " main.Info: Fit iteration         = " << analyseIter    << endl;
    cout << " main.Info: Systematics           = " << systematics    << endl;
    cout << " main.Info: Eta Range             = " << fEtaDown       << " - " <<  fEtaUp << endl;
    cout << " main.Info: Momentum Range        = " << fpDown         << " - " <<  fpUp << endl;
    cout << " main.Info: Centrality            = " << fcentInput      << endl;
    cout << " main.Info: Charge                = " << charge         << endl;
    cout << " ================================================================================= " << endl;
    cout << " main.Info: Centrality (Bin)      = " << fcentInputBin      << endl;
    cout << " main.Info: Eta Range (Bin)       = " << etaRange[0]    << " - " <<  etaRange[1] << endl;
    cout << " main.Info: Momentum Range (Bin)  = " << momRange[0]    << " - " <<  momRange[1] << endl;
    cout << " main.Info: subsample index       = " << subsample      << endl;
    cout << " main.Info: data Tree location  --> " << tmpFileName    << endl;
    cout << " main.Info: fit params location --> " << paramTreeName  << endl;
    cout << " ================================================================================= " << endl;
    
    readFitParams(paramTreeName,analyseIter);
    fgengaus = new TF1("fgengaus",fitFunctionGenGaus,0,1000,5);
    fgengaus->SetNpx(10000);
    
    // Create the TIdentity2D object and start analysis
    TIdentity2D *iden4 = new TIdentity2D(4);      // Set the number of particles to 4
    iden4 -> SetInputDir(inDir);
    iden4 -> SetOutputDir(outDir);
    iden4 -> SetFileName(fileName);
    iden4 -> SetFunctionPointers(EvalFitValue);
    iden4 -> SetLimits(0.,1020.,10.); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
    iden4 -> SetUseSign(charge);
    Long_t nEntries;
    if (isSim){
        cout<<"main.Info: Running for sim"<<endl;
        iden4 -> GetTree(nEntries,treeIdentityMC);
    } else{
        cout<<"main.Info: Running for data"<<endl;
        iden4 -> GetTree(nEntries,treeIdentity);
    }
    
    // Get the centrality bin information using fhCent
    Int_t centBinLimit  = (fpp) ? 0 : fhCent->FindBin(fcentInput+0.025)-1;
    Double_t centrality = fhCent->GetBinCenter(centBinLimit);
    
    // Calculate First moments
    // -------------- DEBUGGING -----------------
    TF1 *fDebug;
    fDebug = (fpp) ? CheckFitParams(centBinLimit+1) : CheckFitParams(centBinLimit) ;
    // ------------------------------------------
    //     
    // Write results into tree
    TFile *outFile = new TFile(Form("TImoments_sim%d_sub%d_sys%d_cent%d_int%d_rap%d_mom%d_sign%d.root",isSim,subsample,systematics,fcentInput, isIntegrated,rapInterval, momInterval,charge),"recreate");
    
    if (test) nEntries = 2000000;
    Int_t usedBins[fnEtaBins][fnCentBins][fnMomBins][2];  // usedBins[fEtaBin][fCentBin][fMomBin][fCharge] = 1;
    
    iden4 -> SetMinPt(0.);
    iden4 -> SetMaxPt(1000.);
    iden4 -> SetMinMom(fpDown);
    iden4 -> SetMaxMom(fpUp);
    
    Int_t bins[3];
    for( Int_t i = 0; i < fnEtaBins; i++ )           // loop over eta bins
        for( Int_t j = 0; j < fnCentBins; j++ )      // loop over cent bins 
            for( Int_t k = 0; k < fnMomBins; k++ )   // loop over momentum bins
                for(Int_t l = 0; l < 2; l++)         // loop over charge combinations
                {
                    usedBins[i][j][k][l] = -1;
                }
                iden4 -> Reset();
            
            // track by track loop --> read all track info 
            for( Int_t i = 0; i < nEntries; i++ )
            {
                if( !iden4 ->  GetEntry(i) ) continue;
                iden4      ->  GetBins( bins );    // reads identity tree and retrives mybin[] info
                fEtaBin  = bins[0];
                fCentBin = bins[1];
                fMomBin  = bins[2];
                fCharge = iden4 -> GetSign();
                //cout<<"testting sign "<<fCharge<<endl;
                if( fCharge == -1) fCharge = 0;
                if( fEtaBin < 0 || fCentBin < 0 || fMomBin < 0 ) 
                {
                    cout<<" main.Error: this should not happen "<< fEtaBin <<"  "<<fCentBin<<"  "<<fMomBin << endl; 
                    continue;
                }
                if( fMomBin  > fnMomBins)  continue;
                if( fCentBin != centBinLimit ) continue;
                if( fMomBin < momRange[0] || fMomBin > momRange[1] ) continue; //no upper momentum cut;;;
                if( fEtaBin < etaRange[0] || fEtaBin > etaRange[1] ) continue;
                usedBins[fEtaBin][fCentBin][fMomBin][fCharge] = 1;
                //                 
                Bool_t isAdd = kFALSE;
                iden4 -> AddEntry(isAdd);
            }
            
            iden4 -> Finalize();
            // First Moments
            el1 = iden4 -> GetMean(electron);
            pi1 = iden4 -> GetMean(pion);
            ka1 = iden4 -> GetMean(kaon);
            pr1 = iden4 -> GetMean(proton);
            
            fCountMean   = iden4 -> GetAverCount(0);
            fCountSecond = iden4 -> GetAverCount(1);
            fCountMix    = iden4 -> GetAverCount(2);
            
            el2 = pi2 = ka2 = pr2 = elpi = elka = elpr = pika = pipr = kapr = -10000;
            el1Int = pi1Int = ka1Int = pr1Int = -1000;
            nnorm = 1;
            
            // Calculate 2. order moments only for full range
            timer.Reset(); timer.Start();
            
            cout << "==================================" << endl;
            cout << " main.Info: calculating integrals " << endl;
            timer.Reset(); timer.Start();
            cout << "==================================" << endl;
            for(Int_t i = 0; i < fnEtaBins; i++)
                for(Int_t j = 0; j < fnCentBins; j++)
                    for(Int_t k = 0; k < fnMomBins; k++)
                        for(Int_t l = 0; l < 2; l++)
                        {
                            if(usedBins[i][j][k][l] != 1) continue;
                            fEtaBin   =  i;
                            fCentBin  =  j;
                            fMomBin   =  k;
                            fCharge = l;
                            iden4  -> AddIntegrals(charge);
                        }
                        iden4 -> CalcMoments();
                    
                    // Second Moments
                    el2 = iden4 -> GetSecondMoment(electron);
                    pi2 = iden4 -> GetSecondMoment(pion);
                    ka2 = iden4 -> GetSecondMoment(kaon);
                    pr2 = iden4 -> GetSecondMoment(proton);
                    
                    // Mixed Moments
                    elpi = iden4 -> GetMixedMoment(electron,pion);
                    elka = iden4 -> GetMixedMoment(electron,kaon);
                    elpr = iden4 -> GetMixedMoment(electron,proton);
                    pika = iden4 -> GetMixedMoment(pion,kaon);
                    pipr = iden4 -> GetMixedMoment(pion,proton);
                    kapr = iden4 -> GetMixedMoment(kaon,proton);
                    
                    //Integrals:
                    el1Int = iden4 -> GetMeanI(electron);
                    pi1Int = iden4 -> GetMeanI(pion);
                    ka1Int = iden4 -> GetMeanI(kaon);
                    pr1Int = iden4 -> GetMeanI(proton);
                    nnorm     = pi1/pi1Int;
                    
                    nEvents   = iden4 -> GetNEvents();
                    cout<<" =============================== first moments =============================== "<<endl;
                    cout<<"events      : "<< nEvents << endl;
                    cout<<"electron    : "<< el1   <<" int: "<< el1Int*nnorm<<endl;
                    cout<<"pion        : "<< pi1   <<" int: "<< pi1Int*nnorm<<endl;
                    cout<<"kaon        : "<< ka1   <<" int: "<< ka1Int*nnorm<<endl;
                    cout<<"proton      : "<< pr1   <<" int: "<< pr1Int*nnorm<<endl;
                    cout<<"fCountMean   : "<<fCountMean<<endl;
                    cout<<"fCountSecond : "<<fCountSecond<<endl;
                    cout<<"electron2   : "<< el2  <<endl;
                    cout<<"pion2       : "<< pi2  <<endl;
                    cout<<"kaon2       : "<< ka2  <<endl;
                    cout<<"proton2     : "<< pr2  <<endl;
                    momTree -> Fill();
                    cout << "====================================" << endl;
                    cout << " main.Info: calculation is finished " << endl;
                    timer.Stop(); timer.Print();
                    cout << "====================================" << endl;
                    
                    outFile -> cd();
                    momTree -> Write();
                    outFile -> Close();
                    
                    delete outFile; //yeni eklave etdim.
                    delete fDebug;
                    delete iden4;
                    //delete moments;
                    return 1;
}


