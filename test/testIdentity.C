#include "iostream"
#include "TIdentity2D.h"
// #include "TIdentityBase.h"
//#include "TTreeStream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1F.h"
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

// ======= Modification Part =============================================================================
Int_t useSign = -100;
Bool_t fpp             = kFALSE;
const Int_t nEtaBins   = 32;///16;
const Int_t nMomBins   = 150;
const Int_t nParticles = 4;
TString treeName       = "fEventTree";
Bool_t test = kFALSE;
const Int_t nCentBins  = 9;
Float_t xCentBins[] = {0, 5,  10,  20, 30, 40, 50, 60, 70, 80};
// =======================================================================================================
Double_t amp[nEtaBins][nCentBins][nMomBins][nParticles][2]; // amp[etaBin][centBin][momBin][particleType], before kurtosis[9][12][1000][6];
Double_t mean[nEtaBins][nCentBins][nMomBins][nParticles][2];
Double_t sigma[nEtaBins][nCentBins][nMomBins][nParticles][2];
Double_t skew[nEtaBins][nCentBins][nMomBins][nParticles][2];
Double_t kurtosis[nEtaBins][nCentBins][nMomBins][nParticles][2];
Int_t etaBin, centBin, momBin;
Int_t idenSign;
// Values to set from outside
TString inDir("/hera/alice/marsland/pFluct/files/analysis/Data/");
TString outDir("./");

// function and the tree object
TF1   *fgengaus = 0;
TTree *dataTree = 0;
Int_t ntest = 1;
Int_t centBinLimit;            // 0-5-10-20-30-40-50-60-70-80 is given as input and converted to bin info (is given from outside)

Float_t getRapidity(Float_t px, Float_t py, Float_t pz, TString pid)
{
    Float_t mp  =  0.938272;
    Float_t mpi =  0.13957;
    Float_t mk  =  0.493677;
    TLorentzVector v;
    Float_t E;
    Float_t mom2 = px*px + py*py + pz*pz;
    if (pid == "proton")
        E = sqrt(mom2 + mp*mp);
    else
        if(pid == "pion")
            E = sqrt(mom2 + mpi*mpi);
        else
            if (pid == "kaon")
                E = sqrt(mom2 + mk*mk);
            else
            { cout << "enter vslid pid for rapidity calculation"<<endl; }
    v.SetPxPyPzE(px, py, pz, E);
    //return v.Rapidity();
    return v.PseudoRapidity();
}

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
TF1 *CheckFitParams(Int_t centdebug){
    
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
    
    TFile f(paramTreeName);
    TTree* tree = (TTree*)f.Get("treeId");
    Int_t nent = tree -> GetEntries();
    
    Int_t sign       = 0;
    Int_t it         = 0;
    Int_t sl         = 0;
    
    Int_t signBin = -100;
    
    
    Int_t myBin[3]  = {0};  // 0; eta, 1;cent, 2;p
    
    Double_t elA  = 0., elM  = 0., elSi = 0., elK = 0., elSk = 0.;
    Double_t piA  = 0., piM  = 0., piSi = 0., piK = 0., piSk = 0.;
    Double_t kaA  = 0., kaM  = 0., kaSi = 0., kaK = 0., kaSk = 0.;
    Double_t prA  = 0., prM  = 0., prSi = 0., prK = 0., prSk = 0.;
    
    tree->SetBranchAddress("sign"   ,&sign);
    tree->SetBranchAddress("it"     ,&it);
    tree->SetBranchAddress("sl"     ,&sl);
    tree->SetBranchAddress("myBin"  ,myBin);
    
    tree->SetBranchAddress("elM" ,&elM);
    tree->SetBranchAddress("piM" ,&piM);
    tree->SetBranchAddress("kaM" ,&kaM);
    tree->SetBranchAddress("prM" ,&prM);
    
    tree->SetBranchAddress("elSi" ,&elSi);
    tree->SetBranchAddress("piSi" ,&piSi);
    tree->SetBranchAddress("kaSi" ,&kaSi);
    tree->SetBranchAddress("prSi" ,&prSi);
    
    tree->SetBranchAddress("elA" ,&elA);
    tree->SetBranchAddress("piA" ,&piA);
    tree->SetBranchAddress("kaA" ,&kaA);
    tree->SetBranchAddress("prA" ,&prA);
    
    tree->SetBranchAddress("elSk" ,&elSk);
    tree->SetBranchAddress("piSk" ,&piSk);
    tree->SetBranchAddress("kaSk" ,&kaSk);
    tree->SetBranchAddress("prSk" ,&prSk);
    
    tree->SetBranchAddress("elK" ,&elK);
    tree->SetBranchAddress("piK" ,&piK);
    tree->SetBranchAddress("kaK" ,&kaK);
    tree->SetBranchAddress("prK" ,&prK);
    
    for(Int_t i = 0; i < nent; ++i)
    {
        // myBin[0] --> Eta, myBin[1]--> Centrality, myBin[2]-->Momentum
        tree -> GetEntry(i);
        //if(sign != useSign ) continue;
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
// -----------------------------------------------------------------------------------------
Double_t EvalFitValue(Int_t i, Double_t x)
{
    
    // i stands for particle type
    fgengaus -> SetParameter(0,amp[etaBin][centBin][momBin][i][idenSign]);
    fgengaus -> SetParameter(1,mean[etaBin][centBin][momBin][i][idenSign]);
    fgengaus -> SetParameter(2,sigma[etaBin][centBin][momBin][i][idenSign]);
    fgengaus -> SetParameter(3,kurtosis[etaBin][centBin][momBin][i][idenSign]);
    fgengaus -> SetParameter(4,skew[etaBin][centBin][momBin][i][idenSign]);
    return fgengaus->Eval(x);
    
}
// -----------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    
    Char_t fileName[255];      //file name of tree
    Char_t fitFileName[255];   // file name for fit functions
    Int_t centInput;
    Int_t subsample;
    Int_t analyseIter;
    Float_t etaDown = -0.8;
    Float_t etaUp   = 0.8;
    Float_t pDown   = 0.2;
    Float_t pUp     = 1.5;
    Float_t ptDown  = 0.;
    Float_t ptUp    = 15;
    
    Int_t ptBin;
    TString selectRapidity = "";
    
    // 4 arguments should be given ./testIden <dEdxListName> <fitFileName> cent subsample
    
    cout<<"NUMBER OF ARGUMENTS "<<argc<<endl;
    
    Int_t ptEffBin  = 0;
    Int_t etaEffBin = 0;
    
    Int_t systematics = -1;
    
    Int_t isIntegrated = 1;
    Int_t rapInterval = 0;
    
    Int_t momInterval = 0;
    
    Int_t isSim = 0;
    
    if(argc == 12)
    {
        sprintf(fileName,"%s",argv[1]);
        sprintf(fitFileName,"%s",argv[2]);
        centInput  = atoi(argv[3]);
        subsample  = atoi(argv[4]);
        analyseIter= atoi(argv[5]);
        //ptEffBin   = atoi(argv[6]) + 1; // provide starting from 0
        //etaEffBin  = atoi(argv[7]) + 1; // provide starting from 0
        //isIntegrated = atoi(argv[6]);
        useSign    = atoi(argv[6]);
        systematics = atoi(argv[7]);
        isIntegrated = atoi(argv[8]);
        rapInterval=atoi(argv[9]);
        momInterval=atoi(argv[10]);
        isSim      = atoi(argv[11]);
        cout<<"read file names from input"<<endl;
    }
    else
    {
        sprintf(fitFileName,"IdentityInput_FitParamTree.root");
        sprintf(fileName,"IdentityInput_DataTree.list");
        centInput  = 0;
        subsample  = 100;
        analyseIter= 6;
        etaDown    = -0.8;
        etaUp      = 0.8;
        pDown      = 0.2;
        pUp        = 1.5;
        cout<<"read default file names"<<endl;
    }
    
    TString simName = "";
    if(isSim)
    {
        simName = "_sim_";
    }
    
    
    if(rapInterval == 0)
    {
        etaDown = -0.8;
        etaUp   =  0.8;
    }

    
    if(rapInterval == 1)
    {
        etaDown = -0.7;
        etaUp   =  0.7;
    }
    
    if(rapInterval == 2)
    {
        etaDown = -0.6;
        etaUp   =  0.6;
    }

    if(rapInterval == 3)
    {
        etaDown = -0.5;
        etaUp   =  0.5;
    }
    
    if(rapInterval == 4)
    {
        etaDown = -0.4;
        etaUp   =  0.4;
    }
    
    if(rapInterval == 5)
    {
        etaDown = -0.3;
        etaUp   =  0.3;
    }
    
    if(rapInterval == 6)
    {
        etaDown = -0.2;
        etaUp   =  0.2;
    }

    if(rapInterval == 7)
    {
        etaDown = -0.1;
        etaUp   =  0.1;
    }
    
    if(rapInterval == 8)
    {
        etaDown = -0.05;
        etaUp   =  0.05;
    }
    
    
    
    Bool_t useMomCuts  = kTRUE; //kFalse if pt cut is used
    
    
    
    Float_t minMom = 0.6;
    Float_t maxMom = 1.5;
    Float_t minPt  = 0.;
    Float_t maxPt  = 15.;
    
    
    
    /////////// MOM  RANGE
    
    if(momInterval == 0 )
    {
        minMom = 0.6;
        maxMom = 1.5;
    }
    
    if(momInterval == 1 )
    {
        minMom = 0.6;
        maxMom = 1.2;
    }
    
    if(momInterval == 2 )
    {
        minMom = 0.6;
        maxMom = 1.;
    }
    
    if(momInterval == 3 )
    {
        minMom = 0.6;
        maxMom = 0.8;
    }
    
    if(momInterval == 4 )
    {
        minMom = 0.2;
        maxMom = 1.5;
    }

    
    
    /////////
    
    

    
    if(!isIntegrated)
    {
        minMom = 0.2;
        maxMom = 3.;
        useMomCuts = kTRUE;
    }
    

    if(useMomCuts)
    {
       pDown = minMom;
       pUp   = maxMom;
    }
    else
    {
        ptDown = minPt;
        ptUp   = maxPt;
    }
    
    
    
    
    TH1F *fhPtot =  new TH1F("fhPtot","Momentum Bins"  ,150   ,0.2 , 3.2 );
    TH1F *fhEta;
    
    if(isSim)
    {
    fhEta = new TH1F("fhEta","Eta bins",8,-0.8,0.8);
    }
    else
    {
    //fhEta = new TH1F("fhEta","Eta bins",32,-0.8,0.8);
    fhEta = new TH1F("fhEta","Eta bins",16,-0.8,0.8);
    }
    
    
    
    Int_t pDownBin =  fhPtot -> FindBin(pDown)    - 1;
    Int_t pUpBin   =  fhPtot -> FindBin(pUp)   - 1;
    
    Int_t momRange[2] = {pDownBin  ,pUpBin};     // 0-367 is whole momentum range each bein is 10 MeV up to 300
    
    Int_t etaDownBin = fhEta -> FindBin(etaDown + 0.0000001)  - 1;
    Int_t etaUpBin   = fhEta -> FindBin(etaUp   - 0.0000001)  - 1;
    
    
    
    Int_t etaRange[2] = {etaDownBin,etaUpBin};
    
    
    cout << " Inputs: " << endl;
    cout << fileName    << endl;
    cout << fitFileName << endl;
    
    cout << "cent        = " <<   centInput   << endl;
    cout << "subsample   = " <<   subsample   << endl;
    cout << "iter        = " <<   analyseIter << endl;
    cout << "etaBin      = " <<   etaEffBin   << endl;
    
    TROOT IdentityMethod("IdentityMethod","compiled identity method");
    
    // prepare histogram to find the centrality bin
    TH1D * hCent = new TH1D( "hCent","Centrality Bins",nCentBins,xCentBins );
    
    // full path to data tree list
    TString tmpFileName(fileName);
    tmpFileName.Prepend(inDir);
    
    // Full path to param tree file
    TString tmpfitFileName(fitFileName);
    TString paramTreeName(inDir);
    paramTreeName.Append(tmpfitFileName);
    
    cout << " ================================================================================= " << endl;
    cout << " Fit iteration         = " << analyseIter    << endl;
    //cout << " Eta Range             = " << etaRange[0]    << " - " <<  etaRange[1] << endl;
    cout << " Momentum Range        = " << momRange[0]    << " - " <<  momRange[1] << endl;
    cout << " dEdx List             = " << fileName       << endl;
    cout << " Fit File              = " << fitFileName    << endl;
    cout << " centrality bin        = " << centInput      << endl;
    cout << " subsample index       = " << subsample       << endl;
    cout << " data Tree location  --> " << tmpFileName    << endl;
    cout << " fit params location --> " << paramTreeName  << endl;
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
    iden4 -> SetUseSign(useSign);
    Long_t nEntries;
    if (isSim)
    {
    iden4 -> GetTree(nEntries,"sim");
    }
    else
    {
        iden4 -> GetTree(nEntries,"DATA");
    }
    
    // Tree for the final results
    Double_t prpiNuDyn=0., prkaNuDyn=0., pikaNuDyn=0., pipiNuDyn=0., kakaNuDyn=0., prprNuDyn=0.;
    Double_t el1=0.,  pi1=0.,  ka1=0.,  pr1=0.,  de1=0.,  mu1=0.;
    Double_t el1Int=0.,  pi1Int=0.,  ka1Int=0.,  pr1Int=0;
    Double_t el2=0.,  pi2=0.,  ka2=0.,  pr2=0.,  de2=0.,  mu2=0.;
    Double_t elpi=0., elka=0., elpr=0., elde=0., elmu=0., pika=0., pipr=0., pide=0.;
    Double_t pimu=0., kapr=0., kade=0., kamu=0., prde=0., prmu=0., demu=0.;
    
    // Get the centrality bin information using hCent
    Int_t centBinLimit  = (fpp) ? 0 : hCent->FindBin(centInput+0.25)-1;
    Double_t centrality = hCent->GetBinCenter(centBinLimit);
    
    // File for the moments output which are dumped into a tree
    
    
    
    // Calculate First moments
    Int_t bins[3];
    
    // -------------- DEBUGGING -----------------
    TF1 *fDebug;
    fDebug = (fpp) ? CheckFitParams(centBinLimit+1) : CheckFitParams(centBinLimit) ;
    // ------------------------------------------
    
    // Write results into tree
    
    
    
    
    
    
    TFile *outFile;
    if(useSign == 0)
        outFile = new TFile(Form("TImoments_cent%d_sub%d_sys%d_int%d_rap%d_mom%d%s.root",centInput, subsample,   systematics, isIntegrated,rapInterval, momInterval, simName.Data()),"recreate");
    else
        if(useSign == 1)
            outFile = new TFile(Form("TImoments_cent%d_sub%d_sys%d_int%d_rap%d_mom%d_pos%s.root",centInput, subsample, systematics, isIntegrated,rapInterval, momInterval, simName.Data()),"recreate");
        else
            outFile = new TFile(Form("TImoments_cent%d_sub%d_sys%d_int%d_rap%d_mom%d_neg%s.root",centInput, subsample, systematics, isIntegrated,rapInterval, momInterval, simName.Data()),"recreate");
    
    
    Double_t nEvents = 0;
    Double_t nnorm   = 0;
    
    
    
    
    // Int_t ptEffBin = 0;
    // Int_t etaEffBin = 0;
    Int_t phiEffBin = 0;
    
    std::string fnameused =   fileName;
    std::string fitnameused = fitFileName;
    
    Double_t countMean;
    Double_t countSecond;
    Double_t countMix;
    
    TTree *momTree = new TTree("momTree","momTree");
    momTree -> Branch("fileName",&fnameused);
    momTree -> Branch("fitFileName",&fitnameused);
    
    momTree -> Branch("countMean",&countMean);
    momTree -> Branch("countSecond",&countSecond);
    momTree -> Branch("countMix",&countMix);
    
    
    momTree -> Branch("sign",&useSign);
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
    momTree -> Branch("ptbin",&ptEffBin);
    momTree -> Branch("etabin",&etaEffBin);
    momTree -> Branch("el1Int",&el1Int);
    momTree -> Branch("pi1Int",&pi1Int);
    momTree -> Branch("ka1Int",&ka1Int);
    momTree -> Branch("pr1Int",&pr1Int);
    momTree -> Branch("nnorm",&nnorm);
    momTree -> Branch("pDown",&pDown);
    momTree -> Branch("pUp",&pUp);
    momTree -> Branch("pDownBin",&pDownBin);
    momTree -> Branch("pUpBin",&pUpBin);
    momTree -> Branch("ptDown",&ptDown);
    momTree -> Branch("ptUp",&ptUp);
    momTree -> Branch("rapInt",&rapInterval);
    momTree -> Branch("etaDown",&etaDown);
    momTree -> Branch("etaUp",&etaUp);
 

    
    
    
    //momTree -> Branch("phibin",&phiEffBin);
    
    TNtuple *nt = new TNtuple("idenNT","idenNt","yp:yk:ypi:pt:phi:wp:wk:wpi:momx:momy:momz");
    
    TH1F *hTempPt = new TH1F("hTempPt","hTempPt",26, 0.2, 1.5); //5 MeV bin
    TH1F*  hTempEta = new TH1F("hTempEta","",16, -0.8, 0.8);
    
    Int_t iptMax  = hTempPt  -> GetXaxis() -> GetNbins() + 1;
    Int_t ietaMax = hTempEta -> GetXaxis() -> GetNbins() + 1;
    
    //nEntries = 20000;
    
    Int_t usedBins[33][11][160][2];
    
   
    if(isIntegrated)
    {
     iptMax  = 1;  //fast check for integrated
     ietaMax = 1;  //fast check for integrated
    }
    
    Int_t iptMin = hTempPt -> FindBin(ptDown); 
    
    if(iptMax == 1)
    {
        iptMin = 1;
    }
    
    iden4 -> SetMinPt(ptDown);
    iden4 -> SetMaxPt(ptUp);
    iden4 -> SetMinMom(pDown);
    iden4 -> SetMaxMom(pUp);
    
    
    
    for (Int_t ipt = iptMin; ipt < iptMax + 1; ipt++) //iptMax is for full range
    {
        ptEffBin  = ipt;
        for(Int_t ieta = 1; ieta < ietaMax + 1; ieta++) //ietaMax is for full range
        {
            etaEffBin = ieta;
            
            if( ( ipt == iptMax && ieta != ietaMax ) || (ipt != iptMax && ieta == ietaMax) ) continue;
            
            
            if ( ipt == iptMax && ieta == ietaMax )
            {
                iptMax  = 1;
                ietaMax = 1;
                
                cout<<"reached max values "<<iptMax<<"  "<<ietaMax<<endl;
                
            }
            
            for( Int_t i = 0; i < 33; i++ )
                for( Int_t j = 0; j < 11; j++ )
                    for( Int_t k = 0; k < 160; k++ )
                        for(Int_t l = 0; l < 2; l++)
                    {
                        usedBins[i][j][k][l] = -1;
                    }
            
            iden4 -> Reset();
            
            for( Int_t i = 0; i < nEntries; i++ )
            {
                if( !iden4 ->  GetEntry(i) ) continue;
                iden4      ->  GetBins( bins );
                ////
                
                etaBin  = bins[0];
                centBin = bins[1];
                momBin  = bins[2];
                
                idenSign = iden4 -> GetSign();
                
                //cout<<"testting sign "<<idenSign<<endl;
                
                if(idenSign == -1) idenSign = 0;
                
                //cout<<"idenSign == "<<idenSign<<endl;
                
                if( etaBin < 0 || centBin < 0 || momBin < 0 ) {cout<<" this should not happen "<< etaBin <<"  "<<centBin<<"  "<<momBin << endl; continue;}
                
//cout<<"bins "<< bins[0]<<"  "<<bins[1]<<"  "<<bins[2]<<endl;
                
                
                //cout<<"etaBinTest "<<etaBin<<endl;
                
               // if( etaBin  > 15  ) continue;
                
               // if(isSim && etaBin > 7) continue;
                
                if( momBin  > 148)  continue;
                
                if( ietaMax != 1 && etaBin != ( etaEffBin -1 ) ) continue;
                
                
                ///// TESTING, WILL BE WRERITTEN
                Float_t mpx      = iden4  ->  GetMomX();
                Float_t mpy      = iden4  ->  GetMomY();
                Float_t mpz      = iden4  ->  GetMomZ();
                Int_t   mcharge  = iden4  ->  GetSign();
                Float_t mpt = sqrt(mpx*mpx + mpy*mpy);
                Float_t mphi = TMath::ATan2(mpy, mpx);
                if( mphi < 0 ) mphi += TMath::TwoPi();
                
                Float_t rapProton = getRapidity(mpx, mpy, mpz, "proton");
                Float_t rapPion   = getRapidity(mpx, mpy, mpz, "pion");
                Float_t rapKaon   = getRapidity(mpx, mpy, mpz, "kaon");
                
                iden4 -> SetEffPion(1.);
                iden4 -> SetEffKaon(1.);
                iden4 -> SetEffProton(1.);
                
                
                
                if( iptMax != 1 && hTempPt -> FindBin(mpt) != ptEffBin ) continue;
                
                
                if( centBin != centBinLimit ) continue;
                
                if( mpt > ptUp || mpt < ptDown) {cout<<"this pt cut should not work  "<<endl; continue;} // pt cut;;
                
                
                if( momBin < momRange[0] ||  momBin > momRange[1] ) continue; //no upper momentum cut;;;
                
                
                if( etaBin < etaRange[0] || etaBin > etaRange[1] ) continue;
                
                
                usedBins[etaBin][centBin][momBin][idenSign] = 1;
                
                iden4 -> AddEntry();
                
                Float_t wwp  = iden4 -> GetwProton();
                Float_t wwk  = iden4 -> GetwKaon();
                Float_t wwpi = iden4 -> GetwPion();
                
                if(i < 200000 && subsample == 0)
                    nt -> Fill(rapProton, rapKaon, rapPion, mpt, mphi, wwp, wwk, wwpi, mpx, mpy, mpz);
            }
            
            iden4 -> Finalize();
            // First Moments
            el1 = iden4 -> GetMean(electron);
            pi1 = iden4 -> GetMean(pion);
            ka1 = iden4 -> GetMean(kaon);
            pr1 = iden4 -> GetMean(proton);
            
            countMean   = iden4 -> GetAverCount(0);
            countSecond = iden4 -> GetAverCount(1);
            countMix    = iden4 -> GetAverCount(2);
            
            el2 = pi2 = ka2 = pr2 = elpi = elka = elpr = pika = pipr = kapr = -10000;
            el1Int = pi1Int = ka1Int = pr1Int = -1000;
            nnorm = 1;
            
            
            cout<<"pt bin: "<<ptEffBin <<"    eta bin: "<<etaEffBin <<"    phi bin: "<<phiEffBin<< endl;
            
            // Calculate 2. order moments only for full range
            if(  iptMax == 1 && ietaMax == 1 )
                
            {
                
                ptEffBin = etaEffBin = 0;
                
                cout<<"new binnings "<<endl;
                cout<<"pt bin: "<<ptEffBin <<"    eta bin: "<<etaEffBin <<"    phi bin: "<<phiEffBin<< endl;
                
                
                
                
                cout<<" calculating integrals " <<endl;
                //Int_t mmm = 0;
                
                for(Int_t i = 0; i < 33; i++)
                    for(Int_t j = 0; j < 11; j++)
                        for(Int_t k = 0; k < 160; k++)
                            for(Int_t l = 0; l < 2; l++)
                        {
                            
                           // mmm++;
                            //cout << " etaBin = " << i           << setw(15) <<  " centBin = "        << k << setw(15) ;
                            //cout << " iter   = " << analyseIter << setw(15) ;
                            //cout << " pBin =   " << j           << setw(25) <<  " centCalculated = " << centBinLimit << endl;
                            
                            if(usedBins[i][j][k][l] != 1) continue;
                            
                            //if(mmm %100 == 0)
                            //{
                              cout<<"still integrals "<< i<<" "<< j <<" "<<k<<"  "<<centInput<<"  "<<l<<endl;
                            //}
                                
                            etaBin   =  i;
                            centBin  =  j;
                            momBin   =  k;
                            idenSign = l;
                            
                            
                            
                            iden4  -> AddIntegrals(useSign);
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
            }

            
            
            
            nEvents   = iden4 -> GetNEvents();
        
            cout<<" =============================== first moments =============================== "<<endl;
            
            
            
            cout<<"events "<< nEvents << endl;
            
            
            
            
            cout<<"electron    : "<< el1   <<" int: "<< el1Int*nnorm<<endl;
            cout<<"pion        : "<< pi1   <<" int: "<< pi1Int*nnorm<<endl;
            cout<<"kaon        : "<< ka1   <<" int: "<< ka1Int*nnorm<<endl;
            cout<<"proton      : "<< pr1   <<" int: "<< pr1Int*nnorm<<endl;
            cout<<"countMean   : "<<countMean<<endl;
            cout<<"countSecond : "<<countSecond<<endl;
            cout<<"electron2   : "<< el2  <<endl;
            cout<<"pion2       : "<< pi2  <<endl;
            cout<<"kaon2       : "<< ka2  <<endl;
            cout<<"proton2     : "<< pr2  <<endl;
            momTree -> Fill();
        } //ptBin
    } //etaBin
    
    outFile -> cd();
    momTree -> Write();
    nt -> Write();
    outFile -> Close();
    
    delete outFile; //yeni eklave etdim.
    delete fDebug;
    delete iden4;
    //delete moments;
    return 1;
}


