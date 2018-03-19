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

using namespace std;
using std::cout;
using std::setw;


TH3F *pionEffPos[9] = {0};
TH3F *protonEffPos[9] = {0};
TH3F *kaonEffPos[9] = {0};


TH3F *pionEffNeg[9] = {0};
TH3F *protonEffNeg[9] = {0};
TH3F *kaonEffNeg[9] = {0};

TFile *fEff[30];

Float_t effPion, effProton, effKaon;

/*
 How to run
 
 Hints
 1) in     TIdentity2D.cc line 111 --> branch names have to be same
 
 cpEbyeMacros
 cd /hera/alice/marsland/pFluct/files/analysis/TIdentity
 make clean; make
 
 cd /hera/alice/marsland/pFluct/files/analysis/TIdentity/example
 cpEbyeMacros
 make clean; make
 
 // alternative to screen session
 nohup ./testIden IdDataTrees_AllStatistics/cent_0/DataTree_0.root         ParamTree.root 0 0  > cent0_10-80pbin.out &
 
 */

enum particles{electron, pion, proton, kaon};

// ======= Modification Part =============================================================================
Int_t useSign = -100;
Bool_t fpp             = kFALSE;
const Int_t nEtaBins   = 16;
const Int_t nMomBins   = 150;
const Int_t nParticles = 4;
Int_t binFactor        = 2;     // 2 for 20Mev mom bin (for MC analysis), 1 for 10Mev bin (for real data analysis)
TString treeName       = "fEventTree";
//   TString treeName = "fEventTreeMC";
//   TString treeName = "dEdxTree";
Bool_t test = kFALSE;
const Int_t nCentBins  = 9;
const Int_t ncentbins = 10; Float_t xCentBins[ncentbins] = {0, 5,  10,  20, 30, 40, 50, 60, 70, 80};
//const Int_t ncentbins = 19; Float_t xCentBins[ncentbins] = {0, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80};


// =======================================================================================================


Double_t amp[nEtaBins][nCentBins][nMomBins][nParticles]; // amp[etaBin][centBin][momBin][particleType], before kurtosis[9][12][1000][6];
Double_t mean[nEtaBins][nCentBins][nMomBins][nParticles];
Double_t sigma[nEtaBins][nCentBins][nMomBins][nParticles];
Double_t skew[nEtaBins][nCentBins][nMomBins][nParticles];
Double_t kurtosis[nEtaBins][nCentBins][nMomBins][nParticles];
Int_t etaBin, centBin, momBin;

// Values to set from outside
TString inDir("/lustre/nyx/alice/users/marsland/pFluct/files/analysis/Data/");
TString outDir("./");


// function and the tree object
TF1   *fgengaus = 0;
TTree *dataTree = 0;
Int_t ntest = 1;
Int_t centBinLimit;            // 0-5-10-20-30-40-50-60-70-80 is given as input and converted to bin info (is given from outside)


void setFitParameters(Int_t i, Int_t eta, Int_t mom, Int_t cent, TF1 *f)
{
    
    cout<<"cent "<<cent<<endl;
    
    f -> SetParameter(0,amp[eta][cent][mom][i]);
    f -> SetParameter(1,mean[eta][cent][mom][i]);
    f -> SetParameter(2,sigma[eta][cent][mom][i]);
    f -> SetParameter(3,kurtosis[eta][cent][mom][i]);
    f -> SetParameter(4,skew[eta][cent][mom][i]);
}


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

void GetEffs(Float_t px, Float_t py, Float_t pz, Int_t charge, Int_t mcent)
{
    Float_t mpt     = sqrt(px*px + py*py);
    Float_t mphi    = TMath::ATan2(py,px);
    if(mphi < 0) mphi += TMath::Pi()*2;
    Float_t mtheta  = TMath::ATan2(mpt,pz);
    //Float_t theta3 = TMath::ACos();
    Float_t meta    = -TMath::Log(TMath::Tan(mtheta/2.));
    
    //cout<<"test "<< mpt << meta << mphi <<endl;
    
    if( charge == -1 )
    {
        //cout<<"negative  "<< mcent <<endl;
        
        //if(!pionEffNeg[mcent]) cout<<"eff histogram is empty"<<endl;
        
        effPion   = pionEffNeg[mcent]   -> GetBinContent(pionEffNeg[mcent]   -> FindBin(mpt, meta, mphi));
        effKaon   = kaonEffNeg[mcent]   -> GetBinContent(kaonEffNeg[mcent]   -> FindBin(mpt, meta, mphi));
        effProton = protonEffNeg[mcent] -> GetBinContent(protonEffNeg[mcent] -> FindBin(mpt, meta, mphi));
    }
    else
        if( charge == 1 )
        {
            effPion   = pionEffPos[mcent]   -> GetBinContent(pionEffPos[mcent]   -> FindBin(mpt, meta, mphi));
            effKaon   = kaonEffPos[mcent]   -> GetBinContent(kaonEffPos[mcent]   -> FindBin(mpt, meta, mphi));
            effProton = protonEffPos[mcent] -> GetBinContent(protonEffPos[mcent] -> FindBin(mpt, meta, mphi));
        }
        else
        {
            effPion = effProton = effKaon = 0.;
        }
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
    g->FixParameter(0,amp[0][centdebug][19][pion]);
    g->FixParameter(1,mean[0][centdebug][19][pion]);
    g->FixParameter(2,sigma[0][centdebug][19][pion]);
    g->FixParameter(3,kurtosis[0][centdebug][19][pion]);
    g->FixParameter(4,skew[0][centdebug][19][pion]);
    
    return g;
    
}
// -----------------------------------------------------------------------------------------
/*
 void GetDeDxTree(TString list){
 
 //
 // Get the dEdxTree (i.e identity tree with branches myBin, mydEdx, ...) either from a root file or from a chain
 // "list" can be either a list of root files containing dEdxTree or a single file
 //
 
 if (list.Contains(".root")) {
 TFile *f = TFile::Open(list);
 dataTree = (TTree*)f->Get(treeName);
 cout << "get dEdxTree from ---------> " << list << endl;
 
 } else {
 TChain *chain = new TChain(treeName);
 TString command="cat ";
 command+=list;
 
 TString pipe(gSystem->GetFromPipe(command.Data()));
 TObjArray *arrLines=pipe.Tokenize("\n");
 TIter nextLine(arrLines);
 TObject *o = 0x0;
 
 // make the chain adding the files
 Int_t filecount = 0;
 while ( o=nextLine() ) {
 cout << " files to be chained -->  " << o->GetName() << endl;
 chain->AddFile(o->GetName());
 filecount++;
 if (test && filecount > ntest) break;
 }
 chain->Lookup();
 dataTree = (TTree*)chain;
 }
 
 }
 */
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
        
        if(sign != useSign ) continue;
        
        if (fpp) myBin[1]=0;
        
        // Analyse only 4 iteration of the iterative fitting procedure
        if (it != analyseIter) continue;
        
        amp[myBin[0]][myBin[1]][myBin[2]][electron] = elA;
        amp[myBin[0]][myBin[1]][myBin[2]][pion]     = piA;
        amp[myBin[0]][myBin[1]][myBin[2]][kaon]     = kaA;
        amp[myBin[0]][myBin[1]][myBin[2]][proton]   = prA;
        
        mean[myBin[0]][myBin[1]][myBin[2]][electron] = elM;
        mean[myBin[0]][myBin[1]][myBin[2]][pion]     = piM;
        mean[myBin[0]][myBin[1]][myBin[2]][kaon]     = kaM;
        mean[myBin[0]][myBin[1]][myBin[2]][proton]   = prM;
        
        sigma[myBin[0]][myBin[1]][myBin[2]][electron] = elSi;
        sigma[myBin[0]][myBin[1]][myBin[2]][pion]     = piSi;
        sigma[myBin[0]][myBin[1]][myBin[2]][kaon]     = kaSi;
        sigma[myBin[0]][myBin[1]][myBin[2]][proton]   = prSi;
        
        skew[myBin[0]][myBin[1]][myBin[2]][electron] = elSk;
        skew[myBin[0]][myBin[1]][myBin[2]][pion]     = piSk;
        skew[myBin[0]][myBin[1]][myBin[2]][kaon]     = kaSk;
        skew[myBin[0]][myBin[1]][myBin[2]][proton]   = prSk;
        
        kurtosis[myBin[0]][myBin[1]][myBin[2]][electron] = elK;
        kurtosis[myBin[0]][myBin[1]][myBin[2]][pion]     = piK;
        kurtosis[myBin[0]][myBin[1]][myBin[2]][kaon]     = kaK;
        kurtosis[myBin[0]][myBin[1]][myBin[2]][proton]   = prK;
        
    }
    
}
// -----------------------------------------------------------------------------------------
Double_t EvalFitValue(Int_t i, Double_t x)
{
    
    // i stands for particle type
    fgengaus -> SetParameter(0,amp[etaBin][centBin][momBin][i]);
    fgengaus -> SetParameter(1,mean[etaBin][centBin][momBin][i]);
    fgengaus -> SetParameter(2,sigma[etaBin][centBin][momBin][i]);
    fgengaus -> SetParameter(3,kurtosis[etaBin][centBin][momBin][i]);
    fgengaus -> SetParameter(4,skew[etaBin][centBin][momBin][i]);
    return fgengaus->Eval(x);
    
}
// -----------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    
    // Pt ranges
    
    Float_t pt_range[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5};
    Int_t  binnumPt = sizeof(pt_range)/sizeof(Float_t) - 1;
    //TH1F*  hTempPt = new TH1F("hTempPt","", binnumPt, pt_range);
    
    TH1F*  hTempPt = new TH1F("hTempPt","", 26, 0.2, 1.5);
    
    TH1F*  hTempEta = new TH1F("hTempEta","",16, -0.8, 0.8);
    TH1F*  hTempPhi = new TH1F("hTempPhi","",25, 0, TMath::TwoPi());
    
    
    // for(Int_t i = 0; i < 27; i++)
    // {
    //     Float_t x = 0.2 + i*0.05;
    //     pt_range[i] = x;
    // }
    //= {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5};
    //
    //Reading eff matrices:
    /*
    TString effDir = "/hera/alice/rustamov/analysis/fluct/TIdentity/example/";
    Int_t effCent[] = {0, 5, 10, 20, 30, 40, 50, 60, 70};
    
    //TFile *fEff = NULL;
    TString pidName[3] = {"pion","kaon","proton"};
    TString effPid;
    Int_t k = 0;
    for(Int_t i = 0; i < 9; i++)
    {
        for(Int_t j = 0; j <3; j++)
        {
            k++;
            effPid = pidName[j];
            
            effDir.Append("eff_");
            effDir.Append(effPid.Data());
            effDir.Append("_");
            effDir += effCent[i];
            effDir.Append("_.root");
            cout<<"eff pion   "<< effDir <<endl;
            
            //if(fEff) delete fEff;
            fEff[k] = new TFile(effDir);
            
            if(j == 0)
            {
                pionEffPos[i] = new TH3F(*(TH3F*)fEff[k] -> Get("effPos"));
                pionEffNeg[i] = new  TH3F(*(TH3F*)fEff[k] -> Get("effNeg"));
                
                
                //cout<<"testtest "<<pionEffNeg[i] -> GetBinContent(pionEffNeg[i] -> FindBin(0.728296, -0.133788, 0.319674))<<endl;
                
                
            }
            else
                if(j == 1)
                {
                    kaonEffPos[i] = new  TH3F(*(TH3F*)fEff[k] -> Get("effPos"));
                    kaonEffNeg[i] = new  TH3F(*(TH3F*)fEff[k] -> Get("effNeg"));
                }
                else
                {
                    protonEffPos[i] = new  TH3F(*(TH3F*)fEff[k] -> Get("effPos"));
                    protonEffNeg[i] = new  TH3F(*(TH3F*)fEff[k] -> Get("effNeg"));
                }
            effDir = "/hera/alice/rustamov/analysis/fluct/TIdentity/example/";
        }
    }
    */
    
    
    
    //gSystem->Exec("cp /hera/alice/marsland/pFluct/files/analysis/TIdentity/example/testIdentity.C .");
    
    //    Read the entries from out
    Char_t fileName[255];      //file name of tree
    Char_t fitFileName[255];   // file name for fit functions
    Int_t centInput;
    Int_t subsample;
    Int_t analyseIter;
    Float_t etaDown;
    Float_t etaUp;
    Float_t pDown = 0.2;
    Float_t pUp   = 1.5;
    Int_t ptBin;
    TString selectRapidity = "";
    
    // 4 arguments should be given ./testIden <dEdxListName> <fitFileName> cent subsample
    
    cout<<"NUMBER OF ARGUMENTS "<<argc<<endl;
    
    Int_t isSim = 0;
    Int_t isIntegrated = 1;
    Int_t rapInterval = 0;
    Int_t systematics = 0;
    Int_t momInterval = 0;
    
    if(argc == 12)
    {
        sprintf(fileName,"%s",argv[1]);
        sprintf(fitFileName,"%s",argv[2]);
        centInput  = atoi(argv[3]);
        subsample  = atoi(argv[4]);
        analyseIter= atoi(argv[5]);
        //etaDown    = atof(argv[6]);
        //etaUp      = atof(argv[7]);
        // pDown      = atof(argv[8]);
        //ptBin        = atoi(argv[8]);
        //pUp        = atof(argv[9]);
        useSign    = atoi(argv[6]);
        //selectRapidity = argv[9];
        
        systematics = atoi(argv[7]);
        isIntegrated = atoi(argv[8]);
        rapInterval=atoi(argv[9]);
        momInterval = atoi(argv[10]);
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
    
    // Convert real eta and mom info to bin information
    //Int_t etaDownBin;
   // Int_t etaUpBin;
    
    etaUp = 0.8;
    etaDown = -0.8;
    
    pDown = 0.3;//0.6;
    pUp   = 0.8;
    
    /*
    if (isSim)
    {
        etaDownBin = Int_t((etaDown*10.)/2.+4.);
        etaUpBin   = Int_t((etaUp*10.)/2.+4.);
    }
    else
    {
        etaDownBin = Int_t((etaDown*10.)+8.);
        etaUpBin   = Int_t((etaUp*10.)+8.);
    }
     
    
    
    
    
    
    Int_t pDownBin   = Int_t((pDown*100.)-20.)/binFactor;
    Int_t pUpBin     = Int_t((pUp*100.)-20.)/binFactor+1;
    */
    // Eta range in bins;
    // p[0,367] --> [0,300] binwidth is 10MeV rest is 100MeV. which mean [0,300]-->[0.2,3.2],
    //          --> very clean window --> [19,20]-[0,39,0.4]
    // eta[0,8] --> [0,8]   binwidth is 0.2, which means [0,8] --> [-0.8,0.8]
    
    
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
    
    Int_t etaDownBin = fhEta -> FindBin(etaDown + 0.0000001)  - 1;
    Int_t etaUpBin   = fhEta -> FindBin(etaUp   - 0.0000001)  - 1;
    
    
    Int_t etaRange[2] = {etaDownBin,etaUpBin};
    Int_t momRange[2] = {pDownBin  ,pUpBin};     // 0-367 is whole momentum range each bein is 10 MeV up to 300
    
    cout << " Inputs: " << endl;
    cout << fileName    << endl;
    cout << fitFileName << endl;
    cout << "cent        = " << centInput   << endl;
    cout << "subsample   = " << subsample   << endl;
    cout << "iter        = " << analyseIter << endl;
    cout << "realEtaDown = " << etaDown     << "   ===   realEtaUp  = " << etaUp      << endl;
    cout << "etaDownBin  = " << etaDownBin  << "   ===   etaUpBin   = " << etaUpBin   << endl;
    cout << "realPDown   = " << pDown       << "   ===   realPUp    = " << pUp        << endl;
    cout << "pDownBin    = " << pDownBin    << "   ===   pUpBin     = " << pUpBin     << endl;
    cout<< "selectRapidity = "<< selectRapidity<<endl;
    
    TROOT IdentityMethod("IdentityMethod","compiled identity method");
    
    // prepare histogram to find the centrality bin
    TH1D * hCent = new TH1D( "hCent","Centrality Bins",ncentbins-1,xCentBins );
    
    // full path to data tree list
    TString tmpFileName(fileName);
    tmpFileName.Prepend(inDir);
    
    // Full path to param tree file
    TString tmpfitFileName(fitFileName);
    TString paramTreeName(inDir);
    paramTreeName.Append(tmpfitFileName);
    
    cout << " ================================================================================= " << endl;
    cout << " Fit iteration         = " << analyseIter    << endl;
    cout << " Eta Range             = " << etaRange[0]    << " - " <<  etaRange[1] << endl;
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
    //GetDeDxTree(tmpFileName);                     // Get the data tree which contains dEdxTree
    //iden4 -> SetInputTree(dataTree);              // Set the data tree in order to be read by TIdentity module
    //iden4 -> SetFileName(fileName);
    iden4 -> SetInputDir(inDir);
    iden4 -> SetOutputDir(outDir);
    iden4 -> SetFileName(fileName);
    iden4 -> SetFunctionPointers(EvalFitValue);
    iden4 -> SetLimits(0.,1020.,1.); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
    iden4 -> SetUseSign(useSign);
    Long_t nEntries;
    if(isSim)
    iden4 -> GetTree(nEntries,"sim");
    else
    iden4 -> GetTree(nEntries);
    
    // Tree for the final results
    Double_t prpiNuDyn=0., prkaNuDyn=0., pikaNuDyn=0., pipiNuDyn=0., kakaNuDyn=0., prprNuDyn=0.;
    Double_t el1=0.,  pi1=0.,  ka1=0.,  pr1=0.,  de1=0.,  mu1=0.;
    Double_t el2=0.,  pi2=0.,  ka2=0.,  pr2=0.,  de2=0.,  mu2=0.;
    Double_t elpi=0., elka=0., elpr=0., elde=0., elmu=0., pika=0., pipr=0., pide=0.;
    Double_t pimu=0., kapr=0., kade=0., kamu=0., prde=0., prmu=0., demu=0.;
    
    // Get the centrality bin information using hCent
    Int_t centBinLimit  = (fpp) ? 0 : hCent->FindBin(centInput+0.25)-1;
    Double_t centrality = hCent->GetBinCenter(centBinLimit);
    
    // File for the moments output which are dumped into a tree
    
    
    
    // Calculate First moments
    Int_t bins[4];
    
    // -------------- DEBUGGING -----------------
    TF1 *fDebug;
    fDebug = (fpp) ? CheckFitParams(centBinLimit+1) : CheckFitParams(centBinLimit) ;
    // ------------------------------------------
    
    // Write results into tree
    TFile *outFile;
    if(useSign == 0)
        outFile = new TFile(Form("TImoments_%d_%d.root",centInput, subsample),"recreate");
    else
        if(useSign == 1)
            outFile = new TFile(Form("TImoments_%d_%d_pos.root",centInput, subsample),"recreate");
        else
            outFile = new TFile(Form("TImoments_%d_%d_neg.root",centInput, subsample),"recreate");
    
    
    Double_t nEvents = 0;
    
    
    
    
    Int_t ptEffBin = 0;
    Int_t etaEffBin = 0;
    Int_t phiEffBin = 0;
    
    TTree *momTree = new TTree("momTree","momTree");
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
    momTree -> Branch("phibin",&phiEffBin);
    
    
    /////momTree -> Branch
    //momTree -> SetBranch("kapr",&pi2);
    //momTree -> SetBranch("kapr",&pi2);
    
            iden4 -> Reset();
    
    
    TH1F *dedxHisto = new TH1F("dedxHisto","dedxHisto",4000,20,1020);
    TH1F *dedxHistoUp = new TH1F("dedxHistoUp","dedxHistoUp",400,0,1000);
    TH1F *dedxHistoDown = new TH1F("dedxHistoDown","dedxHistoDown",400,0,1000);
    
    
    
    TNtuple *nt1 = new TNtuple("nt1","nt1","wpi:wpr:wk");
    TNtuple *nt2 = new TNtuple("nt2","nt2","WPI:WPR:WK");
    
    Int_t fixEtaBin  = 10;
    Int_t fixMomBin  = 10;
    Int_t centFix    = 0;
    
    
    iden4 -> SetMinPt(0.);
    iden4 -> SetMaxPt(15.);
    
    
    
    iden4 -> SetMinMom(pDown);
    iden4 -> SetMaxMom(pUp);
    
   // nEntries = 10000;
    
    Float_t WPI, WPR, WK;
    
    WPI = WPR = WK = 0.;
    
    for(Int_t i = 0; i < nEntries; i++)
                {
                    if(!iden4 ->  GetEntry(i)) continue;
                    iden4     ->  GetBins(bins);
                    ////
                    
                    etaBin  = bins[0];
                    centBin = bins[1];
                    momBin  = bins[2];
                    
                    if( etaBin < 0 || centBin < 0 || momBin < 0 ) {cout<<" this should not happen "<< etaBin <<"  "<<centBin<<"  "<<momBin << endl; continue;}
                    
                    
                    if( momBin < momRange[0] ||  momBin > momRange[1] ) continue;
                    
                    
                    if( etaBin < etaRange[0] || etaBin > etaRange[1] ) continue;
                    
                    
                    
                    //cout<<"centbin == "<<centBin<<endl;
                    
                    if(!centFix) centFix = centBin;
                    
                    if(momBin == fixMomBin && etaBin == fixEtaBin)
                    dedxHisto -> Fill(iden4 -> GetDeDx());
                    
                    if(momBin == fixMomBin + 1 && etaBin == fixEtaBin +1)
                        dedxHistoUp -> Fill(iden4 -> GetDeDx());
                    
                    if(momBin == fixMomBin - 1 && etaBin == fixEtaBin - 1)
                        dedxHistoDown -> Fill(iden4 -> GetDeDx());
                    
                    Bool_t isAdd = kFALSE;
                    
                    //if(momBin == fixMomBin && etaBin == fixEtaBin)
                    {
                       iden4 -> AddEntry(isAdd);
                        
                        Float_t wwp  = iden4 -> GetwProton();
                        Float_t wwk  = iden4 -> GetwKaon();
                        Float_t wwpi = iden4 -> GetwPion();
                        nt1 -> Fill(wwpi, wwp, wwk);
                        
                        if(isAdd)
                        {
                            if(wwpi > -50 && wwp > -50 && wwk > -50 )
                            {
                               WPI += wwpi;
                               WPR += wwp;
                               WK  += wwk;
                            }
                        }
                        else
                        {
                            nt2 -> Fill(WPI, WPR, WK );
                            WPI = WPR = WK  = 0.;
                        }
                    }
                }
                
    TF1 *fpion = new TF1("fpion",fitFunctionGenGaus,0,1000,5);
    fpion->SetNpx(10000);
    setFitParameters(1, fixEtaBin, fixMomBin, centFix,fpion);
    
    TF1 *fproton = new TF1("fproton",fitFunctionGenGaus,0,1000,5);
    fproton -> SetNpx(10000);
    setFitParameters(2, fixEtaBin, fixMomBin, centFix, fproton);
    
    TF1 *fkaon = new TF1("fkaon",fitFunctionGenGaus,0,1000,5);
    fkaon -> SetNpx(10000);
    setFitParameters(3, fixEtaBin, fixMomBin, centFix, fkaon);
    
    TF1 *felectron = new TF1("felectron",fitFunctionGenGaus,0,1000,5);
    felectron -> SetNpx(10000);
    setFitParameters(0, fixEtaBin, fixMomBin, centFix, felectron);
    
    
    dedxHisto -> Write();
    dedxHistoDown -> Write();
    dedxHistoUp -> Write();
    fpion -> Write();
    fproton -> Write();
    fkaon -> Write();
    felectron -> Write();
    nt1 -> Write();
    nt2 -> Write();
    delete fDebug;
    delete iden4;
    //delete moments;
    return 1;
}


