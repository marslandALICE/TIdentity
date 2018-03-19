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
#include "iostream"
#include "TH3F.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "string"

using namespace std;
using std::cout;
using std::setw;

// =======================================================================================================
// Helper Functions
TString   PrintNumInBinary(UInt_t num);
Bool_t    ApplyTreeSelection(UInt_t cut);
Double_t  fitFunctionGenGaus(Double_t *x, Double_t *par);
void      InitializeObjects(Int_t sysSet);
void      readFitParams(TString paramTreeName, Int_t fNthFitIteration);
Double_t  EvalFitValue(Int_t particle, Double_t x);
//
// =======================================================================================================
//
// enums
enum particles{electron=0, pion=1, kaon=2, proton=3};
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
const Int_t fnSignBins     = 3;
const Int_t fnEtaBins      = 16;  // MC: 8, Data:16
const Int_t fnCentBins     = 9;
const Int_t fnParticleBins = 4;   // keep it always 4
const Int_t fnMomBins      = 150;
//
//
// Look up table related
const Int_t nBinsLineShape      = 2000;
Bool_t      fTestMode           = kFALSE;
Bool_t      lookUpTableForLine  = kFALSE;
Int_t       lookUpTableLineMode = 0;     // 0 for TH1D, 1 for TF1
//
//
const Float_t fEtaRangeDown = -0.8;
const Float_t fEtaRangeUp   = 0.8;
const Float_t fMomRangeDown = 0.2;
const Float_t fMomRangeUp   = 3.2;
Float_t xCentBins[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
Int_t xSignBins[] = {-1,1,0}; // coding of signs
//
//
Int_t   fNthFitIteration = 6;
TString treeLineShapes   = "treeId";
TString treeIdentity     = "tracks";   // for data "fIdenTree",    for MC "fIdenTreeMC"  , new data "tracks"
//
//
const Int_t fnCutBins=7;
Int_t fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
/*
    if (sysSet==1)  fCutArr[fnCutBins] = {kCRows60,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    if (sysSet==2)  fCutArr[fnCutBins] = {kCRows100, kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};

    if (sysSet==3)  fCutArr[fnCutBins] = {kCRows80,  kChi2TPC3, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    if (sysSet==4)  fCutArr[fnCutBins] = {kCRows80,  kChi2TPC5, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};

    if (sysSet==5)  fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXYSmall, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    if (sysSet==6)  fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXYLarge, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};

    if (sysSet==7)  fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZSmall, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    if (sysSet==8)  fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZLarge, kEventVertexZ, kClusterRequirementITS, kNewITSCut};

    if (sysSet==9)  fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZSmall, kClusterRequirementITS, kNewITSCut};
    if (sysSet==10) fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZSmall, kClusterRequirementITS, kNewITSCut};

    if (sysSet==11) fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kNewITSCut,             kCRows80};
    if (sysSet==12) fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kCRows80};
*/
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
// From Anar hoca
Double_t nEvents = 0;
Double_t nnorm   = 1.;
Double_t fCountMean;
Double_t fCountSecond;
Double_t fCountMix;
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
Int_t fEtaBin, fCentBin, fMomBin, fSignBin;
UInt_t fCutBit;
Int_t fUsedBins[fnEtaBins][fnCentBins][fnMomBins][fnSignBins];
//
Double_t prpiNuDyn=0., prkaNuDyn=0., pikaNuDyn=0., pipiNuDyn=0., kakaNuDyn=0., prprNuDyn=0.;
Double_t el1=0.,  pi1=0.,  ka1=0.,  pr1=0.,  de1=0.,  mu1=0.;
Double_t el1Int=0.,  pi1Int=0.,  ka1Int=0.,  pr1Int=0;
Double_t el2=0.,  pi2=0.,  ka2=0.,  pr2=0.,  de2=0.,  mu2=0.;
Double_t elpi=0., elka=0., elpr=0., elde=0., elmu=0., pika=0., pipr=0., pide=0.;
Double_t pimu=0., kapr=0., kade=0., kamu=0., prde=0., prmu=0., demu=0.;
//
//
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{

    //
    // 4 arguments should be given ./testIden <dEdxListName> <inputfileNameLineShapes> cent fSubsample
    //

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
    }
    else
    {
        cout<<" main.Error: wrong input list"<<endl;
    }
    // OutputFile
    outFile = new TFile(Form("TIMoments_sub%d_cent_%3.2f_mom_%3.2f_%3.2f_eta_%3.2f_%3.2f.root",fSubsample,fCentInput,fpDownInput,fpUpInput,fEtaDownInput,fEtaUpInput),"recreate");
    //
    // Initialize objects and get the bin information
    InitializeObjects(fSystematic);
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
    //
    TROOT IdentityMethod("IdentityMethod","compiled identity method");
    cout << " ================================================================================= " << endl;
    cout << " main.Info: Inputs: " << endl;
    cout << " main.Info: data Tree             = " << fileNameDataTree       << endl;
    cout << " main.Info: Line Shapes           = " << fileNameLineShapes     << endl;
    cout << " main.Info: Centrality            = " << fCentInput      << endl;
    cout << " main.Info: Subsample index       = " << fSubsample      << endl;
    cout << " main.Info: Eta Range             = " << fEtaDownInput       << " - " <<  fEtaUpInput << endl;
    cout << " main.Info: Momentum Range        = " << fpDownInput         << " - " <<  fpUpInput << endl;
    cout << " main.Info: Fit iteration         = " << fNthFitIteration    << endl;
    cout << " main.Info: Charge                = " << fSignInput         << endl;
    cout << " ================================================================================= " << endl;
    cout << " main.Info: Centrality (Bin)      = " << fCentInputBin      << endl;
    cout << " main.Info: Eta Range (Bin)       = " << etaRange[0]    << " - " <<  etaRange[1] << endl;
    cout << " main.Info: Momentum Range (Bin)  = " << momRange[0]    << " - " <<  momRange[1] << endl;
    cout << " ================================================================================= " << endl;
    //
    //
    readFitParams(fileNameLineShapes,fNthFitIteration);
    fgengaus = new TF1("fgengaus",fitFunctionGenGaus,0,1020,5);
    //
    // Create the TIdentity2D object and start analysis
    TIdentity2D *iden4 = new TIdentity2D(4);      // Set the number of particles to 4
    iden4 -> SetFileName(fileNameDataTree);
    iden4 -> SetFunctionPointers(EvalFitValue);
    iden4 -> SetLimits(0.,1020.,10.); // --> (dEdxMin,dEdxMax,binwidth), if slice histograms are scaled wrt binwidth, then binwidth=1
    iden4 -> SetUseSign(fSignInput);  // pass input sign value to TIdentity module
    Long_t nEntries;
    cout<<" main.Info: Running for data"<<endl;
    iden4 -> GetTree(nEntries,treeIdentity);
    iden4 -> Reset();
    //
    //
    if (fTestMode) nEntries = 10000000;
    Float_t bins[7];
    //
    //     bins[0] = eta;
    //     bins[1] = cent;
    //     bins[2] = ptot;
    //     bins[3] = sign;
    //     bins[4] = cutBit;
    //     bins[5] = cRows;
    //     bins[6] = tpcchi2;
    //
    // track by track loop --> read all track info  and add tracks to the iden4 object
    for( Int_t i = 0; i < nEntries; i++ )
    {
        if( !iden4 ->  GetEntry(i) ) continue;
        iden4      ->  GetBins( bins );    // reads identity tree and retrives mybin[] info
        fEtaBin  = fhEta  -> FindBin(bins[0] + 0.0000001) -1;
        fCentBin = fhCent -> FindBin(bins[1] + 0.0000001) -1;
        fMomBin  = fhPtot -> FindBin(bins[2] + 0.0000001) -1;
        fSignBin = (Int_t)bins[3]+1;  // neg --> 0, neutrol --> 1, pos --> 2
        fCutBit  = (UInt_t)bins[4];
        //
        //
        if (!ApplyTreeSelection(fCutBit)) continue;
        fChi2->Fill(bins[6]);
        fcRows->Fill(bins[5]);
        if (i%1000000==0){
            TString cutBinary = PrintNumInBinary(fCutBit);
            cout << cutBinary << "  "<< bins[0] << "  " << bins[1] << "  " << bins[2] << "  " << bins[3] << "  " << bins[4] << endl;
        }
        //
        //
        if( fCentBin != fCentInputBin ) continue;
        if( fMomBin < momRange[0] || fMomBin > momRange[1] ) continue;
        if( fEtaBin < etaRange[0] || fEtaBin > etaRange[1] ) continue;
        fUsedBins[fEtaBin][fCentBin][fMomBin][fSignBin] = 1;
        if (fSignInput==0) fUsedBins[fEtaBin][fCentBin][fMomBin][fSignBin] = 1;
        //
        Bool_t isAdd = kFALSE;
        iden4 -> AddEntry(isAdd);
    }
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
    // First Moments
    el1 = iden4 -> GetMean(electron);
    pi1 = iden4 -> GetMean(pion);
    ka1 = iden4 -> GetMean(kaon);
    pr1 = iden4 -> GetMean(proton);
    //
    // Second Moments
    el2 = iden4 -> GetSecondMoment(electron);
    pi2 = iden4 -> GetSecondMoment(pion);
    ka2 = iden4 -> GetSecondMoment(kaon);
    pr2 = iden4 -> GetSecondMoment(proton);
    //
    // Mixed Moments
    elpi = iden4 -> GetMixedMoment(electron,pion);
    elka = iden4 -> GetMixedMoment(electron,kaon);
    elpr = iden4 -> GetMixedMoment(electron,proton);
    pika = iden4 -> GetMixedMoment(pion,kaon);
    pipr = iden4 -> GetMixedMoment(pion,proton);
    kapr = iden4 -> GetMixedMoment(kaon,proton);
    //
    //Integrals:
    el1Int = iden4 -> GetMeanI(electron);
    pi1Int = iden4 -> GetMeanI(pion);
    ka1Int = iden4 -> GetMeanI(kaon);
    pr1Int = iden4 -> GetMeanI(proton);
    //
    // Printing
    nnorm     = pi1/pi1Int;
    nEvents   = iden4 -> GetNEvents();
    cout << " =============================== first moments =============================== "<<endl;
    cout << " events      : "<< nEvents << endl;
    cout << " electron    : "<< el1   <<" int: "<< el1Int*nnorm << "  ratio: " << el1/(el1Int*nnorm) << endl;
    cout << " pion        : "<< pi1   <<" int: "<< pi1Int*nnorm << "  ratio: " << pi1/(pi1Int*nnorm) << endl;
    cout << " kaon        : "<< ka1   <<" int: "<< ka1Int*nnorm << "  ratio: " << ka1/(ka1Int*nnorm) << endl;
    cout << " proton      : "<< pr1   <<" int: "<< pr1Int*nnorm << "  ratio: " << pr1/(pr1Int*nnorm) << endl;
    cout << " electron2   : "<< el2  <<endl;
    cout << " pion2       : "<< pi2  <<endl;
    cout << " kaon2       : "<< ka2  <<endl;
    cout << " proton2     : "<< pr2  <<endl;
    momTree -> Fill();
    cout << "====================================" << endl;
    cout << " main.Info: calculation is finished " << endl;
    timer.Stop(); timer.Print();
    cout << "====================================" << endl;
    //
    // Close file and clear memory
    outFile -> cd();
    fChi2->Write("fChi2");
    fcRows->Write("fcRows");
    momTree -> Write();
    outFile -> Close();
    delete outFile; //yeni eklave etdim.
    delete iden4;
    return 1;
}
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------
void readFitParams(TString paramTreeName, Int_t fNthFitIteration)
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

    } else {
        //
        cout << " Read lane shapes from ttree " << endl;
        treeLookUp = (TTree*)fLineShapesLookUpTable->Get(treeLineShapes);
        Int_t nent = treeLookUp -> GetEntries();
        cout << nent <<  "   tree is fine go ahead " << endl;

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
            signBin = sign+1;

            // Analyse only 4 iteration of the iterative fitting procedure
            if (it != fNthFitIteration) continue;
            fAmpArr[myBin[0]][myBin[1]][myBin[2]][electron][signBin]      = elA;
            fAmpArr[myBin[0]][myBin[1]][myBin[2]][pion][signBin]          = piA;
            fAmpArr[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]          = kaA;
            fAmpArr[myBin[0]][myBin[1]][myBin[2]][proton][signBin]        = prA;

            fMeanArr[myBin[0]][myBin[1]][myBin[2]][electron][signBin]     = elM;
            fMeanArr[myBin[0]][myBin[1]][myBin[2]][pion][signBin]         = piM;
            fMeanArr[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]         = kaM;
            fMeanArr[myBin[0]][myBin[1]][myBin[2]][proton][signBin]       = prM;

            fSigmaArr[myBin[0]][myBin[1]][myBin[2]][electron][signBin]    = elSi;
            fSigmaArr[myBin[0]][myBin[1]][myBin[2]][pion][signBin]        = piSi;
            fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]        = kaSi;
            fSigmaArr[myBin[0]][myBin[1]][myBin[2]][proton][signBin]      = prSi;

            fSkewArr[myBin[0]][myBin[1]][myBin[2]][electron][signBin]     = elSk;
            fSkewArr[myBin[0]][myBin[1]][myBin[2]][pion][signBin]         = piSk;
            fSkewArr[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]         = kaSk;
            fSkewArr[myBin[0]][myBin[1]][myBin[2]][proton][signBin]       = prSk;

            fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][electron][signBin] = elK;
            fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][pion][signBin]     = piK;
            fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kaon][signBin]     = kaK;
            fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][proton][signBin]   = prK;

        }
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
void InitializeObjects(Int_t sysSet)
{

    //
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
    // initialize output tree
    //
    momTree = new TTree("momTree","momTree");
    momTree -> Branch("fileNameDataTree",&fileNameDataTree);
    momTree -> Branch("fileNameLineShapes",&fileNameLineShapes);
    momTree -> Branch("fCountMean",&fCountMean);
    momTree -> Branch("fCountSecond",&fCountSecond);
    momTree -> Branch("fCountMix",&fCountMix);
    momTree -> Branch("sign",&fSignInput);
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
    momTree -> Branch("fpDownInput",&fpDownInput);
    momTree -> Branch("fpUpInput",&fpUpInput);
    momTree -> Branch("fpDownBin",&fpDownBin);
    momTree -> Branch("fpUpBin",&fpUpBin);
    momTree -> Branch("fEtaDownInput",&fEtaDownInput);
    momTree -> Branch("fEtaUpInput",&fEtaUpInput);
    //
    //
    const Int_t nMoments = 13;
    TVectorF moments(nMoments);
    TVectorF momentsPos(nMoments);
    TVectorF momentsNeg(nMoments);
    TVectorF momentsCross(nMoments);
    //
    // initialize counters
    for(Int_t i=0;i<nMoments; i++){
        moments[i]=0.;
        momentsPos[i]=0.;
        momentsNeg[i]=0.;
        momentsCross[i]=0.;
    }

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
Bool_t ApplyTreeSelection(UInt_t cut)
{
    //     UInt_t arr[fnCutBins];
    //     for (Int_t i=0;i<fnCutBins;i++){
    //         arr[i] = ((cut >> fCutArr[i]) & 1);
    //     }
    //
    //     return (arr[0]&&arr[1]&&arr[2]&&arr[3]&&arr[4]&&arr[5]&&arr[6]);

    for (Int_t i=0;i<fnCutBins;i++){
        if( ((cut >> fCutArr[i]) & 1) == 0 ) return kFALSE;
    }

}
// -----------------------------------------------------------------------------------------
Double_t fitFunctionGenGaus(Double_t *x, Double_t *par)
{
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
