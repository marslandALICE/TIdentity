#ifndef TIDENTITY2D_H
#define TIDENTITY2D_H
#include "TNamed.h"
#include "TBranch.h"
#include "fstream"
#include "iostream"
#include "TString.h"
#include "TIdentityBase.h"
#include "stdlib.h"
#include "TIdentityFunctions.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <vector>
#include "TString.h"
#include "TH1.h"
#include "TMatrixD.h"
class TFile;
class TTree;
class TH1D;

class TIdentity2D:public TIdentityBase
{
public:
    TIdentity2D();
    TIdentity2D(Int_t size);
    virtual ~TIdentity2D();
    static Double_t GetValue(Double_t *, Double_t *);
    static Double_t GetFunctions(Double_t *, Double_t *);
    static Double_t GetFunctionsMix(Double_t *, Double_t *);
    void InitIden2D(Int_t );
    void CalcMoments();
    void GetMoments();
    Double_t GetSecondMoment(Int_t);
    Double_t GetMixedMoment(Int_t, Int_t);
    Double_t GetNuDyn(Int_t, Int_t);
    Double_t GetMean(Int_t);
    Double_t GetMeanI(Int_t);
    virtual void SetFunctionPointers(fptr);
    virtual void Run();
    void SetFileName(TString  _fileName)   { fileName  = _fileName;  }
    Int_t GetIndex(Int_t, Int_t&, Int_t&);
    Int_t AddEntry();
    Bool_t GetEntry(Int_t);
    void GetTree(Long_t &n, TString idenTreeName);
    TTree *GetTreeFromChain(TString treeList, TString treeName);
    void Finalize();
    void GetBins(const Int_t nExtraBins, Double_t *);

    void InitFunctions();
    //Int_t makeDebug();
    void ResetValues();
    void AddParticles();
    static TIdentityFunctions *functions;

    Double_t GetIntegral(Int_t, Int_t, Int_t);
    Double_t GetIntegralMix(Int_t, Int_t);
    void AddIntegrals(Int_t);
    Double_t GetWI(Int_t, Int_t, Int_t);
    void SetLimits(Float_t, Float_t, Double_t);
    Int_t GetNEvents() const {return fCountVeto;}
    void SetUseSign(Int_t _useSign) {fUseSign = _useSign;}
    void Reset();

    Double_t GetDeDx() { return fMyDeDx; }

    void SetBranchNames(const Int_t tmpNBranches, TString tmpBranchNameArr[])
    {
      fNBranches = tmpNBranches;
      fBranchNames = new TString[fNBranches];
      fBranchVariables = new Float_t[fNBranches];
      for (Int_t i=0; i<fNBranches; i++){
        fBranchNames[i] = tmpBranchNameArr[i];
        fBranchVariables[i] = 0.;
      }
    }

    Double_t GetAverCount(Int_t i) { return fAverCount[i]; }

private:

    Float_t       ffMin;
    Float_t       ffMax;
    Double_t      ffBW;
    TMatrixD     *A;
    //TMatrixD invA;
    Double_t     *B;
    Double_t     *fRecMoments;
    Int_t         fUseSign;
    Int_t         fSign;
    TBranch      *fMyBinBrach;  // just to check which kind of tree is used

    Double_t      WI[10][10];
    Double_t      WI2[10][10];
    Double_t      WIMix[10][10];
    Int_t         fSize_size;
    Int_t         fTSize;
    Int_t         fSizeMatrix;
    Int_t         fTSizeMixed;
    Int_t         fCount;
    Double_t     *fW_sum;
    Double_t     *fAver;
    Double_t     *fAverMixed;
    Double_t     *fAver2;
    Double_t     *fAverI;

    static Int_t  fNParticles;
    static Int_t  fNMixParticles;

    ULong64_t     fPrevEvt;
    ULong64_t     fPrevEvtVeto;
    Int_t         fCountVeto;
    Int_t         fMyBin[3];
    ULong64_t     fEventNum;
    Int_t         fEventNumOldVersion;
    Double_t      fMyDeDx;
    Float_t       fDEdx;
    UInt_t        fCutBit;
    Long_t        fTreeEntries;

    TFile        *TIdentityFile;
    TFile        *debugFile;
    TTree        *TIdentityTree;
    TH1D         *fHistWs[4];
    TH1D         *fHistOmegas[4];
    TString  fileName;
    Char_t TFunctionsName[255];
    TF1 *TFunctions[10];
    TF1 *IFunctions[50][50];
    TF1 *IFunctions2[50][50];
    TF1 *IFunctionsMix[50][50];

    Double_t myIntegral(TF1*);

    std::vector<double> *W;
    std::vector<double> *W2;
    std::vector<double> *Wmixed;

    Int_t fCountPart;
    Int_t fCountPartNeg;
    Int_t fCountPartPos;

    std::vector<double> fCountVec;
    std::vector<double> fCountVec2;
    std::vector<double> fCountVecMix;

    Double_t fAverCount[3];

    Int_t fNBranches;
    TString *fBranchNames;
    Float_t *fBranchVariables;


    ClassDef(TIdentity2D,0)
};
#endif
