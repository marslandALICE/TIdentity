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
    //
    //
    static Double_t GetValue(Double_t *, Double_t *);
    static Double_t GetFunctions(Double_t *, Double_t *);
    static Double_t GetFunctionsMix(Double_t *, Double_t *);
    static TIdentityFunctions *fFunctions;
    //
    //
    void            InitIden2D(Int_t );
    void            CalcMoments();
    void            GetMoments();
    void            SetFileName(TString  _fileName)   { fTIdenFileName  = _fileName;  }
    void            GetTree(Long_t &n, TString idenTreeName);
    void            Finalize();
    void            GetBins(const Int_t nExtraBins, Double_t *);
    void            InitFunctions();
    void            ResetValues();
    void            AddParticles();
    void            AddIntegrals(Int_t);
    void            SetLimits(Float_t, Float_t, Double_t);
    void            SetUseSign(Int_t _useSign) {fUseSign = _useSign;}
    void            Reset();
    //
    //
    virtual void    SetFunctionPointers(fptr);
    virtual void    Run();
    //
    //
    Double_t        GetSecondMoment(Int_t);
    Double_t        GetMixedMoment(Int_t, Int_t);
    Double_t        GetNuDyn(Int_t, Int_t);
    Double_t        GetMean(Int_t);
    Double_t        GetMeanI(Int_t);
    Double_t        myIntegral(TF1*);
    Double_t        GetIntegral(Int_t, Int_t, Int_t);
    Double_t        GetIntegralMix(Int_t, Int_t);
    Double_t        GetWI(Int_t, Int_t, Int_t);
    Double_t        GetDeDx() { return fMyDeDx; }
    Double_t        GetAverCount(Int_t i) { return fAverCount[i]; }
    Int_t           GetIndex(Int_t, Int_t&, Int_t&);
    Int_t           AddEntry();
    Int_t           GetNEvents() const {return fCountVeto;}
    Bool_t          GetEntry(Int_t);
    TTree          *GetTreeFromChain(TString treeList, TString treeName);
    //
    //
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

    TFile        *fTIdentityFile;
    TFile        *fDebugFile;
    TTree        *fTIdentityTree;
    TH1D         *fHistWs[4];
    TH1D         *fHistOmegas[4];
    TString       fTIdenFileName;
    Char_t        fTFunctionsName[255];
    TF1          *fTFunctions[10];
    TF1          *fIFunctions[50][50];
    TF1          *fIFunctions2[50][50];
    TF1          *fIFunctionsMix[50][50];

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
