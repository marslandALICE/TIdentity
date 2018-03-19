#ifndef TIDENTITY2D_H
#define TIDENTITY2D_H
#include "TNamed.h"
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
    void SetInputDir(TString  _inputDir)   { inputDir  = _inputDir;  }
    void SetOutputDir(TString _outputDir)  { outputDir = _outputDir; }
    void SetFileName(TString  _fileName)   { fileName  = _fileName;  }
    Int_t GetIndex(Int_t, Int_t&, Int_t&);
    Int_t AddEntry(Bool_t &);
    Bool_t GetEntry(Int_t);
    void GetTree(Long_t &n, TString idenTreeName = "DATA");
    TTree *GetTreeFromChain(TString treeList, TString treeName);
    void Finalize();
    void GetBins(Int_t *);
    
    Float_t GetMomX() { return momX; }
    Float_t GetMomY() { return momY; }
    Float_t GetMomZ() { return momZ; }
    
    Float_t GetwProton() { return wProton; }
    Float_t GetwKaon()   { return wKaon;   }
    Float_t GetwPion()   { return wPion;   }

    
  
    void SetEffPion(Float_t _effPion )      {effPion = _effPion; }
    void SetEffKaon(Float_t _effKaon )      {effKaon = _effKaon; }
    void SetEffProton(Float_t _effProton )  {effProton = _effProton;}
    
    Int_t GetSign() { return sign; }
    
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
    Int_t GetNEvents() const {return countVeto;}
    void SetUseSign(Int_t _useSign) {useSign = _useSign;}
    void Reset();
    
    Double_t GetDeDx() { return myDeDx; }
    
    void SetMaxPt(Float_t _maxPt) {maxPt = _maxPt;}
    void SetMinPt(Float_t _minPt) {minPt = _minPt;}
    void SetMaxMom(Float_t _maxMom) {maxMom = _maxMom; }
    void SetMinMom(Float_t _minMom) {minMom = _minMom; }
    
    Double_t GetAverCount(Int_t i) { return averCount[i]; }
    
private:
    
    Float_t maxPt, minPt, maxMom, minMom;
    
    Float_t wProton;
    Float_t wKaon;
    Float_t wPion;
    Float_t ffMin, ffMax;
    Double_t ffBW;
    TMatrixD *A;
    //TMatrixD invA;
    Double_t *B;
    Double_t *recMoments;
    Int_t useSign;
    Int_t    sign;
    
    Double_t WI[10][10];
    Double_t WI2[10][10];
    Double_t WIMix[10][10];
    Int_t size_size;
    Int_t TSize;
    Int_t sizeMatrix;
    Int_t TSizeMixed;
    Int_t    count;
    Double_t myMom, myPt;
    Double_t *W_sum;
    Double_t *aver;
    Double_t *averMixed;
    Double_t *aver2;
    Double_t *averI;
    Int_t prevEvt;
    Int_t prevEvtVeto;
    Int_t countVeto;
    Int_t    myBin[3];
    Int_t    evtNum;
    Double_t  myDeDx;
    Float_t  momX;
    Float_t  momY;
    Float_t  momZ;
    
    Float_t effPion;
    Float_t effProton;
    Float_t effKaon;
    
    Long_t    nEntries;
    TFile*   TIdentityFile;
    TFile*   debugFile;
    TTree*   TIdentityTree;
    TH1D*    histoBin;
    TString  inputDir;
    TString  outputDir;
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
    
    Int_t countPart;
    Int_t countPartNeg;
    Int_t countPartPos;
    
    std::vector<double> countVec;
    std::vector<double> countVec2;
    std::vector<double> countVecMix;
    
    Double_t averCount[3];
    
    
    ClassDef(TIdentity2D,0)
};
#endif
