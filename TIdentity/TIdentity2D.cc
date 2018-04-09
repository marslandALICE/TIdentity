//Author A. Rustamov
//TIdentity class for calculation of all second moments
#include "TIdentity2D.h"
#include "TIdentityFunctions.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"

using namespace std;
ClassImp(TIdentity2D)

Int_t TIdentity2D::fNParticles=4;
Int_t TIdentity2D::fNMixParticles=6;

TIdentityFunctions* TIdentity2D::functions = NULL;

TIdentity2D::TIdentity2D()
{
    InitIden2D(4);
}

TIdentity2D::TIdentity2D(Int_t size)
{

    InitIden2D(size);
}

void TIdentity2D::InitIden2D(Int_t size)
{
    cout<<" _________________________________________________________________________________"<<endl;
    cout<<"|                                                                                 |"<<endl;
    cout<<"|  *   ******     *******    *      *    *********   *   *********   *         *  |"<<endl;
    cout<<"|  *   *      *   *          * *    *        *       *       *        *       *   |"<<endl;
    cout<<"|  *   *      *   *          *  *   *        *       *       *         *     *    |"<<endl;
    cout<<"|  *   *      *   * *****    *   *  *        *       *       *          * * *     |"<<endl;
    cout<<"|  *   *      *   *          *    * *        *       *       *             *      |"<<endl;
    cout<<"|  *   *      *   *          *     **        *       *       *            *       |"<<endl;
    cout<<"|  *   ******     ********   *      *        *       *       *         * *        |"<<endl;
    cout<<"|                                                                                 |"<<endl;
    cout<<"| based on:                                                                       |"<<endl;
    cout<<"| 1. PRC 84, 024902 (2011)                                                        |"<<endl;
    cout<<"| 2. PRC 86, 044906 (2012)                                                        |"<<endl;
    cout<<"| see also: PRC 83, 054907 (2011)                                                 |"<<endl;
    cout<<"|                                                                                 |"<<endl;
    cout<<"|contact: a.rustamov@cern.ch                                                      |"<<endl;
    cout<<"|_________________________________________________________________________________|"<<endl;
    cout<<" "<<endl;

    debugFile = new TFile("TIdenDebug.root","recreate");
    for (Int_t i=0;i<size;i++) {
      fHistWs[i]     = new TH1D(Form("hW_%d",i),Form("hW_%d",i),900,0.,30.);
      fHistOmegas[i] = new TH1D(Form("hOmega_%d",i),Form("hOmega_%d",i),500,0.,1.);
    }
    countPart    = 0;
    countPartNeg = 0;
    countPartPos = 0;

    useSign = -1000;
    sign = -1000;
    myBinBrach=0x0;
    countVeto = 0;
    myDeDx=-1000.;
    dEdx=-1000.;

    prevEvtVeto = -1;
    TSize = size;
    TSizeMixed = size*(size-1)/2;
    sizeMatrix = TSize + TSizeMixed;
    W         = new vector<double>     [TSize];
    W2        = new vector<double>     [TSize];
    Wmixed    = new vector<double> [TSizeMixed];
    if(!functions) functions = new TIdentityFunctions();
    W_sum = new Double_t [TSize];
    prevEvt = -1000;
    aver = new Double_t [TSize];
    aver2 = new Double_t [TSize];
    averI = new Double_t [TSize];
    averMixed = new Double_t [TSizeMixed];

    fNParticles = TSize;
    fNMixParticles = TSizeMixed;

    size_size = 3000;
    for(Int_t i = 0; i < 10; i++)
        for(Int_t j = 0; j < 10; j++)
        {
            WI[i][j] = 0.;
            WI2[i][j] = 0.;
            WIMix[i][j] = 0.;
        }
        for(Int_t i = 0; i < TSize; i++ )
        {
            averI[i] = 0;
        }
        A = new TMatrixD(sizeMatrix,sizeMatrix);
        B = new Double_t [sizeMatrix];
        recMoments = new Double_t [sizeMatrix];
        for(Int_t i = 0; i < sizeMatrix; i++)
        {
            recMoments[i] = 0.;
        }
}

void TIdentity2D::Reset()
{
    countVec.clear();
    countVec2.clear();
    countVecMix.clear();
    countPart = countPartNeg = countPartPos = 0;

    for (Int_t i = 0; i < TSize; i++)
    {
        cout<<" TIdentity2D::Reset.Info: resetting vectors"<<endl;
        W[i].clear();
        W2[i].clear();
        aver[i]  = 0;
        aver2[i] = 0;
        averI[i] = 0;
        W_sum[i] = 0;
    }

    for (Int_t i = 0; i < TSizeMixed; i++)
    {
        Wmixed[i].clear();
        averMixed[i] = 0;
    }

    countVeto = 0;
    prevEvtVeto = -1;

    prevEvt = -1000;
    for(Int_t i = 0; i < 10; i++)
        for(Int_t j = 0; j < 10; j++)
        {
            WI[i][j] = 0.;
            WI2[i][j] = 0.;
            WIMix[i][j] = 0.;
        }
        for(Int_t i = 0; i < TSize; i++ )
        {
            averI[i] = 0;
        }

        for(Int_t i = 0; i < sizeMatrix; i++)
        {
            recMoments[i] = 0.;
            B[i] = 0;
        }

}


TIdentity2D::~TIdentity2D()
{
    if(W) {delete [] W; W = 0;}
    if(W2){ delete [] W2; W2 = 0;}
    if(Wmixed) {delete [] Wmixed; Wmixed = 0; }
    if(W_sum) delete [] W_sum;
    if(aver) delete [] aver;
    if(aver2) delete [] aver2;
    if(averI) delete [] averI;
    if(averMixed) delete [] averMixed;
    if(recMoments) delete [] recMoments;
    if(A) delete A;
    if(B) delete [] B;

}

void TIdentity2D::GetTree(Long_t &nent, TString idenTreeName)
{

    if (fileName.Contains(".root")) {
        TIdentityFile = new TFile(fileName);
        cout<<"TIdentity2D::GetTree.Info: We are reading the file "<<fileName<<endl;
        TIdentityTree = (TTree*)TIdentityFile->Get(idenTreeName);
    } else {
        TIdentityTree = GetTreeFromChain(fileName,idenTreeName);
    }
    if (!TIdentityTree) cout << "TIdentity2D::GetTree.Error: tree could not be read" << endl;
    cout<<"TIdentity2D::GetTree.Info: ======================== "<<endl;
    TIdentityTree -> GetListOfBranches()->ls();
    myBinBrach = (TBranch*)TIdentityTree->FindBranch("myBin");
    cout<<"TIdentity2D::GetTree.Info: ======================== "<<endl;
    if (!myBinBrach){ // new version of tree format
      TIdentityTree -> SetBranchAddress("gid"    ,&evtNum);
      TIdentityTree -> SetBranchAddress("dEdx"   ,&dEdx);
      TIdentityTree -> SetBranchAddress("sign"   ,&sign);
      TIdentityTree -> SetBranchAddress("cutBit" ,&cutBit);
      for (Int_t i=0;i<fNBranches;i++){
        TIdentityTree -> SetBranchAddress(fBranchNames[i] ,&fBranchVariables[i]);
      }
    } else { // old version of tree format
      TIdentityTree -> SetBranchAddress("sign",&sign);
      TIdentityTree -> SetBranchAddress("myBin",myBin);
      TIdentityTree -> SetBranchAddress("myDeDx",&myDeDx);
      TIdentityTree -> SetBranchAddress("evtNum",&evtNum);
    }
    nEntries = (Long_t)TIdentityTree -> GetEntries();
    nent = nEntries;
    InitFunctions();
}

TTree *TIdentity2D::GetTreeFromChain(TString treeList, TString treeName)
{
    cout << "TIdentity2D::GetTreeFromChain.Info: Files added to the chain" << endl;
    ifstream fileTmp(treeList);
    char file[255];
    TChain *chain=NULL;
    chain = new TChain(treeName);
    while(fileTmp) {
        fileTmp.getline(file, sizeof(file));  // delim defaults to '\n'
        if(fileTmp) cout << "TIdentity2D::GetTreeFromChain.Info: " << file << endl;
        chain->Add(file);
    }
    fileTmp.close();
    return (TTree*)chain;
}

Double_t TIdentity2D::GetValue(Double_t *xval, Double_t *par)
{
    Double_t xx = xval[0];
    Int_t index = (Int_t)par[0];
    return functions -> GetValue(index, xx);
}

Double_t TIdentity2D::GetFunctions(Double_t *xx, Double_t *par)
{
    Int_t j = (Int_t)par[0];
    Int_t i = (Int_t)par[1];
    Int_t k = (Int_t)par[2];
    Double_t val[fNParticles];
    Double_t sumVal = 0;
    for(Int_t ii = 0; ii < fNParticles; ii++)
    {
        val[ii] = functions -> GetValue(ii, xx[0]);
        sumVal += val[ii];
    }
    Double_t relVal[fNParticles];
    if(sumVal < 1e-15) return 0.;
    for(Int_t m = 0; m < fNParticles; m++)
    {
        relVal[m] = val[m]/sumVal;
        if( k == 2) relVal[m] *= relVal[m];
    }
    return relVal[j]*val[i];
}

Double_t TIdentity2D::GetFunctionsMix(Double_t *xx, Double_t *par)
{
    Int_t j = (Int_t)par[0];
    Int_t k = (Int_t)par[1];
    Double_t val[fNParticles];
    Double_t sumVal = 0;
    for(Int_t ii = 0; ii < fNParticles; ii++)
    {
        val[ii] = functions -> GetValue(ii, xx[0]);
        sumVal += val[ii];
    }
    if(sumVal < 1e-15) return 0.;
    Double_t myVal[fNMixParticles];
    Int_t t = 0;
    for(Int_t m = 0; m < fNParticles-1; m++)
        for(Int_t n = m+1; n < fNParticles; n++)
        {
            myVal[t] = val[m]*val[n]/sumVal/sumVal;
            t++;
        }
        return myVal[j]*val[k];
}

void TIdentity2D::InitFunctions()
{
    for(Int_t i = 0; i < TSize+1; i++)
    {
        sprintf(TFunctionsName,"func[%d]",i);
        TFunctions[i] = new TF1(TFunctionsName,GetValue, ffMin, ffMax,1);
        TFunctions[i] -> SetParameter(0,i);
    }

    for(Int_t i = 0; i < TSize; i++)
        for(Int_t j = 0; j < TSize; j++)
        {
            sprintf(TFunctionsName,"Ifunc[%d][%d]",i,j);
            IFunctions[i][j] = new TF1(TFunctionsName,GetFunctions,ffMin,ffMax,3);
            IFunctions[i][j] -> SetParameter(0,i);
            IFunctions[i][j] -> SetParameter(1,j);
            IFunctions[i][j] -> SetParameter(2,1);

            sprintf(TFunctionsName,"Ifunc2[%d][%d]",i,j);
            IFunctions2[i][j] = new TF1(TFunctionsName,GetFunctions,ffMin,ffMax,3);
            IFunctions2[i][j] -> SetParameter(0,i);
            IFunctions2[i][j] -> SetParameter(1,j);
            IFunctions2[i][j] -> SetParameter(2,2);
        }
        for(Int_t i = 0; i < TSizeMixed; i++)
            for(Int_t j = 0; j < TSize; j++)
            {
                sprintf(TFunctionsName,"IfuncMix[%d][%d]",i,j);
                IFunctionsMix[i][j] = new TF1(TFunctionsName,GetFunctionsMix,ffMin,ffMax,2);
                IFunctionsMix[i][j] -> SetParameter(0,i);
                IFunctionsMix[i][j] -> SetParameter(1,j);
            }
}

void TIdentity2D::SetFunctionPointers(fptr fun)
{
    functions -> funs = fun;
}

void TIdentity2D::Run()
{
    ;
}

Bool_t TIdentity2D::GetEntry(Int_t i)
{
    if(i%5000000 == 0) cout<<" TIdentity2D::GetEntry.Info: track "<<i<<" of "<<nEntries <<endl;
    TIdentityTree -> GetEntry(i);
    if( evtNum != prevEvtVeto && prevEvtVeto >0) {countVeto++;}
    prevEvtVeto = evtNum;
    if ( dEdx>-1 ) myDeDx=dEdx;
    if ( (myDeDx < ffMin || myDeDx > ffMax) && myDeDx > 0 ) return kFALSE;
    // secure the usage of sign=0 which is sum of + and - particles
    if( !(sign == useSign || useSign == 0) ) return kFALSE;
    return kTRUE;
}

Int_t TIdentity2D::AddEntry()
{

    Bool_t isAdd = kTRUE;
    if(evtNum == prevEvt)
    {
        AddParticles();
    }
    else
    {
        if(count != 0)
        {
            countVec.push_back(countPart);
            countVec2.push_back(countPart*countPart);
            countVecMix.push_back(countPartNeg*countPartPos);

            for(Int_t m = 0; m < TSize; m++)
            {
                W[m].push_back(W_sum[m]);
                W2[m].push_back(W_sum[m]*W_sum[m]);
            }
            //
            // Debug hists
            for(Int_t i = 0; i < TSize; i++) fHistWs[i]->Fill(W[i][W[i].size()-1]);
            //
            Int_t t = 0;
            for(Int_t m = 0; m < TSize-1; m++)
                for(Int_t n = m+1; n < TSize; n++)
                {
                    Wmixed[t].push_back(W_sum[m]*W_sum[n]);
                    t++;
                }
                isAdd = kFALSE;
        }
        ResetValues();
        count = 0;
        AddParticles();
    }
    prevEvt = evtNum;
    return 1;
}

void TIdentity2D::ResetValues()
{
    for(Int_t s = 0; s < TSize; s++)
    {
        W_sum[s] = 0.;
    }
    countPart    = 0;
    countPartNeg = 0;
    countPartPos = 0;
}

void TIdentity2D::AddParticles()
{
    Double_t mValue[fNParticles] = {0.};
    Double_t sumValue = 0;

    for(Int_t i = 0; i < TSize; i++)
    {

        if( myDeDx < 0 ) { mValue[i] = 0; count++; continue; }
        mValue[i] = TFunctions[i] -> Eval(myDeDx);
        sumValue += mValue[i];
    }
    //
    // Debug hists for omega values
    for(Int_t i = 0; i < TSize; i++) fHistOmegas[i]->Fill(mValue[i]/sumValue);
    //
    countPart += 1;
    if(sign == -1)
    {
        countPartNeg += 1;
    }
    if(sign == 1)
    {
        countPartPos += 1;
    }

    if(sumValue > 1e-10)
    {
        count++;
        for(Int_t i = 0; i < TSize; i++)
        {
            W_sum[i] += mValue[i]/sumValue; ;
        }
    }

}

void TIdentity2D::Finalize()
{
    cout<<" "<<endl;
    cout<<" TIdentity2D::Finalize.Info: ***************************************************"<<endl;
    cout<<" TIdentity2D::Finalize.Info: ************ number of analyzed events: "<<countVeto+1<<" ******"<<endl;
    cout<<" TIdentity2D::Finalize.Info: ***************************************************"<<endl;
    cout<<" "<<endl;

    for(Int_t m = 0; m < TSize; m++)
    {
        aver[m]  = accumulate(W[m].begin(),  W[m].end(), 0.0)/countVeto;
        aver2[m] = accumulate(W2[m].begin(), W2[m].end(), 0.0)/countVeto;
    }
    Int_t t = 0;
    for(Int_t m = 0; m < TSize-1; m++)
        for(Int_t n = m+1; n < TSize; n++)
        {
            averMixed[t] = accumulate(Wmixed[t].begin(), Wmixed[t].end(), 0.0)/countVeto;
            t++;
        }

        averCount[0] = accumulate(countVec.begin(),  countVec.end(), 0.0)/countVeto;
    averCount[1] = accumulate(countVec2.begin(),  countVec2.end(), 0.0)/countVeto;
    averCount[2] = accumulate(countVecMix.begin(),  countVecMix.end(), 0.0)/countVeto;

    //
    // Write some output to data
    debugFile->cd();
    for (Int_t i=0;i<TSize;i++) fHistWs[i]     ->Write();
    for (Int_t i=0;i<TSize;i++) fHistOmegas[i] ->Write();
    debugFile -> Close();
    delete debugFile;
    //

}

void TIdentity2D::GetBins(const Int_t nExtraBins, Double_t *bins)
{
  // TString branchNames[nBranches]={"eta","cent","ptot","cRows","tpcchi2"};
  if (!myBinBrach){
    bins[0] = evtNum;
    bins[1] = dEdx;
    bins[2] = sign;
    bins[3] = cutBit;
    for (Int_t i=0;i<nExtraBins;i++) bins[i+4] = fBranchVariables[i];
  } else{
    bins[0] = myBin[0];
    bins[1] = myBin[1];
    bins[2] = myBin[2];
    bins[3] = sign;
  }
}

Double_t TIdentity2D::GetIntegral(Int_t i, Int_t j, Int_t k)
{
    if(k == 1)
    {
        return myIntegral(IFunctions[i][j]);
    }
    else
        return
        myIntegral(IFunctions2[i][j]);
}

Double_t TIdentity2D::GetIntegralMix(Int_t i, Int_t j)
{
    return myIntegral(IFunctionsMix[i][j]);
}

void TIdentity2D::SetLimits(Float_t min, Float_t max, Double_t BW)
{
    ffMin = min;
    ffMax = max;
    ffBW = BW;
}

Double_t TIdentity2D::myIntegral(TF1 *Fun)
{

    Double_t xx, sum = 0;
    Double_t step = (ffMax-ffMin)/(2*size_size);
    for(Int_t i = 1; i < 2*size_size+1; i++)
    {
        xx  = ffMin + step*i;
        sum += Fun -> Eval(xx);
    }
    return sum*step;
}

void TIdentity2D::AddIntegrals(Int_t lsign)
{
    if(useSign != lsign) return;

    for(Int_t i = 0; i < TSize; i++)
        for(Int_t j = 0; j < TSize; j++)
        {
            WI[i][j]  += GetIntegral(i,j,1);
            WI2[i][j] += GetIntegral(i,j,2);
        }

        for(Int_t t = 0; t < TSizeMixed; t++)
            for(Int_t j = 0; j < TSize; j++)
            {
                WIMix[t][j] +=  GetIntegralMix(t,j);
            }

            for(Int_t i = 0; i < TSize; i++)
            {
                averI[i] += myIntegral(TFunctions[i])/ffBW;
            }
}

Double_t TIdentity2D::GetWI(Int_t i, Int_t j, Int_t k)
{
    if(k == 1)
        return WI[i][j]/averI[j]/ffBW;
    else
        if(k ==2)
            return
            WI2[i][j]/averI[j]/ffBW;
        else
            return  WIMix[i][j]/averI[j]/ffBW;
}

void TIdentity2D::CalcMoments()
{
    Int_t t  = TSize;
    Int_t tt = TSize;
    for(Int_t i = 0; i < TSize; i++)
    {
        t = TSize;
        for(Int_t j = 0; j < TSize; j++)
        {
            (*A)(i,j) = GetWI(i,j,1)*GetWI(i,j,1);
            if(j == TSize -1) continue;
            for(Int_t k = j+1; k < TSize; k++)
            {
                (*A)(i,t) = 2.*GetWI(i,j,1)*GetWI(i,k,1);
                t++;
            }
        }
    }

    t = TSize;
    for(Int_t i = 0; i < TSize-1; i++)
        for(Int_t j = i+1; j < TSize; j++)
        {
            tt = TSize;
            for(Int_t k = 0; k < TSize; k++)
            {
                (*A)(t,k) = GetWI(i,k,1)*GetWI(j,k,1);
                if(k == TSize-1) continue;
                for(Int_t kk = k+1; kk < TSize; kk++)
                {
                    (*A)(t,tt) = GetWI(i,k,1)*GetWI(j,kk,1) + GetWI(i,kk,1)*GetWI(j,k,1);
                    tt++;
                }
            }
            t++;
        }

        TMatrixD invA = (*A).Invert();

        for(Int_t i = 0; i < TSize; i++)
        {
            B[i] = aver2[i];
            for(Int_t j = 0; j < TSize; j++)
            {
                B[i] -= aver[j]*(GetWI(i,j,2) - GetWI(i,j,1)*GetWI(i,j,1));
            }
        }
        Int_t indA, indB;
        for(Int_t kk = 0; kk < TSizeMixed; kk++)
        {
            B[kk+TSize] = averMixed[kk];
            GetIndex(kk, indA, indB);
            for(Int_t m = 0; m < TSize; m++)
            {
                B[kk+TSize] -= aver[m]*(GetWI(kk,m,0)-GetWI(indA,m,1)*GetWI(indB,m,1));
            }
        }

        for(Int_t k = 0; k < sizeMatrix; k++)
            for(Int_t tt = 0; tt < sizeMatrix; tt++)
            {
                recMoments[k]  += invA(k,tt)*B[tt];
            }
}

void TIdentity2D::GetMoments()
{
    ;
}

Double_t TIdentity2D::GetSecondMoment(Int_t i)
{
    if(i >= sizeMatrix)
    {
        cout<<"TIdentity2D::GetSecondMoment.Info: out of bound"<<endl;
        return -1000.;
    }
    return recMoments[i];
}

Double_t TIdentity2D::GetMixedMoment(Int_t i, Int_t j)
{
    Int_t tmp = i;
    Int_t ind = TSize-1;
    if( i > j) {i = j; j = tmp;}
    for(Int_t k = 0; k < TSize -1; k++)
    {
        for(Int_t kk = k+1; kk < TSize; kk++)
        {
            ind++;
            if (i == k && j == kk) break;
        }
        if(i == k) break;
    }
    return recMoments[ind];
}

Double_t TIdentity2D::GetMean(Int_t i)
{
    if(i >= TSize)
    {
        cout<<"TIdentity2D::GetMean.Info: out of bound"<<endl;
        return -1000.;
    }
    return aver[i];
}

Double_t TIdentity2D::GetMeanI(Int_t i)
{
    if(i >= TSize)
    {
        cout<<"TIdentity2D::GetMeanI.Info: out of bound"<<endl;
        return -1000.;
    }
    return averI[i];
}

Double_t TIdentity2D::GetNuDyn(Int_t i, Int_t j)
{
    Double_t nydyn = GetSecondMoment(i)/GetMean(i)/GetMean(i);
    nydyn += GetSecondMoment(j)/GetMean(j)/GetMean(j);
    nydyn -= 2.*GetMixedMoment(i,j)/GetMean(i)/GetMean(j);
    nydyn -= (1./GetMean(i) + 1./GetMean(j));
    return nydyn;
}

Int_t TIdentity2D::GetIndex(Int_t k, Int_t &a, Int_t &b)
{
    Int_t t = 0;
    for(Int_t m = 0; m < TSize-1; m++)
        for(Int_t n = m+1; n < TSize; n++)
        {
            if(t == k) { a = m; b = n; return 1;}
            t++;
        }
        return 1;
}
