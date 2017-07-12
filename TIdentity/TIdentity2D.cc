//Author A. Rustamov
//TIdentity class for calculation of all second moments
#include "TIdentity2D.h"
#include "TIdentityFunctions.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TF1.h"

using namespace std;
ClassImp(TIdentity2D)

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
 
  countPart    = 0;
  countPartNeg = 0;
  countPartPos = 0;
    
  useSign = -1000;
  sign = -1000;
  countVeto = 0;  
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
        cout<<"resetting vectors"<<endl;
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
    
    
   // useSign = -1000;
   // sign = -1000;
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

void TIdentity2D::GetTree(Long_t &nent, TString MC)
{
  TIdentityFile = new TFile(inputDir+fileName);
  cout<<"we are reading the file "<<inputDir+fileName<<endl;
  //TIdentityTree = (TTree*)TIdentityFile->Get("fEventTree");
  if(MC == "sim")
  {
      cout<<"running for sim"<<endl;
     TIdentityTree = (TTree*)TIdentityFile->Get("fIdenTreeMC");
  }
  else
  {
      cout<<"running for data"<<endl;
     TIdentityTree = (TTree*)TIdentityFile->Get("fIdenTree");
  }
  TIdentityTree -> SetBranchAddress("sign",&sign);
  TIdentityTree -> SetBranchAddress("myBin",myBin);
  TIdentityTree -> SetBranchAddress("myDeDx",&myDeDx);
  TIdentityTree -> SetBranchAddress("evtNum",&evtNum);
  TIdentityTree -> SetBranchAddress("px",&momX);
  TIdentityTree -> SetBranchAddress("py",&momY);
  TIdentityTree -> SetBranchAddress("pz",&momZ);
  nEntries      = (Long_t)TIdentityTree -> GetEntries();
  //histoBin  = new TH1D("histoBin","histoBin",150,ffMin,ffMax);
  //debugFile = new TFile(outputDir+"debugFile.root","recreate");
  nent = nEntries;  
  TIdentityFile -> cd();
  InitFunctions();
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
  Double_t val[4];
  Double_t sumVal = 0;
  for(Int_t ii = 0; ii < 4; ii++)
  {
    val[ii] = functions -> GetValue(ii, xx[0]);
    sumVal += val[ii]; 
  }
  Double_t relVal[4];
  if(sumVal < 1e-15) return 0.;
  for(Int_t m = 0; m < 4; m++)
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
  Double_t val[4];
  Double_t sumVal = 0;
  for(Int_t ii = 0; ii < 4; ii++)
  { 
    val[ii] = functions -> GetValue(ii, xx[0]);
    sumVal += val[ii]; 
  }
  if(sumVal < 1e-15) return 0.;
  Double_t myVal[6];
  Int_t t = 0;
  for(Int_t m = 0; m < 3; m++)
    for(Int_t n = m+1; n < 4; n++)
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
  if(i%5000000 == 0) cout<<"event "<<i<<" of "<<nEntries <<endl;
  TIdentityTree -> GetEntry(i);
    
  /////
    
    Float_t mom = sqrt(momX*momX + momY*momY + momZ*momZ);
    Float_t pt  = sqrt(momX*momX + momY*momY);
    
    if( mom < minMom || mom > maxMom ) return kFALSE;
    if( pt  < minPt  || pt > maxPt )   return kFALSE;
    
  /////
    
    
    
  if( evtNum != prevEvtVeto && prevEvtVeto >0) {countVeto++;}   
  prevEvtVeto = evtNum;
  if ( (myDeDx < ffMin || myDeDx > ffMax) && myDeDx > 0 ) return kFALSE;   
  if(sign != useSign && useSign != 0) return kFALSE;
  return kTRUE;
}

Int_t TIdentity2D::AddEntry(Bool_t &isAdd)
{
  //if(myBin[0] == 1 && myBin[1] == 10 && myBin[2] == 6 && myBin[3] == 3)
    //{
      //histoBin -> Fill(myDeDx);
    //}   
    isAdd = kTRUE;
  if(evtNum == prevEvt) 
    {
        //cout<<"adding1  "<<myDeDx<<"  "<< endl;
      AddParticles();
        //cout<<"added1" << endl;
        
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
        //cout<<"adding"<<endl;
      AddParticles();
        //cout<<"added"<<endl;
        
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
  Double_t mValue[4] = {0.};
  Double_t sumValue = 0;

    //cout<<"before eval "<<myDeDx<<"  "<<myBin[0]<<"  "<<myBin[1]<<"  "<<myBin[2]<<endl;
    
  for(Int_t i = 0; i < TSize; i++)
     {	
	  
      if( myDeDx < 0 ) { mValue[i] = 0; count++; continue; }
      mValue[i] = TFunctions[i] -> Eval(myDeDx);
         //cout<<"i == "<<i <<" "<< endl;
      sumValue += mValue[i];
     }
    
    //cout<<"after eval"<<endl;
    
    wProton = wKaon = wPion = -100.;
    if(sumValue > 1e-10)
    {
     wProton = mValue[2]/sumValue;
     wKaon   = mValue[3]/sumValue;
     wPion   = mValue[1]/sumValue;
    }
    
    countPart += 1;
    if(sign == -1)
    {
        countPartNeg += 1;
    }
    if(sign == 1)
    {
        countPartPos += 1;
    }
    
  //if(sumValue != 0)
    if(sumValue > 1e-10)
    {
      count++;
      for(Int_t i = 0; i < TSize; i++)
	{
	  
        W_sum[i] += mValue[i]/sumValue; //this should be
        
        //cout<<"testing "<< mValue[i]/sumValue <<endl;
      ////// to test
        //if( i == 0 ) W_sum[i] += mValue[i]/sumValue;
        
       
        
        
        
      //////
        
	}
    }
}

void TIdentity2D::Finalize()
{
  cout<<" "<<endl;
  cout<<"************number of analyzed events: "<<countVeto<<"******"<<endl;
  cout<<"***************************************************"<<endl;
  cout<<"***************************************************"<<endl;
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
    
  //CalcMoments(); // Indi deyishdim ki, ayrica 2 - ci momentler hesablansin
    
    
    
 // GetMoments();  bu lazim deyildi onsuz da;
 
 
    
    /*
  ofstream outf("output.txt");
  for(Int_t i = 0; i < TSize; i++)
    {
      outf<<GetMean(i)<<endl;      
    }
  
  for(Int_t i = 0; i < TSize; i++)
    {
      outf<<GetSecondMoment(i)<<endl;
    }
  
  for(Int_t m = 0; m < TSize-1; m++)
    for(Int_t n = m+1; n < TSize; n++)
      {
	outf<<GetMixedMoment(m,n)<<endl;
      }
  */
  //debugFile -> cd();
  //histoBin  -> Write();
  //debugFile -> Close();
}

void TIdentity2D::GetBins(Int_t *bins)
{
  bins[0] = myBin[0];
  bins[1] = myBin[1];
  bins[2] = myBin[2];
  //bins[3] = myBin[3];
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
      B[kk+4] = averMixed[kk];
      GetIndex(kk, indA, indB);
      for(Int_t m = 0; m < TSize; m++)
	{
	  B[kk+4] -= aver[m]*(GetWI(kk,m,0)-GetWI(indA,m,1)*GetWI(indB,m,1));
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
       cout<<"out of bound"<<endl;
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
      cout<<"out of bound"<<endl;
      return -1000.;
    }
  return aver[i];
}

Double_t TIdentity2D::GetMeanI(Int_t i)
{
  if(i >= TSize) 
    {
      cout<<"out of bound"<<endl;
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