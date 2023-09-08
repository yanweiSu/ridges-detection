#include <stdlib.h>
#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  double *InData;
  int NRow,NCol;
  double lambda;

  int *Freq;
  double *FreqOut;
  double *FValOut;
  double *Energy, *FVal;
  double minval, val;
 
  int i,j,k;

  /* Portal to matlab */
  InData = mxGetPr(prhs[0]);
  NRow = mxGetM(prhs[0]);
  NCol = mxGetN(prhs[0]);

  lambda = mxGetScalar(prhs[1]);

  plhs[0] = mxCreateNumericArray(1,&NRow,mxDOUBLE_CLASS,mxREAL);
  FreqOut = mxGetPr(plhs[0]);
  

  /* Main operations start here */
  const double eps = 1e-8;
  double sum = 0;

  Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
  FVal = (double *)mxMalloc(sizeof(double)*NRow*NCol);
  Freq = (int *)mxMalloc(sizeof(int)*NRow);
  for (i=0;i<NRow;++i)
    for (j=0;j<NCol;++j){
      Energy[i*NCol+j] = InData[i+j*NRow];
      sum += Energy[i*NCol+j];
    }

  for (i=0;i<NRow;++i) {
    for (j=0;j<NCol;++j)
      Energy[i*NCol+j] = -log(Energy[i*NCol+j]/sum+eps);
  }

  for (j=0;j<NCol;++j)
    FVal[j] = Energy[j];

  for (i=1;i<NRow;++i) {
    for (j=0;j<NCol;++j)
      FVal[i*NCol+j] = 1e16;
    for (j=0;j<NCol;++j) {
      for (k=0;k<NCol;++k) {
        if (FVal[i*NCol+j] > FVal[(i-1)*NCol+k]+lambda*(k-j)*(k-j))
          FVal[i*NCol+j] = FVal[(i-1)*NCol+k]+lambda*(k-j)*(k-j);
      }
      FVal[i*NCol+j] += Energy[i*NCol+j];
    }
  }

  minval = FVal[(NRow-1)*NCol];
  Freq[NRow-1] = 0;
  for (j=1;j<NCol;++j) {
    if (FVal[(NRow-1)*NCol+j]<minval) {
      minval = FVal[(NRow-1)*NCol+j];
      Freq[NRow-1] = j;
    }
  }

  for (i=NRow-2;i>=0;--i){
    val = FVal[(i+1)*NCol+Freq[i+1]] - Energy[(i+1)*NCol+Freq[i+1]];
    for (j=0;j<NCol;++j)
      if (fabs(val-FVal[i*NCol+j]-lambda*(j-Freq[i+1])*(j-Freq[i+1]))<eps){
        Freq[i] = j;
        break;
      }
  }
          
  for (i=0;i<NRow;++i){
    FreqOut[i] = (double)Freq[i]+1.;
  }

  mxFree(Freq);
  mxFree(Energy);
  mxFree(FVal); 
  return;
}

 
