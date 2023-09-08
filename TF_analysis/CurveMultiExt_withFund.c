#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *InData, *tic, *fundIN;
    int NRow, NCol, Ntic, Nfund;
    double lambda, mu;

    int *Freq, *fundM, *fund;
    double *FreqOut;
    // double *FValOut;
    double *Energy, *FVal;
    double minval, val;

    double multi;
    int bw;
    int i,j,k,l;

    /* Portal to matlab */
    InData = mxGetPr(prhs[0]);
    NRow = mxGetM(prhs[0]);
    NCol = mxGetN(prhs[0]);

    tic = mxGetPr(prhs[1]);
    Ntic = mxGetM(prhs[1]);

    // Fundamental's tic
    fundIN = mxGetPr(prhs[2]);
    Nfund = mxGetM(prhs[2]);
    fund = (int *)mxMalloc(sizeof(int)*(NRow));
    for (i=0; i<NRow; ++i) {
        fund[i] = (int) fundIN[i];
        fund[i] -= 1;   //Matlab starts from 1
        if (fund[i]<0)
            fund[i] = 0;
    }

    lambda = mxGetScalar(prhs[3]);

    mu = mxGetScalar(prhs[4]);

    multi = mxGetScalar(prhs[5]);

    bw = mxGetScalar(prhs[6]);
    
    plhs[0] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);

    FreqOut = mxGetPr(plhs[0]);

    /* Main operations start here */
    if (Nfund != NRow) {
        printf("Error: Size not compatible");
        return;
    }
    const double eps = 1e-8;
    double sum = 0, tmp;

    Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
    FVal = (double *)mxMalloc(sizeof(double)*NRow*(2*bw));
    Freq = (int *)mxMalloc(sizeof(int)*NRow);
    fundM = (int *)mxMalloc(sizeof(int)*(NRow));

    // Build baseLine array
    double freqVal;
    for (i = 0; i<NRow; ++i) {
        freqVal = tic[fund[i]]*multi;
        fundM[i] = 0;
        for (k = Ntic-1; k >= 1; --k) {
            // search for the true indices for the multiple in tfrtic
            if (tic[k]>=freqVal && tic[k-1]<freqVal)
                fundM[i] = k;
        }
    }

    // TFR matrix
    for (i = 0; i < NRow; ++i) {
        for (j = 0; j < NCol; ++j) {
            Energy[i*NCol+j] = InData[i+j*NRow];
            sum += Energy[i*NCol+j];
        }
    }

    for (i = 0;i < NRow;++i) {
        for (j = 0;j < NCol;++j)
            Energy[i*NCol+j] = -log(Energy[i*NCol+j]/sum+eps);
    }

    // Pre-allocation
    // for (i=0; i<NRow; ++i)
    for (j = 0; j < 2*bw; ++j)
        FVal[j] = Energy[j+fundM[0]-bw];

    for (i = 1; i < NRow; ++i) {
        for (j = 0; j < 2*bw; ++j) {
            minval = 1e16;
            for (k = 0; k < 2*bw; ++k) {
                tmp = FVal[(i-1)*(2*bw) + k]\
                 + lambda * ((k+fundM[i-1])-(j+fundM[i])) * ((k+fundM[i-1])-(j+fundM[i]))\
                 + mu * (multi*tic[fund[i-1]]-tic[k+fundM[i-1]-bw]) * (multi*tic[fund[i-1]]-tic[k+fundM[i-1]-bw]);
                if(tmp < minval)
                    minval = tmp;
            }
            FVal[i*(2*bw) + j] = minval + Energy[i*NCol+(j+fundM[i]-bw)];
        }
    }

    minval = FVal[(NRow-1)*(2*bw)];
    Freq[NRow-1] = 0 + fundM[NRow-1] - bw;
    for (j = 0; j < 2*bw; ++j) {
        tmp = FVal[(NRow-1)*(2*bw) + j];
        if (tmp < minval) {
            minval = tmp;
            Freq[NRow-1] = j + fundM[NRow-1] - bw;
        }
    }

    for (i = NRow-2; i >= 0; --i) {
        val = FVal[(i+1)*(2*bw) + (Freq[i+1]-fundM[i+1]+bw)] - Energy[(i+1)*NCol+Freq[i+1]];
        for (j = 0; j < 2*bw; ++j) {
            tmp = FVal[i*(2*bw) + j]\
             + lambda * ((j+fundM[i]-bw)-Freq[i+1]) * ((j+fundM[i]-bw)-Freq[i+1])\
             + mu * (multi*tic[fund[i]]-tic[j+fundM[i]-bw]) * (multi*tic[fund[i]]-tic[j+fundM[i]-bw]);
            tmp = fabs(val-tmp);
            if(tmp < eps) {
                Freq[i] = j + fundM[i] - bw;
                j = 2*bw; // get out of loops
            }
        }
    }

    for (i=0; i<NRow; ++i) {
        FreqOut[i] = (double)Freq[i]+1.;
    }

    mxFree(Freq);
    mxFree(fundM);
    mxFree(Energy);
    mxFree(FVal); 
    return;
}