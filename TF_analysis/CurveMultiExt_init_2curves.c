#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *InData, *tic, *init_double; // TFR, tfrtic*fs, initial point
    double *multi;  // which multiples you want to extract
    double *lambda, *mu;    // smoothness and structure penalties
    double *band;   // The searching band for the fundamental.
    double *band_multi; // The searching bandwidth for the multiples.
    int NRow, NCol, Ntic, Ndim;
    double INIT;    // W.r.t. an initial point or not

    int *Freq0, *Freq1, *fundM, *init;
    double *FreqOut1, *FreqOut0;
    double *Energy, *FVal;
    double minval, val;

    int bw2, bw1;
    int i,j,k,l,j1;

    /* Portal to matlab */
    // TFR matrix: Row is frequency and Column is time
    InData = mxGetPr(prhs[0]);
    NRow = mxGetM(prhs[0]);
    NCol = mxGetN(prhs[0]);

    // tfrsqtic*fs sequence: Give the true freuqency information
    tic = mxGetPr(prhs[1]);
    Ntic = mxGetM(prhs[1]);
    // freuqency resolution
    double fr = tic[1] - tic[0];

    lambda = mxGetPr(prhs[2]);
    double l1 = lambda[0], l2 = lambda[1];

    mu = mxGetPr(prhs[3]);
    double mu1 = mu[0]/fr/fr;

    multi = mxGetPr(prhs[4]);
    double multi1 = multi[0];

    band = mxGetPr(prhs[5]);
    if (band[1] <= band[0]) {
        printf("Invalid searching band!\n");
        return;
    }

    band_multi = mxGetPr(prhs[6]);
    // Check that the multiple's band is within the range of frequecies
    if (band[1]*multi1 + band_multi[0] >= tic[Ntic-1]) {
        band[1] = (tic[Ntic-1] - band_multi[0])/multi1;
        if (band[0] >= band[1]) {
            printf("Out of range!\n");
            return;
        }
    }

    int lowIdx = 0, highIdx = Ntic-1;
    // The tic index corresponding to the true fundamental band
    // They're set to be: lowIdx >= 1 and highIdx <= Ntic-2
    for (k = Ntic-1; k >= 1; --k) {
        // search for the true indices for multiple frequency in tfrtic*fs
        if (tic[k]>=band[0] && tic[k-1]<=band[0])
            lowIdx = k;     // lowIdx cannot be the minimum tic
        if (tic[k]>=band[1] && tic[k-1]<=band[1])
            highIdx = k-1;  // highIdx cannot be the maximum tic
        if ((lowIdx != 0) && (highIdx != Ntic-1))
            break;
    }
    bw1 = highIdx - lowIdx + 1;
    if (bw1<1 || bw1==Ntic) {
        printf("Searching band is too narrow!\n");
        return;
    }
    bw2 = floor(band_multi[0]/fr);

    INIT = mxGetScalar(prhs[7]);

    // The starting frequency point
    init_double = mxGetPr(prhs[8]);
    Ndim = mxGetM(prhs[8]);
    init = (int *)mxMalloc(sizeof(int)*(Ndim));
    for (i = 0; i < Ndim; ++i) {
        init[i] = (int) init_double[i];
        init[i] -= 1;   //Matlab starts from 1
    }
    
    plhs[0] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);

    FreqOut0 = mxGetPr(plhs[0]);
    FreqOut1 = mxGetPr(plhs[1]);

    /* Main operations start here */
    const double eps = 1e-8;
    double sum = 0, tmp;

    printf("tfr_size = %dx%d, band_size = %dx%d, ", NRow, NCol, bw1, bw2);
    printf("Multiple = %.0f\n", multi1);

    Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
    FVal = (double *)mxMalloc(sizeof(double)*NRow*(bw1)*(2*bw2));
    Freq0 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq1 = (int *)mxMalloc(sizeof(int)*NRow);

    // fundM: the tic sequence corresponding to the multiple of the fundamental tic sequence
    fundM = (int *)mxMalloc(sizeof(int)*(2*bw1));

    double freqVal1;
    int tic1;
    for (j1 = 0; j1 < bw1; j1++) {
        freqVal1 = tic[lowIdx+j1]*multi1;
        tic1 = Ntic;
        for (k = lowIdx+j1; k < Ntic-1; k++) {
            if (tic[k]<freqVal1 && tic[k+1]>=freqVal1)
                tic1 = k+1;
        }
        if (tic1==Ntic) { // If we can't find the tic index, end it.
            printf("tic Out of range!\n");
            for (i = 0; i < NRow; ++i) {
                FreqOut0[i] = 1.;
                FreqOut1[i] = 1.;
            }
            return;
        } else {
            fundM[j1*2] = (tic1-bw2 > 0) ? (tic1-bw2) : 0;
            fundM[j1*2+1] = (tic1+bw2 < NCol-1) ? (tic1+bw2) : NCol-1;
        }
    }

    // TFR matrix
    for (i = 0; i < NRow; ++i) {
        for (j = 0; j < NCol; ++j) {
            Energy[i*NCol+j] = InData[i+j*NRow];
            sum += Energy[i*NCol+j];
        }
    }
    // If the energy distribution of the TFR is zero, end the program.
    if (sum < eps) { //if sum == 0
        printf("zero-energy TFR!\n");
        for (i=0; i<NRow; ++i) {
            FreqOut0[i] = 1.;
            FreqOut1[i] = 1.;
        }
        return;
    }

    for (i = 0; i < NRow; ++i) {
        for (j = 0; j < NCol; ++j)
            Energy[i*NCol+j] = -log(Energy[i*NCol+j]/sum + eps);
    }

    // Pre-allocation (i=0)
    for (j1 = 0; j1 < bw1; ++j1) {
        for (int j2 = fundM[j1*2]; j2 < fundM[j1*2+1]; ++j2) {
            FVal[j1*(2*bw2) + (j2-fundM[j1*2])]\
            = Energy[j1+lowIdx] + Energy[j2]\
            + INIT * (l1 * (j1 + lowIdx - init[0]) * (j1 + lowIdx - init[0])\
            + l2 * (j2 - init[1]) * (j2 - init[1])\
            + mu1 * (multi1*tic[j1+lowIdx] - tic[j2]) * (multi1*tic[j1+lowIdx] - tic[j2]));
            // INIT: Indicate if there is an initial point of the curves
            // lambdas: smoothness penalties
            // mus: similarity penalties'
        }
    }

    // Matrix that records the score of each positions
    for (i = 1; i < NRow; ++i) {

    for (j1 = 0; j1 < bw1; ++j1) {
        for (int j2 = fundM[j1*2]; j2 < fundM[j1*2+1]; ++j2) {
            minval = 1e16;
            for (int k1 = 0; k1 < bw1; ++k1) {
                for (int k2 = fundM[k1*2]; k2 < fundM[k1*2+1]; ++k2) {
                    tmp = FVal[(i-1)*(bw1*2*bw2)+k1*(2*bw2)+(k2-fundM[k1*2])] \
                    + l1 * (k1-j1) * (k1-j1) + l2 * (k2-j2) * (k2-j2) \
                    + mu1 * (multi1*tic[k1+lowIdx]-tic[k2]) * (multi1*tic[k1+lowIdx]-tic[k2]);
                    // lambdas: smoothness penalties
                    // mus: similarity penalties
                    if(tmp < minval) // Records the optimal answer at this position
                        minval = tmp;
                }
            }
            FVal[i*(bw1*2*bw2)+j1*(2*bw2)+(j2-fundM[j1*2])] = \
            minval + Energy[i*NCol+(j1+lowIdx)] + Energy[i*NCol+j2];
        }
    }

    }

    // Retrieve the answer back
    minval = FVal[(NRow-1)*(bw1*2*bw2)];
    Freq0[NRow-1] = 0;
    Freq1[NRow-1] = 0;
    for (j1 = 0; j1 < bw1; ++j1) {
        for (int j2 = fundM[j1*2]; j2 < fundM[j1*2+1]; ++j2) {
            tmp = FVal[(NRow-1)*(bw1*2*bw2) + j1*(2*bw2) + (j2-fundM[j1*2])];
            if (tmp < minval) {
                minval = tmp;
                Freq0[NRow-1] = j1;
                Freq1[NRow-1] = j2;
            }
        }
    }

    for (i = NRow-2; i >= 0; --i) {
        val = FVal[(i+1)*(bw1*2*bw2) + Freq0[i+1]*(2*bw2) + (Freq1[i+1]-fundM[Freq0[i+1]*2])] \
        - Energy[(i+1)*NCol+(Freq0[i+1]+lowIdx)] - Energy[(i+1)*NCol+Freq1[i+1]];
        
        for (j1 = 0; j1 < bw1; ++j1) {
            for (int j2 = fundM[j1*2]; j2 < fundM[j1*2+1]; ++j2) {
                tmp = FVal[i*(bw1*2*bw2) + j1*(2*bw2) + (j2-fundM[j1*2])]\
                + l1 * (j1 - Freq0[i+1]) * (j1 - Freq0[i+1]) \
                + l2 * (j2 - Freq1[i+1]) * (j2 - Freq1[i+1]) \
                + mu1 * (multi1*tic[j1+lowIdx] - tic[j2]) * (multi1*tic[j1+lowIdx] - tic[j2]);
                if(fabs(val-tmp) < eps) {
                    // count = 1;
                    Freq0[i] = j1;
                    Freq1[i] = j2;
                    // get out of the loops
                    // j2 += 2*bw2;
                    // j1 = bw1;
                }
            }
        }
    }

    for (i=0; i<NRow; ++i) {
        FreqOut0[i] = (double)Freq0[i]+lowIdx+1.;
        FreqOut1[i] = (double)Freq1[i]+1.;
    }

    mxFree(Freq0);
    mxFree(Freq1);
    mxFree(fundM);
    mxFree(Energy);
    mxFree(FVal); 
    return;
}