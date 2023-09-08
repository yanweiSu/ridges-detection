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

    int *Freq0, *Freq1, *Freq2, *fundM, *init;
    double *FreqOut1, *FreqOut2, *FreqOut0;
    double *Energy, *FVal;
    double minval, val;

    int bw2, bw3, bw1;
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
    double l1 = lambda[0], l2 = lambda[1], l3 = lambda[2];
    mu = mxGetPr(prhs[3]);
    double mu1 = mu[0]/fr/fr, mu2 = mu[1]/fr/fr;

    multi = mxGetPr(prhs[4]);
    double multi1 = multi[0], multi2 = multi[1];

    band = mxGetPr(prhs[5]);
    if (band[1] <= band[0]) {
        printf("Invalid searching band!\n");
        return;
    }

    band_multi = mxGetPr(prhs[6]);
    // Check that the multiple's band is within the range of frequecies
    if (band[1]*multi2 + band_multi[1] >= tic[Ntic-1]) {
        band[1] = (tic[Ntic-1] - band_multi[1])/multi2;
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
            lowIdx = k;
        if (tic[k]>=band[1] && tic[k-1]<=band[1])
            highIdx = k;
        if ((lowIdx != 0) && (highIdx != 0))
            break;
    }
    bw1 = highIdx - lowIdx + 1;
    if (bw1<1 || bw1==Ntic) {
        printf("Searching band is too narrow!\n");
        return;
    }
    bw2 = floor(band_multi[0]/fr);
    bw3 = floor(band_multi[1]/fr);


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
    plhs[2] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);

    FreqOut1 = mxGetPr(plhs[1]);
    FreqOut2 = mxGetPr(plhs[2]);
    FreqOut0 = mxGetPr(plhs[0]);

    /* Main operations start here */
    const double eps = 1e-8;
    double sum = 0, tmp;

    printf("tfr_size = %dx%d, band_size = %dx%dx%d, ", NRow, NCol, bw1, bw2, bw3);
    printf("Multiple = %.0f, %.0f\n", multi1, multi2);

    Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
    FVal = (double *)mxMalloc(sizeof(double)*NRow*(bw1)*(2*bw2)*(2*bw3));
    Freq0 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq1 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq2 = (int *)mxMalloc(sizeof(int)*NRow);

    // fundM: the tic sequence corresponding to the multiple of the fundamental tic sequence
    fundM = (int *)mxMalloc(sizeof(int)*(4*bw1));

    double freqVal1, freqVal2;
    int tic1, tic2;
    for (j1 = 0; j1 < bw1; j1++) {
        freqVal1 = tic[lowIdx+j1]*multi1;
        freqVal2 = tic[lowIdx+j1]*multi2;

        tic1 = Ntic;
        tic2 = Ntic;
        for (k = lowIdx+j1; k < Ntic-1; k++) {
            if (tic[k]<freqVal1 && tic[k+1]>=freqVal1)
                tic1 = k+1;
            if (tic[k]<freqVal2 && tic[k+1]>=freqVal2)
                tic2 = k+1;
        }
        if (tic1==Ntic || tic2==Ntic) {
            printf("Out of range!\n");
            for (i=0; i<NRow; ++i) {
                FreqOut0[i] = 1.;
                FreqOut1[i] = 1.;
                FreqOut2[i] = 1.;
            }
            return;
        } else {
            fundM[j1*4] = (tic1-bw2 > 0) ? (tic1-bw2) : 0;
            fundM[j1*4+1] = (tic1+bw2 < NCol-1) ? (tic1+bw2) : NCol-1;
            fundM[j1*4+2] = (tic2-bw3 > 0) ? (tic2-bw3) : 0;
            fundM[j1*4+3] = (tic2+bw3 < NCol-1) ? (tic2+bw3) : NCol-1;
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
        printf("zero TFR!\n");
        for (i=0; i<NRow; ++i) {
            FreqOut0[i] = 1.;
            FreqOut1[i] = 1.;
            FreqOut2[i] = 1.;
        }
        return;
    }

    for (i = 0; i < NRow; ++i) {
        for (j = 0; j < NCol; ++j)
            Energy[i*NCol+j] = -log(Energy[i*NCol+j]/sum + eps);
    }

    // Pre-allocation (i=0)
    for (j1 = 0; j1 < bw1; ++j1) {
        for (int j2 = fundM[j1*4]; j2 < fundM[j1*4+1]; ++j2) {
            for (int j3 = fundM[j1*4+2]; j3 < fundM[j1*4+3]; ++j3) {
                FVal[j1*(2*bw2*2*bw3) + (j2-fundM[j1*4])*(2*bw3) + (j3-fundM[j1*4+2])]\
                 = Energy[j1+lowIdx] + Energy[j2] + Energy[j3]\
                 + INIT * (l1 * (j1 + lowIdx - init[0]) * (j1 + lowIdx - init[0])\
                 + l2 * (j2 - init[1]) * (j2 - init[1])\
                 + l3 * (j3 - init[2]) * (j3 - init[2])\
                 + mu1 * (multi1*tic[j1+lowIdx] - tic[j2]) * (multi1*tic[j1+lowIdx] - tic[j2])\
                 + mu2 * (multi2*tic[j1+lowIdx] - tic[j3]) * (multi2*tic[j1+lowIdx] - tic[j3]));
                // INIT: Indicate if there is an initial point of the curves
                // lambdas: smoothness penalties
                // mus: similarity penalties'
           }
        }
    }   


    // Matrix that records the score of each positions
    for (i = 1; i < NRow; ++i) {

    for (j1 = 0; j1 < bw1; ++j1) {
        for (int j2 = fundM[j1*4]; j2 < fundM[j1*4+1]; ++j2) {
            for (int j3 = fundM[j1*4+2]; j3 < fundM[j1*4+3]; ++j3) {
                minval = 1e16;
                for (int k1 = 0; k1 < bw1; ++k1) {
                    for (int k2 = fundM[k1*4]; k2 < fundM[k1*4+1]; ++k2) {
                        for (int k3 = fundM[k1*4+2]; k3 < fundM[k1*4+3]; ++k3) {
                            tmp = FVal[(i-1)*(bw1*2*bw2*2*bw3)+k1*(2*bw2*2*bw3)+(k2-fundM[k1*4])*(2*bw3)+(k3-fundM[k1*4+2])]\
                             + l1 * (k1-j1) * (k1-j1)\
                             + l2 * (k2-j2) * (k2-j2)\
                             + l3 * (k3-j3) * (k3-j3)\
                             + mu1 * (multi1*tic[k1+lowIdx]-tic[k2]) * (multi1*tic[k1+lowIdx]-tic[k2])\
                             + mu2 * (multi2*tic[k1+lowIdx]-tic[k3]) * (multi2*tic[k1+lowIdx]-tic[k3]);
                            // lambdas: smoothness penalties
                            // mus: similarity penalties
                            if(tmp < minval) // Records the optimal answer at this position
                                minval = tmp;
                        }
                    }
                }
                FVal[i*(bw1*2*bw2*2*bw3)+j1*(2*bw2*2*bw3)+(j2-fundM[j1*4])*(2*bw3)+(j3-fundM[j1*4+2])] = \
                minval + Energy[i*NCol+(j1+lowIdx)] + Energy[i*NCol+j2] + Energy[i*NCol+j3];
           }
        }
    }

    }

    // Retrieve the answer back
    minval = FVal[(NRow-1)*(bw1*2*bw2*2*bw3)];
    Freq0[NRow-1] = 0;
    Freq1[NRow-1] = 0;
    Freq2[NRow-1] = 0;

    for (j1 = 0; j1 < bw1; ++j1) {
        for (int j2 = fundM[j1*4]; j2 < fundM[j1*4+1]; ++j2) {
            for (int j3 = fundM[j1*4+2]; j3 < fundM[j1*4+3]; ++j3) {
                tmp = FVal[(NRow-1)*(bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + (j2-fundM[j1*4])*(2*bw3) + (j3-fundM[j1*4+2])];
                if (tmp < minval) {
                    minval = tmp;
                    Freq0[NRow-1] = j1;
                    Freq1[NRow-1] = j2;
                    Freq2[NRow-1] = j3;
                }
            }
        }
    }

    for (i = NRow-2; i >= 0; --i) {
        val = FVal[(i+1)*(bw1*2*bw2*2*bw3) + Freq0[i+1]*(2*bw2*2*bw3)\
         + (Freq1[i+1]-fundM[Freq0[i+1]*4])*(2*bw3) + (Freq2[i+1]-fundM[Freq0[i+1]*4+2])]\
         - Energy[(i+1)*NCol+(Freq0[i+1]+lowIdx)] - Energy[(i+1)*NCol+Freq1[i+1]] - Energy[(i+1)*NCol+Freq2[i+1]];
        
        for (j1 = 0; j1 < bw1; ++j1) {
            for (int j2 = fundM[j1*4]; j2 < fundM[j1*4+1]; ++j2) {
                for (int j3 = fundM[j1*4+2]; j3 < fundM[j1*4+3]; ++j3) {
                    tmp = FVal[i*(bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + (j2-fundM[j1*4])*(2*bw3) + (j3-fundM[j1*4+2])]\
                     + l1 * (j1 - Freq0[i+1]) * (j1 - Freq0[i+1])\
                     + l2 * (j2 - Freq1[i+1]) * (j2 - Freq1[i+1])\
                     + l3 * (j3 - Freq2[i+1]) * (j3 - Freq2[i+1])\
                     + mu1 * (multi1*tic[j1+lowIdx] - tic[j2]) * (multi1*tic[j1+lowIdx] - tic[j2])\
                     + mu2 * (multi2*tic[j1+lowIdx] - tic[j3]) * (multi2*tic[j1+lowIdx] - tic[j3]);
                    if(fabs(val-tmp) < eps) {
                        // count = 1;
                        Freq0[i] = j1;
                        Freq1[i] = j2;
                        Freq2[i] = j3;
                        // get out of the loops
                        // j3 += 2*bw3;
                        // j2 += 2*bw2;
                        // j1 = bw1;
                    }
                }
            }
        }
    }

    for (i=0; i<NRow; ++i) {
        FreqOut0[i] = (double)Freq0[i]+lowIdx+1.;
        FreqOut1[i] = (double)Freq1[i]+1.;
        FreqOut2[i] = (double)Freq2[i]+1.;
    }

    mxFree(Freq0);
    mxFree(Freq1);
    mxFree(Freq2);
    mxFree(fundM);
    mxFree(Energy);
    mxFree(FVal); 
    return;
}