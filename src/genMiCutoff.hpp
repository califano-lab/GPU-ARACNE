#ifndef GENMICUTOFF_HPP
#define GENMICUTOFF_HPP

#include <cuda.h>
#include "util.hpp"
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <cmath>
#define PAIRS 100000
#define OUTLIERS 100
#define FITTING 10000

/*
 * this constructs the null model of the aracne computation
 */


// generate random matrix
// top half for the first element of a pair
// lower half for the second element of a pair

float computeMiThreshold(unsigned int nSamples, float pValue, unsigned int seed)
{
    // generate a random permuted matrix
    Matrix<int> *h_randomMatrix = new Matrix<int>(PAIRS * 2, nSamples);
    h_randomMatrix->permute(seed);

    // load the random permuted matrix onto GPU
    unsigned int *d_randomMatrix;
    HANDLE_ERROR (cudaMalloc((void **)&d_randomMatrix, h_randomMatrix->size()));
    HANDLE_ERROR (cudaMemcpy((void *)d_randomMatrix, (void *)h_randomMatrix->memAddr(), 
                h_randomMatrix->size(), cudaMemcpyHostToDevice));
    
    // declare the result variable and 
    // call miAP() to compute all the mutual information values in the d_miResult
    float *d_miResult;
    miAP(d_randomMatrix, PAIRS, nSamples, &d_miResult);
    
    // on device sort of this 100000 element array
    // first need to wrap the device pointers
    thrust::device_ptr<float> w_miResult(d_miResult);
    thrust::stable_sort(w_miResult, w_miResult + PAIRS);
    
    // copy data back to host
    float *h_miResult = new float[PAIRS];
    HANDLE_ERROR (cudaMemcpy((void *)h_miResult, (void *)d_miResult, 
                PAIRS * sizeof (float), cudaMemcpyDeviceToHost) );
    HANDLE_ERROR (cudaDeviceSynchronize());    
    // copy (PAIRS - OUTLIERS - FITTING) to (PAIRS - OUTLIERS)
    int firstPointToFit = PAIRS - OUTLIERS - FITTING;
    float *X_logPValues = new float[FITTING];
    float *Y_miValues = new float[FITTING];
    for (int i = 0; i < FITTING; i++){
        X_logPValues[i] = log((float)(PAIRS - (firstPointToFit + i)) / (float)PAIRS);
        Y_miValues[i] = h_miResult[firstPointToFit + i];
    }

    // linear regression of X_logPValues with Y_miValues
    // let's hard code it
    float X_bar = 0;
    float Y_bar = 0;
    float X_square_bar = 0;
    float XY_bar = 0;
    for (int i = 0; i < FITTING; i++){
        X_bar += X_logPValues[i];
        X_square_bar += X_logPValues[i] * X_logPValues[i];
        Y_bar += Y_miValues[i];
        XY_bar += X_logPValues[i] * Y_miValues[i];
    }
    X_bar = X_bar / FITTING;
    Y_bar = Y_bar / FITTING;
    XY_bar = XY_bar / FITTING;
    X_square_bar = X_square_bar / FITTING;
    float beta = (XY_bar - X_bar * Y_bar) / (X_square_bar - X_bar * X_bar);
    float alpha = Y_bar - beta * X_bar;

    // return the mutual information threshold based on the pValue given
    HANDLE_ERROR(cudaFree(d_miResult));
    HANDLE_ERROR(cudaFree(d_randomMatrix));
    delete h_randomMatrix;
    delete[] h_miResult;
    delete[] X_logPValues;
    delete[] Y_miValues;
    return log(pValue) * beta + alpha;
}

#endif
