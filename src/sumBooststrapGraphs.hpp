#ifndef PRUNEGRAPH_H
#define PRUNEGRAPH_H
#include <cuda.h>
#include "util.hpp"

// the caller launches 2000 * 2000 * 20000 threads
// the kernel is implemented here
// first dimension: TF1
// second dimension: TF2
// third dimension: all genes
// decision array should be initialized to true


template<typename T>
__global__
void label(T *d_rawGraph, unsigned int nRows, unsigned int nCols, bool *d_decision, unsigned int *d_TFGeneIdx)
{
    unsigned int TF1 = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int TF2 = blockIdx.y * blockDim.y + threadIdx.y;
    if (TF1 >= nRows || TF2 >= nRows) return;
    T val0 = d_rawGraph[TF1 * nCols + d_TFGeneIdx[TF2]];   
    if  (val0 < 0) return;
    unsigned int geneIdx = blockIdx.z * blockDim.z + threadIdx.z;
    T val1 = d_rawGraph[TF1 * nCols + geneIdx];
    T val2 = d_rawGraph[TF2 * nCols + geneIdx];  
    if (val1 < 0 || val2 < 0) return;
    
    if (val0 < val1 && val0 < val2) d_decision[TF1 * nCols + d_TFGeneIdx[TF2]] = false;
    if (val1 < val0 && val1 < val2) d_decision[TF1 * nCols + geneIdx] = false;
    if (val2 < val0 && val2 < val1) d_decision[TF2 * nCols + geneIdx] = false;
}

template<typename T>
__global__
void sumBooststrapGraphs(T *d_sumMiGraph, *d_sumCountGraph, *d_currGraph, unsigned int nRows, unsigned int nCols, bool *d_decision)
{
    unsigned int rowIdx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int colIdx = blockIdx.y * blockDim.y + threadIdx.y;
    if (rowIdx >= nRows || colIdx >= nCols) return;
    if (d_decision[rowIdx * nCols + colIdx])
    { 
        d_sumCountGraph[rowIdx * nCols + colIdx] += 1;
        d_sumMiGraph[rowIdx * nCols + colIdx] += d_currGraph[rowIdx * nCols + colIdx];
    }
}

template<typename T>
__host__
void graphMatrixTo3cols(T *d_rawGraph, unsigned int nRows, unsigned int nCols, unsigned int *d_TFGeneIdx)
{
    bool *d_;
    HANDLE_ERROR(cudaMalloc((void **)&d_decision, sizeof(bool) * nRows * nCols));
    // set memory to -1 initializes all values to true
    HANDLE_ERROR(cudaMemset((void *)d_decision, -1, sizeof(bool) * nRows * nCols));
    // sequentially launche two kernels
    dim3 blockDimLabel(8, 8, 16);
    dim3 gridDimLabel(ceil(nRows / 8.0), ceil(nRows / 8.0), ceil(nCols / 16.0));
    label<<<gridDimLabel, blockDimLabel>>>(d_rawGraph, nRows, nCols, d_decision, d_TFGeneIdx);
    HANDLE_ERROR(cudaGetLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());
    dim3 blockDimCut(16, 64, 1);
    dim3 gridDimCut(ceil(nRows / 16.0), ceil(nCols / 64.0), 1);
    cut<<<gridDimCut, blockDimCut>>>(d_rawGraph, nRows, nCols, d_decision);
    HANDLE_ERROR(cudaGetLastError());
    HANDLE_ERROR(cudaDeviceSynchronize());
    cudaFree(d_decision);
}
#endif
