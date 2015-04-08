#ifndef PRUNEGRAPH_H
#define PRUNEGRAPH_H
#include <cuda.h>

// the caller launches 2000 * 2000 * 20000 threads
// the kernel is implemented here
// first dimension: TF1
// second dimension: TF2
// third dimension: all genes
// decision array should be initialized to true


//           val0
//TF1 <--------------> TF2
//     \            /
//      \          /
//       \        /
//  val1  \      /  val2
//         \    /
//          \  /
//          Gene

template<typename T>
__global__
void label(T *miData, unsigned int nRows, unsigned int nCols, bool *decision)
{
    unsigned int TF1 = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int TF2 = blockIdx.y * blockDim.y + threadIdx.y;
    T val0 = miData[TF1 * nRows + TF2];   
    if  (val0 < 0) return;
    unsigned int geneIdx = blockIdx.z * blockDim.z + threadIdx.z;
    T val1 = miData[TF1 * nRows + geneIdx];
    T val2 = miData[TF2 * nRows + geneIdx]  
    if (val1 < 0 || val2 < 0) return;
    
    if (val0 < val1 && val0 < val2) decision[TF1 * nRows + TF2] = false;
    if (val1 < val0 && val1 < val2) decision[TF1 * nRows + geneIdx] = false;
    if (val2 < val0 && val2 < val1) decision[TF2 * nRows + geneIdx] = false;

    // handle symmetry issues
}

template<typename T>
__global__
void cut(T *miData, unsigned int nRows, unsigned int nCols, bool *decision){}

template<typename T>
__host__
void pruneGraph(T *rawGraph, unsigned int nRows, unsigned int nCols, unsigned int *d_TFGeneIdx)
{
    bool *decsion = new bool[nTFs * nGenes];
    // initialization here
    dim3 blockDim(8, 8, 16);
    dim3 gridDim(ceil(nRows / 8.0), ceil(nRows / 8.0), ceil(nCols / 16.0));
    label<<<gridDim, blockDim>>>(rawGraph, nRows, nCols, decision);
    dim3 blockDim(16, 64, 1);
    dim3 gridDim(ceil(nRows / 16.0), ceil(nCols / 64.0), 1);
    cut<<<gridDim, blockDim>>>(rawGraph, nRows, nCols, decision);
    delete decision;
}
#endif
