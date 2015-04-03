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
void DPI(T *miData, unsigned int nRows, unsigned int nCols, bool *decision)
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

#endif
