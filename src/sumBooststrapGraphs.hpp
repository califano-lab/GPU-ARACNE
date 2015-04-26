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
__host__
void matTo3cols(T *graphMiMat, *graphCountMat, unsigned int nRows, unsigned int nCols, T cutValue )
{
    T tfidx[];
    T geneidx[];
    T mi[];
    T count[];
    for (int i = 0; i < nRows; ++i )
    {
      for (int j = 0; i< nCols; ++j )
      {
	id = i * nCols + j;
	mi = graphMat[ i * nCols + j ];
	if( mi > cutValue )
	{
	  outGraph[i * 4 ] = {i ,j , mi, graphCountMat [id] };
	} 
      }
    }
    // 
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
    bool *d_graph3cols;
    HANDLE_ERROR(cudaMalloc((void **)&d_graph3cols, sizeof(float) * nRows * nCols )); // TODO need to calculate the size 
    // set memory to -1 initializes all values to true
    dim3 bDim(8, 8, 16);
    dim3 gDim(ceil(nRows /(1.0 * bDim.x )), ceil(nRows /(1.0 * bDim.y)), ceil(nCols /(1.0 * bDim.z) ));
    float cutValue = 0.0;
    matTo3cols_d <<< gDim, bDim >>>(d_rawGraph, nRows, nCols, cutValue );

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
