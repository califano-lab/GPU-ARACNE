#ifndef miAP_HPP
#define miAP_HPP
#include "Cube.hpp"
#include "util.hpp"
#include <cstdio>
#include <cstdlib>

template <typename T>
__global__ 
void computeMi(T *d_rankMatrix, unsigned int nTFs, unsigned int nGenes, unsigned int nSamples, 
        unsigned int *d_TFGenesIdx, T **d_rawGraph)
{
    unsigned int TFIdx = blockIdx.x;
    unsigned int geneIdx = blockIdx.y;
    // need pick two rows of the d_rankMatrix
    // first row: d_TFGenesIdx[TFIdx]
    // second row: geneIdx
    __shared__ Cube *cubeArray = new Cube[nSamples / 4 + 1];
    unsigned int head = 0;
    unsigned int tail = 0;
    Cube cube;
    cube.totalCount = nSamples;
    cube.upper = ......................
    dim3 blockDim(256, 1, 1);
    dim3 gridDim(ceil(nSamples / 256.0), 1, 1);
    // launch another kernel here


}





template <typename T>
void miAP(T *d_rankMatrix, unsigned int nTFs, unsigned int nGenes, unsigned int nSamples, 
        unsigned int *d_TFGeneIdx, T **d_rawGraph)
{
    dim3 gridDim(nTFs, nGenes, 1);
    dim3 blockDim(1, 1, 1);

    HANDLE_ERROR( cudaMalloc((void **)d_rawGraph, sizeof(T) * nTFs * nGenes) );
    computeMi<<<blockDim, gridDim>>>(d_rankMatrix, nTFs, nGenes, nSamples, d_TFGeneIdx, d_mi);
    HANDLE_ERROR( cudaGetLastError() );
    HANDLE_ERROR( cudaDeviceSynchronize() );
}
    
#endif 
