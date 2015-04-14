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
    // first row (define it X): d_TFGenesIdx[TFIdx]
    // second row (define it Y): geneIdx
    __shared__ Cube *cubeArray = new Cube[nSamples / 4 + 1];
    unsigned int head = 0;
    unsigned int tail = 1;
    Cube cube;
    cube.totalCount = nSamples;
    cube.upper = round(nSamples - 1);
    cube.lower = 0;
    cube.left = 0;
    cube.right = round(nSamples - 1);
    cubeArray[head] = cube;
    __shared__ unsigned int upperRight = 0;
    __shared__ unsigned int upperLeft = 0;
    __shared__ unsigned int lowerLeft = 0 ; 
    __shared__ unsigned int lowerRight = 0;
    do{
        unsigned int colIdx = threadIdx.x;
        unsigned int middleX = (cubeArray[head].left + cubeArray[head].right) / 2;
        unsigned int middleY = (cubeArray[head].upper + cubeArray[head].lower) / 2;
        unsigned int coordX = round(d_rankMatrix[d_TFGenesIdx[TFIdx] * nSamples + colIdx]);
        unsigned int coordY = round(d_rankMatrix[geneIdx * nSamples + colIdx]);
        
        if (cubeArray[head].left <= coordX && coordX < middleX){
            if (cubeArray[head].lower <= coordY && coordY < middleY){
                atomicAdd(&lowerLeft, 1);
            } else {
                atomicAdd(&upperLeft, 1);
            }
        } else {
            if (cubeArray[head].lower <= coordY && coordY < middleY){
                atomicAdd(&lowerRight, 1);
            } else {
                atomicAdd(&upperRight, 1);
            }
        }
        __syncthreads();

        // designate one thread to compute this how to optimize?
        if (threadIdx.x == 0){
            if (chiSeq(upperRight, upperLeft, lowerLeft, lowerRight)) {// significant
                cubeArray[tail].totalCount = upperRight;
                cubeArray[tail].upper 
                    ........ ...........
                cubeArray[tail + 1].totalCount = upperLeft;
                cubeArray[tail + 2].totalCount = lowerLeft;
                cubeArray[tail + 3].totalCount = lowerRight;
                tail = tail + 4;
                head = head + 1;            
            }
        }
    } while(head < tail);
}





template <typename T>
void miAP(T *d_rankMatrix, unsigned int nTFs, unsigned int nGenes, unsigned int nSamples, 
        unsigned int *d_TFGeneIdx, T **d_rawGraph)
{
    dim3 gridDim(nTFs, nGenes, 1);
    dim3 blockDim(nSamples, 1, 1);

    HANDLE_ERROR( cudaMalloc((void **)d_rawGraph, sizeof(T) * nTFs * nGenes) );
    computeMi<<<blockDim, gridDim>>>(d_rankMatrix, nTFs, nGenes, nSamples, d_TFGeneIdx, d_mi);
    HANDLE_ERROR( cudaGetLastError() );
    HANDLE_ERROR( cudaDeviceSynchronize() );
}
    
#endif 
