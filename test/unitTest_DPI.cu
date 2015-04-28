// basic viable unit of execution.
#include "Matrix.hpp"
#include "pruneGraph.hpp"
#include <cuda.h>
#include <cstdlib>
#include <cstdio>
#include "util.hpp"

// argument list: 
// <TFfile> <datafile> <nTFs> <nGenes> <nSamples> <nBootstraps>
int main(int argc, char *argv[])
{
   
    // DPI to prune network
    Matrix<float> *h_miValue = new Matrix<float>(3,3);
    h_miValue->element(0,0) = 0.7;
    h_miValue->element(0,1) = 0.2;
    h_miValue->element(0,2) = 0.1;
    h_miValue->element(1,0) = 0.2;
    h_miValue->element(1,1) = 0.7;
    h_miValue->element(1,2) = 0.3;
    h_miValue->element(2,0) = 0.1;
    h_miValue->element(2,1) = 0.3;
    h_miValue->element(2,2) = 0.7;
    unsigned int nTFs = 3;
    unsigned int nGenes = 3;
    unsigned int *h_TFGeneIdx = new unsigned int[3];
    h_TFGeneIdx[0] = 0;
    h_TFGeneIdx[1] = 1;
    h_TFGeneIdx[2] = 2;

    h_miValue->print();
    float *d_miValue;
    unsigned int *d_TFGeneIdx;
    HANDLE_ERROR(cudaMalloc((void **)&d_miValue, h_miValue->size()));
    HANDLE_ERROR(cudaMemcpy((void *)d_miValue, (void *)h_miValue->memAddr(), h_miValue->size(), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMalloc((void **)&d_TFGeneIdx, sizeof (unsigned int) * nTFs));
    HANDLE_ERROR(cudaMemcpy((void *)d_TFGeneIdx, (void *)h_TFGeneIdx, sizeof(unsigned int) * nTFs, cudaMemcpyHostToDevice));
    pruneGraph(d_miValue, nTFs, nGenes, d_TFGeneIdx);
    HANDLE_ERROR(cudaMemcpy((void *)h_miValue->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaDeviceSynchronize());
    h_miValue->print();
    delete h_miValue;
    delete h_TFGeneIdx;
    HANDLE_ERROR(cudaFree(d_miValue));
    HANDLE_ERROR(cudaFree(d_TFGeneIdx));
    return 0;
}
