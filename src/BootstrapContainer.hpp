#ifndef BOOTSTRAPCONTAINER_HPP
#define BOOTSTRAPCONTAINER_HPP
#include "util.hpp"
#include <cuda.h>

__global__
void aggregate(float *d_miContainer, unsigned int *d_miCounter, float *d_miValue, unsigned int nTFs, unsigned int nGenes)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nTFs * nGenes) return;
    atomicAdd(d_miContainer + idx, d_miValue[idx]);
    atomicAdd(&d_counterContainer + idx, 1);
}

class BootstrapContainer
{
public:
    BootstrapContainer(unsigned int nTFs, unsigned int nGenes)
    {
        HANDLE_ERROR(cudaMalloc((void **)&d_miContainer, sizeof(float) * nTFs * nGenes));
        HANDLE_ERROR(cudaMemset((void *)d_miContainer, (float)0, sizeof(float) * nTFs * nGenes));
        HANDLE_ERROR(cudaMalloc((void **)&d_countContainer, sizeof(unsigned int) * nTFs * nGenes));
        HANDLE_ERROR(cudaMemset((void *)d_countContainer, 0, sizeof(unsigned int) * nTFs * nGenes));
    }

    ~BootstrapContainer()
    {
        HANDLE_ERROR(cudaFree(d_miContainer));
        HANDLE_ERROR(cudaFree(d_countContainer));
    }

    void addToContainer(float *d_miValue)
    {
        dim3 blockDim(1024, 1, 1);
        dim3 gridDim(ceil(nTFs * nGenes / 1024.0), 1, 1);
        aggregate<<<gridDim, blockDim>>>(d_miContainer, d_miCounter, d_miValue, nTFs, nGenes);
        HANDLE_ERROR(cudaGetLastError());
        HANDLE_ERROR(cudaDeviceSynchronize());
    }

    const float *miValueMemAddr() 
    {
        return d_miContainer;
    }

    void condenseGraph(float lambda)
    {
        // ask Jing
        // this function uses the count matrix and Poisson distribution to determine the final edges
        // the final result is saved back to d_miContainer
        // TODO
    }

private:
    float *d_miContainer;
    unsigned int *d_countContainer;
    unsigned int nTFs;
    unsigned int nGenes;
}

#endif
