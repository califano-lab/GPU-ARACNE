#ifndef BOOTSTRAPCONTAINER_HPP
#define BOOTSTRAPCONTAINER_HPP
#include "util.hpp"
#include <cuda.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <cmath>

__global__
void aggregate(float *d_miContainer, unsigned int *d_countContainer, 
        float *d_miValue, unsigned int nTFs, unsigned int nGenes)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nTFs * nGenes) return;
    atomicAdd(d_miContainer + idx, d_miValue[idx] * (d_miValue[idx] > 0));
    atomicAdd(d_countContainer + idx, 1 * (d_miValue[idx] > 0));
}

__global__
void filter(float *d_miContainer, unsigned int *d_countContainer, 
        int threshold, unsigned int nTFs, unsigned int nGenes)
{
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nTFs * nGenes) return;
    d_miContainer[idx] = (d_miContainer[idx] / d_countContainer[idx]) * (d_countContainer[idx] >= threshold);
}

class BootstrapContainer
{
public:
    BootstrapContainer(unsigned int nTFs, unsigned int nGenes)
    {
        this->nTFs = nTFs;
        this->nGenes = nGenes;
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
        aggregate<<<gridDim, blockDim>>>(d_miContainer, d_countContainer, d_miValue, nTFs, nGenes);
        HANDLE_ERROR(cudaGetLastError());
        HANDLE_ERROR(cudaDeviceSynchronize());
    }

    const float *miValueMemAddr() 
    {
        return d_miContainer;
    }

    void condenseGraph(float pValue)
    {
        // this function uses the count matrix and Poisson distribution to determine the final edges
        // the final result is saved back to d_miContainer
        // thrust to reduce (compute sum)
        thrust::device_ptr<unsigned int> w_countContainer (d_countContainer);
        int totalCount = thrust::reduce(w_countContainer, w_countContainer + nTFs * nGenes, 
                (unsigned int)0, thrust::plus<unsigned int>());
        // compute the average count of the graph, which is also called lambda
        float lambda = (float)totalCount / (float)(nTFs * nGenes); 
        // let's hard code Poisson equation
        float currentTerm = 1 / exp(lambda);
        float cd = currentTerm; // cumulative density (just one value not a function)
        int k = 0;
        for (k = 1; k < totalCount; k++){
            if (cd > 1.0 - pValue){
                break;
            } else {
                currentTerm = currentTerm * lambda / k;
                cd += currentTerm;
            }
        }

        // apply threshold k
        dim3 blockDim(1024, 1, 1);
        dim3 gridDim(ceil(nTFs * nGenes / 1024.0), 1, 1);
        filter<<<gridDim, blockDim>>>(d_miContainer, d_countContainer, k, nTFs, nGenes);
        HANDLE_ERROR(cudaGetLastError());
        HANDLE_ERROR(cudaDeviceSynchronize());
    }

private:
    float *d_miContainer;
    unsigned int *d_countContainer;
    unsigned int nTFs;
    unsigned int nGenes;
};

#endif
