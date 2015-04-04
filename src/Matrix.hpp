#ifndef MATRIX_H
#define MATRIX_H
#include <cstdlib>
#include <cuda.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#define INIT_VALUE -1
#define __CUDA__ __host__ __device__
//#define __CUDA__ 
//matrix implemented in 1-D array

__global__ void rankRow(unsigned int *d_data, unsigned int nRows, unsigned int nCols)
{
    unsigned int rowIdx = blockDim.x * blockIdx.x + threadIdx.x;
    // this is slow, think to optimize
    unsigned int *helperArray = new unsigned int[nCols];
    for (int i = 0; i < nCols; i++){
        helperArray[i] = i;
    }
    thrust::stable_sort_by_key(d_data+rowIdx*nCols + d_data+rowIdx*nCols + nCols, helperArray);
    for (int i = 0; i < nCols; i++){
        d_data[rowIdx * nCols + helperArray[i]] = i;
    }
}

class Matrix
{
private:
    float *data;
    unsigned int nRows;
    unsigned int nCols;
   
public:
    // constructor square matrix
    __CUDA__ Matrix(unsigned int size);

    // special constructor called in kernel function
    __CUDA__ Matrix(unsigned int m, unsigned int n, float *data);

    // constructor general 2-D matrix
    __CUDA__ Matrix(unsigned int m, unsigned int n);

    // copy constructor
    __CUDA__ Matrix(const Matrix& rhs);
    
    // overloading assignment operator
    __CUDA__ Matrix& operator = (const Matrix& rhs);
   
    // destructor
    __CUDA__ ~Matrix();

    // check interaction
    __CUDA__ bool hasCorrelation(unsigned int geneIdx1, unsigned int geneIdx2)
    {
        return !(data[geneIdx1 * nCols + geneIdx2] < 0);
    }

    __CUDA__ void setValue(unsigned int geneIdx1, unsigned int geneIdx2, float val)
    {
        data[geneIdx1 * nCols + geneIdx2] = val;
    }

    // rank each row
    __host__ Matrix *getRankMatrix()
    {
        unsigned int *d_data;
        cudaMalloc((void **)&d_data, size());
        cudaMemcpy((void *)d_data, (void *)data, cudaMemcpyHostToDevice);
        unsigned int threadsPerBlock = 1024;
        rankRow<<<ceil(nRows / (1.0*threadsPerBlock)), threadsPerBlock>>>
            (d_data, nRows, nCols);
        return new Matrix(nRows, nCols, d_data);
    }
    
    __CUDA__ float& element(unsigned int geneIdx1, unsigned int geneIdx2)
    {
        return data[geneIdx1 * nCols + geneIdx2];
    }

    const __CUDA__ float* memAddr()
    {
        return data;
    }

    __CUDA__ size_t size()
    {
        return (size_t)(sizeof(float) * nRows * nCols);
    }

    __CUDA__ unsigned int getNumRows()
    {
        return nRows;
    }

    __CUDA__ unsigned int getNumCols()
    {
        return nCols;
    }

    __host__ void print();
};


#endif
