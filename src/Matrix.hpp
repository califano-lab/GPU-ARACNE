#ifndef MATRIX_H
#define MATRIX_H
#include <cuda.h>
#include <thrust/sort.h>
#define INIT_VALUE -1
#define __CUDA__ __host__ __device__
//#define __CUDA__ 
//matrix implemented in 1-D array
__global__ void rankRow(float *dataStart, unsigned int nRows, unsigned int nCols);

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
    __host__ void rank()
    {
        unsigned int threadsPerBlock = 512;
        unsigned int numBlocks = ceil(nRows / (1.0 * threadsPerBlock));
        // it will change data (mutating)
        rankRow<<<numBlocks, threadsPerBlock>>>(data, nRows, nCols);
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

};


#endif
