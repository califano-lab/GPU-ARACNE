#ifndef MATRIX_H
#define MATRIX_H
#include <cuda.h>
#define INIT_VALUE -1
#define __CUDA__ __host__ __device__
//#define __CUDA__ 
//matrix implemented in 1-D array
template <typename T>
class Matrix
{
private:
    T *data;
    unsigned int nRows;
    unsigned int nCols;
    __CUDA__ rankRow(T *dataStart)
    {
        unsigned int rowIdx = blockDim.x * blockIdx.x + threadIdx.x;
        if (rowIdx >= nRows) return;
        T *arrayStart = dataStart + rowIdx * nCols 

        // work here
        //
        //
        //
        // work here
    }

public:
    // constructor square matrix
    __CUDA__ Matrix(unsigned int size)
    {
        nRows = size;
        nCols = size;
        data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nCols + j] = INIT_VALUE;
            }
        }
    }

    // special constructor called in kernel function
    __CUDA__ Matrix(unsigned int m, unsigned int n, T *data)
    {
        nRows = m; 
        nCols = n;
        // use shallow copy of the pointer
        // data should be on device already
        this->data = data;
    }

    // constructor general 2-D matrix
    __CUDA__ Matrix(unsigned int m, unsigned int n)
    {
        nRows = m;
        nCols = n;
        data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nCols + j] = INIT_VALUE;
            }
        }
    }

    // copy constructor
    __CUDA__ Matrix(const Matrix<T>& rhs)
    {
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nCols + j] = INIT_VALUE;
            }
        }
    }

    // overloading assignment operator
    __CUDA__ Matrix<T>& operator = (const Matrix<T>& rhs)
    {
        if (this->data == rhs.data) return *this;
        delete[] data;
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nCols + j] = INIT_VALUE;
            }
        }
        return *this;
    }

    // destructor
    __CUDA__ ~Matrix()
    {
        delete[] data;
    }

    // check interaction
    __CUDA__ bool hasCorrelation(unsigned int geneIdx1, unsigned int geneIdx2)
    {
        return !(data[geneIdx1 * nCols + geneIdx2] < 0);
    }

    __CUDA__ void setValue(unsigned int geneIdx1, unsigned int geneIdx2, T val)
    {
        data[geneIdx1 * nCols + geneIdx2] = val;
    }

    // rank each row
    __CUDA__ void rank()
    {
        threadsPerBlock = 512;
        numBlocks = ceil(nRows / (1.0 * threadsPerBlock));
        // it will change data (mutating)
        rankRow<<<numBlocks, threadsPerBlock>>>(data);
    } 

    __CUDA__ T& element(unsigned int geneIdx1, unsigned int geneIdx2)
    {
        return data[geneIdx1 * nCols + geneIdx2];
    }

    const __CUDA__ T* memAddr()
    {
        return data;
    }

    __CUDA__ size_t size()
    {
        return (size_t)(sizeof(T) * nRows * nCols);
    }

    __CUDA__ unsigned int nRows()
    {
        return nRows;
    }

    __CUDA__ unsigned int nCols()
    {
        return nCols;
    }

};


#endif
