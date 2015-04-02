#ifndef MATRIX_H
#define MATRIX_H
#include <cuda.h>
#define INIT_VALUE -1
#define __CUDA__ __host__ __device__
//#define __CUDA__ 
//matrix implemented in 1-D array
template <typename T>
class MiMatrix
{
private:
    T *data;
    unsigned int nRows;
    unsigned int nCols;

public:
    // constructor square matrix
    __CUDA__ MiMatrix(unsigned int size)
    {
        nRows = size;
        nCols = size;
        data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nRows + j] = INIT_VALUE;
            }
        }
    }

    // constructor general 2-D matrix
    __CUDA__ MiMatrix(unsigned int m, unsigned int n)
    {
        nRows = m;
        nCols = n;
        data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nRows + j] = INIT_VALUE;
            }
        }
    }

    // copy constructor
    __CUDA__ MiMatrix(const MiMatrix<T>& rhs)
    {
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nRows + j] = INIT_VALUE;
            }
        }
    }

    // overloading assignment operator
    __CUDA__ MiMatrix<T>& operator = (const MiMatrix<T>& rhs)
    {
        if (this->data == rhs.data) return *this;
        delete[] data;
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nRows + j] = INIT_VALUE;
            }
        }
        return *this;
    }

    // destructor
    __CUDA__ ~MiMatrix()
    {
        delete[] data;
    }

    // check interaction
    __CUDA__ bool hasCorrelation(unsigned int geneIdx1, unsigned int geneIdx2)
    {
        return !(data[geneIdx1 * nRows + geneIdx2] < 0);
    }

    __CUDA__ void setValue(unsigned int geneIdx1, unsigned int geneIdx2, T val)
    {
        data[geneIdx1 * nRows + geneIdx2] = val;
    }

    __CUDA__ T& element(unsigned int geneIdx1, unsigned int geneIdx2)
    {
        return data[geneIdx1 * nRows + geneIdx2];
    }

    const __CUDA__ T* memAddr()
    {
        return data;
    }

    __CUDA__ size_t size()
    {
        return (size_t)(sizeof(T) * nRows * nCols);
    }

};


#endif
