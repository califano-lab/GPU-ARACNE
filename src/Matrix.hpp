#ifndef MATRIX_H
#define MATRIX_H
#include <cstdlib>
#include <cuda.h>
#include <thrust/sort.h>
#define INIT_VALUE -1
#define __CUDA__ __host__ __device__
//#define __CUDA__ 
//matrix implemented in 1-D array

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
        // all run on CPU with thrust facility on GPU
        Matrix *helperMatrix = new Matrix(nRows, nCols);
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                helperMatrix->element(i, j) = j;
            }
        }
        for (int i = 0; i < nRows; i++){
            thrust::stable_sort_by_key(this->data + i * nCols, 
                    this->data + i * nCols + nCols, helperMatrix->data + i * nCols);
        }
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                this->element(i, helperMatrix->element(i, j)) = j;
            }
        }
        delete helperMatrix;
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
