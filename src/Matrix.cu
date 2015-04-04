#include "Matrix.hpp"
#define INIT_VALUE -1
#define __CUDA__ __host__ __device__
//#define __CUDA__ 
//matrix implemented in 1-D array

__global__ void rankRow(float *dataStart, unsigned int nRows, unsigned int nCols)
{
    unsigned int rowIdx = blockDim.x * blockIdx.x + threadIdx.x;
    if (rowIdx >= nRows) return;
    float *arrayStart = dataStart + rowIdx * nCols; 
    float *arrayEnd = arrayStart + nCols;

        // work here
        //
        //
        //
        // work here
}

// constructor square matrix
__CUDA__ Matrix::Matrix(unsigned int size)
{
    nRows = size;
    nCols = size;
    data = new float[nRows * nCols];
    for (int i = 0; i < nRows; i++){
        for (int j = 0; j < nCols; j++){
            data[i * nCols + j] = INIT_VALUE;
        }
    }
}

// special constructor called in kernel function
__CUDA__ Matrix::Matrix(unsigned int m, unsigned int n, float *data)
{
    nRows = m; 
    nCols = n;
    // use shallow copy of the pointer
    // data should be on device already
    this->data = data;
}

// constructor general 2-D matrix
__CUDA__ Matrix::Matrix(unsigned int m, unsigned int n)
{
    nRows = m;
    nCols = n;
    data = new float[nRows * nCols];
    for (int i = 0; i < nRows; i++){
        for (int j = 0; j < nCols; j++){
            data[i * nCols + j] = INIT_VALUE;
        }
    }
}

// copy constructor
__CUDA__ Matrix::Matrix(const Matrix& rhs)
{
    nRows = rhs.nRows;
    nCols = rhs.nCols;
    data = new float[nRows * nCols];
    for (int i = 0; i < nRows; i++){
        for (int j = 0; j < nCols; j++){
            data[i * nCols + j] = INIT_VALUE;
        }
    }
}

    // overloading assignment operator
__CUDA__ Matrix& Matrix::operator = (const Matrix& rhs)
{
    if (this->data == rhs.data) return *this;
    delete[] data;
    nRows = rhs.nRows;
    nCols = rhs.nCols;
    data = new float[nRows * nCols];
    for (int i = 0; i < nRows; i++){
        for (int j = 0; j < nCols; j++){
            data[i * nCols + j] = INIT_VALUE;
        }
    }
    return *this;
}

// destructor
__CUDA__ Matrix::~Matrix()
{
    delete[] data;
}
