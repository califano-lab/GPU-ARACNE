#include "Matrix.hpp"
#define INIT_VALUE -1
#define __CUDA__ __host__ __device__
//#define __CUDA__ 
//matrix implemented in 1-D array

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
__host__ Matrix& Matrix::operator = (const Matrix& rhs)
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
    int status;
    cuPointerGetAttribute((void *)&status, CU_POINTER_ATTRIBUTE_MEMORY_TYPE, data);
    if (status == CU_MEMORYTYPE_HOST) 
        delete[] data;
    else
        cudaFree(data);
}

__host__ void Matrix::print()
{
    for (int i = 0; i < nRows; i++){
        for (int j = 0; j < nCols; j++){
            std::cout << this->element(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
