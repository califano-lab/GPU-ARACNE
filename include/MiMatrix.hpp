#ifndef MATRIX_H
#define MATRIX_H
#include <cuda.h>
#define INIT_VALUE -1
#define __CUDA__ __host__ __device__
//#define __CUDA__ 
template <typename T>
class MiMatrix
{
private:
    T **data;
    unsigned int nRows;
    unsigned int nCols;

public:
    // constructor square matrix
    __CUDA__ MiMatrix(unsigned int size)
    {
        nRows = size;
        nCols = size;
        data = new T*[size];
        for (int i = 0; i < size; i++){
            data[i] = new T[size];
            for (int j = 0; j < size; j++){
                data[i][j] = INIT_VALUE;
            }
        }
    }

    // constructor general 2-D matrix
    __CUDA__ MiMatrix(unsigned int m, unsigned int n)
    {
        nRows = m;
        nCols = n;
        data = new T*[m];
        for (int i = 0; i < m; i++){
            data[i] = new T[n];
            for (int j = 0; j < n; j++){
                data[i][j] = INIT_VALUE;
            }
        }
    }

    // copy constructor
    __CUDA__ MiMatrix(const MiMatrix<T>& rhs)
    {
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        data = new T*[nRows];
        for (int i = 0; i < nRows; i++){
            data[i] = new T[nCols];
            for (int j = 0; j < nCols; j++){
                data[i][j] = INIT_VALUE;
            }
        }
    }

    // overloading assignment operator
    __CUDA__ MiMatrix<T>& operator = (const MiMatrix<T>& rhs)
    {
        if (this->data == rhs.data) return *this;
        for (int i = 0; i < nRows; i++){
            delete[] data[i];
        }
        delete[] data;
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        data = new T*[nRows];
        for (int i = 0; i < nRows; i++){
            data[i] = new T[nCols];
            for (int j = 0; j < nCols; j++){
                data[i][j] = INIT_VALUE;
            }
        }
        return *this;
    }

    // destructor
    __CUDA__ ~MiMatrix()
    {
        for (int i = 0; i < nRows; i++){
            delete[] data[i];
        }
        delete[] data;
    }

    // check interaction
    __CUDA__ bool hasCorrelation(unsigned int geneIdx1, unsigned int geneIdx2)
    {
        return !(data[geneIdx1][geneIdx2] < 0);
    }

};


#endif
