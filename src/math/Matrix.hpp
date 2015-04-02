#ifndef MATRIX_H
#define MATRIX_H
#include <cuda.h>
#define INIT_VALUE -1

template <typename T>
class MiMatrix
{
private:
    std::vector<std::vector<T>> data;
    std::map<unsigned int, unsigned int> gene2Idx;
    unsigned int nRows;
    unsigned int nCols;

public:
    // constructor square matrix
    __host__ __device__ MiMatrix(unsigned int size);
    // constructor general 2-D matrix
    __host__ __device__ MiMatrix(unsigned int m, unsigned int n);
    // copy constructor
    __host__ __device__ MiMatrix(const MiMatrix<T>& rhs);
    // overloading assignment operator
    __host__ __device__ MiMatrix<T>& operator = (const MiMatrix<T>& rhs);
    // destructor
    __host__ __device__ ~MiMatrix();
    // check interaction
    __host__ __device__ bool hasCorrelation(unsigned int geneIdx1, unsigned int geneIdx2);

}


#endif
