#ifndef MATRIX_HPP
#define MATRIX_HPP
#include "util.hpp"
#include <cstdlib>

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

#include <cuda.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

#define INIT_VALUE -1
#define __CUDA__ __host__ __device__

//#define __CUDA__ 
//matrix implemented in 1-D array

template <typename T>
__global__ void rankRow(T *d_data, unsigned int tabletSize, unsigned int nCols, unsigned int *d_rank)
{
    unsigned int tabIdx = blockDim.x * blockIdx.x + threadIdx.x;
    if (tabIdx >= tabletSize) return;
    // this is slow, think to optimize
    unsigned int *helperArray = new unsigned int[nCols];
    for (int i = 0; i < nCols; i++){
        helperArray[i] = i;
    }

    // rank elements in each row
    /*
     * this chunk of code can replace the next chunk of code the system support CUDA dynamic parallelism
    thrust::stable_sort_by_key(thrust::seq, d_data+tabIdx*nCols, d_data+tabIdx*nCols + nCols, helperArray);
    */

    for (int i = 0; i < nCols-1; i++){
        for (int j = i+1; j < nCols; j++) {
            if (*(d_data+tabIdx*nCols+i) > *(d_data+tabIdx*nCols+j)) {
                int temp = *(d_data+tabIdx*nCols+i);
    		*(d_data+tabIdx*nCols+i) = *(d_data+tabIdx*nCols+j);
		*(d_data+tabIdx*nCols+j) = temp;
                temp = helperArray[i];
                helperArray[i] = helperArray[j];
                helperArray[j] = temp;
	    }
	}
    }
    for (int i = 0; i < nCols; i++){
        d_rank[tabIdx * nCols + helperArray[i]] = i;
    }

    delete[] helperArray;
}

template <typename T>
class Matrix
{
private:
    T *data;
    unsigned int nRows;
    unsigned int nCols;
    // type of data: HOST or DEVICE
   
public:
    // constructor square matrix
    __host__ Matrix(unsigned int size)
    {
        nRows = size;
        nCols = size;
        HANDLE_ERROR(cudaHostAlloc((void **)&data, size * size * sizeof(T), 0));
        //data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nCols + j] = (T)INIT_VALUE;
            }
        }
    }

    // special constructor called in kernel function
    __host__ Matrix(unsigned int m, unsigned int n, T *data)
    {
        nRows = m; 
        nCols = n;
        // use shallow copy of the pointer
        // data should be on device already
        this->data = data;
    }

    // constructor general 2-D matrix
    __host__ Matrix(unsigned int m, unsigned int n)
    {
        nRows = m;
        nCols = n;
        HANDLE_ERROR(cudaHostAlloc((void **)&data, nRows * nCols * sizeof(T), cudaHostAllocDefault));
        //data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nCols + j] = (T)INIT_VALUE;
            }
        }
    }

    // copy constructor
    __host__ Matrix(const Matrix& rhs)
    {
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        HANDLE_ERROR(cudaHostAlloc((void **)&data, nRows * nCols * sizeof(T), cudaHostAllocDefault));
        //data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nCols + j] = rhs.data[i * nCols + j];
            }
        }
    }
    
    // overloading assignment operator
    __host__ Matrix& operator = (const Matrix<T>& rhs)
    {
        if (this->data == rhs.data) return *this;
        HANDLE_ERROR(cudaFreeHost(data));
        nRows = rhs.nRows;
        nCols = rhs.nCols;
        HANDLE_ERROR(cudaHostAlloc((void **)&data, nRows * nCols * sizeof(T), cudaHostAllocDefault));
        //data = new T[nRows * nCols];
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                data[i * nCols + j] = rhs.data[i * nCols + j];
            }
        }
        return *this;
    }

    // destructor
    __host__ ~Matrix()
    {
        HANDLE_ERROR(cudaFreeHost(data));
        //delete[] data;
    }

    // check interaction
    __host__ bool hasCorrelation(unsigned int geneIdx1, unsigned int geneIdx2)
    {
        return !(data[geneIdx1 * nCols + geneIdx2] < 0);
    }

    __host__ Matrix<T> *getCols( int *colsIdx, unsigned int numCols ) 
    {
	// if ( numCols == 0 ) return; 
        Matrix<T> *subMatrix = new Matrix<T>(nRows,  numCols);
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < numCols; ++j ){
                subMatrix->setValue( i, j, data[ i * nCols + colsIdx[j] ] );
            }
        }
        return subMatrix;
    } 

    __host__ Matrix<T> *bootstrapMatrix( ) 
    {
        srand(time(NULL));
        Matrix<T> *bsMatrix = new Matrix<T>(nRows, nCols);
	int bsj; 
        for (int j = 0; j < nCols; ++j ){
	    bsj = rand() % (int)nCols;
	    for (int i = 0; i < nRows; i++){
		// generate random indexs with replicates
                // bsMatrix->setValue( i, j, data[ i * nCols + bsj ] );
                // direct handle mem address more efficient
                bsMatrix->data[i * nCols + j] = data[i * nCols + bsj];
            }
        }
        return bsMatrix;
    } 
 
   
    __host__ void permute(unsigned int seed) 
    {
        std::srand(seed);
	std::vector<int> baseSeq(nCols); 
	for (int i = 0; i < nCols; i++)
	  baseSeq[i] = i;

	for (int i = 0; i < nRows; i++) 
	{
	    std::random_shuffle(baseSeq.begin(), baseSeq.end());
	    for (int j = 0; j < nCols; j++) 
	    {
                data[i * nCols + j] = baseSeq[j];
	    }
        }
    }

    __host__ void setValue(unsigned int geneIdx1, unsigned int geneIdx2, T val)
    {
        data[geneIdx1 * nCols + geneIdx2] = val;
    }

    // rank each row
    __host__ unsigned int *getRankMatrix()
    {
        T *d_data;
        cudaMalloc((void **)&d_data, size());
        HANDLE_ERROR(cudaMemcpy((void *)d_data, (void *)data, size(), cudaMemcpyHostToDevice));
        unsigned int *d_rank;
        cudaMalloc((void **)&d_rank, nRows * nCols * sizeof(int));
        unsigned int tabletSize = 128;
        unsigned int nRows_cpy = nRows;
        unsigned int realTabletSize; 
        for (int i = 0; i < ceil(nRows / (1.0*tabletSize)); i++){
            if (nRows_cpy >= tabletSize){
                realTabletSize = tabletSize;
                nRows_cpy = nRows_cpy - tabletSize;
            } else {
                realTabletSize = nRows_cpy;
            }

            dim3 threadsPerBlock(128, 1, 1);
            dim3 blocksPerGrid((unsigned int)ceil(realTabletSize/(1.0 * 128)), 1, 1);
            rankRow<<<blocksPerGrid, threadsPerBlock>>>
                (d_data + i * tabletSize * nCols, realTabletSize, nCols, d_rank + i * tabletSize * nCols);
            HANDLE_ERROR(cudaGetLastError());
            HANDLE_ERROR(cudaDeviceSynchronize());
        }
        HANDLE_ERROR(cudaFree(d_data));
        return d_rank;
    }
    
    __host__ T& element(unsigned int geneIdx1, unsigned int geneIdx2)
    {
        return data[geneIdx1 * nCols + geneIdx2];
    }

    const __host__ T* memAddr()
    {
        return data;
    }

    __host__ size_t size()
    {
        return (size_t)(sizeof(T) * nRows * nCols);
    }

    __host__ unsigned int getNumRows()
    {
        return nRows;
    }

    __host__ unsigned int getNumCols()
    {
        return nCols;
    }
    
    __host__ void print()
    {
        for (int i = 0; i < nRows; i++){
            for (int j = 0; j < nCols; j++){
                std::cout << this->element(i, j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    __host__ void printHead()
    {
        for (int i = 0; i < 5; i++){
            for (int j = 0; j < 5; j++){
                std::cout << this->element(i, j) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

};

#endif
