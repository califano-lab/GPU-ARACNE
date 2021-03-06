#include "miAP.hpp"

// Chi Square statistical test
__device__ 
bool chiSeq(unsigned int nSamples, unsigned int a, unsigned int b, unsigned int c, unsigned int d){
    if (a + b + c + d == nSamples) return true; // first cube
    if (a + b + c + d < 8) return false;
    
    float expected = (a + b + c + d) / 4.0;
    float testStat = ((a - expected) * (a - expected) + (b - expected) * (b - expected) + 
        (c - expected) * (c - expected) + (d - expected) * (d - expected)) / expected;
    if (testStat < CHI_THRESHOLD) return false; // not significant
    else return true;
}

// called by experiment miAP()
__global__ 
void computeMi(unsigned int *d_rankMatrix, unsigned int nTFs, unsigned int nGenes, unsigned int nSamples, 
        unsigned int *d_TFGenesIdx, float *d_rawGraph, float miThreshold)
{
    unsigned int TFIdx = blockIdx.x;
    unsigned int geneIdx = blockIdx.y;
    // need pick two rows of the d_rankMatrix
    // first row (define it X): d_TFGenesIdx[TFIdx]
    // second row (define it Y): geneIdx
    extern __shared__ Cube cubeArray[];
    unsigned int head = 0;
    unsigned int tail = 1;
    Cube cube;
    cube.totalCount = nSamples;
    cube.upper = nSamples - 1;
    cube.lower = 0;
    cube.left = 0;
    cube.right = nSamples - 1;
    cubeArray[head] = cube;
    __shared__ unsigned int upperRight;
    __shared__ unsigned int upperLeft;
    __shared__ unsigned int lowerLeft; 
    __shared__ unsigned int lowerRight;
    float miValue = 0;
    __syncthreads();
    do{
        upperRight = 0;
        upperLeft = 0;
        lowerLeft = 0;
        lowerRight = 0;
        unsigned int colIdx = threadIdx.x; // at most 1024 samples
        unsigned int middleX = (cubeArray[head].left + cubeArray[head].right) / 2;
        unsigned int middleY = (cubeArray[head].upper + cubeArray[head].lower) / 2;
        unsigned int coordX = (d_rankMatrix[d_TFGenesIdx[TFIdx] * nSamples + colIdx]);
        unsigned int coordY = (d_rankMatrix[geneIdx * nSamples + colIdx]);

        // optimize to reduce thread divergence
        atomicAdd(&lowerLeft, cubeArray[head].left <= coordX & coordX <= middleX & 
                cubeArray[head].lower <= coordY & coordY <= middleY);
        atomicAdd(&upperLeft, cubeArray[head].left <= coordX & coordX <= middleX & 
                middleY < coordY & coordY <= cubeArray[head].upper);
        atomicAdd(&lowerRight, middleX < coordX & coordX <= cubeArray[head].right &
                cubeArray[head].lower <= coordY & coordY <= middleY);
        atomicAdd(&upperRight, middleX < coordX & coordX <= cubeArray[head].right &
                middleY < coordY & coordY <= cubeArray[head].upper);
        __syncthreads();
        
        if (chiSeq(nSamples, upperRight, upperLeft, lowerLeft, lowerRight)) {// significant
            cubeArray[tail].totalCount = upperRight;
            cubeArray[tail].upper = cubeArray[head].upper;
            cubeArray[tail].lower = middleY + 1;
            cubeArray[tail].left = middleX + 1;
            cubeArray[tail].right = cubeArray[head].right;
            
            cubeArray[tail + 1].totalCount = upperLeft;
            cubeArray[tail + 1].upper = cubeArray[head].upper;
            cubeArray[tail + 1].lower = middleY + 1;
            cubeArray[tail + 1].left = cubeArray[head].left;
            cubeArray[tail + 1].right = middleX;
            
            cubeArray[tail + 2].totalCount = lowerLeft;
            cubeArray[tail + 2].upper = middleY;
            cubeArray[tail + 2].lower = cubeArray[head].lower;
            cubeArray[tail + 2].left = cubeArray[head].left;
            cubeArray[tail + 2].right = middleX;
            
            cubeArray[tail + 3].totalCount = lowerRight;
            cubeArray[tail + 3].upper = middleY;
            cubeArray[tail + 3].lower = cubeArray[head].lower;
            cubeArray[tail + 3].left = middleX + 1;
            cubeArray[tail + 3].right = cubeArray[head].right;
            // no locking needed
            // head and tail on register        
            tail = tail + 4;
            head = head + 1;
        } else {
            // compute the actual mutual information value here
            unsigned int countX = cubeArray[head].right - cubeArray[head].left + 1;
            unsigned int countY = cubeArray[head].upper - cubeArray[head].lower + 1;
            float logRight;
            if (cubeArray[head].totalCount == 0)
                logRight = 1;
            else
                logRight = (float)cubeArray[head].totalCount / (float)(countX * countY) ;
            // miValue and head on register
            // no locking needed
            miValue += (float)cubeArray[head].totalCount *  log(logRight);
            head = head + 1;
        }
        __syncthreads();
    } while(head < tail);
    d_rawGraph[TFIdx * nGenes + geneIdx] = 
        (miValue / nSamples + log((float)nSamples)  > miThreshold) * (miValue / nSamples + log((float)nSamples)); 
}

// called by null model miAP()
__global__
void computeMi(unsigned int *d_randomMatrix, unsigned int nPairs, unsigned int nSamples, float *d_miResult)
{
    unsigned int rowIdx = blockIdx.x;
    // need pick two rows 
    // first row (define it X): rowIdx
    // second row (define it Y): rowIdx * 2 
    extern __shared__ Cube cubeArray[];
    unsigned int head = 0;
    unsigned int tail = 1;
    Cube cube;
    cube.totalCount = nSamples;
    cube.upper = nSamples - 1;
    cube.lower = 0;
    cube.left = 0;
    cube.right = nSamples - 1;
    cubeArray[head] = cube;
    __shared__ unsigned int upperRight;
    __shared__ unsigned int upperLeft;
    __shared__ unsigned int lowerLeft; 
    __shared__ unsigned int lowerRight;
    float miValue = 0;
    __syncthreads();
    do{
        upperRight = 0;
        upperLeft = 0;
        lowerLeft = 0;
        lowerRight = 0;
        unsigned int colIdx = threadIdx.x;
        unsigned int middleX = (cubeArray[head].left + cubeArray[head].right) / 2;
        unsigned int middleY = (cubeArray[head].upper + cubeArray[head].lower) / 2;
        unsigned int coordX = d_randomMatrix[rowIdx * nSamples + colIdx];
        unsigned int coordY = d_randomMatrix[rowIdx * 2 * nSamples + colIdx];
        
        // optimize to reduce thread divergence
        atomicAdd(&lowerLeft, cubeArray[head].left <= coordX & coordX <= middleX & 
                cubeArray[head].lower <= coordY & coordY <= middleY);
        atomicAdd(&upperLeft, cubeArray[head].left <= coordX & coordX <= middleX & 
                middleY < coordY & coordY <= cubeArray[head].upper);
        atomicAdd(&lowerRight, middleX < coordX & coordX <= cubeArray[head].right &
                cubeArray[head].lower <= coordY & coordY <= middleY);
        atomicAdd(&upperRight, middleX < coordX & coordX <= cubeArray[head].right &
                middleY < coordY & coordY <= cubeArray[head].upper);
        __syncthreads();
      
        if (chiSeq(nSamples, upperRight, upperLeft, lowerLeft, lowerRight)) {// significant
            cubeArray[tail].totalCount = upperRight;
            cubeArray[tail].upper = cubeArray[head].upper;
            cubeArray[tail].lower = middleY;
            cubeArray[tail].left = middleX;
            cubeArray[tail].right = cubeArray[head].right;
            
            cubeArray[tail + 1].totalCount = upperLeft;
            cubeArray[tail + 1].upper = cubeArray[head].upper;
            cubeArray[tail + 1].lower = middleY;
            cubeArray[tail + 1].left = cubeArray[head].left;
            cubeArray[tail + 1].right = middleX - 1;
            
            cubeArray[tail + 2].totalCount = lowerLeft;
            cubeArray[tail + 2].upper = middleY - 1;
            cubeArray[tail + 2].lower = cubeArray[head].lower;
            cubeArray[tail + 2].left = cubeArray[head].left;
            cubeArray[tail + 2].right = middleX - 1;
            
            cubeArray[tail + 3].totalCount = lowerRight;
            cubeArray[tail + 3].upper = middleY - 1;
            cubeArray[tail + 3].lower = cubeArray[head].lower;
            cubeArray[tail + 3].left = middleX;
            cubeArray[tail + 3].right = cubeArray[head].right;
            // no locking needed
            // head and tail on register        
            tail = tail + 4;
            head = head + 1;            
        } else {
            // compute the actual mutual information value here
            unsigned int countX = cubeArray[head].right - cubeArray[head].left + 1;
            unsigned int countY = cubeArray[head].upper - cubeArray[head].lower + 1;
            float logRight;
            if (cubeArray[head].totalCount == 0)
                logRight = 1;
            else
                logRight = (float)cubeArray[head].totalCount / (float)(countX * countY);
            // miValue and head on register
            // no locking needed
            miValue += (float)cubeArray[head].totalCount *  log(logRight);
            head = head + 1;
        }
        __syncthreads();
    } while(head < tail);
    d_miResult[rowIdx] = miValue /nSamples + log((float)nSamples) ;
}

// this function computes the mutual information given the ranked experimental matrix
// X and Y are mapped faciliated by d_TFGenesIdx (see details in source code)
// function called during actual MI computation 
void miAP(unsigned int *d_rankMatrix, unsigned int nTFs, unsigned int nGenes, unsigned int nSamples, 
        unsigned int *d_TFGeneIdx, float **d_rawGraph, float miThreshold)
{
    dim3 gridDim(nTFs, nGenes, 1);
    dim3 blockDim(nSamples, 1, 1);

    HANDLE_ERROR( cudaMalloc((void **)d_rawGraph, sizeof(float) * nTFs * nGenes) );
    computeMi<<<gridDim, blockDim, nSamples * sizeof(Cube)>>>
        (d_rankMatrix, nTFs, nGenes, nSamples, d_TFGeneIdx, *d_rawGraph, miThreshold);
    HANDLE_ERROR( cudaGetLastError() );
    HANDLE_ERROR( cudaDeviceSynchronize() );
}
 
// this function computes the mutual information given random matrix
// the first half of the matrix is regarded as X and the second half Y
// function called during null model computation
// input d_randomMatix should be 2*nPairs by nSamples matrix
// d_miResult, just an address of a pointer to T is fine
void miAP(unsigned int *d_randomMatrix, unsigned int nPairs, unsigned int nSamples, float **d_miResult)
{
    dim3 gridDim(nPairs, 1, 1);
    dim3 blockDim(nSamples, 1, 1);

    HANDLE_ERROR( cudaMalloc((void **)d_miResult, sizeof(float) * nPairs) );
    computeMi<<<gridDim, blockDim, nSamples * sizeof(Cube)>>>(d_randomMatrix, nPairs, nSamples, *d_miResult);
    HANDLE_ERROR( cudaGetLastError() );
    HANDLE_ERROR( cudaDeviceSynchronize() );
}

