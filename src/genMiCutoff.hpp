//
//  genMiCutoff.h
//  aracneGPU
//
//  Created by jing on 3/24/15.
//  Copyright (c) 2015 jing. All rights reserved.
//

// #ifndef __aracneGPU__genMiCutoff__
// #define __aracneGPU__genMiCutoff__

// compile with: nvcc -arch=sm_20 -lcurand -o t89 t89.cu
#include <stdio.h>
#include <curand.h>
#include <curand_kernel.h>



#define GPU_CHECKERROR( err ) (gpuCheckError( err, __FILE__, __LINE__ ))
static void gpuCheckError( cudaError_t err,
                          const char *file,
                          int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
               file, line );
        exit( EXIT_FAILURE );
    }
}


__device__ 
float getnextrand(curandState *state){

  return (float)(curand_uniform(state));
}

__device__ 
int getnextrandscaled(curandState *state, int scale){

  return (int) scale * getnextrand(state);
}


__global__ 
void initCurand(curandState *state, unsigned long seed){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(seed, 0, 0, &state[idx]);
}

__global__ 
void testrand(curandState *state, int *a1, int *a2){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;

    a1[idx] = getnextrandscaled(&state[idx], SCALE);
    a2[idx] = getnextrandscaled(&state[idx], SCALE);
}

int main() {


	SCALE = 49000; 
	DSIZE = 5000; 
	nTPB  = 256; 

    int *h_a1, *h_a2, *d_a1, *d_a2;

    GPU_CHECKERROR( curandState *devState );

    h_a1 = (int *)malloc( DSIZE * sizeof(int));
    if (h_a1 == 0) { printf("malloc fail\n"); return 1; }
    h_a2 = (int *)malloc( DSIZE * sizeof(int));
    if (h_a2 == 0) {printf("malloc fail\n"); return 1;}

    GPU_CHECKERROR( cudaMalloc((void**)&d_a1, DSIZE * sizeof(int)) );
    GPU_CHECKERROR( cudaMalloc((void**)&d_a2, DSIZE * sizeof(int)) );
    GPU_CHECKERROR( cudaMalloc((void**)&devState, DSIZE * sizeof(curandState)) );
    
    GPU_CHECKERROR( initCurand<<< (DSIZE+nTPB-1)/nTPB, nTPB >>>(devState, 1) );
    GPU_CHECKERROR( cudaDeviceSynchronize() );

    GPU_CHECKERROR( testrand<<< (DSIZE+nTPB-1)/nTPB, nTPB >>>(devState, d_a1, d_a2) );
     GPU_CHECKERROR( cudaDeviceSynchronize() );

     GPU_CHECKERROR( cudaMemcpy(h_a1, d_a1, DSIZE*sizeof(int), cudaMemcpyDeviceToHost) );
     GPU_CHECKERROR( cudaMemcpy(h_a2, d_a2, DSIZE*sizeof(int), cudaMemcpyDeviceToHost) );
     cudaCheckErrors("cudamemcpy");
     printf("1st returned random value is %d\n", h_a1[0]);
     printf("2nd returned random value is %d\n", h_a2[0]);

     for (int i=1; i< DSIZE; i++){
       if (h_a1[i] != h_a1[0]) {
         printf("mismatch on 1st value at %d, val = %d\n", i, h_a1[i]);
         return 1;
         }
       if (h_a2[i] != h_a2[0]) {
         printf("mismatch on 2nd value at %d, val = %d\n", i, h_a2[i]);
         return 1;
         }
       }
     printf("thread values match!\n");

}

// #endif /* defined(__aracneGPU__genMiCutoff__) */
