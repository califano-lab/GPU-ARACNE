#ifndef UTIL_HPP
#define UTIL_HPP
#include <cuda.h>
#include <cstdio>
#define HANDLE_ERROR( err ) (gpuCheckError( err, __FILE__, __LINE__ ))
static void gpuCheckError(cudaError_t err, const char *file, int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

#define HOST 0
#define DEVICE 1

#endif
