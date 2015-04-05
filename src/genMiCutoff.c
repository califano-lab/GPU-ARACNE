//
//  genMiCutoff.cpp
//  aracneGPU
//
//  Created by jing on 3/24/15.
//  Copyright (c) 2015 jing. All rights reserved.
//

#include <iostream>
#include <stdlib.h> 
#include <stdio.h>
#include <time.h> 

#include "genMiCutoff.hpp"

//__device__ __host__ 
int* genRandIntArray ( int array[], const int N )
{
 // init random seed 
   // struct timeval tv;
   // gettimeofday( &tv, NULL);
   // int usec = tv.tv_usec ; 
   // srand(usec);

   printf("Total: %d\n", N);
     
   int *a ;
   a = ( int * ) malloc( N * sizeof( int ) ) ;
   for ( int i=0; i< N; ++i)
   {
     a[i] = i + 1;
   }

   for ( int i = N - 1; i > 0; --i)
   {
     size_t j = (unsigned int ) ( drand48() * (i+ 1));
     int t = a[j] ; 
     a[j] = a[i]; 
     a[i] = t; 
   }
   for ( int i = 0; i < N; ++i ) 
   {
     printf("%d\t",  a[i]) ;
   }
   printf( "\n");
   return a;
}

float calMIcutoff ( float * mi, const int  NTOPPERCENT, const float pCut )
{
  float coeff = [2];
  return coeff;
}

__global__
float calMIap_d (int *X, int *Y, const int N){
  //
  
 if ( gridIdx.x > N | gridIdx.y > N) return; 

 // init shared mem;
   
}

float calMIcutCoeff (const int Nsmp, const int Nperm) {
    // generate  permutation arrays
    int *h_X;
    h_X = (int * ) malloc( sizeof (int) * Nsmp);
    h_X = genRandIntArray( h_X, Nsmp );

    int *h_Y;
    h_Y = (int * ) malloc( sizeof (int) * Nsmp);
    h_Y = genRandIntArray( h_Y, Nsmp );
   
//    // Grid size and block Size
//    int numBlockx = 1; // number of pairs 
//    int numBlocky = 1;
//
//    int tileSize = 5;
//    int numThreadXPerBlock = Nsmp / tileSize ; 
//    int numThreadYPerBlock = Nsmp / tileSize ; 
//
//    // allocate memory, and copy the data over:
//
//    int *h_mi;  
//    h_mi = (float *) malloc( sizeof( float) * numBlockx * numBlocky);  
//    
//    int *d_X,  *d_Y;
//    int *d_mi;
//     
//    GPU_CHECKERROR( cudaMalloc ((void **) &d_X, Nsmp * sizeof (int)) );
//    GPU_CHECKERROR( cudaMalloc ((void **) &d_Y, Nsmp * sizeof (int)) );
//      
//    GPU_CHECKERROR( cudaMemcpy(d_X, h_X, Nsmp * sizeof(int), cudaMemcpyHostToDevice) ) ;
//    GPU_CHECKERROR( cudaMemcpy(d_Y, h_Y, Nsmp * sizeof(int), cudaMemcpyHostToDevice) ) ;
//        
//    // launch the kernel:
//
//    calMIap_d <<< (numBlockx, numBlocky, 1), (numThreadXPerBlock, numThreadYPerBlock, 1) >>> ( float * d_mi, int * d_X, int * d_Y, int Nsmp);
//    
//
//    // copy results back
//    GPU_CHECKERROR( cudaDeviceSynchronize() );
//    GPU_CHECKERROR( cudaMemcpy(h_mi,d_mi, numBlockx * numBlocky * sizeof(float), cudaMemcpyDeviceToHost) ) ;
//
//    // Postprocessing to generate cutoff
//    free(h_X); 
//    free(h_Y); 
//    free(h_mi); 
//    cudafree(d_X);
//    cudafree(d_Y);
//    cudafree(d_mi);

    return 0.0;
}

