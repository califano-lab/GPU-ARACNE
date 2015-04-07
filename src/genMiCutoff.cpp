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
//#include <curand.h>
//#include <curand_kernel.h>
#include <time.h> 

#include "genMiCutoff.hpp"

using namespace std;

#define NPERM 100000 // number of randomization
#define NTOPPERCENT 10 // percentage of top MI's used to fitting

__device__ __host__ 
int* genRandIntArray ( int array[], const int &N )
{
 // init random seed 
   struct timeval tv;
   gettimeofday( &tv, NULL);
   int usec = tv.tv_usec ; 
   srand48(usec);

   printf("%d", N);
     
   int *a ;
   a = malloc( N * sizeof( int ) ) ;
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
   return a;
}


int main() {
    
    // generate  permutation arrays
    const int Nsmp = 5;
    int *X = malloc( sizeof (int) * Nsmp);
    X = genRandIntArray( X, Nsmp );
    free(X); 
    return 0;
   }

