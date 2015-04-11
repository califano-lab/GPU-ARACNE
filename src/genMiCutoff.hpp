//  Created by jing on 3/24/15.
//  Copyright (c) 2015 jing. All rights reserved.

// #ifndef __aracneGPU__genMiCutoff__
// #define __aracneGPU__genMiCutoff__

__host__
int *genRandIntArray_h ( int *b, const int N , const int M )
{
   // helper function to generate random permuation of N samples init random seed 
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

   for ( int i = 0; i < M; ++i ) 
   {
     b[i] = a[i];
   }

   return b ;
}


// #endif /* defined(__aracneGPU__genMiCutoff__) */
