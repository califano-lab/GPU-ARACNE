//
//  genMiCutoff.cpp
//  aracneGPU
//

#include <iostream>
#include <stdlib.h> 
#include <stdio.h>
#include <time.h> 
#include <math.h> 
#include <conio.h>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

void genMIcutoff( float *dMat, float pValue, float *miCutoff ) 
{
  // host function to control miCutoff calculation 

  // 1. git sample size and ngene
  const int Ngene = dMat.nrow();
  const int Nsmp = dMat.ncol();
  const int NsubSmp = 15;

  // 2. init sample proption
  int subSmp[ NsubSmp ]; 
  for ( int i = 1; i < Nsubsmp; ++i )
  {
    subSmp[i-1] = (int) (0.3 + (i - 1) * 0.05) * Nsmp  ; 
  }

  const int numRep = 3;
  Nperm = 100000; 
  
  // init array to hold result  on host
  float *alpha_h; 
  *alpha_h = (float *) malloc( sizeof(float) * NsubSmp * numRep ) ; 

  float *beta_h;
  *beta_h = (float *) malloc( sizeof(float) * NsubSmp * numRep ) ; 

  // add some random nosie to data  kernel, leave data on kernel 

  // copy alpha_h, beta_h to device 
  // do permutation for each subsampling 
  float *dmat_sub_d;
  cudamalloc( *dmat_sub_d, sizeof(dMat_h) ) ; // use the maximum size   

  // for each sub sample and each replicate calculate threshold and coefficient 
  // run for 14 * 3 times 
  for ( int isubSmp = 0; isubSmp < NsubSmp ; ++isubSmp ) {
     for ( int inumRep  = 0; inumRep < numRep ; ++ inumRep ) {
        // take a subset of dMat,  
        // and generate a similiar size of array with random permuated integer 
        genRandIntArray_h( *dmat_sub_h, *subSmp[isubSmp] ) ; 

        // overload the dmat_sub_d pointer every time launch a size >= previous_size
        cudaMemcpy( dmat_sub_d, dmat_sub_h, 0, Host2Device );      

        // launch kernel 
	extraplotaMICutoff_d <<< ( 1,1,1), (1,1,1) >>> ( *dmat_sub_d, Nperm, *coef_d ); 

	// copy result back
	__cudasyncDevice(); 
	cudaMemcpy( coef_h, coef_d, Device2Host) ;
   }
  }

  // linear fitting 
    float coef[2]; 
    coef = subSmp ~ mean(beta_h) ;  // linear regression 

    float a = mean( alpha_h ) ; // overall alpha 
    float b = coef[0]; 
    float c = coef[1]; 

  // calcuate MIcutoff  
    miCutoff = ( log( pValue )  - a ) / ( b* n + c ) ;

    return; 

}

// local function 
void extrapolateMICutoff_d( float *dMat, const int Nperm, float *nullMI ) 
{
  // input dMat is random permuted integers with Nsmp_sub <= Nsmp Columns 
  // and Ngene rows;
  const int Ngene = dMat.getNumRows() ;
  const int Nsmp = dMat.getNumColumns() ;
  // random take 2 index 

  int randomI; 
  int randomJ; 

  // prepare random datadata 
  X = dMat[ randomI ]; 
  Y = dMat[ randomJ ];


}

__global__
void genCoef_subSmp ( float *nullMIArray , float *coef, const int Nperm, float pPvalue ) 
{
  // run use one block 
  // sort 
  const int N = Nperm;  
  thrust::stable_sort_by_key(thrust::seq , nullMIArray , nullMIArray + Nperm);
  

  // generate density
  const int Npoint = Nperm / 100 ; 
  int tdx = threadIdx.x + blockIdx.x * blockDim.x;
 
  // init cumulative  
  float cumf1[Npoint] ;   
  float cumf2[Npoint] ;   
  bool cumf2Bool[Npoint] ;   
  if (tdx < Npoint )  {
      cumf1[ tdx + i ] = 0.0 ; 
      cumf2[ tdx + i ] = 0.0 ; 
      cumf2Bool[ tdx + i ] = true ; 
    }
  }
  __syncthreads();

  // take Npoint  
  int cumf2BoolSum = 0; 
  if( tdx < Npoint ) {
    cumf1[ tdx ] = nullMIArray[ tdx * 100  ];
    float temp = (Nperm - tdx * 100) /(1.0 * Nperm); 
    cumf2[ tdx ] = temp ;
    bool tempB = (temp == 0 ); 
    cumf2Bool[ tdx ] =  tempB ;  
    cumf2BoolSum += (int) tempB ;   
  }
  __syncthreads();

  const int Nfit = 100;
  if( tdx < Npoint ) {
    cumf1[ tdx ]  = nullMIArray[ tdx * 100 ]; 
    cumf2[ tdx ]  = nullMIArray[ tdx * 100 ];  
         
  }
  
}

void genRandIntArray_h ( int *a, const int N )
{
   // helper function to generate random permuation of N samples init random seed 

   printf("Total: %d\n", N);
     
   // int *a ;
   // a = ( int * ) malloc( N * sizeof( int ) ) ;
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
   // for ( int i = 0; i < N; ++i ) 
   //{
    // printf("%d\t",  a[i]) ;
  // }
   //printf( "\n");
   //return a;
   return;
}
