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

// main function to  miCutoff calculation 
__host__
float genMIcutoff( float *dMat, float pValue ) 
{

  // get sample size and ngene
  const int nGene = dMat.nRows();
  const int nSample = dMat.nCols();
  const int nRep = 3;
  const int nPerm = 100000; 

  // init sample propotion
  const int nSubsmp = 15;
  int subSmp[ nSubsmp ]; 
  for ( int i = 1; i <= nSubsmp; ++i )
  {
    subSmp[i-1] = (int) (0.3 + (i - 1) * 0.05) * nSample  ; 
  }

  // < OPTIONAL > add some random nosie to data  kernel, leave data on kernel 

  // init array to hold result  on host
  // copy alpha_h, beta_h to device 

  // for each subsmp size, calculate alpha, beta on GPU independently 

  
  // do permutation for each subsampling 
  // for each sub sample and each replicate calculate threshold and coefficient 
  // init for one loop
  Matrix<T> *dMatSub_h = new Matrix<T>( nGene, nSubsmp ) ;
  Matrix<T> *rMatSub_h = new Matrix<T>( nGene, nSubsmp ) ;

  int *rMatSub_d; 
  HANDLE_ERROR( cudaMalloc( (void **) &rMatSub_d, sizeof(int) * nSample * nGene) );

  
  // define block dim, grid dim 
  dim3 blockSize( 8, 8, 16 );
  dim3 gridSize( ceil(nRows / 8.0), ceil(nRows / 8.0), ceil(nCols / 16.0) );
  int *tempRandIdxArray; // local random subsampling index
  tempRandIdxArray = ( int * ) malloc( sizeof(int) * nSample );


  // run for 14 * 3 times, consider to use streams later , TODO
  float *coef_h;  
  *coef_h  = (float *) malloc( sizeof(float) * nSubsmp * numRep * 2 ) ; 
  float *coef_d ;
  HANDLE_ERROR( cudaMalloc( (void **) &coef_d,   sizeof(float) * 2) );
  HANDLE_ERROR( cudaMemset( (void *)coef_d, 1.0, sizeof(float) * 2) );

  // HANDLE_ERROR( cudaMalloc( (void **) &alpha_d, sizeof(float) * nSubsmp * nRep) );

  for ( int iSubsmp = 0; iSubsmp < nSubsmp ; ++iSubsmp ) {
     for ( int iRep  = 0; iRep < nRep ; ++iRep ) {
        // take a subset of dMat
        // overload the dmat_sub_d pointer every time launch a size >= previous_size
	 genRandIntArray_h( tempRandIdxArray, nSample, nSubsmp );
	 dMatSub_h = dMat.getCols( tempRandIdxArray );

        // and generate a similiar size of array with random permuated integer 
         rMatSub_h = dMatSub_h.genRandPermMatrix(); 
	 HANDLE_ERROR( cudaMemcpy( (void *) dMatSub_d, (void *) dMatSub_h->memAddr(), dMatSub_h->size(), cudaMemcpyHostToDevice)  ) ; 
	
        // launch kernel 
	 extraplotaMICutoff_d <<< gridSize, blockSize >>> (dMatSub_d, nPerm, coef_d );

	 HANDLE_ERROR(cudaGetLastError());
    	 HANDLE_ERROR(cudaDeviceSynchronize());
	 __cudasyncDevice(); 
	// copy result back
	 HANDLE_ERROR( cudaMemcpy( coef_h[ iSubsmp * nSubsmp + iRep ] , coef_d, sizeof(float) * 2, cudaMemcpyDeviceToHost) );
	
   }
  }

  // linear fitting 
    float coef[2]; 
    coef[0] = 1.0;
    coef[1] = 1.0;
    //coef = subSmp ~ mean(beta_h) ;  // linear regression 

    float a = mean( alpha_h ) ; // overall alpha 
    float b = coef[0]; 
    float c = coef[1]; 

  // calcuate MIcutoff  
    float miCutoff; 
    miCutoff = ( log( pValue )  - a ) / ( b* n + c ) ;

    return miCutoff; 

}

// local function 
// void extrapolateMICutoff_d( float *dMat, const int Nperm, float *nullMI ) 
// {
//   // input dMat is random permuted integers with nSubsmp <= nSample Columns 
//   // and nGene rows;
//   const int nGene = dMat.getNumRows() ;
//   const int nSample = dMat.getNumColumns() ;
//   // random take 2 index 
// 
//   int randomI; 
//   int randomJ; 
// 
//   // prepare random datadata 
//   X = dMat[ randomI ]; 
//   Y = dMat[ randomJ ];
// 
// }

//__global__
//void genCoef_subSmp ( float *nullMIArray , float *coef, const int Nperm, float pPvalue ) 
//{
//  // run use one block 
//  // sort 
//  const int N = Nperm;  
//  thrust::stable_sort_by_key(thrust::seq , nullMIArray , nullMIArray + Nperm);
//  
//
//  // generate density
//  const int Npoint = Nperm / 100 ; 
//  int tdx = threadIdx.x + blockIdx.x * blockDim.x;
// 
//  // init cumulative  
//  float cumf1[Npoint] ;   
//  float cumf2[Npoint] ;   
//  bool cumf2Bool[Npoint] ;   
//  if (tdx < Npoint )  {
//      cumf1[ tdx + i ] = 0.0 ; 
//      cumf2[ tdx + i ] = 0.0 ; 
//      cumf2Bool[ tdx + i ] = true ; 
//    }
//  }
//  __syncthreads();
//
//  // take Npoint  
//  int cumf2BoolSum = 0; 
//  if( tdx < Npoint ) {
//    cumf1[ tdx ] = nullMIArray[ tdx * 100  ];
//    float temp = (Nperm - tdx * 100) /(1.0 * Nperm); 
//    cumf2[ tdx ] = temp ;
//    bool tempB = (temp == 0 ); 
//    cumf2Bool[ tdx ] =  tempB ;  
//    cumf2BoolSum += (int) tempB ;   
//  }
//  __syncthreads();
//
//  const int Nfit = 100;
//  if( tdx < Npoint ) {
//    cumf1[ tdx ]  = nullMIArray[ tdx * 100 ]; 
//    cumf2[ tdx ]  = nullMIArray[ tdx * 100 ];  
//         
//  }
//  
//}
//


