#ifndef SUMBOOSTSTRAPGRAPHS_H
#define SUMBOOSTSTRAPGRAPHS_H
#include <cuda.h>
#include <random>
#include "util.hpp"
#define POISSONCUT 0.05
// the caller launches 2000 * 2000 * 20000 threads
// the kernel is implemented here
// first dimension: TF1
// second dimension: TF2
// third dimension: all genes
// decision array should be initialized to true


template<typename T>
__global__
void sumBooststrapGraphs(T *d_sumMiGraph, int *d_sumCountGraph, T *d_currGraph, unsigned int nRows, unsigned int nCols)
{
  // sum the bootstrapping mi value, and edge information 
    unsigned int rowIdx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int colIdx = blockIdx.y * blockDim.y + threadIdx.y;
    if (rowIdx >= nRows || colIdx >= nCols) return;
    if (!d_currGraph[rowIdx * nCols + colIdx] < 0)
    { 
        d_sumCountGraph[rowIdx * nCols + colIdx] += 1;
        d_sumMiGraph[rowIdx * nCols + colIdx] += d_currGraph[rowIdx * nCols + colIdx];
    }
}

__global__
void calMean4Poisson_d( T *d_sumMiGraph, T *d_sumCountGraph, unsigned int nRows, unsigned int nCols, long totalEdges, long totalOccurence)
{
  // to get mean mi for each edge, average edges of all bootstraping
  // call after all bootstraps finished 
    unsigned int rowIdx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int colIdx = blockIdx.y * blockDim.y + threadIdx.y;

    if (rowIdx >= nRows || colIdx >= nCols) return;

    long matIdx = rowIdx * nCols + colIdx ;
    if (!d_currGraph[ matIdx ] < 0)
    { 
        d_sumMiGraph[ matIdx ] = d_sumMiGraph[ matIdx ] / d_sumCountGraph[ matIdx ];
        atomicAdd(totalOccurence, d_sumCountGraph[ matIdx ] ); //TODO: think of better ways to to this ;
        atomicAdd(totalEdge, 1 );
    }

}

template<typename T>
__host__
void poissonIntegrate(T *bsMiGraph, int *bsCountGraph, long totalOccurence, long totalEdges,  unsigned int nRows, unsigned int nCols, Matrix<string> *TFList, Matrix<string> *geneLabels)
{
   // host function used to calculate poisson model for 100 bootstraps  
   // default poisson pvalue set to 0.05
  
   // init poisson cdf library;
   float meanEdges = totalOccurence / ( 1.0 * totalEdges );
   map <int, float > *pcdf;
   poissonLibrary( meanEdges, totalOccurence, &pcdf);

   float mi; int id; 
   long tempOccurence; 
   float tempPvalue;  
   printf("TF\tGene\tMI\tpoissonPvalue\n");
   for (int i = 0; i < nRows; ++i )
   {
     for (int j = 0; i< nCols; ++j )
     {
       id = i * nCols + j;
       mi = bsMiGraph->element(i , j );
       tempOccurence = bsCountGraph->element( i, j ) ; 
       tempPvalue = pcdf.find(tempOccurence)->second;
       if( tempPvalue > POISSONCUT )
       { 
         printf("%s\t%s\t%f\t%f\n", TFList->element(i,0) ,geneLabels->element(j,0) , mi, tempPvalue  );
       } 
     }
   }
   printf("\n### ----------  YEAH ----------- #####\n");
   printf("### ---- YOU REACH THE END ----- #####\n");
   printf("### --------- GO TO BED -------- #####\n");
}

__host__ 
void poissonLibrary ( float mean, long totalOccurence, map <int, float> *pcdf ) 
{
  // generate a library of poisson mean <-> P(X > k) probability  mapping 
  // return a map object, used as pcdf.find(X)->second to get a float p-value
    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<> d( mean );
    std::map<int, int> hist;
    for( int n = 0; n < totalOccurence + 1; ++n ) 
    {
      ++hist[d(gen)];
    }

    float temp;
    for (int i = 0; i < hist.size(); i++){
      for ( int j = 0; j <= i; j++) 
      {
	temp += hist.find(i)->second;
      }
      pcdf[i] = 1 - temp / (1.0 * totalOccurence);
    }
}

#endif
