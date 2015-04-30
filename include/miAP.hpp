#ifndef miAP_HPP
#define miAP_HPP
#include "Cube.hpp"
#include "util.hpp"
#include <cuda.h>
#include <cstdio>
#include <cstdlib>
#define CHI_THRESHOLD 7.815 // 95%

/* 
 * Chi Square statistical test
 * if all four quards together add up to less than 8, subdivision stops and this function returns false.
 * for the first iteration, full size square should always return true
 *
 */
__device__ 
bool chiSeq(unsigned int nSamples, unsigned int a, unsigned int b, unsigned int c, unsigned int d);

/* 
 * called by experiment miAP()
 *
 */
__global__ 
void computeMi(unsigned int *d_rankMatrix, unsigned int nTFs, unsigned int nGenes, unsigned int nSamples, 
        unsigned int *d_TFGenesIdx, float *d_rawGraph, float miThreshold);

/* 
 * called by null model miAP()
 *
 */
__global__
void computeMi(unsigned int *d_randomMatrix, unsigned int nPairs, unsigned int nSamples, float *d_miResult);

/* 
 * this function computes the mutual information given the ranked experimental matrix
 * X and Y are mapped faciliated by d_TFGenesIdx (see details in source code)
 * function called during actual MI computation 
 *
 */
void miAP(unsigned int *d_rankMatrix, unsigned int nTFs, unsigned int nGenes, unsigned int nSamples, 
        unsigned int *d_TFGeneIdx, float **d_rawGraph, float miThreshold);

/* 
 * this function computes the mutual information given random matrix
 * the first half of the matrix is regarded as X and the second half Y
 * function called during null model computation
 * input d_randomMatix should be 2*nPairs by nSamples matrix
 * d_miResult, just an address of a pointer to T is fine
 *
 */
void miAP(unsigned int *d_randomMatrix, unsigned int nPairs, unsigned int nSamples, float **d_miResult);
   
#endif 
