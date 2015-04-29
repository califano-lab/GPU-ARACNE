#ifndef PRUNEGRAPH_HPP
#define PRUNEGRAPH_HPP
#include <cuda.h>
#include "util.hpp"

/*
 *
 * the caller launches nTFs * nTFs * nGenes threads
 * the kernel is implemented here
 * first dimension: TF1
 * second dimension: TF2
 * third dimension: all genes
 * decision array should be initialized to false (0)
 *
 */
__global__
void label(float *d_rawGraph, unsigned int nRows, unsigned int nCols, int *d_decision, unsigned int *d_TFGeneIdx);


/*
 *  simple cut() function: set the edge value to be -1 when the edge is cut
 */
__global__
void cut(float *d_rawGraph, unsigned int nRows, unsigned int nCols, int *d_decision);


/*
 * the following function sequentially calls two kernels label() and cut() to label the edges to be delelted
 * and then cut the edges to be deleted
 *
 */
__host__
void pruneGraph(float *d_rawGraph, unsigned int nRows, unsigned int nCols, unsigned int *d_TFGeneIdx);

#endif
