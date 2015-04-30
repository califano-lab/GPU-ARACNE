#ifndef INPUTOUTPUT_HPP
#define INPUTOUTPUT_HPP

#include"Matrix.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <tr1/unordered_map>
#include "util.hpp"

/* 
 * load input data into Matrix.
 * if mat == NULL, then only the first column is loaded into the matrix
 * for example, the transcription factor list may be loaded in this way
 * nRows: number of genes usually
 * nCols: number of samples usually
 *
 */
void loadMatrix(Matrix<float> **mat, std::vector<std::string> *geneLabels,
        const char* filename, unsigned int nRows, unsigned int nCols);

/*
 * create mapping between transcription factor index to general gene index
 *
 */
    
void createMapping(unsigned int **d_TFGeneIdx, std::vector<std::string> *TFList, 
        std::vector<std::string> *geneLabels, unsigned int nTFs, unsigned int nGenes);
#endif
