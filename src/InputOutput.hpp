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
using namespace std;

void loadMatrix(Matrix<float> **mat, vector<string> *geneLabels,
        const char* filename, unsigned int nRows, unsigned int nCols)
{
    if (mat != NULL)
        *mat = new Matrix<float>(nRows, nCols);
    //geneLabels = new vector<string>(nRows);
    ifstream fs;
    fs.open(filename);
    if (!fs.is_open()) {
        perror("File open failed.");
        exit(1);
    }
    string line;
    string token;
    int rowIdx = 0;
    int colIdx = 0;
    bool skip = true;
    while(getline(fs, line)) {
        colIdx = 0;
        istringstream lineStream(line);
        skip = true;
        while(getline(lineStream, token, ',')) {
            if (skip){
                (*geneLabels)[rowIdx] = token;
                skip = false;
                continue;
            }
            char *end;
            (*mat)->element(rowIdx, colIdx) = strtof(token.c_str(), &end);
            colIdx++;
        }
        rowIdx++;
    }
    fs.close();
}

void createMapping(unsigned int **d_TFGeneIdx, vector<string> *TFList, 
        vector<string> *geneLabels, unsigned int nTFs, unsigned int nGenes)
{
    unsigned int *h_TFGeneIdx = new unsigned int[nTFs];
    tr1::unordered_map<string, unsigned int> genePool;
    for (unsigned int i = 0; i < nGenes; i++){
        genePool[(*geneLabels)[i]] = i;
    }
    for (unsigned int i = 0; i < nTFs; i++){
        h_TFGeneIdx[i] = genePool[(*TFList)[i]];
    }
    HANDLE_ERROR (cudaMalloc((void **)d_TFGeneIdx, sizeof(unsigned int) * nTFs));
    HANDLE_ERROR (cudaMemcpy((void *)*d_TFGeneIdx, h_TFGeneIdx, sizeof(unsigned int) *nTFs, cudaMemcpyHostToDevice));
    HANDLE_ERROR (cudaDeviceSynchronize());
    delete[] h_TFGeneIdx;
}
#endif
