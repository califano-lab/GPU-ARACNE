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

void loadMatrix(Matrix<float> **mat, Matrix<string> **geneLabels,
        const char* filename, unsigned int nRows, unsigned int nCols)
{
    if (mat != NULL)
        *mat = new Matrix<float>(nRows, nCols);
    *geneLabels = new Matrix<string>(nRows, 1);
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
                (*geneLabels)->element(rowIdx, 0) = token;
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

void createMapping(unsigned int **d_TFGeneIdx, Matrix<string> *TFList, 
        Matrix<string> *geneLabels, unsigned int nTFs, unsigned int nGenes)
{
    unsigned int *h_TFGeneIdx = new unsigned int[nTFs];
    tr1::unordered_map<string, unsigned int> genePool;
    for (unsigned int i = 0; i < nGenes; i++){
        genePool[geneLabels->element(i, 0)] = i;
    }
    for (unsigned int i = 0; i < nTFs; i++){
        h_TFGeneIdx[i] = genePool[TFList->element(i, 0)];
    }
#ifdef TEST
    for (int i = 0; i < nTFs; i++){
        cout << h_TFGeneIdx[i] << endl;
    }
#endif
    HANDLE_ERROR (cudaMalloc((void **)d_TFGeneIdx, sizeof(unsigned int) * nTFs));
    HANDLE_ERROR (cudaMemcpy((void *)*d_TFGeneIdx, h_TFGeneIdx, sizeof(unsigned int) *nTFs, cudaMemcpyHostToDevice));
    HANDLE_ERROR (cudaDeviceSynchronize());
}
#endif
