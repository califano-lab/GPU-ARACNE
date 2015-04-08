// basic viable unit of execution.
#define TEST
#include "InputOutput.hpp"
#include "Matrix.hpp"
#include "pruneGraph.hpp"
#include <cstdlib>
#include <cstdio>

// argument list: 
// <TFfile> <datafile> <nTFs> <nGenes> <nSamples> <nBootstraps>
int main(int argc, char *argv[])
{
    // argument check
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0] 
            << " <TFfile> <datafile> <nTFs> <nGenes> <nSamples> <nBootstraps>" << std::endl;  
        exit(1);
    } 
    char *TFFilename = argv[1];
    char *dataFilename = argv[2];
    unsigned int nTFs = atoi(argv[3]);
    unsigned int nGenes = atoi(argv[4]);
    unsigned int nSamples = atoi(argv[5]);
    unsigned int nBootstraps = atoi(argv[6]);

    // import transcription factor list
    Matrix<std::string> *TFList;
    loadMatrix(NULL, &TFList, TFFilename, nTFs, 1);
    
    // import data
    Matrix<float> *dataMat;
    Matrix<std::string> *geneLabels;
    loadMatrix(&dataMat, &geneLabels, dataFilename, nGenes, nSamples);
    
    // create TF to genes index mapping
    unsigned int *d_TFGeneIdx;
    createMapping(&d_TFGeneIdx, TFList, geneLabels, nTFs, nGenes);

    // rank data
    float *d_rankMat = dataMat->getRankMatrix();
#ifdef TEST
    geneLabels->print();
    TFList->print();
    dataMat->print();
#endif    
    delete dataMat;
    // calculate MIcutoff
    // at this pint d_rankMat is a nGenes * nSamples matrix with rank
    // this array is already in the GPU


    // build network
    // the output of this part should be nTFs * nGenes matrix stored in a plain 1-D array
    float *rawGraph = NULL; 
    
    // DPI to prune network
    pruneGraph(rawGraph, nTFs, nGenes, d_TFGeneIdx);
    
    // output data

    // cleanup
    delete TFList;
    delete geneLabels;
    cudaFree(d_rankMat);
    cudaFree(d_TFGeneIdx);
    return 0;
}
