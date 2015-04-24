// basic viable unit of execution.
//#define TEST
#include "InputOutput.hpp"
#include "Matrix.hpp"
#include "pruneGraph.hpp"
#include "miAP.hpp"
#include "genMiCutoff.hpp"
#include <cstdlib>
#include <cstdio>

// argument list: 
// <TFfile> <datafile> <nTFs> <nGenes> <nSamples> <nBootstraps>
int main(int argc, char *argv[])
{
    // argument check
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] 
            << " <TFfile> <datafile> <nTFs> <nGenes> <nSamples> <nBootstraps> <pValue>" << std::endl;  
        exit(1);
    } 
    char *TFFilename = argv[1];
    char *dataFilename = argv[2];
    unsigned int nTFs = atoi(argv[3]);
    unsigned int nGenes = atoi(argv[4]);
    unsigned int nSamples = atoi(argv[5]);
    unsigned int nBootstraps = atoi(argv[6]);
    float pValue = atof(argv[7]);

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
    Matrix<float> *h_ranked = new Matrix<float>(nGenes, nSamples);
    cudaMemcpy((void *)h_ranked->memAddr(), (void *)d_rankMat, h_ranked->size(), cudaMemcpyDeviceToHost);
    h_ranked->print();
    delete h_ranked;
#endif    
    delete dataMat;

    // calculate MIcutoff
    // at this pint d_rankMat is a nGenes * nSamples matrix with rank
    // this array is already in the GPU
    
    unsigned int seed = 1;
//    float miThreshold = computeMiThreshold(nSamples, pValue, seed);
    float miThreshold = computeMiThreshold(1000, 0.00000001, seed);
    std::cout << miThreshold << std::endl;
    // build network
    // the output of this part should be nTFs * nGenes matrix stored in a plain 1-D array 
    float *d_miValue;
    miAP(d_rankMat, nTFs, nGenes, nSamples, d_TFGeneIdx, &d_miValue, miThreshold);

#ifdef TEST
    Matrix<float> *h_miValue = new Matrix<float>(nTFs, nGenes);
    cudaMemcpy((void *)h_miValue->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost);
    h_miValue->print();
    delete h_miValue;
#endif
    
    // DPI to prune network
    pruneGraph(d_miValue, nTFs, nGenes, d_TFGeneIdx);

#ifdef TEST
    Matrix<float> *h_miValue_pruned = new Matrix<float>(nTFs, nGenes);
    cudaMemcpy((void *)h_miValue_pruned->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost);
    h_miValue_pruned->print();
    delete h_miValue_pruned;
#endif
    // output data

    // cleanup
    delete TFList;
    delete geneLabels;
    cudaFree(d_rankMat);
    cudaFree(d_TFGeneIdx);
    cudaFree(d_miValue);
    return 0;
}
