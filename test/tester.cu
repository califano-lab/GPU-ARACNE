// basic viable unit of execution.
#define TEST
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
    std::vector<std::string> TFList(nTFs);
    loadMatrix(NULL, &TFList, TFFilename, nTFs, 1);
    // import data
    Matrix<float> *dataMat;
    std::vector<std::string> geneLabels(nGenes);
    loadMatrix(&dataMat, &geneLabels, dataFilename, nGenes, nSamples);
    // create TF to genes index mapping
    unsigned int *d_TFGeneIdx;
    createMapping(&d_TFGeneIdx, &TFList, &geneLabels, nTFs, nGenes);
    // rank data
    unsigned int *d_rankMat = dataMat->getRankMatrix();

#ifdef TEST
    for (int i = 0; i < TFList.size(); i++){
        std::cout<< TFList[i] << std::endl;
    }
    for (int i = 0; i < geneLabels.size(); i++){
        std::cout << geneLabels[i] << std::endl;
    }
    dataMat->print();
    Matrix<unsigned int> *h_ranked = new Matrix<unsigned int>(nGenes, nSamples);
    cudaMemcpy((void *)h_ranked->memAddr(), (void *)d_rankMat, h_ranked->size(), cudaMemcpyDeviceToHost);
    h_ranked->print();
    delete h_ranked;
    unsigned int *h_TFGeneIdx = new unsigned int [nTFs];
    cudaMemcpy((void *)h_TFGeneIdx, (void *)d_TFGeneIdx, sizeof(unsigned int) * nTFs, cudaMemcpyDeviceToHost);
    std::cout << "TFGenesIdx: " << std::endl;
    for (int i = 0; i < nTFs; i++){
        std::cout<< h_TFGeneIdx[i] << std::endl;
    }
#endif    
    delete dataMat;

    // calculate MIcutoff
    // at this pint d_rankMat is a nGenes * nSamples matrix with rank
    // this array is already in the GPU
    
    unsigned int seed = 1;
    float miThreshold = computeMiThreshold(nSamples, pValue, seed);
//    float miThreshold = computeMiThreshold(1000, 0.00000001, seed);
#ifdef TEST
    std::cout << miThreshold << std::endl;
#endif
    // build network
    // the output of this part should be nTFs * nGenes matrix stored in a plain 1-D array 
    float *d_miValue;
    miAP(d_rankMat, nTFs, nGenes, nSamples, d_TFGeneIdx, &d_miValue, (float)0.12);

#ifdef TEST
    Matrix<float> *h_miValue = new Matrix<float>(nTFs, nGenes);
    cudaMemcpy((void *)h_miValue->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost);
    HANDLE_ERROR(cudaDeviceSynchronize());
    h_miValue->print();
    delete h_miValue;
#endif
    
    // DPI to prune network
    pruneGraph(d_miValue, nTFs, nGenes, d_TFGeneIdx);

#ifdef TEST
    Matrix<float> *h_miValue_pruned = new Matrix<float>(nTFs, nGenes);
    cudaMemcpy((void *)h_miValue_pruned->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost);
    HANDLE_ERROR(cudaDeviceSynchronize());
    h_miValue_pruned->print();
    delete h_miValue_pruned;
#endif
    // output data

    // cleanup
    cudaFree(d_rankMat);
    cudaFree(d_TFGeneIdx);
    cudaFree(d_miValue);
    return 0;
}
