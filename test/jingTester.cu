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
// <TFfile> <datafile> <nTFs> <nGenes> <nSamples> <nBootstraps> <pValue>
int main(int argc, char *argv[])
{
    //
    // argument check
    //
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

    // 
    // import data  
    //
    Matrix<std::string> *TFList;
    loadMatrix(NULL, &TFList, TFFilename, nTFs, 1);
    
    Matrix<float> *dataMat;
    Matrix<std::string> *geneLabels;
    loadMatrix(&dataMat, &geneLabels, dataFilename, nGenes, nSamples);
    
    // 
    // calculate MIcutoff
    // 
    unsigned int seed = 1;
    float miThreshold = computeMiThreshold(nSamples, pValue, seed);
    printf( "Finish miThreshold ... %f \n ", miThreshold) ;


    // 
    // do  bootstraps on  orignal matrix
    // 
    unsigned int *d_TFGeneIdx;
    createMapping(&d_TFGeneIdx, TFList, geneLabels, nTFs, nGenes);

    float *d_bsMat;
    float *d_rankMat;
    float *d_miValue;
    for (ibs = 0; i < nBootstraps; i ++ ) {

      d_bsMat = dataMat->bootstrapMatrix();
      d_rankMat = d_bsMat->getRankMatrix();

      Matrix<float> *h_ranked = new Matrix<float>(nGenes, nSamples);
      cudaMemcpy((void *)h_ranked->memAddr(), (void *)d_rankMat, h_ranked->size(), cudaMemcpyDeviceToHost);

      miAP(d_rankMat, nTFs, nGenes, nSamples, d_TFGeneIdx, &d_miValue, miThreshold);

      Matrix<float> *h_miValue = new Matrix<float>(nTFs, nGenes);
      cudaMemcpy((void *)h_miValue->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost);
      HANDLE_ERROR(cudaDeviceSynchronize());
      h_miValue->print();
      delete h_miValue;
#end  if
      
      // DPI to prune network
      pruneGraph(d_miValue, nTFs, nGenes, d_TFGeneIdx);

    }

delete h_ranked;
delete dataMat;
#ifdef TEST
    Matrix<float> *h_miValue_pruned = new Matrix<float>(nTFs, nGenes);
    cudaMemcpy((void *)h_miValue_pruned->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost);
    HANDLE_ERROR(cudaDeviceSynchronize());
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
