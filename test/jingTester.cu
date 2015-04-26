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
    printf( "Step 1 miThreshold ...  ") ;
    float miThreshold = computeMiThreshold(nSamples, pValue, seed);
    printf( " Done %f \n", miThreshold) ;

    // 
    // do  bootstraps on  orignal matrix
    // 

    printf( "Step 2 Bootstrapping ... "); 
    unsigned int *d_TFGeneIdx;
    createMapping(&d_TFGeneIdx, TFList, geneLabels, nTFs, nGenes);

    Matrix<float> *d_bsMat  = new Matrix<float>(nGenes, nSamples);

    d_bsMat = dataMat->bootstrapMatrix();
    printf( "Done\n"); 

    Matrix<float> *h_ranked = new Matrix<float>(nGenes, nSamples);
    float *d_rankMat;
    float *d_miValue;
    float *d_bsMiValue;
    int *d_bsMiCount;
    nBootstraps = 1;
    for (int ibs = 0; ibs < nBootstraps; ibs++ ) {

      d_bsMat = dataMat->bootstrapMatrix();
      d_rankMat = d_bsMat->getRankMatrix();
#ifdef TEST
      printf("Original data -----------------\n");
      //dataMat->printHead();
      printf("Bootstrap data-----------------\n");
      //d_bsMat->printHead();
      printf("bs rank data-----------------\n");
      cudaMemcpy((void *)h_ranked->memAddr(), (void *)d_rankMat, h_ranked->size(), cudaMemcpyDeviceToHost);
      HANDLE_ERROR(cudaDeviceSynchronize());
      h_ranked->printHead();
#endif

      printf( "AdaP ... "); 
      miAP(d_rankMat, nTFs, nGenes, nSamples, d_TFGeneIdx, &d_miValue, miThreshold);
      printf( " \n"); 
      // DPI 

      printf( "DPI ..."); 
      pruneGraph(d_miValue, nTFs, nGenes, d_TFGeneIdx);
      printf( "\n"); 

      // consolidate
      printf( "consolidate ..."); 
      printf( " \n"); 
      //Matrix<float> *h_miValue = new Matrix<float>(nTFs, nGenes);
      //cudaMemcpy((void *)h_miValue->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost);
      //HANDLE_ERROR(cudaDeviceSynchronize());
      //h_miValue->print();
      //delete h_miValue;
      //  Matrix<float> *h_miValue_pruned = new Matrix<float>(nTFs, nGenes);
      // cudaMemcpy((void *)h_miValue_pruned->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost);
      // HANDLE_ERROR(cudaDeviceSynchronize());
      // h_miValue_pruned->print();
      // delete h_miValue_pruned;
     
    }

//
//ship d_bsMiValue and d_bsMiCount back, fit poisson model 
//
    Matrix<float> *h_bsMiValue = new Matrix<float>(nTFs, nGenes);
    Matrix<int>   *h_bsMiCount = new Matrix<int>(nTFs, nGenes);
    
    cudaMemcpy((void *)h_bsMiValue->memAddr(), (void *)d_bsMiValue, h_bsMiValue->size(), cudaMemcpyDeviceToHost);
    cudaMemcpy((void *)h_bsMiCount->memAddr(), (void *)d_bsMiCount, h_bsMiCount->size(), cudaMemcpyDeviceToHost);
    HANDLE_ERROR(cudaDeviceSynchronize());
    delete h_bsMiValue;
    delete h_bsMiCount;
     
    delete h_ranked;
    delete dataMat;

    // cleanup
    delete TFList;
    delete geneLabels;
    cudaFree(d_rankMat);
    cudaFree(d_TFGeneIdx);
    cudaFree(d_miValue);
    return 0;
}
