// basic viable unit of execution.
#define TEST
#include "InputOutput.hpp"
#include "Matrix.hpp"
#include "pruneGraph.hpp"
#include "sumBooststrapGraphs.hpp"
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
    // TODO: add stream control to all  bootstraps 
    // 

    printf( "Step 2 Bootstrapping ... "); 
    unsigned int *d_TFGeneIdx;
    createMapping(&d_TFGeneIdx, TFList, geneLabels, nTFs, nGenes);

    Matrix<float> *d_bsMat  = new Matrix<float>(nGenes, nSamples);
    Matrix<float> *h_ranked = new Matrix<float>(nGenes, nSamples);

    nBootstraps = 2; // for test 
    float *d_rankMat;
    float *d_miValue;
    float *d_bsMiValue; 
    int   *d_bsMiCount;

    dim3 bDim(32, 32);
    dim3 gDim( ceil(nGenes /(1.0 * bDim.x) ), ceil(nSamples / (1.0 * bDim.y) ) ) ;

    for (int ibs = 0; ibs < nBootstraps; ibs++ ) {

      printf("Bootstrap %d .... ", ibs);
      d_bsMat = dataMat->bootstrapMatrix();
      d_rankMat = d_bsMat->getRankMatrix();
      printf(" Done\n");
#ifdef TEST
      //printf("Original data -----------------\n");
      //dataMat->printHead();
      //printf("Bootstrap data-----------------\n");
      //d_bsMat->printHead();
      //printf("bs rank data-----------------\n");
      //cudaMemcpy((void *)h_ranked->memAddr(), (void *)d_rankMat, h_ranked->size(), cudaMemcpyDeviceToHost);
      //HANDLE_ERROR(cudaDeviceSynchronize());
      //h_ranked->printHead();
#endif

      printf( "AdaP ... "); 
      miAP(d_rankMat, nTFs, nGenes, nSamples, d_TFGeneIdx, &d_miValue, miThreshold);
      printf( " \n"); 

      // DPI 
      printf( "DPI ..."); 
      pruneGraph(d_miValue, nTFs, nGenes, d_TFGeneIdx);
      printf( "\n"); 
      
      //TODO: ship data from next stream

      // merge bs graphs 
      printf( "merge bs graphs ..."); 
      sumBooststrapGraphs <<< bDim, gDim >>>(d_bsMiValue, *d_bsMiCount, *d_miValue, nGenes, nSamples);
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
    printf( "Count edges and occurences of all bs graphs ..."); 

    Matrix<float> *h_bsMiValue = new Matrix<float>(nTFs, nGenes);
    Matrix<int>   *h_bsMiCount = new Matrix<int>(nTFs, nGenes);
    
    cudaMemcpy((void *)h_bsMiValue->memAddr(), (void *)d_bsMiValue, h_bsMiValue->size(), cudaMemcpyDeviceToHost);
    cudaMemcpy((void *)h_bsMiCount->memAddr(), (void *)d_bsMiCount, h_bsMiCount->size(), cudaMemcpyDeviceToHost);
    HANDLE_ERROR(cudaDeviceSynchronize());

    calMean4Poisson_d<<< bDim, gDim>>>( T *d_sumMiGraph, T *d_sumCountGraph, unsigned int nRows, unsigned int nCols, long totalEdges, long totalOccurence)

    printf( "\n");

    //
    // integrating bootstrapping mi and count
    //
    // TODO: current output to stdout, add output to file option
    poissonIntegrate( h_bsMiGraph, h_bsCountGraph, totalOccurence, totalEdges,  nRows, nCols, &TFList, &geneLabels)

    //
    // cleanup
    //
    delete h_bsMiValue;
    delete h_bsMiCount;
    delete h_ranked;

    delete dataMat;
    delete TFList;
    delete geneLabels;
    cudaFree(d_rankMat);
    cudaFree(d_TFGeneIdx);
    cudaFree(d_miValue);
    cudaFree(d_bsMiValue);
    cudaFree(d_bsMiCount);
    return 0;
}
