// basic viable unit of execution.
#include "InputOutput.hpp"
#include "Matrix.hpp"
#include "genMiCutoff.hpp"
#include <cstdlib>
#include <cstdio>

// argument list: 
// < datafile > < nGenes > < nSamples > < nBootstraps > < pValueThreshold > 
int main(int argc, char *argv[])
{
    // argument check
    if (argc != 8) {
        std::cerr << "Usage: " << argv[0] 
            << " <TFfile> <datafile> <nTFs> <nGenes> <nSamples> <nBootstraps> <pvalue>" << std::endl;  
        exit(1);
    } 
    char *TFFilename = argv[1];
    char *dataFilename = argv[2];
    unsigned int nTFs = atoi(argv[3]);
    unsigned int nGenes = atoi(argv[4]);
    unsigned int nSamples = atoi(argv[5]);
    unsigned int nBootstraps = atoi(argv[6]);
    long pvalue = atol( argv[7] ); 

    // import transcription factor list
    Matrix<std::string> *TFList;
    loadMatrix(NULL, &TFList, TFFilename, nTFs, 1);
    
    // import data
    Matrix<float> *dataMat;
    Matrix<std::string> *geneLabels;
    loadMatrix(&dataMat, &geneLabels, dataFilename, nGenes, nSamples);
    
    // rank data
    
    // calculate MIcutoff
    std::srand( unsigned ( std::time(0) ) );	
    const int nSample = 20;
    const int nSubsmp = 10;
    int *tempRandIdxArray; // local random subsampling index
    tempRandIdxArray = ( int * ) malloc( sizeof(int) * nSample );
    genRandIntArray_h( tempRandIdxArray, nSample, nSubsmp );

    // at this pint d_rankMat is a nGenes * nSamples matrix with rank
    // this array is already in the GPU


    // build network
    
    // dpi to prune network
    
    // output data
    return 0;
}
