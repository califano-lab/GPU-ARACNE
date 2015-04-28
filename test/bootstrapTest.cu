// basic viable unit of execution.
//#define TEST
#include "InputOutput.hpp"
#include "Matrix.hpp"
#include "pruneGraph.hpp"
#include "miAP.hpp"
#include "genMiCutoff.hpp"
#include "BootstrapContainer.hpp"
#include <cstdlib>
#include <cstdio>
#define POISSON_P_VALUE 0.05
#define SMALL_LIMIT 0.00001

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
    std::cerr << "### Data Import Completed ###" << std::endl;

    // create TF to genes index mapping
    unsigned int *d_TFGeneIdx;
    createMapping(&d_TFGeneIdx, &TFList, &geneLabels, nTFs, nGenes);
    std::cerr << "### Mapping Created ###" << std::endl;
#ifdef TEST
    for (int i = 0; i < TFList.size(); i++){
        std::cout << TFList[i] << std::endl;
    }
    for (int i = 0; i < geneLabels.size(); i++){
        std::cout << geneLabels[i] << std::endl;
    }
    dataMat->printHead();
    unsigned int *h_TFGeneIdx = new unsigned int [nTFs];
    cudaMemcpy((void *)h_TFGeneIdx, (void *)d_TFGeneIdx, sizeof(unsigned int) * nTFs, cudaMemcpyDeviceToHost);
    std::cout << "TFGenesIdx: " << std::endl;
    for (int i = 0; i < nTFs; i++){
        std::cout<< h_TFGeneIdx[i] << std::endl;
    }
#endif
 
    // calculate MIcutoff
    unsigned int seed = 1;
    std::cerr << "### Constructing Null Model ###" << std::endl;
    float miThreshold = computeMiThreshold(nSamples, pValue, seed);
    std::cerr << "### Null Model Constructed ###" << std::endl;
#ifdef TEST
    std::cout << miThreshold << std::endl;
#endif

    std::cerr << "### Start Bootstrapping ###" << std::endl;
    // preparing containter
    BootstrapContainer bContainer(nTFs, nGenes);

    for (int iteration = 0; iteration < nBootstraps; iteration++){
        // draw random samples with repeats from original sample pool
        Matrix<float> *subsample = dataMat->bootstrapMatrix();

        // rank data
        unsigned int *d_rankMat = subsample->getRankMatrix();
#ifdef TEST
        Matrix<unsigned int> *h_ranked = new Matrix<unsigned int>(nGenes, nSamples);
        cudaMemcpy((void *)h_ranked->memAddr(), (void *)d_rankMat, h_ranked->size(), cudaMemcpyDeviceToHost);
        h_ranked->printHead();
        delete h_ranked;
#endif  
        delete subsample;
        
        // build network
        // the output of this part should be nTFs * nGenes matrix stored in a plain 1-D array 
        float *d_miValue;
        miAP(d_rankMat, nTFs, nGenes, nSamples, d_TFGeneIdx, &d_miValue, miThreshold);
#ifdef TEST
        Matrix<float> *h_miValue = new Matrix<float>(nTFs, nGenes);
        cudaMemcpy((void *)h_miValue->memAddr(), (void *)d_miValue, h_miValue->size(), cudaMemcpyDeviceToHost);
        HANDLE_ERROR(cudaDeviceSynchronize());
        h_miValue->printHead();
        delete h_miValue;
#endif
        HANDLE_ERROR(cudaFree(d_rankMat));

        // Data processing inequality
        pruneGraph(d_miValue, nTFs, nGenes, d_TFGeneIdx);
#ifdef TEST
        Matrix<float> *h_miValue_pruned = new Matrix<float>(nTFs, nGenes);
        cudaMemcpy((void *)h_miValue_pruned->memAddr(), (void *)d_miValue, h_miValue_pruned->size(), cudaMemcpyDeviceToHost);
        HANDLE_ERROR(cudaDeviceSynchronize());
        h_miValue_pruned->printHead();
        delete h_miValue_pruned;
#endif

        // Data aggregation
        bContainer.addToContainer(d_miValue);
        HANDLE_ERROR(cudaFree(d_miValue));
    }
    std::cerr << "### Bootstrapping Finished ###" << std::endl;
    
    std::cerr << "### Condensing Graph ###" << std::endl;
    // condense the graph based on Poisson distribution
    bContainer.condenseGraph(POISSON_P_VALUE);

    // delete the original data matrix after bootstrapping 
    delete dataMat;

    // output data
    std::cerr << "### Writing Output ###" << std::endl;
    Matrix<float> *h_miValue_final = new Matrix<float>(nTFs, nGenes);
    cudaMemcpy((void *)h_miValue_pruned->memAddr(), (void *)bContainer.miValueMemAddr(), 
            h_miValue_final->size(), cudaMemcpyDeviceToHost);
    HANDLE_ERROR(cudaDeviceSynchronize());
    for (int i = 0; i < nTFs; i++){
        for (int j = 0; j < nGenes; j++){
            if (h_miValue_final->element(i, j) > SMALL_LIMIT)
                std::cout << TFList[i] << " " << geneLabels[j] << " " 
                    << h_miValue_final->element(i, j) << std::endl;
        }
    }

    // clean-up
    delete h_miValue_final;
    cudaFree(d_TFGeneIdx);
    std::cerr << "### Done! ###" << std::endl;
    return 0;
}
