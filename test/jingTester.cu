#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Matrix.hpp"
#include "genMiCutoff.hpp"
#include "genMiCutoff.h"
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2) {
        perror("sample number and permutation number!");
        exit(1);
    }
    int nSamples = atoi(argv[1]);
    int nperm = atoi(argv[2]); // this will change to fix value after testing  

    srand(time(NULL));

    Matrix<float> *h_mat = new Matrix<float>(nGenes, nSamples);
    for (int i = 1; i < nGenes; i++){
        for (int j = 0; j < nSamples; j++){
            h_mat->element(i,j) = (rand() % 100) / 20.0;
        }
    }
    h_mat->print();

    cout << "start calculating null model mi" << endl;
    clock_t start = clock();

    float h_micut = getRankMatrix(10,2); 

    clock_t stop = clock();
    
    return 0;
}
