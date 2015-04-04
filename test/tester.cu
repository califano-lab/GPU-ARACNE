#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Matrix.hpp"
using namespace std;
int main(int argc, char *argv[])
{
    int nGenes = atoi(argv[1]);
    int nSamples = atoi(argv[2]);
    srand(time(NULL));
    Matrix *mat = new Matrix(nGenes, nSamples);
    for (int i = 1; i < nGenes; i++){
        for (int j = 0; j < nSamples; j++){
            mat->element(i,j) = (rand() % 100) / 20.0;
        }
    }
    cout << "start ranking" << endl;
    //mat->print();
    clock_t start = clock();
    Matrix *d_mat = mat->getRankMatrix();
    clock_t stop = clock();
    //mat->print();
    cout << "done ranking" << endl;
    cout << "Time taken: " << (float)(stop - start)/CLOCKS_PER_SEC << endl;
    delete(mat);
    delete(d_mat);
    
    return 0;
}
