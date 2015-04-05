#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Matrix.hpp"
using namespace std;
int main(int argc, char *argv[])
{
    if (argc != 3) {
        perror("need two more arguments!");
        exit(1);
    }
    int nGenes = atoi(argv[1]);
    int nSamples = atoi(argv[2]);
    srand(time(NULL));
    Matrix<float> *h_mat = new Matrix<float>(nGenes, nSamples);
    for (int i = 1; i < nGenes; i++){
        for (int j = 0; j < nSamples; j++){
            h_mat->element(i,j) = (rand() % 100) / 20.0;
        }
    }
    h_mat->print();
    cout << "start ranking" << endl;
    clock_t start = clock();
    float *d_arr = h_mat->getRankMatrix();
    clock_t stop = clock();
    float *h_arr = new float[h_mat->size()];
    cudaMemcpy((void *)h_arr, (void *)d_arr, h_mat->size(), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    Matrix<float> *h_ret = new Matrix<float>(nGenes, nSamples, h_arr);

    h_ret->print();
    cout << "done ranking" << endl;
    cout << "Time taken: " << (float)(stop - start)/CLOCKS_PER_SEC << endl;
    delete(h_mat);
    cudaFree(d_arr);
    delete(h_ret);
    
    return 0;
}
