#include <iostream>
#include <cstdlib>
#include <ctime>
//#include "Matrix.hpp"
//#include "genMiCutoff.hpp"
#include "genMiCutoff.hpp"
#include "genMiCutoff.cpp"
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 3) {
        perror("sample number and permutation number!");
        exit(1);
    }
    const int Nsmp = atoi(argv[1]);
    const int Nperm = atoi(argv[2]); // this will change to fix value after testing  

    srand(time(NULL));

    cout << "start calculating null model mi .... " << endl;
    clock_t start = clock();

    float h_micut;
    h_micut  = calMIcutCoeff( Nsmp, Nperm ); 

    clock_t stop = clock();
    
    return 0;
}
