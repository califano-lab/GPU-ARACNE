#include <iostream>
#include "MiMatrix.hpp"
int main()
{
    MiMatrix<float> *mat = new MiMatrix<float>(2000, 20000);
    if (mat->hasCorrelation(1,1)){
        std::cout << "has correlation" << std::endl;
    } else {
        std::cout << "no correlation" << std::endl;
    }
    delete(mat);
    
    return 0;
}
