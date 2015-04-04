#include <iostream>
#include "Matrix.hpp"
int main()
{
    Matrix *mat = new Matrix(2000, 20000);
    mat->element(1,1) = 10.;
    if (mat->hasCorrelation(1,1)){
        std::cout << "has correlation" << std::endl;
    } else {
        std::cout << "no correlation" << std::endl;
    }
    std::cout << "data location = " << mat->memAddr() << std::endl;
    std::cout << "data size = " << mat->size() << std::endl;
    delete(mat);
    
    return 0;
}
