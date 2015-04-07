#ifndef INPUTOUTPUT_HPP
#define INPUTOUTPUT_HPP

#include"Matrix.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
using namespace std;

void loadMatrix(Matrix<float> **mat, 
        const char* filename, unsigned int nRows, unsigned int nCols)
{
    *mat = new Matrix<float>(nRows, nCols);
    ifstream fs;
    fs.open(filename);
    if (!fs.is_open()) {
        perror("File open failed.");
        exit(1);
    }
    string line;
    string token;
    int rowIdx = 0;
    int colIdx = 0;
    while(getline(fs, line)) {
        colIdx = 0;
        istringstream lineStream(line);
        while(getline(lineStream, token, ',')) {
            char *end;
            (*mat)->element(rowIdx, colIdx) = strtof(token.c_str(), &end);
            colIdx++;
        }
        rowIdx++;
    }
    fs.close();
}


#endif
