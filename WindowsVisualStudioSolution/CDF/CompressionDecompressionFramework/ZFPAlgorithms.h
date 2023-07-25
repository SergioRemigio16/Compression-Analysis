#ifndef _ZFPALGORITHMS_H_
#define _ZFPALGORITHMS_H_

#include <zfp.h>
#include <vector>
#include <complex>
#include <iostream>
#include <chrono>
#include <fftw3.h>
#include <numeric>
#include <cmath>

namespace ZFPAlgorithms {

    unsigned char* compressMatrix(double* originalData, int x, int y, int z, double rate, int& size);
    double* decompressMatrix(unsigned char* buffer, int bufferSize);
    /*
    double* decompressMatrixFixedRate(CompressionResult& result);
    double* decompressMatrixAccuracy(CompressionResult& result);
    */
}

#endif // _ZFPALGORITHMS_H_ 
