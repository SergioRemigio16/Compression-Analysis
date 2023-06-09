#ifndef _ZFPALGORITHMS_H_
#define _ZFPALGORITHMS_H_

#include <zfp.h>
#include <vector>
#include <complex>
#include <iostream>
#include <chrono>
#include <complex>
#include <zfp.h>
#include <fftw3.h>
#include <vector>
#include <numeric>
#include <cmath>

namespace ZFPAlgorithms {

    struct CompressionResult {
        int x;
        int y;
        int z;
        /// <summary>
        /// Holds the compressed data
        /// </summary>
        std::vector<unsigned char> buffer;
        size_t bufsize;
        double precision;
        uint rate;
        CompressionResult() : x(0), y(0), z(0), bufsize(0), rate(0), precision(0.0) {}
    };

    struct CompressionResultFFT {
        int x, y, z; // Original dimensions
        int compressSize; // Number of strongest frequencies we keep
        std::vector<int> indices; // Indices of the kept frequencies
        std::vector<std::complex<double>> frequencies; // Complex values of the kept frequencies
    };

    CompressionResult compressMatrixFixedRate(double*& originalData, int x, int y, int z, uint rate);
    CompressionResult compressMatrixAccuracy(double*& originalData, int x, int y, int z, double precision);
    double* decompressMatrixFixedRate(CompressionResult& result);
    double* decompressMatrixAccuracy(CompressionResult& result);

    CompressionResultFFT compressMatrixFFT(const double*& originalMatrix, int x, int y, int z, int compressSize);
    double* decompressMatrixFFT(const CompressionResultFFT& compressionResult);

    /// <summary>
	/// Calculate byte size of compression result
	/// </summary>
	/// <param name="compressionResult"></param>
	/// <returns></returns>
	size_t calculateDecompressedDataBytes(const ZFPAlgorithms::CompressionResult& compressionResult);
}

#endif // _ZFPALGORITHMS_H_ 
