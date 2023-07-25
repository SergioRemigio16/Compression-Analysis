#include "TimingExperiment.h"
#include "ZFPAlgorithms.h"

#include <zfp.h>
#include <vector>
#include <complex>
#include <iostream>
#include <chrono>
#include <fftw3.h>
#include <numeric>
#include <cmath>

int main() {
    runTimingExperiment();
}


/*
int main() {
    int x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    double rate = 20; // Compression ratio. 0.9 means keeping 90% of the orignial data

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);

    // Compressing the data
    int byteStreamSize;
    unsigned char* byteStream = ZFPAlgorithms::compressMatrix(originalMatrix, x, y, z, rate, byteStreamSize);

    // Decompressing the received data
    double* decompressedMatrix = ZFPAlgorithms::decompressMatrix(byteStream, byteStreamSize);

    // Printing size of original and serialized compressed data
    std::cout << "Size of original matrix (bytes): " << (x * y * z) * sizeof(double) << "\n";
    std::cout << "Size of serialized compressed data (bytes): " << byteStreamSize << "\n";

    // Printing various types of error between original and decompressed data
    Utilities::printError(originalMatrix, decompressedMatrix, x, y, z);

    // Freeing the memory
    delete[] originalMatrix;
    delete[] decompressedMatrix;
    delete[] byteStream;

    return 0;
}
*/
