#pragma once
class FFTAlgorithms
{
};



/*

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h>

// Define typedef for convenience
typedef std::vector<std::complex<double>> ComplexVec;

int main() {
    size_t x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    double compressionRatio = 0.1; // Compression ratio. 0.5 means keeping half of the frequency components.

    // Define original matrix
    double* originalMatrix = new double[x * y * z];

    // Initialize the original matrix with some data
    for (size_t i = 0; i < x * y * z; ++i) {
        originalMatrix[i] = i;
    }

    // Print original matrix
    std::cout << "Original Matrix:\n";
    for (size_t i = 0; i < x * y * z; ++i) {
        std::cout << originalMatrix[i] << ' ';
    }
    std::cout << "\n";

    // Perform FFT
    ComplexVec complexMatrix(x * y * z);
    fftw_plan p = fftw_plan_dft_r2c_3d(x, y, z, originalMatrix,
        reinterpret_cast<fftw_complex*>(complexMatrix.data()), FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    // Compress the data by zeroing out the least significant components
    size_t compressedSize = static_cast<size_t>(x * y * z * compressionRatio);
    for (size_t i = compressedSize; i < x * y * z; ++i) {
        complexMatrix[i] = 0;
    }

    // Perform inverse FFT
    double* decompressedMatrix = new double[x * y * z];
    fftw_plan q = fftw_plan_dft_c2r_3d(x, y, z, reinterpret_cast<fftw_complex*>(complexMatrix.data()),
        decompressedMatrix, FFTW_ESTIMATE);
    fftw_execute(q);
    fftw_destroy_plan(q);

    // FFTW's backwards transform does not normalize the result, so we have to do it manually
    for (size_t i = 0; i < x * y * z; ++i) {
        decompressedMatrix[i] /= (x * y * z);
    }

    // Print decompressed matrix
    std::cout << "Decompressed Matrix:\n";
    for (size_t i = 0; i < x * y * z; ++i) {
        std::cout << decompressedMatrix[i] << ' ';
    }
    std::cout << "\n";

    // Print size of original and compressed matrices in bytes
    std::cout << "Size of original matrix (bytes): " << (x * y * z) * sizeof(double) << "\n";
    std::cout << "Size of compressed matrix (bytes): " << compressedSize * sizeof(std::complex<double>) << "\n";

    // Compute and print Mean Squared Error
    double mse = 0;
    for (size_t i = 0; i < x * y * z; ++i) {
        double error = originalMatrix[i] - decompressedMatrix[i];
        mse += error * error;
    }
    mse /= (x * y * z);
    std::cout << "Mean Squared Error: " << mse << "\n";

    // Deallocate memory
    delete[] originalMatrix;
    delete[] decompressedMatrix;

    return 0;
}
*/

