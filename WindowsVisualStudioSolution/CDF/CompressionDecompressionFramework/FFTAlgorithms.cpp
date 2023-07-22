#include "FFTAlgorithms.h"

#include <iostream>
#include <cmath>
#include <fftw3.h>
#include <algorithm>
#include <cstring>  

/*
int main() {
    size_t x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    double compressionRatio = 0.6; // Compression ratio. 0.5 means keeping half of the frequency components.

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);

    // Compressing the data
    std::vector<FrequencyComponent> compressedData = compressData(originalMatrix, x, y, z, compressionRatio);

    // Serializing the compressed data for sending via MPI

    size_t serializedSize;
    unsigned char* serializedCompressedData = serializeCompressedData(compressedData, serializedSize);

    // Receiving the serialized data via MPI and deserializing it
    // Note: In this code, we simply reuse the serialized data for demonstration purposes. In a real application,
    // the data would be sent and received via MPI.
    std::vector<FrequencyComponent> receivedCompressedData = deserializeCompressedData(serializedCompressedData, serializedSize);

    // Decompressing the received data
    double* decompressedMatrix = decompressData(receivedCompressedData, x, y, z);

    // Printing comparison of original and decompressed data
    Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);

    // Printing size of original and serialized compressed data
    std::cout << "Size of original matrix (bytes): " << (x * y * z) * sizeof(double) << "\n";
    std::cout << "Size of serialized compressed data (bytes): " << serializedSize << "\n";

    // Printing various types of error between original and decompressed data
    printError(originalMatrix, decompressedMatrix, x, y, z);

    // Freeing the memory
    delete[] serializedCompressedData;
    delete[] originalMatrix;
    delete[] decompressedMatrix;

    return 0;
}

*/


void printError(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z) {
	// Computing and printing Mean Squared Error, Mean Absolute Error,
    // Maximum Absolute Error, Root Mean Squared Error and Peak Signal to Noise Ratio.
    double mse = 0;
    for (size_t i = 0; i < x * y * z; ++i) {
        double error = originalMatrix[i] - decompressedMatrix[i];
        mse += error * error;
    }
    mse /= (x * y * z);

    std::cout << "Mean Squared Error: " << mse << "\n";
    // Compute and print Mean Absolute Error
    double mae = 0;
    for (size_t i = 0; i < x * y * z; ++i) {
        double error = std::abs(originalMatrix[i] - decompressedMatrix[i]);
        mae += error;
    }
    mae /= (x * y * z);
    std::cout << "Mean Absolute Error: " << mae << "\n";

    // Compute and print Maximum Absolute Error
    double maxError = 0;
    for (size_t i = 0; i < x * y * z; ++i) {
        double error = std::abs(originalMatrix[i] - decompressedMatrix[i]);
        maxError = std::max(maxError, error);
    }
    std::cout << "Max Absolute Error: " << maxError << "\n";

    // Compute and print Root Mean Squared Error
    double rmse = sqrt(mse);
    std::cout << "Root Mean Squared Error: " << rmse << "\n";

    // Compute and print Peak Signal to Noise Ratio (PSNR)
    double maxOriginalValue = *std::max_element(originalMatrix, originalMatrix + (x * y * z));
    double psnr = 20 * log10(maxOriginalValue / rmse);
    std::cout << "Peak Signal to Noise Ratio (in dB): " << psnr << "\n";

}

bool compareMagnitude(const MagnitudeIndexPair& a, const MagnitudeIndexPair& b) {
    return a.first < b.first;
}

unsigned char* serializeCompressedData(const std::vector<FrequencyComponent>& compressedData, size_t& size) {
	// Calculating the size of the data in bytes
    size = sizeof(FrequencyComponent) * compressedData.size();
    // Allocating an unsigned char buffer of sufficient size
    unsigned char* buffer = new unsigned char[size];
    // Copying the data from the vector to the buffer
    std::memcpy(buffer, compressedData.data(), size);
    // Returning the buffer
    return buffer;
}

std::vector<FrequencyComponent> deserializeCompressedData(const unsigned char* buffer, size_t size) {
    // Calculating the number of elements in the buffer
    std::vector<FrequencyComponent> compressedData(size / sizeof(FrequencyComponent));
    // Copying the data from the buffer to the vector
    std::memcpy(compressedData.data(), buffer, size);
    // Returning the vector
    return compressedData;
}

double* decompressData(const std::vector<FrequencyComponent>& compressedData, size_t x, size_t y, size_t z) {
    // Initialize complexMatrix with zeros
    ComplexVec complexMatrix(x * y * z, 0.0);

    // Fill in the significant frequency components from compressedData
    for (const auto& component : compressedData) {
        complexMatrix[component.index] = component.value;
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

    return decompressedMatrix;

}

std::vector<FrequencyComponent> compressData(double* originalMatrix, size_t x, size_t y, size_t z, double compressionRatio) {
    // Performing FFT on the original matrix
    ComplexVec complexMatrix(x * y * z);
    fftw_plan p = fftw_plan_dft_r2c_3d(x, y, z, originalMatrix,
        reinterpret_cast<fftw_complex*>(complexMatrix.data()), FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

    // Get magnitudes and corresponding indices
    std::vector<MagnitudeIndexPair> magnitudes;
    for (size_t i = 0; i < x * y * z; ++i) {
        magnitudes.push_back(MagnitudeIndexPair(std::abs(complexMatrix[i]), i));
    }

    // Sort the magnitudes
    std::sort(magnitudes.begin(), magnitudes.end(), compareMagnitude);

    // Store only the significant frequency components in compressedData
    size_t significantSize = static_cast<size_t>((1.0 - compressionRatio) * x * y * z);
    std::vector<FrequencyComponent> compressedData;
    for (size_t i = x * y * z - significantSize; i < x * y * z; ++i) {
        compressedData.push_back({ complexMatrix[magnitudes[i].second], magnitudes[i].second });
    }

    return compressedData;
}



