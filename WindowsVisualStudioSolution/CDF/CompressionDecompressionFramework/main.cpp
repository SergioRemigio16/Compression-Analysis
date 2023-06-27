#include "TimingExperiment.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <algorithm>

void printError(const double* originalMatrix, const double* decompressedMatrix, int x, int y, int z) {
    // Compute and print Mean Squared Error
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



// Define typedef for convenience
typedef std::vector<std::complex<double>> ComplexVec;

// Define a pair to store magnitude and index
typedef std::pair<double, size_t> MagnitudeIndexPair;

// New struct to hold frequency component and its index
struct FrequencyComponent {
    std::complex<double> value;
    size_t index;
};

bool compareMagnitude(const MagnitudeIndexPair& a, const MagnitudeIndexPair& b) {
    return a.first < b.first;
}

int main() {
    size_t x = 3, y = 7, z = 7; // Dimensions of matrix. Modify as needed.
    double compressionRatio = 0.9; // Compression ratio. 0.5 means keeping half of the frequency components.

    // Define original matrix
    double* originalMatrix = Utilities::createMatrixWave(x, y, z, 1, 3.1415, 0, 1.0, 1.0, 7);
    // Perform FFT
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

    // Define your compressed data as a vector of FrequencyComponent
    std::vector<FrequencyComponent> compressedData;

    // Store only the significant frequency components in compressedData
    size_t significantSize = static_cast<size_t>((1.0 - compressionRatio) * x * y * z);
    for (size_t i = x * y * z - significantSize; i < x * y * z; ++i) {
        compressedData.push_back({ complexMatrix[magnitudes[i].second], magnitudes[i].second });
    }

    // To decompress, zero out complexMatrix and then fill in only the significant frequency components from compressedData
    std::fill(complexMatrix.begin(), complexMatrix.end(), 0.0);
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

    Utilities::printComparison(originalMatrix, decompressedMatrix, x, y, z);

    // Print size of original and compressed matrices in bytes
    std::cout << "Size of original matrix (bytes): " << (x * y * z) * sizeof(double) << "\n";
    std::cout << "Size of compressed matrix (bytes): " << compressedData.size() * sizeof(FrequencyComponent) << "\n";

    printError(originalMatrix, decompressedMatrix, x, y, z);

    // Deallocate memory
    delete[] originalMatrix;
    delete[] decompressedMatrix;

    return 0;
}





/*
int main() {
	return runTimingExperiment();
}
*/


