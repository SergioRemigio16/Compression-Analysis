#include "Utilities.h"
#include <cmath>
#include <chrono>
#include <iostream>

double* Utilities::createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase) {
    double* matrix = new double[x * y * z];
    for (int i = 0; i < x * y * z; i++) {
        double waveValue = amplitude * sin(frequency * i + phase);
        matrix[i] = waveValue;
    }
    return matrix;
}

double* Utilities::createMatrixRandom(int x, int y, int z, double minVal, double maxVal) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(minVal, maxVal);

    double* matrix = new double[x * y * z];
    for (int i = 0; i < x * y * z; i++) {
        matrix[i] = dis(gen);
    }
    return matrix;
}

double Utilities::measureTimingOverhead() {
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> overhead = end - start;
    return overhead.count();
}

void Utilities::printMatrix(const double* matrix, int x, int y, int z) {
    for (int i = 0; i < z; i++) {
        std::cout << "z = " << i << ":\n";
        for (int j = 0; j < y; j++) {
            for (int k = 0; k < x; k++) {
                std::cout << matrix[i * x * y + j * x + k] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
}

double Utilities::calculateMSE(const double* originalData, const double* decompressedData, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        double diff = originalData[i] - decompressedData[i];
        sum += diff * diff;
    }
    return sum / size;
}
