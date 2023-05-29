#include "Utilities.h"
#include <cmath>
#include <chrono>
#include <iostream>

// Generate matrix data based on a wave distribution.
// Input: The dimensions (x, y, z), frequency, amplitude and phase of the wave.
// Output: A vector representing the wave-distributed matrix.
std::vector<double> Utilities::createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase) {
    std::vector<double> matrix(x * y * z);
    for (int i = 0; i < x * y * z; i++) {
        double waveValue = amplitude * sin(frequency * i + phase);
        matrix[i] = waveValue;
    }
    return matrix;
}

// Generate matrix data based on a random distribution.
// Input: The dimensions (x, y, z), and minimum and maximum value of the random distribution.
// Output: A vector representing the randomly distributed matrix.
std::vector<double> Utilities::createMatrixRandom(int x, int y, int z, double minVal, double maxVal) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(minVal, maxVal);

    std::vector<double> matrix(x * y * z);
    for (int i = 0; i < x * y * z; i++) {
        matrix[i] = dis(gen);
    }
    return matrix;
}

// Measure the overhead of timing process.
// Output: The overhead of timing in seconds.
double Utilities::measureTimingOverhead() {
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> overhead = end - start;
    return overhead.count();
}

// Print the provided matrix to the console.
// Input: The matrix data, and its dimensions (x, y, z).
void Utilities::printMatrix(const std::vector<double>& matrix, int x, int y, int z) {
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
