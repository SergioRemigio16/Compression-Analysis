#include "Utilities.h"
#include <cmath>

std::vector<double> Utilities::createMatrixWave(int x, int y, int z, double frequency, double amplitude, double phase) {
    std::vector<double> matrix(x * y * z);
    for (int i = 0; i < x * y * z; i++) {
        double waveValue = amplitude * sin(frequency * i + phase);
        matrix[i] = waveValue;
    }
    return matrix;
}

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

