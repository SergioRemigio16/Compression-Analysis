#include "Utilities.h"
#include <cmath>

std::vector<float> Utilities::createMatrixWave(int x, int y, int z, float frequency, float amplitude, float phase) {
    std::vector<float> matrix(x * y * z);
    for (int i = 0; i < x * y * z; i++) {
        float waveValue = amplitude * sin(frequency * i + phase);
        matrix[i] = waveValue;
    }
    return matrix;
}

std::vector<float> Utilities::createMatrixRandom(int x, int y, int z, float minVal, float maxVal) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(minVal, maxVal);

    std::vector<float> matrix(x * y * z);
    for (int i = 0; i < x * y * z; i++) {
        matrix[i] = dis(gen);
    }
    return matrix;
}

